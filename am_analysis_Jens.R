library(plyr); library(boot); library(ggplot2); library(gridExtra)

#Temporary Change 1

####1. Read in data####
d=read.csv("df_export_18spp_v4.txt") #Takes ~1 minute
spp=read.csv("spgrp codes.csv")
####2. A little EDA (Optional)####
#What proportion of data points have mortality?
length(which(d$mort.bin==1))/length(d$mort.bin) #About 4%
#plot(Pnorm~Psd,data=d) #note the very strong correlation between mean and standard deviation precipitation
#plot(Tnorm~Tsd,data=d)

####3. Process data####
#Set up storage variables
p.list=list()
B.Pz1=data.frame("Estimate"=rep(NA,length(unique(d$spgroup))),"Std.Error"=NA,row.names=unique(d$spgroup))
for(i in unique(d$spgroup)){#Analyze every species separately
#spg=unique(d$spgroup)[1] #Pick your species group; can eventually be embedded in For loop.
d1=d[d$spgroup==i,] #Subset data to species group in question
d1$cell.id=paste(d1$pos.x,d1$pos.y,sep=".") #Identify the cell as unique combo of x and y position.
#hist(table(d1$cell.id)) #See how many cells have all 7 years
outliers=which(d1$mort.tpa > (mean(d1$mort.tpa)+sd(d1$mort.tpa)*100)) #Find outliers with really high mortality
#CHECKME: Is there a way to do this by looking for "natural breaks" in the tail of the distribution?
if(!!length(outliers)){# If there are outliers to remove, remove them
  d1=d1[-outliers,]}

which.cells=#This code chunk identifies those cells that have the full 7 years worth of data for the given species.
  dimnames(table(d1$cell.id))[[1]][which(as.vector(table(d1$cell.id))==7)]
#CHECKME: WHy is the length of which.cells so much smaller than what the histogram (above) says it should be?

d1= #Subset the data to those cells with all 7 years
  d1[d1$cell.id %in% which.cells,]

tmp= #Calculate mean newly-dead tpa for each cell across all 7 years.
  ddply(d1,.(cell.id),summarize,mort.mean=mean(mort.tpa))
d1=d1[order(d1$cell.id),] #Sort the data frame by cell.id
d1$mort.rel= #For each year, calculate mortality relative to the average mortality in that pixel (as a difference)
  d1$mort.tpa-tmp[rep(row.names(tmp),each=7),"mort.mean"]
#CHECKME: It's a problem that cells that had no mortality in 6 out of 7 years and lots of mortality in 1 out of 7 years both count as extreme examples, but count independently.

####4. Analyze mortality####
#d1$mort.rel is the y variable to use. It captures anomalously large mortality events for a given pixel. This accounts for the non-independence of repeated measures at a single pixel. The other thing we could add would be a spatial random effects term based on pos.x and pos.y, which would account for the non-independence of anomalously high mortality events during the same year at neighboring pixels (for example). In other words, this would account for the fact that there should be correlation in anomalously high mortality events among pixels that are closer together, rather than farther apart (because some mortality polygons contained several pixels within that one "mortality event")

#hist(d1$mort.rel,breaks=200)  
#hist(d1$mort.rel[-which(d1$mort.rel==0)],breaks=200)  

#Hypothesis 1a: Do drier sites have higher mean mortality across all years?
d1a= #Subset data into unique cells with 7 years worth of data, and calculate average mortality (tpa) per year
  ddply(d1,.(cell.id),summarize,
        mort.mean_tpa=mean(mort.tpa),mort.bin=max(mort.bin),
        Tnorm=mean(Tnorm),Pnorm=mean(Pnorm),
        Tsd=mean(Tsd),Psd=mean(Psd))
#hist(d1a$mort.mean)
m.bin <- glm(mort.bin ~ Tnorm * Pnorm, data = d1a, family = binomial(link = logit))
#summary(m.bin)
m.gam <- glm(mort.mean_tpa ~ Tnorm*Pnorm, data = subset(d1a, mort.bin == 1), family = Gamma(link = log))
#summary(m.gam)
d1a$mort.bin.pred=inv.logit(predict.glm(m.bin,newdata=d1a))
d1a$mort.cond_tpa.pred=d1a$mort.bin.pred*exp(predict.glm(m.gam,newdata=d1a))

#Plot Hypothesis 1a
ggplot(d1a) + 
  geom_point(aes(x=Tnorm,y=Pnorm,colour=mort.cond_tpa.pred)) + 
  scale_colour_gradientn(colours=rev(rainbow(3))) + 
  theme_bw()

#Hypothesis 1b: Do drier sites have higher mortality anomalies?
d1b= #Subset data into unique cells with 7 years worth of data, and calculate average mortality (tpa) per year
  ddply(d1,.(cell.id),summarize,
        mort.max_anom=max(mort.rel),mort.bin=max(mort.bin),
        Tnorm=mean(Tnorm),Pnorm=mean(Pnorm),
        Tsd=mean(Tsd),Psd=mean(Psd))
m.bin <- glm(mort.bin ~ Tnorm * Pnorm, data = d1b, family = binomial(link = logit))
#summary(m.bin)
m.gam <- glm(mort.max_anom ~ Tnorm*Pnorm, data = subset(d1b, mort.bin == 1), family = Gamma(link = log))
#summary(m.gam)
d1b$mort.bin.pred=inv.logit(predict.glm(m.bin,newdata=d1b))
d1b$mort.cond_tpa_anom.pred=d1b$mort.bin.pred*exp(predict.glm(m.gam,newdata=d1b))

#Plot Hypothesis 1b
# ggplot(d1b) + 
#   geom_point(aes(x=Tnorm,y=Pnorm,colour=mort.cond_tpa_anom.pred)) + 
#   scale_colour_gradientn(colours=rev(rainbow(3))) + 
#   theme_bw()
#Gives exactly the same results as Hypothesis 1a, because the cells with higher mortality anomalies also have higher mortality averaged across all years (which was tested in H1a)

#Hypothesis 1c: Which lag is the best predictor of mortality anomaly?
#Note: Using full dataset here; not using a hurdle model because distribution of mortality anomaly is continuous and contains both positive and negative values. It is still zero-inflated though.

m.0=lm(mort.rel~Tz0*Pz0, data=d1)
m.1=lm(mort.rel~Tz0*Pz0 + Tz1*Pz1, data=d1)
m.2=lm(mort.rel~Tz0*Pz0 + Tz1*Pz1 + Tz2*Pz2, data=d1)
m.3=lm(mort.rel~Tz0*Pz0 + Tz1*Pz1 + Tz2*Pz2 + Tz3*Pz3, data=d1)
m.4=lm(mort.rel~Tz0*Pz0 + Tz1*Pz1 + Tz2*Pz2 + Tz3*Pz3 + Tz4*Pz4, data=d1)
AIC(m.0,m.1,m.2,m.3,m.4)
#The more previous years you add, the better the fit gets, but the biggest gain in AIC is from just adding 1 year of lag, so let's go with m.1
m.1.only=lm(mort.rel~Tz1*Pz1,data=d1)
m.2.only=lm(mort.rel~Tz2*Pz2,data=d1)
m.3.only=lm(mort.rel~Tz3*Pz3,data=d1)
m.4.only=lm(mort.rel~Tz4*Pz4,data=d1)
AIC(m.0, m.1, m.1.only,m.2.only,m.3.only,m.4.only)
#Z scores from 1 year ago are a better predictor than z scores from the past year! Or than z scores from >2 years ago
#summary(m.0);summary(m.1);summary(m.1.only)
d1$mort.rel.pred=inv.logit(predict.glm(m.1.only,newdata=d1))
d1$mort.rel.pred=predict.glm(m.1.only,newdata=d1)
#CHECKME: I don't think I need inv.logit here
#CHECKME: Alternative to the above code is m.1.only$fitted (?)

#Plot Hypothesis 1c
#Plot for best model (1 year prior, T and P z scores)
spp_id=spp[spp$spgroup==i,"scientific.name"]
p1a=
  ggplot(d1) + 
  geom_point(aes(x=Tz1,y=Pz1,colour=mort.rel.pred)) + 
  scale_colour_gradientn(colours=rev(rainbow(3))) + 
  theme_bw()+labs(title="Model predictions from 1 year lag; \n showing last year's z scores")+
  annotate("text",label=as.character(spp_id),1.5,1.5)
p1b=
  ggplot(d1) + 
  geom_point(aes(x=Tz0,y=Pz0,colour=mort.rel.pred)) + 
  scale_colour_gradientn(colours=rev(rainbow(3))) + 
  theme_bw()+labs(title="Model predictions from 1 year lag; \n showing this year's z scores")
p.list[[i]]=grid.arrange(p1a,p1b,ncol=1)
B.Pz1[as.character(i),c(1,2)]=summary(m.1.only)$coefficients[3,c(1,2)]
}

#CHECKME 1:CHeck error messages- why are we getting fitted probabilities 0 or 1? Why is that a problem?
#CHECkME 2: Can we build in spatial random effects? How? SPAMM?. DO THIS FIRST, MIGHT FIX OTHER PROBLEMS
#CHECKME 3: Look at deficit, aet, and snowpack as additional predictore (/alternative predictors)

B.Pz1$Spgrp=rownames(B.Pz1)
#Precip effects, species ranks
#B.Pz1=B.Pz1[rev(order(B.Pz1$Estimate)),]
B.Pz1$Spgrp=factor(B.Pz1$Spgrp)
B.Pz1$Spgrp=reorder(B.Pz1$Spgrp,-B.Pz1$Estimate)
ggplot(B.Pz1, aes(Spgrp,Estimate)) + 
geom_point() + 
geom_errorbar(aes(ymin=Estimate-2*Std.Error,ymax=Estimate+2*Std.Error))+
theme_bw()
