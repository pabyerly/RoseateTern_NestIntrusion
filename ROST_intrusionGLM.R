library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(car)
library(AICcmodavg)
library(effects)
library(lattice)  
library(agricolae)

rost=read.csv("rost_intrusion_site.csv")
rost_ind=read.csv("rost_intrusion_ind.csv")
rost_master=read.csv("ROST_MasterData_Field.csv")

#nest success by island and year
success_summary=rost_master%>%
  group_by(island, year) %>%
  summarise(mean=mean(Fledge, na.rm=TRUE),
            SD=sd(Fledge, na.rm=TRUE), SE=SD/(sqrt(118)))

#predation rate by island and year
pred_summary=rost_master%>%
  group_by(island, year) %>%
  summarise(mean=mean(pred, na.rm=TRUE),
            SD=sd(pred, na.rm=TRUE), SE=SD/(sqrt(118)))

#intrusion rate by island and year
int_summary_year=rost%>%
  group_by(Island, Year) %>%
  summarise(int_mean=mean(n_intrusions),
            SD=sd(n_intrusions, na.rm=TRUE), SE=SD/(sqrt(118)))
            

#mean intrusion rate by island
int_summary=rost%>%
  group_by(Island) %>%
  summarise(int_mean=mean(intrusion_rate),
            SE=((sd(intrusion_rate, na.rm=TRUE))/sqrt(118)))

##########################################################################################

#effects of intrusion on predation and nest survival 
#1. model effect of intrusion rate on nest survival
mod1=lmer(mean_success ~ 1+(1|Island), data=rost)
mod2=lmer(mean_success ~ early_intrusion+(1|Island), data=rost)
anova(mod1, mod2)
delta=AICc(mod1, return.K = FALSE)-AICc(mod2, return.K = FALSE)
delta
#model output from lmertest: for mod 2, intrusion rate
contest1D(mod2, c(0, 1), confint=TRUE)

#2. model effect of intrusion rate on nest survival
mod1=lmer(mean_pred ~ 1+(1|Island), data=rost)
mod2=lmer(mean_pred ~ early_intrusion+(1|Island), data=rost)
anova(mod1, mod2)
AICc(mod1, return.K = FALSE)
AICc(mod2, return.K = FALSE)
delta=AICc(mod1, return.K = FALSE)-AICc(mod2, return.K = FALSE)
delta
#model output from lmertest: for mod 2, intrusion rate
contest1D(mod2, c(0, 1), confint=TRUE)

####################################################################################################################

#mixed effects models 

#scale large covariate
size <- scale(rost$colony_n, center = TRUE, scale = TRUE)
#rename others for model fit
n_intrusions <- rost$n_intrusions
place <- rost$place
dens<-rost$dens
height<-rost$height_cover
cover<-rost$per_cover
lagu<-rost$n_lagu
tern<-rost$n_terns

full=glmer(n_intrusions~size+place+dens+height+cover+lagu+tern+(1|Island)+offset(log(hours)), family="poisson", data=rost)
car::vif(full)

#intrusion rate 
#create model list 
mod=list()
#null and full model 
mod[[1]] <-glmer(n_intrusions~1+(1|Island)+offset(log(hours)), family="poisson", data=rost)
mod[[2]] <-glmer(n_intrusions~size+place+dens+height+cover+lagu+tern+(1|Island)+offset(log(hours)), family="poisson", data=rost)
#hypothesis: colony size only
mod[[3]] <-glmer(n_intrusions~size+(1|Island)+offset(log(hours)), family="poisson", data=rost)
#hypotheses: biological (conspecifics, preds, other terns)
mod[[4]] <-glmer(n_intrusions~lagu+(1|Island)+offset(log(hours)), family="poisson", data=rost)
mod[[5]] <-glmer(n_intrusions~size+lagu+(1|Island)+offset(log(hours)), family="poisson", data=rost)
mod[[6]] <-glmer(n_intrusions~lagu+tern+(1|Island)+offset(log(hours)), family="poisson", data=rost)
#hypotheses: site covariates
mod[[7]] <-glmer(n_intrusions~place+dens+(1|Island)+offset(log(hours)), family="poisson", data=rost)
mod[[8]] <-glmer(n_intrusions~cover+height+place+(1|Island)+offset(log(hours)), family="poisson", data=rost)
mod[[9]] <-glmer(n_intrusions~cover+height+dens+(1|Island)+offset(log(hours)), family="poisson", data=rost)
#hypotheses: biological + site covariates + colony size 
mod[[10]] <-glmer(n_intrusions~size+place+dens+(1|Island)+offset(log(hours)), family="poisson", data=rost)
mod[[11]] <-glmer(n_intrusions~size+cover+dens+(1|Island)+offset(log(hours)), family="poisson", data=rost)
mod[[12]] <-glmer(n_intrusions~size+lagu+dens+(1|Island)+offset(log(hours)), family="poisson", data=rost)
mod[[13]] <-glmer(n_intrusions~size+lagu+tern+dens+(1|Island)+offset(log(hours)), family="poisson", data=rost)

#create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)
#top model 
mod <-glmer(n_intrusions~colony_n+cover+dens+(1|Island), family="poisson", data=rost)

#effects: https://strengejacke.github.io/ggeffects/
mod_predict1 = ggpredict(mod, terms="colony_n")
plot1=ggplot(mod_predict, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0, 8) +
  labs(title ="a.", 
       x = "Colony size", y = "Intrusions/Hour") + 
  pubtheme
plot1

mod_predict2 = ggpredict(mod, terms="cover")
plot2=ggplot(mod_predict2, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0, 8) +
  labs(title ="b.", 
       x = "% Nest Cover", y = "Intrusions/Hour") + 
  pubtheme
plot2

mod_predict3 = ggpredict(mod, terms="dens")
plot3=ggplot(mod_predict3, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0, 50) +
  labs(title ="c.", 
       x = "Nest Site Density", y = "Intrusions/Hour") + 
  pubtheme
plot3

#########################################################################################################

#individual behavior at nest
rost_ind=read.csv("rost_intrusion_ind.csv")
rost=read.csv("rost_intrusion_site.csv")

rost_nest=merge(rost_ind, rost, by= c("camera"))

#mean of behavior: stays
mean(rost_nest$stay.defend)
#standard deviation (binomial formula)
sd=sqrt((1-0.38)^2/141)

#mean of behavior: leaves
mean(rost_nest$leaves)
#standard error
sd=sqrt((1-0.62)^2/141)


#summarize by pred type 
defense=rost_nest%>%
  group_by(pred_type) %>%
  summarise(stay.defend=mean(stay.defend),
            SD=sd(stay.defend))

#summarize by pred type 
fleeing=rost_nest%>%
  group_by(pred_type) %>%
  summarise(leaves=mean(leaves),
            SD=sd(leaves))

response=as.factor(rost_nest$response)
#INDIVIDUAL BEHAVIOR: social effects on leaves nest when intruder comes 
#create model list
defends=rost_nest$stay.defend
#scale large covariates
size <- scale(rost_nest$colony_n, center = TRUE, scale = TRUE)
julian <- scale(rost_nest$julian, center = TRUE, scale = TRUE)
n_intrusions <- rost_nest$n_intrusions
place <- rost_nest$place
dens<-rost_nest$dens
height<-rost_nest$height_cover
cover<-rost_nest$per_cover
lagu<-rost_nest$n_lagu
tern<-rost_nest$n_terns

#defends/stays
#covariates: n, place, dens, lagu, tern, period, cover
mod=list()
#null and full model 
mod[[1]] <-glmer(defends~1+(1|nest), family="binomial", data=rost_nest)
mod[[2]] <-glmer(defends~size+place+dens+lagu+tern+period+(1|nest), family="binomial", data=rost_nest)
mod[[3]] <-glmer(defends~size+(1|nest), family="binomial", data=rost_nest)
mod[[4]] <-glmer(defends~size+lagu+tern+dens+(1|nest)+(1|Island), family="binomial", data=rost_nest)
mod[[5]] <-glmer(defends~place+dens+(1|nest), family="binomial", data=rost_nest)
mod[[6]] <-glmer(defends~dens+(1|nest), family="binomial", data=rost_nest)
mod[[7]] <-glmer(defends~size+lagu+tern+(1|nest), family="binomial", data=rost_nest)
mod[[8]] <-glmer(defends~period+(1|nest), family="binomial", data=rost_nest)
mod[[9]] <-glmer(defends~size+place+(1|nest), family="binomial", data=rost_nest)
mod[[10]] <-glmer(defends~size+place+period+(1|nest), family="binomial", data=rost_nest)
mod[[11]] <-glmer(defends~size+lagu+tern+period+(1|nest), family="binomial", data=rost_nest)
mod[[12]] <-glmer(defends~size+dens+(1|nest), family="binomial", data=rost_nest)

#create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

mod2 <-glmer(defends~size+dens+(1|nest), family="binomial", data=rost_nest)
mod_predict4 = ggpredict(mod2, terms="size")
plot4=ggplot(mod_predict4, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0, 50) +
  labs(title ="c.", 
       x = "Colony Size", y = "Stays/Defends") + 
  pubtheme
plot4

#leaves
mod=list()
#null and full model 
mod[[1]] <-glmer(leaves~1+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[2]] <-glmer(leaves~size+place+dens+lagu+tern+time_light+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[3]] <-glmer(leaves~size+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[4]] <-glmer(leaves~size+lagu+tern+dens+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[5]] <-glmer(leaves~place+dens+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[6]] <-glmer(leaves~dens+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[7]] <-glmer(leaves~size+lagu+tern+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[8]] <-glmer(leaves~period+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[9]] <-glmer(leaves~size+place+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[10]] <-glmer(leaves~size+place+period+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[11]] <-glmer(leaves~size+lagu+tern+period+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[12]] <-glmer(leaves~size+tern+dens+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod[[13]] <-glmer(leaves~size+dens+time_light+(1|nest), family=binomial(link="logit"), data=rost_nest)

#create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(mod), sep = " ")
##generate AICc table
aictab(cand.set = mod, modnames = Modnames, sort = TRUE)
#round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = mod, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

mod3 <-glmer(leaves~colony_n+tern+dens+(1|nest), family=binomial(link="logit"), data=rost_nest)
mod_predict5 = ggpredict(mod3, terms="colony_n")
plot4=ggplot(mod_predict5, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0, 2) +
  labs(title ="a.", 
       x = "Colony Size", y = "Leaves") + 
  pubtheme
plot4

mod_predict6 = ggpredict(mod3, terms="dens")
plot5=ggplot(mod_predict6, aes(x, predicted)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  ylim(0, 2) +
  labs(title ="b.", 
       x = "Nest Site Density", y = "Leaves") + 
  pubtheme
plot5
