#==================== Light interception and LAI ==============
rm(list=ls(all=T))
library(nlme);library(multcomp);library(emmeans);library(ggpubr);library(ggeffects)
data <- read.csv(file="canopy.csv", header=T,sep=",",dec=".");dim(data);head(data)
str(data)
data <- data[, c(1:15)]
data$year<-as.factor(data$year);data$composition<-as.factor(data$composition)
data$plot<-as.factor(data$plot);data$sub<-as.factor(data$sub)
data$trt<-as.factor(data$trt);data$eu<-as.factor(data$eu)
data$eu_f<-as.factor(data$eu_f);data$eu_ff<-as.factor(data$eu_ff);data$period<-as.factor(data$period)
data <-na.omit(data)

fit <- nlme(
  li ~ 1 - exp(-b * lai),
  fixed = b ~ trt*composition,
  random = b ~ 1 | year/period/eu,
  correlation = corCompSymm(form = ~ time | year/period/eu),
  #correlation = corAR1(form = ~ time | eu),
  data = data,
  control = nlmeControl(msMaxIter = 1000),  # Increase maximum number of iterations
  start = c(b = rep(0.6, nlevels(data$trt)*nlevels(data$composition)))  # Provide an initial value for b for each Treatment level
)
summary(fit);plot(fit)
warp.lsm <- lsmeans(fit, ~ composition | trt, param = "b")
estimate=cld(warp.lsm, adjust="None")
summary(emmeans(fit, pairwise ~ composition | trt, param = "b"), infer=TRUE)


my_function <- function(x, b) {
  1 - exp(-b * x)
}
x1 <- seq(0, 10, length.out = 206);y1 <- my_function(x1, estimate[6,3]);y2 <- my_function(x1, estimate[2,3])
x2 <- seq(0, 10, length.out = 205);y3 <- my_function(x2, estimate[5,3]);y4 <- my_function(x2, estimate[3,3])
x3 <- seq(0, 10, length.out = 204);y5 <- my_function(x3, estimate[4,3]);y6 <- my_function(x3, estimate[1,3])

data$group <- paste(data$trt, data$composition, sep = "_")

plot_1=data[data$trt=="control",]
plot_2=data[data$trt=="cut",]
plot_3=data[data$trt=="cut+nutrients",]

plot1 <- ggplot(plot_1, aes(x = lai, y = li, color = factor(composition))) +
  geom_line(aes(x = x1, y = y2), color = "red",size=1) +
  geom_line(aes(x = x1, y = y1), color = "blue", size=1) +
  scale_color_manual(values = c("blue", "red")) +
  geom_point(aes(), alpha=0.2, show.legend=F) +  
  geom_hline(yintercept = .9, linetype = "dashed", color = "black", alpha=0.4)+
  labs(y = "Light interception (% incident PAR)",
       x = "Leaf area index")+
  theme_bw()+xlim(0, 8)+
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", alpha=0.4)

plot2 <- ggplot(plot_2, aes(x = lai, y = li, color = factor(composition))) +
  geom_line(aes(x = x2, y = y4), color = "red",size=1) +
  geom_line(aes(x = x2, y = y3), color = "blue", size=1) +
  scale_color_manual(values = c("blue", "red")) +
  geom_point(aes(), alpha=0.2, show.legend=F) +  
  geom_hline(yintercept = .9, linetype = "dashed", color = "black", alpha=0.4)+
  labs(y = "Light interception (% incident PAR)",
       x = "Leaf area index")+
  theme_bw()+xlim(0, 8)+
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", alpha=0.4)

plot3 <- ggplot(plot_3, aes(x = lai, y = li, color = factor(composition))) +
  geom_line(aes(x = x3, y = y6), color = "red",size=1) +
  geom_line(aes(x = x3, y = y5), color = "blue", size=1) +
  scale_color_manual(values = c("blue", "red")) +
  geom_point(aes(), alpha=0.2, show.legend=F) +  
  geom_hline(yintercept = .9, linetype = "dashed", color = "black", alpha=0.4)+
  labs(y = "Light interception (% incident PAR)",
       x = "Leaf area index")+
  theme_bw()+xlim(0, 8)+
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", alpha=0.4)

ggarrange(plot1, plot2, plot3, nrow = 1, align = "hv")

#============ NOW CANOPY ======================

rm(list=ls(all=T))
library(nlme);library(multcomp)
data <- read.csv(file="canopy.csv", header=T,sep=",",dec=".");dim(data);head(data)
str(data)
data <- data[, c(1:15)]
data$year<-as.factor(data$year);data$composition<-as.factor(data$composition)
data$plot<-as.factor(data$plot);data$sub<-as.factor(data$sub)
data$trt<-as.factor(data$trt);data$eu<-as.factor(data$eu)
data$eu_f<-as.factor(data$eu_f);data$eu_ff<-as.factor(data$eu_ff);data <-na.omit(data)
data$rpm_mm<-data$rpm_mm/10

fit0 <- lme(lai ~ composition*trt*rpm_mm,
            random = ~1 | year/period/eu,
            correlation = corCompSymm(form = ~ time | year/period/eu),
            data = data,
            control = lmeControl(msMaxIter = 10000)) 
summary(fit0)
anova(fit0)
cld(emtrends(fit0, pairwise~ trt * composition, var = "rpm_mm"), adjust='none')
emmeans(fit0, ~ trt * composition*rpm_mm, at = list(rpm_mm = 0))
summary(emtrends(fit0, pairwise ~ composition | trt, var = "rpm_mm"), infer=TRUE)

pr = ggpredict(fit0, terms = c("rpm_mm", "composition", "trt"))
biomass = plot(pr,add.data = T, jitter = T, dot.alpha = .4,
               dot.size = 1, show.title = F, colors = c("blue", "red")) +
  theme_bw()
labs(x = "Canopy height (cm)",
     y = "Leaf area index")+
  abline(h = 4, col = "lightgray", lty = 3)
biomass+xlim(0,25)+ylim(0,8)+
  geom_hline(yintercept = 4, linetype = "dashed", color = "black", alpha=0.4)+
  labs(x = "Canopy height (cm)",
       y = "Leaf area index")

#============ Now same canopy attributes for the grazing trial  ===========
rm(list=ls(all=T))
library(nlme);library(multcomp);library(emmeans);library(ggpubr);library(ggeffects)
data <- read.csv(file="grazing_lai.csv", header=T,sep=",",dec=".");dim(data);head(data)
str(data)
#data <- data[, c(1:15)]
data$year<-as.factor(data$year);data$period<-as.factor(data$period)
data$clipping<-as.factor(data$clipping);data$Plot<-as.factor(data$Plot)
data$pasture<-as.factor(data$pasture);data$Cage<-as.factor(data$Cage)
data$site<-as.factor(data$site)
data <-na.omit(data)

fit <- nlme(
  li ~ 1 - exp(-b * lai),
  fixed = b ~ pasture,
  random = b ~ 1 | year/clipping/site,
  correlation = corAR1(form = ~ time | year/clipping/site),
  #correlation = corAR1(form = ~ time | eu),
  data = data,
  control = nlmeControl(msMaxIter = 1000),  # Increase maximum number of iterations
  start = c(b = rep(0.6, nlevels(data$pasture)))  # Provide an initial value for b for each Treatment level
)
summary(fit);plot(fit)
warp.lsm <- lsmeans(fit, pairwise ~ pasture, param = "b")
estimate=cld(warp.lsm, adjust="None")


my_function <- function(x, b) {
  1 - exp(-b * x)
}
x1 <- seq(0, 10, length.out = 888);y1 <- my_function(x1, estimate[2,2]);y2 <- my_function(x1, estimate[1,2])


plot1 <- ggplot(data, aes(x = lai, y = li, color = factor(pasture))) +
  geom_point(aes(), alpha=0.2, show.legend=F) +  
  geom_line(aes(x = x1, y = y2), color = "red",size=1) +
  geom_line(aes(x = x1, y = y1), color = "blue", size=1) +
  scale_color_manual(values = c("blue", "red")) +
  
  geom_hline(yintercept = .9, linetype = "dashed", color = "black", alpha=0.4)+
  labs(y = "Light interception (% incident PAR)",
       x = "Leaf area index")+
  theme_bw()+xlim(0, 8)+
  geom_vline(xintercept = 3.8, linetype = "dashed", color = "black", alpha=0.4)
plot1


#============ NOW CANOPY ======================
str(data)
fit0 <- lme(lai ~ pasture*cm,
            random = ~1 | year/clipping/site,
            correlation = corCompSymm(form = ~ time | year/clipping/site),
            data = data,
            control = lmeControl(msMaxIter = 10000)) 
summary(fit0);plot(fit0)
anova(fit0)
cld(emtrends(fit0, pairwise ~ pasture|cm, var = "cm"), adjust='none')
emmeans(fit0, ~ pasture*cm, at = list(cm = 0))



pr = ggpredict(fit0, terms = c("cm", "pasture"))
biomass = plot(pr,add.data = T, jitter = T, dot.alpha = .4,
               dot.size = 1, show.title = F, colors = c("blue", "red")) +
  theme_bw()
labs(x = "Canopy height (cm)",
     y = "Leaf area index")+
  abline(h = 4, col = "lightgray", lty = 3)
biomass+xlim(0,25)+ylim(0,8)+
  geom_hline(yintercept = 3.8, linetype = "dashed", color = "black", alpha=0.4)+
  labs(x = "Canopy height (cm)",
       y = "Leaf area index")
