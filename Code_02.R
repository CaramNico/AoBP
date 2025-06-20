library(nlme);library(MuMIn); library(ggeffects); library(nlme);library(glmmTMB)
library(nortest);library(DHARMa)library(ggplot2);library(lmerTest)
# ========= STUDY 1 =================
rm(list=ls())

# Read excel from Github. Sheet = "Study 1 RUE" ========
x$year = as.factor(x$year);x$subplot = as.factor(x$subplot)
x$date = as.factor(x$date)
x$composition = as.factor(x$composition);x$group_rp = as.factor(x$group_rp)
x$plot = as.factor(x$plot);x$sub = as.factor(x$sub);
x$treatment = as.factor(x$treatment);x$trt = as.factor(x$trt)
x$suff = as.factor(x$suff)

# ==== ANOVA ========
#x <- x[, c(1:13,22)]
x<-na.omit(x)
str(x)
fit <- lme(rue ~ treatment*composition*year,
           random = ~ 1 | plot/sub,
           weights = varIdent(form = ~ 1 | treatment*year),
           correlation = corCAR1(form = ~ clipping | plot/sub),
           data = x)

AIC(fit)
plot(fit);summary(fit)
anova(fit)
#corMatrix(fit$modelStruct$corStruct)
cld(lsmeans(fit, ~ composition), adjust="None")
emmeans(fit,pairwise ~ composition, adjust = "none")


pred <- ggpredict(fit, terms = c("treatment", "composition", "year")) 
rue= plot(pred,dot.size = 2, show.title = F, colors = c("blue", "red"))+theme_bw()+
  labs(x = "Treatment",
       y = "RUE (g MJ-1)")
rue

# ======== Linear and Non-linear regressions =========
year1=x[x$year=="2021",]
year2=x[x$year=="2022",]

lp_95 <- nlrq(rue ~ SSplin(rp*100, Asym, mid, scal), 
              data = x, tau = 0.95);summary(lp_95);AIC(lp_95)
lp_95_Y1 <- nlrq(rue ~ SSplin(rp*100, Asym, mid, scal), 
                 data = year1, tau = 0.95);summary(lp_95_Y1);AIC(lp_50_Y1)
lp_95_Y2 <- nlrq(rue ~ SSplin(rp*100, Asym, mid, scal), 
                 data = year2, tau = 0.95);summary(lp_95_Y2);AIC(lp_95_Y2)
lp_50_mod2 <- nlsLM(rue ~ SSplin(rp*100, a, xs, b), data = x, 
                    start = list(a = 0.67462, xs = 0.31000, b = -0.87393 ))
summary(lp_50_mod2);AIC(lp_50_mod2)

lp_50_Y1 <- nlsLM(rue ~ SSplin(rp*100, a, xs, b), data = year1, 
                  start = list(a = 0.67462, xs = 0.31000, b = -0.87393 ))
summary(lp_50_Y1);AIC(lp_50_Y1)

lp_50_Y2 <- nlsLM(rue ~ SSplin(rp*100, a, xs, b), data = year2, 
                  start = list(a = 0.67462, xs = 0.31000, b = -0.87393 ))

summary(lp_50_Y2);AIC(lp_50_Y2)
year2$rp <-year2$rp*100
lm_10 <- rq(rue ~ rp, data = x, tau = 0.1);summary(lm_10);AIC(lm_10)
lm_10_Y1 <- rq(rue ~ rp, data = year1, tau = 0.1);summary(lm_10_Y1);AIC(lm_10_Y1)
lm_10_Y2 <- rq(rue ~ rp, data = year2, tau = 0.1);summary(lm_10_Y2);AIC(lm_10_Y2)

plot1 <- ggplot(data = x, aes(x = rp*100, y = rue)) + 
  geom_point(aes(color=n_conc, shape = treatment)) + 
  scale_colour_gradient(low = "blue", high = "red")+
  geom_line(aes(y = fitted(lp_95)))+
  geom_line(data = df_line, aes(x = rp*100, y = predicted_lp), color = "grey",alpha=0.5) +
  geom_line(data = year2, aes(x = rp , y = predict(lp_95_Y2)), color = "grey", alpha=0.5,linetype = 'dashed') +
  geom_line(aes(y = fitted(lp_50_mod2)))+
  geom_line(data = df_line2, aes(x = rp*100, y = predicted_lp2), color = "grey",alpha=0.5) +
  geom_line(data = year2, aes(x = rp , y = predict(lp_50_Y2)), color = "grey",alpha=0.5,linetype = 'dashed') +
  geom_line(aes(y = fitted(lm_10)))+
  geom_line(data = year1, aes(x = rp*100 , y = predict(lm_10_Y1)), color = "grey",alpha=0.5) +
  geom_line(data = year2, aes(x = rp , y = predict(lm_10_Y2)), color = "grey", alpha=0.5, linetype = 'dashed') +
  #geom_line(aes(y = fitted(lp_80)), linetype = 'longdash')+
  theme_bw()+
  labs(x = "Legume proportion (%)",
       y = "RUE (g MJ-1)")
plot1

# Read excel from Github. Sheet = "Study 2 RUE" ========

str(x)
fit <- lme(rue ~ past*year,
           random = ~ 1 | clipping,
           weights = varIdent(form = ~ 1 | past*year),
           #correlation = corCAR1(form = ~ clipping | cage),
           data = x)
fit <- lm(rue ~ past*year, data = x);AIC(fit2)
AIC(fit);plot(fit);summary(fit)

means <- lsmeans(fit, ~ past);means
estimate=cld(means, adjust="None");estimate

pred <- ggpredict(fit, terms = c("year", "past")) 
rue= plot(pred,dot.size = 2, show.title = T, colors = c("blue", "red"))+theme_bw()+
  labs(x = "Year",
       y = "RUE (g MJ-1)")
rue

year1=x[x$year=="2022",]
year2=x[x$year=="2023",]

lm_95 <- rq(rue ~ rp + I(rp^2), data = x, tau = 0.95);summary(lm_95, se="boot");AIC(lm_95)
lm_95_Y1 <- rq(rue ~ rp + I(rp^2), data = year1, tau = 0.95);summary(lm_95_Y1, se="boot");AIC(lm_95_Y1)
#lm_95_Y1lm <- rq(rue ~ rp, data = year1, tau = 0.95);summary(lm_95_Y1lm, se="boot");AIC(lm_95_Y1lm)
lm_95_Y2 <- rq(rue ~ rp + I(rp^2), data = year2, tau = 0.95);summary(lm_95_Y2, se="boot");AIC(lm_95_Y2)


#======= 50th Quantile ==========

lm_50 <- rq(rue ~ rp + I(rp^2), data = x);summary(lm_50, se="boot");AIC(lm_50)
lm_50_Y1 <- rq(rue ~ rp + I(rp^2), data = year1);summary(lm_50_Y1, se="boot");AIC(lm_50_Y1)
#lm_95_Y1lm <- rq(rue ~ rp, data = year1, tau = 0.95);summary(lm_95_Y1lm, se="boot");AIC(lm_95_Y1lm)
lm_50_Y2 <- rq(rue ~ rp + I(rp^2), data = year2);summary(lm_50_Y2, se="boot");AIC(lm_50_Y2)



#======= 10th Quantile ==========

lm_10 <- rq(rue ~ rp + I(rp^2), data = x, tau = 0.1);summary(lm_10, se='boot');AIC(lm_10)
lm_10_Y1 <- rq(rue ~ rp, data = year1, tau = 0.1);summary(lm_10_Y1, se='boot');AIC(lm_10_Y1)
lm_10_Y2 <- rq(rue ~ rp + I(rp^2), data = year2, tau = 0.1);summary(lm_10_Y2, se='boot');AIC(lm_10_Y2)


# ========== Now plot ===============
library(ggplot2)
str(x)
plot1 <- ggplot(data = x, aes(x = rp, y = rue)) + 
  geom_point(aes(color=n_con, shape = year)) + 
  scale_colour_gradient(low = "blue", high = "red")+
  geom_line(aes(y = fitted(lm_95)))+
  geom_line(data = year1, aes(x = rp, y = predict(lm_95_Y1)), color = "grey", alpha=0.5) +
  geom_line(data = year2, aes(x = rp, y = predict(lm_95_Y2)), color = "grey", alpha=0.5,linetype = 'dashed') +
  geom_line(aes(y = fitted(lm_50)))+
  geom_line(data = year1, aes(x = rp, y = predict(lm_50_Y1)), color = "grey", alpha=0.5) +
  geom_line(data = year2, aes(x = rp, y = predict(lm_50_Y2)), color = "grey",alpha=0.5,linetype = 'dashed') +
  geom_line(aes(y = fitted(lm_10)))+
  geom_line(data = year1, aes(x = rp, y = predict(lm_10_Y1)), color = "grey",alpha=0.5) +
  geom_line(data = year2, aes(x = rp, y = predict(lm_10_Y2)), color = "grey", alpha=0.5, linetype = 'dashed') +
  #geom_line(aes(y = fitted(lp_80)), linetype = 'longdash')+
  theme_bw()+
  labs(x = "Legume proportion (%)",
       y = "RUE (g MJ-1)")
plot1

# NOW LOOCV ==============


# Initialize variables to store RMSE values and actual vs predicted values
rmse <- numeric()
actual_vs_predicted <- data.frame(Actual = numeric(), Predicted = numeric())
for (i in 1:nrow(x)) {
  test_data <- x[i, ]
  train_data <- x[-i, ]
  #model <- lm(rue ~ rp, data = train_data)
  #model <- lm(rue ~ rp + I(rp^2), data = train_data)
  #model <- nlsLM(rue ~ SSplin(rp_perc, Asym, mid, scal), data = train_data)
  model <- nlsLM(rue ~ SSpquad(rp_perc, Asym, mid, scal, c), data = train_data)
  predictions <- predict(model, newdata = test_data)
  rmse[i] <- sqrt(mean((test_data$rue - predictions)^2))
  actual_vs_predicted <- rbind(actual_vs_predicted, 
                               data.frame(Actual = test_data$rue, Predicted = predictions))
}

average_rmse <- round(mean(rmse), 4);average_rmse

# = Replace train and test data for each study, year and period =========
# = Same for the quantiles below 

# E.g., 95th quantile ========
for (i in 1:nrow(x)) {
  test_data <- x[i, ]
  train_data <- x[-i, ]
  #model <- rq(rue ~ rp_perc, data = train_data, tau=0.95)
  #model <- rq(rue ~ rp_perc + I(rp_perc^2), data = train_data, tau=0.95)
  #model <- nlrq(rue ~ SSplin(rp_perc, Asym, mid, scal), data = train_data, tau = 0.95)
  model <- nlrq(rue ~ SSpquad(rp_perc, Asym, mid, scal, c), data = x, tau = 0.95)
  predictions <- predict(model, newdata = test_data)
  rmse[i] <- sqrt(mean((test_data$rue - predictions)^2))
  actual_vs_predicted <- rbind(actual_vs_predicted, 
                               data.frame(Actual = test_data$rue, Predicted = predictions))
}

average_rmse <- round(mean(rmse), 4);average_rmse












