# Libraries
library(mgcv)
library(fields)
library(dplyr)
library(glmnet)
library(MuMIn)
library(lawstat)
library(car)
library(ggplot2)

# Read File
covid_cases <- read.csv("CovidCases.csv")

head(covid_cases)
# ------------------------------------------------------------------------------------------------------------------------
# Q1
# Model a Poisson Initially
poi_mdl <- glm(Deaths ~ . -Country -Confirmed, offset=log(Confirmed),
              data=covid_cases, family=poisson)
summary(poi_mdl)

# Checking for Overdispersion by modelling a quasipoisson
q_poi_mdl <-  glm(Deaths ~ . -Country -Confirmed, offset=log(Confirmed),
              data=covid_cases, family=quasipoisson)
summary(q_poi_mdl)

confint(q_poi_mdl)
# ------------------------------------------------------------------------------------------------------------------------
# Q2
# filter the data to be 10 or more deaths
filtered_data <- covid_cases %>% filter(Deaths >= 10)
filtered_data

# refit the quasi poisson on the filtered data
filt_poi_mdl <- glm(Deaths ~ . -Country - Confirmed,
                offset=log(Confirmed), data=filtered_data, family=quasipoisson)
summary(filt_poi_mdl)
# ------------------------------------------------------------------------------------------------------------------------
# Q3
# Take the covariates to check correlation against
covariates <- c('PopDensity',
                'MedianAge','UrbanPop','Bed','Lung','HealthExp','GDP')
# Correlation
pairs(subset(filtered_data, select=covariates),
      upper.panel=NULL, pch=19, cex=0.3)

# Residual Plots
residualPlots(filt_poi_mdl,
              type="pearson",
              smooth=list(smoother=gamLine, col="#377eb8"),
              fitted=FALSE,
              col="grey",
              pch=19,
              cex=0.3)

# ACF autocorrelation plots
acf(residuals(filt_poi_mdl, type="pearson"), main="filt_poi_mdl")
# Runs test
runs.test(residuals(filt_poi_mdl, type="pearson"))
# Overall Residual plots
residualPlots(filt_poi_mdl, ~1, fitted=TRUE)

# mean Variance relationship
plot(fitted(mdl), residuals(mdl, type="pearson"), xlab="Mean fitted values",
     ylab="Variance of Pearson residuals", pch=19)
# ------------------------------------------------------------------------------------------------------------------------
# Q4
# Making a poisson model to get Liklihood
filt_poi_mdl_normal <- glm(Deaths ~ . -Country -Confirmed, offset=log(Confirmed), data=filtered_data, family=poisson)

options(na.action = "na.fail")
dre = dredge(filt_poi_mdl_normal, rank="QAIC", chat=summary(filt_poi_mdl)$dispersion)
head(dre,5)
# ------------------------------------------------------------------------------------------------------------------------
# Q5 Dredge with Seed
set.seed(42)
xmatrix <- model.matrix( ~ scale(PopDensity) +scale(MedianAge)
              +scale(UrbanPop) +scale(Bed)+ scale(Lung)+ scale(HealthExp)+ scale(GDP) -1 , data=filtered_data)

head(xmatrix)

# Calculate Lasso and Cross validation Lasso
LASSO <- glmnet(xmatrix, filtered_data$Deaths, family="poisson",
                offset=log(filtered_data$Confirmed), alpha=1)
cvLASSO <- cv.glmnet(xmatrix, filtered_data$Deaths,
  family="poisson", offset=log(filtered_data$Confirmed), alpha=1, nfolds=10)

# Get the minimum
log(cvLASSO$lambda.min)
# Plot
par(mfrow=c(1, 2))
plot(LASSO, xvar="lambda")
abline(v=log(cvLASSO$lambda.min), lwd=4, lty=2)
plot(cvLASSO)
abline(v=log(cvLASSO$lambda.min), lwd=4, lty=2)
# Get Coefficients
coef(cvLASSO, cvLASSO$lambda.min)
# ------------------------------------------------------------------------------------------------------------------------
# Q6
# Model for k=5
k<- 5
PRS_model_5 <- mgcv::gam(Deaths ~ s(PopDensity, k=k) +s(MedianAge, k=k)
              +s(UrbanPop, k=k) +s(Bed, k=k)+ s(Lung, k=k)+ s(HealthExp, k=k)+ s(GDP, k=k),
              family=quasipoisson, data=filtered_data, offset=log(Confirmed))

summary(PRS_model_5)
PRS_model_5$sp

# Model for k=10
k<- 10
PRS_model_10 <- mgcv::gam(Deaths ~ s(PopDensity, k=k) +s(MedianAge, k=k)
              +s(UrbanPop, k=k) +s(Bed, k=k)+ s(Lung, k=k)+ s(HealthExp, k=k)+ s(GDP, k=k),
              family=quasipoisson, data=filtered_data, offset=log(Confirmed))

summary(PRS_model_10)
PRS_model_10$sp

# Plot Partial Residuals for k=5 and k=10
par(mfrow=c(3,3))
plot(PRS_model_5, shade=T, residuals=T, ylim=c(-10,10))

par(mfrow=c(3,3))
plot(PRS_model_10, shade=T, residuals=T, ylim=c(-10,10))
# ---------------------------------------------------------------------------------------------------------------
# Q7
covid_time = read.csv("CovidConfirmedTime.csv")
covid_time

# reorder data by maximum day value for each country
max_data <- covid_time %>% group_by(Country) %>% top_n(1, Day) %>% arrange(desc(Confirmed))
max_data

# Plot bar chart
ggplot(max_data, aes(x=reorder(Country, Confirmed), y=Confirmed)) +
  geom_bar(stat = "identity", fill="steelblue") +
  labs(y = "Confirmed Count", x = "Country") +
   ggtitle("Maximum Confirmed Count") +
   theme(plot.title = element_text(hjust = 0.5)) + coord_flip()

# Plot line chart
list_of_countries <- c("United Kingdom", "Korea, South",
                        "US", "Italy", "Norway", "Austria", "Spain", "Sweden",
                      "Portugal", "Germany", "Australia", "Iran")

subset(covid_time, Country %in% list_of_countries) %>%
  ggplot( aes(x=Day, y=Confirmed, group=Country, color=Country)) +
    geom_line()
# ------------------------------------------------------------------------------------------------------------------------
# Q9 Plot
average_country <- aggregate(covid_time$Confirmed, list(covid_time$Day), mean)
names(average_country) <- c("Day", "Confirmed")

ggplot()+ geom_line(data=subset(covid_time, Country %in% c("Germany", "United Kingdom")),
                    mapping=aes(x=Day, y =Confirmed, group=Country, color=Country)) +
                    geom_line(data = average_country, mapping=aes(x =Day,  y = Confirmed, col = "Average Country"))

# Get the fatlity rate for countries
ft_rate_country <- covid_cases %>% filter(Country %in% c("Germany", "United Kingdom"))
ft_rate_country$rate = ft_rate_country$Deaths / ft_rate_country$Confirmed
ft_rate_country

# ------------------------------------------------------------------------------------------------------------------------
# Q12
# Reading in the fixed and random effects tables
random_data <- read.csv("Random.csv")
fixed_data <- read.csv("Fixed.csv")

# Examing parts of the table
random_data %>% filter(Effect== "Day")
head(random_data %>% filter(Effect== "Intercept") %>% arrange(desc(abs(Estimate))),10)
head(random_data %>% filter(Effect== "Day") %>% arrange(desc(abs(Estimate))),10)


# Plot the histogram of slopes with the fixed effect added.
ggplot(data=subset(random_data, Effect %in% c("Day")) ,
          aes(x=Estimate+subset(fixed_data, Effect %in% c("Day"))$Estimate)) +
          geom_histogram(color="black",fill="white", bins=30)+
          labs(y = "Frequency", x = "Log of Slope")

# Make the predictions from the data by building intercepts and days
# for 49 days
days <- seq(0,49)
# Get the top 3 countries by absolute slope deviation
top3 <- head(random_data %>% filter(Effect %in% c("Day") ) %>% arrange(desc(abs(Estimate))),3)
top3
# Generate the "average country" data
norm_int <- subset(fixed_data, Effect == "Intercept")$Estimate
norm_day <- subset(fixed_data, Effect == "Day")$Estimate
norm_values <- days * norm_day + norm_int

# generate data for 3 countries
list_obj <- list()
c_names <- c()
for (i in top3$Subject){
  c_names <- c(c_names, i)
  out_list <- c()
  intercept <- norm_int + random_data %>% filter(Subject == i) %>% subset(Effect == "Intercept") %>% select(Estimate)
  day_value <- norm_day + random_data %>% filter(Subject == i) %>% subset(Effect == "Day") %>% select(Estimate)
  values <- days * day_value$Estimate + intercept$Estimate

  df <- data.frame(cbind(days, values))
  names(df) <- c("days", "values")
  list_obj[[i]] <- df
}

list_obj
c_names

# PLot on the normal scale
ggplot() + geom_line(data=data.frame(cbind(days, norm_values)),
                    mapping=aes(x=days, y=exp(norm_values), color="Average Country")) +
          geom_line(data=list_obj[[c_names[1]]],
                              mapping=aes(x=days, y=exp(values), color=paste(c_names[1], "Fitted", sep=" "))) +
          geom_line(data=list_obj[[c_names[2]]],
                              mapping=aes(x=days, y=exp(values), color=paste(c_names[2], "Fitted", sep=" "))) +
          geom_line(data=list_obj[[c_names[3]]],
                              mapping=aes(x=days, y=exp(values), color=paste(c_names[3], "Fitted", sep=" "))) +
          geom_line(data=subset(covid_time, Country == "Turkey"),
                              mapping=aes(x=Day, y=Confirmed, color=paste(c_names[1], "Observed", sep=" "))) +
          geom_line(data=subset(covid_time, Country == "Korea, South"),
                              mapping=aes(x=Day, y=Confirmed, color=paste(c_names[2], "Observed", sep=" "))) +
          geom_line(data=subset(covid_time, Country == "US"),
                              mapping=aes(x=Day, y=Confirmed, color=paste(c_names[3], "Observed", sep=" "))) +
                              labs(y = "Confirmed Cases", x = "Days", color="Country")

# Plot on the Log Scale
ggplot() + geom_line(data=data.frame(cbind(days, norm_values)),
                    mapping=aes(x=days, y=(norm_values), color="Average Country")) +
          geom_line(data=list_obj[[c_names[1]]],
                              mapping=aes(x=days, y=(values), color=paste(c_names[1], "Fitted", sep=" "))) +
          geom_line(data=list_obj[[c_names[2]]],
                              mapping=aes(x=days, y=(values), color=paste(c_names[2], "Fitted", sep=" "))) +
          geom_line(data=list_obj[[c_names[3]]],
                              mapping=aes(x=days, y=(values), color=paste(c_names[3], "Fitted", sep=" "))) +
          geom_line(data=subset(covid_time, Country == "Turkey"),
                              mapping=aes(x=Day, y=log(Confirmed), color=paste(c_names[1], "Observed", sep=" "))) +
          geom_line(data=subset(covid_time, Country == "Korea, South"),
                              mapping=aes(x=Day, y=log(Confirmed), color=paste(c_names[2], "Observed", sep=" "))) +
          geom_line(data=subset(covid_time, Country == "US"),
                              mapping=aes(x=Day, y=log(Confirmed), color=paste(c_names[3], "Observed", sep=" "))) +
                              labs(y = "Confirmed Cases", x = "Days", color="Country")
