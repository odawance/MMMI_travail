setwd("D:/Téléchargements/DAWANCE_Oceane_MMMI_travail")

par(mfrow = c(1,1))

## Daily report

daily_report <- read.table("Daily_number_of_new_reported_COVID-19_cases_and_deaths_worldwide.csv", header=TRUE, sep=",")
daily_report <- daily_report[daily_report$countriesAndTerritories == "Italy",]
daily_report <- daily_report[daily_report$year == "2020",]
daily_report <- as.data.frame(apply(daily_report, 2, rev))
View(daily_report)

barplot(height = as.numeric(daily_report$cases), names.arg = daily_report$dateRep, col = "green", 
        xlab = "Date", ylab = "Nombre de cas", main = "Cas de COVID-19 par jour en Italie en 2020")
barplot(height = as.numeric(daily_report$deaths), names.arg = daily_report$dateRep, col = "red", 
        xlab = "Date", ylab = "Nombre de deces", main = "Deces lies au COVID-19 par jour en Italie en 2020")

par(mfrow = c(1,2))

plot(cumsum(daily_report$cases), type = "l", xlab = "Nombre de jours depuis le 01 janvier 2020", ylab = "Nombre total de cas", main = "Nombre total de cas de COVID-19 en Italie en 2020", col = "green")
plot(cumsum(daily_report$deaths), type = "l", xlab = "Nombre de jours depuis le 01 janvier 2020", ylab = "Nombre total de deces", main = "Nombre total de deces lies au COVID-19 en Italie en 2020", col = "red")

## Cases with age

cases_with_age <- read.table("Data_on_the_14-day_age-specific_notification_rate_of_new_COVID-19_cases.csv", header=TRUE, sep=",")
cases_with_age <- cases_with_age[cases_with_age$ï..country == "Italy",]
cases_with_age <- cases_with_age[substr(cases_with_age$year_week,1,4) == "2020",]
View(cases_with_age)

age_groups <- unique(cases_with_age$age_group)

subset_data <- subset(cases_with_age, age_group == age_groups[1])
plot(seq_along(subset_data$year_week), subset_data$new_cases, type = "b", 
     xlab = "Semaine", ylab = "Nouveaux cas",ylim=c(0,max(cases_with_age$new_cases,na.rm = TRUE)), 
     main = "Nouveaux cas de COVID-19 par semaine en Italie en 2020", col=1)
for (i in 2:length(age_groups)){
  subset_data <- subset(cases_with_age, age_group == age_groups[i])
  lines(seq_along(subset_data$year_week), subset_data$new_cases, type = "b", col=i)
}
legend("topleft", legend = c(age_groups[1],age_groups[2],age_groups[3],age_groups[4],age_groups[5],age_groups[6]), 
       col = c(1,2,3,4,5,6), pch = c(1,1,1,1,1,1))

subset_data <- subset(cases_with_age, age_group == age_groups[1])
plot(seq_along(subset_data$year_week), subset_data$new_cases, type = "b", 
     xlab = "Semaine", ylab = "Nouveaux cas",ylim=c(0,1000), 
     main = "Nouveaux cas de COVID-19 par semaine en Italie en 2020", col=1)
for (i in 2:length(age_groups)){
  subset_data <- subset(cases_with_age, age_group == age_groups[i])
  lines(seq_along(subset_data$year_week), subset_data$new_cases, type = "b", col=i)
}
legend("topright", legend = c(age_groups[1],age_groups[2],age_groups[3],age_groups[4],age_groups[5],age_groups[6]), 
       col = c(1,2,3,4,5,6), pch = c(1,1,1,1,1,1))

## Hospitalizations

hospitalizations <- read.table("covid-hospitalizations.csv", header=TRUE, sep=",")
hospitalizations <- hospitalizations[hospitalizations$entity == "Italy",]
hospitalizations <- hospitalizations[substr(hospitalizations$date,1,4) == "2020",]
View(hospitalizations)

ICU_occupancy <- hospitalizations[hospitalizations$indicator == "Daily ICU occupancy",]
ICU_occupancy <- c(ICU_occupancy[,5])
plot(1:length(ICU_occupancy), ICU_occupancy, col="orange", xlab="Nombre de jours depuis le 24 fevrier 2020", 
     ylab="Nombre de personnes en soins intensifs", main="Occupation journaliere des USI en Italie en 2020", )

hospital_occupancy <- hospitalizations[hospitalizations$indicator == "Daily hospital occupancy",]
hospital_occupancy <- c(hospital_occupancy[,5])
plot(1:length(hospital_occupancy), hospital_occupancy, col="blue", xlab="Nombre de jours depuis le 24 fevrier 2020", 
     ylab="Nombre de personnes hospitalisees", main="Occupation journaliere des hopitaux en Italie en 2020", )

