## Read and prepare for the analyses the data downloaded from 
## Dipartimento della Protezione Civile
## COVID-19 Italia - Monitoraggio della situazione 
## http://arcg.is/C1unv
## Source: https://github.com/pcm-dpc/COVID-19


# national data ---------------------------------------------------

COVID19 = fread("data/dpc-covid19-ita-andamento-nazionale.csv")
# consider only molecular swabs and corresponding positive cases
COVID19$date = as.Date(COVID19$data)
COVID19[date >= "2021-01-15", "tamponi"] = 
  COVID19[date >= "2021-01-15", "tamponi_test_molecolare"]
COVID19[date >= "2021-01-15", "totale_casi"] = 
  COVID19[date >= "2021-01-15", "totale_positivi_test_molecolare"]
#
COVID19[, positives := { x = c(NA, diff(totale_casi)); ifelse(x < 0, NA, x) } ]
COVID19[, swabs  := { x = c(NA, diff(tamponi)); ifelse(x < 0, NA, x) } ]
# correct for some errors in the recorded data
COVID19[date == "2020-12-17", "swabs"] <- 185320
COVID19[date == "2020-06-19", "positives"] <- 251
COVID19[date == "2021-03-22", "positives"] <- 13846
COVID19[date == "2021-03-22", "positives"] <- 13846
#
COVID19[, y := positives/swabs]
COVID19 = COVID19[,c("date", "positives", "swabs", "y")]
# set dummy variable for days after sundays and holidays to consider 
# limited testing capacity
Italian_holidays = as.Date(
  c("2020-01-01", "2020-01-06", "2020-04-12", "2020-04-13", "2020-04-25",
    "2020-05-01", "2020-06-02", "2020-08-15", "2020-11-01", "2020-12-08", 
    "2020-12-25", "2020-12-26",
    "2021-01-01", "2021-01-06", "2021-04-04", "2021-04-05", "2021-04-25", 
    "2021-05-01", "2021-06-02", "2021-08-15", "2021-11-01", "2021-12-08", 
    "2021-12-25", "2021-12-26"))
COVID19[, weekend := (weekdays(date) %in% c("Sunday", "Monday") |
                      date %in% (Italian_holidays+1))]

# regional data ---------------------------------------------------

COVID19_REGIONI = fread("data/dpc-covid19-ita-regioni.csv")
COVID19_REGIONI$date = as.Date(COVID19_REGIONI$data)
REGIONE = unique(COVID19_REGIONI$denominazione_regione)
COVID19_REGIONS = vector(mode = "list", length = length(REGIONE))
names(COVID19_REGIONS) = REGIONE
for(i in 1:length(REGIONE))
{
  DT = copy(COVID19_REGIONI)
  DT = DT[denominazione_regione == REGIONE[i],]
  # consider only molecular swabs and corresponding positive cases
  DT[date >= "2021-01-15", "tamponi"] = 
    DT[date >= "2021-01-15", "tamponi_test_molecolare"]
  DT[date >= "2021-01-15", "totale_casi"] = 
    DT[date >= "2021-01-15", "totale_positivi_test_molecolare"]
  #
  DT[, positives := { x = c(NA, diff(totale_casi)); ifelse(x < 0, NA, x) } ]
  DT[, swabs     := { x = c(NA, diff(tamponi)); ifelse(x < 0, NA, x) } ]
  # correct for some errors in the recorded data
  if(REGIONE[i] == "Umbria") 
    DT[date == "2021-02-08"]$swabs <- 754
  #
  DT[, y := positives/swabs]
  DT = DT[,c("date", "positives", "swabs", "y")]
  DT[, weekend := (weekdays(date) %in% c("Sunday", "Monday") |
                   date %in% (Italian_holidays+1))]
  COVID19_REGIONS[[i]] = DT
}
rm(COVID19_REGIONI)

Rt_REGIONI = fread("data/dpc-covid19-ita-regioni.csv")
Rt_REGIONI = fread("data/Rt_regions.csv")
Rt_REGIONI$date = as.Date(Rt_REGIONI$date)

