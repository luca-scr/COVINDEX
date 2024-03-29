# EPICENTRO: Dati della Sorveglianza integrata COVID-19 in Italia
# https://www.epicentro.iss.it/coronavirus/sars-cov-2-sorveglianza-dati
# https://www.epicentro.iss.it/coronavirus/open-data/calcolo_rt_italia.zip

library(EpiEstim)

EpiEstim_Rt_estimate <- function(date_of_analysis = Sys.Date(), file)
{

  ## parametri dell'intervallo seriale stimati da dati di contact tracing lombardi
  shape.stimato <- 1.87
  rate.stimato <- 0.28
  
  ## massimo numero di giorni dell'intervallo seriale
  N <- 300
  
  ## definisco la distribuzione dell'intervallo seriale
  intervallo.seriale <- dgamma(0:N, shape=shape.stimato, rate=rate.stimato) 
  
  ## normalizzo la distribuzione dell'intervallo seriale in modo che la somma faccia 1
  SI <- (intervallo.seriale/sum(intervallo.seriale)) 
  
  ## leggo la curva epidemica da un file con 3 colonne separate da spazi: data, numero di casi trasmessi localmente, numero di casi importati
  
  # LS:begin 
  if(missing(file))
  {
    # scarico file zip e leggo file csv in automatico
    url = "https://www.epicentro.iss.it/coronavirus/open-data/calcolo_rt_italia.zip"
    file = tempfile()
    download.file(url, file)
  }
  # oppure leggo file dati
  files = unzip(file, list=TRUE, junkpaths = TRUE)$Name
  # files = files[-agrep("__MACOSX", files)]
  f = files[agrep("curva_epidemica_Italia", files)]
  curva.epidemica <- read.table(unz(file, f))
  curva.epidemica[,1] <- as.Date(curva.epidemica[,1])
  names(curva.epidemica) <- c("dates", "local", "imported") ## assegno i nomi richiesti dal pacchetto EpiEstim
  date_of_analysis <- as.Date(date_of_analysis)
  curva.epidemica <- subset(curva.epidemica, dates <= date_of_analysis)
  # LS:end
  
## calcolo la stima di R applicando la funzione estimate_R del pacchetto EpiEstim
  stima <- estimate_R(incid = curva.epidemica, 
                      method = "non_parametric_si", 
                      config = make_config(list(si_distr = SI, 
                                                n1 = 10000, 
                                                mcmc_control = make_mcmc_control(thin=1, burnin=1000000))))

  ###################
  ### Attenzione! ###
  ###################
  ## La stima e' calcolata su tutta la curva epidemica specificata, ma il
  ##  pacchetto non puo' tenere conto dei ritardi di inserimento nel dato.
  ## Le stime di Rt varieranno man mano che vengono inseriti nuovi casi con 
  ## data di inizio sintomi indietro nel tempo.
  ## Per questo motivo ISS considera valide le stime fino a 14 giorni prima 
  ## della data in cui viene effettuata la stima.
  ## Questo ritardo puo' cambiare nel tempo
  
  ## estraggo i risultati di interesse
  R.medio <- stima$R$`Mean(R)` ## valore medio
  R.lowerCI <- stima$R$`Quantile.0.025(R)` ## estremo inferiore dell'intervallo di confidenza
  R.upperCI <- stima$R$`Quantile.0.975(R)` ## estremo superiore dell'intervallo di confidenza

  ## estraggo le date di riferimento per la stima di R
  ## la data rappresenta il giorno centrale intorno a cui e' calcolata la media
  ##  mobile di Rt (con finestra di ampiezza pari a una settimana)
  sel.date <- stima$R[, "t_end"]
  date <- curva.epidemica[sel.date,1] 

  out <- data.frame(date = as.Date(date), 
                    # incidence = curva.epidemica[curva.epidemica$dates %in% as.Date(date),"local"],
                    Rt = R.medio,
                    lower = R.lowerCI,
                    upper = R.upperCI)
  out <- subset(out, date <= (max(as.Date(date))-14))
  return(out)
}
