# R code to reproduce the analysis of the paper
# Scrucca L. (2021) A COVINDEX based on a GAM beta regression model 
#   with an application to the COVID-19 pandemic in Italy

source("covindex.R")
load_install_package("ggplot2")
load_install_package("ggthemes")
load_install_package("ggrepel")
load_install_package("gridExtra")
load_install_package("kableExtra")

# Data ------------------------------------------------------------

source("covindex_data.R")
head(COVID19)
head(COVID19_REGIONS)

# Italy -----------------------------------------------------------

# Starting date of analysis
COVID19 = COVID19[date >= "2020-03-01",]
# Ending date of analysis
COVID19 = COVID19[date <= "2021-06-30",]
COVID19$t = 1:nrow(COVID19)

ggplot(COVID19, aes(x = date, y = y, size = swabs)) +
  geom_point(alpha = 0.5) +
  scale_size(range = c(.01, 3), name="Swabs:  ",
             # labels =  scales::unit_format(unit = "K", scale = 1e-3)) +
             labels =  scales::label_number()) +
  labs(x = "", y = "Test positive rate") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels = monthyear_labels) +
  scale_y_continuous(lim = c(0,0.35),
                     labels = scales::label_percent(accuracy = 1), 
                     sec.axis = dup_axis(name = ""),
                     minor_breaks = function(y) minor_breaks(y,2)) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", vjust = 1, size = 10),
        legend.pos = "bottom")

# Fit the model
mod = covindex_gam_betareg(
  y = COVID19$y, 
  t = COVID19$t,
  x = data.frame(weekend = COVID19$weekend),
  wts = COVID19$swabs/mean(COVID19$swabs, na.rm=TRUE))
summary(mod)

ggplot(as.data.frame(mod[c("k", "df", "AIC")]), aes(x = k, y = AIC)) +
  geom_line() + 
  geom_point(col = "white", size = 3) + 
  geom_point() + 
  scale_x_continuous(name = "k", breaks = function(lim) pretty(lim, n=9),
                     sec.axis = sec_axis(trans = approxfun(crit$bs, crit$df), 
                                         breaks = function(lim) pretty(lim, n=7),
                                         name = "EDF")) +
  geom_vline(xintercept = mod$k_opt, lty = 2) +
  theme_bw()

plot(mod, residuals = TRUE, shade = TRUE, rug = FALSE)
qq.gam(mod, rep = 100, level = 0.95, type = "deviance", pch = 1)

par(mfrow = c(1,2))
acf(residuals(mod, type = "deviance"), lag.max = 21)
pacf(residuals(mod, type = "deviance"), lag.max = 21, ylim = c(-0.5, 0.5))

# Predictions
pred = simpred_covindex_gam_betareg(
  mod,
  newdata = data.frame(t = COVID19$t, 
                       weekend = FALSE) )
CI_y = cbind(COVID19, 
             list2DF(pred[c("mu", 
                            "lower_credint", "upper_credint", 
                            "lower_predint", "upper_predint")]))

ggplot(CI_y, aes(x = date, y = mu)) +
  # model fit
  geom_ribbon(aes(ymin = lower_credint, ymax = upper_credint),
              alpha = 1, fill = "grey50") +
  geom_ribbon(aes(ymin = lower_predint, ymax = upper_predint),
              alpha = 0.5, fill = "grey50") +
  geom_line(col = "dodgerblue2", lwd = 1) +
  # data
  geom_point(aes(y = y), size = 0.5) +
  #
  labs(y= "Test positive rate", x = "") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels = monthyear_labels) +
  scale_y_continuous(lim = c(0,0.34), 
                     breaks = seq(0,0.3,by=0.05),
                     expand = expansion(mult = 0, add = 0),
                     labels = scales::label_percent(accuracy = 1), 
                     sec.axis = dup_axis(name = ""),
                     minor_breaks = function(y) minor_breaks(y,2)) +
  theme_minimal()

covindex_sim = apply(pred$mu_sim, 2, covindex)
credint = apply(covindex_sim, 1, quantile, 
                prob = c(0.025, 0.975), na.rm = TRUE)
CI_covindex = cbind(COVID19, 
                    covindex = covindex(pred$mu),
                    lower_credint = credint[1,], 
                    upper_credint = credint[2,])

ggplot(CI_covindex, aes(x = date, y = covindex)) +
  geom_hline(yintercept = 1, lty = 2) +
  # model fit
  geom_ribbon(aes(ymin = lower_credint, ymax = upper_credint),
              alpha = 1, fill = "grey") +
  geom_line(col = "dodgerblue2", lwd = 1) +
  #
  labs(y = "COVINDEX", x = "") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels =  monthyear_labels) +
  scale_y_continuous(limits = c(0.5,2), expand=c(0,0),
                     trans = "log10", 
                     breaks = seq(0.5, 2, by=0.1),
                     sec.axis = dup_axis(name=""),
                     minor_breaks = function(y) minor_breaks(y,1)) +
  theme_minimal()

COVINDEX_lag = data.frame("date" = COVID19$date, 
                          "covindex1"  = covindex(pred$mu, lag = 1),
                          "covindex3"  = covindex(pred$mu, lag = 3),
                          "covindex5"  = covindex(pred$mu, lag = 5),
                          "covindex7"  = covindex(pred$mu, lag = 7),
                          "covindex10" = covindex(pred$mu, lag = 10),
                          "covindex14" = covindex(pred$mu, lag = 14),
                          check.names = FALSE)  

cols = c("black", "#D7191C", "#FDAE61", "dodgerblue2", "#A6D96A", "#1A9641")
matplot(COVINDEX_lag$date, COVINDEX_lag[,-1], type = "l",
        lty = 1, col = cols, log = "y", xaxt = "n",
        xlab = "", ylab = "COVINDEX")
axis.Date(side = 1, at = seq(min(COVINDEX_lag$date), max(COVINDEX_lag$date), by="1 mon"), format="%m-%y")

		
fit = COVID19[, c("date", "t", "y")]
fit[, mu := pred$mu]
fit[, covindex := covindex(fit$mu)]
xlim <- range(0, max(fit$mu,na.rm=TRUE)*1.05)
ylim <- range(0.5, 2, range(fit$covindex,na.rm=TRUE)*c(0.95,1.05))
df = copy(fit)
df = df[mday(date) == 1  | date == "2020-03-08",] 
df[, label := format(date, format = "%b %y")]
df[, nudge_x := c(NA, 0.01, 0.02, 0, 0.03, 0.02, 0, 0.02, 0.02, 
                  0.01, 0, -0.01, -0.01, -0.01, 0.02, -0.01, 0)]
df[, nudge_y := c(NA, 0.05, -0.05, -0.05, -0.01, 0, 0.1, 0, 0, 
                  0.05, -0.05, 0.05, 0.05, 0.05, 0, 0.02, 0.03)]
quadcols = c("red3", "darkorange", "gold2")
quadcols = sapply(quadcols, function(...) adjustcolor(..., alpha = 0.4))

ggplot(fit, aes(x = mu, y = covindex)) +
  coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) + 
  annotate("rect", fill = quadcols[3],
           xmin = xlim[1], xmax = 0.05, 
           ymin = ylim[1], ymax = 1) +
  annotate("rect", fill = quadcols[2],
           xmin = 0.05, xmax = xlim[2], 
           ymin = ylim[1], ymax = 1) +
  annotate("rect", fill = quadcols[2],
           xmin = xlim[1], xmax = 0.05, 
           ymin = 1, ymax = ylim[2]) +
  annotate("rect", fill = quadcols[1],
           xmin = 0.05, xmax = xlim[2], 
           ymin = 1, ymax = ylim[2]) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_vline(xintercept = 0.05, lty = 2) +
  geom_path(col = "grey30", size = 1, lineend = "round") +
  geom_point(aes(x = tail(mu,1), y = tail(covindex,1)), size = 2) +
  geom_text_repel(data = df, aes(x = mu, y = covindex, label = label), 
                  nudge_x = df$nudge_x,
                  nudge_y = df$nudge_y) +
  scale_x_continuous(breaks = seq(0, max(xlim), by = 0.05),
                     labels = scales::label_percent(accuracy = 1), 
                     minor_breaks = function(x) minor_breaks(x,2)) +
  scale_y_continuous(trans = "log10", 
                     breaks = seq(floor(min(ylim)*10)/10,
                                  ceiling(max(ylim)*10)/10, 
                                  by=0.1)) +
  labs(x = "Test Positive Rate", y = "COVINDEX") + 
  theme_minimal()


# Empirical coverage ----------------------------------------------

# Starting date of analysis 2020-09-01
DT = copy(COVID19[date >= as.Date("2020-09-01")])
# Forecast from 2021-02-01 to 2021-05-31
forecast = list(start = as.Date("2021-02-01"), 
                end = as.Date("2021-05-31"))
forecast$ndays = as.numeric(forecast$end - forecast$start + 1)
# for time horizon of 
forecast$horizon = 1:14

COVERAGE = matrix(NA, nrow = forecast$ndays, ncol = length(forecast$horizon))
for(i in 1:forecast$ndays)
{
  # fit model
  DT_train = DT[date < (forecast$start+i-1),]
  mod = covindex_gam_betareg(
    y = DT_train$y, 
    t = DT_train$t, 
    x = data.frame(weekend = DT_train$weekend),
    wts = DT_train$swabs/mean(DT_train$swabs, na.rm=TRUE))
  # predict at horizon
  DT_test = COVID19[date %in% (forecast$start + (i-1) + (forecast$horizon-1))]
  pred = simpred_covindex_gam_betareg(
    mod, 
    newdata = data.frame(t = DT_test$t, 
                         weekend = DT_test$weekend))
  DT_test = cbind(DT_test, list2DF(pred[c("mu", 
                                          "lower_credint", "upper_credint", 
                                          "lower_predint", "upper_predint")]))
  COVERAGE[i,] = (DT_test[, y >= lower_predint & y <= upper_predint])
}

colnames(COVERAGE) = forecast$horizon
colMeans(COVERAGE)

dt = data.table(
  date = rep(seq(forecast$start, forecast$end, by = 1)-1, times = ncol(COVERAGE)),
  horizon = rep(as.factor(forecast$horizon), each = nrow(COVERAGE)), 
  value = c(COVERAGE))

plot1 = 
  ggplot(dt, aes(x = date, y = horizon, fill = value)) + 
  geom_tile(colour = "darkgrey") +
  scale_fill_manual(values = c(adjustcolor("red", alpha.f = 0.5), "white")) +
  labs(x = "", y = "Prediction time horizon (days)") +
  scale_x_date(date_breaks = "1 month",
               expand = c(0.01,0.01),
               date_labels = monthyear_labels) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot2 = 
  ggplot(data.frame(coverage = colMeans(COVERAGE), horizon = forecast$horizon),
         aes(x = coverage, y = horizon)) +
  geom_point(pch = 15) +
  geom_vline(xintercept = 0.95, lty = 2) +
  scale_x_continuous(name = "Coverage", 
                     lim = c(0.5,1), 
                     labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(breaks = forecast$horizon, 
                     name = NULL,
                     position = "right") +
  theme_minimal()

grid.arrange(plot1, plot2, nrow = 1, ncol = 2, widths = c(3,1))


# Italian regions -------------------------------------------------

REGIONS = names(COVID19_REGIONS)
START = "2020-09-01"
DATES = c("2020-11-01", "2020-12-01", "2021-02-01", "2021-03-01")
out = vector(mode = "list", length = length(REGIONS))
names(out) = REGIONS
for(i in 1:length(REGIONS))
{
  DT = COVID19_REGIONS[[i]]
  DT = DT[date >= START,]
  DT$t = 1:nrow(DT)
  mod = covindex_gam_betareg(
    y = DT$y, 
    t = DT$t, 
    x = data.frame(weekend = DT$weekend),
    wts = DT$swabs/mean(DT$swabs, na.rm=TRUE) )
  pred = simpred_covindex_gam_betareg(
    mod, 
    newdata = data.frame(t = DT$t, 
                         weekend = FALSE) )
  COVINDEX = cbind(DT, mu = pred$mu, covindex = covindex(pred$mu))
  sel = sapply(DATES, function(d) which(COVINDEX$date == d))
  out[[i]] = COVINDEX[sel, c("date", "mu", "covindex")]
}

quadcols = c("red3", "darkorange", "gold2")
quadcols = sapply(quadcols, function(...) adjustcolor(..., alpha = 0.4))
xlim = range(0, 0.25, range(sapply(out, function(x) range(x$mu,na.rm=TRUE))))
ylim = range(0.5, 2, range(sapply(out, function(x) range(x$covindex,na.rm=TRUE))))
plot = vector(mode = "list", length = length(DATES))
for(i in seq(DATES))
{
  DT = lapply(out, function(x) x[date == DATES[i],])
  DT = data.table(region = REGIONS, 
                  mu = sapply(DT, function(x) x$mu),
                  covindex = sapply(DT, function(x) x$covindex))
  plot[[i]] = 
    ggplot(DT, aes(x = mu, y = covindex, label = region)) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) + 
    annotate("rect", fill = quadcols[3],
             xmin = xlim[1], xmax = 0.05, 
             ymin = ylim[1], ymax = 1) +
    annotate("rect", fill = quadcols[2],
             xmin = 0.05, xmax = xlim[2], 
             ymin = ylim[1], ymax = 1) +
    annotate("rect", fill = quadcols[2],
             xmin = xlim[1], xmax = 0.05, 
             ymin = 1, ymax = ylim[2]) +
    annotate("rect", fill = quadcols[1],
             xmin = 0.05, xmax = xlim[2], 
             ymin = 1, ymax = ylim[2]) +
    geom_point() +
    geom_text_repel() +
    geom_hline(yintercept = 1, lty = 2) +
    geom_vline(xintercept = 0.05, lty = 2) +
    scale_y_continuous(trans = "log10", 
                       breaks = c(0.5, 0.75, 1, 1.5, 2)) +
    ggtitle("", subtitle = format(as.Date(DATES[i]), "%d-%m-%Y")) +
    labs(x = "Test Positive Rate", y = "COVINDEX") + 
    theme_minimal()
}  
grid.arrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
             nrow = 2, ncol = 2)


# Comparison of COVINDEX with Rt ----------------------------------

## Italy ----

source("EpiEstim_Rt_estimate.R")
date_of_analysis = as.Date("2021-03-28")
date_of_analysis = as.Date("2021-06-30")
Rt = EpiEstim_Rt_estimate(date_of_analysis, 
                          file = "data/calcolo_rt_italia.zip")

DT = COVID19[date <= date_of_analysis,]
mod = covindex_gam_betareg(
  y = DT$y, 
  t = DT$t, 
  x = data.frame(weekend = DT$weekend),
  wts = DT$swabs/mean(DT$swabs, na.rm=TRUE) )
pred = simpred_covindex_gam_betareg(
  mod, 
  newdata = data.frame(t = DT$t,
                       weekend = FALSE) )
covindex_sim = apply(pred$mu_sim, 2, covindex)
credint = apply(covindex_sim, 1, quantile, 
                prob = c(0.025, 0.975), na.rm = TRUE)
covindex_pred = data.table(date = DT$date, 
                           covindex = covindex(pred$mu),
                           lower_credint = credint[1,], 
                           upper_credint = credint[2,])

plot1 = 
  ggplot(data = Rt, aes(x = date, y = Rt)) +
  geom_hline(yintercept = 1, lty = 2) +
  # covindex
  geom_ribbon(data = covindex_pred,
              aes(x = date, ymin = lower_credint, ymax = upper_credint),
              alpha = 0.8, fill = "grey", inherit.aes = FALSE) +
  geom_line(data = covindex_pred,
            aes(x = date, y = covindex, color = "COVINDEX"), size = 1) +
  # Rt
  geom_line(aes(color = "Rt"), lty = 3, size = 1) +
   geom_line(data = Rt, aes(x = date, y = Rt, color = "Rt"), size = 1) +
  #
  labs(title = "", x = "", y = "", color = "Legend") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels = monthyear_labels,
               limits = c(as.Date("2020-03-01"), NA)) +
  scale_y_continuous(limits = c(0.5,2), expand=c(0,0),
                     trans = "log10", 
                     breaks = seq(0.5, 2, by=0.1),
                     sec.axis = dup_axis(),
                     minor_breaks = function(y) minor_breaks(y,1)) +
  scale_color_manual(name = "",
                     values = c("COVINDEX" = "dodgerblue2", 
                                "Rt" = "red3"),
                     labels = c("COVINDEX", "Rt")) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.margin = margin(c(0, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 20, unit = "pt")))


date_of_analysis = as.Date("2020-12-05")
Rt = EpiEstim_Rt_estimate(date_of_analysis, 
                          file = "data/calcolo_rt_italia.zip")
Rt_post = EpiEstim_Rt_estimate(date_of_analysis+7, 
                          file = "data/calcolo_rt_italia.zip")

DT = COVID19[date <= date_of_analysis,]
mod = covindex_gam_betareg(y = DT$y, 
                           t = DT$t, 
                           x = data.frame(weekend = DT$weekend),
                           wts = DT$swabs/mean(DT$swabs, na.rm=TRUE))
pred = simpred_covindex_gam_betareg(mod, 
                                    newdata = data.frame(t = DT$t,
                                                         weekend = FALSE))
covindex_sim = apply(pred$mu_sim, 2, covindex)
credint = apply(covindex_sim, 1, quantile, 
                prob = c(0.025, 0.975), na.rm = TRUE)
covindex_pred = data.table(date = DT$date, 
                           covindex = covindex(pred$mu),
                           lower_credint = credint[1,], 
                           upper_credint = credint[2,])

plot2 = ggplot(data = Rt_post, aes(x = date, y = Rt)) +
  geom_hline(yintercept = 1, lty = 2) +
  # covindex
  geom_ribbon(data = covindex_pred,
              aes(x = date, ymin = lower_credint, ymax = upper_credint),
              alpha = 0.8, fill = "grey", inherit.aes = FALSE) +
  geom_line(data = covindex_pred,
            aes(x = date, y = covindex, color = "COVINDEX"), size = 1) +
  # Rt
  geom_line(aes(color = "Rt"), lty = 3, size = 1) +
  geom_line(data = Rt, aes(x = date, y = Rt, color = "Rt"), size = 1) +
  #
  labs(title = "", x = "", y = "", color = "Legend") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels = monthyear_labels,
               limits = c(as.Date("2020-09-01"), NA)) +
  scale_y_continuous(limits = c(0.5,2), expand=c(0,0),
                     trans = "log10", 
                     breaks = seq(0.5, 2, by=0.1),
                     # sec.axis = dup_axis(),
                     minor_breaks = function(y) minor_breaks(y,1)) +
  scale_color_manual(name = "",
                     values = c("COVINDEX" = "dodgerblue2", 
                                "Rt" = "red3"),
                     labels = c("COVINDEX", "Rt")) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(legend.position = "none",
        legend.margin = margin(c(0, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 20, unit = "pt")))

date_of_analysis = as.Date("2021-02-28")
Rt = EpiEstim_Rt_estimate(date_of_analysis, 
                          file = "data/calcolo_rt_italia.zip")
Rt_post = EpiEstim_Rt_estimate(date_of_analysis+7, 
                          file = "data/calcolo_rt_italia.zip")

DT = COVID19[date <= date_of_analysis,]
mod = covindex_gam_betareg(y = DT$y, 
                           t = DT$t, 
                           x = data.frame(weekend = DT$weekend),
                           wts = DT$swabs/mean(DT$swabs, na.rm=TRUE))
pred = simpred_covindex_gam_betareg(mod, 
                                    newdata = data.frame(t = DT$t,
                                                         weekend = FALSE))
covindex_sim = apply(pred$mu_sim, 2, covindex)
credint = apply(covindex_sim, 1, quantile, 
                prob = c(0.025, 0.975), na.rm = TRUE)
covindex_pred = data.table(date = DT$date, 
                           covindex = covindex(pred$mu),
                           lower_credint = credint[1,], 
                           upper_credint = credint[2,])

plot3 = ggplot(data = Rt_post, aes(x = date, y = Rt)) +
  geom_hline(yintercept = 1, lty = 2) +
  # covindex
  geom_ribbon(data = covindex_pred,
              aes(x = date, ymin = lower_credint, ymax = upper_credint),
              alpha = 0.8, fill = "grey", inherit.aes = FALSE) +
  geom_line(data = covindex_pred,
            aes(x = date, y = covindex, color = "COVINDEX"), size = 1) +
  # Rt
  geom_line(aes(color = "Rt"), lty = 3, size = 1) +
  geom_line(data = Rt, aes(x = date, y = Rt, color = "Rt"), size = 1) +
  #
  labs(title = "", x = "", y = "", color = "Legend") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels = monthyear_labels,
               limits = c(as.Date("2021-01-01"), NA)) +
  scale_y_continuous(limits = c(0.5,2), expand=c(0,0),
                     trans = "log10", 
                     breaks = seq(0.5, 2, by=0.1),
                     position = "right",
                     # sec.axis = dup_axis(),
                     minor_breaks = function(y) minor_breaks(y,1)) +
  scale_color_manual(name = "",
                     values = c("COVINDEX" = "dodgerblue2", 
                                "Rt" = "red2"),
                     labels = c("COVINDEX", "Rt")) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(legend.position = "none",
        legend.margin = margin(c(0, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 20, unit = "pt")))
				
grid.arrange(plot1, plot2, plot3,
             layout_matrix = matrix(c(1,2,1,3), nrow = 2, ncol = 2))


## Lombardia ----

START = "2020-09-01"
END   = "2021-06-30"
DT = COVID19_REGIONS$Lombardia
DT = DT[date >= START,]
DT = DT[date <= END,]
DT$t = 1:nrow(DT)
DT

mod = covindex_gam_betareg(
  y = DT$y, 
  t = DT$t, 
  x = data.frame(weekend = DT$weekend),
  wts = DT$swabs/mean(DT$swabs, na.rm=TRUE))
summary(mod)

pred = simpred_covindex_gam_betareg(
  mod,
  newdata = data.frame(t = DT$t, 
                       weekend = FALSE) )
CI_y = cbind(DT, list2DF(pred[c("mu", 
                                "lower_credint", "upper_credint", 
                                "lower_predint", "upper_predint")]))

ggplot(CI_y, aes(x = date, y = mu)) +
  # model fit
  geom_ribbon(aes(ymin = lower_credint, ymax = upper_credint),
              alpha = 1, fill = "grey50") +
  geom_ribbon(aes(ymin = lower_predint, ymax = upper_predint),
              alpha = 0.5, fill = "grey50") +
  geom_line(col = "dodgerblue2", lwd = 1) +
  # data
  geom_point(aes(y = y), size = 0.5) +
  #
  labs(y= "Test positive rate", x = "") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels = monthyear_labels) +
  scale_y_continuous(lim = c(0,0.34), 
                     breaks = seq(0,0.3,by=0.05),
                     expand = expansion(mult = 0, add = 0),
                     labels = scales::label_percent(accuracy = 1), 
                     sec.axis = dup_axis(name = ""),
                     minor_breaks = function(y) minor_breaks(y,2)) +
  theme_minimal()

covindex_sim = apply(pred$mu_sim, 2, covindex)
credint = apply(covindex_sim, 1, quantile, 
                prob = c(0.025, 0.975), na.rm = TRUE)
CI_covindex = cbind(DT, 
                    covindex = covindex(pred$mu),
                    lower_credint = credint[1,], 
                    upper_credint = credint[2,])

ggplot(CI_covindex, aes(x = date, y = covindex)) +
  geom_hline(yintercept = 1, lty = 2) +
  # model fit
  geom_ribbon(aes(ymin = lower_credint, ymax = upper_credint),
              alpha = 1, fill = "grey") +
  geom_line(aes(color = "COVINDEX"), size = 1) +
  geom_point(data = Rt_REGIONI[Regione == "Lombardia",],
             aes(x = date, y = Rt, color = "Rt")) +
  geom_line(data = Rt_REGIONI[Regione == "Lombardia",],
            aes(x = date, y = Rt, color = "Rt")) + 
  #
  scale_color_manual(name = "",
                     values = c("COVINDEX" = "dodgerblue2", 
                                "Rt" = "red3"),
                     labels = c("COVINDEX", "Rt")) +
  labs(y = "COVINDEX", x = "") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels =  monthyear_labels) +
  scale_y_continuous(trans = "log10", 
                     breaks = seq(0.5, 2.5, by=0.1),
                     sec.axis = dup_axis(name=""),
                     minor_breaks = function(y) minor_breaks(y,1)) +
  coord_cartesian(ylim = c(0.5, 2.5), expand = FALSE) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.margin = margin(c(0, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 20, unit = "pt")))

## Umbria ----

START = "2020-09-01"
END   = "2021-06-30"
DT = COVID19_REGIONS$Umbria
DT = DT[date >= START,]
DT = DT[date <= END,]
DT$t = 1:nrow(DT)
DT

mod = covindex_gam_betareg(
  y = DT$y, 
  t = DT$t, 
  x = data.frame(weekend = DT$weekend),
  wts = DT$swabs/mean(DT$swabs, na.rm=TRUE))
summary(mod)

pred = simpred_covindex_gam_betareg(
  mod,
  newdata = data.frame(t = DT$t, 
                       weekend = FALSE) )
CI_y = cbind(DT, list2DF(pred[c("mu", 
                                "lower_credint", "upper_credint", 
                                "lower_predint", "upper_predint")]))

ggplot(CI_y, aes(x = date, y = mu)) +
  # model fit
  geom_ribbon(aes(ymin = lower_credint, ymax = upper_credint),
              alpha = 1, fill = "grey50") +
  geom_ribbon(aes(ymin = lower_predint, ymax = upper_predint),
              alpha = 0.5, fill = "grey50") +
  geom_line(col = "dodgerblue2", lwd = 1) +
  # data
  geom_point(aes(y = y), size = 0.5) +
  #
  labs(y= "Test positive rate", x = "") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels = monthyear_labels) +
  scale_y_continuous(lim = c(0,0.34), 
                     breaks = seq(0,0.3,by=0.05),
                     expand = expansion(mult = 0, add = 0),
                     labels = scales::label_percent(accuracy = 1), 
                     sec.axis = dup_axis(name = ""),
                     minor_breaks = function(y) minor_breaks(y,2)) +
  theme_minimal()

covindex_sim = apply(pred$mu_sim, 2, covindex)
credint = apply(covindex_sim, 1, quantile, 
                prob = c(0.025, 0.975), na.rm = TRUE)
CI_covindex = cbind(DT, 
                    covindex = covindex(pred$mu),
                    lower_credint = credint[1,], 
                    upper_credint = credint[2,])

ggplot(CI_covindex, aes(x = date, y = covindex)) +
  geom_hline(yintercept = 1, lty = 2) +
  # model fit
  geom_ribbon(aes(ymin = lower_credint, ymax = upper_credint),
              alpha = 1, fill = "grey") +
  geom_line(aes(color = "COVINDEX"), size = 1) +
  geom_point(data = Rt_REGIONI[Regione == "Umbria",],
             aes(x = date, y = Rt, color = "Rt")) +
  geom_line(data = Rt_REGIONI[Regione == "Umbria",],
            aes(x = date, y = Rt, color = "Rt")) + 
  #
  scale_color_manual(name = "",
                     values = c("COVINDEX" = "dodgerblue2", 
                                "Rt" = "red3"),
                     labels = c("COVINDEX", "Rt")) +
  labs(y = "COVINDEX", x = "") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels =  monthyear_labels) +
  scale_y_continuous(trans = "log10", 
                     breaks = seq(0.5, 2.5, by=0.1),
                     sec.axis = dup_axis(name=""),
                     minor_breaks = function(y) minor_breaks(y,1)) +
  coord_cartesian(ylim = c(0.5, 2), expand = FALSE) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.margin = margin(c(0, 0, 0, 0)),
        legend.text = element_text(margin = margin(r = 20, unit = "pt")))
