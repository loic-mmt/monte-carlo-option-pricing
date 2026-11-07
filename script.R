
data <- read.csv("data/processed/selected_contracts.csv")
set.seed(123)

display <- data.frame(
  Echeance = data$expiry,
  Strike = data$K,
  Jours = round(365 * data$T, 2),
  `K/S0` = round(data$moneyness, 4),
  `IV (%)` = round(100 * data$sigma, 2),
  `Mid (points)` = round(data$midpoint, 2),
  check.names = FALSE
)

print(display, row.names = FALSE)


black_scholes <- function(spot_0, r, q, sigma, T, N = 1000) {
  Z <- rnorm(N)
  S_T <- spot_0 * exp((r - q - sigma^2 / 2) * T + sigma * sqrt(T) * Z)
  return(S_T)
}

calculate_Eq <- function(strike, spot) {
  mean(pmax(strike - spot, 0))
}

theoric_price <- function(r, T, Eq) {
  exp(-r * T) * Eq
}

analytic_price <- function(S_0, K, r, q, T, d1, d2) {
  V_BS <- K * exp(-r * T) * pnorm(-d2) - S_0 * exp(-q * T) * pnorm(-d1)
  return(V_BS)
}

calculate_d1 <- function(S_0, K, r, q, sigma, T) {
  denominator <- sigma * sqrt(T)
  numerator <- log(S_0 / K) + (r - q + (sigma^2 / 2)) * T
  d1 <- numerator / denominator
  return(d1)
}

calculate_d2 <- function(d1, sigma, T) {
  d2 <- d1 - sigma * sqrt(T)
  return(d2)
}

ic_percentile <- function(gains, B = 1000) {
  prix_bootstrap <- numeric(B)
  for (b in seq_len(B)) {
    moyennes <- numeric(length(gains))
    for (s in seq_along(gains)) {
      x <- gains[[s]]
      moyennes[s] <- mean(sample(x, size = length(x), replace = TRUE))
    }
    prix_bootstrap[b] <- mean(moyennes)
  }
  quantile(prix_bootstrap, c(0.025, 0.975))
}

nb_IC <- 100
B_bootstrap <- 1000


p <- data[1, ]

d1 <- calculate_d1(p$S0, p$K, p$r, p$q, p$sigma, p$T)
d2 <- calculate_d2(d1, p$sigma, p$T)
V_BT <- analytic_price(p$S0, p$K, p$r, p$q, p$T, d1, d2)


N <- c(1000, 5000, 20000, 100000)
list_V_0 <- numeric(1000)
RMSE <- numeric(length(N))
mean_V_0 <- numeric(length(N))
VAR <- numeric(length(N))
IC_mc <- matrix(NA_real_, nrow = nb_IC, ncol = 2)


for (n in seq_along(N)) {
  for (i in seq_len(1000)) {
    S_T <- black_scholes(p$S0, p$r, p$q, p$sigma, p$T, N[n])
    Eq <- calculate_Eq(p$K, S_T)
    V_0 <- theoric_price(p$r, p$T, Eq)

    list_V_0[i] <- V_0
    if (N[n] == 100000 && i <= nb_IC) {
      gains_mc <- exp(-p$r * p$T) * pmax(p$K - S_T, 0)
      IC_mc[i, ] <- ic_percentile(list(gains_mc), B = B_bootstrap)
    }
  }
  RMSE[n] <- sqrt(mean((list_V_0 - V_BT)^2))
  mean_V_0[n] <- mean(list_V_0)
  VAR[n] <- var(list_V_0)
}

options(scipen = 999)
chart <- data.frame(
  V_BT = V_BT,
  V_0 = mean_V_0,
  RMSE = RMSE,
  VAR = VAR,
  N = N,
  check.names = FALSE
)
prix_mc <- list_V_0


black_scholes_anti <- function(S0, K, r, q, sigma, T, N = 1000, retour_gains = FALSE) {
  stopifnot(N %% 2 == 0)
  Z <- rnorm(N / 2)
  Z_anti <- -Z
  S_T <- S0 * exp((r - q - sigma^2 / 2) * T + sigma * sqrt(T) * Z)
  S_T_anti <- S0 * exp((r - q - sigma^2 / 2) * T + sigma * sqrt(T) * Z_anti)

  gains <- exp(-r * T) * (pmax(K - S_T, 0) + pmax(K - S_T_anti, 0)) / 2
  if (retour_gains) return(gains)
  return(mean(gains))
}

N <- c(1000, 5000, 20000, 100000)
list_V_0 <- numeric(1000)
RMSE <- numeric(length(N))
mean_V_0 <- numeric(length(N))
VAR <- numeric(length(N))
IC_anti <- matrix(NA_real_, nrow = nb_IC, ncol = 2)

for (n in seq_along(N)) {
  for (i in seq_len(1000)) {
    gains_anti <- black_scholes_anti(
      p$S0, p$K, p$r, p$q, p$sigma, p$T, N[n], retour_gains = TRUE
    )
    list_V_0[i] <- mean(gains_anti)
    if (N[n] == 100000 && i <= nb_IC) {
      IC_anti[i, ] <- ic_percentile(list(gains_anti), B = B_bootstrap)
    }
  }
  RMSE[n] <- sqrt(mean((list_V_0 - V_BT)^2))
  mean_V_0[n] <- mean(list_V_0)
  VAR[n] <- var(list_V_0)
}

options(scipen = 999)
chart2 <- data.frame(
  V_BT = V_BT,
  V_0 = mean_V_0,
  RMSE = RMSE,
  VAR = VAR,
  N = N,
  check.names = FALSE
)
prix_anti <- list_V_0


black_scholes_strat <- function(S0, K, r, q, sigma, T, N = 1000, n_strata = 5, retour_gains = FALSE) {
  stopifnot(N %% n_strata == 0)

  gains_strates <- vector("list", n_strata)

  for (s in seq_len(n_strata)) {
    U <- runif(N / n_strata, min = (s - 1) / n_strata, max = s / n_strata)

    Z <- qnorm(U)
    S_T <- S0 * exp((r - q - sigma^2 / 2) * T + sigma * sqrt(T) * Z)

    gains_strates[[s]] <- exp(-r * T) * pmax(K - S_T, 0)
  }

  if (retour_gains) return(gains_strates)
  mean(sapply(gains_strates, mean))
}

N <- c(1000, 5000, 20000, 100000)
list_V_0 <- numeric(1000)
RMSE <- numeric(length(N))
mean_V_0 <- numeric(length(N))
VAR <- numeric(length(N))
IC_strat <- matrix(NA_real_, nrow = nb_IC, ncol = 2)

for (n in seq_along(N)) {
  for (i in seq_len(1000)) {
    gains_strat <- black_scholes_strat(
      p$S0, p$K, p$r, p$q, p$sigma, p$T, N[n], retour_gains = TRUE
    )
    list_V_0[i] <- mean(sapply(gains_strat, mean))
    if (N[n] == 100000 && i <= nb_IC) {
      IC_strat[i, ] <- ic_percentile(gains_strat, B = B_bootstrap)
    }
  }
  RMSE[n] <- sqrt(mean((list_V_0 - V_BT)^2))
  mean_V_0[n] <- mean(list_V_0)
  VAR[n] <- var(list_V_0)
}

options(scipen = 999)
chart3 <- data.frame(
  V_BT = V_BT,
  V_0 = mean_V_0,
  RMSE = RMSE,
  VAR = VAR,
  N = N,
  check.names = FALSE
)
prix_strat <- list_V_0


comparaison <- rbind(
  cbind(Methode = "Monte-Carlo", chart),
  cbind(Methode = "Antithétiques", chart2),
  cbind(Methode = "Stratification", chart3)
)

print(comparaison, row.names = FALSE)


anciens_parametres <- par(c("mfrow", "mar", "mgp", "cex"))
par(mfrow = c(1, 1), mar = c(3, 3.5, 0.8, 0.5),
    mgp = c(1.9, 0.6, 0), cex = 0.9)

reference_RMSE <- chart$RMSE[1] * sqrt(chart$N[1] / chart$N)
plot(chart$N, chart$RMSE, type = "b", log = "xy", xaxt = "n",
     xlab = "Nombre d'évaluations du payoff (N)", ylab = "RMSE (points)",
     ylim = range(chart$RMSE, chart2$RMSE, chart3$RMSE, reference_RMSE),
     col = "steelblue4", pch = 16, lwd = 2)
axis(1, at = chart$N, labels = format(chart$N, scientific = FALSE, trim = TRUE))
lines(chart2$N, chart2$RMSE, type = "b", col = "darkorange2", pch = 17, lwd = 2)
lines(chart3$N, chart3$RMSE, type = "b", col = "darkgreen", pch = 15, lwd = 2)
lines(chart$N, reference_RMSE, col = "grey40", lty = 2, lwd = 1.5)
legend("topright",
       legend = expression("Monte-Carlo", "Antithétiques", "Stratification", N^{-1/2}),
       col = c("steelblue4", "darkorange2", "darkgreen", "grey40"),
       pch = c(16, 17, 15, NA), lty = c(1, 1, 1, 2),
       lwd = c(2, 2, 2, 1.5), bty = "n", cex = 0.85)
par(anciens_parametres)


prix_methodes <- list(
  "Monte-Carlo" = prix_mc,
  "Antithétiques" = prix_anti,
  "Stratification" = prix_strat
)
IC_par_methode <- list(
  "Monte-Carlo" = IC_mc,
  "Antithétiques" = IC_anti,
  "Stratification" = IC_strat
)
bilan_IC <- data.frame(
  Methode = names(IC_par_methode),
  Couverture_pct = sapply(IC_par_methode, function(IC) {
    100 * mean(IC[, 1] <= V_BT & V_BT <= IC[, 2])
  }),
  Largeur_moyenne = sapply(IC_par_methode, function(IC) {
    mean(IC[, 2] - IC[, 1])
  }),
  row.names = NULL
)
estimation_IC <- sapply(prix_methodes, function(prix) prix[1])
IC_bootstrap <- lapply(IC_par_methode, function(IC) IC[1, ])


limites_x <- range(unlist(prix_methodes), V_BT, estimation_IC, unlist(IC_bootstrap))
coupures <- seq(limites_x[1], limites_x[2], length.out = 41)
histogrammes <- lapply(prix_methodes, hist, breaks = coupures, plot = FALSE)

anciens_parametres <- par(c("mfrow", "mar", "mgp", "cex", "xaxs", "yaxs"))
layout(rbind(1:3, c(4, 4, 4)), heights = c(1, 0.12))

par(mar = c(2.8, 2.8, 2, 0.35), mgp = c(1.65, 0.45, 0), cex = 0.95, xaxs = "i", yaxs = "i")
for (nom in names(prix_methodes)) {
  limite_y <- max(histogrammes[[nom]]$density) * 1.08
  plot(histogrammes[[nom]], freq = FALSE,
       main = nom, xlab = "Prix estimé (points)", ylab = "Densité",
       xlim = limites_x, ylim = c(0, limite_y), col = "grey85", border = "white")
  abline(v = V_BT, col = "red", lwd = 2)
  abline(v = estimation_IC[nom], col = "blue", lwd = 1.5)
  abline(v = IC_bootstrap[[nom]], col = "blue", lty = 2, lwd = 1.5)
  box(bty = "l")
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("V_BT", "Estimation retenue", "IC bootstrap (95 %)"),
       col = c("red", "blue", "blue"), lty = c(1, 1, 2), lwd = c(2, 1.5, 1.5),
       horiz = TRUE, bty = "n", cex = 0.9)
par(anciens_parametres)


print(bilan_IC, row.names = FALSE)


N_comparaison <- 100000
nb_repetitions <- 1000
chart4 <- data.frame()

for (s in seq_len(nrow(data))) {
  contrat <- data[s, ]

  d1 <- calculate_d1(contrat$S0, contrat$K, contrat$r, contrat$q, contrat$sigma, contrat$T)
  d2 <- calculate_d2(d1, contrat$sigma, contrat$T)
  V_BT_contrat <- analytic_price(contrat$S0, contrat$K, contrat$r, contrat$q, contrat$T, d1, d2)

  prix_contrats <- matrix(
    NA_real_, nrow = nb_repetitions, ncol = 3,
    dimnames = list(NULL, c("Monte-Carlo", "Antithétiques", "Stratification"))
  )

  for (i in seq_len(nb_repetitions)) {
    S_T <- black_scholes(
      contrat$S0, contrat$r, contrat$q, contrat$sigma, contrat$T,
      N = N_comparaison
    )
    prix_contrats[i, "Monte-Carlo"] <- theoric_price(
      contrat$r, contrat$T, calculate_Eq(contrat$K, S_T)
    )
    prix_contrats[i, "Antithétiques"] <- black_scholes_anti(
      contrat$S0, contrat$K, contrat$r, contrat$q, contrat$sigma, contrat$T,
      N = N_comparaison
    )
    prix_contrats[i, "Stratification"] <- black_scholes_strat(
      contrat$S0, contrat$K, contrat$r, contrat$q, contrat$sigma, contrat$T,
      N = N_comparaison
    )
  }

  for (methode in colnames(prix_contrats)) {
    prix <- prix_contrats[, methode]
    rmse <- sqrt(mean((prix - V_BT_contrat)^2))
    chart4 <- rbind(chart4, data.frame(
      Echeance = contrat$expiry,
      Strike = contrat$K,
      Methode = methode,
      V_BT = V_BT_contrat,
      V_0 = mean(prix),
      RMSE = rmse,
      RMSE_relatif_pct = 100 * rmse / V_BT_contrat,
      VAR = var(prix)
    ))
  }
}

print(chart4, row.names = FALSE)
