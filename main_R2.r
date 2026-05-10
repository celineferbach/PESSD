# ── Librairies ────────────────────────────────────────────────────────────────
library(plm)
library(dplyr)
library(systemfit)


# ── Chargement des données ────────────────────────────────────────────────────
# panelV2.csv est généré par main.ipynb (cellule d'export)
# Colonnes : Country, Year, dep_sante, pib_ph, part_65,
#            practiciens, lits, tx_deces_2ans, chomage, dep_sante_lag
panel <- read.csv("panelV2.csv", sep = ",")
head(panel)

# Renommer pour plus de clarté
colnames(panel) <- c("pays", "annee", "dep_sante", "pib_ph", "part_65",
                     "practiciens", "lits", "tx_deces_2ans", "chomage", "dep_sante_lag")

# Déclarer la structure panel
panel <- pdata.frame(panel, index = c("pays", "annee"))
panel <- subset(panel, select = -c(dep_sante_lag))
head(panel)


# ══════════════════════════════════════════════════════════════════════════════
# 1. Blundell-Bond (System GMM) — modèle de base
# ══════════════════════════════════════════════════════════════════════════════
# transformation = "ld" → level and difference = System GMM (Blundell-Bond)

bb_model <- pgmm(
  dep_sante ~ plm::lag(dep_sante, 1) + plm::lag(dep_sante, 2) + plm::lag(dep_sante, 3)
            + part_65        # exogène
            + pib_ph         # endogène
            + chomage        # exogène
            + tx_deces_2ans  # endogène

            | lag(dep_sante, 4:5)
            + lag(pib_ph, 4:5)
            + lag(tx_deces_2ans, 4:5)
            + chomage
            + part_65,

  data           = panel,
  effect         = "individual",
  model          = "twosteps",
  transformation = "ld",
  collapse       = TRUE
)

summary(bb_model, robust = TRUE)

# Tests AR
ar1 <- mtest(bb_model, order = 1); print(ar1)
ar2 <- mtest(bb_model, order = 2); print(ar2)
ar3 <- mtest(bb_model, order = 3); print(ar3)
ar4 <- mtest(bb_model, order = 4); print(ar4)


# ══════════════════════════════════════════════════════════════════════════════
# 2. Tests de stationnarité
# ══════════════════════════════════════════════════════════════════════════════

purtest(panel$dep_sante,    test = "madwu", exo = "intercept", lags = 1)
purtest(panel$dep_sante,    test = "ips",   exo = "intercept", lags = 1)
purtest(panel[, "pib_ph"],       test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "part_65"],      test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "tx_deces_2ans"],test = "madwu", exo = "intercept", lags = 1)

panel <- panel %>%
  mutate(
    log_pib_ph  = log(pib_ph),
    d_pib_ph    = c(NA, diff(pib_ph)),
    d_part_65   = c(NA, diff(part_65)),
    log_part_65 = log(part_65)
  )

purtest(panel[, "log_pib_ph"],  test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "d_pib_ph"],    test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "d_part_65"],   test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "log_part_65"], test = "madwu", exo = "intercept", lags = 1)


# ══════════════════════════════════════════════════════════════════════════════
# 3. Blundell-Bond v5 — variables transformées (logs + différences premières)
# ══════════════════════════════════════════════════════════════════════════════

panel <- panel %>%
  mutate(
    log_dep_sante = log(dep_sante),
    log_chomage   = log(chomage),
    d_pib_ph      = c(NA, diff(pib_ph)),
    d_part_65     = c(NA, diff(part_65))
  )

panel_v5 <- subset(panel, select = c("log_dep_sante", "log_chomage",
                                     "d_pib_ph", "d_part_65", "tx_deces_2ans"))

bb_model_v5 <- pgmm(
  log_dep_sante ~ plm::lag(log_dep_sante, 1)
                + d_part_65
                + d_pib_ph
                + log_chomage
                + tx_deces_2ans

                | plm::lag(log_dep_sante, 3:5)
                + plm::lag(d_pib_ph, 2:3)
                + plm::lag(tx_deces_2ans, 2:3)
                + log_chomage
                + d_part_65,

  data           = panel_v5,
  effect         = "individual",
  model          = "twosteps",
  transformation = "ld",
  collapse       = TRUE
)

summary(bb_model_v5, robust = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# 4. Arellano-Bond (Difference GMM)
# ══════════════════════════════════════════════════════════════════════════════

ab_model <- pgmm(
  log_dep_sante ~ plm::lag(log_dep_sante, 1)
                + d_part_65
                + d_pib_ph
                + tx_deces_2ans

                | plm::lag(log_dep_sante, 2:3)
                + plm::lag(d_pib_ph, 2:3)
                + plm::lag(tx_deces_2ans, 2:3)
                + d_part_65,

  data           = panel_v5,
  effect         = "individual",
  model          = "twosteps",
  transformation = "d",
  collapse       = TRUE
)

summary(ab_model, robust = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# 5. Blundell-Bond en logs complets (avec practiciens et lits)
# ══════════════════════════════════════════════════════════════════════════════

# Recharger panel complet avec practiciens et lits
panel_full <- read.csv("panelV2.csv", sep = ",")
colnames(panel_full) <- c("pays", "annee", "dep_sante", "pib_ph", "part_65",
                          "practiciens", "lits", "tx_deces_2ans", "chomage", "dep_sante_lag")
panel_full <- pdata.frame(panel_full, index = c("pays", "annee"))
panel_full <- subset(panel_full, select = -c(dep_sante_lag))

panel_log <- panel_full %>%
  mutate(
    log_dep_sante   = log(dep_sante),
    log_pib_ph      = log(pib_ph),
    log_practiciens = log(practiciens),
    log_lits        = log(lits),
    log_part_65     = log(part_65)
  )

bb_model_log <- pgmm(
  log_dep_sante ~ plm::lag(log_dep_sante, 1)
                + log_part_65
                + log_pib_ph
                + log_practiciens
                + log_lits
                + tx_deces_2ans

                | plm::lag(log_dep_sante, 2:3)
                + plm::lag(log_pib_ph, 2:3)
                + plm::lag(tx_deces_2ans, 2:3)
                + log_practiciens
                + log_lits
                + log_part_65,

  data           = panel_log,
  effect         = "twoways",
  model          = "twosteps",
  transformation = "ld",
  collapse       = TRUE
)

summary(bb_model_log, robust = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# 6. 3SLS — système d'équations simultanées
# ══════════════════════════════════════════════════════════════════════════════

demean_twoway <- function(x, id, time) {
  grand_mean <- mean(x, na.rm = TRUE)
  mean_id    <- ave(x, id,   FUN = function(z) mean(z, na.rm = TRUE))
  mean_time  <- ave(x, time, FUN = function(z) mean(z, na.rm = TRUE))
  x - mean_id - mean_time + grand_mean
}

vars_3sls <- c("dep_sante", "pib_ph", "part_65", "practiciens", "lits", "tx_deces_2ans")

panel_dm2 <- as.data.frame(panel_full)
panel_dm2[, vars_3sls] <- lapply(
  panel_dm2[, vars_3sls],
  function(x) demean_twoway(x, id = panel_dm2$pays, time = panel_dm2$annee)
)

panel_dm2 <- pdata.frame(panel_dm2, index = c("pays", "annee"))
panel_dm2 <- panel_dm2[c("dep_sante", "pib_ph", "part_65",
                          "practiciens", "lits", "tx_deces_2ans")]
colnames(panel_dm2)[colnames(panel_dm2) == "tx_deces_2ans"] <- "tx_deces"
colnames(panel_dm2)[colnames(panel_dm2) == "part_65"]       <- "vieil"

eq_dep_sante <- as.formula("dep_sante ~ vieil + pib_ph + practiciens + lits + tx_deces")
eq_pib       <- as.formula("pib_ph ~ dep_sante + vieil + practiciens + lits")

system_eq <- list(sante = eq_dep_sante, pib = eq_pib)
instruments <- as.formula("~ vieil + lits + practiciens + tx_deces")

model_3sls <- systemfit(
  system_eq,
  method = "3SLS",
  inst   = instruments,
  data   = panel_dm2
)

summary(model_3sls)
