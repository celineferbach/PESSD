#install packages
install.packages("languageserver")
install.packages("httpgd")
install.packages("radian")
install.packages(c("plm", "pdynmc"))
install.packages("systemfit")

#libraries
library(plm)
library(dplyr)
library(systemfit)


#chargement des données
panel <- read.csv("panelV2.csv", sep=",")
head(panel)


################################################## Régression ##################################################################

#on renomme les colonnes pour plus de clarté
colnames(panel) <- c("pays", "annee", "dep_sante", "pib_ph", "part_65",
                  "practiciens", "lits", "tx_deces_2ans", "chomage", "dep_sante_lag")

#déclarer la structure en panel des données
panel <- pdata.frame(panel, index = c("pays", "annee"))
panel <- subset(panel, select = -c(pays, annee, dep_sante_lag))
head(panel)


# Estimation Blundell-Bond (System GMM) avec pgmm()
# pgmm() implémente Arellano-Bond (difference GMM) ET Blundell-Bond (system GMM)
# - "DPD" = dynamic panel data
# - effect = "individual" pour les effets fixes individuels
# - transformation = "ld" → "level and difference" = System GMM (Blundell-Bond)
#   (transformation = "d" seul = Arellano-Bond / difference GMM uniquement)

bb_model <- pgmm(
  dep_sante ~ plm::lag(dep_sante, 1) + plm::lag(dep_sante, 2) + plm::lag(dep_sante, 3)
             + part_65          # exogène
             + pib_ph           # endogène
             + chomage          # exogène
             + lits
             + practiciens
             + tx_deces_2ans    # endogène (hypothèse nouvelle)

             | lag(dep_sante, 4:5)     # instruments pour le lag de dep_sante
             + lag(pib_ph, 4:5)        # instruments pour pib_ph (endogène)
             + lag(tx_deces_2ans, 4:5) # instruments pour tx_deces_2ans (endogène)
             + lits
             + practiciens
             + chomage                 # exogène → instrument pour lui-même
             + part_65,                # exogène → instrument pour lui-même

  data           = panel,
  effect         = "individual",
  model          = "twosteps",
  transformation = "ld",
  collapse = TRUE
)

summary(bb_model, robust = TRUE)

# Calculer les tests AR(1), AR(2), AR(3) et AR(4) manuellement
ar1 <- mtest(bb_model, order = 1)
ar2 <- mtest(bb_model, order = 2)
ar3 <- mtest(bb_model, order = 3)
ar4 <- mtest(bb_model, order = 4)

# Afficher les résultats
print(ar1)
print(ar2)
print(ar3)
print(ar4)



#Vérification de la stationnarité de dep_sante

# Test de Fisher (combine des ADF individuels sur chaque pays)
purtest(panel$dep_sante, test = "madwu", exo = "intercept", lags = 1)

# Ou test IPS
purtest(panel$dep_sante, test = "ips", exo = "intercept", lags = 1)

purtest(panel[, "pib_ph"],       test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "part_65"],      test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "tx_deces_2ans"], test = "madwu", exo = "intercept", lags = 1)

panel <- panel %>%
  mutate(
    log_pib_ph  = log(pib_ph),      # option 1 : log (souvent suffit à stationnariser)
    d_pib_ph    = c(NA, diff(pib_ph)) # option 2 : différence première
  )

# Vérifier si le log stationnarise
purtest(panel[, "log_pib_ph"], test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "d_pib_ph"], test = "madwu", exo = "intercept", lags = 1)

panel <- panel %>%
  mutate(
    # Option 1 : différence première (taux de vieillissement annuel)
    d_part_65 = c(NA, diff(part_65)),
    
    # Option 2 : log
    log_part_65 = log(part_65)
  )

purtest(panel[, "d_part_65"],   test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "log_part_65"], test = "madwu", exo = "intercept", lags = 1)


panel <- panel %>%
  mutate(
    log_dep_sante   = log(dep_sante),
    log_chomage = log(chomage),
    #log_lits        = log(lits),
    d_pib_ph        = c(NA, diff(pib_ph)),
    d_part_65       = c(NA, diff(part_65))
  )

panel <- subset(panel, select = c("log_dep_sante", "log_chomage", "d_pib_ph", "d_part_65", "tx_deces_2ans"))

bb_model_v5 <- pgmm(
  log_dep_sante ~ plm::lag(log_dep_sante, 1)
                + d_part_65             # différence première, exogène
                + d_pib_ph              # différence première, endogène
                + log_chomage           # exogène
                + tx_deces_2ans         # endogène, stationnaire

                | plm::lag(log_dep_sante, 3:5)
                + plm::lag(d_pib_ph, 2:3)       # instruments pour d_pib_ph
                + plm::lag(tx_deces_2ans, 2:3)  # instruments pour tx_deces_2ans
                + log_chomage
                + d_part_65,                     # exogène → instrument pour elle-même

  data           = panel,
  effect         = "individual",
  model          = "twosteps",
  transformation = "ld",
  collapse       = TRUE
)

summary(bb_model_v5, robust = TRUE)


################################################# Arellano-Bond #####################################

# Passer transformation = "d" au lieu de "ld"
# Moins efficace mais plus conservateur
ab_model <- pgmm(
  log_dep_sante ~ plm::lag(log_dep_sante, 1)
                + d_part_65
                + d_pib_ph
                + practiciens
                + lits
                + tx_deces_2ans

                | plm::lag(log_dep_sante, 2:3)
                + plm::lag(d_pib_ph, 2:3)
                + plm::lag(tx_deces_2ans, 2:3)
                + d_part_65
                + practiciens
                + lits,

  data           = panel,
  effect         = "individual",
  model          = "twosteps",
  transformation = "d",   # ← difference GMM uniquement
  collapse       = TRUE
)

summary(ab_model, robust = TRUE)



################################################ Avec des logs ? #########################################################

# Créer les variables log dans le data.frame AVANT pdata.frame
panel_raw <- read.csv("panelV2.csv", sep = ",")
colnames(panel_raw) <- c("index", "pays", "annee", "dep_sante", "pib_ph", "part_65",
                  "practiciens", "lits", "tx_deces_2ans", "chomage", "dep_sante_lag")

panel_raw <- pdata.frame(panel, index = c("pays", "annee"))
panel_raw <- subset(panel, select = -c(index, chomage, dep_sante_lag))
head(panel_raw)

panel_log <- panel_raw %>%
  mutate(
    log_dep_sante  = log(dep_sante),
    log_pib_ph     = log(pib_ph),
    log_practiciens = log(practiciens),
    log_lits       = log(lits),
    log_part_65    = log(part_65)   # optionnel
  )


# Modèle en log
bb_model_log <- pgmm(
  log_dep_sante ~ plm::lag(log_dep_sante, 1)
                + log_part_65
                + log_pib_ph
                + log_practiciens      # exogène
                + log_lits             # exogène
                + tx_deces_2ans        # endogène, pas logifié

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














######################################## 3SLS ###########################################

# Déméanisation two-way (individus + années)
demean_twoway <- function(x, id, time) {
  grand_mean  <- mean(x, na.rm = TRUE)
  mean_id     <- ave(x, id,   FUN = function(z) mean(z, na.rm = TRUE))
  mean_time   <- ave(x, time, FUN = function(z) mean(z, na.rm = TRUE))
  x - mean_id - mean_time + grand_mean
}

vars <- c("dep_sante", "pib_ph", "part_65", "practiciens", "lits", "tx_deces_2ans")

panel_dm2 <- panel
panel_dm2[, vars] <- lapply(
  panel[, vars],
  function(x) demean_twoway(x, id = panel$pays, time = panel$annee)
)

panel_dm2 <- pdata.frame(panel_dm2, index = c("pays", "annee"))
colnames(panel_dm2)[colnames(panel_dm2) == "tx_deces_2ans"] <- "tx_deces"
colnames(panel_dm2)[colnames(panel_dm2) == "part_65"] <- "vieil"
panel_dm2 <- panel_dm2[c("dep_sante", "pib_ph", "vieil", "practiciens", "lits", "tx_deces")]

# Définition du système d'équations simultanées ─────────────────────────
# On essaie avec les variables endogènes suivantes : pib et part_65

eq_dep_sante <- as.formula("dep_sante ~ vieil + pib_ph + practiciens + lits + tx_deces")
eq_pib    <- as.formula("pib_ph ~ dep_sante + vieil + practiciens + lits")

system <- list(
  sante = eq_dep_sante,
  pib    = eq_pib
)

# Les instruments sont toutes les variables exogènes de l'ensemble du système
instruments <- as.formula("~ vieil + lits + practiciens + tx_deces")

# Estimation 3SLS
model_3sls <- systemfit(
  system,
  method = "3SLS",
  inst   = instruments,
  data   = panel_dm2
)

summary(model_3sls)