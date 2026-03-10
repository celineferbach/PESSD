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
head(panel)   # voir les premières lignes


################################################## Régression ##################################################################

#on renomme les colonnes pour plus de clarté
colnames(panel) <- c("pays", "annee", "dep_sante", "pib_ph", "part_65",
                  "practiciens", "lits", "tx_deces_2ans", "dep_sante_lag")

#déclarer la structure en panel des données
panel <- pdata.frame(panel, index = c("pays", "annee"))


# Estimation Blundell-Bond (System GMM) avec pgmm()
# pgmm() implémente Arellano-Bond (difference GMM) ET Blundell-Bond (system GMM)
# - "DPD" = dynamic panel data
# - effect = "individual" pour les effets fixes individuels
# - transformation = "ld" → "level and difference" = System GMM (Blundell-Bond)
#   (transformation = "d" seul = Arellano-Bond / difference GMM uniquement)

bb_model <- pgmm(
  dep_sante ~ lag(dep_sante, 1)    # dynamique
             + part_65             # exogène (structure démographique, peu endogène)
             + pib_ph              # endogène (corrélation bidirectionnelle avec santé)
             + practiciens         # endogène (l'offre répond aux dépenses)
             + lits                # endogène (idem)
             + tx_deces_2ans     # endogène (causalité inverse possible)

             # ── Instruments ──────────────────────────────────────────────
             | lag(dep_sante, 2:4)  # instr. pour la dépendante retardée
             + lag(pib_ph, 2:4)      # instr. pour gdp_pc (endogène)
             + lag(practiciens, 2:4)  # instr. pour physicians (endogène)
             + lag(lits, 2:4)        # instr. pour beds (endogène)
             + tx_deces_2ans        # exogene ?
             + part_65,             # exogène → instrument pour elle-même (niveau)

  data           = panel,
  effect         = "individual",
  model          = "twosteps",
  transformation = "ld"
)

summary(bb_model, robust = TRUE)



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