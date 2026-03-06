#install packages
install.packages("languageserver")
install.packages("httpgd")
install.packages("radian")
install.packages(c("plm", "pdynmc"))

#libraries
library(plm)
library(dplyr)


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
