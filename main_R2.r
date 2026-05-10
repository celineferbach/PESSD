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
panel <- read.csv("panelV2.csv", sep = ",")
head(panel)

#on renomme pour plus de clarté
colnames(panel) <- c("pays", "annee", "dep_sante", "pib_ph", "part_65",
                     "practiciens", "lits", "tx_deces_2ans", "chomage", "dep_sante_lag")

# Déclarer la structure panel
panel <- pdata.frame(panel, index = c("pays", "annee"))
panel <- subset(panel, select = -c(dep_sante_lag))
head(panel)


# ══════════════════════════════════════════════════════════════════════════════
# Blundell-Bond (System GMM) 
# ══════════════════════════════════════════════════════════════════════════════
# pgmm() implémente Arellano-Bond (difference GMM) ET Blundell-Bond (system GMM)
# - "DPD" = dynamic panel data
# - effect = "individual" pour les effets fixes individuels
# - transformation = "ld" → "level and difference" = System GMM (Blundell-Bond)
#   (transformation = "d" seul = Arellano-Bond / difference GMM uniquement)

bb_model <- pgmm(
  dep_sante ~ plm::lag(dep_sante, 1) + plm::lag(dep_sante, 2) + plm::lag(dep_sante, 3)
            + part_65        # exogène
            + pib_ph         # endogène
            + chomage        # exogène
            + tx_deces_2ans  # endogène

             | lag(dep_sante, 4:5)     # instruments pour le lag de dep_sante
             + lag(pib_ph, 4:5)        # instruments pour pib_ph (endogène)
             + lag(tx_deces_2ans, 4:5) # instruments pour tx_deces_2ans (endogène)
             + chomage                 # exogène → instrument pour lui-même
             + part_65,                # exogène → instrument pour lui-même

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



#Tests de stationnarité

# Test de Fisher (combine des ADF individuels sur chaque pays)
purtest(panel$dep_sante,    test = "madwu", exo = "intercept", lags = 1)


purtest(panel$dep_sante,    test = "ips",   exo = "intercept", lags = 1)
purtest(panel[, "pib_ph"],       test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "part_65"],      test = "madwu", exo = "intercept", lags = 1)
purtest(panel[, "tx_deces_2ans"],test = "madwu", exo = "intercept", lags = 1)

