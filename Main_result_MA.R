library(devtools, quietly = TRUE)
#devtools::install_github("jdwilson4/NetSurv")
library(blockmodels)
library(matrixStats)
library(NetSurv, quietly = TRUE)
library(Matrix, quietly = TRUE)
library(Rlab, quietly = TRUE)
library(igraph)
library(igraphdata)
library(nett)
library(randnet)
library(PRIMME)
library(mclust)
# zusätzlich RTOOL4.4 und riverplot+ seperat installieren um dynsbm 0.7 installieren zu können (, da es archiviert wurde)
library(riverplot)
library(dynsbm)
library(kableExtra)
library(webshot)
library(readr)
library(irlba)
library(gridExtra)
library(grid)
library(knitr)
library(quantmod)
library(tidyverse)

clustersigniert<-function (A, K, tau = 1, lap = FALSE, nstart = 30, iter.max = 100) 
{
  avg.d <- mean(colSums(abs(A)))
  A.tau <- abs(A) + tau * avg.d/nrow(abs(A))
  if (!lap) {
    SVD <- irlba(A.tau, nu = K, nv = K)
  }
  else {
    d.tau <- colSums(A.tau)
    L.tau <- diag(1/sqrt(d.tau)) %*% A.tau %*% diag(1/sqrt(d.tau))
    SVD <- irlba(L.tau, nu = K, nv = K)
  }
  km <- kmeans(SVD$v[, 1:K], centers = K, nstart = nstart, 
               iter.max = iter.max)
  return(list(cluster = km$cluster, loss = km$tot.withinss))
}
clustering_Guerses <- function(dynamic_networksdc, K, max_iter = 30, tau = 1, lap = FALSE, nstart = 30, iter.max = 100) {
  T <- length(dynamic_networksdc)
  n <- nrow(dynamic_networksdc[[1]])
  
  # Initialisierung
  a <- 1
  b <- 1
  c <- 1
  
  temp <- dynamic_networksdc
  
  # Berechnung von sum_neu
  neu <- list()
  for (t in 2:T) {
    neu[[t]] <- a * (temp[[t]] - temp[[t - 1]] * temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  # Persistente Kanten
  pers <- list()
  for (t in 2:T) {
    pers[[t]] <- b * (temp[[t - 1]] * temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  
  # Surviving status
  sum_surv <- c * (temp[[1]] - temp[[T]])
  
  # Berechnung von Wij
  Wij <- sum_pers + sum_neu + sum_surv
  
  # Initialisierung für Iterationen
  differenz <- list()
  abc <- list()
  abc[[1]] <- c(0, 0, 0)
  abc[[2]] <- c(a, b, c)
  
  for (ü in 3:max_iter) {
    differenz[[ü]] <- abc[[ü - 1]] - abc[[ü - 2]]
    ap <- sqrt(sum(differenz[[ü]]^2))  # Euklidische Norm
    
    if (ap > 0.0001) {
      # Normalized Spectral Clustering
      spec <- clustersigniert(Wij, K, tau = tau, lap = lap, nstart = nstart, iter.max = iter.max)
      community_assignments_estdc <- spec$cluster
      
      # Schätzung der Markov Dynamik
      counts_within <- matrix(0, 2, 2)  # Zeilen: vorheriger Zustand, Spalten: nächster Zustand
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:T) {
        prev_adj <- dynamic_networksdc[[t - 1]]
        curr_adj <- dynamic_networksdc[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_estdc[i] == community_assignments_estdc[j]) {
                counts_within[prev_adj[i, j] + 1, curr_adj[i, j] + 1] <- counts_within[prev_adj[i, j] + 1, curr_adj[i, j] + 1] + 1
              } else {
                counts_between[prev_adj[i, j] + 1, curr_adj[i, j] + 1] <- counts_between[prev_adj[i, j] + 1, curr_adj[i, j] + 1] + 1
              }
            }
          }
        }
      }
      
      P_within_est <- counts_within / rowSums(counts_within)
      P_between_est <- counts_between / rowSums(counts_between)
      
      # Berechnung neuer a, b, c Werte
      P00 <- P_within_est[1, 1] 
      P01 <- P_within_est[1, 2]
      P10 <- P_within_est[2, 1]
      P11 <- P_within_est[2, 2]
      Q00 <- P_between_est[1, 1]
      Q01 <- P_between_est[1, 2]
      Q10 <- P_between_est[2, 1]
      Q11 <- P_between_est[2, 2]
      
      a_est <- log((P10 * P01 * Q00^2) / (Q10 * Q01 * P00^2))
      b_est <- log((P11 * Q00) / (P00 * Q11))
      c_est <- (1 / (T - 1)) * log((P10 * Q00) / (Q10 * P00))
      
      Wij <- b_est * sum_pers + a_est * sum_neu + c_est * sum_surv
      abc[[ü]] <- c(a_est, b_est, c_est)
    } else {
      z <- community_assignments_estdc
      break
    }
  }
  
  return(list(community_assignments = z, P_within_est = P_within_est, P_between_est = P_between_est))
}


dax_tickers <- c(
  "ADS.DE", "AIR.DE", "ALV.DE", "BAS.DE", "BAYN.DE", "BEI.DE", 
  "BMW.DE", "BNR.DE", "CBK.DE", "CON.DE", "1COV.DE", "DHL.DE", 
  "DTG.DE", "DB1.DE", "DBK.DE", "DTE.DE", "EOAN.DE", "ENR.DE", 
  "FME.DE", "FRE.DE", "HEI.DE", "HEN3.DE", "HNR1.DE", "IFX.DE", 
  "MBG.DE", "MRK.DE", "MTX.DE", "MUV2.DE", "P911.DE", "PAH3.DE", 
  "RWE.DE", "RHM.DE", "SAP.DE", "SIE.DE", "SY1.DE", "ENR.DE", 
  "VOW3.DE", "VNA.DE", "QIA.DE", "ZAL.DE"
)

# Parameter
start_date <- as.Date("2024-01-01")
num_weeks <- 7
threshold <- 0.4

# Liste zur Speicherung der Adjazenzmatrizen
A_dax_list <- list()

# Automatisierung über die ersten 10 Wochen
for (week in 1:num_weeks) {
  # Zeitraum berechnen
  week_start <- start_date + (week - 1) * 7
  week_end <- week_start + 6
  
  # Daten laden
  getSymbols(dax_tickers, src = "yahoo", from = week_start, to = week_end, auto.assign = TRUE)
  
  # Schlusskurse zusammenführen und tägliche Renditen berechnen
  prices <- do.call(merge, lapply(dax_tickers, function(x) Cl(get(x))))
  
  # Sicherstellen, dass es genug Daten gibt
  if (nrow(prices) > 1) {
    returns <- na.omit(diff(log(prices)))
    
    # Korrelationen berechnen
    correlation_matrix <- cor(returns)
    
    # Kantenliste basierend auf Korrelationen über dem Schwellenwert erstellen
    edges <- data.frame(from = character(), to = character(), weight = numeric(), stringsAsFactors = FALSE)
    for (i in 1:(ncol(correlation_matrix)-1)) {
      for (j in (i+1):ncol(correlation_matrix)) {
        if ((correlation_matrix[i, j]) > threshold) {
          edges <- rbind(edges, data.frame(from = colnames(correlation_matrix)[i], to = colnames(correlation_matrix)[j], weight = correlation_matrix[i, j]))
        }
      }
    }
    
    # Graph erstellen und Adjazenzmatrix speichern
    gdax <- graph_from_data_frame(edges, directed = FALSE)
    A_dax_list[[week]] <- as.matrix(as_adjacency_matrix(gdax))
  } else {
    message(paste("Nicht genug Daten für Woche", week))
  }
}

M1<-A_dax_list[[1]]
M2<-A_dax_list[[2]]
M3<-A_dax_list[[3]]
M4<-A_dax_list[[4]]
M5<-A_dax_list[[5]]
M6<-A_dax_list[[6]]
M7<-A_dax_list[[7]]
vemdax <-BM_bernoulli_multiplex("SBM", list(M1,M2,M3,M4,M5,M6,M7))
vemdax$estimate()
which.max(vemdax$ICL)

Clust<-clustering_Guerses(A_dax_list,2)

membership<-Clust$community_assignments
palette <- rainbow(max(membership))
igraph_networksdax <-
  lapply(A_dax_list, function(adj_matrix) {
    graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  })
par(mfrow = c(2, 3)) # Adjust based on the number of time steps
#alle Snapshots in einem Bild:
set.seed(17)
for (t in 1:5) {
  plot(
    igraph_networksdax[[t]],
    main = paste("Zeitpunkt", t),
    vertex.size = 5,
    vertex.label = NA,
    vertex.color = membership
  )
}
