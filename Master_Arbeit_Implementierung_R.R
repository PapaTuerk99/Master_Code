############################################################
# Implementierungen meiner Masterarbeit (und ein paar extras)
############################################################
# Alle benötigten Pakete
# install.packages("devtools")
# install.packages("blockmodels")
# install.packages("matrixStats")
# install.packages("Matrix")
# install.packages("Rlab")
# install.packages("igraph")
# install.packages("igraphdata")
# install.packages("nett")
# install.packages("randnet")
# install.packages("PRIMME")
# install.packages("mclust")
# install.packages("riverplot")
# install.packages("dynsbm")
# install.packages("kableExtra")
# install.packages("webshot")
# install.packages("readr")
# install.packages("irlba")
# installiere NetSurv
# install.packages("quantmod")
# install.packages("tidyverse")

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

############################################################
# Kapitel 2.6 Beispiele
# 2.6.1
# R-Paket randnet:
set.seed(17)
dtf <-
  BlockModel.Gen(
    8,
    60,
    Pi = cbind(17, 25, 18),
    K = 3,
    beta = 0.10,
    rho = 0.9,
    simple = FALSE,
    power = FALSE
  )
Af <- Matrix(dtf$A)
gf <- graph.adjacency (Af , mode = "undirected")
plot(gf, vertex.color = dtf$g) 
#SBM:
set.seed(17)
sbm <-
  BlockModel.Gen(
    8,
    60,
    Pi = cbind(17, 25, 18),
    K = 3,
    beta = 0.1,
    rho = 0,
    simple = FALSE,
    power = FALSE
  )
Asbm <- Matrix(sbm$A)
g <- graph.adjacency (Asbm , mode = "undirected")
plot(g, vertex.color = sbm$g, )
# Pi ist Blockgrößen; 8 ist expected degree; K ist die Anzahl an Blöcken; rho=0 bedueted SBM;
# power False bedeuted U[0.2;1] falls rho>0.
# beta ist die Propotion zu p_out/p_in
##############################
# Beispiel in Abschnitt 2.6.2
# Poisson-DCSBM ohne Schleife (Funktion von netSurv übernommen):
Poisson_DCSBM <-function (n,
            k = 2,
            P,
            sizes = c(round(n / 2), n - round(n / 2)),
            random.community.assignment = c(FALSE, TRUE),
            community.labels = NULL,
            delta = rep(0, k),
            edge.list = c(FALSE, TRUE))
  {
    if (sum(sizes) != n) {
      stop("argument sizes must sum to n")
    }
    if (length(sizes) != k) {
      stop("argument sizes must be of length k")
    }
    if (!is.null(community.labels)) {
      Membership <- as.numeric(as.factor(community.labels))
      k <- length(unique(Membership))
      Y <- matrix(rep(0, n * k), ncol = k)
      for (i in 1:k) {
        Y[which(Membership == i), i] = 1
      }
    }
    if (is.null(community.labels)) {
      Y <- matrix(rep(0, n * k), ncol = k)
      index <- list()
      possible <- 1:n
      Membership <- rep(1, n)
      random.community.assignment <- random.community.assignment[1]
      if (random.community.assignment == TRUE) {
        for (i in 1:k) {
          index[[i]] <- sample(possible, sizes[i])
          Membership[index[[i]]] = i
          possible <- setdiff(possible, index[[i]])
        }
      }
      if (random.community.assignment == FALSE) {
        for (i in 1:k) {
          index[[i]] <- possible[1:sizes[i]]
          Membership[index[[i]]] = i
          possible <- setdiff(possible, index[[i]])
        }
      }
      for (i in 1:k) {
        Y[index[[i]], i] = 1
      }
    }
    edge.list <- edge.list[1]
    if (length(which(delta > 1)) > 0) {
      stop("argument delta parameters must be less than or equal to 1")
    }
    thetas <- rep(1, n)
    for (i in 1:k) {
      thetas[Membership == i] <- runif(sizes[i], min = 1 -
                                         delta[i], max = 1 + delta[i])
    }
    theta.matrix <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        theta.matrix[i, j] <- thetas[i] * thetas[j]
      }
    }
    expected.A <- Y %*% P %*% t(Y)
    expected.A <- theta.matrix * expected.A
    temp <- matrix(rpois(n ^ 2, matrix(expected.A, ncol = 1)),
                   ncol = n)
    Adj <- matrix(0, ncol = n, nrow = n)
    Adj[upper.tri(Adj)] <- temp[upper.tri(temp)]
    Adj <- Adj + t(Adj)
    diag(Adj) <- 0
    Adj <- Matrix(Adj, sparse = TRUE)
    if (edge.list == TRUE) {
      return(list(
        Adjacency = as.vector(Adj),
        Thetas = thetas,
        Membership = Membership
      ))
    }
    if (edge.list == FALSE) {
      return(list(
        Adjacency = Adj,
        Thetas = thetas,
        Membership = Membership
      ))
    }
  }
set.seed(17)
net <- Poisson_DCSBM(
    n = 50,
    k = 2,
    P = cbind(c(0.3, 0.10), c(0.10, 0.30)),
    sizes = c(22, 28),
    random.community.assignment = FALSE,
    delta = c(0.2, 0.7),
    edge.list = TRUE
  )
net
mat <-
  matrix(net$Adjacency, nrow = 50, ncol = 50) #Vektor zur Matrix machen
gi <-
  graph.adjacency(mat, mode = "undirected") #Adjazenzmatrix generieren zum plotten
plot(gi, vertex.color = net$Membership) #Plotten mit Membership Färbung



################################################################
# Abschnitt 2.6.3:
# temporales Bernoulli SBM mit Markov-Dynamik für statisches z:
# kleiner Graph mit T=6 n=20 K=2 p_in=0,3 und p_out=0.1.Markov-Dynamik ist P_within und P_between
# Parameter:
n <- 20 # Knotenanzahl
k <- 2  # Anzahl der Blöcke
T <- 6 # Anzahl der Snapshots
# Statische community Zugehörigkeit:
set.seed(17)
community_assignments <- sample(1:k, n, replace = TRUE)
# Markov Matrix:
# z_i=z_j community:
P_within <- matrix(c(0.9, 0.1,
                     0.1, 0.9),
                   nrow = 2, ncol = 2)
# Verschiedene communities
P_between <- matrix(c(0.95, 0.05,
                      0.05, 0.95),
                    nrow = 2, ncol = 2)
# homogenes Bernoulli SBM mit p_in=0,3 ; p_out=0,1
adjacency_matrix <- matrix(0, n, n)
set.seed(17)
for (i in 1:n) {
  for (j in 1:i) {
    if (i != j) {
      if (community_assignments[i] == community_assignments[j]) {
        prob <- 0.3  # also p_in
      } else {
        prob <- 0.1  # p_out
      }
      adjacency_matrix[i, j] <- rbinom(1, 1, prob)
      adjacency_matrix[j, i] <-
        adjacency_matrix[i, j] #Symmetrie garantieren
    }
  }
}
# Listen zum speichern:
dynamic_networks <- list()
dynamic_networks[[1]] <- adjacency_matrix
set.seed(17)
for (t in 2:T) {
  previous_adjacency_matrix <- dynamic_networks[[t - 1]]
  current_adjacency_matrix <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:i) {
      if (i != j) {
        if (community_assignments[i] == community_assignments[j]) {
          transition_probs <- P_within[previous_adjacency_matrix[i, j] + 1,]
        } else {
          transition_probs <- P_between[previous_adjacency_matrix[i, j] + 1,]
        }
        current_state <- sample(0:1, 1, prob = transition_probs)
        current_adjacency_matrix[i, j] <- current_state
        current_adjacency_matrix[j, i] <- current_state
      }
    }
  }
  dynamic_networks[[t]] <- current_adjacency_matrix
}
igraph_networks <- lapply(dynamic_networks, function(adj_matrix) {
  graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
})
par(mfrow = c(2, 3)) # Adjust anhand von T
# alle Snapshots in einem Schaubild:
for (t in 1:T) {
  plot(
    igraph_networks[[t]],
    main = paste("Zeitpunkt", t),
    vertex.size = 5,
    vertex.label = NA,
    vertex.color = community_assignments
  )
}


###############################
# temporales Bernoulli DC-SBM mit Markov-Dynamik für statisches z:
# kleiner Graph mit T=6 n=20 K=2 p_in=0,3 und p_out=0.1.Markov-Dynamik ist P_within und P_between
# Zusätzlich Theta uniform zwischen 0.6 und 1.4
n <- 20
K <- 2
T <- 6
# Community-Zuordnungen (Beispiel: zufällige Zuordnung)
set.seed(23)
community_assignmentsdc <- sample(1:K, n, replace = TRUE)

# Theta-Werte für jeden Knoten (Beispiel: zufällige Werte)
set.seed(23)
theta <- runif(n, 0.5, 1.5)
# Markov-Übergangsmatrizen
P_within <- matrix(c(0.9, 0.1,
                     0.1, 0.9),
                   nrow = 2, ncol = 2)
# Verschiedene communities
P_between <- matrix(c(0.95, 0.05,
                      0.05, 0.95),
                    nrow = 2, ncol = 2)
# Initiale Adjazenzmatrix
set.seed(23)
initial_adj <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      if (community_assignmentsdc[i] == community_assignmentsdc[j]) {
        initial_adj[i, j] <- rbinom(1, 1, theta[i] * theta[j] * 0.3) #p_in
      } else {
        initial_adj[i, j] <- rbinom(1, 1, theta[i] * theta[j] * 0.1) #p_out
      }
      initial_adj[j, i] <- initial_adj[i, j]
    }
  }
}
# Liste der dynamischen Netzwerke
dynamic_networksdc <- list()
dynamic_networksdc[[1]] <- initial_adj

# Simulation der dynamischen Netzwerke
set.seed(23)
for (t in 2:T) {
  prev_adj <- dynamic_networksdc[[t - 1]]
  curr_adj <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        if (community_assignmentsdc[i] == community_assignmentsdc[j]) {
          if (prev_adj[i, j] == 1) {
            prob <- P_within[2, 2]
          } else {
            prob <- P_within[1, 2]
            prob <- prob * theta[i] * theta[j]
          }
          if (prob > 1) {
            prob <- 1
          }
          curr_adj[i, j] <- rbinom(1, 1, prob)
        } else {
          if (prev_adj[i, j] == 1) {
            prob <- P_between[2, 2]
          } else {
            prob <- P_between[1, 2]
            prob <- prob * theta[i] * theta[j]
          }
          if (prob > 1) {
            prob <- 1
          }
          curr_adj[i, j] <- rbinom(1, 1, prob)
        }
        curr_adj[j, i] <- curr_adj[i, j]
      }
    }
  }
  dynamic_networksdc[[t]] <- curr_adj
}
igraph_networksdc <-
  lapply(dynamic_networksdc, function(adj_matrix) {
    graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  })
par(mfrow = c(2, 3)) # Adjust based on the number of time steps
#alle Snapshots in einem Bild:
for (t in 1:T) {
  plot(
    igraph_networksdc[[t]],
    main = paste("Zeitpunkt", t),
    vertex.size = 5,
    vertex.label = NA,
    vertex.color = community_assignmentsdc
  )
}

expected_degree_from_adj_matrices <- function(adj_matrices, node) {
  # Initialisieren der Liste zur Speicherung der Grade
  degrees <- c()
  
  # Iteriere über jede Adjazenzmatrix
  for (adj_matrix in adj_matrices) {
    # Berechne den Grad des Knotens in dieser Adjazenzmatrix
    degree <- sum(adj_matrix[node,])
    
    # Füge den Grad zur Liste hinzu
    degrees <- c(degrees, degree)
  }
  
  # Berechne den Durchschnitt der Grade
  expected_deg <- mean(degrees)
  return(expected_deg)
}
exp_sbm <- list()
for (r in 1:20) {
  expect <- expected_degree_from_adj_matrices(dynamic_networks, r)
  exp_sbm[r] <- expect
}
exp_dcsbm <- list()
for (r in 1:20) {
  expect <- expected_degree_from_adj_matrices(dynamic_networksdc, r)
  exp_dcsbm[r] <- expect
}
exp_sbm_vec <- unlist(exp_sbm)
exp_dcsbm_vec <- unlist(exp_dcsbm)
# empirische Varianz und Durchschnitt:
mean(exp_sbm_vec)
var(exp_sbm_vec)
mean(exp_dcsbm_vec)
var(exp_dcsbm_vec)
mean(theta)
##########################################################
# Beispiele aus dem Abschnitt 3.1
# Bethe Hessian K schätzen wie im Paper mit t=5 !!!
# Übernommen von BHMC.estimate, damit man t nach belieben anpassen kann
BH <- function (A, K.max = 15)
{
  if (K.max <= 2)
    K.max <- 2
  d <- colSums(A)
  n <- nrow(A)
  I <- as(diag(rep(1, n)), "dgCMatrix")
  D <- as(diag(d), "dgCMatrix")
  r <- sqrt(mean(d))
  BH <- (r ^ 2 - 1) * I - r * A + D
  rho <- sort(eigs_sym(BH, K.max, which = "SA")$values)
  diff <- rho[2:K.max] - 5 * rho[1:(K.max - 1)]
  return(list(K = max(which(diff > 0)), values = rho))
}
#Test mit Poisson DC-SBM (2.6.2) aus Bsp 3.3:
BH(mat, K.max = 9)

#Non backtracking matrix für (2.6.2):
get.backtrack <- function(g) {
  n <- ecount(g)
  g <- as.directed(g)
  h <- line.graph(g)
  P <- rbind(1:n, (n + 1):(2 * n))
  P <- cbind(P, P[c(2, 1), ])
  h <- h - edges(E(h, P))
  get.adjacency(h)
}
x <- get.backtrack(gi)
e <- eigen(x, only.values = TRUE)$values
spektralradius <- sqrt(max(abs(e)))
di <- e - spektralradius
filtered_vector <- Re(di[Re(di) > 0])
K <- length(filtered_vector)
show(K)

########################################################
# Abschnitt 3.5:
# Modelwahl: Homogenes Poisson DC-SBM erstellen mit n=500, (220,280) Blockgröße, w=(0,3 0,1)
# Beispiel 3.15:
set.seed(17)
nets <- Poisson_DCSBM(
    n = 500,
    k = 2,
    P = cbind(c(0.3, 0.10), c(0.10, 0.30)),
    sizes = c(220, 280),
    random.community.assignment = FALSE,
    delta = c(0.5, 0.5),
    edge.list = TRUE
  )
#Test Normierung für Identifizierbarkeit des Modells Abweichung:
for (i in 1:500) {
  if (nets$Membership[i] == 1) {
    nets$Thetas[i] <- nets$Thetas[i] / 220
  } else{
    nets$Thetas[i] <- nets$Thetas[i] / 280
  }
}
matt <- matrix(nets$Adjacency, nrow = 500, ncol = 500)
#graphmatt <- graph.adjacency(matt, mode ="undirected")
#plot(graphmatt)
#Bethe Hessian für t=7 um K_max zu schätzen:
# t=7 um zu zeigen, dass die spektrale Lücke groß ist.
BHmax <- function (A, K.max = 15)
{
  if (K.max <= 2)
    K.max <- 2
  d <- colSums(A)
  n <- nrow(A)
  I <- as(diag(rep(1, n)), "dgCMatrix")
  D <- as(diag(d), "dgCMatrix")
  r <- sqrt(mean(d))
  BH <- (r ^ 2 - 1) * I - r * A + D
  rho <- sort(eigs_sym(BH, K.max, which = "SA")$values)
  diff <- rho[2:K.max] - 7 * rho[1:(K.max - 1)]
  return(list(K = max(which(diff > 0)), values = rho))
}
K_max <- BHmax(matt, K.max = 10)
K_max <- K_max$K
# K wird als 2 geschätzt
# Nun ECV Methode aus dem R-Paket randnet:
# set.seed(17)
# Selection<-ECV.block(matt, max.K=K_max, cv = NULL, B = 10, holdout.p = 0.1, tau = 0, dc.est = 2, kappa = NULL)
# Selection$l2.model
# DC-SBM mit K=2 als Resultat
# Test für K_max=5
set.seed(17)
Select <-ECV.block(
    matt,
    max.K = 5,
    cv = NULL,
    B = 10,
    holdout.p = 0.1,
    tau = 0,
    dc.est = 2,
    kappa = NULL
  )
Select$l2
Select$dc.l2
Select$l2.model
################
# Beispiel 3.17:
# Robustheit: 20 mal testen ob die Methode das richtige Ergebnisse liefert:
# Hinweis: Es kann ein paar Minuten dauern
ergrob <- list()
set.seed(17)
for (x in 1:20) {
  Status <-ECV.block(
      matt,
      max.K = 5,
      cv = NULL,
      B = 10,
      holdout.p = 0.1,
      tau = 0,
      dc.est = 2,
      kappa = NULL
    )
  ergrob[[x]] <- Status$l2.model
}
print(ergrob)
################
# Beispiel 3.18:
# Karate Club laden:
karate_club <- graph.famous("Zachary")
# plot.igraph(karate_club, main = "Zachary Karate Club Network")
# wahre Spaltung:
ground_truth <-
  c(1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    1,
    1,
    1,
    1,
    2,
    2,
    1,
    1,
    2,
    1,
    2,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2)
colors <- ifelse(ground_truth == 1, "green", "red")
# Plot mit wahren Communitys:
set.seed(17)
plot(karate_club, vertex.color = colors, main = "Karate Klub mit den echten Blöcken")
# ECV bei kleinem Graphen:
# Adjazenzmatrix erzeugen:
adj_matrix_karate <-
  as.matrix(get.adjacency(karate_club, sparse = FALSE))
ergebnisse2 <- list()
set.seed(17)
for (y in 1:20) {
  Status2 <-ECV.block(
      adj_matrix_karate,
      max.K = 4,
      cv = NULL,
      B = 10,
      holdout.p = 0.3,
      tau = 0,
      dc.est = 2,
      kappa = NULL
    )
  ergebnisse2[[y]] <- Status2$l2.model
}
print(ergebnisse2)
# Zum Vergleich, Bethe Hessian Schätzer erkennt echtes K:
BH(adj_matrix_karate)
####################
# Beispiel 3.19:
#schwächen Analyse für zu kleine synthetische Graphen: (n=100 aber p=0,2)
ergebnisse<-list()
ergebnissetotal<-list()
set.seed(17)
for(q in 1:20){
netstest <-Poisson_DCSBM(n = 100, k = 2, P = cbind(c(0.3, 0.10), c(0.10, 0.30)),
                         sizes = c(42, 58), random.community.assignment = FALSE,
                         delta = c(0.5, 0.5), edge.list = TRUE)
for(i in 1:100){
  if(netstest$Membership[i]==1){
    netstest$Thetas[i]<-netstest$Thetas[i]/42
  } else{netstest$Thetas[i]<-netstest$Thetas[i]/58}
}  
mattest <- matrix(netstest$Adjacency, nrow = 100, ncol=100)
for(x in 1:20){
Status <- ECV.block(mattest, max.K=5, cv = NULL, B = 10, holdout.p = 0.2, tau = 0, dc.est = 2, kappa = NULL)
ergebnisse[[x]]<-Status$l2.model
}
ergebnissetotal[[q]]<-ergebnisse
}
print(ergebnissetotal)
# Häufigstes Modell anzeigen:
häufigstes_element <- function(liste) {
  häufigkeiten <- table(unlist(liste))
  häufigstes <- names(häufigkeiten)[which.max(häufigkeiten)]
  return(häufigstes)
}
# Das häufigste Element in jeder Liste der Liste finden
häufigste_ergebnisse <- sapply(ergebnissetotal, häufigstes_element)
# Ausgabe der häufigsten Ergebnisse:
print(häufigste_ergebnisse)



##################################################################################
#Kapitel 4.1.1:
# Datenbeispiel für nicht evolvierenden Graphen:
# Kerninghan-Lin EM ist in phython zu finden
# fremde Datenwerk network (Dolphin laden) von https://networks.skewed.de/net/dolphins heruntergeladen
# Working Dirctory einstellen!!!!
g<-read.graph("dolphin.gml",format=c("gml"))
adj_matrix <- as.matrix(get.adjacency(g))
# Anzahl der Blöcke K schätzen mittels ICL:
set.seed(17)
vemdolph <- BM_bernoulli("SBM", adj_matrix )
set.seed(17)
vemdolph$estimate()
which.max(vemdolph$ICL)
#vemdolph <- BM_poisson("SBM", adj_matrix )
#vemdolph$estimate()
#which.max(vemdolph$ICL)
BH(adj_matrix, 6) #2 
# ECV für p=0,3 :
#set.seed(17)
#ED<-ECV.block(adj_matrix, max.K=5, cv = NULL, B = 10, holdout.p = 0.3, tau = 0, dc.est = 2, kappa = NULL)
#ED$l2.model
# ECV versagt hier beispielsweise
# Bild mit Graph und Graph mit ground truth Färbung
par(mfrow = c(1, 2))
set.seed(17)
plot(g, vertex.size = 5, vertex.label = NA, main="Ursprungsgraph")
# manuell die wahres z einpflegen
ground_truth_clusters <- c(
  1, 2, 1, 3, 4, 2, 2, 2, 3, 2, 1, 4, 3, 2, 3, 4, 3, 2,
  4, 2, 3, 4, 2, 4, 4, 2, 2, 2, 1, 4, 1, 2, 2, 3, 3, 4, 
  3, 3, 3, 3, 3, 2, 1, 3, 3, 4, 3, 1, 2, 3, 3, 4, 3, 3,
  2, 4, 2, 2, 3, 3, 2, 3
)
palette <- rainbow(max(ground_truth_clusters))
set.seed(17)
plot(g,
     vertex.color = palette[ground_truth_clusters],
     vertex.size = 5,
     vertex.label = NA,
     edge.arrow.size = 0.5,
     main = "Graph mit wahrer Färbung")
#Schaubild mit den ganzen Clustering erzeugen:
#V-EM:
par(mfrow=c(1,4))
sbm_dolph_hat <- apply(vemdolph$memberships[[4]]$Z, 1, which.max)
adj_plot<-graph.adjacency(adj_matrix , mode ="undirected")
palette <- rainbow(max(sbm_dolph_hat))
VEMdolph<- adjustedRandIndex(sbm_dolph_hat, ground_truth_clusters)
set.seed(17)
plot(g,
     vertex.color = palette[sbm_dolph_hat],
     vertex.size = 5,
     vertex.label = NA,
     edge.arrow.size = 0.5,
     main = "V-EM Clustering")
#normalized Laplace Clustering:
set.seed(17)
j<-reg.SP(adj_matrix, 4, tau = 0.5, lap = TRUE,nstart=30,iter.max=100) 
membership<-j$cluster
palette <- rainbow(max(membership))
Lapdolph<- adjustedRandIndex(membership, ground_truth_clusters)
set.seed(17)
plot(g,
     vertex.color = palette[membership],
     vertex.size = 5,
     vertex.label = NA,
     edge.arrow.size = 0.5,
     main = "Laplace Clustering")
#Louvain Algorithmus zur Modularity Maximierung:
set.seed(17)
communities <- cluster_louvain(g,NULL,1)
membership <- membership(communities)
palette <- rainbow(max(membership))
Loudolph<- adjustedRandIndex(membership, ground_truth_clusters)
set.seed(17)
plot(g,
     vertex.color = palette[membership],
     vertex.size = 5,
     vertex.label = NA,
     edge.arrow.size = 0.5,
     main = "Clustering mit Louvain-Algorithmus")
#Leiden im Vergleich:
set.seed(17)
communities<-cluster_leiden(g,objective_function = "modularity")
membership <- membership(communities)
palette <- rainbow(max(membership))
Leidolph<- adjustedRandIndex(membership, ground_truth_clusters)
set.seed(17)
plot(g,
     vertex.color = palette[membership],
     vertex.size = 5,
     vertex.label = NA,
     edge.arrow.size = 0.5,
     main = "Clustering mit Leiden-Algorithmus")
#Tabelle erzeugen:
#Übersicht der Genauigkeiten: (VEM statt V-EM, da sonst Fehlermeldung kommt)
df <- data.frame(Louvain = Loudolph , Laplace = Lapdolph, Leiden = Leidolph, VEM = VEMdolph)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

# Speichern als PDF in der working directory:
pdf_file <- "tabelle_dolph.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()

#ARI vom Kerninghan-Lin (siehe Phython:)
#Manuelle Übertragung der Community Erkennungen: (1: Dunkel orange; 2: hellblau; 3: dunkelblau; 4: helles orange) 
Newman_sbm<- c(3, 1, 3, 3, 4, 2, 3, 1, 3, 4, 3,
               4, 4, 2, 2, 2, 1, 2, 4, 1, 3, 2,
               4, 2, 2, 4, 4, 3, 3, 2, 3, 4, 2,
               2, 2, 4, 2, 2, 2, 4, 3, 1, 3, 2,
               3, 2, 4, 3, 4, 4, 1, 2, 4, 3, 1,
               4, 4, 2, 4, 4, 4, 4)

ARINewmanSBM<-adjustedRandIndex(Newman_sbm, ground_truth_clusters) # 0.063
Newman_dcsbm<- c(4, 3, 4, 4, 2, 3, 3, 3, 3, 3, 4,
                 2, 2, 3, 2, 2, 2, 3, 2, 3, 4, 4,
                 3, 4, 4, 3, 3, 3, 1, 2, 1, 3, 3,
                 4, 4, 3, 3, 4, 2, 3, 2, 3, 4, 4,
                 3, 4, 2, 4, 3, 4, 2, 4, 4, 4, 3,
                 4, 3, 3, 3, 4, 3, 3)
ARINewmanDCSBM<-adjustedRandIndex(Newman_dcsbm, ground_truth_clusters) # 0.325

######################################################################################
# Abschnitt 4.1.2:
# synthetische Beispiele:
# nicht evolvierender Fall:
# simpler homogener bernoulli SBM mit statischen Communitys:
# Anzahl der Knoten in jeder Community festlegen:
community_sizes <- c(17, 25, 18)
# Gesamtanzahl der Knoten
num_nodes <- sum(community_sizes)
# Wahrscheinlichkeitsmatrix für Verbindungen innerhalb und zwischen den Communities
p_matrix <- matrix(c(0.3, 0.05, 0.05,
                     0.05, 0.3, 0.05,
                     0.05, 0.05, 0.3), 
                   nrow = 3, byrow = TRUE)

# Erzeuge einen zufälligen SBM-Graphen
# set.seed(17)
# sbm_graph <- sample_sbm(num_nodes, pref.matrix = p_matrix, block.sizes = community_sizes)
# plot(sbm_graph)
# Erstelle die Adjazenzmatrix des Graphen
# adjacency_matrix <- as.matrix(as_adjacency_matrix(sbm_graph))

# Definiere die Community-Zugehörigkeiten der Knoten
communities <- rep(1:length(community_sizes), community_sizes)

# Plotte den vereinfachten Graphen und färbe die Knoten nach den Communities
#plot(sbm_graph,
#     vertex.color = rainbow(length(community_sizes))[communities],
#     vertex.size = 10,
#    edge.arrow.size = 0.5,
#    main = "Vereinfachter Stochastisches Blockmodell (SBM) Netzwerk mit festen Communities")

# Nun folge von diesen Matrizen generieren und adjaznzmatrizen in einer Liste speichern:
temp<-list()
AnzahlSnap<-10
set.seed(17)
for(t in 1:AnzahlSnap){
g<-sample_sbm(num_nodes, pref.matrix = p_matrix, block.sizes = community_sizes)
mat<-as.matrix(as_adjacency_matrix(g))
temp[[t]] <- mat
}

#Blockanzahl schätzen ICL-Kriterium Blockmodell Paket und V-EM Clustering
#for (i in 1:length(temp)) {
# assign(paste0("M", i), temp[[i]])}
# einzelne SBM Testen:
#set.seed(17)
#test<-BM_bernoulli("SBM",M1)
#set.seed(17)
#test$estimate()
#sbmtest_Z_hat <- apply(test$memberships[[3]]$Z, 1, which.max)
#adjustedRandIndex(sbmtest_Z_hat, communities)


#set.seed(17)
#my_model <- BM_bernoulli_multiplex("SBM",list(M1,M2,M3,M4,M5))
#set.seed(17)
#my_model$estimate()
#which.max(my_model$ICL)

 
# Laplace Clustering der Daten und berechnung des ARI:
ariLap<-list()
set.seed(17)
for(t in 1:AnzahlSnap){
  stemp<-temp[[t]]
  h<-reg.SP(stemp, 3, tau = 0.5, lap = TRUE,nstart=30,iter.max=100)   
  h<-h$cluster
  ariLap[[t]] <- adjustedRandIndex(communities, h)
  }
# Durchschnittliche ARI:
#arivlaplace<-as.vector(ariLap)
meanvekLap <- sapply(ariLap, mean)
Laplace<-mean(meanvekLap)

#Analog für Louvain:
arilou<-list()
set.seed(17)
for(t in 1:AnzahlSnap){
  stemp<-temp[[t]]
  graph <- graph_from_adjacency_matrix(stemp, mode = "undirected")
   h<-cluster_louvain(graph, NULL, 0.5)   
  membership1 <- membership(h)
  arilou[[t]] <- adjustedRandIndex(communities, membership1)
}
#arivlouvain<-as.vector(arilou)
meanvekLou <- sapply(arilou, mean)
Louvain<-mean(meanvekLou)
# Leiden alg
arilei<-list()
set.seed(17)
for(t in 1:AnzahlSnap){
  stemp<-temp[[t]]
  graph <- graph_from_adjacency_matrix(stemp, mode = "undirected")
  h<-cluster_leiden(graph, objective_function = "modularity")   
  membership1 <- membership(h)
  arilei[[t]] <- adjustedRandIndex(communities, membership1)
}
#arivleiden<-as.vector(arilei)
meanvek <- sapply(arilei, mean)
Leiden<-mean(meanvek)
###
# Übersicht der einzelnen ARI erstellen in Form einer Tabelle:
# Für die Übersicht runden wir auf die 5. Nachkomma Stelle
# Erstellen eines DataFrames
df <- data.frame(Louvain = meanvekLou , Laplace = meanvekLap, Leiden = meanvek)
# Runden der Werte auf die fünfte Dezimalstelle
df <- round(df, 5)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
pdf_file <- "tabelle_sbm.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()



############
# Übersicht der Durchschnitte (SBM)
df <- data.frame(Louvain = Louvain , Laplace = Laplace, Leiden = Leiden)
df <- round(df, 5)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
pdf_file <- "tabelle_sbm_mean.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()


##########################################################
# DC-SBM via sample_dcsbm aus paket nett:
Adcsbm<-list()
AnzahlSnap<-10
set.seed(17)
for(w in 1:AnzahlSnap){
k<-3
theta_range <- c(0.5, 1.5)
n<-60
p_in<- 0.3
p_out<- 0.05
community_sizesdc <- c(17, 25, 18)
  # Generiere Blockzuweisungen fixiert aber
  block_assignmentsdc<-rep(1:length(community_sizesdc), community_sizesdc)
  # Generiere Theta-Werte gleichmäßig verteilt im angegebenen Bereich
  theta <- runif(n, min = theta_range[1], max = theta_range[2])
  
  # Erstelle die Blockwahrscheinlichkeitsmatrix
  B <- matrix(p_out, nrow = k, ncol = k)
  diag(B) <- p_in
  
  # Erzeuge die Adjazenzmatrix mit sample_dcsbm
  A <- sample_dcsbm(block_assignmentsdc, B, theta)
  Aad<-as.matrix(A)
  Adcsbm[[w]]<-Aad
}
#Laplace Clustering der Daten und berechnung des ARI:
ariLapdc<-list()
set.seed(17)
for(t in 1:AnzahlSnap){
  stemp<-Adcsbm[[t]]
  h<-reg.SP(stemp, 3, tau = 0.5, lap = TRUE,nstart=30,iter.max=100)   
  h<-h$cluster
  ariLapdc[[t]] <- adjustedRandIndex(block_assignmentsdc, h)
}
# Durchschnittliche ARI:
# arivlaplacedc<-as.vector(ariLapdc)
meanvekLapdc <- sapply(ariLapdc, mean)
Laplacedc<-mean(meanvekLapdc)
#Analog für Louvain:
ariloudc<-list()
set.seed(17)
for(t in 1:AnzahlSnap){
  stemp<-Adcsbm[[t]]
  graph <- graph_from_adjacency_matrix(stemp, mode = "undirected")
  h<-cluster_louvain(graph, NULL, 0.5)   
  membership1 <- membership(h)
  ariloudc[[t]] <- adjustedRandIndex(block_assignmentsdc, membership1)
}
#arivlouvaindc<-as.vector(ariloudc)
meanvekLoudc <- sapply(ariloudc, mean)
Louvaindc<-mean(meanvekLoudc)
# Leiden alg
arileidc<-list()
set.seed(17)
for(t in 1:AnzahlSnap){
  stemp<-Adcsbm[[t]]
  graph <- graph_from_adjacency_matrix(stemp, mode = "undirected")
  h<-cluster_leiden(graph, objective_function = "modularity")   
  membership1 <- membership(h)
  arileidc[[t]] <- adjustedRandIndex(block_assignmentsdc, membership1)
}
#arivleidendc<-as.vector(arileidc)
meanvekdc <- sapply(arileidc, mean)
Leidendc<-mean(meanvekdc)

# Übersicht erstellen der einzelnen ARIs in Form einer Tabelle:
# Erstellen eines DataFrames
df <- data.frame(Louvain  = meanvekLoudc , Laplace = meanvekLapdc, Leiden = meanvekdc)
df <- round(df, 5)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_dcsbm.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()
###
# Übersicht erstellen der einzelnen ARIs in Form einer Tabelle:
# Erstellen eines DataFrames
df <- data.frame(Louvain  = Louvaindc , Laplace = Laplacedc, Leiden = Leidendc)
df <- round(df, 5)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_dcsbm_mean.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()
########################
# Bethe-Hessian Schätzer für unsere Beispiele:
# SBM:
BH_sbm<-list()
for(t in 1:AnzahlSnap){
  stemp<-temp[[t]]
  numb<-BH(stemp, 7)
  BH_sbm[[t]]<-numb$K
}
# DC-SBM:
BH_dcsbm<-list()
for(t in 1:AnzahlSnap){
  stemp<-Adcsbm[[t]]
  numb<-BH(stemp, 7)
  BH_dcsbm[[t]]<-numb$K
  }


##################################################################################################
# Abschnitt 4.2
# Satz 3.38 als eine Funktion um alle Bedingungen zu überprufen
# matrix1 ist P matrix 2 ist Q:
pruefe_ungleichungen <- function(matrix1, matrix2) {
  if (matrix1[1, 1] >= matrix2[1, 1]) {
    stop("Fehler: Bedingung a) nicht erfüllt")
  }
  
  if (matrix2[1, 1]*matrix1[2,2] < matrix2[2, 2]*matrix1[1,1]) {
    stop("Fehler: Bedingung b) nicht erfüllt")
  }
  
  if (matrix1[2, 1] <= matrix2[2, 1]) {
    stop("Fehler: Bedingung c) nicht erfüllt")
  }
  
  return("Alle Ungleichungen sind erfüllt.")
}
P_within <- matrix(c(0.8, 0.2,  
                     0.2, 0.8), 
                   nrow = 2, ncol = 2)
P_between <- matrix(c(0.95, 0.1,  
                      0.05, 0.9), 
                    nrow = 2, ncol = 2)
pruefe_ungleichungen(P_within,P_between)
# Unser folgendes Beispiel erfüllt alle Bedingungen, wodurch W_ij ground truth positive Einträge besitzt

#######################################################
# Unser Algorithmus mit Louvain, spektral Relaxation und Leiden. Gamma wird nur grob selbst als 1 festgelegt:
# Unser Modell: 10 Realisierungen Markov hom. SBM und Wahl 0.001 als Abbruchwert und 0.05 Holdout bei falschen negativen b:
# Wir untersuchen zusätzlich die Anzahl an Schritten bis die Iteration abbricht.
# Hinweis: c ist gamma; b ist beta; und a ist alpha 
#######################################################
euklidische_norm <- function(v) {
  sqrt(sum(v^2))
}
# Alle Listen erstellen:
Genauigkeit<-list()
Schritte<-list()
Genauigkeitlou<-list()
Schrittelou<-list()
Genauigkeitlei<-list()
Schrittelei<-list()
Genauigkeitsumme<-list()
set.seed(23)
for(ä in 1:10){
  # Parameters
  n <- 100 # number of nodes
  k <- 3  # number of communities
  T <- 10 # number of time steps
  # Static community assignments (uniform generiert)
  community_assignments <- sample(1:k, n, replace = TRUE)
  # Markov transition Matrix:
  # z_i=z_j community
  P_within <- matrix(c(0.8, 0.2,  
                       0.2, 0.8), 
                     nrow = 2, ncol = 2)
  
  # Matrix Q:
  P_between <- matrix(c(0.95, 0.1,  
                        0.05, 0.9), 
                      nrow = 2, ncol = 2)
  adjacency_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:i) {
      if (i != j) {
        if (community_assignments[i] == community_assignments[j]) {
          prob <- 0.1  # p_in
        } else {
          prob <- 0.05 # p_out
        }
        adjacency_matrix[i, j] <- rbinom(1, 1, prob)
        adjacency_matrix[j, i] <- adjacency_matrix[i, j]
      }
    }
  }
  # Liste zum Speichern:
  dynamic_networks <- list()
  dynamic_networks[[1]] <- adjacency_matrix
  for (t in 2:T) {
    previous_adjacency_matrix <- dynamic_networks[[t - 1]]
    current_adjacency_matrix <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:i) {
        if (i != j) {
          if (community_assignments[i] == community_assignments[j]) {
            transition_probs <- P_within[previous_adjacency_matrix[i, j] + 1, ]
          } else {
            transition_probs <- P_between[previous_adjacency_matrix[i, j] + 1, ]
          }
          current_state <- sample(0:1, 1, prob = transition_probs)
          current_adjacency_matrix[i, j] <- current_state
          current_adjacency_matrix[j, i] <- current_state
        }
      }
    }
    
    dynamic_networks[[t]] <- current_adjacency_matrix
  }
  
  #Initialisierung a,b,c=1
  temp<-list()
  temp<-dynamic_networks
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  #persistente Kanten W_i,j
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  #surviving status:
  c<-1
  sum_surv<-c*(temp[[1]]-temp[[T]])
  Wij<-sum_pers+sum_neu+sum_surv
  
  differenz<-list()
  abc<-list()
  abc[[1]]<-c(0,0,0)
  abc[[2]]<-c(a,b,c)
  for (ü in 3:30){
    differenz[[ü]] <- abc[[ü-1]] - abc[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    if(ap>0.0001){
      #normalized spectral clustering:
      spec<-reg.SP(Wij, k, tau = 0.5, lap = TRUE,nstart=30,iter.max=100)
      community_assignments_est<-spec$cluster
      
      #Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)  # rows: previous state, columns: next state
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networks)) {
        prev_adj <- dynamic_networks[[t - 1]]
        curr_adj <- dynamic_networks[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_est[i] == community_assignments_est[j]) {
                counts_within[prev_adj[i, j] + 1, curr_adj[i, j] + 1] <- counts_within[prev_adj[i, j] + 1, curr_adj[i, j] + 1] + 1
              } else {
                counts_between[prev_adj[i, j] + 1, curr_adj[i, j] + 1] <- counts_between[prev_adj[i, j] + 1, curr_adj[i, j] + 1] + 1
              }
            }
          }
        }
      }
      
      # Schätzung von P
      P_within_est <- counts_within / rowSums(counts_within)
      
      # Schätzung von Q
      P_between_est <- counts_between / rowSums(counts_between)
      
      # Ausgabe der geschätzten Markov Dynamik:
      P_within_est
      P_between_est
      # Binden und neue Gewichte a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2))
      b_est<-log((P11*Q00)/(P00*Q11))
      c_est<-(1/(T-1))*log((P10*Q00)/(Q10*P00))
      if(-0.05<b_est & b_est<0){
        b<-0.05
      }
      Wij<-b_est*sum_pers+ a_est*sum_neu+c_est*sum_surv
      abc[[ü]]<-c(a_est, b_est, c_est)
    }
    else {
      z<-community_assignments_est 
      break
    }
  }
  Genauigkeit[[ä]]<-adjustedRandIndex(z,community_assignments)
  Schritte[[ä]]=length(abc)-1 
  ##################
  # Louvain
  # Initialisierung a,b,c=1
  temp<-list()
  temp<-dynamic_networks
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  #persistente Kanten W_i,j
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  #surviving status:
  c<-1
  sum_surv<-c*(temp[[1]]-temp[[T]])
  Wij<-sum_pers+sum_neu+sum_surv
  
  differenz<-list()
  abc<-list()
  abc[[1]]<-c(0,0,0)
  abc[[2]]<-c(a,b,c)
  # Zusätzlich Gamma selber fixieren für a;b=1 
  gamma<-1
  for (ü in 3:30){
    differenz[[ü]] <- abc[[ü-1]] - abc[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    if(ap>0.0001){
      #Louvain:
      GrWij<-graph_from_adjacency_matrix(Wij, mode="undirected")
      lou<-cluster_louvain(GrWij,NULL, gamma )
      community_assignments_est<-lou$membership
      
      #Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)  # rows: previous state, columns: next state
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networks)) {
        prev_adj <- dynamic_networks[[t - 1]]
        curr_adj <- dynamic_networks[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_est[i] == community_assignments_est[j]) {
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
      
      # Ausgabe der geschätzten Markov Dynamiken: 
      P_within_est
      P_between_est
      # Binden und neue a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2))
      b_est<-log((P11*Q00)/(P00*Q11))
      c_est<-(1/(T-1))*log((P10*Q00)/(Q10*P00))
      if(-0.05<b_est & b_est<0){
        b_est<-0.05
      }
      if(-0.1<c_est & c_est<0){
        c_est<-0.05
      }
      
      Wij<-b_est*sum_pers+ a_est*sum_neu+c_est*sum_surv
      gamma<-1
      abc[[ü]]<-c(a_est, b_est, c_est)
    }
    else {
      z<-community_assignments_est 
      
      break
    }
  }
  Genauigkeitlou[[ä]]<-adjustedRandIndex(z,community_assignments)
  Schrittelou[[ä]]=length(abc)-1 
  #################
  # Leiden wird nun getestet
  # Initialisierung a,b,c=1
  temp<-list()
  temp<-dynamic_networks
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  #persistente Kanten W_i,j
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  #surviving status:
  c<-1
  sum_surv<-c*(temp[[1]]-temp[[T]])
  Wij<-sum_pers+sum_neu+sum_surv
  
  differenz<-list()
  abc<-list()
  abc[[1]]<-c(0,0,0)
  abc[[2]]<-c(a,b,c)
  # Zusätzlich Gamma selber fixieren für a;b=1 
  gamma<-1
  for (ü in 3:30){
    differenz[[ü]] <- abc[[ü-1]] - abc[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    if(ap>0.0001){
      #Leiden:
      GrWij<-graph_from_adjacency_matrix(Wij, mode="undirected")
      lou<-cluster_leiden(GrWij, objective_function = "modularity",NULL, gamma )
      community_assignments_est<-lou$membership
      
      #Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)  # rows: previous state, columns: next state
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networks)) {
        prev_adj <- dynamic_networks[[t - 1]]
        curr_adj <- dynamic_networks[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_est[i] == community_assignments_est[j]) {
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
      
      # Ausgabe der geschätzten Markov-Dynamiken:
      P_within_est
      P_between_est
      #Binden und neue a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2))
      b_est<-log((P11*Q00)/(P00*Q11))
      c_est<-(1/(T-1))*log((P10*Q00)/(Q10*P00))
      if(-0.05<b_est & b_est<0){
        b_est<-0.05
      }
      if(-0.1<c_est & c_est<0){
        c_est<-0.05
      }
      Wij<-b_est*sum_pers+ a_est*sum_neu+c_est*sum_surv
      abc[[ü]]<-c(a_est, b_est, c_est)
    }
    else {
      z<-community_assignments_est 
      
      break
    }
  }
  Genauigkeitlei[[ä]]<-adjustedRandIndex(z,community_assignments)
  Schrittelei[[ä]]=length(abc)-1 
  ########################   
  # einfach aggregiertes Clustering     
  temp<-list()
  temp<-dynamic_networks
  sum_temp <- Reduce("+", temp)
  spec<-reg.SP(sum_temp, k, tau = 0.5, lap = TRUE,nstart=30,iter.max=100)
  community_assignments_est<-spec$cluster
  n<-adjustedRandIndex(community_assignments_est, community_assignments)
  Genauigkeitsumme[[ä]]<-n
  
}   
# Bild erstellen von den ARI Werten (Durchschnitte hier unnötig): 
Genau_Lap <- sapply(Genauigkeit, mean)
# Genau_Lap_mean<-mean(Genau_Lap)
Genau_Lou <- sapply(Genauigkeitlou, mean)
# Genau_Lou_mean<-mean(Genau_Lou)
Genau_Lei <- sapply(Genauigkeitlei, mean)
# Genau_Lei_mean<-mean(Genau_Lei)
Genau_simp<- sapply(Genauigkeitsumme, mean)
# Genau_simp_mean<- mean(Genau_simp)
df <- data.frame(Louvain  = Genau_Lou , Laplace = Genau_Lap, Leiden = Genau_Lei, Summe = Genau_simp)
df <- round(df, 5)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_temp_sbm.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()

# Anzahl der Schritte :
Schritte_Lap <- sapply(Schritte, mean)
Schritte_Lap_mean<-mean(Schritte_Lap)
Schritte_Lou <- sapply(Schrittelou, mean)
Schritte_Lou_mean<-mean(Schritte_Lou)
Schritte_Lei <- sapply(Schrittelei, mean)
Schritte_Lei_mean<-mean(Schritte_Lei)
# Tabelle mit den Schritten erstellen:
df <- data.frame(Louvain  = Schritte_Lou , Laplace = Schritte_Lap, Leiden = Schritte_Lei)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_SBM_Schritte.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()


################################################################################
# Test für extrem geringe Anzahl an Snapshots. (T=3)
# Unser Modell 10 Realisierungen Markov hom. SBM und Wahl 0.001 als Abbruchwert und 0.05 Holdout bei negativen aber nahem b:
# Funktion für Abbruchbedingung:
euklidische_norm <- function(v) {
  sqrt(sum(v^2))
}
# Alle Listen zum speichern erstellen:
Genauigkeit<-list()
Schritte<-list()
Genauigkeitlou<-list()
Schrittelou<-list()
Genauigkeitlei<-list()
Schrittelei<-list()
Genauigkeitsumme<-list()
set.seed(23)
for(ä in 1:10){
# Parameter festlegen
  n <- 100 # number of nodes
  k <- 3  # number of communities
  T <- 4 # number of time steps
  
  # Statische community assignments (uniform generiert)
  community_assignments <- sample(1:k, n, replace = TRUE)
  
  # Markov transition Matrix:
  # z_i=z_j community
  P_within <- matrix(c(0.8, 0.2,  
                       0.2, 0.8), 
                     nrow = 2, ncol = 2)
  
  # Different communities
  P_between <- matrix(c(0.95, 0.1,  
                        0.05, 0.9), 
                      nrow = 2, ncol = 2)
  # homogenes Bernoulli SBM mit p_in=0,3 ; p_out=0,1
  adjacency_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:i) {
      if (i != j) {
        if (community_assignments[i] == community_assignments[j]) {
          prob <- 0.3  # p_in
        } else {
          prob <- 0.1 # p_out
        }
        adjacency_matrix[i, j] <- rbinom(1, 1, prob)
        adjacency_matrix[j, i] <- adjacency_matrix[i, j]
      }
    }
  }
  # Initialize a list to store the adjacency matrices for each time step
  dynamic_networks <- list()
  dynamic_networks[[1]] <- adjacency_matrix
  
  for (t in 2:T) {
    previous_adjacency_matrix <- dynamic_networks[[t - 1]]
    current_adjacency_matrix <- matrix(0, n, n)
    
    for (i in 1:n) {
      for (j in 1:i) {
        if (i != j) {
          if (community_assignments[i] == community_assignments[j]) {
            transition_probs <- P_within[previous_adjacency_matrix[i, j] + 1, ]
          } else {
            transition_probs <- P_between[previous_adjacency_matrix[i, j] + 1, ]
          }
          current_state <- sample(0:1, 1, prob = transition_probs)
          current_adjacency_matrix[i, j] <- current_state
          current_adjacency_matrix[j, i] <- current_state
        }
      }
    }
    
    dynamic_networks[[t]] <- current_adjacency_matrix
  }
  #Initialisierung a,b,c=1
  temp<-list()
  temp<-dynamic_networks
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  #persistente Kanten W_i,j:
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  #surviving status:
  c<-1
  sum_surv<-c*(temp[[1]]-temp[[T]])
  Wij<-sum_pers+sum_neu+sum_surv # Wij ohne Gewichtung erstellt
  # technische Erstellung der Abbruchbedingung
  differenz<-list()
  abc<-list()
  abc[[1]]<-c(0,0,0)
  abc[[2]]<-c(a,b,c)
  for (ü in 3:30){
    differenz[[ü]] <- abc[[ü-1]] - abc[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    # Abbruchbedingung hier festgelegt
    if(ap>0.0001){
      #normalized spectral clustering:
      spec<-reg.SP(Wij, k, tau = 0.5, lap = TRUE,nstart=30,iter.max=100)
      community_assignments_est<-spec$cluster
      
      #Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)  # rows: previous state, columns: next state
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networks)) {
        prev_adj <- dynamic_networks[[t - 1]]
        curr_adj <- dynamic_networks[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_est[i] == community_assignments_est[j]) {
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
      
      # Angabe der geschätzten Markov-Dynamiken:
      P_within_est
      P_between_est
      # Binden und neue a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2))
      b_est<-log((P11*Q00)/(P00*Q11))
      c_est<-(1/(T-1))*log((P10*Q00)/(Q10*P00))
      Wij<-b_est*sum_pers+ a_est*sum_neu+c_est*sum_surv
      abc[[ü]]<-c(a_est, b_est, c_est)
    }
    else {
      z<-community_assignments_est 
      break
    }
  }
  Genauigkeit[[ä]]<-adjustedRandIndex(z,community_assignments)
  Schritte[[ä]]=length(abc)-1 
  ####################################################
  # Louvain
  # Initialisierung a,b,c=1:
  temp<-list()
  temp<-dynamic_networks
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  #persistente Kanten W_i,j
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  #surviving status:
  c<-1
  sum_surv<-c*(temp[[1]]-temp[[T]])
  Wij<-sum_pers+sum_neu+sum_surv
  
  differenz<-list()
  abc<-list()
  abc[[1]]<-c(0,0,0)
  abc[[2]]<-c(a,b,c)
  # Zusätzlich Gamma selber fixieren für a;b=1. Bei uns gamma=1. 
  gamma<-1
  for (ü in 3:30){
    differenz[[ü]] <- abc[[ü-1]] - abc[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    if(ap>0.0001){
      # Louvain:
      GrWij<-graph_from_adjacency_matrix(Wij, mode="undirected")
      lou<-cluster_louvain(GrWij,NULL, gamma )
      community_assignments_est<-lou$membership
      
      # Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)  
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networks)) {
        prev_adj <- dynamic_networks[[t - 1]]
        curr_adj <- dynamic_networks[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_est[i] == community_assignments_est[j]) {
                counts_within[prev_adj[i, j] + 1, curr_adj[i, j] + 1] <- counts_within[prev_adj[i, j] + 1, curr_adj[i, j] + 1] + 1
              } else {
                counts_between[prev_adj[i, j] + 1, curr_adj[i, j] + 1] <- counts_between[prev_adj[i, j] + 1, curr_adj[i, j] + 1] + 1
              }
            }
          }
        }
      }
      
      # Estimate transition probabilities for in-community edges
      P_within_est <- counts_within / rowSums(counts_within)
      
      # Estimate transition probabilities for between-community edges
      P_between_est <- counts_between / rowSums(counts_between)
      
      # Print the estimated Markov transition matrices
      P_within_est
      P_between_est
      #Binden und neue a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2))
      b_est<-log((P11*Q00)/(P00*Q11))
      c_est<-(1/(T-1))*log((P10*Q00)/(Q10*P00))
      Wij<-b_est*sum_pers+ a_est*sum_neu+c_est*sum_surv
      gamma<-1
      abc[[ü]]<-c(a_est, b_est, c_est)
    }
    else {
      z<-community_assignments_est 
      
      break
    }
  }
  z<-community_assignments_est
  Genauigkeitlou[[ä]]<-adjustedRandIndex(z,community_assignments)
  Schrittelou[[ä]]=length(abc)-1 
  #########################################################
  # Leiden
  # Initialisierung a,b,c=1
  temp<-list()
  temp<-dynamic_networks
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  #persistente Kanten W_i,j
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  #surviving status:
  c<-1
  sum_surv<-c*(temp[[1]]-temp[[T]])
  Wij<-sum_pers+sum_neu+sum_surv
  
  differenz<-list()
  abc<-list()
  abc[[1]]<-c(0,0,0)
  abc[[2]]<-c(a,b,c)
  # Zusätzlich Gamma selber schätzen =1 für a;b=1 
  
  gamma<-1
  for (ü in 3:30){
    differenz[[ü]] <- abc[[ü-1]] - abc[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    if(ap>0.0001){
      #Leiden:
      GrWij<-graph_from_adjacency_matrix(Wij, mode="undirected")
      lou<-cluster_leiden(GrWij, objective_function = "modularity",NULL, gamma )
      community_assignments_est<-lou$membership
      
      #Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)  # rows: previous state, columns: next state
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networks)) {
        prev_adj <- dynamic_networks[[t - 1]]
        curr_adj <- dynamic_networks[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_est[i] == community_assignments_est[j]) {
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
      
      # Angabe der geschätzten Markov-Dynamik
      P_within_est
      P_between_est
      # Binden und neue a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2))
      b_est<-log((P11*Q00)/(P00*Q11))
      c_est<-(1/(T-1))*log((P10*Q00)/(Q10*P00))
      Wij<-b_est*sum_pers+ a_est*sum_neu+c_est*sum_surv
      abc[[ü]]<-c(a_est, b_est, c_est)
    }
    else {
      z<-community_assignments_est 
      
      break
    }
  }
   z<-community_assignments_est
   Genauigkeitlei[[ä]]<-adjustedRandIndex(z,community_assignments)
  Schrittelei[[ä]]=length(abc)-1 
  ######################################################################  
  # simple aggregiertes Clustering     
  temp<-list()
  temp<-dynamic_networks
  sum_temp <- Reduce("+", temp)
  spec<-reg.SP(sum_temp, k, tau = 0.5, lap = TRUE,nstart=30,iter.max=100)
  community_assignments_est<-spec$cluster
  n<-adjustedRandIndex(community_assignments_est, community_assignments)
  Genauigkeitsumme[[ä]]<-n
  
} # nach einem durchgang wird ein neuer synthetischer Graph erstellt und ales wiederholt  


# Bilder erzeugen für die einzelnen ARI und Durchschnitte:
Genau_Lap <- sapply(Genauigkeit, mean)
Genau_Lap_mean<-mean(Genau_Lap)
Genau_Lou <- sapply(Genauigkeitlou, mean)
Genau_Lou_mean<-mean(Genau_Lou)
Genau_Lei <- sapply(Genauigkeitlei, mean)
Genau_Lei_mean<-mean(Genau_Lei)
Genau_simp<- sapply(Genauigkeitsumme, mean)
Genau_simp_mean<- mean(Genau_simp)
# Zuerst Tabelle mit den entsprechenden ARI Werten:
df <- data.frame(Louvain  = Genau_Lou , Laplace = Genau_Lap, Leiden = Genau_Lei, Summe= Genau_simp)
df <- round(df, 5)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_ari_evol_t_klein.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()

# Übersicht ARI Durchschnitte: 
df <- data.frame(Louvain  = Genau_Lou_mean , Laplace = Genau_Lap_mean, Leiden = Genau_Lei_mean, Summe= Genau_simp_mean)
df <- round(df, 5)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_ari_mean_evol_t_klein.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()

# Anzahl der Schritte in einer Tabelle darstellen:
Schritte_Lap <- sapply(Schritte, mean)
Schritte_Lap_mean<-mean(Schritte_Lap)
Schritte_Lou <- sapply(Schrittelou, mean)
Schritte_Lou_mean<-mean(Schritte_Lou)
Schritte_Lei <- sapply(Schrittelei, mean)
Schritte_Lei_mean<-mean(Schritte_Lei)
# Tabelle mit den Schritten erstellen:
df <- data.frame(Louvain  = Schritte_Lou , Laplace = Schritte_Lap, Leiden = Schritte_Lei)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_evol_tklein_Schritte.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()







#######################################################################
# Unsere Beobachtung für den Vorschlag aus 3.9:
# Von den Authoren den Algorithmus testen:
# seperat getestet:
Genauigkeitpaper<-list()
Schrittepaper<-list()
set.seed(23)
for(ä in 1:10){
  
  # Parameters
  n <- 100 # number of nodes
  k <- 3  # number of communities
  T <- 10 # number of time steps
  
  # Static community assignments
  community_assignments <- sample(1:k, n, replace = TRUE)
  
  # Markov transition Matrix:
  # z_i=z_j community
  P_within <- matrix(c(0.8, 0.2,
                       0.2, 0.8), 
                     nrow = 2, ncol = 2)
  
  # Different communities
  P_between <- matrix(c(0.95, 0.1,  
                        0.05, 0.9), 
                      nrow = 2, ncol = 2)
  # homogenes Bernoulli SBM mit p_in=0,3 ; p_out=0,1
  adjacency_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:i) {
      if (i != j) {
        if (community_assignments[i] == community_assignments[j]) {
          prob <- 0.1  # Initial probability of edge within the same community
        } else {
          prob <- 0.05 # Initial probability of edge between different communities
        }
        adjacency_matrix[i, j] <- rbinom(1, 1, prob)
        adjacency_matrix[j, i] <- adjacency_matrix[i, j]
      }
    }
  }
  # Initialize a list to store the adjacency matrices for each time step
  dynamic_networks <- list()
  dynamic_networks[[1]] <- adjacency_matrix
  
  for (t in 2:T) {
    previous_adjacency_matrix <- dynamic_networks[[t - 1]]
    current_adjacency_matrix <- matrix(0, n, n)
    
    for (i in 1:n) {
      for (j in 1:i) {
        if (i != j) {
          if (community_assignments[i] == community_assignments[j]) {
            transition_probs <- P_within[previous_adjacency_matrix[i, j] + 1, ]
          } else {
            transition_probs <- P_between[previous_adjacency_matrix[i, j] + 1, ]
          }
          current_state <- sample(0:1, 1, prob = transition_probs)
          current_adjacency_matrix[i, j] <- current_state
          current_adjacency_matrix[j, i] <- current_state
        }
      }
    }
    
    dynamic_networks[[t]] <- current_adjacency_matrix
  }
  
  #Initialisierung a,b,c=1
  temp<-list()
  temp<-dynamic_networks
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  #persistente Kanten W_i,j
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  Wij<-sum_pers+sum_neu
  
  differenz<-list()
  ab<-list()
  ab[[1]]<-c(0,0)
  ab[[2]]<-c(a,b)
  for (ü in 3:30){
    differenz[[ü]] <- ab[[ü-1]] - ab[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    if(ap>0.0001){
      #normalized spectral clustering:
      spec<-reg.SP(Wij, k, tau = 0, lap = TRUE,nstart=30,iter.max=100)
      community_assignments_est<-spec$cluster
      
      #Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)  # rows: previous state, columns: next state
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networks)) {
        prev_adj <- dynamic_networks[[t - 1]]
        curr_adj <- dynamic_networks[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_est[i] == community_assignments_est[j]) {
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
      
      # Ausgabe der geschätzten Markov-Dynamik:
      P_within_est
      P_between_est
      # Binden und neue a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01)/(Q10*Q01))
      b_est<-log((P11)/(Q11))
      if(-0.05<b_est & b_est<0){
        b<-0.05
      }
      Wij<-b_est*sum_pers+ a_est*sum_neu
      ab[[ü]]<-c(a_est, b_est)
    }
    else {
      z<-community_assignments_est 
      
      break
    }
  }
  Genauigkeitpaper[[ä]]<-adjustedRandIndex(z,community_assignments)
  Schrittepaper[[ä]]=length(ab)-1 
}
# Funktioniert obwohl negatives b (beta):
Genauigkeitpaper
Schrittepaper
# Jeder Diagonaleintrag von D größer 0 ? 
# Wir testen es aus für die letze Realisierung:
compute_degree_matrix <- function(weight_matrix) {
  # Symmetrie Test
  if (!isSymmetric(weight_matrix)) {
    stop("Weight matrix must be symmetric.")
  }
  degrees <- rowSums(weight_matrix)
  degree_matrix <- diag(degrees)
  
  return(degree_matrix)
}
D<-compute_degree_matrix(Wij)
all(D>=0) 
#Ja -> Idee: signiertes Clustering verwenden um immer clustern zu können

########################################
# Algorithmus aus Kapitel 3.9
# spektrales Clustering über signierte Laplace Matrix:
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
# Unser Modell aus 3.9 für 20 Realisierungen Markov hom. SBM und Wahl 0.001 als Abbruchwert und 0.05 Holdout bei negativen aber nahem b:
# Wir testen ein extrem ungünstiges Szenario:
P_within <- matrix(c(0.7, 0.4,  # From 0 to 0 or 1
                     0.3, 0.6), # From 1 to 0 or 1
                   nrow = 2, ncol = 2)

# Different communities
P_between <- matrix(c(0.3, 0.1,  # From 0 to 0 or 1
                      0.7, 0.9), # From 1 to 0 or 1
                    nrow = 2, ncol = 2)
T<-3
# ground truth Werte für alpha und beta für unseren Algorithmus:
P00 <- P_within_est[1, 1] 
P01<- P_within_est[1, 2]
P10<-P_within_est[2, 1]
P11<-P_within_est[2, 2]
Q00<-P_between_est[1,1]
Q01<-P_between_est[1,2]
Q10<-P_between_est[2,1]
Q11<-P_between_est[2,2]
a_gt<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2))
b_gt<-log((P11*Q00)/(P00*Q11))
c_gt<-(1/(T-1))*log((P10*Q00)/(Q10*P00))
# ground truth Werte für alpha und beta für den Algorithmus aus dem Paper:
a_gt_aut<-log((P10*P01)/(Q10*Q01))
b_gt_aut<-log((P11)/(Q11))
# Liste erzeugen:
Genauigkeit_erweiterung<-list()
Schritte_erweiterung<-list()
Parameter<- list()
set.seed(23)
for(ä in 1:10){
  # Parameters
  n <- 100 # number of nodes
  k <- 3  # number of communities
  T <- 3 # number of time steps
  
  # Static community assignments
  community_assignments <- sample(1:k, n, replace = TRUE)
  
  # Markov transition Matrix:
  # z_i=z_j community
  P_within <- matrix(c(0.7, 0.4,  # From 0 to 0 or 1
                       0.3, 0.6), # From 1 to 0 or 1
                     nrow = 2, ncol = 2)
  
  # Different communities
  P_between <- matrix(c(0.3, 0.1,  # From 0 to 0 or 1
                        0.7, 0.9), # From 1 to 0 or 1
                      nrow = 2, ncol = 2)
  # homogenes Bernoulli SBM mit p_in=0,3 ; p_out=0,1
  adjacency_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:i) {
      if (i != j) {
        if (community_assignments[i] == community_assignments[j]) {
          prob <- 0.02  # Initial probability of edge within the same community
        } else {
          prob <- 0.01  # Initial probability of edge between different communities
        }
        adjacency_matrix[i, j] <- rbinom(1, 1, prob)
        adjacency_matrix[j, i] <- adjacency_matrix[i, j]
      }
    }
  }
  # Initialize a list to store the adjacency matrices for each time step
  dynamic_networks <- list()
  dynamic_networks[[1]] <- adjacency_matrix
  
  for (t in 2:T) {
    previous_adjacency_matrix <- dynamic_networks[[t - 1]]
    current_adjacency_matrix <- matrix(0, n, n)
    
    for (i in 1:n) {
      for (j in 1:i) {
        if (i != j) {
          if (community_assignments[i] == community_assignments[j]) {
            transition_probs <- P_within[previous_adjacency_matrix[i, j] + 1, ]
          } else {
            transition_probs <- P_between[previous_adjacency_matrix[i, j] + 1, ]
          }
          current_state <- sample(0:1, 1, prob = transition_probs)
          current_adjacency_matrix[i, j] <- current_state
          current_adjacency_matrix[j, i] <- current_state
        }
      }
    }
    
    dynamic_networks[[t]] <- current_adjacency_matrix
  }
  
  #Initialisierung a,b,c=1
  temp<-list()
  temp<-dynamic_networks
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  #persistente Kanten W_i,j
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  #surviving status:
  c<-1
  sum_surv<-c*(temp[[1]]-temp[[T]])
  Wij<-sum_pers+sum_neu+sum_surv
  
  differenz<-list()
  abc<-list()
  abc[[1]]<-c(0,0,0)
  abc[[2]]<-c(a,b,c)
  for (ü in 3:30){
    differenz[[ü]] <- abc[[ü-1]] - abc[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    if(ap>0.0001){
      #normalized spectral clustering:
      spec<-clustersigniert(Wij, k, tau=0, lap=TRUE,nstart=30,iter.max=100)
      community_assignments_est<-spec$cluster
      
      #Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networks)) {
        prev_adj <- dynamic_networks[[t - 1]]
        curr_adj <- dynamic_networks[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_est[i] == community_assignments_est[j]) {
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
      
      # Angabe der geschätzten Markov-Dynamik:
      P_within_est
      P_between_est
      # Binden und neue a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2))
      b_est<-log((P11*Q00)/(P00*Q11))
      c_est<-(1/(T-1))*log((P10*Q00)/(Q10*P00))
      Wij<-b_est*sum_pers+ a_est*sum_neu+c_est*sum_surv
      abc[[ü]]<-c(a_est, b_est, c_est)
    }
    else {
      z<-community_assignments_est 
      
      break
    }
  }
  Genauigkeit_erweiterung[[ä]]<-adjustedRandIndex(z,community_assignments)
  Schritte_erweiterung[[ä]]=length(abc)-1 
  Parameter[[ä]]<- abc
  }
# Man sieht exakte Wiederherstellung und 3 bis 4 Schritte nur:
Genauigkeit_erweiterung
Schritte_erweiterung
# Die letzetn Parameter Werte filtern:
Parameter_letzte<-list()
for(ü in 1:10){
  p<-Schritte_erweiterung[[ü]]
  Parameter_letzte[[ü]]<-Parameter[[ü]][[p]]
}
# Aufteilen
alp<-list()
bet<-list()
varphi<-list()
for(ü in 1:10){
 alp[[ü]]<-Parameter_letzte[[ü]][1]
 bet[[ü]]<-Parameter_letzte[[ü]][2]
 varphi[[ü]]<-Parameter_letzte[[ü]][3]
}
alp<-sapply(alp, mean)
bet<-sapply(bet, mean)
varphi<-sapply(varphi, mean)
# Tabelle erstellen:
df <- data.frame(Alpha= alp, Beta=bet, Phi=varphi)
df<- round(df, 5)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_Parameter_Schritte.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()
# Sehr gute Ergebnisse, die vorherigen Varianten haben nicht funktioniert wegen a,b<0.
# Nur wenige Schritte notwendig !!!
##################################################################
# Unseren Algorithmus nun probieren für das evolvierende DC-SBM:
# Vergleich von unserer Erweiterung/ Algorithmus der Authoren und simples Clustering:
# Bedingungen Prüfen: 
P_within <- matrix(c(0.3, 0.3,  
                     0.7, 0.7), 
                   nrow = 2, ncol = 2)

# Different communities
P_between <- matrix(c(0.95, 0.35,  # From 0 to 0 or 1
                      0.05, 0.65), # From 1 to 0 or 1
                    nrow = 2, ncol = 2)
P00 <- P_within[1, 1] 
P01<- P_within[1, 2]
P10<-P_within[2, 1]
P11<-P_within[2, 2]
Q00<-P_between[1,1]
Q01<-P_between[1,2]
Q10<-P_between[2,1]
Q11<-P_between[2,2]
a_aut_gt_dc<-log((P10*P01)/(Q10*Q01)) # 2.48
b_aut_gt_dc<-log((P11)/(Q11))         # 0.074
# Algorithmen der Authoren ist Wohldefiniert.
a_gt_dc<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2)) # 4.790
b_gt_dc<-log((P11*Q00)/(P00*Q11))             # 1.22
T<-5 
c_gt_dc<-(1/(T-1))*log((P10*Q00)/(Q10*P00))   # 0.26
# Alle Gewichte positiv, aber:
pruefe_ungleichungen(P_within,P_between)
# c) aus 3.42 nicht erfüllt

################################################################
# Falls nicht zuvor geladen:
euklidische_norm <- function(v) {
  sqrt(sum(v^2))
}
# temporales DC-SBM: 
# 10 Realisierungen Markov hom. DC-SBM und Wahl 0.001 als Abbruchwert:
Genauigkeitdcsbm<-list()
Schrittedcsbm<-list()
Genauigkeitsummedc<-list()
Genauigkeitpaper_ver<-list()
Schrittepaper_ver<-list()
set.seed(23)
for(ä in 1:10){
  # Parameter:
  n <- 300
  K <- 17
  T <- 5
  
  community_assignmentsdc <- sample(1:K, n, replace = TRUE)
  # Theta-Werte für jeden Knoten (Beispiel: zufällige Werte)
  theta <- runif(n, 0.2, 1.8) # Theta uniform von 0.2 bis 1.8
  P_within <- matrix(c(0.3, 0.3,  
                       0.7, 0.7), 
                     nrow = 2, ncol = 2)
  
  # Different communities
  P_between <- matrix(c(0.95, 0.35,  
                        0.05, 0.65), 
                      nrow = 2, ncol = 2)
  # Initiale Adjazenzmatrix
  initial_adj <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(i!=j){
        if (community_assignmentsdc[i] == community_assignmentsdc[j]) {
          initial_adj[i, j] <- rbinom(1, 1, theta[i] * theta[j] * 0.05)
        } else {
          initial_adj[i, j] <- rbinom(1, 1, theta[i] * theta[j] * 0.02)
        }
        initial_adj[j, i] <- initial_adj[i, j]
      }
    }
  }
 
  dynamic_networksdc <- list()
  dynamic_networksdc[[1]] <- initial_adj
  
  # Simulation der dynamischen Netzwerke
  for (t in 2:T) {
    prev_adj <- dynamic_networksdc[[t - 1]]
    curr_adj <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if(i!=j){
          if (community_assignmentsdc[i] == community_assignmentsdc[j]) {
            if (prev_adj[i, j] == 1) {
              prob <- P_within[2, 2]
            } else {
              prob <- P_within[1, 2]
              prob <- prob * theta[i] * theta[j]
            }
            if(prob>1){
              prob<-1
            }
            curr_adj[i, j] <- rbinom(1, 1, prob)
          } else {
            if (prev_adj[i, j] == 1) {
              prob <- P_between[2, 2]
            } else {
              prob <- P_between[1, 2]
              prob <- prob * theta[i] * theta[j]
            }
             # Absicherung:
            if(prob>1){
              prob<-1
            }
            curr_adj[i, j] <- rbinom(1, 1, prob)
          }
          curr_adj[j, i] <- curr_adj[i, j]
        }
      }
    }
    dynamic_networksdc[[t]] <- curr_adj
  }

  # Initialisierung a,b,c=1
  temp<-list()
  temp<-dynamic_networksdc
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  # persistente Kanten W_i,j
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  # surviving status:
  c<-1
  sum_surv<-c*(temp[[1]]-temp[[T]])
  Wij<-sum_pers+sum_neu+sum_surv
  
  differenz<-list()
  abc<-list()
  abc[[1]]<-c(0,0,0)
  abc[[2]]<-c(a,b,c)
  for (ü in 3:30){
    differenz[[ü]] <- abc[[ü-1]] - abc[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    if(ap>0.0001){
      #normalized spectral clustering:
      spec<-clustersigniert(Wij, K, tau = 0, lap = TRUE,nstart=30,iter.max=100)
      community_assignments_estdc<-spec$cluster
      
      #Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)  # rows: previous state, columns: next state
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networksdc)) {
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
      
      # Ausgabe von P und Q geschätzt:
      P_within_est
      P_between_est
      #Binden und neue a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01*Q00^2)/(Q10*Q01*P00^2))
      b_est<-log((P11*Q00)/(P00*Q11))
      c_est<-(1/(T-1))*log((P10*Q00)/(Q10*P00))
      Wij<-b_est*sum_pers+ a_est*sum_neu+c_est*sum_surv
      abc[[ü]]<-c(a_est, b_est, c_est)
    }
    else {
      z<-community_assignments_estdc 
      
      break
    }
  }
  z<-community_assignments_estdc
  Genauigkeitdcsbm[[ä]]<-adjustedRandIndex(z,community_assignmentsdc)
  Schrittedcsbm[[ä]]=length(abc)-1 
  
  temp<-list()
  dynamic_networks<-dynamic_networksdc
  temp<-dynamic_networks
  a<-1
  neu<-list()
  for(t in 2:T){
    neu[[t]]<-a*(temp[[t]]-temp[[t-1]]*temp[[t]])
  }
  sum_neu <- Reduce(`+`, neu[2:T])
  
  #persistente Kanten W_i,j
  b<-1
  pers<-list()
  for(t in 2:T){
    pers[[t]]<-b*(temp[[t-1]]*temp[[t]])
  }
  sum_pers <- Reduce(`+`, pers[2:T])
  Wij<-sum_pers+sum_neu
  
  differenz<-list()
  ab<-list()
  ab[[1]]<-c(0,0)
  ab[[2]]<-c(a,b)
  for (ü in 3:30){
    differenz[[ü]] <- ab[[ü-1]] - ab[[ü - 2]]
    ap<-euklidische_norm(differenz[[ü]])
    if(ap>0.0001){
      #normalized spectral clustering:
      spec<-reg.SP(Wij, K, tau = 0, lap = TRUE,nstart=30,iter.max=100)
      community_assignments_est<-spec$cluster
      
      #Schätzer für die Markov Dynamik von den geschätzten Clustering
      counts_within <- matrix(0, 2, 2)  # rows: previous state, columns: next state
      counts_between <- matrix(0, 2, 2)
      
      for (t in 2:length(dynamic_networks)) {
        prev_adj <- dynamic_networks[[t - 1]]
        curr_adj <- dynamic_networks[[t]]
        
        for (i in 1:n) {
          for (j in 1:i) {
            if (i != j) {
              if (community_assignments_est[i] == community_assignments_est[j]) {
                counts_within[prev_adj[i, j] + 1, curr_adj[i, j] + 1] <- counts_within[prev_adj[i, j] + 1, curr_adj[i, j] + 1] + 1
              } else {
                counts_between[prev_adj[i, j] + 1, curr_adj[i, j] + 1] <- counts_between[prev_adj[i, j] + 1, curr_adj[i, j] + 1] + 1
              }
            }
          }
        }
      }
      
      # Estimate transition probabilities for in-community edges
      P_within_est <- counts_within / rowSums(counts_within)
      
      # Estimate transition probabilities for between-community edges
      P_between_est <- counts_between / rowSums(counts_between)
      
      # Print the estimated Markov transition matrices
      P_within_est
      P_between_est
      #Binden und neue a,b,c berechnen:
      P00 <- P_within_est[1, 1] 
      P01<- P_within_est[1, 2]
      P10<-P_within_est[2, 1]
      P11<-P_within_est[2, 2]
      Q00<-P_between_est[1,1]
      Q01<-P_between_est[1,2]
      Q10<-P_between_est[2,1]
      Q11<-P_between_est[2,2]
      a_est<-log((P10*P01)/(Q10*Q01))
      b_est<-log((P11)/(Q11))
      if(-0.05<b_est & b_est<0){
        b<-0.05
      }
      Wij<-b_est*sum_pers+ a_est*sum_neu
      ab[[ü]]<-c(a_est, b_est)
    }
    else {
      z<-community_assignments_est 
      break
    }
  }
  z<-community_assignments_est
  Genauigkeitpaper_ver[[ä]]<-adjustedRandIndex(z,community_assignmentsdc)
  Schrittepaper_ver[[ä]]=length(ab)-1 
# simples Clustering:
  temp<-list()
  temp<-dynamic_networksdc
  sum_temp <- Reduce("+", temp)
  spec<-reg.SP(sum_temp, K, tau = 0.5, lap = TRUE,nstart=30,iter.max=100)
  community_assignments_estdc<-spec$cluster
  q<-adjustedRandIndex(community_assignments_estdc, community_assignmentsdc)
  Genauigkeitsummedc[[ä]]<-q
 }
# Tabellen erstellen und als PDF speichern:
Genauigkeitsummedc <- sapply(Genauigkeitsummedc, mean)
Genauigkeitdcsbm <- sapply(Genauigkeitdcsbm, mean)
Genauigkeitpaper_ver <- sapply(Genauigkeitpaper_ver, mean)
# Tabelle mit den Schritten erstellen:
df <- data.frame( Vorschlag= Genauigkeitdcsbm , Paper = Genauigkeitpaper_ver, Summe = Genauigkeitsummedc)
df<-round(df,5)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_ari_dcsbm.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()
# Anzahl der Schritte:
Schrittedcsbm <- sapply(Schrittedcsbm, mean)
Schrittepaper_ver <- sapply(Schrittepaper_ver, mean)
df <- data.frame( Vorschlag= Schrittedcsbm , Paper = Schrittepaper_ver)
# Installiere das benötigte PhantomJS für webshot:
webshot::install_phantomjs()
# Erstellen der Tabelle und Speichern als Bild
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  save_kable("tabelle.html")
# Konvertieren der HTML-Datei in ein Bild
webshot("tabelle.html", "tabelle.png")
# Die Tabelle in der Konsole anzeigen:
kable(df, format = "html", table.attr = "style='width:30%;'") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# Als PDF speichern:
pdf_file <- "tabelle_schritte_dcsbm.pdf"
pdf(pdf_file, width = 8, height = 6)
grid.table(df)
dev.off()
#########################################################################