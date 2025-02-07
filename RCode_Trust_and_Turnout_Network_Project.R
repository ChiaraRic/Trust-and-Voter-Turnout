# load required packages
library(tidyverse)
library(bootnet)
library(psychonetrics)
library(qgraph)
library(psych)
library(NetworkComparisonTest)
library(mgm)
library(dplyr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path <- "C:/Users/chiar/OneDrive/Desktop/Network_Project/"

# load and prepare data
dat <- read.csv("data.csv", sep= ",")

dat_subset <- dat %>%
  select(ppltrst, 
         pplfair, 
         pplhlp, 
         trstprl, 
         trstplt, 
         trstlgl, 
         trstplc, 
         trstprt, 
         trstep, 
         trstun,
         vote)


dat_subset <- dat_subset %>%
  mutate_all(na_if, 99) %>%
  mutate_all(na_if, 88) %>%
  mutate_all(na_if, 77) %>%
  mutate_all(na_if, 66) %>%
  mutate(vote = "is.na<-" (vote, vote > 2))

# determine type for all variables
# g = continuous, c = categorical (including binary)
type <- c(rep("g", 10), "c")


# estimate network
est_Net <- estimateNetwork(dat_subset,
                default = "mgm", 
                type = type,
                tuning = 0, 
                missing = "listwise", 
                criterion = "EBIC", 
                rule = "AND")

# edge weight matrix
estEdges <- est_Net$graph
dimnames(estEdges) <- list(est_Net$labels, est_Net$labels)
estEdges
write.csv(estEdges, file = paste0(path,"WeightsMatrix.csv"))

# set maximum for plotting to maximum of estimated network
maximum_weight <- max(estEdges)

# Define variables and their groups for color coding
groups <- c(rep("social", 3), rep("political", 7), "vote") # Example for 5 nodes
group_colors <- c("political" = "#C6DDF0", "social" = "#CFF8B4", "vote" = "white")
node_colors <- group_colors[groups]


# plot network
png(filename = paste0(path, "EstimatedNetwork.png"), 
    res=300, 
    width = 13, 
    height = 12.5, 
    units = "cm")
plot(est_Net, 
     layout = "spring", 
     maximum = maximum_weight, 
     color = node_colors)  
dev.off()

# centrality plot
png(filename = paste0(path, "CentralityPlot.png"), 
    res=300, 
    width = 20, 
    height = 12.5, 
    units = "cm")
centralityPlot(est_Net, 
               scale = "raw0", 
               labels = names(dat_subset), 
               include = c("Strength",
                           "Closeness", 
                           "Betweenness"))
dev.off()

# strength only plot
png(filename = paste0(path, "StrengthPlot.png"), 
    res=300, 
    width = 20, 
    height = 12.5, 
    units = "cm")
centralityPlot(est_Net, 
               scale = "raw0", 
               labels = names(dat_subset), 
               include = "Strength")
dev.off()

# test edge weight accuracy
set.seed(2024)
boot_edgeAccuracy <- bootnet(est_Net, 
              nBoots = 1000, 
              type = "nonparametric",
              nCores = parallel::detectCores())


png(filename = paste0(path, "BootstrappingEdgeAccuracy.png"), 
    res=300, 
    width = 20, 
    height = 20, 
    units = "cm")
plot(boot_edgeAccuracy, 
     order = "sample", 
     plot = "interval", 
     split0 = TRUE
) 
dev.off()


# case-dropping bootstrap to test stability of centrality estimates
set.seed(2024)
boot_stability <- bootnet(est_Net, 
                          nBoots = 1000, 
                          default = "mgm", 
                          type = "case", 
                          nCores = parallel::detectCores(), 
                          statistics = c("Strength",
                                         "Closeness", 
                                         "Betweenness")
                          )
png(filename = paste0(path, "BootstrappingCentrality.png"), 
    res=300, 
    width = 20, 
    height = 20, 
    units = "cm")
plot(boot_stability, statistics = c("Strength",
                                  "Closeness", 
                                  "Betweenness"))
dev.off()

# calculate correlation stability (CS) coefficient
corStability(boot_stability)

# calculate small world index
smallworldIndex(qgraph(est_Net$graph))

# average shortest path length (APL), 
# which is quantified by taking the average of all shortest path lengths (geodesic distance) between each pair of nodes

### sensitivity analysis ----

# count non-voters
n_non_voters <- dat_subset %>%
  filter(vote == 2) %>%
  nrow()

# get all voter and non-voter ID's
dat_subset2 <- dat_subset
dat_subset2$id <- c(1:nrow(dat_subset))

non_voters <- dat_subset2 %>%
  filter(vote == 2) %>%
  pull(id)

voters <- dat_subset2 %>%
  filter(vote == 1) %>%
  pull(id)

# set layout
Layout <- averageLayout(est_Net$graph, est_Net$graph)

# matrix to save ID's of sampled voters 
samples_voters <- matrix(NA, nrow = 213, ncol = 10)

# sample voters 
for (samples in 1:10){
  selected_voters <- sample(x = voters, size = n_non_voters, replace = FALSE)
  samples_voters[,samples] <- selected_voters
}  

calculate_edge_weights <- function(network_all, net_sub){
  estimated_edges <- net_sub$graph
  dimnames(estimated_edges) <- list(network_all$labels, network_all$labels)
  return(estimated_edges)
}


# matrix for NCT results
NCT_tests <- matrix(NA, nrow = 2, ncol = 10)

for (i in 1:10){
sample_voters <- samples_voters[, i]  
sample_ids <- c(non_voters, sample_voters)
sub_sample <- dat_subset %>%
    slice(sample_ids)

# estimate networks
sub_net <- estimateNetwork(sub_sample,
                             default = "mgm", 
                             type = type,
                             tuning = 0, 
                             missing = "listwise", 
                             criterion = "EBIC", 
                             rule = "AND") 
edge_weights <- calculate_edge_weights(network_all = est_Net, net_sub = sub_net)
write.csv(edge_weights, file = paste0(path,"WeightsMatrix_Sub_net", i, ".csv"))

# plot networks
png(filename = paste0(path, "Sub_Net", i, ".png"), 
    res=300, 
    width = 20, 
    height = 20, 
    units = "cm")
plot(sub_net, 
     layout = Layout, 
     maximum = maximum_weight)
dev.off()

# NCT 
NCT_sensitivity <- NCT(est_Net, 
                       sub_net,
                       it=100,
                       edges = "all", 
                       p.adjust.methods = "none")
# Global Strength test
NCT_tests[1, i] <- NCT_sensitivity$glstrinv.pval
# Omnibus test
NCT_tests[2, i] <- NCT_sensitivity$nwinv.pval

# prep for next round
assign(paste0("sub_net", i), sub_net)
}

# save NCT results 
NCT_tests
dimnames(NCT_tests) <- list(
  c("Global_Strength_Test",
    "Omnibus_Test"), 
  paste0("Run", 1:10))

write.csv(NCT_tests, file = paste0(path,"NCT_test.csv"))












