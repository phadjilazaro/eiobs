# config:
# ######
# remotes::install_github("DvP17/mrio")
library(mrio)
library(tidyr)
library(readxl)
library(tidyverse)
library(ggplot2)
library(networkD3)
library(dplyr)
library(hrbrthemes)
library(webshot)
library(gridExtra)
library(fmsb)
library(RColorBrewer)
library(rgdal)
library(broom)
library(maps)
library(viridis)
library(rworldmap)
library(ggalluvial)
library(pheatmap)
library(ggplotify)
library(grid)
library(gridExtra)
# ######

setwd("C:/Users/paulh/Dropbox/Documents/These/Fond/Main/P1.SA&macroTransition/iii.ApplyingSAIO/Biophysical_Promise/Open Access/Data")


# functions: 
# ######
cascade_network <- function(graph, node, q, depth = NA, L_curr) {            
  # graph = DC.graph 
  # node = ChocSect
  # q = 0.09
  # depth=NA 
  
  if (is.na(depth)) {
    # Find the maximum possible depth (one vertex per layer)
    depth <- length(V(graph))
  }
  cn.graph <- make_empty_graph(n = 0, directed = TRUE) + vertex(node)
  V(cn.graph)[node]$layer <- 1
  V(cn.graph)[node]$input_loss <- 1
  for (l in 1:depth) {
    print(l)
    new <- 0
    # Loop over all vertices in the previous layer
    for (v in V(cn.graph)[layer == l]$name) {
      q.val <- quantile(E(graph)[.from(V(graph)[v])]$weight, 1-q)
      v.elist <- E(graph)[.from(V(graph)[v]) & weight >= q.val]
      for (e in v.elist[order(v.elist$weight, decreasing=TRUE)]) {
        v.end <- ends(graph, e)[,2]
        if (!v.end %in% V(cn.graph)$name) {
          cn.graph <- cn.graph + vertex(v.end, layer = l+1, input_loss=0)
          new <- new + 1
        }
        if (!v.end %in% V(cn.graph)[layer <= l]$name){
          cn.graph <- cn.graph + edge(V(cn.graph)[v], V(cn.graph)[v.end], weight = E(graph)[e]$weight*V(cn.graph)[v]$input_loss)
          V(cn.graph)[v.end]$input_loss<-V(cn.graph)[v.end]$input_loss+L_curr[which(exio_indus_codes$IndustryTypeCode == v.end),which(exio_indus_codes$IndustryTypeCode == v)]/diag(L_curr)[which(exio_indus_codes$IndustryTypeCode == v.end)]
        }
      }
    }
    # Break if no more nodes to assign -- that is, nothing "new"
    if (!new) break
  }
  return(cn.graph)
}

create_beau_network <- function(cn) {
  dataNT <- toVisNetworkData(cn)
  dataNT$nodes$"shape" <- "dot"
  dataNT$nodes$"style" <- "wedged"
  #dataNT$nodes$"shape" <- "diamond" # 
  #dataNT$nodes$"shape" <- "box" #if you want nodes's names in boxes
  dataNT$nodes$"physics" <- FALSE
  dataNT$nodes$font.size <- 40
  dataNT$nodes$size <- 22
  dataNT$nodes$"color" <- "black"
  
  dataNT$edges$"color" <- "black"
  dataNT$edges$font.size <- 30
  # dataNT$edges$arrow.size <- 10
  # dataNT$edges <- dataNT$edges[,c(1,2,9,3:8)]
  # dataNT$edges$label <- dataNT$edges$weight
  dataNT$edges$value <- dataNT$edges$weight
  dataNT$edges <- dataNT$edges %>% rename(label=weight)
  dataNT$edges$label <- as.character(round(dataNT$edges$label*100,5))
  
  
  return(
    visNetwork(nodes = dataNT$nodes, edges = dataNT$edges, main = paste("Main backward output connections (", ChocSect,")")) %>%
      visEdges(arrows = "to") %>% visLayout(hierarchical = FALSE) %>% visPhysics(enabled = TRUE)
  )
}
# ######


# # I/ Import and format raw data
# #################################
# #################################

# # I.1/ BACH (Balance sheet) data (ds.bach_formatted)
# ################
ds_nace_names <- read.csv("NACE_REV2_20201021_190552.csv", sep=";", stringsAsFactors = F) %>% 
  select(Code, ADDEDnomination, Parent, Description)

ds.bach <- read.csv("20230630.csv",stringsAsFactors = F,sep=";", skip=1)
# remove false sectors
ds.bach <- ds.bach %>% filter(!sector %in% c("Zc", "Z0", "Mc"))

bach_regions <- unique(ds.bach$country)

ds.bach_formatted <- ds.bach %>% 
  filter(size=="0", sample==0, nchar(ds.bach$sector)>1) %>% 
  select(country,year,sector,total_assets,turnover,gross_value_added,nb_firms,employees, 
         eval(grep("wm",colnames(ds.bach),value=T))) %>%
  mutate_at(vars((starts_with("A") & ends_with("wm")) | (starts_with("L") & ends_with("wm")) | (starts_with("E") & ends_with("wm"))),
            .funs = list(~ total_assets/100 * .x))  %>%
  mutate_at(vars((starts_with("I") & ends_with("wm")) | (starts_with("N") & ends_with("wm"))),
            .funs = list(~ total_assets/100 * .x)) 

base_bach_cols <- c("country", "year", "sector", "turnover", "total_assets", "gross_value_added", "nb_firms", "employees")
bach_liability_cols <- c("Lp_wm", "L1_wm", "L2_wm", "L31_wm", "L32_wm", "L4_wm", "L5_wm", "L6_wm")
bach_equity_cols <- c("E_wm")

# ################

# # I.2/ Exiobase (EIO) data (F_sc_regionEindus ; L)
# ################
load("exiobase_regions.Rdata")
load("exiobase_industries.Rdata")
load("exiobase_region_groups.Rdata")
year <- 2019
type = "ixi"
exio_indus_codes <- data.frame(read_excel("types_version2.2.2.xlsx", sheet=4))
n_industries <- 163
exio_region_codes_full <- data.frame(read_excel("Classification_jiec12715-sup-0009-suppmat-9.xlsx", sheet=5))
exio_region_codes <- exio_region_codes_full[-1,]
rm(exio_region_codes_full)
regions <- exio_region_codes$DESIRE.code
n_region <- length(regions)

path <- c(
  paste0("IOT_", year, "_", type, "/A.txt"),
  paste0("IOT_", year, "_", type, "/Y.txt"),
  paste0("IOT_", year, "_", type, "/satellite/F.txt"),
  paste0("IOT_", year, "_", type, "/satellite/F_hh.txt")
)

n <- ifelse(type == "ixi", 7989, 9802)

# read matrices
A <- as.matrix(data.table::fread(path[1], select = 3:n, skip = 3,
                                 header = F))
FD <- as.matrix(data.table::fread(path[2], select = 3:345, skip = 3,
                                  header = F))
Q <- data.matrix(data.table::fread(path[3], select = 2:(n-1), skip = 2,
                                   header = F), rownames.force = NA)
# Remove rows that have been added in the stressor "now-casted" database relative to the original EXIOBASE ("Energy Carrier": rows 1105 to 1113)
Q <- Q[-(1),][-(1105:1113),]

# ... leontief invese and total output
I <- diag(ncol(A))
L <- solve(I - A)
X <- L %*% FD
# Z <- A %*% diag(rowSums(X))

xout <- as.matrix(rowSums(X))
totalinput <- t(xout)

E <- t(mrio::cf_exio$cf_cc) %*% Q

# CABERNARD ALGO links
E_per_x <- E / totalinput
E_per_x[which(is.nan(E_per_x))] <- 0 # remove NaNs
E_per_x[which(is.infinite(E_per_x))] <- 0 # remove Infinites
E_per_x[which(E_per_x < 0)] <- 0 # remove Negatives


diag_Epx <- diag(as.vector(E_per_x[1,]))
diag_FD  <- diag(as.vector(rowSums(FD)))
scEmat   <- diag_Epx %*% L %*% diag_FD

exio_region.indus <- rep(regions, each=n_industries)
exio_region.indus <- paste0(exio_region.indus, rep(exio_indus_codes$IndustryTypeCode, n_region))
F_sc_regionEindus_all <- as.data.frame(scEmat)
colnames(F_sc_regionEindus_all) <- exio_region.indus

F_sc_regionEindus <- F_sc_regionEindus_all[,which(substr(colnames(F_sc_regionEindus_all), 1, 2) %in% bach_regions)]
F_sc_regionEindus <- t(F_sc_regionEindus)

# ################


# # I.3/ NACE - Exiobase Concordance table (lst_concordance_ExioNace)
# ################ 

ds <- data.frame(read_excel("Correspondance_exiobase2_to_NACE2008_v2 - Changed.xlsx"))
#### COMMENT.0: the concordance table use Exiobase.2 category names, but it is the same sectors than in Exiobase.3. 
ds$NACE.Code.A <- substr(ds$NACE.Code.A, 1, 3)
ds <- ds[,1:4]
#### COMMENT.1: I40.2.a and I40.2.b are differentiated in the concordance table, but not in the Exiobase data. So, I merge them as "i40.2"
ds[which(ds$IndustryTypeCode=="I40.2.a"), ] <- c("i40.2", "Manufacture and Distribution of gas", "D35", "Manufacture and Distribution of gas")
ds[which(ds$IndustryTypeCode=="I40.2.b"), ] <- c("i40.2", "Manufacture and Distribution of gas", "D35", "Manufacture and Distribution of gas")
#### COMMENT.2: BACH data includes sectors "M701" and "M702" when the concordance table only includes "M70". So, I add the two categories "M701" and "M702" in the correspondance table and allocate them to the same exio-sector than "M70" (then, their respectve turnover will serve to weight the allocation of env. footprint)
ds <- rbind(ds, ds[which(ds$NACE.Code.A=="M70"),])
ds$NACE.Code.A[nrow(ds)] <- "M701"
ds <- rbind(ds, ds[which(ds$NACE.Code.A=="M70"),])
ds$NACE.Code.A[nrow(ds)] <- "M702"
#### COMMENT.2bis: same issue with "K642" -> "K64"
ds <- rbind(ds, ds[which(ds$NACE.Code.A=="K64"),])
ds$NACE.Code.A[nrow(ds)] <- "K642"

ds_concordance_ExioNace <- unique(ds)

# Buid the EXIO to NACE list
##
lst_concordance_ExioNace <- vector(mode = "list", length = dim(exio_indus_codes)[1])
names(lst_concordance_ExioNace) <- exio_indus_codes$IndustryTypeCode
for (i in 1:length(lst_concordance_ExioNace)){
  for (j in 1:dim(ds_concordance_ExioNace)[1]){
    if (names(lst_concordance_ExioNace)[i] == ds_concordance_ExioNace$IndustryTypeCode[j] & !(ds_concordance_ExioNace$NACE.Code.A[j] %in% lst_concordance_ExioNace[[i]])){
      lst_concordance_ExioNace[[i]] <- append(lst_concordance_ExioNace[[i]], ds_concordance_ExioNace$NACE.Code.A[j])
    }
  }
}

# Buid the NACE to EXIO list
###
nace_indus_code <- unique(ds_concordance_ExioNace$NACE.Code.A)
lst_concordance_NaceExio <- vector(mode = "list", length = length(nace_indus_code))
names(lst_concordance_NaceExio) <- nace_indus_code
for (i in 1:length(lst_concordance_NaceExio)){
  for (j in 1:dim(ds_concordance_ExioNace)[1]){
    if (names(lst_concordance_NaceExio)[i] == ds_concordance_ExioNace$NACE.Code.A[j] & !(ds_concordance_ExioNace$IndustryTypeCode[j] %in% lst_concordance_NaceExio[[i]])){
      lst_concordance_NaceExio[[i]] <- append(lst_concordance_NaceExio[[i]], ds_concordance_ExioNace$IndustryTypeCode[j])
    }
  }
}
# ################ 


# # II. Plot results 
# #################################
# #################################


# II.1. Responsibility case study (V.B.1)
# ################

# II.1.A Prepare dataset with exio stressors by bach sectors
# #######

year <- 2019
yr=year
bach_countries <- unique(ds.bach_formatted$country)

# Adjust the targeted regions to match BACH and EXIOBASE
lst_set_region <- vector(mode = "list", length = length(bach_countries))
names(lst_set_region) <- bach_countries
for (i in 1:length(lst_set_region)){
  lst_set_region[[i]] <- c((1:163)+163*(i-1))
}
lst_bach_F_regionEindus <- vector(mode = "list", length = length(bach_countries))
names(lst_bach_F_regionEindus) <- bach_countries
lst_bach_F_regionEindus_sub <- vector(mode = "list", length = length(c("FR", "DE", "ES", "IT", "PL")))
names(lst_bach_F_regionEindus) <- c("FR", "DE", "ES", "IT", "PL")


# Build the dataset with exio stressors by bach sectors
ctry = "FR"

# Select turnover of the country-year
sec.turnover <- ds.bach_formatted %>% select(country, year, sector, turnover) %>% filter(country %in% ctry ,year==yr)
sec.turnover$turnover[which(is.na(sec.turnover$turnover))] <- 0 # remove NaNs

ctryNum <- which(ctry == bach_countries)

# Get the stressor matrix by exio sectors
exio_F_sc_regionEindus <- F_sc_regionEindus[lst_set_region[[ctryNum]], ]

# Set up dataframe
ds_bach_F_i <- data.frame(row.names = c("sector", exio_region.indus), stringsAsFactors = F)
for (i in 1:dim(sec.turnover)[1]){
  ds_bach_F_i <- cbind(ds_bach_F_i, 0, stringsAsFactors = FALSE)
  ds_bach_F_i[1,i] <- sec.turnover$sector[i]
  colnames(ds_bach_F_i)[i] <- sec.turnover$country[i]
}
# Fill dataframe with Exio stressors by bach sectors
for (i in 1:dim(ds_bach_F_i)[2]){
  for (j in 1:length(names(lst_concordance_ExioNace))){
    if (as.character(ds_bach_F_i[1,i]) %in% lst_concordance_ExioNace[[j]] & length(lst_concordance_ExioNace[[j]])==1) {
      ds_bach_F_i[2:nrow(ds_bach_F_i),i] <- as.numeric(ds_bach_F_i[2:nrow(ds_bach_F_i),i]) + 
        as.numeric(exio_F_sc_regionEindus[j,1:ncol(exio_F_sc_regionEindus)])
    }
    if (as.character(ds_bach_F_i[1,i]) %in% lst_concordance_ExioNace[[j]] & length(lst_concordance_ExioNace[[j]])>1) {
      Yshare <- sec.turnover$turnover[i]/sum(sec.turnover$turnover[which(sec.turnover$sector %in% lst_concordance_ExioNace[[j]])], na.rm = T)
      ds_bach_F_i[2:nrow(ds_bach_F_i),i] <- as.numeric(ds_bach_F_i[2:nrow(ds_bach_F_i),i]) + 
        as.numeric(exio_F_sc_regionEindus[j,1:ncol(exio_F_sc_regionEindus)])*Yshare
    }
  }
}
ds_bach_F_i <- as.data.frame(t(ds_bach_F_i))
ds_bach_F_i[,-1] <- as.data.frame(lapply(ds_bach_F_i[,-1], function(x) as.numeric(as.character(x))))
ds_bach_F_i$country <- ctry
print(ctry)

ds_bach_F_i_sub <- ds_bach_F_i
ds_bach_F_i_sub$sector <- substr(ds_bach_F_i_sub$sector, 1, 1)
ds_bach_F_i_sub <- ds_bach_F_i_sub %>% group_by(sector, country) %>% summarise_all(list(sum))
colnames(ds_bach_F_i_sub)[3:ncol(ds_bach_F_i_sub)] <- 
  paste0(exiobase_regions$ID_region_group, "_", rep(exiobase_industries$ID_industry_group, n_region))
ds_bach_F_i_sub_num <- ds_bach_F_i_sub[,-c(1,2)]
ds_bach_F_i_sub_num <- t(ds_bach_F_i_sub_num)
ds_bach_F_i_sub_num <- as.data.frame(ds_bach_F_i_sub_num)
ds_bach_F_i_sub_num$sect_coun <- sub("\\..*", "", rownames(ds_bach_F_i_sub_num)) # before the point
ds_bach_F_i_sub_num <- ds_bach_F_i_sub_num  %>% group_by(sect_coun) %>% summarise_all(list(sum))
sect_coun <- ds_bach_F_i_sub_num$sect_coun
ds_bach_F_i_sub_num <- as.data.frame(t(ds_bach_F_i_sub_num[,-1]))
colnames(ds_bach_F_i_sub_num) <- sect_coun
ds_bach_F_i_sub <- cbind(ds_bach_F_i_sub[,c(1,2)],ds_bach_F_i_sub_num)
lst_bach_F_regionEindus_sub[[ctry]] <- ds_bach_F_i_sub

# Embedded footprint in liabilities 
ds.bach_formatted_sub <- ds.bach_formatted[,1:49]
ds.bach_formatted_sub$sector <- substr(ds.bach_formatted_sub$sector,1,1)
ds.bach_formatted_sub <- ds.bach_formatted_sub %>% group_by(sector, country, year) %>% summarise_all(list(sum))
df.liabs_fp <- ds.bach_formatted_sub %>%
  select(all_of(base_bach_cols), all_of(bach_liability_cols) , all_of(bach_equity_cols)) %>%
  filter(year==yr) %>%
  filter(country =="FR") %>% # !!
  gather(key = "variable", value = "value", -c(country, year, sector, total_assets)) %>%
  left_join(lst_bach_F_regionEindus_sub[["FR"]], by = c("sector", "country")) 

ds.embedded_fp <-  df.liabs_fp
ds.embedded_fp[,7:ncol(ds.embedded_fp)] <- ds.embedded_fp[,7:ncol(ds.embedded_fp)] * ds.embedded_fp$value / ds.embedded_fp$total_assets
ds.embedded_fp <- ds.embedded_fp %>% 
  select(-value, -total_assets) %>% 
  rename(instrument = variable) %>% 
  filter(instrument %in% c("L1_wm", "L2_wm", "L31_wm", "L32_wm", "Lp_wm", "L4_wm", "L5_wm", "L6_wm", "L6_wm", "E_wm"))
ds.embedded_fp <- ds.embedded_fp %>% gather(key = "region.indus", value = "value", -c(country, year, sector, instrument))
ds.embedded_fp$region.indus <- sub("\\..*", "", ds.embedded_fp$region.indus)
ds.embedded_fp$region <- substr(sub("_.*", "", ds.embedded_fp$region.indus), 2,4)
ds.embedded_fp$indus <- sub('.+_(.+)', '\\1', ds.embedded_fp$region.indus)
exiobase_region_groups$ID_region_group <- as.character(exiobase_region_groups$ID_region_group) 
ds.embedded_fp <- ds.embedded_fp %>% left_join(exiobase_region_groups %>% 
                                                 rename(region = ID_region_group), by = c("region")) 
ds.embedded_fp <- ds.embedded_fp %>% select(-region) %>% rename(region = Exiobase_region_group)
ds.embedded_fp <- ds.embedded_fp %>% select(-region.indus) %>% group_by(sector, country, year, instrument, indus, region) %>% summarise_all(funs(sum), na.rm = TRUE)

# #######


# II.1.B. Plot the sankey 
# #######

ds.embedded_fp_financial <- ds.embedded_fp %>% 
  filter(!is.na(value)) %>% filter(!value==0)
ds.embedded_fp_financial$country <- paste0(ds.embedded_fp_financial$country, "_F")
ds.embedded_fp_financial$footprint <- "GHG"
ds.embedded_fp_financial$owner <- "Bank"
ds.embedded_fp_financial$owner[which(ds.embedded_fp_financial$instrument %in% c("L32_wm", "Lp_wm", "L4_wm", "L5_wm", "L6_wm"))] <- "Other non-financial"
ds.embedded_fp_financial$owner[which(ds.embedded_fp_financial$instrument %in% c("L1_wm", "L31_wm"))] <- "Other financial"
ds.embedded_fp_financial$owner[which(ds.embedded_fp_financial$instrument %in% c("E_wm"))] <- "Equity shareholder"

ds.embedded_fp_financial <- ds.embedded_fp_financial %>% ungroup() %>% 
  mutate(country = "France") %>% 
  mutate(instrument = sub("_wm.*", "", instrument))
my_palette <- brewer.pal(name="Set3",n=12)[-1]
ggplot(ds.embedded_fp_financial %>% filter(country == "France") %>% 
         select(-year) %>% 
         filter(!is.na(value)),
       aes(y = value,
           axis1 = footprint, axis2 = region, axis3 = indus, axis4 = country, axis5 = sector, axis6 = instrument, axis7 = owner,
           fill = instrument)) +
  geom_flow() +
  geom_alluvium(aes(fill=instrument)) +
  scale_x_discrete(limits = c("Footprint", "Region", "Indus", "Country", "Funded Sector", "Liability Instrument", "Asset Owner")) +
  geom_stratum(alpha = .2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
  xlab("GHG") + 
  theme_ipsum() + # theme_minimal() +
  scale_fill_manual(values = my_palette) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# #######



# II.2. Exposure case study (V.B.2)
# ################

# II.1.A. Prepare datasets with the shock effects and the liabilities issued
# #######

# II.1.A.a. The Economic effects dataset (for the network plot)
# ####

# Shock on Exio sectors 
L_FR <- matrix(0, nrow = dim(L)[1], ncol = dim(L)[2])
nstart_FR = 1631
nend_FR = nstart_FR + 162 

L_FR[nstart_FR:nend_FR,nstart_FR:nend_FR] <- L[nstart_FR:nend_FR,nstart_FR:nend_FR]
I <- diag(x = 1, nrow=dim(L)[1], ncol = dim(L)[1])
L_FR_ind <- L_FR - I

Chosen_export_shock <- 1
ChocSect <- "i23.2"
chocSect_num <- which(exio_indus_codes$IndustryTypeCode==ChocSect) + (nstart_FR-1)
n = 163 * 49 
oneVecChoc <- c(rep(0, n)) 
oneVecChoc[chocSect_num] <- 1

ChocSect_sect_stranding.Y <- Chosen_export_shock*L_FR_ind %*% matrix(oneVecChoc, nrow = n)
effect_FR <- ChocSect_sect_stranding.Y[nstart_FR:nend_FR]
df_effect_FR <- exio_indus_codes %>% select(IndustryTypeCode, IndustryTypeName) %>% 
  mutate(effect = effect_FR)


# Transform Exio effects into NACE sectors
df_effect_FR$nace_sector <- NA
for (exio_i in df_effect_FR$IndustryTypeCode) {
  for (nace_i in names(lst_concordance_NaceExio)) {
    if (exio_i %in% lst_concordance_NaceExio[[nace_i]]) {
      df_effect_FR$nace_sector[which(df_effect_FR$IndustryTypeCode == exio_i)] <- nace_i
    }
  }
}
df_NACE_names <- ds_nace_names %>% ungroup() %>% filter(ADDEDnomination != "") %>% distinct() %>%
  select(ADDEDnomination, Description) %>%
  rename(nace_sector = ADDEDnomination,
         nace_name = Description)

df_effect_FR_NACE <-  df_NACE_names %>% 
  left_join(df_effect_FR %>% 
              group_by(nace_sector) %>% 
              summarise_at(vars(effect), funs(sum(.))) %>% 
              ungroup()) %>% 
  mutate(effect = ifelse(is.na(effect), 0, effect))

# ####

# II.1.A.b. The Liability dataset (for the pie charts)
# ####

df.liabs_fp <- ds.bach_formatted %>% 
  select(all_of(base_bach_cols), all_of(bach_liability_cols) , all_of(bach_equity_cols)) %>% 
  filter(year==2018) %>%
  rename(L2_Credit_Institutions = L2_wm,
         L4_Trade_payables = L4_wm,
         Lp_Provisions = Lp_wm,
         L31_Other_financial_creditors = L31_wm,
         L32_Other_non_financial_creditors = L32_wm, 
         E_Equities = E_wm,
         L1_Bonds_and_obligations = L1_wm,
         L5_Payments_received_on_account_for_orders = L5_wm,
         L6_Deferred_liabilities = L6_wm)
# ####

# II.1.B. Plot the network and the pie charts (both are then visually combined in PowerPoint) 
# #######

# Network plot 
varContent <- L_FR[nstart_FR:nend_FR,nstart_FR:nend_FR]
rownames(varContent) <- exio_indus_codes$IndustryTypeCode
library(igraph)
DC.graph <- simplify(graph_from_adjacency_matrix(t(varContent), mode = "directed", weighted = TRUE))
cn <- cascade_network(DC.graph, ChocSect, q = 0.01380, depth=NA, varContent)
create_beau_network(cn)


# Financial pie charts 
df_liabs_FR <- df.liabs_fp %>% filter(country=="FR") %>% 
  select(-employees, -nb_firms, -gross_value_added, -turnover) %>% 
  gather(key = liab_type, value = value, -c("country", "year", "sector", "total_assets"))

transport_block <- c("H52", "H50", "H49", "G45")
manuf_block <- c("C30","C26","C32","C25","C24")
trade_block <- c("G46", "G47")
finance_block <- c("M69", "M70", "M71", "M73", "M74", "M75", "N78", "N81", "N82", "N80", "Q87", "M701", "M702",
                   "H53", "J59", "J60", "J61", 
                   "K64", "K642", "L68")
refinery_sector <- c("C19")

write.csv2(data.frame(refinery_sector, transport_block, trade_block, finance_block, manuf_block))

my_palette <- brewer.pal(name="Paired",n=12) 

block <- "transport_block"
#  "transport_block", "manuf_block", "trade_block", "finance_block", "refinery_sector"
tot_assets <- format(sum(df_liabs_FR %>% filter(sector %in% get(block)) %>% distinct(total_assets), na.rm = T), big.mark   = " ")
df_liabs_FR %>% filter(sector %in% get(block)) %>% 
  ggplot(., aes(x="", y=value, fill=liab_type))+
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_palette) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(axis.text.x=element_blank(),
        axis.title = element_blank()) +
  annotate("text", x=-4, y=10, colour = "white",
           label= paste0(tot_assets," \u20ac"), size = 7) 

# ####



