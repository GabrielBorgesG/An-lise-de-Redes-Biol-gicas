install.packages("data.table")
install.packages("igraph")
install.packages("dplyr")
install.packages("magrittr")
install.packages("BiocManager")
install.packages("tibble")
library(data.table)
library(igraph)
library(dplyr)
library(magrittr)
library(BiocManager)
library(tibble)

getwd()

# IL11
il11_arq <- read.delim("/home/bioinformatica/Gabriel_Borges/Rstudio/IL11.tsv", stringsAsFactor = F)
il11_graph <- graph_from_data_frame(il11_arq, directed = FALSE, vertices = NULL)
il11_graph
plot(il11_graph)

# Inflammation
inflam_arq <- read.delim("/home/bioinformatica/Gabriel_Borges/Rstudio/Inflam.tsv", stringsAsFactor = F)
inflam_graph <- graph_from_data_frame(inflam_arq, directed = FALSE, vertices = NULL)
inflam_graph
plot(inflam_graph)

# Autophagy
autop_arq <- read.delim("/home/bioinformatica/Gabriel_Borges/Rstudio/Autop.tsv", stringsAsFactor = F)
autop_graph <- graph_from_data_frame(autop_arq, directed = FALSE, vertices = NULL)
autop_graph
plot(autop_graph)

# Unindo as tabelas
total_arq <- bind_rows(il11_arq, inflam_arq, autop_arq)
total_graph <- graph_from_data_frame(total_arq, directed = FALSE, vertices = NULL)
total_graph
plot(total_graph)

# Centralidades

## Centralidade Grau
df_grau <- as.data.frame(degree(total_graph))
View(df_grau)
setDT(df_grau, keep.rownames = "Nodo")
colnames(df_grau)[2] ="Grau"
Grau_médio = mean(df_grau$Grau)
Hubs = df_grau %>% filter(df_grau$Grau > Grau_médio)
Não_hubs = df_grau %>% filter(df_grau$Grau < Grau_médio)

## Centralidade de intermediação
df_interm <- as.data.frame(betweenness(total_graph))
setDT(df_interm, keep.rownames = "Nodos")
colnames(df_interm)[2] ="Intermedialidade_<bet>"
Intermedialidade_média = mean(df_interm$"Intermedialidade_<bet>")

Gargalo = df_interm %>% filter(df_interm$"Intermedialidade_<bet>" > Bet_médio)
Não_gargalo = df_interm %>% filter(df_interm$"Intermedialidade_<bet>" < Bet_médio)

## Centralidade por Proximidade
df_prox <- closeness(total_graph, mode = "all", weights = E(total_graph)$weight)
# colnames(df_prox)[1] = "Proximidade" # Não funcionou
df_prox
Proximidade_média <- mean(df_prox)

## Centralidade Eigenvector
df_autovetor <- eigen_centrality(total_graph, directed = FALSE, scale = TRUE, weights = NULL)
df_autovetor
autovetor_médio <- mean(df_autovetor[["vector"]])

# Gráficos de dispersão (betweenness x node degree; betweenness x eigenvector)

  # Intermedialidade x Grau
    df_grau_interm = cbind(df_grau, df_interm)
      # Apagando a coluna repetida Nodes
      colnames(df_grau_interm)[3] ="Apagar"
      df_grau_interm$"Apagar"<-NULL
    plot(df_grau_interm$Grau, df_grau_interm$"Intermedialidade_<bet>", xlab="Grau", ylab="Intermedialidade")
    abline(h = Intermedialidade_média) # h (horizontal)
    abline(v = Grau_médio) # v (vertical)
    plot(df_grau_interm$"Intermedialidade_<bet>", df_grau_interm$Grau, xlab="Intermedialidade", ylab="Grau")
    abline(v = Intermedialidade_média)
    abline(h = Grau_médio)
    
  # Intermedialidade x Eigenvector
    df_autovetor_interm = cbind(df_interm, df_autovetor[["vector"]])
    colnames(df_autovetor_interm)[3] ="Autovetor"
    plot(df_autovetor_interm$Autovetor, df_autovetor_interm$"Intermedialidade_<bet>", xlab="Autovetor", ylab="Intermedialidade")
    abline(h = Intermedialidade_média)
    abline(v = autovetor_médio)
    plot(df_autovetor_interm$"Intermedialidade_<bet>", df_autovetor_interm$Autovetor, xlab="Intermedialidade", ylab="Autovetor")
    abline(v = Intermedialidade_média)
    abline(h = autovetor_médio)
    
  # Intermedialidade x Proximidade
    df_prox_interm = cbind(df_interm, df_prox)
    colnames(df_prox_interm)[3] ="Proximidade"
    plot(df_prox_interm$"Intermedialidade_<bet>", df_prox_interm$Proximidade, xlab="Intermedialidade", ylab="Proximidade")
    abline(h = Proximidade_média)
    abline(v = Intermedialidade_média)
    
  # Grau x Autovetor
    df_grau_autovetor = cbind(df_grau, df_autovetor[["vector"]])
    colnames(df_grau_autovetor)[3] ="Autovetor"
    plot(df_grau_autovetor$Autovetor, df_grau_autovetor$Grau, xlab="Autovetor", ylab="Grau")
    abline(h = Grau_médio)
    abline(v = autovetor_médio)
    plot(df_grau_autovetor$Grau, df_grau_autovetor$Autovetor, xlab="Grau", ylab="Autovetor")
    abline(v = Grau_médio)
    abline(h = autovetor_médio)
    
  # Grau x Proximidade
    df_grau_prox = cbind(df_grau, df_prox)
    colnames(df_grau_prox)[3] ="Proximidade"
    plot(df_grau_prox$Grau, df_grau_prox$Proximidade, xlab="Grau", ylab="Proximidade")
    abline(h = Proximidade_média)
    abline(v = Grau_médio)
    plot(df_grau_prox$Proximidade, df_grau_prox$Grau, xlab="Proximidade", ylab="Grau")
    abline(v = Proximidade_média)
    abline(h = Grau_médio)
    
  # Proximidade x Autovetor
    df_prox_autovetor = cbind(df_prox, df_autovetor[["vector"]])
    colnames(df_prox_autovetor)[1] ="Proximidade"
    colnames(df_prox_autovetor)[2] ="Autovetor"
    df_prox_autovetor <- as.data.frame(df_prox_autovetor)
    df_prox_autovetor <- df_prox_autovetor %>% rownames_to_column(var = "rownames")
    colnames(df_prox_autovetor)[1] ="Nodos"
    plot(df_prox_autovetor$Autovetor, df_prox_autovetor$Proximidade, xlab="Autovetor", ylab="Proximidade")
    abline(h = Proximidade_média)
    abline(v = autovetor_médio)
    plot(df_prox_autovetor$Proximidade, df_prox_autovetor$Autovetor, xlab="Proximidade", ylab="Autovetor")
    abline(v = Proximidade_média)
    abline(h = autovetor_médio)
    
## Comunidades
total_graph_prox <- walktrap.community(total_graph, weights = E(total_graph)$weight, steps = 4, merges = TRUE, modularity = FALSE)
total_graph_prox
print(membership(total_graph_prox))
plot(total_graph_prox, total_graph)

par(oma=c(3,3,3,3)) # all sides have 3 lines of space
par(mar=c(5,4,4,2) + 0.1) # mar=c(b,l,t,r)


    # Visualizar cada comunidade
    for(i in 1:15){
      print(total_graph_prox[i])
    }
    
    # Nomes dos nodos
    nomes_nodos <- as.data.frame(total_graph_prox$names)
    # Número de nodos em cada comunidade
    comunidade_do_nodo <- as.data.frame(total_graph_prox$membership)
    # Juntando os dois
    Comunidades <- cbind(nomes_nodos, comunidade_do_nodo)
    colnames(Comunidades)[1] = "Nodo"
    colnames(Comunidades)[2] = "ID_comunidade"
    
    # Número de membros dentro de cada comunidade
    numero_de_membros <- sizes(total_graph_prox)
    numero_de_membros
    plot(numero_de_membros)
    
# IL11 é parte de uma comunidade grande?
    
  # Separando as comunidades:
    
    Cluster_1 <- subset(Comunidades, ID_comunidade == "1")
    Cluster_2 <- subset(Comunidades, ID_comunidade == "2")
    Cluster_3 <- subset(Comunidades, ID_comunidade == "3")
    Cluster_4 <- subset(Comunidades, ID_comunidade == "4")
    Cluster_5 <- subset(Comunidades, ID_comunidade == "5")
    Cluster_6 <- subset(Comunidades, ID_comunidade == "6")
    Cluster_7 <- subset(Comunidades, ID_comunidade == "7")
    Cluster_8 <- subset(Comunidades, ID_comunidade == "8")
    Cluster_9 <- subset(Comunidades, ID_comunidade == "9")
    Cluster_10 <- subset(Comunidades, ID_comunidade == "10")
    Cluster_11 <- subset(Comunidades, ID_comunidade == "11")
    Cluster_12 <- subset(Comunidades, ID_comunidade == "12")
    Cluster_13 <- subset(Comunidades, ID_comunidade == "13")
    Cluster_14 <- subset(Comunidades, ID_comunidade == "14")
    Cluster_15 <- subset(Comunidades, ID_comunidade == "15")
    
  # Sim, da comunidade 4, com 66 nodos
    
# ClusterProfiler para avaliar os processos biológicos associados com as comunidades
# das quais a IL11 faz parte

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Hs.eg.db", dependencies = TRUE)
options(BioC_mirror = "https://bioconductor.org")
library(clusterProfiler)
library(org.Hs.eg.db)

# Enriquecimento por ontologia gênica
for(i in 1:15){
  PB <- enrichGO(gene = Cluster_1, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                 pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.2)
  head(PB)
  barplot(PB, showCategory=10)
  dotplot(PB, showCategory=10)
  cnetplot(PB, categorySize="pvalue", foldChange=genes)
  emapplot(PB, title = "Enrichment Map")
}

gene_symbols <- Cluster_1
converted_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, columns = "ENTREZID", keytype = "SYMBOL")

# enrichGO: Realiza análise de enriquecimento GO
# gene: Um vetor de IDs de genes para análise
# OrgDb: Especifica o banco de dados de anotações (ex.: org.Hs.eg.db para humanos)
# ont: Define a ontologia GO a ser usada (BP, MF, CC)
# pAdjustMethod: Método para ajuste de múltiplos testes (ex.: Benjamini-Hochberg)
# pvalueCutoff e qvalueCutoff: Definem os valores de corte para significância dos resultados