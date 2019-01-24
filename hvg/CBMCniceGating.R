
# Load data ---------------------------------------------------------------

citeDIR <- "~/sds/sd17l002/p/scRNAseq_datasets/CITEseq_NatMethods_2017/data/"
rna_file <- paste0(citeDIR, "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz")
adt_file <- paste0(citeDIR, "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz")

rawC <- as.matrix(read.csv(gzfile(rna_file), row.names = 1))
prot <- as.matrix(read.csv(gzfile(adt_file), row.names = 1))

# exclude mouse cells:
is_mouse <- colSums(rawC[grepl("MOUSE", rownames(rawC)),]) / colSums(rawC) > .1
rawC <- rawC[grepl("HUMAN", rownames(rawC)), ! is_mouse]
prot <- prot[, ! is_mouse]




# Normalize / PCA / UMAP---------------------------------------------------------------


normC <- log1p(t( t(rawC) / colSums(rawC)))
# simplifying Seurat's code to compute Centered Log Ratio (CLR):
norm_prot <- apply(prot, 1, function(x) {
  log1p( (x) /
           (exp(sum(log1p((x)[x > 0]), na.rm = TRUE)/length(x + 1))) ) })


library(umap)
library(irlba)
Uprot <- umap(norm_prot)
pca   <- irlba::prcomp_irlba(t(normC), n = 50)
Urna <- umap(pca$x)





# Seurat run (vignette) ---------------------------------------------------
library(Seurat)
# From Seurat's multimodal vignette:
# (more var.genes than default, no regression, 13 PCs)
cbmc <- CreateSeuratObject(raw.data = rawC)
cbmc <- MakeSparse(cbmc)
cbmc <- NormalizeData(cbmc, display.progress = FALSE)
cbmc <- FindVariableGenes(cbmc, do.plot = FALSE, y.cutoff = 0.5)
cbmc <- ScaleData(cbmc, display.progress = FALSE, display.progress = FALSE)
cbmc <- RunPCA(cbmc, pcs.print = 0)
# cbmc <- FindClusters(cbmc, dims.use = 1:13, print.output = FALSE)
cbmc <- RunTSNE(cbmc, dims.use = 1:13)
cbmc <- RunUMAP(cbmc, dims.use = 1:13)











# Define groundtruth cell types -------------------------------------------

# This is informed by my own PBMC-based knowledge plus diverse paper searches.
# I feel like ideally I would double-check this with a biologists who knows
# CBMCs well, or confirm it with a public FACS dataset from flowrepository.

# The goal here is to identify cell types with NO false-positives, so I kick
# out a few cells too many rather than labeling a cell with the wrong cell type.

# I positively select erythrocytes using mRNA markers and then kick out all doublets /
# multiplets / other cells using protein markers.
# 

# For excluding erythrocytes and platelets / megakaryocytes, there are
# no good protein marker available, we use RNA for it:
not_ery <- normC["HUMAN_HBB", ]  < .05  &
  normC["HUMAN_HBG2", ] < .05 &
  normC["HUMAN_HBA1", ] < .04 
not_platelet <- normC["HUMAN_GP9", ] < .005 &
  normC["HUMAN_PF4", ] < .005 &
  normC["HUMAN_PPBP", ] < .005 

# Note that this might also exclude doublets and multiplets containing
# erythrocytes, so I'll define high-confidence erys further below.



# ...cell type definitions (high-confidence) -----------------------------------------



# T cells 
# Surprising: CD4 T cells seem to express CD14!
# FYI: The largest isoform CD45RA is expressed on naïve T cells. Activated and
#      memory T lymphocytes express the shortest CD45 isoform, CD45RO, which
#      lacks RA, RB, and RC exons. (Krzywinska, PLoS One 2016)
Tcell <- not_ery       &
    not_platelet  &
    norm_prot[, "CD19"] < 1.5 &
    norm_prot[, "CD56"] < 1.1 &
    norm_prot[, "CD11c"] < 1  &
    norm_prot[, "CD3" ] > 1   &
    (
      (norm_prot[, "CD4"] > 1  &  norm_prot[, "CD8"] < 1 )  |
        (norm_prot[, "CD4"] < .5  &  norm_prot[, "CD8"] > 3 )
    )  &
    norm_prot[, "CD34" ] < 1   # informed by UMAP on norm_prot
T4      <- Tcell & norm_prot[, "CD4"] > 1  &  norm_prot[, "CD8"] < 1 
T8      <- Tcell & norm_prot[, "CD4"] < .5  &  norm_prot[, "CD8"] > 3 


 
Bcell <- not_ery       &
    not_platelet  &
    norm_prot[, "CD19"] > 2 &
    norm_prot[, "CD56"] < 1.1 &
    norm_prot[, "CD11c"] < 1  &
    norm_prot[, "CD3"]   < .5

HSC <- not_ery       &
# https://www.stemcell.com/human-hematopoietic-stem-and-progenitor-cell-phenotyping-panels.html :
    not_platelet  &
    norm_prot[, "CD19"] < 1.5 &
    norm_prot[, "CD56"] < 1.1 &
    norm_prot[, "CD11c"] < 1  &
    norm_prot[, "CD34"] > 1.5 &
    # surprising that CD11c- HSCs have very high colsums...
    norm_prot[, "CD3"] < .5  # informed by UMAP on norm_prot

NK0 <- not_ery & not_platelet &
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4839597/ (Krzywinska, PLoS One 2016):
      norm_prot[, "CD56"] > 1  &
      norm_prot[, "CD3"]  < .5 &
      norm_prot[, "CD45RA"] > 1 &
      norm_prot[, "CD19"] < 1.5 &
      norm_prot[, "CD16"] > 1
      # 20 % CD11c+ NK cells is normal in adults (http://www.jimmunol.org/content/177/8/5659.long, Aranami 2006),
      # but I find they systematically have larger colSums and are CD14+, so I'm suspicious and
      # keep them separate:
NK_CD11cHigh <- NK0 & norm_prot[, "CD11c"] > 1.5  
NK           <- NK0 & norm_prot[, "CD11c"] < 1.2


# I am quite confident these are pure erys, because they are grouped in
# UMAP etc.
ery <- !not_ery &
  norm_prot[, "CD19"] < 1 &
  norm_prot[, "CD3"]  < .5  &
  norm_prot[, "CD56"] < 1 &
  norm_prot[, "CD14"]  < .5  &
  norm_prot[, "CD16"]  < 1  &
  norm_prot[, "CD34"]  < 1  



M14 <- not_ery & not_platelet &
# negative for CD3, CD8, CD19, CD56 and CD4 (https://onlinelibrary.wiley.com/doi/pdf/10.1002/cyto.b.20609)
  norm_prot[, "CD14"] > 1 &
  norm_prot[, "CD16"]  < 1  &
  norm_prot[, "CD19"] < 1.5 &
  norm_prot[, "CD56"] < 1 &
  norm_prot[, "CD3"]  < .5  &
  norm_prot[, "CD4"]  < 1  &
  norm_prot[, "CD8"]  < 1.5  &
  norm_prot[, "CD34"]  < 1.3 
  # once you find papers supporting it, also require CD45RA < .5? (CD11c > 1 is already.)




# I can't seem to confidently define intermediate Monos, ncMonos and DCs.
# To benchmark Mono14 clusters/classifications, however, I do include cells
# that certainly are not classical monocytes, and might be one of these:
notM14_perhapsDC  <-
  normC["HUMAN_HLA-DRA", ] > quantile(normC["HUMAN_HLA-DRA", ],
                                      probs = 1 - 100/length(normC["HUMAN_HLA-DRA",])) &
  norm_prot[, "CD14"] < .3 &
  norm_prot[, "CD19"] < 1 
  

notM14_perhapsM16 <- not_ery & not_platelet &
  norm_prot[, "CD14"] > .5 &
  norm_prot[, "CD11c"] > 1  &
  norm_prot[, "CD16"] > 1.2  &
  norm_prot[, "CD19"] < 1.5 &
  norm_prot[, "CD56"] < 1 &
  norm_prot[, "CD3"]  < .5  &
  norm_prot[, "CD4"]  < 1  &
  norm_prot[, "CD8"]  < 1.5  &
  norm_prot[, "CD34"]  < 1.3 



# ... uncertain cell type definitions -------------------------------------


T_naive <- Tcell & norm_prot[, "CD45RA"] > 1.5
T_mem   <- Tcell & norm_prot[, "CD45RA"] < .5
 

platelet1 <- !not_platelet &
  norm_prot[, "CD19"] < 1 &
  norm_prot[, "CD3"]  < .5  &
  norm_prot[, "CD56"] < 1 &
  norm_prot[, "CD14"]  < .8  &
  norm_prot[, "CD16"]  < 1  &
  norm_prot[, "CD34"]  < 1  &
  norm_prot[, "CD11c"] < 1  
platelet2 <- !not_platelet &
  norm_prot[, "CD19"] > 1 &
  norm_prot[, "CD3"]  > .5  &
  norm_prot[, "CD56"] > 1 &
  norm_prot[, "CD14"] >  .8  &
  norm_prot[, "CD16"] >  1  &
  norm_prot[, "CD34"] >  1  &
  norm_prot[, "CD11c"] >  1  





# DCs https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1594912/ (Breitling, Infect Immun 2006),
# also Sorg, Blood 1999:
# HLA-DRA+
# negative for CD3, CD14, CD16, CD19, CD20, and CD56
# cDCs are CD1C+, pDCs are CLEC4C+ 
oldDC  <- not_ery       &
    not_platelet  &
    norm_prot[, "CD19"] < 1.5 &
    norm_prot[, "CD56"] < 1.1 &
    norm_prot[, "CD3"]  < .5  &
    norm_prot[, "CD16"]  < .5  &
    norm_prot[, "CD14"]  < .3  &
    norm_prot[, "CD34"]  < 1  
# we have enriched for DC markers:
table(oldDC, CD1C=normC["HUMAN_CD1C", ]!=0)    
table(oldDC, HLA_DR=normC["HUMAN_HLA-DRA", ]!=0)    

# ...doublets ----------------------------------------------------------------

# I want to discriminate 'multiplets' from 'simple doublets':
# multiplets are commonly observed in the CITEseq data and have > 2 markers coexpressed,
#            often all of them
# doublets have 2 markers coexpressed but not the others
doublet_TB <- norm_prot[, "CD3"] > 1 &
              norm_prot[, "CD19"] > 2  &
              norm_prot[, "CD34"] < 1.5 &
              norm_prot[, "CD56"] < 1.5 &
              norm_prot[, "CD11c"] < 1
# I was slightly worried that CD3+CD19+ double-positive cells might actually
# be a normal thing in chord blood, but quick online searches shows that people
# routinely sort T cells by depleting B cells via CD19 - see for example:
# Nomura, Experimental Hematology 2001
# (https://www.sciencedirect.com/science/article/pii/S0301472X01006890?via%3Dihub)

doublet_TNK <- norm_prot[, "CD3"] > 1 &
              norm_prot[, "CD19"] < 1.5  &
              norm_prot[, "CD34"] < 1.5 &
              norm_prot[, "CD56"] > 1.5 &
              norm_prot[, "CD11c"] < 1
doublet_T48 <- norm_prot[, "CD4"] > 1 &
               norm_prot[, "CD8"] > 2 &
               norm_prot[, "CD19"] < 1.3 &
               norm_prot[, "CD56"] < 1 &
               norm_prot[, "CD11c"] < 1      

doublet_T4Mono <- norm_prot[, "CD3"] > 1 &
  # not entirely sure about this one. CD8 T cells are reported to
  # sometimes express CD11c (in disease). In chord blood, however,
  # I have heard no mention of it, so perhaps these doublets are valid.
                  norm_prot[, "CD11c"] > 1 &
                  norm_prot[, "CD8"] < 1.5 &
                  norm_prot[, "CD19"]< 1.5




# multiplets, as opposed to 'simple doublets', express MORE than 2 mutually-exclusive
# markers:
multiplets  <- rowSums(
               cbind(norm_prot[, "CD3"] > .5,
                     norm_prot[, "CD19"] > 2,
                     norm_prot[, "CD56"] > 1.5,
                     norm_prot[, "CD34"] > 1.5)) > 2







# ...Groundtruth Summary ----------------------------------------------------------


# High-confidence classes:  7x celltype,
#                           2x anti-celltype (notM14),
#                           4x doublets, 
#                           1x multiplets

# as dataframe:
groundtruth <- data.frame(
                 ery = ery,
                 #Tcell = Tcell, 
                 T4=T4, 
                 T8 = T8, 
                 doublet_T48=doublet_T48,
                 Bcell=Bcell,
                 HSC=HSC,
                 NK = NK,
                 #NK_CD11cHigh=NK_CD11cHigh,
                 #mono=mono,
                 M14=M14,
                 notM14_perhapsDC=notM14_perhapsDC,
                 notM14_perhapsM16=notM14_perhapsM16,
                 doublet_T4Mono=doublet_T4Mono,
                 doublet_TB = doublet_TB,
                 doublet_TNK = doublet_TNK,
                 multiplets = multiplets
                 )
# collapsed into a single vector:
groundtruth_vector <- rep("undefined", nrow(groundtruth))
setCT <- function(vector, ct) {vector[groundtruth[, ct]] <- ct; return(vector)}
for(ct in colnames(groundtruth)){
  groundtruth_vector <<- setCT(groundtruth_vector, ct)
}

groundtruth_vector <- factor(groundtruth_vector,
                             levels =c(
                                 "HSC",
                                 "ery",
                                 "Bcell",
                                 "doublet_TB",
                                 "T4",   "doublet_T48",
                                 "T8",      
                                 "doublet_TNK",
                                 "doublet_T4Mono",
                                 "NK",
                                 "M14",
                                
                                 "notM14_perhapsDC",
                                 "notM14_perhapsM16",
                                 
                                 "multiplets",
                                 "undefined"
                                 ))




ggplotScale_groundtruth <- scale_color_manual(values=c(
  "T4"               = "#B2DF8A",
  "T8"               = "#33A02C",
  "doublet_T48"       = "red",
  "Bcell"            = "#AE017E",
  "HSC"              = "black",
  "NK"               = "#F16913",
  "ery"              = "darkred",
  "M14"              = "#4292C6",
  "notM14_perhapsDC" = "#08519C",
  "notM14_perhapsM16"= "#08519C",
  "doublet_TB"       = "red",
  "doublet_TNK"      = "red",
  "doublet_T4Mono"   = "red",
  "multiplets"       ="red",
  "undefined"        = adjustcolor("grey", alpha.f = .3)))
  
  
  


# ...summary (upsetR) --------------------------------------------------------
library(UpSetR)
upset( groundtruth * 1,
      nsets = 18
)







# Shortcut Functions ---------------------------------------------------------------
# functions:
   U <- function(ids, umap=Uprot, clr = "black"){
     plot(umap$layout, pch=20, col = adjustcolor("grey", alpha.f = .3))
     # drop=FALSE necessary when length(ids)is 1:
     points(umap$layout[ids,, drop=FALSE], pch=20, col = clr)
   }
   gate <- function(x = "CD3", y = "CD19", ids=1:100, lowCol="grey", highCol="red"){
    plot(norm_prot[, x], norm_prot[, y], pch=20, col = lowCol, xlab = x, ylab=y)
    points(norm_prot[ids, x, drop=FALSE], norm_prot[ids, y, drop=FALSE], pch=20, col = highCol) 
   }
   prot_colSums <- colSums(prot)
   rna_colSums <- colSums(rawC)
   G <- function(ids, clr = scales::muted("red")){
     par(mfrow = c(2,3))
     gate("CD3", "CD19", ids, highCol = clr)
     gate("CD14", "CD16", ids, highCol = clr)
     gate("CD4", "CD8", ids, highCol = clr)
     gate("CD34", "CD56", ids, highCol = clr)
     gate("CD45RA", "CD11c", ids, highCol = clr)
     plot(prot_colSums, rna_colSums, log = "xy", pch=20, col = adjustcolor("grey", alpha.f = .4))
     points(prot_colSums[ids], rna_colSums[ids], pch=20, col = clr)
     # plot(normC["HUMAN_HLA-DRA", ], normC["HUMAN_MKI67",], pch=20, col = "grey")
     # points(normC["HUMAN_HLA-DRA", ids], normC["HUMAN_MKI67", ids], pch=20, col = clr)
     par(mfrow = c(1,1))}

plotDF <- function(seurat){ data.frame( seurat@dr$tsne@cell.embeddings,
                            seurat@dr$umap@cell.embeddings,
                            seurat@meta.data,
                            TrueClass=groundtruth_vector,
                            groundtruth
)}
S <- function(ids, clr = "black"){
  df <- plotDF(cbmc)
  par(mfrow=c(1,2))
  plot(df$tSNE_1, df$tSNE_2, pch=20, col = adjustcolor("grey", alpha.f = .3),
       main = "Seurat (vignette)")
  points(df$tSNE_1[ids,drop=F], df$tSNE_2[ids, drop=F], pch=20, col = clr)
  plot(df$UMAP1, df$UMAP2, pch=20, col = adjustcolor("grey", alpha.f = .3))
  points(df$UMAP1[ids,drop=F], df$UMAP2[ids, drop=F], pch=20, col = clr)
  par(mfrow=c(1,1))
}
  
  

V <- function(gene = "CD3E"){
  gene <- paste0("HUMAN_", gene)
 p <-
  ggplot(data.frame(Groundtruth = groundtruth_vector, Expression = normC[gene,]),
       aes(Groundtruth, Expression, col = Groundtruth))+
  geom_jitter() + ggtitle(gene)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
   ggplotScale_groundtruth + guides(col=FALSE)
 return(p)
}






# Helper functions for classification -------------------------------------

tids <- function(gene = "CD3E", n = 20){
  # top IDs: from the cells expressing this gene at all, get the indices of the n highest.
  expr <- normC[paste0("HUMAN_", gene), ] 
  expressors <- sum(expr > 0)
  if(n > expressors){warning(paste0("Max n reached; there are only ",
                                    expressors,
                                    " cells expressing ", gene))
    n <- expressors}
  order(expr, decreasing = T)[1:n]
}
removeOverlaps <- function( vectorList = list(A=c("A", "B", "C"), B = c("X", "Y", "B"))) {
  # from each vector in vectorList, the elements occuring in any other vector 
  # are removed.
  upset(fromList(vectorList), nsets = length(vectorList))
  multipleOccurences <- unlist(vectorList)[ duplicated(unlist(vectorList))]
  multipleOccurences <- unique(multipleOccurences)
  vectorList <- lapply(vectorList, function(v) { v <- v[! v %in% multipleOccurences]})
  return(vectorList)
}





# Savepoint ---------------------------------------------------------------

#save.image(file = "/home/felix/PhDother/scAnalysis/sc_methods/hvg/savepoint/start_CBMCniceGating.rda")
load(file.path("/home/felix/PhDother/scAnalysis/sc_methods/hvg/savepoint/",
               "start_CBMCniceGating.rda"))
library(umap)
library(glmnet)
library(irlba)
library(Seurat)
library(tidyverse)
library(UpSetR)
library(pheatmap)
library(RColorBrewer)









# Confirm B/T doublets ----------------------------------------------------



genesB <- c(  
  "HUMAN_CD19",
  "HUMAN_MS4A1",
  "HUMAN_CD79A",
  "HUMAN_CD79B",
  "HUMAN_IGHM",
  "HUMAN_IGHD"
  )

genesT <- c(
  "HUMAN_TRAC",
  "HUMAN_TRAT1",
  "HUMAN_TRBC2",
  "HUMAN_CD3E",
  "HUMAN_CD3G",
  "HUMAN_CD3D"
)


# expression levels:
t(normC[genesT, ]) %>% as.data.frame %>% rownames_to_column(var = "Cell") %>%
as_tibble() %>% mutate(Groundtruth = groundtruth_vector) %>%
  gather("Gene", "Expression", -Groundtruth, -Cell) %>%
  filter(Groundtruth %in% c("Bcell", "T4", "T8", "doublet_TB")) %>%
ggplot()+geom_jitter(aes(Groundtruth, Expression, col = Groundtruth))+
  facet_wrap(~Gene, scales = "free_y")

# percent of expressing cells:
t(normC[genesT, ]) %>% as.data.frame %>% rownames_to_column(var = "Cell") %>%
  as_tibble() %>% mutate(Groundtruth = groundtruth_vector) %>%
  gather("Gene", "Expression", -Groundtruth, -Cell) %>%
  filter(Groundtruth %in% c("Bcell", "T4", "T8", "doublet_TB")) %>%
  group_by(Groundtruth, Gene) %>%
  summarize(percent_expressing = sum(Expression>0)/length(Expression)) %>%
  ggplot()+geom_point(aes(Groundtruth, percent_expressing, col=Groundtruth))+
  facet_wrap(~Gene, scales = "free_y")







# Seurat, B/T doublets  --------------------------------------------------------

TSNEPlot(cbmc)
DimPlot(cbmc, reduction.use = "umap")

# B/T-cell doublets fail spectacularly:
  ggplot(plotDF(cbmc)) +
  geom_point(aes(tSNE_1, tSNE_2, col = doublet_TB)) +
  scale_color_manual(values = c('TRUE'="black", 'FALSE'=adjustcolor("grey", alpha.f = .2)))
# all at once:
  ggplot(plotDF(cbmc)) +
  geom_point(aes(tSNE_1, tSNE_2, col = TrueClass)) + ggplotScale_groundtruth


  


  
  
  
  
  
  

# multinomial playground -------------------------------------------------------------


# We envision a workflow where the researcher selects cells based on marker genes
# and informed by UMAP/tSNE. The below method is as simple as it gets:

# select training cells ('goldstandard', 'labelled data'):
trainIDs <- list(
  T = tids("CD3E", 300),
  B = tids("CD79B", 150),
  M14 = tids("CD14", 200),
  CD16 = tids("FCGR3A", 100),
  NK = tids("NCAM1", 40),
  E =tids("HBB", 50),
  HSC = tids("CD34", 20),
  Plt = tids("GP9", 30))
  
trainIDs <- removeOverlaps(trainIDs)


Labels <- unlist(lapply(names(trainIDs),
                         function(class) rep(class, length(trainIDs[[class]])) ))
IDsForLabels <- unlist(trainIDs)

table(Truth=groundtruth_vector[IDsForLabels], Train=Labels)
# Amazing how many we get right with very simple means of marker expression!
# Only CD16 does not work, and we'll notice this below.

# train model:
cv<-cv.glmnet(x = pca$x[IDsForLabels, ],
          y = Labels,
          family = "multinomial",
          alpha = .5)






# Score jitters

pred <- as.data.frame(
          predict(cv,
                  pca$x,
                  type ="response")[,,1]
  )
p <- ggplot(data.frame(Groundtruth = groundtruth_vector, pred)) + 
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     ggplotScale_groundtruth
p + geom_jitter(aes(Groundtruth, T, col=Groundtruth), width = .3)
#








# classification scores:
  pred <- as.data.frame(
          predict(cv,
                  pca$x,
                  type ="response")[,,1]
  )
# for plotting we massage this:
  rownames(pred) <- colnames(normC)
  pred$Truth <- groundtruth_vector
  pred <- droplevels(pred[pred$Truth != "undefined",])
  

  
# make heatmap
pheatmap(t(pred[order(pred$Truth), colnames(pred)!="Truth"]), cluster_rows = F, cluster_cols = F,
           annotation_col= pred[order(pred$Truth), "Truth", drop=FALSE],
           color = rev(RColorBrewer::brewer.pal(8, "RdBu")),
           breaks=0:8/8,
           show_rownames = TRUE, show_colnames = FALSE)
# overall impression is good!
  
    
# doublets are handled well:  
plot(pred$B, pred$T, col = case_when(pred$Truth == "doublet_TB" ~ "red",
                                     pred$Truth %in% c("T4", "T8") ~ "#B2DF8A",
                                     pred$Truth == "Bcell" ~  "black",
                                     TRUE ~ adjustcolor("grey", alpha.f = .4)),
     pch=20); legend(
       "topright",
       legend = c("T", "B", "doublet_TB"), pch=20, col = c("#B2DF8A",1, "red"),
       inset=.1)
plot(pred$M14, pred$T, col = ifelse(pred$Truth == "doublet_T4Mono", "red",  adjustcolor("grey", alpha.f = .4)), pch=20,
     main = "Doublets of T4 and Monocyte")
  





  
  
# Number of correct classifications:  
  
# maximal response decides class:
res <- table(
    apply(pred[, colnames(pred)!="Truth"],
          1,
          function(x) colnames(pred)[which(x==max(x))]),
    groundtruth_vector[groundtruth_vector!="undefined", drop=T]
)# t(apply(res,1,function(x) {round(100*x/sum(x))}))  


# response > .9
res <- table(apply(pred[, colnames(pred)!="Truth"],
                   1,
                   function(x) {
                     x <- x > .9
                     ifelse(length(which(x))==0, NA, colnames(pred)[which(x)])
                   } ), groundtruth_vector[groundtruth_vector!="undefined", drop=T], useNA = "always")
# t(apply(res,1,function(x) {round(100*x/sum(x))}))  
plot(pred$CD16, pred$NK, col = case_when(groundtruth_vector =="NK" ~ "blue",
                                         groundtruth_vector == "notM14_perhapsM16" ~ "green",
                                         TRUE ~ adjustcolor("grey", alpha.f = .2)),
     pch = 20)
# NK cells have the problem that classes "NK" and "CD16" are both trained on
# groundtruth NK cells (because selection of training data is poor).
# Hence the bad performance, I suppose.


res

# Seurat clustering:
# res <- as.matrix(table(cbmc@meta.data$res.0.8[groundtruth_vector!="undefined"], groundtruth_vector[groundtruth_vector!="undefined", drop=T]))



pheatmap(t(apply(res,1,function(x) {round(100*x/sum(x), 1)}))  ,
         breaks = -1:100, cluster_cols = F, cluster_rows = F,
         color = c("grey30", colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(101)),
         display_numbers=T, number_format = "%.1f",
         fontsize_number = 12, number_color = 1)


         
 




































# Downsample, CD8/CD4 T cells, ... ----------------------------------------




# downsample:
rawC <- matrix(rbinom(rawC, rawC, .8), ncol=ncol(rawC),
               dimnames = list(rownames(rawC), colnames(rawC)))

