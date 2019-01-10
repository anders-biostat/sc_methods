




# exclude erythrocytes. No good protein marker available, we use RNA for it:
not_ery <- normC["HUMAN_HBB", ]  < .05  &
  normC["HUMAN_HBG2", ] < .05 &
  normC["HUMAN_HBA1", ] < .04 
# same for platelets / megakaryocytes:
not_platelet <- normC["HUMAN_GP9", ] < .005 &
  normC["HUMAN_PF4", ] < .005 &
  normC["HUMAN_PPBP", ] < .005 


Bcells <- norm_prot[, "CD19"] > 2 & not_ery & not_platelet    #  &  norm_prot[, "CD3"] < .5
Tcells <- norm_prot[, "CD3"] > 1  & not_ery & not_platelet    #  &norm_prot[, "CD19"] < 1.5 
HSCs   <- norm_prot[, "CD34"] > 2 & not_ery & not_platelet 
# monocyte exclusion as in paper https://onlinelibrary.wiley.com/doi/pdf/10.1002/cyto.a.20859
NKcells<- norm_prot[, "CD14"] < .4 & norm_prot[, "CD56"] > 1 & not_ery & not_platelet 
DCs    <- norm_prot[, "CD14"] < .4 & normC["HUMAN_HLA-DRA", ] > .005
Mono14 <- norm_prot[, "CD14"] > .6 & normC["HUMAN_HLA-DRA", ] > .001 & norm_prot[, "CD16"] < 1
Mono16 <- norm_prot[, "CD14"] > .6 & normC["HUMAN_HLA-DRA", ] > .001 & norm_prot[, "CD16"] > 1.5





#        T O   D O :
#
#
#        exclude multiplets
#
#




# functions:
gate <- function(x = "CD3", y = "CD19", ids=1:100, lowCol="grey", highCol="red"){
 plot(norm_prot[, x], norm_prot[, y], pch=20, col = lowCol, xlab = x, ylab=y)
 points(norm_prot[ids, x], norm_prot[ids, y], pch=20, col = highCol) 
}
G <- function(ids, clr = scales::muted("red")){
  par(mfrow = c(2,3))
  gate("CD3", "CD19", ids, highCol = clr)
  gate("CD14", "CD16", ids, highCol = clr)
  gate("CD4", "CD8", ids, highCol = clr)
  gate("CD34", "CD56", ids, highCol = clr)
  gate("CD45RA", "CD11c", ids, highCol = clr)
  plot(normC["HUMAN_HLA-DRA", ], normC["HUMAN_MKI67",], pch=20, col = "grey")
  points(normC["HUMAN_HLA-DRA", ids], normC["HUMAN_MKI67", ids], pch=20, col = clr)
  par(mfrow = c(1,1))}

G(Tcells)

