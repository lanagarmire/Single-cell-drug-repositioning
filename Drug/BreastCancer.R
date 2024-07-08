#Download datasets GSE113197 and GSE123926 from GEO before running this script.
#Dataset link: 
#GSE145926: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113197
#GSE158055: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123926

#Load Library
library('Asgard')
library('Seurat')

setwd("./")
#Load normal sample Ind5 from GSE113197 dataset
celltype<-read.table(file="Normal_celltype.txt",header = T,check.names=FALSE)
data<-read.table(file="GSM3099847_Ind5_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype2<-subset(celltype,sample=="Ind5" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype2))
data<-data[,common]
Epithelial2 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype2,cell=colnames(data),type="Normal"))

#Load normal sample Ind6 from GSE113197 dataset
data<-read.table(file="GSM3099848_Ind6_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype3<-subset(celltype,sample=="Ind6" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype3))
data<-data[,common]
Epithelial3 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype3,cell=colnames(data),type="Normal"))

#Load normal sample Ind7 from GSE113197 dataset
data<-read.table(file="GSM3099849_Ind7_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype4<-subset(celltype,sample=="Ind7" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype4))
data<-data[,common]
Epithelial4 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype4,cell=colnames(data),type="Normal"))

#Load cancer sample PDX110 from GSE123926 dataset
TNBC_PDX.data<- Read10X(data.dir = "GSM3516947_PDX110")
TNBC.PDX2 <- CreateSeuratObject(counts = TNBC_PDX.data, project = "TNBC", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(TNBC_PDX.data),cell=colnames(TNBC_PDX.data),sample="PDX-110",type="TNBC.PDX"))

#Load cancer sample PDX322 from GSE123926 dataset
TNBC_PDX.data<- Read10X(data.dir = "GSM3516948_PDX322")
TNBC.PDX3 <- CreateSeuratObject(counts = TNBC_PDX.data, project = "TNBC", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(TNBC_PDX.data),cell=colnames(TNBC_PDX.data),sample="PDX-332",type="TNBC.PDX"))

#single-cell alignment
SC.list<-list(TNBC.PDX2=TNBC.PDX2,TNBC.PDX3=TNBC.PDX3,Epithelial2=Epithelial2,Epithelial3=Epithelial3,Epithelial4=Epithelial4)
CellCycle=TRUE #Set it TRUE if you want to do Cell Cycle Regression
anchor.features=2000

for (i in 1:length(SC.list)) {
  SC.list[[i]] <- NormalizeData(SC.list[[i]], verbose = FALSE)
  SC.list[[i]] <- FindVariableFeatures(SC.list[[i]], selection.method = "vst",
                                       nfeatures = anchor.features, verbose = FALSE)
}
SC.anchors <- FindIntegrationAnchors(object.list = SC.list,anchor.features = anchor.features, dims = 1:15)
SC.integrated <- IntegrateData(anchorset = SC.anchors, dims = 1:15)
DefaultAssay(SC.integrated) <- "integrated"
if(CellCycle){
  ##Cell Cycle Regression
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  SC.integrated <- CellCycleScoring(SC.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SC.integrated <- ScaleData(SC.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SC.integrated))
  SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
}else{
  ##Run the standard workflow for visualization and clustering
  SC.integrated <- ScaleData(SC.integrated, verbose = FALSE)
  SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
}
##t-SNE and Clustering
SC.integrated <- RunUMAP(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindNeighbors(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindClusters(SC.integrated, algorithm = 1, resolution = 0.4)

##Cell Type Annotation, set by.CellType=TRUE if you want to annotate cell  type.
by.CellType=FALSE
if(by.CellType == TRUE){
  library('SingleR')
  data <- as.matrix(SC.integrated@assays$RNA@data)
  hpca.se <- HumanPrimaryCellAtlasData()
  pred.hpca <- SingleR(test = data, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
  cell.label <- data.frame(row.names = row.names(pred.hpca),celltype=pred.hpca$labels)
  if(length(SC.integrated@meta.data$celltype)>0){
    SC.integrated@meta.data$celltype <- cell.label$celltype
  }else{
    SC.integrated@meta.data <- cbind(SC.integrated@meta.data,cell.label)
  }
  new.cells <- data.frame()
  for(i in unique(SC.integrated$seurat_clusters)){
    sub.data <- subset(SC.integrated,seurat_clusters==i)
    temp <- table(sub.data@meta.data$celltype)
    best.cell <- names(which(temp==temp[which.max(temp)]))
    cells.temp <- data.frame(cell.id=row.names(sub.data@meta.data),celltype=best.cell)
    new.cells <- rbind(new.cells,cells.temp)
  }
  cell.meta <- SC.integrated@meta.data
  cell.id <- rownames(cell.meta)
  row.names(new.cells) <- new.cells[,1]
  new.cells <- new.cells[cell.id,]
  SC.integrated@meta.data$celltype <- new.cells$celltype
}else{
  SC.integrated@meta.data$celltype <- paste0("C",as.numeric(SC.integrated@meta.data$seurat_clusters))
}


#Change sample names
sample<-SC.integrated@meta.data$sample
sample[which(sample=="Ind5")]<-"Normal1"
sample[which(sample=="Ind6")]<-"Normal2"
sample[which(sample=="Ind7")]<-"Normal3"
SC.integrated@meta.data$sample<-sample

#View alignment results
pdf(file = "BC_plot.pdf",width = 12,height = 3)
DimPlot(SC.integrated, reduction = "umap", split.by = "sample",group.by = "celltype")
dev.off()
#Save data
saveRDS(SC.integrated,file="TNBC_SCdata.rds")

#Get differential genes from limma
SC.integrated=readRDS("TNBC_SCdata.rds")
Case=c("PDX-110","PDX-332")
Control=c("Normal1","Normal2","Normal3")
library('limma')
DefaultAssay(SC.integrated) <- "RNA"
set.seed(123456)
min.cells=3 # The minimum number of cells for a cell type. A cell type is omitted if it has less cells than the minimum number.
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.integrated@meta.data$celltype)){
  Idents(SC.integrated) <- "celltype"
  c_cells <- subset(SC.integrated, celltype == i)
  Idents(c_cells) <- "type"
  Samples=c_cells@meta.data
  Controlsample <- row.names(subset(Samples,sample %in% Control))
  Casesample <- row.names(subset(Samples,sample %in% Case))
  if(length(Controlsample)>min.cells & length(Casesample)>min.cells){
    expr <- as.matrix(c_cells@assays$RNA@data)
    new_expr <- as.matrix(expr[,c(Casesample,Controlsample)])
    new_sample <- data.frame(Samples=c(Casesample,Controlsample),type=c(rep("Case",length(Casesample)),rep("Control",length(Controlsample))))
    row.names(new_sample) <- paste(new_sample$Samples,row.names(new_sample),sep="_")
    expr <- new_expr
    bad <- which(rowSums(expr>0)<3)
    if(length(bad)>0){
      expr <- expr[-bad,]
    }
    mm <- model.matrix(~0 + type, data = new_sample)
    fit <- lmFit(expr, mm)
    contr <- makeContrasts(typeCase - typeControl, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrasts = contr)
    tmp <- eBayes(tmp)
    C_data <- topTable(tmp, sort.by = "P",n = nrow(tmp))
    C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$t,adj.P.Val=C_data$adj.P.Val,P.Value=C_data$P.Value)
    Gene.list[[i]] <- C_data_for_drug
    C_names <- c(C_names,i)
  }
}
names(Gene.list) <- C_names
#Save data
saveRDS(Gene.list,file="TNBC_genelist_limma.rds")

#Get differential genes from Seurat
SC.integrated=readRDS("TNBC_SCdata.rds")
library('Seurat')
DefaultAssay(SC.integrated) <- "RNA"
set.seed(123456)
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.integrated@meta.data$celltype)){
  Idents(SC.integrated) <- "celltype"
  c_cells <- subset(SC.integrated, celltype == i)
  Idents(c_cells) <- "type"
  C_data <- FindMarkers(c_cells, ident.1 = "TNBC.PDX", ident.2 = "Normal")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_logFC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
#Save data
saveRDS(Gene.list,file="TNBC_genelist_seurat.rds")

#Get differential genes from Seurat DESeq2
SC.integrated=readRDS("TNBC_SCdata.rds")
library('Seurat')
DefaultAssay(SC.integrated) <- "RNA"
set.seed(123456)
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.integrated@meta.data$celltype)){
  Idents(SC.integrated) <- "celltype"
  c_cells <- subset(SC.integrated, celltype == i)
  Idents(c_cells) <- "type"
  C_data <- FindMarkers(c_cells, ident.1 = "TNBC.PDX", ident.2 = "Normal", test.use = "DESeq2")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_logFC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
#Save data
saveRDS(Gene.list,file="TNBC_genelist_DESeq2.rds")

#Get differential genes from EdgeR
library('edgeR')
SC.integrated=readRDS("TNBC_SCdata.rds")
Case=c("PDX-110","PDX-332")
Control=c("Normal1","Normal2","Normal3")
DefaultAssay(SC.integrated) <- "RNA"
set.seed(123456)
min.cells=3 # The minimum number of cells for a cell type. A cell type is omitted if it has less cells than the minimum number.
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.integrated@meta.data$celltype)){
  Idents(SC.integrated) <- "celltype"
  c_cells <- subset(SC.integrated, celltype == i)
  Idents(c_cells) <- "type"
  Samples=c_cells@meta.data
  Controlsample <- row.names(subset(Samples,sample %in% Control))
  Casesample <- row.names(subset(Samples,sample %in% Case))
  if(length(Controlsample)>min.cells & length(Casesample)>min.cells){
    expr <- as.matrix(c_cells@assays$RNA@data)
    new_expr <- as.matrix(expr[,c(Casesample,Controlsample)])
    new_sample <- data.frame(Samples=c(Casesample,Controlsample),type=c(rep("Case",length(Casesample)),rep("Control",length(Controlsample))))
    row.names(new_sample) <- paste(new_sample$Samples,row.names(new_sample),sep="_")
    expr <- new_expr
    bad <- which(rowSums(expr>0)<3)
    expr <- expr[-bad,]
    group <- new_sample$type
    dge <- DGEList(counts=expr, group=group)
    group_edgeR <- factor(group,levels = c("Control","Case"))
    design <- model.matrix(~ group_edgeR)
    dge <- estimateDisp(dge, design = design)
    fit <- glmFit(dge, design)
    res <- glmLRT(fit)
    C_data <- res$table
    C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$logFC,adj.P.Val=p.adjust(C_data$PValue,method = "BH"),P.Value=C_data$PValue)
    Gene.list[[i]] <- C_data_for_drug
    C_names <- c(C_names,i)
  }
}
names(Gene.list) <- C_names
#Save data
saveRDS(Gene.list,file="TNBC_genelist_edgeR.rds")

#Drug repurposing
for(i in c("limma","edgeR","DESeq2","seurat")){
Gene.list<-readRDS(paste0("TNBC_genelist_",i,".rds"))
my_gene_info<-read.table(file="DrugReference/breast_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="DrugReference/breast_drug_info.txt",sep="\t",header = T,quote = "")
cmap.ref.profiles<-GetDrugRef(drug.response.path = 'DrugReference/breast_rankMatrix.txt',
                                 probe.to.genes = my_gene_info, drug.info = my_drug_info)
Drug.ident.res<-GetDrug(gene.data = Gene.list, drug.ref.profiles = cmap.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type = "FDA")
saveRDS(Drug.ident.res,file=(paste0("TNBC_drugs_FDA_",i,".rds")))
}

#Drug score
library(cmapR)
SC.integrated<-readRDS("TNBC_SCdata.rds")
Gene.list<-readRDS("TNBC_genelist_limma.rds")
Drug.ident.res<-readRDS("TNBC_drugs_limma.rds")
GSE92742.gctx.path="GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
Tissue="breast"
for(i in c("limma","edgeR","DESeq2","seurat")){
  Gene.list<-readRDS(paste0("TNBC_genelist_",i,".rds"))
  Drug.ident.res<-readRDS(paste0("TNBC_drugs_FDA_",i,".rds"))  
  Drug.score<-DrugScore(SC.integrated=SC.integrated,
                     Gene.data=Gene.list,
                     Cell.type=NULL,
                     Drug.data=Drug.ident.res,
                     FDA.drug.only=T,
                     Case=Case,
                     Tissue=Tissue,
                     GSE92742.gctx=GSE92742.gctx.path,
                     GSE70138.gctx=GSE70138.gctx.path)
   saveRDS(Drug.score,file=paste0("TNBC_drugscore_FDA_",i,".rds"))
}

#Drug repurposing using drugs/compounds
Gene.list<-readRDS("TNBC_genelist_limma.rds")
my_gene_info<-read.table(file="DrugReference/breast_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="DrugReference/breast_drug_info.txt",sep="\t",header = T,quote = "")
cmap.ref.profiles<-GetDrugRef(drug.response.path = 'DrugReference/breast_rankMatrix.txt',
                                 probe.to.genes = my_gene_info, drug.info = my_drug_info)
Drug.ident.res = GetDrug(gene.data = Gene.list, drug.ref.profiles = cmap.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type="all")
saveRDS(Drug.ident.res,file="TNBC_drugs_limma_all.rds")

#Drug/compound score
library(cmapR)
SC.integrated<-readRDS("TNBC_SCdata.rds")
Gene.list<-readRDS("TNBC_genelist_limma.rds")
Drug.ident.res<-readRDS("TNBC_drugs_limma_all.rds")
GSE92742.gctx.path="GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
Tissue="breast"

cell_metadata <- SC.integrated@meta.data
cell_metadata$cluster <- SC.integrated@meta.data$celltype

Drug.score<-DrugScore(cell_metadata, cluster_degs = Gene.list, 
                        cluster_drugs = Drug.ident.res, tissue = "breast", 
                        case = Case, gse92742_gctx_path = GSE92742.gctx.path, 
                        gse70138_gctx_path = GSE70138.gctx.path)
saveRDS(Drug.score,file="TNBC_drugscore_all.rds")

#Drug combination
library(cmapR)
SC.data<-readRDS("TNBC_SCdata.rds")
Gene.list<-readRDS("TNBC_genelist_limma.rds")
Drug.ident.res<-readRDS("TNBC_drugs_limma.rds")
GSE92742.gctx.path="GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
Tissue="breast"
Case=c("PDX-110","PDX-332")
Drug.combinations<-DrugCombination(SC.integrated=SC.data,
                      Gene.data=Gene.list,
                      Drug.data=Drug.ident.res,
                      Drug.FDR=0.05,
                      FDA.drug.only=TRUE,
                      Combined.drugs=2,
                      Case=Case,
                      Tissue="breast",
                      GSE92742.gctx=GSE92742.gctx.path,
                      GSE70138.gctx=GSE70138.gctx.path)
saveRDS(Drug.combinations,file="TNBC_drugcombinations.rds")
#Select mono-drugs
Drug.ident.res<-readRDS("TNBC_drugs_limma.rds")
Final.drugs<-TopDrug(SC.integrated=SC.integrated,
                   Drug.data=Drug.ident.res,
                   Drug.FDR=0.05,
                   FDA.drug.only=TRUE,
                   Case=Case
)
saveRDS(Final.drugs,file="TNBC_selected_drugs_limma.rds")
#Select drug combinations
Final.combinations<-TopCombination(Drug.combination=Drug.combinations,
                   Combination.FDR=0.05,
                   Min.combination.score=1
)
saveRDS(Final.combinations,file="TNBC_selected_drugcombinations.rds")
Final.combinations<-readRDS("TNBC_selected_drugcombinations.rds")

#Personalized drug combination
SC.data<-readRDS("TNBC_SCdata.rds")
Gene.list<-readRDS("TNBC_genelist.rds")
Drug.ident.res<-readRDS("TNBC_drugs.rds")
GSE92742.gctx.path="GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
Tissue="breast"
Case="PDX-110"
Drug.combinations<-DrugCombination(SC.integrated=SC.data,
                                   Gene.data=Gene.list,
                                   Drug.data=Drug.ident.res,
                                   Drug.FDR=0.05,
                                   FDA.drug.only=TRUE,
                                   Combined.drugs=2,
                                   Case=Case,
                                   Tissue="breast",
                                   GSE92742.gctx=GSE92742.gctx.path,
                                   GSE70138.gctx=GSE70138.gctx.path)
PDX1.combinations<-TopCombination(Drug.combination=Drug.combinations,
                                   Combination.FDR=0.05,
                                   Min.combination.score=1
)
Case="PDX-332"
Drug.combinations<-DrugCombination(SC.integrated=SC.data,
                                   Gene.data=Gene.list,
                                   Drug.data=Drug.ident.res,
                                   Drug.FDR=0.05,
                                   FDA.drug.only=TRUE,
                                   Combined.drugs=2,
                                   Case=Case,
                                   Tissue="breast",
                                   GSE92742.gctx=GSE92742.gctx.path,
                                   GSE70138.gctx=GSE70138.gctx.path)
PDX2.combinations<-TopCombination(Drug.combination=Drug.combinations,
                                  Combination.FDR=0.05,
                                  Min.combination.score=1
)
