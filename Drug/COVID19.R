#Download dataset GSE145926 andd GSE158055 from GEO before running this script.
#Dataset link: 
#GSE145926: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926
#GSE158055: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158055

#Load library
library('Asgard')
library('Seurat')

setwd("./")

#Load data

data<-Read10X_h5(file="GSM4339771_C143_filtered_feature_bc_matrix.h5")
C143 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C143",type="Severe.COVID19",outcome="Deceased"))
C143[["percent.mt"]] <- PercentageFeatureSet(C143, pattern = "^MT-")
C143 <- subset(C143, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
saveRDS(C143,file="C143.rds")

data<-Read10X_h5(file="GSM4339773_C145_filtered_feature_bc_matrix.h5")
C145 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C145",type="Severe.COVID19",outcome="Cured"))
C145[["percent.mt"]] <- PercentageFeatureSet(C145, pattern = "^MT-")
C145 <- subset(C145, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
saveRDS(C145,file="C145.rds")

data<-Read10X_h5(file="GSM4339774_C146_filtered_feature_bc_matrix.h5")
C146 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C146",type="Severe.COVID19",outcome="Deceased"))
C146[["percent.mt"]] <- PercentageFeatureSet(C146, pattern = "^MT-")
C146 <- subset(C146, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
saveRDS(C146,file="C146.rds")

data<-Read10X_h5(file="GSM4475051_C148_filtered_feature_bc_matrix.h5")
C148 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C148",type="Severe.COVID19",outcome="Cured"))
C148[["percent.mt"]] <- PercentageFeatureSet(C148, pattern = "^MT-")
C148 <- subset(C148, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
saveRDS(C148,file="C148.rds")

data<-Read10X_h5(file="GSM4475052_C149_filtered_feature_bc_matrix.h5")
C149 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C149",type="Severe.COVID19",outcome="Cured"))
C149[["percent.mt"]] <- PercentageFeatureSet(C149, pattern = "^MT-")
C149 <- subset(C149, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
saveRDS(C149,file="C149.rds")

data<-Read10X_h5(file="GSM4475053_C152_filtered_feature_bc_matrix.h5")
C152 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C152",type="Severe.COVID19",outcome="Cured"))
C152[["percent.mt"]] <- PercentageFeatureSet(C152, pattern = "^MT-")
C152 <- subset(C152, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
saveRDS(C152,file="C152.rds")

Cells <- read.csv("GSE158055_COVID19/GSE158055_cell_annotation.csv")
Cells <- subset(Cells, sampleID %in% c("S-S086-1","S-S085-1","S-S087-1","S-S088-1","S-S089-1","S-S090-1","S-S006","S-S008","S-S009"))
data<-Read10X(data.dir = "GSE158055_COVID19/selected_samples",gene.column = 1)
Z.data <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample=Cells$sampleID,type="Severe.COVID19",outcome="Cured"))
Z.data[["percent.mt"]] <- PercentageFeatureSet(Z.data, pattern = "^MT-")
Z.data <- subset(Z.data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
Z.data@meta.data$outcome[which(Z.data@meta.data$sample %in% c("S-S086-1","S-S087-1"))]="Deceased"


#Merge data
Severe.COVID19<-merge(x=C143,y=c(C145,C146,C148,C149,C152))

#single-cell alignment and annotation
SC.list<-list(Z.data=Z.data,Severe.COVID19=Severe.COVID19)
CellCycle=TRUE #Set it TRUE if you want to do Cell Cycle Regression
anchor.features=2000
by.CellType=TRUE

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
saveRDS(SC.integrated,file="COVID_SCdata.rds")
SC.integrated=readRDS('COVID_SCdata.rds')

##Cell Type Annotation, set by.CellType=TRUE if you want to annotate cell  type.
if(by.CellType == TRUE){
  library('SingleR')
  data <- SC.integrated@assays$RNA@data
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
#Save data
saveRDS(SC.integrated,file="COVID_SCdata_annotated.rds")

#Creat the vector of sample names
SC.data=readRDS('COVID_SCdata_annotated.rds')
Deceased.Severe<-c("COVID19.C143","COVID19.C146","S-S086-1","S-S087-1")
Cured.Severe<-setdiff(unique(SC.data@meta.data$sample),Deceased.Severe)

#Get differential genes from limma
library('limma')
DefaultAssay(SC.data) <- "RNA"
set.seed(123456)
Case=Deceased.Severe
Control=Cured.Severe
min.cells=3 # The minimum number of cells for a cell type. A cell type is omitted if it has less cells than the minimum number.
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.data@meta.data$celltype)){
  Idents(SC.data) <- "celltype"
  c_cells <- subset(SC.data, celltype == i)
  Idents(c_cells) <- "outcome"
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
saveRDS(Gene.list,file="COVID19_genelist_limma.rds")

#Get differential genes from Seurat
SC.data=readRDS("COVID_SCdata_annotated.rds")
library('Seurat')
DefaultAssay(SC.data) <- "RNA"
set.seed(123456)
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.data@meta.data$celltype)){
  Idents(SC.data) <- "celltype"
  c_cells <- subset(SC.data, celltype == i)
  Idents(c_cells) <- "outcome"
  C_data <- FindMarkers(c_cells, ident.1 = "Deceased", ident.2 = "Cured")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_logFC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val)
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
#Save data
saveRDS(Gene.list,file="COVID19_genelist_seurat.rds")

#Get differential genes from Seurat DESeq2
SC.data=readRDS("COVID_SCdata_annotated.rds")
library('Seurat')
DefaultAssay(SC.data) <- "RNA"
set.seed(123456)
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.data@meta.data$celltype)){
  Idents(SC.data) <- "celltype"
  c_cells <- subset(SC.data, celltype == i)
  Idents(c_cells) <- "outcome"
  C_data <- FindMarkers(c_cells, ident.1 = "Deceased", ident.2 = "Cured", test.use = "DESeq2") #DESeq2 reports an error here that due to the sparsity of single-cell data
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_logFC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val)
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
#Save data
saveRDS(Gene.list,file="COVID19_genelist_DESeq2.rds")

#Get differential genes from EdgeR
library('edgeR')
SC.data=readRDS("COVID_SCdata_annotated.rds")
Case=Deceased.Severe
Control=Cured.Severe
DefaultAssay(SC.data) <- "RNA"
set.seed(123456)
min.cells=3 # The minimum number of cells for a cell type. A cell type is omitted if it has less cells than the minimum number.
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.data@meta.data$celltype)){
  Idents(SC.data) <- "celltype"
  c_cells <- subset(SC.data, celltype == i)
  Idents(c_cells) <- "outcome"
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
saveRDS(Gene.list,file="COVID19_genelist_edgeR.rds")

#Drug repurposing using FDA-approved drugs
my_gene_info<-read.table(file="DrugReference/lung_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="DrugReference/lung_drug_info.txt",sep="\t",header = T,quote = "")
cmap.ref.profiles = GetDrugRef(drug.response.path = 'DrugReference/lung_rankMatrix.txt',
                                 probe.to.genes = my_gene_info, drug.info = my_drug_info)
for(i in c("limma","edgeR","seurat")){
  Gene.list<-readRDS(paste0("COVID19_genelist_",i,".rds"))
  Drug.ident.res<-GetDrug(gene.data = Gene.list, drug.ref.profiles = cmap.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type = "FDA")
  saveRDS(Drug.ident.res,file=(paste0("COVID19_drugs_",i,"_FDA.rds")))
}
  
#Select FDA-approved mono-drugs
SC.data=readRDS('COVID_SCdata_annotated.rds')
Drug.ident.res<-readRDS("COVID19_drugs_limma_FDA.rds")
#Drug.ident.res<-readRDS("Windows/COVID19_FDA_drugs-windows-DC.rds")
Final.drugs<-TopDrug(SC.integrated=SC.data,
                     Drug.data=Drug.ident.res,
                     Drug.FDR=0.05,
                     FDA.drug.only=TRUE,
                     Case=Deceased.Severe
)

#Save data
saveRDS(Final.drugs,file="COVID19_selected_FDA_drugs.rds")

#Drug score
library(cmapR)
SC.data<-readRDS("COVID_SCdata_annotated.rds")
Gene.list<-readRDS("COVID19_genelist_limma.rds")
Drug.ident.res<-readRDS("COVID19_drugs_limma_FDA.rds")
GSE92742.gctx.path="GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
Tissue="lung"
Cell.type=c("Neutrophils","T_cells","Monocyte","NK_cell")
Deceased.Severe<-c("COVID19.C143","COVID19.C146","S-S086-1","S-S087-1")

for(i in c("limma","edgeR","seurat")){
  Gene.list<-readRDS(paste0("COVID19_genelist_",i,".rds"))
  Drug.ident.res<-readRDS(paste0("COVID19_drugs_",i,"_FDA.rds"))  
  Drug.score<-DrugScore(SC.integrated=SC.data,
                        Gene.data=Gene.list,
                        Cell.type=Cell.type,
                        Drug.data=Drug.ident.res,
                        FDA.drug.only=T,
                        Case=Deceased.Severe,
                        Tissue=Tissue,
                        GSE92742.gctx=GSE92742.gctx.path,
                        GSE70138.gctx=GSE70138.gctx.path)
  saveRDS(Drug.score,file=paste0("COVID19_drugscore_FDA_",i,".rds"))
}


#Drug repurposing using drugs/compounds
Gene.list<-readRDS("COVID19_genelist_limma.rds")
my_gene_info<-read.table(file="DrugReference/lung_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="DrugReference/lung_drug_info.txt",sep="\t",header = T,quote = "")
cmap.ref.profiles = GetDrugRef(drug.response.path = 'DrugReference/lung_rankMatrix.txt',
                               probe.to.genes = my_gene_info, drug.info = my_drug_info)
Drug.ident.res = GetDrug(gene.data = Gene.list, drug.ref.profiles = cmap.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type="all")
saveRDS(Drug.ident.res,file="COVID19_drugs_all.rds")

#Drug/compound score
Gene.list<-readRDS(paste0("COVID19_genelist_limma.rds"))
Drug.ident.res<-readRDS(paste0("COVID19_drugs_all.rds"))  
Drug.score<-DrugScore(SC.integrated=SC.data,
                      Gene.data=Gene.list,
                      Cell.type=Cell.type,
                      Drug.data=Drug.ident.res,
                      FDA.drug.only=T,
                      Case=Deceased.Severe,
                      Tissue=Tissue,
                      GSE92742.gctx=GSE92742.gctx.path,
                      GSE70138.gctx=GSE70138.gctx.path)
saveRDS(Drug.score,file=paste0("COVID19_all_drugscore_top4celltypes.rds"))
