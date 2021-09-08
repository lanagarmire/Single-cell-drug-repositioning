#Download datasets GSE132509 from GEO before running this script.

library('Asgard')
library('Seurat')

setwd("./")

#Load dataset
data <- Read10X(data.dir = "PBMMC_1")
PBMMC_1 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Normal",patient="PBMMC_1",sample="PBMMC_1"))
data <- Read10X(data.dir = "PBMMC_2")
PBMMC_2 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Normal",patient="PBMMC_2",sample="PBMMC_2"))
data <- Read10X(data.dir = "PBMMC_3")
PBMMC_3 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Normal",patient="PBMMC_3",sample="PBMMC_3"))
data <- Read10X(data.dir = "ETV6_RUNX1_1")
ETV6_RUNX1_1 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Pre-B ETV6-RUNX1 ALL",patient="ETV6_RUNX1_1",sample="ETV6_RUNX1_1"))
data <- Read10X(data.dir = "ETV6_RUNX1_2")
ETV6_RUNX1_2 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Pre-B ETV6-RUNX1 ALL",patient="ETV6_RUNX1_2",sample="ETV6_RUNX1_2"))
data <- Read10X(data.dir = "ETV6_RUNX1_3")
ETV6_RUNX1_3 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Pre-B ETV6-RUNX1 ALL",patient="ETV6_RUNX1_3",sample="ETV6_RUNX1_3"))
data <- Read10X(data.dir = "ETV6_RUNX1_4")
ETV6_RUNX1_4 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Pre-B ETV6-RUNX1 ALL",patient="ETV6_RUNX1_4",sample="ETV6_RUNX1_4"))
data <- Read10X(data.dir = "HHD_1")
HHD_1 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Pre-B HHD ALL",patient="HHD_1",sample="HHD_1"))
data <- Read10X(data.dir = "HHD_2")
HHD_2 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Pre-B HHD ALL",patient="HHD_2",sample="HHD_2"))
data <- Read10X(data.dir = "PRE_T_1")
PRE_T_1 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Pre-T ALL",patient="PRE_T_1",sample="PRE_T_1"))
data <- Read10X(data.dir = "PRE_T_2")
PRE_T_2 <- CreateSeuratObject(counts = data, project = "Leukemia", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),type="Pre-T ALL",patient="PRE_T_2",sample="PRE_T_2"))


#single-cell alignment
case.type="Pre-T ALL"
Case=c("PRE_T_1","PRE_T_2")
target.cell="T_cells"
Control=c("PBMMC_1","PBMMC_2","PBMMC_3")

SC.list<-list(PRE_T_1=PRE_T_1,PRE_T_2=PRE_T_2,PBMMC_1=PBMMC_1,PBMMC_2=PBMMC_2,PBMMC_3=PBMMC_3)
CellCycle=TRUE #Set it TRUE if you want to do Cell Cycle Regression
anchor.features=2000 
by.CellType=TRUE #Set it TRUE if you want to annotate cell  type
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
#View alignment results
pdf(file = "TALL_plot_all.pdf",width = 8,height = 3.5)
DimPlot(SC.integrated, reduction = "umap", split.by = "type",group.by = "celltype",label = T)
dev.off()
#Save data
saveRDS(SC.integrated,file="TALL_SCdata_all.rds")

#SC.integrated=readRDS("TALL_SCdata_all.rds")
#Select target cells
SC.integrated<-subset(SC.integrated,celltype %in% target.cell)
temp.table=data.frame(row.names=unique(paste0("C",as.numeric(SC.integrated@meta.data$seurat_clusters))),new=paste0("C",1:length(unique(SC.integrated@meta.data$seurat_clusters))))
new.cluster=temp.table[paste0("C",as.numeric(SC.integrated@meta.data$seurat_clusters)),]
SC.integrated@meta.data$celltype <- new.cluster

#View alignment results
pdf(file = "TALL_plot_Tcell.pdf",width = 8,height = 3.5)
DimPlot(SC.integrated, reduction = "umap", split.by = "type",group.by = "celltype",label = T)
dev.off()

#Save data
saveRDS(SC.integrated,file="TALL_SCdata_Tcell.rds")


#Cell type changes
library('ggplot2')
library('cowplot')
data<-as.data.frame(SC.integrated@meta.data)
celltypes<-unique(data$celltype)
final.table<-data.frame()
for(i in c("Normal",case.type)){
  sub.data<-subset(data,type==i)
  severity<-unique(sub.data$type)
  celltype.freq<-round(100*table(sub.data$celltype)/nrow(sub.data),2)
  celltype.freq<-celltype.freq[celltypes]
  final.table<-rbind(final.table,celltype.freq)
}
colnames(final.table)<-celltypes
row.names(final.table)<-c("Normal",case.type)

#Get differential genes from limma
SC.integrated=readRDS("TALL_SCdata_Tcell.rds")
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
    expr <- c_cells@assays$RNA@data
    new_expr <- expr[,c(Casesample,Controlsample)]
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
saveRDS(Gene.list,file="TALL_genelist_limma.rds")

#Get differential genes from Seurat
SC.integrated=readRDS("TALL_SCdata_Tcell.rds")
library('Seurat')
DefaultAssay(SC.integrated) <- "RNA"
set.seed(123456)
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.integrated@meta.data$celltype)){
  Idents(SC.integrated) <- "celltype"
  c_cells <- subset(SC.integrated, celltype == i)
  Idents(c_cells) <- "type"
  C_data <- FindMarkers(c_cells, ident.1 = case.type, ident.2 = "Normal")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_logFC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val)
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
#Save data
saveRDS(Gene.list,file="TALL_genelist_seurat.rds")

#Get differential genes from Seurat DESeq2
SC.integrated=readRDS("TALL_SCdata_Tcell.rds")
library('Seurat')
DefaultAssay(SC.integrated) <- "RNA"
set.seed(123456)
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.integrated@meta.data$celltype)){
  Idents(SC.integrated) <- "celltype"
  c_cells <- subset(SC.integrated, celltype == i)
  Idents(c_cells) <- "type"
  C_data <- FindMarkers(c_cells, ident.1 = case.type, ident.2 = "Normal", test.use = "DESeq2")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_logFC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val)
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
#Save data
saveRDS(Gene.list,file="TALL_genelist_DESeq2.rds")

#Get differential genes from EdgeR
library('edgeR')
SC.integrated=readRDS("TALL_SCdata_Tcell.rds")
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
saveRDS(Gene.list,file="TALL_genelist_edgeR.rds")

#Drug repurposing

for(i in c("limma","edgeR","DESeq2","seurat")){
Gene.list<-readRDS(paste0("TALL_genelist_",i,".rds"))
my_gene_info<-read.table(file="DrugReference/haematopoietic-and-lymphoid-tissue_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="DrugReference/haematopoietic-and-lymphoid-tissue_drug_info.txt",sep="\t",header = T,quote = "")
cmap.ref.profiles<-GetDrugRef(drug.response.path = 'DrugReference/haematopoietic-and-lymphoid-tissue_rankMatrix.txt',
                                 probe.to.genes = my_gene_info, drug.info = my_drug_info)
Drug.ident.res<-GetDrug(gene.data = Gene.list, drug.ref.profiles = cmap.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type = "FDA")
saveRDS(Drug.ident.res,file=(paste0("TALL_drugs_FDA_",i,".rds")))
}

#Drug score
GSE92742.gctx.path="GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
Tissue="haematopoietic and lymphoid tissue"
for(i in c("limma","edgeR","DESeq2","seurat")){
  Gene.list<-readRDS(paste0("TALL_genelist_",i,".rds"))
  Drug.ident.res<-readRDS(paste0("TALL_drugs_FDA_",i,".rds"))  
  Drug.score<-DrugScore(SC.integrated=SC.integrated,
                     Gene.data=Gene.list,
                     Cell.type=NULL,
                     Drug.data=Drug.ident.res,
                     FDA.drug.only=T,
                     Case=Case,
                     Tissue=Tissue,
                     GSE92742.gctx=GSE92742.gctx.path,
                     GSE70138.gctx=GSE70138.gctx.path)
   saveRDS(Drug.score,file=paste0("TALL_drugscore_FDA_",i,".rds"))
}

#Drug repurposing using drugs/compounds
Gene.list<-readRDS("TALL_genelist_limma.rds")
my_gene_info<-read.table(file="DrugReference/haematopoietic-and-lymphoid-tissue_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="DrugReference/haematopoietic-and-lymphoid-tissue_drug_info.txt",sep="\t",header = T,quote = "")
cmap.ref.profiles<-GetDrugRef(drug.response.path = 'DrugReference/haematopoietic-and-lymphoid-tissue_rankMatrix.txt',
                              probe.to.genes = my_gene_info, drug.info = my_drug_info)
Drug.ident.res = GetDrug(gene.data = Gene.list, drug.ref.profiles = cmap.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type="all")
saveRDS(Drug.ident.res,file="TALL_drugs_limma_all.rds")

#Drug/compound score
Gene.list<-readRDS(paste0("TALL_genelist_limma.rds"))
Drug.ident.res<-readRDS(paste0("TALL_drugs_limma_all.rds"))  
Drug.score<-DrugScore(SC.integrated=SC.integrated,
                      Gene.data=Gene.list,
                      Cell.type=NULL,
                      Drug.data=Drug.ident.res,
                      FDA.drug.only=T,
                      Case=Case,
                      Tissue=Tissue,
                      GSE92742.gctx=GSE92742.gctx.path,
                      GSE70138.gctx=GSE70138.gctx.path)
saveRDS(Drug.score,file=paste0("TALL_drugscore_all.rds"))


