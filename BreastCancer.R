#Download datasets GSE113197 and GSE123926 from GEO before running this script.

library('Asgard')
library('Seurat')
setwd("Your_local_path/")

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
SC.data<-SCalignment(SC.list)

#Change sample names
sample<-SC.data@meta.data$sample
sample[which(sample=="Ind5")]<-"Normal1"
sample[which(sample=="Ind6")]<-"Normal2"
sample[which(sample=="Ind7")]<-"Normal3"
SC.data@meta.data$sample<-sample

#View alignment results
DimPlot(SC.data, reduction = "umap", split.by = "sample", label = TRUE)

#Save data
saveRDS(SC.data,file="TNBC_SCdata.rds")

#Get differential genes
Case=c("PDX-110","PDX-332")
Control=c("Normal1","Normal2","Normal3")
Gene.list<-GetGene(SC.integrated=SC.data,Case=Case,Control=Control,min.cells=3)
#Save data
saveRDS(Gene.list,file="TNBC_genelist.rds")


#Drug repurposing
Gene.list<-readRDS("TNBC_genelist.rds")
my_gene_info<-read.table(file="Your_loacal_path_for_drug_reference/breast_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="Your_loacal_path_for_drug_reference/breast_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = 'Your_loacal_path_for_drug_reference/breast_rankMatrix.txt',
                                 probe.to.genes = my_gene_info, drug.info = my_drug_info)
Drug.ident.res = GetDrug(gene.data = Gene.list, drug.ref.profiles = drug.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type = "all")
saveRDS(Drug.ident.res,file="TNBC_drugs.rds")

#Drug combination
SC.data<-readRDS("TNBC_SCdata.rds")
Gene.list<-readRDS("TNBC_genelist.rds")
Drug.ident.res<-readRDS("TNBC_drugs.rds")
GSE92742.gctx.path="Your_local_path/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="Your_local_path/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
Tissue="breast"
Drug.combinations<-DrugCombination(SC.integrated=SC.data,
                      Gene.data=Gene.list,
                      Drug.data=Drug.ident.res,
                      Drug.FDR=0.1,
                      FDA.drug.only=TRUE,
                      Combined.drugs=2,
                      Case=Case,
                      Tissue="breast",
                      GSE92742.gctx=GSE92742.gctx.path,
                      GSE70138.gctx=GSE70138.gctx.path)
saveRDS(Drug.combinations,file="TNBC_drugcombinations.rds")

#Select mono-drugs
Final.drugs<-TopDrug(SC.integrated=SC.data,
                   Drug.data=Drug.ident.res,
                   Drug.FDR=0.1,
                   FDA.drug.only=TRUE,
                   Case=Case
)
saveRDS(Final.drugs,file="TNBC_selected_drugs.rds")
#Select drug combinations
Final.combinations<-TopCombination(Drug.combination=Drug.combinations,
                   Combination.FDR=0.1,
                   Min.combination.score=1
)
saveRDS(Final.combinations,file="TNBC_selected_drugcombinations.rds")
