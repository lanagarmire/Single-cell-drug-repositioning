#Download dataset GSE145926 from GEO before running this script.

library('Asgard')

setwd("Your_local_path/")

#Load data
data<-Read10X_h5(file="GSM4339769_C141_filtered_feature_bc_matrix.h5")
C141 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C141",type="Mild.COVID19",outcome="Cured"))
C141[["percent.mt"]] <- PercentageFeatureSet(C141, pattern = "^MT-")
C141 <- subset(C141, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4339770_C142_filtered_feature_bc_matrix.h5")
C142 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C142",type="Mild.COVID19",outcome="Cured"))
C142[["percent.mt"]] <- PercentageFeatureSet(C142, pattern = "^MT-")
C142 <- subset(C142, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4339771_C143_filtered_feature_bc_matrix.h5")
C143 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C143",type="Severe.COVID19",outcome="Death"))
C143[["percent.mt"]] <- PercentageFeatureSet(C143, pattern = "^MT-")
C143 <- subset(C143, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4339772_C144_filtered_feature_bc_matrix.h5")
C144 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C144",type="Mild.COVID19",outcome="Cured"))
C144[["percent.mt"]] <- PercentageFeatureSet(C144, pattern = "^MT-")
C144 <- subset(C144, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4339773_C145_filtered_feature_bc_matrix.h5")
C145 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C145",type="Severe.COVID19",outcome="Cured"))
C145[["percent.mt"]] <- PercentageFeatureSet(C145, pattern = "^MT-")
C145 <- subset(C145, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4339774_C146_filtered_feature_bc_matrix.h5")
C146 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C146",type="Severe.COVID19",outcome="Death"))
C146[["percent.mt"]] <- PercentageFeatureSet(C146, pattern = "^MT-")
C146 <- subset(C146, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4475048_C51_filtered_feature_bc_matrix.h5")
C51 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="Healthy.C51",type="Healthy",outcome="Healthy"))
C51[["percent.mt"]] <- PercentageFeatureSet(C51, pattern = "^MT-")
C51 <- subset(C51, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4475049_C52_filtered_feature_bc_matrix.h5")
C52 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="Healthy.C52",type="Healthy",outcome="Healthy"))
C52[["percent.mt"]] <- PercentageFeatureSet(C52, pattern = "^MT-")
C52 <- subset(C52, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4475050_C100_filtered_feature_bc_matrix.h5")
C100 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="Healthy.C100",type="Healthy",outcome="Healthy"))
C100[["percent.mt"]] <- PercentageFeatureSet(C100, pattern = "^MT-")
C100 <- subset(C100, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4475051_C148_filtered_feature_bc_matrix.h5")
C148 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C148",type="Severe.COVID19",outcome="Cured"))
C148[["percent.mt"]] <- PercentageFeatureSet(C148, pattern = "^MT-")
C148 <- subset(C148, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4475052_C149_filtered_feature_bc_matrix.h5")
C149 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C149",type="Severe.COVID19",outcome="Cured"))
C149[["percent.mt"]] <- PercentageFeatureSet(C149, pattern = "^MT-")
C149 <- subset(C149, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

data<-Read10X_h5(file="GSM4475053_C152_filtered_feature_bc_matrix.h5")
C152 <- CreateSeuratObject(counts = data, project = "COVID19", min.cells = 3, min.features = 200,meta.data=data.frame(row.names=colnames(data),cell=colnames(data),sample="COVID19.C152",type="Severe.COVID19",outcome="Cured"))
C152[["percent.mt"]] <- PercentageFeatureSet(C152, pattern = "^MT-")
C152 <- subset(C152, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

#Merge data
Mild.COVID19<-merge(x=C141,y=c(C142,C144))
Severe.COVID19<-merge(x=C143,y=c(C145,C146,C148,C149,C152))
Healthy<-merge(x=C51,y=c(C52,C100))

#single-cell alignment and annotation
SC.list<-list(Mild.COVID19=Mild.COVID19,Severe.COVID19=Severe.COVID19,Healthy=Healthy)
SC.data<-SCalignment(SC.list,by.CellType=TRUE)

#Save data
saveRDS(SC.data,file="COVID19_SCdata.rds")

#Creat the vector of sample names
Mild.COVID19<-c("COVID19.C141","COVID19.C142","COVID19.C144")
Severe.COVID19<-c("COVID19.C143","COVID19.C145","COVID19.C146","COVID19.C148","COVID19.C149","COVID19.C152")
Healthy<-c("Healthy.C51","Healthy.C52","Healthy.C100")

#Get differential genes
Gene.list<-GetGene(SC.integrated=SC.data,Case=Severe.COVID19,Control=Mild.COVID19,min.cells=3)

#Save data
saveRDS(Gene.list,file="COVID19_genelist.rds")

#Drug repurposing using FDA-approved drugs
Gene.list<-readRDS("COVID19_genelist.rds")
my_gene_info<-read.table(file="Your_loacal_path_for_drug_reference/lung_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="Your_loacal_path_for_drug_reference/lung_drug_info.txt",sep="\t",header = T,quote = "")
cmap.ref.profiles = GetDrugRef(drug.response.path = 'Your_loacal_path_for_drug_reference/lung_rankMatrix.txt',
                                 probe.to.genes = my_gene_info, drug.info = my_drug_info)
Drug.ident.res = GetDrug(gene.data = Gene.list, drug.ref.profiles = cmap.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type="FDA")
saveRDS(Drug.ident.res,file="COVID19_FDA_drugs.rds")

#Select FDA-approved mono-drugs
SC.data<-readRDS("COVID19_SCdata.rds")
Drug.ident.res<-readRDS("COVID19_FDA_drugs.rds")
Final.drugs<-TopDrug(SC.integrated=SC.data,
                     Drug.data=Drug.ident.res,
                     Drug.FDR=0.1,
                     FDA.drug.only=TRUE,
                     Case=Severe.COVID19
)

#Save data
saveRDS(Final.drugs,file="COVID19_selected_FDA_drugs.rds")

#Drug repurposing using compounds
Gene.list<-readRDS("COVID19_genelist.rds")
my_gene_info<-read.table(file="Your_loacal_path_for_drug_reference/lung_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="Your_loacal_path_for_drug_reference/lung_drug_info.txt",sep="\t",header = T,quote = "")
cmap.ref.profiles = GetDrugRef(drug.response.path = 'Your_loacal_path_for_drug_reference/lung_rankMatrix.txt',
                               probe.to.genes = my_gene_info, drug.info = my_drug_info)
Drug.ident.res = GetDrug(gene.data = Gene.list, drug.ref.profiles = cmap.ref.profiles, repurposing.unit = "drug", connectivity = "negative", drug.type="compounds")
saveRDS(Drug.ident.res,file="COVID19_drugs.rds")

#Select mono-compounds
SC.data<-readRDS("COVID19_SCdata.rds")
Drug.ident.res<-readRDS("COVID19_compounds.rds")
Final.drugs<-TopDrug(SC.integrated=SC.data,
                     Drug.data=Drug.ident.res,
                     Drug.FDR=0.1,
                     FDA.drug.only=FALSE,
                     Case=Severe.COVID19
)

#Save data
saveRDS(Final.drugs,file="COVID19_selected_compounds.rds")
