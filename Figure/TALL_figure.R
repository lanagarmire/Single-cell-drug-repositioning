#Download datasets GSE92742 and GSE70138 from GEO before running this script.
library('Asgard')
library('Seurat')

setwd("./")

#Load data
SC.data<-readRDS("TALL_SCdata_all.rds")
#View alignment results
pdf(file = "Fig5A.pdf",width = 6,height = 3)
DimPlot(SC.data, reduction = "umap", split.by = "type",group.by = "celltype")
dev.off()

#Load data
SC.data<-readRDS("TALL_SCdata_Tcell.rds")
SC.data@meta.data$type[which(SC.data@meta.data$type=="Normal")]="Normal T cells"
#View alignment results
pdf(file = "Fig5B.pdf",width = 6,height = 3)
DimPlot(SC.data, reduction = "umap", split.by = "type",group.by = "celltype")
dev.off()


#Pathway analysis
GL<-readRDS("TALL_genelist_limma.rds")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)  
Pathway.table=data.frame()
for(i in paste0("C",c(1:4))){
  Gene.table<-data.frame()
  limma.data<-GL[[i]]
  limma.data<-subset(limma.data,P.Value<0.05 & abs(score)>1)
  limma.data<-limma.data %>% top_n(n = 1000, wt = abs(score))
  temp.data<-data.frame(Celltype=i,Gene=row.names(limma.data),Value=limma.data$score)
  Gene.table<-rbind(Gene.table,temp.data)
  row.names(Gene.table)=Gene.table$Gene
  #Pathway
  original_gene_list <- Gene.table$Value
  names(original_gene_list) <- Gene.table$Gene
  gene_list<-na.omit(original_gene_list)
  ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
  dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
  df2 = Gene.table[row.names(Gene.table) %in% dedup_ids$SYMBOL,]
  df2$Y = dedup_ids$ENTREZID
  kegg_gene_list <- df2$Value
  names(kegg_gene_list) <- df2$Y
  kegg_gene_list<-na.omit(kegg_gene_list)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  GSE <- enrichKEGG(gene = df2$Y,
                    pvalueCutoff = 0.05,
                    organism = 'hsa')
  Pathway=subset(GSE@result,pvalue<0.05)
  sig.pos=grep("signaling pathway",Pathway$Description)
  a.pos=grep("Apoptosis",Pathway$Description)
  c.pos=grep("Cell cycle",Pathway$Description)
  rm.pos1=grep("B cell",Pathway$Description)
  sig.pos=setdiff(sig.pos,c(rm.pos1))
  Pathway=Pathway[c(sig.pos,a.pos,c.pos),]
  Temp=data.frame(Pathway=Pathway$Description,Adjusted.Pvalue=Pathway$p.adjust,Cluster=i)
  Pathway.table=rbind(Pathway.table,Temp)
}
P.s=table(Pathway.table$Pathway)
Pnames=names(P.s[which(P.s>1)])
Pathway.table=subset(Pathway.table,Pathway%in%Pnames)
pdf(file = "Fig5_marker_bubble.pdf",width = 5,height = 3)
ggplot(Pathway.table,aes(x=Pathway, y=Cluster, size=-log10(Adjusted.Pvalue),colour = 'red')) +
  geom_point(alpha=1) +
  scale_size(range = c(1, 5), name="-log10(Adjusted.Pvalue)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ guides(color=FALSE)+ theme(legend.position="top")
dev.off()

#Drug score plot
SC.data<-readRDS("TALL_SCdata_Tcell.rds")
Gene.list<-readRDS("TALL_genelist_limma.rds")
Drug.ident.res<-readRDS("TALL_drugs_limma.rds")
Case=c("PRE_T_1","PRE_T_2")
Tissue="haematopoietic and lymphoid tissue"
GSE92742.gctx.path="GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
Combined.Drug.score<-DrugScore(SC.integrated=SC.data,
                      Gene.data=Gene.list,
                      Cell.type=NULL,
                      Drug.data=Drug.ident.res,
                      FDA.drug.only=T,
                      Case=Case,
                      Tissue=Tissue,
                      GSE92742.gctx=GSE92742.gctx.path,
                      GSE70138.gctx=GSE70138.gctx.path)
Score.list<-data.frame(Patient="Overall",Drug=row.names(Combined.Drug.score),DrugScore=Combined.Drug.score$Drug.therapeutic.score,Pvalue=Combined.Drug.score$P.value,FDR=Combined.Drug.score$FDR)

j=0
for(i in Case){
  j=j+1
  Drug.score<-DrugScore(SC.integrated=SC.data,
                        Gene.data=Gene.list,
                        Cell.type=NULL,
                        Drug.data=Drug.ident.res,
                        FDA.drug.only=T,
                        Case=i,
                        Tissue=Tissue,
                        GSE92742.gctx=GSE92742.gctx.path,
                        GSE70138.gctx=GSE70138.gctx.path)
  Temp=data.frame(Patient=paste0("Patient",j),Drug=row.names(Drug.score),DrugScore=Drug.score$Drug.therapeutic.score,Pvalue=Drug.score$P.value,FDR=Drug.score$FDR)
  Score.list=rbind(Score.list,Temp)
}
library(Hmisc)
Score.list$Drug<-capitalize(Score.list$Drug)
sig.list<-subset(Score.list,DrugScore>quantile(Score.list$DrugScore, 0.9,na.rm=T))
Drug.order<-row.names(Combined.Drug.score)[order(Combined.Drug.score$Drug.therapeutic.score,decreasing = F)]
Drug.order<-capitalize(Drug.order)
sig.list$Drug<-factor(sig.list$Drug,levels = Drug.order)
library('ggplot2')
pdf(file = "Drug_individual_TALL.pdf",width = 3,height = 3)
ggplot(sig.list,aes(x=Patient, y=Drug, size=DrugScore, color=Patient)) +
  geom_point(alpha=1) +
  scale_size(name="DrugScore")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+ guides(color=FALSE)
dev.off()

