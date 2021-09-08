#Download dataset GSE92742 and GSE70138 from GEO before running this script.

library('Asgard')
library('Seurat')

setwd("./")

#View alignment results
SC.data=readRDS('COVID_SCdata_annotated.rds')
pdf(file = "COVID19_immune_cc_removed_S_celltype.pdf",width = 7,height = 3.5)
DimPlot(SC.data, reduction = "umap", split.by = "outcome",group.by ="celltype",label = T,pt.size=1)+NoLegend()
dev.off()

#Cell type changes
library('ggplot2')
library('cowplot')
data<-as.data.frame(SC.data@meta.data)
celltypes<-unique(data$celltype)
final.table<-data.frame()
for(i in c("Cured","Deceased")){
  sub.data<-subset(data,outcome==i)
  severity<-unique(sub.data$type)
  celltype.freq<-round(100*table(sub.data$celltype)/nrow(sub.data),2)
  celltype.freq<-celltype.freq[celltypes]
  final.table<-rbind(final.table,celltype.freq)
}
colnames(final.table)<-celltypes
row.names(final.table)<-c("Cured","Deceased")

##Proportion plot
plot.data<-data.frame()
for(i in 1:length(colnames(final.table))){
  plot.temp<-data.frame(Cell.Type=colnames(final.table)[i], Cluster.Size=final.table[,i], Group=row.names(final.table))
  plot.data<-rbind(plot.data,plot.temp)
}
plot.data$Group = factor(plot.data$Group, levels=c("Deceased","Cured"))
plot.data$Cell.Type = factor(plot.data$Cell.Type, levels=c("Monocyte","Macrophage","T_cells","Neutrophils","Epithelial_cells","B_cell","NK_cell"))
p1<-ggplot(plot.data, aes(Cell.Type, Cluster.Size, fill=Group)) +
  labs(x="", y = "Proportion of cells")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1))
new.data<-log10(final.table["Deceased",]/final.table["Cured",])
plot.data<-data.frame()
for(i in 1:length(colnames(final.table))){
  plot.temp<-data.frame(Cell.Type=colnames(new.data)[i], Value=new.data[,i])
  plot.data<-rbind(plot.data,plot.temp)
}
plot.data$Cell.Type = factor(plot.data$Cell.Type, levels=unique(plot.data$Cell.Type[order(plot.data$Value,decreasing = T)]))
p2<-ggplot(plot.data, aes(Cell.Type,Value, fill=Value)) +
  labs(x="", y = "Fold change of cluster size\n(Deceased VS Cured, Log10)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1))
pdf(file = "Celltype_change.pdf",width = 7,height = 3.5)
plot_grid(p1,p2,ncol=2)
dev.off()

#Pathway analysis
GL<-readRDS("COVID19_genelist_limma.rds")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(dplyr)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)  
Pathway.table=data.frame()
for(i in c("Neutrophils","T_cells","Monocyte","NK_cell")){
  Gene.table<-data.frame()
  limma.data<-GL[[i]]
  limma.data<-subset(limma.data,adj.P.Val<0.05& abs(score)>1)
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
  Pathway=subset(GSE@result,p.adjust<0.05)
  key.pathways=c("Coronavirus disease - COVID-19",
                 "Toll-like receptor signaling pathway",
                 "NF-kappa B signaling pathway",
                 "Chemokine signaling pathway",
                 "TNF signaling pathway",
                 "IL-17 signaling pathway",
                 "JAK-STAT signaling pathway",
                 "T cell receptor signaling pathway")
  Pathway=subset(Pathway,Description %in% key.pathways)
  Temp=data.frame(Pathway=Pathway$Description,Adjusted.Pvalue=Pathway$p.adjust,Cluster=i)
  Pathway.table=rbind(Pathway.table,Temp)
}
Pathway.table$Cluster = factor(Pathway.table$Cluster, levels=c("Monocyte","T_cells","Neutrophils","NK_cell"))
pdf(file = "Fig6_marker_bubble.pdf",width = 7,height = 3.5)
ggplot(Pathway.table,aes(x=Pathway, y=Cluster, size=-log10(Adjusted.Pvalue),colour = 'red')) +
  geom_point(alpha=1) + guides(color=FALSE)+ theme(legend.position="top")+
  scale_size(range = c(1, 5), name="-log10(Adjusted.Pvalue)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1))
dev.off()

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
Drug.score<-DrugScore(SC.integrated=SC.data,
                      Gene.data=Gene.list,
                      Cell.type=Cell.type,
                      Drug.data=Drug.ident.res,
                      FDA.drug.only=T,
                      Case=Deceased.Severe,
                      Tissue=Tissue,
                      GSE92742.gctx=GSE92742.gctx.path,
                      GSE70138.gctx=GSE70138.gctx.path)
saveRDS(Drug.score,file="COVID19_FDA_drugscore_top4celltypes.rds")
Score.list<-data.frame(Patient="Overall",Drug=row.names(Drug.score),DrugScore=Drug.score$Drug.therapeutic.score,Pvalue=Drug.score$P.value,FDR=Drug.score$FDR)

Combined.Drug.score=Drug.score

j=0
for(i in Deceased.Severe){
  j=j+1
  Drug.score<-DrugScore(SC.integrated=SC.data,
                        Gene.data=Gene.list,
                        Cell.type=Cell.type,
                        Drug.data=Drug.ident.res,
                        FDA.drug.only=T,
                        Case=i,
                        Tissue=Tissue,
                        GSE92742.gctx=GSE92742.gctx.path,
                        GSE70138.gctx=GSE70138.gctx.path)
  Temp=data.frame(Patient=paste0("Patient",j),Drug=row.names(Drug.score),DrugScore=Drug.score$Drug.therapeutic.score,Pvalue=Drug.score$P.value,FDR=Drug.score$FDR)
  Score.list=rbind(Score.list,Temp)
}
Score.list$Drug<-capitalize(Score.list$Drug)
sig.list<-subset(Score.list,FDR<0.05 & DrugScore>quantile(Score.list$DrugScore, 0.9,na.rm=T))
Drug.order<-row.names(Combined.Drug.score)[order(Combined.Drug.score$Drug.therapeutic.score,decreasing = T)]
library(R.utils)
Drug.order<-capitalize(Drug.order)
sig.list$Drug<-factor(sig.list$Drug,levels = Drug.order)
library('ggplot2')
pdf(file = "Drug_individual.pdf",width = 7,height = 3)
ggplot(sig.list,aes(x=Drug, y=Patient, size=DrugScore, color=Patient)) +
  geom_point(alpha=1) +
  scale_size(name="DrugScore")+ guides(color=FALSE)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1))
dev.off()

