#Download datasets GSE92742 and GSE70138 from GEO before running this script.
library('Asgard')
library('Seurat')

setwd("./")

#Load data
SC.data<-readRDS("TNBC_SCdata.rds")
#View alignment results
pdf(file = "Fig4A.pdf",width = 6,height = 3)
DimPlot(SC.data, reduction = "umap", split.by = "type",group.by = "celltype")
dev.off()

#Pathway analysis
GL<-readRDS("TNBC_genelist_limma.rds")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)  
Pathway.table=data.frame()
for(i in paste0("C",c(1:8))){
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
  sig.pos=grep("signaling pathway",Pathway$Description)
  a.pos=grep("Apoptosis",Pathway$Description)
  c.pos=grep("Cell cycle",Pathway$Description)
  Pathway=Pathway[c(sig.pos,a.pos,c.pos),]
  Temp=data.frame(Pathway=Pathway$Description,Adjusted.Pvalue=Pathway$p.adjust,Cluster=i)
  Pathway.table=rbind(Pathway.table,Temp)
}
pdf(file = "Fig4_marker_bubble.pdf",width = 6,height = 3)
ggplot(Pathway.table,aes(x=Cluster, y=Pathway, size=-log10(Adjusted.Pvalue),colour = 'red')) +
  geom_point(alpha=1) + guides(color=FALSE)+ theme(legend.position="top")+
  scale_size(range = c(1, 5), name="-log10(Adjusted.Pvalue)")
dev.off()

#Drug score plot
Gene.list<-readRDS("TNBC_genelist_limma.rds")
Drug.ident.res<-readRDS("TNBC_drugs_limma.rds")
Case<-c("PDX-110","PDX-332")
Tissue="breast"
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
  Temp=data.frame(Patient=i,Drug=row.names(Drug.score),DrugScore=Drug.score$Drug.therapeutic.score,Pvalue=Drug.score$P.value,FDR=Drug.score$FDR)
  Score.list=rbind(Score.list,Temp)
}
library(Hmisc)
Score.list$Drug<-capitalize(Score.list$Drug)
sig.list<-subset(Score.list,FDR<0.05 & DrugScore>quantile(Score.list$DrugScore, 0.99))
Drug.order<-row.names(Combined.Drug.score)[order(Combined.Drug.score$Drug.therapeutic.score,decreasing = F)]
Drug.order<-capitalize(Drug.order)
sig.list$Drug<-factor(sig.list$Drug,levels = Drug.order)
library('ggplot2')
pdf(file = "Drug_individual_TNBC.pdf",width = 3.5,height = 5)
ggplot(sig.list,aes(x=Patient, y=Drug, size=DrugScore, color=Patient)) +
  geom_point(alpha=1) +
  scale_size(name="DrugScore")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0))
dev.off()

#Drug combination score plot
SC.data<-readRDS("TNBC_SCdata.rds")
Gene.list<-readRDS("TNBC_genelist_limma.rds")
Drug.ident.res<-readRDS("TNBC_drugs_limma.rds")
Case<-c("PDX-110","PDX-332")
Tissue="breast"
GSE92742.gctx.path="GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"
Combined.Drug.score<-DrugCombination(SC.integrated=SC.data,
                                     Gene.data=Gene.list,
                                     Drug.data=Drug.ident.res,
                                     Drug.FDR=0.05,
                                     FDA.drug.only=TRUE,
                                     Combined.drugs=2,
                                     Case=Case,
                                     Tissue=Tissue,
                                     GSE92742.gctx=GSE92742.gctx.path,
                                     GSE70138.gctx=GSE70138.gctx.path)
Score.list<-data.frame(Patient="Overall",
                       Drug=paste0(capitalize(Combined.Drug.score$Drug1),"+",capitalize(Combined.Drug.score$Drug2)),
                       DrugScore=Combined.Drug.score$Combination.therapeutic.score,
                       Pvalue=Combined.Drug.score$P.value,
                       FDR=Combined.Drug.score$FDR)

for(i in Case){
  Drug.score<-DrugCombination(SC.integrated=SC.data,
                              Gene.data=Gene.list,
                              Drug.data=Drug.ident.res,
                              Drug.FDR=0.05,
                              FDA.drug.only=TRUE,
                              Combined.drugs=2,
                              Case=i,
                              Tissue=Tissue,
                              GSE92742.gctx=GSE92742.gctx.path,
                              GSE70138.gctx=GSE70138.gctx.path)
  Temp=data.frame(Patient=i,
                  Drug=paste0(capitalize(Drug.score$Drug1),"+",capitalize(Drug.score$Drug2)),
                  DrugScore=Drug.score$Combination.therapeutic.score,
                  Pvalue=Drug.score$P.value,
                  FDR=Drug.score$FDR)
  Score.list=rbind(Score.list,Temp)
}
sig.list<-subset(Score.list,FDR<0.05 & DrugScore>quantile(Score.list$DrugScore, 0.8))
Drug.order<-sig.list$Drug[order(sig.list$DrugScore,decreasing = F)]
Drug.order<-unique(Drug.order)
sig.list$Drug<-factor(sig.list$Drug,levels = Drug.order)
library('ggplot2')
pdf(file = "Drug_combination_individual_TNBC.pdf",width = 4,height = 4)
ggplot(sig.list,aes(x=Patient, y=Drug, size=DrugScore, color=Patient)) +
  geom_point(alpha=1) +
  scale_size(name="DrugScore")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0))
dev.off()
