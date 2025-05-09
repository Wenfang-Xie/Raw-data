library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)

dir <- dir("./")
samples_name = c("GSM6622299", "GSM6622300", "GSM6622301")
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- readRDS(dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}
names(scRNAlist) <- samples_name
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
##############################################################################################################

scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<10000&percent.mito<10)
scRNA = SCTransform(scRNA, vars.to.regress="percent.mito", verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
colnames(scRNA@meta.data)[1] = "Sample"

scRNA = RunHarmony(scRNA, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
###################
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:30, reduction="harmony")
colors = c("#8BE0A4", "#C9DB74", "#FE88B1", "#B497E7", "#DCB0F2", "#87C55F", "#D3B484", "#F6CF71", "#9EB9F3", "#F89C74", "#66C5CC","midnightblue","yellowgreen")
###############################################################
custome_theme_1=theme(axis.text.y=element_text(family="Times",face="plain")
                      ,axis.text.x=element_text(family="Times",face="plain")
                      ,plot.title = element_text(hjust = 0.5,family="Times",face="plain")
                      ,axis.title.x=element_text(family="Times",face="plain")
                      ,axis.title.y=element_text(family="Times",face="plain")
                      ,legend.title = element_text(family="Times",face="plain")
                      ,legend.text = element_text(family="Times",face="plain"))

sc_umap = DimPlot(scRNA,cols=colors,group.by='Sample',
                  reduction="umap",
                  label = F, 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
sc_umap
# LabelClusters(sc_umap,id = 'cell_type',family = 'Times',size = 6,fontface = 'bold',color = 'red')
sc_umap_Sample = DimPlot(scRNA,cols=colors,group.by='Sample',
                         reduction="umap",
                         label = "F", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

sc_umap_Sample
p1=sc_umap+custome_theme_1
p2=sc_umap_Sample+custome_theme_1
p1=ggarrange(AugmentPlot(p1,dpi = 300), ncol = 1, nrow = 1,
             legend ='right',legend.grob=get_legend(p1),
             font.label = list(size = 12, face = "bold",family ='Times'))
ggsave(filename = 'sample_umap.pdf',plot = p1,he=8,width = 9)
#########################################################################
p = VlnPlot(scRNA, features=c("nFeature_RNA"), cols=colors, pt.size=0)+NoLegend()+theme(text=element_text(family="Times"))
ggsave("./violin_nFeature_RNA.pdf", p, width=6, height=6)
p = VlnPlot(scRNA, features=c("percent.mito"), cols=colors, pt.size=0)+NoLegend()+theme(text=element_text(family="Times"))
ggsave("./violin_percent_mito.pdf", p, width=6, height=6)

scRNA <- FindClusters(scRNA, resolution = seq(0.1,0.5,by=0.05))
clustree(scRNA)
mydata <- FindClusters(scRNA, resolution=0.4)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25)
write.table(markers, "subcluster_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")

cell_label = c("CD8+T Cell","Melanocyte","Th Cell","Glial cell","Melanocyte","Melanocyte","B Cell","Endothelial cell",
               "Glial cell","Fibroblast","Smooth muscle Cell","Glial cell","Monocyte")
#C0:CD8+T Cell:CD8A	CD8B
#C1,4,5:Melanocyte:MLANA	MITF	PMEL
#C2:Th Cell:IL32	IL7R	CD69
#C3,8,11:Glial cell :S100A6	S100B	S100A10
#C6:B Cell:IGHG1	IGHA1
#C7:Endothelial cell:PECAM1	EGFL7
#C9:Fibroblast:LUM DCN
#C10:Smooth muscle Cell:TAGLN ACTA2
#C12:Monocyte:LYZ	S100A9	S100A8
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
mydata = subset(mydata, cell_type %in% c("CD8+T Cell","Melanocyte","Th Cell","Glial cell","Melanocyte","Melanocyte","B Cell","Endothelial cell",
                                         "Glial cell","Fibroblast","Smooth muscle Cell","Glial cell","Monocyte"))
###############################################################
custome_theme_1=theme(axis.text.y=element_text(family="Times",face="plain")
                      ,axis.text.x=element_text(family="Times",face="plain")
                      ,plot.title = element_text(hjust = 0.5,family="Times",face="plain")
                      ,axis.title.x=element_text(family="Times",face="plain")
                      ,axis.title.y=element_text(family="Times",face="plain")
                      ,legend.title = element_text(family="Times",face="plain")
                      ,legend.text = element_text(family="Times",face="plain"))
sc_umap = DimPlot(mydata,cols=colors,group.by='cell_type',
                  reduction="umap",
                  label = F, 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
sc_umap
# LabelClusters(sc_umap,id = 'cell_type',family = 'Times',size = 6,fontface = 'bold',color = 'red')
sc_umap_Sample = DimPlot(mydata,cols=colors,group.by='cell_type',
                         reduction="umap",
                         label = "F", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

sc_umap_Sample
p1=sc_umap+custome_theme_1
p2=sc_umap_Sample+custome_theme_1
p1=ggarrange(AugmentPlot(p1,dpi = 300), ncol = 1, nrow = 1,
             legend ='right',legend.grob=get_legend(p1),
             font.label = list(size = 12, face = "bold",family ='Times'))
ggsave(filename = 'cell_umap.pdf',plot = p1,he=8,width = 9)
#########################################################################
saveRDS(mydata, "scRNA_raw.rds")
genes = c("CD8A","CD8B","MLANA","MITF","PMEL","IL32","IL7R","CD69","S100A6","S100B",
          "S100A10","IGHG1","IGHA1","PECAM1","EGFL7","LUM","DCN","TAGLN","ACTA2","LYZ","S100A9","S100A8")
p = DotPlot(mydata, features=genes)+coord_flip()+scale_color_distiller(palette="RdYlBu")+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("marker_dotplot.pdf", p, width=7, height=6)

p = DoHeatmap(mydata, features=genes, group.colors=colors, label=FALSE)+scale_fill_gradientn(colors=c("white", "snow", "firebrick3"))+theme(text=element_text(family="Times"))
ggsave("marker_heatmap.pdf", p, width=10, height=5)

gene2 = c("GZMA", "GSDMB")
p = DotPlot(mydata, features=gene2)+coord_flip()+scale_color_distiller(palette="RdYlBu")+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("pyroptosis_dotplot.pdf", p, width=7, height=6)

p = DoHeatmap(mydata, features=gene2, group.colors=colors, label=FALSE)+scale_fill_gradientn(colors=c("white", "snow", "firebrick3"))+theme(text=element_text(family="Times"))
ggsave("pyroptosis_heatmap.pdf", p, width=10, height=5)

GZMA<-FeaturePlot(mydata, features=c("GZMA"), cols=c("lightgray", "red"))+NoLegend()
ggsave("./GZMA.pdf", GZMA, width=6, height=6)
GSDMB<-FeaturePlot(mydata, features=c("GSDMB"), cols=c("lightgray", "red"))+NoLegend()
ggsave("./GSDMB.pdf", GSDMB, width=6, height=6)
GZMAv<-VlnPlot(mydata, features=c("GZMA"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())
ggsave("./GZMA_Vln.pdf", GZMAv, width=6, height=6)
GSDMBv<-VlnPlot(mydata, features=c("GSDMB"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())
ggsave("./GSDMB_Vln.pdf", GSDMBv, width=6, height=6)



cell_DEGs.list=split(x=markers,f=markers$cluster)
cell_DEGs.list=sapply(cell_DEGs.list, function(x){subset(x,select='gene',drop=TRUE)})
library(topGO)
library(clusterProfiler)
library(Rgraphviz)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(ggpubr)
subcluster.name=names(cell_DEGs.list)
subcluster.name=subcluster.name[order(subcluster.name)]
subcluster.name=c("0","1","4","5")
p=list()
for (s in subcluster.name) {
  erich.go.BP = enrichGO(gene =  cell_DEGs.list[[s]],OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.05)
  erich.go.BP.res=erich.go.BP@result
  write.csv(erich.go.BP.res,paste0('./',s,'_enrichment.csv'))
  erich.go.BP.res=erich.go.BP.res[erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
  erich.go.BP.res=erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust)
  p[[which(subcluster.name==s)]]=ggplot(data = erich.go.BP.res,
                                        mapping = aes(x=Count,y=reorder(Description,Count), fill = -log10(p.adjust))) +
    geom_bar(stat="identity")+ theme_bw()+
    scale_fill_gradient(low="grey80",high =c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02")[which(subcluster.name==s)])+
    labs(y=NULL,title = "subcluster ")+ggtitle(s)+
    theme(text = element_text(family = 'Times',size = 14))+
    scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
  
}
length(p)
subcluster.enrichment.fig <- plot_grid(plotlist = p, ncol = 3, nrow = 2)
pdf('subcluster_enrichment.pdf',height = 12,width = 18,onefile = F)
subcluster.enrichment.fig
dev.off()

library(msigdbr) 
library(fgsea)
library(Seurat)


mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
#mdb_kegg = mdb_c2 [grep("^KEGG",mdb_c2 $gs_name),]
fgsea_sets<- mdb_c2 %>% split(x = .$gene_symbol, f = .$gs_name)

cmarkers <- FindMarkers(mydata,ident.1 = "Melanocyte", only.pos = TRUE, 
                        min.pct = 0.1, logfc.threshold = 0.25)
head(cmarkers)

cmarkers$genes = rownames(cmarkers)
cluster0.genes<- cmarkers %>% arrange(desc(avg_logFC)) %>% dplyr::select(genes,avg_logFC)
ranks<- deframe(cluster0.genes)

#fgsea
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
term<-fgseaRes$pathway
scores <- fgseaRes$ES
pvalues <- fgseaRes$pval
result_df <- data.frame(term,scores, pvalues)
write.csv(result_df, "gsea_Melanocyte.csv", row.names = FALSE)
a=plotEnrichment(fgsea_sets[["REACTOME_PYROPTOSIS"]],
                 ranks) + labs(title="REACTOME_PYROPTOSIS")

ggsave("REACTOME_PYROPTOSIS_Melanocyte.pdf", a, width=8, height=5)


library(CellChat)
cellchat = createCellChat(object=mydata, group.by="cell_type", meta=mydata@meta.data)
groupSize = as.numeric(table(cellchat@idents))
CellChatDB = CellChatDB.human
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search="Cell-Cell Contact")
cellchat@DB <- CellChatDB.use

cellchat = subsetData(cellchat)
future::plan("multisession", workers=4)

cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)  
cellchat = projectData(cellchat, PPI.human)

#####################################################################################

cellchat = computeCommunProb(cellchat, raw.use=FALSE, population.size=TRUE)

cellchat = filterCommunication(cellchat, min.cells=10)

cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)
#####################################################################################

b=netVisual_bubble(cellchat, sources.use=c(2,3,4,5,6,7,8,9), targets.use=c(1), color.heatmap="Spectral")+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=15, hjust=0.5, face="bold", size=10), axis.text.y=element_text(face="bold", size=10))
ggsave("cell_chat.pdf", b, width=12, height=15)