options(stringsAsFactors = F)

library(AnnoProbe)
library(GEOquery)  
library(ggplot2) 
library(ggpubr) 
library(ggrepel)
library(reshape2)
library(patchwork)
library(stringr)
library(affy)
library(dplyr)
library(affyQCReport)
library(affyPLM)
library(clusterProfiler)
library(org.Hs.eg.db)
library(FactoMineR)
library(factoextra)  
library(pheatmap)
library(GSVA)
library(msigdbr)
library(xlsx)
library(ConsensusClusterPlus)
library(estimate)
library(ChAMP)
library(systemPipeR)
library(UpSetR)
library(msigdbr)
library(VennDiagram)
library(enrichplot)



call <- data.frame(msigdbr())
call <- call[call$gs_cat %in% c("H","C2","C5"),c(3,5,13)]
call <- call[grep("UP$",call$gs_name,invert = T),]
call <- call[grep("DN$",call$gs_name,invert = T),]

c_go_correlation <- data.frame(msigdbr())
c_go_correlation <- c_go_correlation[c_go_correlation$gs_cat %in% c("C5"),c(2,3)] %>% unique()
colnames(c_go_correlation) <- c("ONTOLOGY","ID")
c_go_correlation$ONTOLOGY <- gsub(pattern = "GO:",replacement = "",c_go_correlation$ONTOLOGY)



probe_anno <- read.table("~/zjw/20210927retinoblastoma/GPL28718_HGU133Plus2_Hs_ENTREZG_desc.annot.txt")
colnames(probe_anno) <- c("probe_ID","ENTREZID")
probe_anno <- merge(probe_anno,bitr(geneID = probe_anno[,2],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db),by="ENTREZID")



 

## 28
immune_cell_signature <- read.table("~/zjw/20210927retinoblastoma/28_immune_gene_list.tsv",sep = "\t",header = 1)








immune_cell_signature_list <- list()
for(i in 1:length(unique(immune_cell_signature$Cell.type))){
  immune_cell_signature_list[[i]] <- immune_cell_signature$Metagene[immune_cell_signature$Cell.type== (unique(immune_cell_signature$Cell.type)[i])]
}
names(immune_cell_signature_list)<- unique(immune_cell_signature$Cell.type)










risk_factor <- read.xlsx("~/zjw/20210927retinoblastoma/Supplementary_data1_from_Liu_et_al.xlsx",sheetIndex = 1)[1:102,1:14]
risk_factor <- risk_factor[,c(1,4,7,8,9,10,11,12)]
risk_factor$Germline.RB1.mutation[risk_factor$Germline.RB1.mutation=="Yes/Mosaicism"] <- "yes"
risk_factor[risk_factor=="NA"] <- NA
risk_factor$Diameter..mm. <- as.numeric(risk_factor$Diameter..mm.)








gset <- getGEO("GSE58780", GSEMatrix =TRUE, getGPL=FALSE)
expression_matrix <- exprs(gset[[1]]) 
pd <- pData(gset[[1]])[,c(2,1)]
pd$sample.ID <- pd$title %>% 
  gsub(pattern = "_REP2",replacement = "") %>% 
  gsub(pattern = "_REP1",replacement = "") %>%
  strsplit(split = " ") %>%
  lapply(function(x) x[1])%>%
  unlist()





expression_matrix[1:4,1:4] 
expression_matrix_symbol <- data.frame(expression_matrix)
expression_matrix_symbol$probe_ID <- rownames(expression_matrix_symbol)
expression_matrix_symbol <- merge(expression_matrix_symbol,probe_anno,by="probe_ID")
rownames(expression_matrix_symbol) <- expression_matrix_symbol$SYMBOL
expression_matrix_symbol <- expression_matrix_symbol[,!colnames(expression_matrix_symbol) %in% colnames(probe_anno)]
expression_matrix_symbol <- as.matrix(expression_matrix_symbol)
expression_matrix_symbol <- expression_matrix_symbol[,c(1:59)]
# expression_matrix_symbol <- expression_matrix_symbol[,1:63]
# expression_matrix_symbol <- expression_matrix_symbol[,1:63 %>% setdiff(c(5,7,8,28))]


## ssgsea
gsva_matrix <- gsva(expr = as.matrix(expression_matrix_symbol), gset.idx.list = immune_cell_signature_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)




## heatmap for cluster  
## ward.D", "ward.D2", "single", "complete", "average" "mcquitty"  "median"  "centroid" 
p <- pheatmap(gsva_matrix,
              clustering_method = "ward.D",
              scale = "row",
              cellheight = 10,
              ccellwidth = 10,
              fontsize = 10)
# p <- pheatmap(cor(gsva_matrix),clustering_method = "ward.D")
# p$tree_col
clusters <- cutree(p$tree_col, k = 2)
temp <- names(clusters)
clusters <- as.character(clusters)
table(clusters)
clusters <- paste0("C",clusters)
names(clusters) <- temp
clusters
clusters_df <- data.frame(group=clusters,
                          geo_accession=names(clusters))




annotation_col = data.frame(cluster=clusters)
rownames(annotation_col)<-colnames(gsva_matrix)

ann_colors = list(
  # Time = c("white", "firebrick"),
  cluster = c("High-ICI" = "#00468B", "Low-ICI" = "#ED0000")
  # GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)

p <- pheatmap(gsva_matrix,
              clustering_method = "ward.D",
              scale = "row",show_colnames = F,
              # border_color = NA,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              cellheight = 10,
              ccellwidth = 10,
              fontsize = 10)



original_classify <- read.xlsx("/home/zjw/zjw/20210927retinoblastoma/Supplementary_data1_from_Liu_et_al.xlsx",sheetIndex = 5)
colnames(original_classify) <- original_classify[3,]
original_classify <- original_classify[c(-1,-2,-3),c(1,2,3)] %>% na.omit()
colnames(original_classify)[c(1,3)] <- c("sample.ID","original.group")

original_classify_new <- merge(clusters_df[,c("sample.ID","group")],original_classify,by="sample.ID")




## PCA
res.pca <- PCA(t(gsva_matrix), graph = FALSE)
# fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
# fviz_pca_var(res.pca, col.var = "black")
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = clusters, # color by groups
             palette = c("#00468B","#ED0000"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")





clusters_df <- merge(clusters_df,pd,by="geo_accession")
clusters_df <- merge(clusters_df,risk_factor,by="sample.ID")
clusters_df$Diameter_group <- as.character(clusters_df$Diameter..mm.)
clusters_df$Diameter_group[clusters_df$Diameter..mm.>3.98 & clusters_df$Diameter..mm.<=6.67] <- "3.98-6.67"
clusters_df$Diameter_group[clusters_df$Diameter..mm.>6.67 & clusters_df$Diameter..mm.<=9.33] <- "6.67–9.33"
clusters_df$Diameter_group[clusters_df$Diameter..mm.>9.33 & clusters_df$Diameter..mm.<=12] <- "9.33–12"
clusters_df$Diameter_group[clusters_df$Diameter..mm.>12   & clusters_df$Diameter..mm.<=14.7] <- "12–14.7"
clusters_df$Diameter_group[clusters_df$Diameter..mm.>14.7 & clusters_df$Diameter..mm.<=17.3] <- "14.7–17.3"
clusters_df$Diameter_group[clusters_df$Diameter..mm.>17.3 & clusters_df$Diameter..mm.<=20] <- "17.3–20"

clusters_df$Diameter_group <- factor(clusters_df$Diameter_group,levels = c("3.98-6.67","6.67–9.33","9.33–12","12–14.7","14.7–17.3","17.3–20"))

clusters_df$Germline.RB1.mutation[clusters_df$Germline.RB1.mutation=="yes"] <- "RB1 mutated"
clusters_df$Germline.RB1.mutation[clusters_df$Germline.RB1.mutation=="no"] <- "RB1 not-mutated"
clusters_df$Germline.RB1.mutation <- factor(clusters_df$Germline.RB1.mutation,levels = c("RB1 mutated","RB1 not-mutated"))

clusters_df$Laterality <- factor(clusters_df$Laterality,levels = c("bilateral","unilateral"))

clusters_df$Growth <- factor(clusters_df$Growth,levels = c("endophytic","exophytic","mixed"))

clusters_df$Necrosis[clusters_df$Necrosis=="necrosis_present"] <- "present"
clusters_df$Necrosis <- factor(clusters_df$Necrosis,levels = c("present","none"))

clusters_df$Optic.nerve.invasion <- factor(clusters_df$Optic.nerve.invasion,levels = c("none","prelaminar","intralaminar","postlaminar"))

clusters_df$Choroid.and.sclera.invasion[clusters_df$Choroid.and.sclera.invasion=="choroid_minimal"] <- "minimal"
clusters_df$Choroid.and.sclera.invasion[clusters_df$Choroid.and.sclera.invasion=="choroid_deep"] <- "deep"
clusters_df$Choroid.and.sclera.invasion[clusters_df$Choroid.and.sclera.invasion=="choroid_extended"] <- "extended"
clusters_df$Choroid.and.sclera.invasion[clusters_df$Choroid.and.sclera.invasion=="choroid_extended_and_sclera"] <- "sclera incasion"

clusters_df$Choroid.and.sclera.invasion <- factor(clusters_df$Choroid.and.sclera.invasion,levels = c("none","minimal","deep","extended","sclera incasion"))




a <- na.omit(clusters_df[,c(3,12)])
a <- plyr::count(a)

chisq.test(matrix(c(0,3,6,5,2,1,4,1,21,10),byrow = F,ncol = 2))

a <- na.omit(clusters_df[,c(3,12)])
a <- plyr::count(a)
p1 <- ggplot(a,aes(x = group, y = freq,fill = Diameter_group)) +
  geom_bar(position = "fill",stat = "identity",width = 0.8,color="black") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x="",y="Proportion",
       title = paste("germline",as.character(chisq.test(matrix(a$freq,byrow = F,nrow = 2))$p.value),sep = "   ")) +
  scale_fill_manual(values = as.character(jdb_palette("brewer_blue",type="continuous"))[c(1:5)*200]) 



a <- na.omit(clusters_df[,c(3,5)])
a <- plyr::count(a)
p2 <- ggplot(a,aes(x = group, y = freq,fill = Germline.RB1.mutation)) +
  geom_bar(position = "fill",stat = "identity",width = 0.8,color="black") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x="",y="Proportion",
       title = paste("germline",as.character(chisq.test(matrix(a$freq,byrow = F,nrow = 2))$p.value),sep = "   ")) +
  scale_fill_manual(values = c("black","grey"))




a <- na.omit(clusters_df[,c(3,6)])
a <- plyr::count(a)
p3 <- ggplot(a,aes(x = group, y = freq,fill = Laterality)) +
  geom_bar(position = "fill",stat = "identity",width = 0.8,color="black") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x="",y="Proportion",
       title = paste("Laterality",as.character(fisher.test(matrix(a$freq,byrow = T,nrow = 2))$p.value),sep = "   "))+
  scale_fill_manual(values = c("#000080","#BBFFFF"))





a <- na.omit(clusters_df[,c(3,7)])
a <- plyr::count(a)
p4 <- ggplot(a,aes(x = group, y = freq,fill = Growth)) +
  geom_bar(position = "fill",stat = "identity",width = 0.8,color="black") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x="",y="Proportion",
       title = paste("Growth",as.character(fisher.test(matrix(a$freq,byrow = T,nrow = 2))$p.value),sep = "   ")) +
  scale_fill_manual(values = c("purple","yellow","#00CD66"))




a <- na.omit(clusters_df[,c(3,9)])
a <- plyr::count(a)
p5 <- ggplot(a,aes(x = group, y = freq,fill = Necrosis)) +
  geom_bar(position = "fill",stat = "identity",width = 0.8,color="black") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x="",y="Proportion",
       title = paste("Necrosis",as.character(fisher.test(matrix(a$freq,byrow = T,nrow = 2))$p.value),sep = "   "))+
  scale_fill_manual(values = c("#8B814C","grey"))






a <- na.omit(clusters_df[,c(3,10)])
a <- plyr::count(a)
p6 <- ggplot(a,aes(x = group, y = freq,fill = Optic.nerve.invasion)) +
  geom_bar(position = "fill",stat = "identity",width = 0.8,color="black") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x="",y="Proportion",
       title = paste("Optic.nerve.invasion",as.character(fisher.test(matrix(a$freq,byrow = T,nrow = 2))$p.value),sep = "   "))+
  scale_fill_manual(values = as.character(jdb_palette("ocean_green",type="continuous"))[c(5:8)*125]) 




a <- na.omit(clusters_df[,c(3,11)])
a <- plyr::count(a)
p7 <- ggplot(a,aes(x = group, y = freq,fill = Choroid.and.sclera.invasion)) +
  geom_bar(position = "fill",stat = "identity",width = 0.8,color="black") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x="",y="Proportion",
       title = paste("Choroid.and.sclera.invasion",as.character(fisher.test(matrix(a$freq,byrow = T,nrow = 2))$p.value),sep = "   ")) +
  scale_fill_manual(values = as.character(jdb_palette("solar_glare",type="continuous"))[c(550,650,750,850,950)]) 


ggarrange(plotlist = list(p1,p2,p3,p4,p5,p6,p7))




box <- melt(gsva_matrix)
colnames(box) <- c("type","geo_accession","NES")
box <- merge(box,data.frame(group=clusters,
                            geo_accession=names(clusters)),by="geo_accession")
box$type <- factor(box$type,levels = unique(immune_cell_signature$Cell.type))
ggboxplot(box, x = "type", y = "NES", fill = "group" ) +
  stat_compare_means(aes(group=group), label = "p.signif",method = "t.test")+
  labs(y="NES",x="") +
  scale_fill_manual(values = c("#00468B","#ED0000")) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=45,hjust = 1,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  guides(fill=FALSE)



p <- pheatmap(cor(t(gsva_matrix)),
              clustering_method = "ward.D",
              cluster_rows = T,
              cluster_cols = T,
              cellwidth=10,
              cellheight=10,
              fontsize=10,
              bk = unique(c(seq(0,1, length=100))),
              scale = "none")







## ESTIMATE
write.table(as.data.frame(expression_matrix_symbol),"/home/zjw/zjw/20210927retinoblastoma/expression_matrix_symbol.txt",quote = F,row.names = T,col.names = T,sep = "\t")
filterCommonGenes(input.f="/home/zjw/zjw/20210927retinoblastoma/expression_matrix_symbol.txt",#
                  output.f="/home/zjw/zjw/20210927retinoblastoma/RB_9699genes.gct",#
                  id="GeneSymbol")#
estimateScore(input.ds = "/home/zjw/zjw/20210927retinoblastoma/RB_9699genes.gct",
              output.ds="/home/zjw/zjw/20210927retinoblastoma/RB_estimate_score.gct", 
              platform="affymetrix")
plotPurity(scores="/home/zjw/zjw/20210927retinoblastoma/RB_estimate_score.gct", samples="GSM5121283", 
           platform="affymetrix")
scores <- read.table("/home/zjw/zjw/20210927retinoblastoma/RB_estimate_score.gct",skip = 2,header = T)[,-1]
rownames(scores) <- scores[,1]
scores <- t(scores)[-1,]
scores <- melt(scores)
colnames(scores)[1] <- "geo_accession"
scores$value <- as.numeric(scores$value)
scores <- merge(scores,clusters_df[,2:3],by="geo_accession")

my_comparisons <- list(c("High-ICI", "Low-ICI"))
p1 <- ggboxplot(data = scores[scores$Var2=="StromalScore",],x = "group",y = "value",ylab = "Stromal Score",fill = "group",xlab = "") +
  scale_fill_manual(values = c("#00468B","#ED0000")) +
  stat_compare_means(comparisons=my_comparisons, label = "p.signif",method = "t.test") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  guides(fill=FALSE)

p2 <- ggboxplot(data = scores[scores$Var2=="ImmuneScore",],x = "group",y = "value",ylab = "Immune Score",fill = "group",xlab = "") +
  scale_fill_manual(values = c("#00468B","#ED0000")) +
  stat_compare_means(comparisons=my_comparisons, label = "p.signif",method = "t.test") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  guides(fill=FALSE)

p3 <- ggboxplot(data = scores[scores$Var2=="ESTIMATEScore",],x = "group",y = "value",ylab = "ESTIMATE Score",fill = "group",xlab = "") +
  scale_fill_manual(values = c("#00468B","#ED0000")) +
  stat_compare_means(comparisons=my_comparisons, label = "p.signif",method = "t.test") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  guides(fill=FALSE)

p4 <- ggboxplot(data = scores[scores$Var2=="TumorPurity",],x = "group",y = "value",ylab = "Tumor Purity",fill = "group",xlab = "") +
  scale_fill_manual(values = c("#00468B","#ED0000")) +
  stat_compare_means(comparisons=my_comparisons, label = "p.signif",method = "t.test") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=0,hjust = 0.5,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  guides(fill=FALSE)

ggarrange(plotlist = list(p1,p2,p3,p4),nrow = 1)







## Limma

design <- model.matrix(~0+factor(clusters %>% gsub(pattern = "-",replacement = "_"),levels = c("Low_ICI","High_ICI")))
colnames(design) <- levels(factor(clusters %>% gsub(pattern = "-",replacement = "_")))
rownames(design) <- colnames(expression_matrix_symbol)
design

contrast.matrix <- makeContrasts(paste0(unique(clusters %>% gsub(pattern = "-",replacement = "_")),collapse = "-"),levels = design)
contrast.matrix ##



##step1
fit <- lmFit(expression_matrix_symbol,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)




volcano <- nrDEG

volcano$abs <- abs(volcano$logFC)

volcano$gene <- rownames(volcano)

volcano$adj.P.Val <- -log2(volcano$adj.P.Val)


volcano$threshold <- as.factor(ifelse(volcano$adj.P.Val > (-log2(0.05)) & abs(volcano$logFC) >= log2(1.5), ifelse(volcano$logFC> log2(1.5) ,'Low-ICI up','High-ICI up'),'No significance'))


# volcano <- volcano[!volcano$gene %in% as.character(as.data.frame(table(volcano$gene))[!as.data.frame(table(volcano$gene))[,2]==1,1]),]

volcano$threshold <- factor(volcano$threshold,levels = c('High-ICI up','Low-ICI up','No significance'))


ggplot(data = volcano, aes(x = logFC, y = adj.P.Val, color=threshold)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","#ED0000","grey")) +
  xlim(c(-2.6, 2.6)) +
  #ylim(c(0, 13))+
  geom_vline(xintercept=c(-log2(1.5),log2(1.5)),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=4,col="black",lwd=0.5) +
  labs(x="log2 (fold change)",y="-log2 (adjusted p-value)") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11)) +
  geom_text_repel(
    data = volcano[volcano$gene %in% c("TOP2A","MKI67","CENPE","CENPF","TTL",
                                       "PDC","OPN1SW","HLA-DRA","C1QB","IGKC"),] ,
    aes(label = gene),
    size = 3.5,
    box.padding = unit(3, "lines"),
    point.padding = unit(0, "lines"), segment.color = "black", show.legend = FALSE ,
    force = T
  )


expression_matrix_symbol_DEG <- expression_matrix_symbol[rownames(expression_matrix_symbol) %in% volcano$gene[volcano$threshold %in% c("High-ICI up","Low-ICI up")],]
expression_matrix_symbol_DEG <- expression_matrix_symbol_DEG[,c(names(clusters)[clusters=="High-ICI"],names(clusters)[clusters=="Low-ICI"])]
# expression_matrix_symbol_DEG <- log2(expression_matrix_symbol_DEG+1)

# pheatmap(expression_matrix_symbol_DEG)

p <- pheatmap(expression_matrix_symbol_DEG,
              clustering_method = "ward.D2",
              scale = "row",
              cluster_cols = F,
              treeheight_row = 40,
              show_rownames = F,
              show_colnames = F,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              color = colorRampPalette(as.character(jdb_palette("brewer_celsius")))(200),
              # cellheight = 10,
              # ccellwidth = 10,
              fontsize = 10)


clusters_row <- cutree(p$tree_row, k = 3)
temp <- names(clusters_row)
clusters_row <- as.character(clusters_row)
table(clusters_row)
clusters_row <- paste0("Gene cluster ",clusters_row)
names(clusters_row) <- temp
annotation_row <- data.frame(Gene_cluster=clusters_row)

p <- pheatmap(expression_matrix_symbol_DEG,
              clustering_method = "ward.D2",
              scale = "row",
              cluster_cols = F,
              treeheight_row = 40,
              show_rownames = F,
              show_colnames = F,
              annotation_col = annotation_col,
              annotation_row = annotation_row,
              annotation_colors = ann_colors,
              color = colorRampPalette(as.character(jdb_palette("brewer_celsius")))(200),
              # cellheight = 10,
              # ccellwidth = 10,
              fontsize = 10)

go_genechuster1 <- enrichGO(rownames(annotation_row)[annotation_row$Gene_cluster=="Gene cluster 1"],OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "ALL")

a <- data.frame(go_up)

dotplot(go_genechuster1,
        title='Top5 GO terms of Gene cluster 1',
        showCategory=5,split='ONTOLOGY')+ 
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))


go_genechuster2 <- enrichGO(rownames(annotation_row)[annotation_row$Gene_cluster=="Gene cluster 2"],OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "ALL")

a <- data.frame(go_up)

dotplot(go_genechuster2,
        title='Top5 GO terms of Gene cluster 2',
        showCategory=5,split='ONTOLOGY')+ 
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))


go_genechuster3 <- enrichGO(rownames(annotation_row)[annotation_row$Gene_cluster=="Gene cluster 3"],OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "ALL")

a <- data.frame(go_up)

dotplot(go_genechuster3,
        title='Top5 GO terms of Gene cluster 3',
        showCategory=5,split='ONTOLOGY')+ 
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))





glist <- volcano$logFC
names(glist) <- as.character(volcano$gene)
glist <- sort(glist,decreasing = T)
glist <- glist[!duplicated(names(glist))]


gsea <- GSEA(glist, TERM2GENE = call[,c(1,2)],
             #TERM2NAME = call_new[,c(1,3)],
             pvalueCutoff = 0.05,pAdjustMethod = "none",maxGSSize = 10000,minGSSize = 0,exponent = 1,nPerm = 5000)

gseaplot2(gsea,geneSetID = c("WP_SENESCENCE_AND_AUTOPHAGY_IN_CANCER"),subplots = c(1,2),pvalue_table = F,title = "WP: senescence and autophagy in cancer") +
  gseaplot2(gsea,geneSetID = c("GO_RETINAL_CONE_CELL_DIFFERENTIATION"),subplots = c(1,2),pvalue_table = F,title = "GO: retinal cone cell differentiation") +
  gseaplot2(gsea,geneSetID = c("WP_RETINOBLASTOMA_GENE_IN_CANCER"),subplots = c(1,2),pvalue_table = F,title = "WP: retinoblastoma gene in cancer") 




aa <- data.frame(gsea)














load("/home/zjw/zjw/20210927retinoblastoma/GSE58783_Methylation_processed_data.RData")

pvalue_data <- processed_data[,grep("PValue",colnames(processed_data))]
processed_data <- processed_data[,grep("Beta",colnames(processed_data))]
colnames(pvalue_data) <- colnames(processed_data)




meth_sample <- read.csv("/home/zjw/zjw/20210927retinoblastoma/meth_sample.csv",header = F,col.names = c("Sample_Name","sample.ID"))[1:70,]
meth_sample$sample.ID <- meth_sample$sample.ID %>%
  strsplit(split = " ") %>%
  lapply(function(x) x[1]) %>%
  unlist() %>%
  gsub(pattern = "_REP1",replacement = "")
meth_sample <- meth_sample[meth_sample$sample.ID %in% clusters_df$sample.ID,]
meth_sample <- meth_sample[order(meth_sample$sample.ID),]
meth_sample <- merge(meth_sample,clusters_df[,c("sample.ID","group")],by="sample.ID")
meth_sample <- meth_sample[order(meth_sample$Sample_Name),]



processed_data <- processed_data[,colnames(processed_data) %in% meth_sample$Sample_Name]
pvalue_data <- pvalue_data[,colnames(pvalue_data) %in% meth_sample$Sample_Name]
processed_data <- processed_data[,order(colnames(processed_data))]
pvalue_data <- pvalue_data[,order(colnames(pvalue_data))]


myFilter <- champ.filter(beta = processed_data,
                         #---原始数据中的β矩阵，甲基化程度:<未-部分-完全> → <0.2 - 0.6>
                         pd = meth_sample,
                         #---pd(phenotype data)，来自Sample_sheet.csv文件，记录了每个样本的信息，如表型、编号等
                         intensity=NULL,
                         #---强度值矩阵
                         Meth=NULL,
                         #---甲基化矩阵
                         UnMeth=NULL,
                         #---非甲基化矩阵
                         detP=pvalue_data,
                         #---Detected P value identifies failed positions defined as both the methylated and unmethylated channel reporting background signal levels.每一个芯片上的数据的可信度
                         #---detP计算：https://blog.csdn.net/Joshua_HIT/article/details/72868459
                         beadcount=NULL,
                         #---Beadcount information for Green and Red Channal, need for filterBeads.(default = NULL) form ChAMP
                         #---number of beads used to summarize the red and green channel from minfi
                         autoimpute=TRUE,
                         #---参见champ.process()
                         filterDetP=TRUE,
                         #---是否过滤detP较大的探针
                         detPcut=0.01,
                         #---设定的detP阈值，值越小越good
                         SampleCutoff=0.1,
                         #---每个样本的failed detP value 占比超过SampleCutoff=0.1时该sample会被丢弃
                         ProbeCutoff=0,
                         #---在移除failed sample后，每个探针中failed detP value 占比超过ProbeCutoff时该探针会被丢弃
                         filterBeads=TRUE,
                         #---每个探针结合了序列的bead少于3个时会被设为NA
                         beadCutoff=0.05,
                         #---设定的beadcount frequency阈值
                         filterNoCG = TRUE,
                         #---是否过滤非CpG位点的探针
                         filterSNPs = TRUE,
                         #---是否过滤SNP位点附近的探针
                         filterMultiHit = TRUE,
                         #---是否过滤映射到基因组多个位置的探针
                         filterXY = TRUE,
                         #---是否过滤检测性染色体上的位点的探针
                         population = NULL,
                         #---筛选特定的人群，如"AFR"(Africa),人群列表见http://www.internationalgenome.org/category/population/ 
                         fixOutlier = TRUE,
                         #---是否修改异常值，当该参数为T时，若beta值小于0，则会用最小正数代替该值，若beta值大于1，则会用最大整数代替该值
                         arraytype = "450K")










myDMP <- champ.DMP(beta = myFilter$beta,
                   pheno = myFilter$pd$group %>% gsub(pattern = "-",replacement = "_"),
                   arraytype = "450K",
                   adjPVal = 1)



myDMP <- myDMP[[1]]



myDMP$DMP <- "No significance"
myDMP$DMP[myDMP$adj.P.Val <= 0.05 & myDMP$deltaBeta <= -0.1] <- "High-ICI hypermethylated probes"
myDMP$DMP[myDMP$adj.P.Val <= 0.05 & myDMP$deltaBeta >=  0.1] <- "Low-ICI hypermethylated probes"
myDMP$logP <- (-log2(myDMP$adj.P.Val))
myDMP$gene <- as.character(myDMP$gene)

myDMP$DMP <- as.character(myDMP$DMP)
myDMP$cgi <- as.character(myDMP$cgi)
myDMP$UCSC_CpG_Islands_Name <- as.character(myDMP$UCSC_CpG_Islands_Name)

myDMP$cgi_new <- myDMP$cgi
myDMP$temp <- as.character(myDMP$DMP)
myDMP$temp <- myDMP$UCSC_CpG_Islands_Name %>% strsplit("-") %>% lapply(function(x) x[2])%>% unlist()

myDMP$cgi_new[myDMP$Strand=="F" & myDMP$cgi=="shore" & myDMP$MAPINFO<myDMP$temp] <- "N shore"
myDMP$cgi_new[myDMP$Strand=="F" & myDMP$cgi=="shore" & myDMP$MAPINFO>myDMP$temp] <- "S shore"
myDMP$cgi_new[myDMP$Strand=="R" & myDMP$cgi=="shore" & myDMP$MAPINFO<myDMP$temp] <- "S shore"
myDMP$cgi_new[myDMP$Strand=="R" & myDMP$cgi=="shore" & myDMP$MAPINFO>myDMP$temp] <- "N shore"

myDMP$cgi_new[myDMP$Strand=="F" & myDMP$cgi=="shelf" & myDMP$MAPINFO<myDMP$temp] <- "N shelf"
myDMP$cgi_new[myDMP$Strand=="F" & myDMP$cgi=="shelf" & myDMP$MAPINFO>myDMP$temp] <- "S shelf"
myDMP$cgi_new[myDMP$Strand=="R" & myDMP$cgi=="shelf" & myDMP$MAPINFO<myDMP$temp] <- "S shelf"
myDMP$cgi_new[myDMP$Strand=="R" & myDMP$cgi=="shelf" & myDMP$MAPINFO>myDMP$temp] <- "N shelf"



MP_processed_data_new <- processed_data[rownames(processed_data) %in% rownames(myDMP[myDMP$DMP %in% c("High-ICI hypermethylated probes","C2 hypermethylated probes"),]),
                                        c(meth_sample$Sample_Name[meth_sample$group=="High-ICI"],meth_sample$Sample_Name[meth_sample$group=="C2"])]


p <- pheatmap(MP_processed_data_new,
              clustering_method = "ward.D",
              scale = "row",
              cluster_cols = F,
              treeheight_row = 0,
              show_rownames = F,
              show_colnames = F,
              annotation_col = data.frame(cluster=meth_sample$group,row.names = meth_sample$Sample_Name),
              annotation_colors = ann_colors,
              # annotation_row = rbind(data.frame(type=rep("High-ICI hypermethylated probes",nrow(myDMP[myDMP$DMP==c("High-ICI hypermethylated probes"),])),row.names = rownames(myDMP[myDMP$DMP==c("High-ICI hypermethylated probes"),])),
              #                        data.frame(type=rep("C2 hypermethylated probes",nrow(myDMP[myDMP$DMP==c("C2 hypermethylated probes"),])),row.names = rownames(myDMP[myDMP$DMP==c("C2 hypermethylated probes"),]))),
              color = colorRampPalette(as.character(jdb_palette("brewer_celsius")))(50),
              # cellheight = 10,
              # ccellwidth = 10,
              fontsize = 10)


C1h_cgi_df <- as.data.frame(myDMP[myDMP$DMP=="High-ICI hypermethylated probes",]$cgi_new %>% table())
colnames(C1h_cgi_df) <- c("cgi","Freq")
C1h_cgi_df$cgi <- factor(C1h_cgi_df$cgi,levels = c("N shelf","N shore","island","S shore","S shelf","opensea"))
C1h_cgi_df[,1] <- paste(C1h_cgi_df[,1]," (",round(C1h_cgi_df[,2]/sum(C1h_cgi_df[,2])*100,2),"%)",sep = "")
C1h_cgi_df$cgi <- factor(C1h_cgi_df$cgi, levels=c(grep("N shelf",C1h_cgi_df[,1],value = T),
                                                  grep("N shore",C1h_cgi_df[,1],value = T),
                                                  grep("island",C1h_cgi_df[,1],value = T),
                                                  grep("S shore",C1h_cgi_df[,1],value = T),
                                                  grep("S shelf",C1h_cgi_df[,1],value = T),
                                                  grep("opensea",C1h_cgi_df[,1],value = T)), ordered=TRUE)
p1 <- ggplot(C1h_cgi_df, aes(x = 1, weight = Freq, fill = cgi)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL) +
  scale_x_continuous(breaks = NULL) 
  # labs(title = "C1 hypermethylated probes")

p2 <- ggbarplot(data = C1h_cgi_df,x = "cgi", y = "Freq",fill = "cgi") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1)) +
  scale_x_discrete(labels=c("N shelf","N shore","island","S shore","S shelf","opensea")) +
  labs(x="",y="Number of probe") +
  geom_text(mapping = aes(label=Freq),vjust = 0, nudge_y = 0.7)+
  guides(fill=FALSE,color=FALSE)

C2h_cgi_df <- as.data.frame(myDMP[myDMP$DMP=="Low-ICI hypermethylated probes",]$cgi_new %>% table())
colnames(C2h_cgi_df) <- c("cgi","Freq")
C2h_cgi_df$cgi <- factor(C2h_cgi_df$cgi,levels = c("N shelf","N shore","island","S shore","S shelf","opensea"))
C2h_cgi_df[,1] <- paste(C2h_cgi_df[,1]," (",round(C2h_cgi_df[,2]/sum(C2h_cgi_df[,2])*100,2),"%)",sep = "")
C2h_cgi_df$cgi <- factor(C2h_cgi_df$cgi, levels=c(grep("N shelf",C2h_cgi_df[,1],value = T),
                                                  grep("N shore",C2h_cgi_df[,1],value = T),
                                                  grep("island",C2h_cgi_df[,1],value = T),
                                                  grep("S shore",C2h_cgi_df[,1],value = T),
                                                  grep("S shelf",C2h_cgi_df[,1],value = T),
                                                  grep("opensea",C2h_cgi_df[,1],value = T)), ordered=TRUE)
p3 <- ggplot(C2h_cgi_df, aes(x = 1, weight = Freq, fill = cgi)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL) +
  scale_x_continuous(breaks = NULL) 
  # labs(title = "C2 hypermethylated probes")

p4 <- ggbarplot(data = C2h_cgi_df,x = "cgi", y = "Freq",fill = "cgi") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1)) +
  scale_x_discrete(labels=c("N shelf","N shore","island","S shore","S shelf","opensea")) +
  labs(x="",y="Number of probe") +
  geom_text(mapping = aes(label=Freq),vjust = 0, nudge_y = 0.7)+
  guides(fill=FALSE,color=FALSE)


p2+p1+p4+p3



rbind(myDMP[myDMP$DMP %in% c("High-ICI hypermethylated probes","Low-ICI hypermethylated probes"),],
      myDMP[myDMP$DMP %in% c("No significance"),][sample(1:400000,size = 400000),])



ggplot(data = myDMP, aes(x = deltaBeta, y = logP, color=DMP)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","#ED0000","grey")) +
  xlim(c(-0.45, 0.45)) +
  #ylim(c(0, 13))+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=4,col="black",lwd=0.5) +
  labs(x="Δβ",y="-log2 (adjusted p-value)") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11))



myDMP_C1h <- myDMP[myDMP$DMP=="High-ICI hypermethylated probes",]
myDMP_C1h <- myDMP_C1h[!myDMP_C1h$gene=="",]

upset_C1h <- data.frame(unique(myDMP_C1h$gene),0,0,0,0,0,0)
colnames(upset_C1h) <- c("gene","TSS1500","TSS200","1stExon","5'UTR","Body","3'UTR")


for (i in 1:length(upset_C1h$gene)) {
  temp <- myDMP_C1h[myDMP_C1h$gene==upset_C1h$gene[i],]$feature %>% unique()
  for (j in temp) {
    eval(parse(text=paste0('upset_C1h$`',j,'`[',i,'] <- 1')))
  }
}

upset(upset_C1h, 
      nsets = 6, 
      nintersects = 100, 
      sets = c("Body","TSS1500","5'UTR","3'UTR","TSS200","1stExon"),
      mb.ratio = c(0.55, 0.45),
      order.by = c("freq", "degree"), 
      decreasing = c(TRUE,FALSE),
      main.bar.color = rev(colorRampPalette(as.character(jdb_palette("brewer_orange")))(29))[1:29],
      sets.bar.color = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3", 'brown2'),
      text.scale = c(1.5, 1.2, 1.2, 1, 1.4, 1.4), # ytitle, ylabel, xtitle, xlabel, sets, number,
      sets.x.label = "Number of probe"
      )


C1h_feature_df <- as.data.frame(myDMP[myDMP$DMP=="High-ICI hypermethylated probes",]$feature %>% table())
C1h_feature_df <- C1h_feature_df[!C1h_feature_df$.=="IGR",]
colnames(C1h_feature_df) <- c("feature","Freq")
C1h_feature_df$feature <- factor(C1h_feature_df$feature,levels = c("Body","TSS1500","5'UTR","3'UTR","TSS200","1stExon"))
C1h_feature_df[,1] <- paste(C1h_feature_df[,1]," (",round(C1h_feature_df[,2]/sum(C1h_feature_df[,2])*100,2),"%)",sep = "")
C1h_feature_df$feature <- factor(C1h_feature_df$feature, levels=c(grep("Body",C1h_feature_df[,1],value = T),
                                                                  grep("TSS1500",C1h_feature_df[,1],value = T),
                                                                  grep("5'UTR",C1h_feature_df[,1],value = T),
                                                                  grep("3'UTR",C1h_feature_df[,1],value = T),
                                                                  grep("TSS200",C1h_feature_df[,1],value = T),
                                                                  grep("1stExon",C1h_feature_df[,1],value = T)), ordered=TRUE)

ggplot(C1h_feature_df, aes(x = 1, weight = Freq, fill = feature)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL) +
  scale_x_continuous(breaks = NULL)  +
  scale_fill_manual(values = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3", 'brown2'))



myDMP_C2h <- myDMP[myDMP$DMP=="Low-ICI hypermethylated probes",]
myDMP_C2h <- myDMP_C2h[!myDMP_C2h$gene=="",]

upset_C2h <- data.frame(unique(myDMP_C2h$gene),0,0,0,0,0,0)
colnames(upset_C2h) <- c("gene","TSS1500","TSS200","1stExon","5'UTR","Body","3'UTR")


for (i in 1:length(upset_C2h$gene)) {
  temp <- myDMP_C2h[myDMP_C2h$gene==upset_C2h$gene[i],]$feature %>% unique()
  for (j in temp) {
    eval(parse(text=paste0('upset_C2h$`',j,'`[',i,'] <- 1')))
  }
}

upset(upset_C2h, 
      nsets = 6, 
      nintersects = 100, 
      sets = c("Body","5'UTR","TSS1500","3'UTR","TSS200","1stExon"),
      mb.ratio = c(0.55, 0.45),
      order.by = c("freq", "degree"), 
      decreasing = c(TRUE,FALSE),
      main.bar.color = rev(colorRampPalette(as.character(jdb_palette("brewer_orange")))(12))[1:15],
      sets.bar.color = c("dodgerblue", "darkorange1", "goldenrod1", "orchid3", "seagreen3", "brown2"),
      text.scale = c(1.5, 1.2, 1.2, 1, 1.4, 1.4), # ytitle, ylabel, xtitle, xlabel, sets, number,
      sets.x.label = "Number of probe"
)


C2h_feature_df <- as.data.frame(myDMP[myDMP$DMP=="Low-ICI hypermethylated probes",]$feature %>% table())
C2h_feature_df <- C2h_feature_df[!C2h_feature_df$.=="IGR",]
colnames(C2h_feature_df) <- c("feature","Freq")
C2h_feature_df$feature <- factor(C2h_feature_df$feature,levels = c("Body","5'UTR","TSS1500","TSS200","3'UTR","1stExon"))
C2h_feature_df[,1] <- paste(C2h_feature_df[,1]," (",round(C2h_feature_df[,2]/sum(C2h_feature_df[,2])*100,2),"%)",sep = "")
C2h_feature_df$feature <- factor(C2h_feature_df$feature, levels=c(grep("Body",C2h_feature_df[,1],value = T),
                                                                  grep("5'UTR",C2h_feature_df[,1],value = T),
                                                                  grep("TSS1500",C2h_feature_df[,1],value = T),
                                                                  grep("TSS200",C2h_feature_df[,1],value = T),
                                                                  grep("3'UTR",C2h_feature_df[,1],value = T),
                                                                  grep("1stExon",C2h_feature_df[,1],value = T)), ordered=TRUE)

ggplot(C2h_feature_df, aes(x = 1, weight = Freq, fill = feature)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  labs(x=NULL,y=NULL) +
  scale_x_continuous(breaks = NULL)  +
  scale_fill_manual(values = c("dodgerblue", "darkorange1", "goldenrod1", "orchid3", "seagreen3", "brown2"))




temp <- as.data.frame(table(as.character(myDMP$gene)))
temp$Var1 <- as.character(temp$Var1)
temp <- temp[order(temp$Freq,decreasing = T),]
temp <- temp[!temp$Var1=="",]

volcano_meth <- myDMP[myDMP$gene %in% temp$Var1[temp$Freq==1],c("adj.P.Val","deltaBeta","gene","DMP")]

for (i in temp$Var1[temp$Freq>1]) {
  print(i)
  temp_temp <- myDMP[myDMP$gene==i,]
  if (temp_temp$DMP %>% table() %>% length()==1) {
    volcano_meth <- rbind(volcano_meth,
                          data.frame(adj.P.Val=mean(temp_temp$adj.P.Val),
                                     deltaBeta=mean(temp_temp$deltaBeta),
                                     gene=i,
                                     DMP=temp_temp$DMP[1]
                                  ))
  } else {
    if (
      ("High-ICI hypermethylated probes" %in% (temp_temp$DMP %>% table() %>% names())) 
      & 
      ("Low-ICI hypermethylated probes" %in% (temp_temp$DMP %>% table() %>% names()))
      ) {
      
    } else {
      temp_temp <- temp_temp[!temp_temp$DMP=="No significance",]
      volcano_meth <- rbind(volcano_meth,
                            data.frame(adj.P.Val=mean(temp_temp$adj.P.Val),
                                       deltaBeta=mean(temp_temp$deltaBeta),
                                       gene=i,
                                       DMP=temp_temp$DMP[1]
                         ))
    }
  }
}



volcano_meth$DMP <- volcano_meth$DMP %>% gsub(pattern = "probes",replacement = "gene")

volcano_meth$logP <- (-log2(volcano_meth$adj.P.Val))

ggplot(data = volcano_meth, aes(x = deltaBeta, y = logP, color=DMP)) +
  geom_point(alpha=0.9, size=1.3) +
  scale_color_manual(values=c("#00468B","#ED0000","grey")) +
  xlim(c(-0.4, 0.4)) +
  #ylim(c(0, 13))+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(0.05),lty=4,col="black",lwd=0.5) +
  labs(x="Δβ",y="-log2 (adjusted p-value)") +
  theme_classic() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  # theme(panel.border = element_rect(linetype=1,color="black",size=2))+
  theme(axis.line.x =element_line(linetype=1,color="black",size=1),
        axis.line.y = element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        #legend.position = c(0.8,0.7),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11)) +
  geom_text_repel(
    data = volcano_meth[volcano_meth$gene %in% c("CD83","HLA-DOA","IRF4","DOK3","CXCR1"),] ,
    aes(label = gene),
    size = 3.5,
    box.padding = unit(3, "lines"),
    point.padding = unit(0, "lines"), segment.color = "black", show.legend = FALSE ,
    force = T
  )


vennPlot(overLapper(list("High-ICI upregulated genes"=volcano$gene[volcano$threshold=="High-ICI up"],
                         "Low-ICI upregulated genes"=volcano$gene[volcano$threshold=="Low-ICI up"],
                         "High-ICI hypermethylated genes"=volcano_meth$gene[volcano_meth$DMP=="High-ICI hypermethylated gene"] %>% as.character(),
                         "Low-ICI hypermethylated genes"=volcano_meth$gene[volcano_meth$DMP=="Low-ICI hypermethylated gene"] %>% as.character()), type="vennsets"
                    )
         )


write.table(data.frame(c(intersect(volcano$gene[volcano$threshold=="Low-ICI up"],volcano_meth$gene[volcano_meth$DMP=="High-ICI hypermethylated gene"] %>% as.character()),
                         intersect(volcano$gene[volcano$threshold=="High-ICI up"],volcano_meth$gene[volcano_meth$DMP=="Low-ICI hypermethylated gene"] %>% as.character()))),"/home/zjw/zjw/20210927retinoblastoma/BOTH_for_PPI.txt",quote = F,row.names = F,col.names = F,sep = "\t")





enrich_C1u <- enricher(volcano[volcano$threshold=="High-ICI up",] %>% rownames(), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
enrich_C2u <- enricher(volcano[volcano$threshold=="Low-ICI up",] %>% rownames(), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
enrich_C1h <- enricher(volcano_meth$gene[volcano_meth$DMP=="High-ICI hypermethylated gene"] %>% as.character(), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
enrich_C2h <- enricher(volcano_meth$gene[volcano_meth$DMP=="Low-ICI hypermethylated gene"] %>% as.character(), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
enrich_C1uC2h <- enricher(intersect(volcano$gene[volcano$threshold=="High-ICI up"],volcano_meth$gene[volcano_meth$DMP=="Low-ICI hypermethylated gene"] %>% as.character()), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
enrich_C1hC2u <- enricher(intersect(volcano$gene[volcano$threshold=="Low-ICI up"],volcano_meth$gene[volcano_meth$DMP=="High-ICI hypermethylated gene"] %>% as.character()), TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
enrich_G1 <- enricher(rownames(annotation_row)[annotation_row$Gene_cluster=="Gene cluster 1"], TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
enrich_G2 <- enricher(rownames(annotation_row)[annotation_row$Gene_cluster=="Gene cluster 2"], TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)
enrich_G3 <- enricher(rownames(annotation_row)[annotation_row$Gene_cluster=="Gene cluster 3"], TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)

enrich_temp <- enricher(c("ALCAM","EGF","FARP1","PROM1","ALK","KIRREL3","MDK","BIRC5",
                          "OIP5","SKP2","SMC4","CDCA2","NCAPD2","KIFC1","CDC20","LMNB1",
                          "LRP8","CEP76","RSPO1","SEMA6A","ZHX2","FLT1","LYN","KCNH2",
                          "GNB4","PDC","RGS16","GUCA1B","ROM1","GUCA1C","TRIM59"), 
                        TERM2GENE = call[,c(1,2)],pvalueCutoff = 0.05,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0)


a <- data.frame(enrich_C1uC2h)




temp <- enrich_C1u
temp@result <- temp@result[grep("^GO",temp@result$ID),]
temp@result <- merge(temp@result,c_go_correlation,by="ID")
temp@result <- temp@result[order(temp@result$p.adjust,decreasing = F),]
temp@result$Description <- gsub(pattern = "GO_",replacement = "",temp@result$Description)
dotplot(temp,title='Top 5 terms of C1 upregulated DEGs', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")

temp <- enrich_C1u
temp@result <- temp@result[grep("^KEGG|^HALLMARK|^WP",temp@result$ID),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "KEGG_|HALLMARK_|WP_",replacement = "",temp@result$Description)
dotplot(temp,title='Top 5 terms of C1 upregulated DEGs', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")



temp <- enrich_C2u
temp@result <- temp@result[grep("^GO",temp@result$ID),]
temp@result <- merge(temp@result,c_go_correlation,by="ID")
temp@result <- temp@result[order(temp@result$p.adjust,decreasing = F),]
temp@result$Description <- gsub(pattern = "GO_",replacement = "",temp@result$Description)
dotplot(temp,title='Top 5 terms of C2 upregulated DEGs', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")

temp <- enrich_C2u
temp@result <- temp@result[grep("^KEGG|^HALLMARK|^WP",temp@result$ID),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "KEGG_|HALLMARK_|WP_",replacement = "",temp@result$Description)
dotplot(temp,title='Top 5 terms of C2 upregulated DEGs', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")



temp <- enrich_C1h
temp@result <- temp@result[grep("^GO",temp@result$ID),]
temp@result <- merge(temp@result,c_go_correlation,by="ID")
temp@result <- temp@result[order(temp@result$p.adjust,decreasing = F),]
temp@result$Description <- gsub(pattern = "GO_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 GO terms of C1 hypermethylated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))

temp <- enrich_C1h
temp@result <- temp@result[grep("^KEGG|^HALLMARK|^REACTOME",temp@result$ID),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "KEGG_|HALLMARK_|REACTOME_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 terms of C1 hypermethylated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))

temp <- enrich_C1h
temp@result <- temp@result[temp@result$ID %in% c("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
                                                 "GO_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
                                                 "GO_T_CELL_DIFFERENTIATION",
                                                 "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                                                 "GO_T_CELL_ACTIVATION",
                                                 "WP_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
                                                 "WP_DEVELOPMENT_OF_PULMONARY_DENDRITIC_CELLS_AND_MACROPHAGE_SUBSETS",
                                                 "GO_DENDRITIC_CELL_CHEMOTAXIS",
                                                 "GO_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION",
                                                 "GO_DENDRITIC_CELL_DIFFERENTIATION",
                                                 "GO_DENDRITIC_CELL_MIGRATION"),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "KEGG_|HALLMARK_|REACTOME_|WP_|REACTOME_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Terms related to immunologi', showCategory=10) +
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))



temp <- enrich_C2h
temp@result <- temp@result[grep("^GO",temp@result$ID),]
temp@result <- merge(temp@result,c_go_correlation,by="ID")
temp@result <- temp@result[order(temp@result$p.adjust,decreasing = F),]
temp@result$Description <- gsub(pattern = "GO_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 GO terms of C2 hypermethylated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))

temp <- enrich_C2h
temp@result <- temp@result[grep("^KEGG|^REACTOME|^WP",temp@result$ID),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "KEGG_|REACTOME_|WP_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 terms of C2 hypermethylated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))






temp <- enrich_C1uC2h
temp@result <- temp@result[grep("^GO",temp@result$ID),]
temp@result <- merge(temp@result,c_go_correlation,by="ID")
temp@result <- temp@result[order(temp@result$p.adjust,decreasing = F),]
temp@result$Description <- gsub(pattern = "GO_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 GO terms of C1 hypomethylated upregulated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))

temp <- enrich_C1uC2h
temp@result <- temp@result[grep("^REACTOME|^PID",temp@result$ID),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "	REACTOME_|PID_",replacement = "",temp@result$Description)
dotplot(temp,title='Top 5 terms of C1 hypomethylated upregulated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))



temp <- enrich_C1hC2u
temp@result <- temp@result[grep("^GO",temp@result$ID),]
temp@result <- merge(temp@result,c_go_correlation,by="ID")
temp@result <- temp@result[order(temp@result$p.adjust,decreasing = F),]
temp@result$Description <- gsub(pattern = "GO_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 GO terms of C2 hypomethylated upregulated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))

temp <- enrich_C1hC2u
temp@result <- temp@result[grep("^REACTOME|^PID",temp@result$ID),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "	REACTOME_|PID_",replacement = "",temp@result$Description)
dotplot(temp,title='Top 5 terms of C2 hypomethylated upregulated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))







a <- data.frame(enrich_G3)

temp <- enrich_G1
temp@result <- temp@result[grep("^KEGG|^REACTOME|^HALLMARK",temp@result$ID),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "REACTOME_|KEGG_|HALLMARK_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 terms of gene cluster 1', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))



temp <- enrich_G2
temp@result <- temp@result[grep("^KEGG|^REACTOME|^HALLMARK",temp@result$ID),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "REACTOME_|KEGG_|HALLMARK_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 terms of gene cluster 2', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))




temp <- enrich_G3
temp@result <- temp@result[grep("^PID|^REACTOME|^HALLMARK",temp@result$ID),]
temp@result <- temp@result[grep("REACTOME_COOPERATION_OF_PDCL_PHLP1_AND_TRIC_CCT_IN_G_PROTEIN_BETA_FOLDING",temp@result$Description,invert = T),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "REACTOME_|PID_|HALLMARK_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 terms of gene cluster 3', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free") +
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))


temp <- enrich_temp
temp@result <- temp@result[grep("^GO",temp@result$ID),]
temp@result <- merge(temp@result,c_go_correlation,by="ID")
temp@result <- temp@result[order(temp@result$p.adjust,decreasing = F),]
temp@result$Description <- gsub(pattern = "GO_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 GO terms of C1 hypermethylated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))

temp <- enrich_C1h
temp@result <- temp@result[grep("^KEGG|^HALLMARK|^REACTOME",temp@result$ID),]
temp@result$ONTOLOGY <- temp@result$ID %>% strsplit("_") %>% lapply(function(x) x[1])%>% unlist()
temp@result$Description <- gsub(pattern = "KEGG_|HALLMARK_|REACTOME_",replacement = "",temp@result$Description)
temp@result$Description <- gsub(pattern = "_",replacement = " ",temp@result$Description) %>% tolower()
dotplot(temp,title='Top 5 terms of C1 hypermethylated genes', showCategory=5,split='ONTOLOGY') +
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(axis.text.x = element_text(face = "plain",size = 13,color="black"),
        axis.text.y = element_text(face = "plain",size = 13,color="black"))

