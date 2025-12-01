# nolint start: line_length_linter

ccol<-rev(colorRampPalette(brewer.pal(n=11, name="BrBG"))(80))
ccol2<-rev(colorRampPalette(brewer.pal(n=11, name="RdYlBu"))(50))

pheatmap(degs_plot[,1:10], show_rownames=T, labels_row=degs_plot$genename,
         show_colnams = T, angle_col= "0",
         cellheight = 6, cellwidth = 8, fontsize_row=5, cluster_cols=T, scale='row', fontsize_col=8,
         cluster_rows = T, color=ccol2, cutree_rows = 2, cutree_cols = 2, labels_col = genets,
         treeheight_col = 10, treeheight_row = 10, legend=T,
         filename="./figures/heatmap_DEGs_annotonly.png", width=8, height=15)

## heat map of only glutamine-associated DEGs 

glut_degs<-droplevels(as.data.frame(degs[grepl("glut", degs$genename, ignore.case = TRUE),]))
pheatmap(glut_degs[,1:10], 
         show_rownames=T, labels_row=glut_degs$genename,fontsize_row=5, scale='row',cluster_rows = T, 
         show_colnams = T, angle_col= "0",cluster_cols=T, fontsize_col=8, labels_col = genets,
         cellheight = 6, cellwidth = 8,  
         color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 5, treeheight_row = 5, legend=T,
         filename="./figures/glut_DEGs.png", width=5, height=5)



## heat map of only acly-coA and other nutrient transporter-assocociated DEGs 
sug_lip_degs<-droplevels(as.data.frame(degs[grepl("acyl-|acyltrans|fatty|solute|lipid|cholesterol|apolipo|myo-inositol", degs$genename, ignore.case = TRUE),]))
pheatmap(sug_lip_degs[,1:10], 
         show_rownames=T, labels_row=sug_lip_degs$genename,fontsize_row=5, scale='row',cluster_rows = T, 
         show_colnams = T, angle_col= "0",cluster_cols=T, fontsize_col=8, labels_col = genets,
         cellheight = 6, cellwidth = 8,  
         color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 5, treeheight_row = 5, legend=T,
         filename="./figures/sug_lip_DEGs.png", width=5, height=5)


# PCA's = Evaluar estadísticamente si Genet y SymbiontState explican la variación entre muestras observada en el PCA.

pca = prcomp(t(assay(rlogMF)), center = TRUE, scale. = FALSE)
adonis2(pca$x~Genet+SymbiontState,data=coldata, method = 'eu')

fullpca<-as.data.frame(pca$x)
fullpca$Sample<-row.names(fullpca)
fullpca<-cbind(fullpca, coldata)

(a<-ggplot(fullpca, aes(PC1, PC2, shape=Genet)) +
  geom_hline(yintercept =0)+geom_vline(xintercept =0) +
  geom_point(aes(colour=SymbiontState,fill=SymbiontState, stroke=1.1), size=7)+
  xlab(paste0("PC1 (44.1%)")) + # from summary(pca)
  ylab(paste0("PC2 (15.4%)")) + 
  theme_bw()+
  scale_shape_manual(values=c(21,22,23,24,25), guide=guide_legend(override.aes = list(size = 4)))+
  scale_fill_manual(values=c("#80CDC1", "#BF812D"))+
  scale_colour_manual(values=c("#003C30", "#543005"),
                      labels=c("Aposymbiotic", "Symbiotic"), 
                      guide=guide_legend(override.aes = list(size = 4,colour = c("#80CDC1", "#BF812D"))))+
  annotate("text", x=30, y=40, label=paste("italic(Adonis) * ' p'[Genet]<0.01"), parse=TRUE)+
  annotate("text", x=30, y=45, label=paste("italic(Adonis) * ' p'[Branch]>0.2"), parse=TRUE)+
  theme(panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour="black", fill=NA, size=1.5), 
        axis.text = element_text(face="bold", colour="black", size=8),
        legend.position = "bottom",
        legend.background = element_rect(colour="grey"),
        legend.text = element_text(face="bold", colour="black", size=8),
        axis.title = element_text(face="bold", colour = "black", size=10),
        legend.title = element_text(face="bold", colour="black", size=10))+
  labs(shape="Colony", colour="Branch")+ guides(fill=FALSE))


# PCA of DEGs
rlogMF_deg<-rlogMF[rownames(conds),]
pca_degs = prcomp(t(assay(rlogMF_deg)), center = TRUE, scale. = FALSE)

adonis2(pca_degs$x~Genet+SymbiontState,data=coldata,  method = 'eu')

pca_degs_df<-as.data.frame(pca_degs$x)
pca_degs_df$Sample<-row.names(pca_degs_df)
pca_degs_df<-cbind(pca_degs_df, coldata)


(b<-ggplot(pca_degs_df, aes(PC1, PC2)) +
  geom_hline(yintercept =0)+geom_vline(xintercept =0) +
  geom_point(aes(color=SymbiontState, fill=SymbiontState, shape=Genet, stroke=1.1), size=7)+
  stat_ellipse(geom="polygon", aes(fill=SymbiontState), alpha=0.3)+
  xlab(paste0("PC1 (63.5%)")) + # from summary(pca_degs)
  ylab(paste0("PC2 (16.5%)")) +
  theme_minimal()+
  scale_shape_manual(values=c(21,22,23,24,25))+
  scale_fill_manual(values=c("#80CDC1", "#BF812D"))+
  scale_colour_manual(values=c("#003C30", "#543005"),
                      labels=c("Aposymbiotic", "Symbiotic"), 
                      guide=guide_legend(override.aes = list(size = 4,colour = c("#80CDC1", "#BF812D"))))+
  annotate("text", x=13, y=13.5, label=paste("italic(Adonis) * ' p'[Genet]<0.01"), parse=TRUE)+
  annotate("text", x=13, y=12, label=paste("italic(Adonis) * ' p'[Branch]<0.01"), parse=TRUE)+
  theme(panel.grid.minor=element_blank(),
        panel.border = element_rect(colour="black", fill=NA, size=1.5),
        axis.text = element_text(face="bold", colour="black", size=8),
        legend.position = "bottom",
        legend.background = element_rect(colour="grey"),
        legend.text = element_text(face="bold", colour="black", size=8),
        axis.title = element_text(face="bold", colour = "black", size=10),
        legend.title = element_text(face="bold", colour="black", size=10))+
  labs(shape="Colony", colour="Branch")+guides(fill=FALSE)+
  guides(shape = guide_legend(override.aes = list(size = 4)),
         colour = guide_legend(override.aes = list(size = 4))))


ggarrange(a, b, ncol=2, nrow=1, common.legend = TRUE, legend = "bottom", labels="AUTO")+
  theme(plot.margin = margin(0,0,2,0))
ggsave("./figures/Fig2_new.png", dpi=300, height=6, width=10, units="in")

# nolint end: line_length_linter