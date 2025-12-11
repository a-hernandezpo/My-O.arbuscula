# nolint start: line_length_linter


## Volcano plot con TODA la información solicitada:

# Contar genes en cada categoría
up_count <- sum(df_volc$threshold == "Up")
down_count <- sum(df_volc$threshold == "Down")
ns_count <- sum(df_volc$threshold == "NS")

# Preparar datos para etiquetar TODOS los genes significativos
# Genes UP (rojos) - todos
genes_up <- df_volc %>% 
  filter(threshold == "Up") %>%
  mutate(label_color = "red")

# Genes DOWN (azules) - todos
genes_down <- df_volc %>% 
  filter(threshold == "Down") %>%
  mutate(label_color = "blue")

# Combinar todos los genes significativos
genes_sig <- bind_rows(genes_up, genes_down)

# Crear el volcano plot con TODA la información
p <- ggplot(df_volc, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c("Down" = "blue", "Up" = "red", "NS" = "grey"),
    name = "Tipo de regulación:",
    labels = c(
      paste0("Down (", down_count, " genes)"),
      paste0("NS (", ns_count, " genes)"),
      paste0("Up (", up_count, " genes)")
    )
  ) +

  # Agregar etiquetas para TODOS los genes UP (rojos)
  geom_text_repel(
    data = genes_up,
    aes(label = genename),
    color = "darkred",
    size = 2.8,
    max.overlaps = 100,  # Aumentar para mostrar más etiquetas
    box.padding = 0.35,
    point.padding = 0.3,
    segment.color = "darkred",
    segment.alpha = 0.5,
    min.segment.length = 0.1,
    force = 2
  ) +
  
  # Agregar etiquetas para TODOS los genes DOWN (azules)
  geom_text_repel(
    data = genes_down,
    aes(label = genename),
    color = "darkblue",
    size = 2.8,
    max.overlaps = 100,
    box.padding = 0.35,
    point.padding = 0.3,
    segment.color = "darkblue",
    segment.alpha = 0.5,
    min.segment.length = 0.1,
    force = 2
  ) +
  
  # Líneas de corte finas
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black", alpha = 0.5) +
  
  # Etiquetas de corte FUERA del gráfico
  annotate("text", x = -2, y = 0, label = "-2", vjust = 2, size = 5) +
  annotate("text", x = 2, y = 0, label = "2", vjust = 2, size = 5) +
  annotate("text", x = min(df_volc$log2FoldChange), 
           y = -log10(0.1), label = "0.1", hjust = -0.2, vjust = -1, size = 5) +
  
  # Anotaciones con conteos en el gráfico
  annotate("text", 
           x = min(df_volc$log2FoldChange) * 0.95, 
           y = max(-log10(df_volc$padj)) * 0.98,
           label = paste("DOWN:", down_count, "genes"),
           color = "blue", size = 4.5, fontface = "bold", hjust = 0) +
  
  annotate("text", 
           x = max(df_volc$log2FoldChange) * 0.95, 
           y = max(-log10(df_volc$padj)) * 0.98,
           label = paste("UP:", up_count, "genes"),
           color = "red", size = 4.5, fontface = "bold", hjust = 1) +
  
  # Título con resumen completo
  labs(
    #title = "GS Volcano plot: Sym vs Apo",
    subtitle = paste("Genes totales: ", nrow(df_volc), 
                     " | Significativos: ", up_count + down_count,
                     " (", round((up_count + down_count)/nrow(df_volc)*100, 1), "%)", sep = ""),
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    #plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5, margin = margin(b = 15), position = "bottom"),
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    #panel.grid.minor = element_blank(),
    
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),  # tamaño título eje X
    axis.title.y = element_text(size = 20, face = "bold"),  # tamaño título eje Y

  ) +
  
  # Ajustar límites para mejor visualización
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1))

# Mostrar el gráfico
print(p)

# Guardar con alta resolución (ajustar tamaño según número de genes)
width_size <- 12 + (nrow(genes_sig) / 50)  # Ajustar ancho según número de genes
height_size <- 10 + (nrow(genes_sig) / 100)  # Ajustar alto según número de genes

ggsave("./GS figures/GS Volcano_Sym_vs_Apo_COMPLETO.png", 
       plot = p, 
       width = min(width_size, 20),  # Máximo 20 pulgadas
       height = min(height_size, 16),  # Máximo 16 pulgadas
       dpi = 300,
       bg = "white")


# Opción: Guardar también versión PDF para mejor escalado
ggsave("./GS figures/GS Volcano_Sym_vs_Apo_COMPLETO.pdf", 
       plot = p, 
       width = 16, 
       height = 12,
       device = cairo_pdf)

# Mostrar resumen detallado en consola

cat("Total genes analizados:", nrow(df_volc), "\n") #25428
cat("Genes UP-regulados (rojos):", up_count, "\n") #71
cat("Genes DOWN-regulados (azules):", down_count, "\n") #7
cat("Genes no significativos (gris):", ns_count, "\n") #25350
cat("Porcentaje significativos:", round((up_count + down_count)/nrow(df_volc)*100, 1), "%\n") #0.3%

if (up_count > 0) {
  cat("\nTop 5 genes UP-regulados:\n")
  print(head(df_volc[df_volc$threshold == "Up", ] %>% 
               arrange(padj) %>% 
               select(gene, log2FoldChange, padj), 5))
}

if (down_count > 0) {
  cat("\nTop 5 genes DOWN-regulados:\n")
  print(head(df_volc[df_volc$threshold == "Down", ] %>% 
               arrange(padj) %>% 
               select(gene, log2FoldChange, padj), 5))
}





ccol<-rev(colorRampPalette(brewer.pal(n=11, name="BrBG"))(80))
ccol2<-rev(colorRampPalette(brewer.pal(n=11, name="RdYlBu"))(50))

pheatmap(degs_plot[,1:10], show_rownames=T, labels_row=degs_plot$genename,
         show_colnams = T, angle_col= "0",
         cellheight = 6, cellwidth = 8, fontsize_row=5, cluster_cols=T, scale='row', fontsize_col=8,
         cluster_rows = T, color=ccol2, cutree_rows = 2, cutree_cols = 2, labels_col = genets,
         treeheight_col = 10, treeheight_row = 10, legend=T,
         filename="./GS figures/heatmap_DEGs_annotonly.png", width=8, height=15)

## heat map of only glutamine-associated DEGs 

glut_degs<-droplevels(as.data.frame(degs[grepl("glut", degs$genename, ignore.case = TRUE),]))
pheatmap(glut_degs[,1:10], 
         show_rownames=T, labels_row=glut_degs$genename,fontsize_row=5, scale='row',cluster_rows = T, 
         show_colnams = T, angle_col= "0",cluster_cols=T, fontsize_col=8, labels_col = genets,
         cellheight = 6, cellwidth = 8,  
         color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 5, treeheight_row = 5, legend=T,
         filename="./GS figures/glut_DEGs.png", width=5, height=5)



## heat map of only acly-coA and other nutrient transporter-assocociated DEGs 
sug_lip_degs<-droplevels(as.data.frame(degs[grepl("acyl-|acyltrans|fatty|solute|lipid|cholesterol|apolipo|myo-inositol", degs$genename, ignore.case = TRUE),]))
pheatmap(sug_lip_degs[,1:10], 
         show_rownames=T, labels_row=sug_lip_degs$genename,fontsize_row=5, scale='row',cluster_rows = T, 
         show_colnams = T, angle_col= "0",cluster_cols=T, fontsize_col=8, labels_col = genets,
         cellheight = 6, cellwidth = 8,  
         color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 5, treeheight_row = 5, legend=T,
         filename="./GS figures/sug_lip_DEGs.png", width=5, height=5)


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
ggsave("./GS figures/Fig2_new.png", dpi=300, height=6, width=10, units="in")

# nolint end: line_length_linter