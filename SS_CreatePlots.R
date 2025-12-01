# nolint start: line_length_linter

ccol<-rev(colorRampPalette(brewer.pal(n=11, name="BrBG"))(80))
ccol2<-rev(colorRampPalette(brewer.pal(n=11, name="RdYlBu"))(50))

pheatmap(degss_plot[,1:10], show_rownames=T, labels_row=degss_plot$genename,
         show_colnams = T, angle_col= "0",
         cellheight = 6, cellwidth = 8, fontsize_row=5, cluster_cols=T, scale='row', fontsize_col=8,
         cluster_rows = T, color=ccol2, cutree_rows = 2, cutree_cols = 2, labels_col = genetss,
         treeheight_col = 10, treeheight_row = 10, legend=T,
         filename="./SS figures/SS heatmap_DEGs_annotonly.png", width=8, height=15)

## heat map of only glutamine-associated DEGs 

glut_degss<-droplevels(as.data.frame(degss[grepl("glut", degss$genename, ignore.case = TRUE),]))
pheatmap(glut_degss[,1:10], 
         show_rownames=T, labels_row=glut_degss$genename,fontsize_row=5, scale='row',cluster_rows = T, 
         show_colnams = T, angle_col= "0",cluster_cols=T, fontsize_col=8, labels_col = genetss,
         cellheight = 6, cellwidth = 8,  
         color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 5, treeheight_row = 5, legend=T,
         filename="./SS figures/SS glut_DEGs.png", width=5, height=5)



## heat map of only acly-coA and other nutrient transporter-assocociated DEGs 
sug_lip_degss<-droplevels(as.data.frame(degss[grepl("acyl-|acyltrans|fatty|solute|lipid|cholesterol|apolipo|myo-inositol", degss$genename, ignore.case = TRUE),]))
pheatmap(sug_lip_degss[,1:10], 
         show_rownames=T, labels_row=sug_lip_degss$genename,fontsize_row=5, scale='row',cluster_rows = T, 
         show_colnams = T, angle_col= "0",cluster_cols=T, fontsize_col=8, labels_col = genetss,
         cellheight = 6, cellwidth = 8,  
         color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 5, treeheight_row = 5, legend=T,
         filename="./SS figures/SS sug_lip_DEGs.png", width=5, height=5)

# nolint end: line_length_linter