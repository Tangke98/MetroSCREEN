
#' Do the heatmap to show the difference.

#' @param object A seurat object.
#' @param cluster Info to cluster the cells.
#' @param feature Feature to show.
#' @param width Width of the figure.
#' @param height Height of the figure.
#  @export

doheatmap_feature=function(object,cluster,feature,width,height,cols){
    mat <- as.matrix(GetAssayData(object, slot = "counts"))
    object$cell_type=as.factor(object@meta.data[,cluster])
    cluster_info <- sort(object$cell_type)
    filtered_feature <- feature[feature %in% rownames(mat)]
    mat <-mat[match(filtered_feature, rownames(mat)), names(cluster_info)]
    # col <- brewer.pal(12, "Set3")
    # col_use=col[1:length(levels(cluster_info))]
    # names(col_use) <- levels(cluster_info)
     top_anno=HeatmapAnnotation(Type=cluster_info,
            col=list(Type=cols))
    col_heatmap <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)
    options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
    p=Heatmap(as.matrix(mat),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_split = cluster_info,
        top_annotation = top_anno,
        column_title = NULL,col=c('#E0F3F8','white','#DE77AE')
       ) 
    print(p)
}

