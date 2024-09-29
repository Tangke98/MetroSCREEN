
#' Do the heatmap to show the difference.

#' @param object A matrix, the row is the MetaModule and the column is the sample/metacell, the value is the MetaModule score.
#' @param sample_info A vector including ample information for the sample/metacell.
#' @param feature Feature to show.
#' @param width Width of the figure.
#' @param height Height of the figure.
#' @export

doheatmap_feature=function(object,sample_info,feature,width,height,cols){
    cluster_info <- sort(sample_info)
    filtered_feature <- feature[feature %in% rownames(object)]
    mat <-object[match(filtered_feature, rownames(object)), names(cluster_info)]
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
        column_title = NULL,col=c('#E0F3F8','white','#DE77AE'),
        heatmap_legend_param = list(legend_direction = "horizontal",title = "MetaModule score")
    )
    draw(p, heatmap_legend_side = "bottom")
}


#' Find the differentially enriched MetaModule.

#' @param object A matrix, the row is the MetaModule and the column is the sample/metacell, the value is the MetaModule score.
#' @param sample_info A vector including ample information for the sample/metacell.
#' @param sample_tech scRNA or bulk.
#' @export

FindAllMarkers_MetaModule=function(object,sample_info,sample_tech){
    if (sample_tech=='scRNA'){
        object.seurat <- CreateSeuratObject(counts = object, project = "metacell", min.cells = 0, min.features = 0)
        object.seurat$state=sample_info
        object.seurat@active.ident=as.factor(object.seurat$state)
        MetaModule.markers <- FindAllMarkers(object.seurat, only.pos = TRUE) 
    }
    if (sample_tech=='bulk'){
        object.seurat <- CreateSeuratObject(counts = object, project = "metacell", min.cells = 0, min.features = 0)
        object.seurat$state=sample_info
        object.seurat@active.ident=as.factor(object.seurat$state)
        MetaModule.markers <- FindAllMarkers(object.seurat, only.pos = TRUE,test.use='t') 
    }
    return(MetaModule.markers)
}
