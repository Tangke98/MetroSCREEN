#' Calculate the MetaRegulon score for each MetaModule in cell subtypes/metabolic subgtypes/sample conditions.
#' 
#' @param object A Seurat object
#' @param feature The feature used to cluster cells, typically chosen from the column names of the meta.data in the Seurat object 
#' @param state The cell subtype, metabolic subtype, or sample condition
#' @param interested_MM The name of the metabolic reaction or interested gene set
#' @param MM_list A list of metabolic reactions or interested gene sets
#' @param markers Marker genes enriched in the corresponding state
#' @param lisa_file The path to the LISA results for the marker genes
#' @param ligand_target_matrix The path to the ligand-target matrix (provided)
#' @param lr_network The ligand-receptor network matrix (provided)
#' @param sample_tech The technology used for sequencing
#' @param output_path The path to save the MetaRegulon results
#' @param RP_path The RP results (provided)
#' @param file_name The file name for the MetaRegulon results
#' 
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import nichenetr
#' @import ggpubr
#' @import GSVA
#' @import ComplexHeatmap
#' @import limma
#' @import RColorBrewer
#' @import rPref
#' @import kpcalg
#' @import igraph
#' @import RobustRankAggreg
#' @import parallel
#' @import ggraph
#' @import reshape2
#' 
#' @export

cal_MetaRegulon=function(object,feature,state,interested_MM,MM_list,markers,lisa_file,ligand_target_matrix,lr_network,sample_tech,output_path,RP_path,file_name){
    RP_score=readRDS(RP_path)
    tmp_dir <- file.path(output_path,file_name)
    if (!dir.exists(tmp_dir)) {
        dir.create(tmp_dir)
    }
    # if (sample_tech=='scRNA'){
    if (length(object@assays)>1) {
    # if 'integrated' exists，set it as the default
        DefaultAssay(object) <- "integrated"
        data=GetAssayData(object)
        data_use=data[,rownames(object@meta.data[object[[feature]]==state,])]
    } else {
        data=GetAssayData(object)
        data_use=data[,rownames(object@meta.data[object[[feature]]==state,])]
    }
    # }
    # if (sample_tech=='bulk'){
    #     data_use=object
    # }

    file_res=paste0(tmp_dir,'/',file_name,".rds") ## get the state data to do the following analysis
    if (!file.exists(file_res)){
        saveRDS(data_use,file_res)
    }
    cl <- makeCluster(3) ## only three functions,set to 3
    clusterExport(cl,
                  list("output_path", "file_name", "state", "file_res", "markers",'lisa_file',"MM_list", "interested_MM",'RP_score','ligand_target_matrix','lr_network',
                       "get_cor",'cal_tf_score','cal_target_score',
                       "gene_gene_interaction", "ligand_target_interaction",'tr_target_interaction'),envir = environment())
    clusterEvalQ(cl, {
        suppressMessages({
            libraries <- c("Seurat", "dplyr", "tidyr", "ggplot2",
                        "nichenetr", "ggpubr", "GSVA", 
                        "ComplexHeatmap", "limma", "RColorBrewer",
                        "rPref", "igraph", "RobustRankAggreg", "parallel",
                        "MetroSCREEN", "reshape2")
            lapply(libraries, library, character.only = TRUE)
            # source('/fs/home/tangke/metabolism/tool/MetroSCREEN/R/kpc.r')
        })
    })
    tasks <- list(gene_gene_task = list(fun = 'gene_gene_interaction', args = list(output_path,file_name,file_res,MM_list,interested_MM)),
                  ligand_target_task = list(fun = 'ligand_target_interaction', args = list(output_path,file_name,state,file_res,markers,MM_list,interested_MM,ligand_target_matrix,lr_network)),
                  tr_target_task = list(fun = 'tr_target_interaction', args = list(output_path,file_name,lisa_file,file_res,RP_score,MM_list,interested_MM))
                 )
    results <- parLapply(cl, tasks, function(task) {
        tryCatch({
            do.call(task$fun, task$args)
            return(paste("Task completed successfully: ", task$fun, "; Result: ", capture.output(result)))
        }, error = function(e) {
            return(paste("Error in function", task$fun, ":", e$message))
        })
    })

    stopCluster(cl) 
    
    tr_cor=readRDS(paste0(output_path,file_name,"/",file_name,":tr_activity_cor.rds"))
    ligand_cor=readRDS(paste0(output_path,file_name,"/",file_name,":lr_activity_cor.rds"))
    gene_cor=readRDS(paste0(output_path,file_name,"/",file_name,":gg_activity_cor.rds"))
    tr_activity=readRDS(paste0(output_path,file_name,"/",file_name,":tr_activity.rds"))
    lr_activity=readRDS(paste0(output_path,file_name,"/",file_name,":lr_activity.rds"))

    data=t(as.matrix(readRDS(file_res)))
    repref_each_res=lapply(as.list(interested_MM),repref_each,tr_cor=tr_cor,tr_activity=tr_activity,lr_activity=lr_activity,ligand_cor=ligand_cor,gene_cor=gene_cor,MM_list=MM_list,data=data,output_path=output_path,file_name=file_name)
}

#' cal_tf_score
#' 
#' @param tf 
#' @param data 
#' 
#' @export

cal_tf_score=function(tf,data){
    tf_orig=data[rownames(data) %in% tf,] ## lisa res 
    tf_orig_zscore=apply(tf_orig,1,function(x){(x-mean(x))/sd(x)}) ## normalize
    cols_with_nan <- apply(tf_orig_zscore, 2, function(x) any(is.nan(x))) ## remove na
    tf_orig_zscore_filtered <- tf_orig_zscore[, !cols_with_nan]
    tf_orig_zscore_filtered[tf_orig_zscore_filtered>4]=4 
    tf_orig_zscore_filtered[tf_orig_zscore_filtered<(-4)]=(-4)            
    return(tf_orig_zscore_filtered)
}

#' cal_target_score
#' 
#' @param lisa_res_in_metroscreen A seurat object
#' @param factor Feature used to cluster cells, usually use the colnames of the meta.data of the seurat object 
#' @param data State of the interested cell types pr metabolic state
#' @param RP_score
#' 
#' @export

cal_target_score=function(lisa_res_in_metroscreen,factor,data,RP_score){
    sample_id=lisa_res_in_metroscreen[lisa_res_in_metroscreen$factor==factor,'sample_id'][1] 
    if (nchar(sample_id)>2){
        dir_id=substring(sample_id,1,3)
    }else{
        dir_id='000'
    }
    
    tf_target=RP_score[,c('gene',paste0(dir_id,'.',sample_id,'_gene_score.txt'))]
    target=tf_target[tf_target[,2]>2,'gene']
    
    if (length(target)<500){
        tf_target_use_order=tf_target[order(tf_target[,2],decreasing = TRUE),]
        tf_target_use_order=tf_target_use_order[tf_target_use_order[,2]>0,]

        if (nrow(tf_target_use_order>=500)){
            target=tf_target_use_order[1:500,'gene']
        }else{
            target=tf_target_use_order[,'gene']
        }
    }
    
    target_orig=data[rownames(data) %in% unique(target),] 
    target_orig_zscore=apply(target_orig,1,function(x){(x-mean(x))/sd(x)})
    cols_with_nan <- apply(target_orig_zscore, 2, function(x) any(is.nan(x)))
    target_orig_zscore_filtered <- target_orig_zscore[, !cols_with_nan] 
    target_orig_zscore_filtered_mean=apply(target_orig_zscore_filtered,1,mean)                       
    target_orig_mean=as.data.frame(target_orig_zscore_filtered_mean)
                           
    target_orig_mean[target_orig_mean>4]=4 
    target_orig_mean[target_orig_mean<(-4)]=(-4)      
                           
    colnames(target_orig_mean)=factor
    return(target_orig_mean)
}

#' get_cor
#' 
#' @param cor
#' @param MM_list 
#' @param MM_index 
#' @param class
#' 
#' @export

get_cor=function(cor,MM_list,MM_index,class){
    cor=cor[,colnames(cor) %in% MM_list[[MM_index]]]
    if(is.matrix(cor)){
        cor <- cor[, colSums(is.na(cor)) != nrow(cor)]
        if(is.matrix(cor)){
            cor_max=as.data.frame(apply(na.omit(cor),1,max))
            cor_max$MR=rownames(cor_max)
        }
        else{
            cor_max=as.data.frame(na.omit(cor))
            cor_max$MR=rownames(cor_max)
        }
    } else{
        cor_max=as.data.frame(na.omit(cor))
        cor_max$MR=rownames(cor_max)
    }
    colnames(cor_max)[1] = switch(class,
                                   'gene' = 'Gene_Gene_interaction',
                                   'ligand' = 'Ligand_Receptor_interaction',
                                   'tr' = 'TR_Target_interaction')
    cor_max = cor_max[!rownames(cor_max) %in% MM_list[[MM_index]],]
    threshold = ifelse(class == 'gene', 0.3, 0)
    cor_max = cor_max[cor_max[[1]] > threshold, ]
    return(cor_max)
    
}

#' gene_gene_interaction
#' 
#' @param output_path 
#' @param file_name 
#' @param file_res 
#' @param MM_list
#' @param interested_MM
#' 
#' @export

gene_gene_interaction=function(output_path,file_name,file_res,MM_list,interested_MM){
    
    ## get the genome-wide correlation
    gg_cor_res=paste0(output_path,file_name,"/",file_name,":gg_activity_cor.rds")
    if (!file.exists(gg_cor_res)){
        data=readRDS(file_res)
        data=t(as.matrix(data))
        result=data %>%
            apply(2,function(x){10*(2**x - 1)}) %>%
            apply(2,function(x){log2(mean(x) + 1)}) 
        expressed_genes <- names(result[order(result, decreasing = TRUE)][1:3000]) ## screen the expressed genes
        data_exp=data[,expressed_genes]
        data_mg=data[rownames(data_exp),colnames(data) %in% unique(unlist(MM_list))]
        cor=cor(data_exp,data_mg)
        get_cor_res=lapply(as.list(interested_MM),get_cor,cor=cor,MM_list=MM_list,class='gene')
        names(get_cor_res)=interested_MM
        saveRDS(get_cor_res,gg_cor_res)    
    }
}

#' ligand_target_interaction
#' 
#' @param output_path 
#' @param file_name 
#' @param state 
#' @param file_res
#' @param markers
#' @param MM_list
#' @param interested_MM
#' @param ligand_target_matrix
#' @param lr_network
#' 
#' @export

ligand_target_interaction=function(output_path,file_name,state,file_res,markers,MM_list,interested_MM,ligand_target_matrix,lr_network){
    
    ## get the ligand activity
    ligand_target_matrix = readRDS(ligand_target_matrix)
    lr_network = readRDS(lr_network)
    
    lr_activity_res=paste0(output_path,file_name,"/",file_name,":lr_activity.rds")
    if (!file.exists(lr_activity_res)){
        marker_state=markers[markers$cluster==state,'gene']
        data=t(as.matrix(readRDS(file_res)))
        data_use=data[,intersect(colnames(data),marker_state)]
        expressed_genes <-colnames(data_use)
        
        ligands = lr_network$from %>% unique()
        expressed_ligands = intersect(ligands,unique(markers$gene))
        receptors = lr_network$to %>% unique()
        expressed_receptors = intersect(receptors,marker_state) 
        
        potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% .$from %>% unique()
        background_expressed_genes = expressed_genes %>% .[. %in% rownames(ligand_target_matrix)]
        
        expression_scaled = data_use %>% .[,background_expressed_genes] %>% scale_quantile()
        
        ligand_activities = predict_single_cell_ligand_activities(cell_ids = rownames(data_use), 
                                                            expression_scaled = expression_scaled, 
                                                            ligand_target_matrix = ligand_target_matrix, 
                                                            potential_ligands = potential_ligands)
        normalized_ligand_activities = as.data.frame(normalize_single_cell_ligand_activities(ligand_activities))
        rownames(normalized_ligand_activities)=normalized_ligand_activities[,1]
        normalized_ligand_activities=normalized_ligand_activities[,-1]
        saveRDS(normalized_ligand_activities,lr_activity_res)
    }

    ## get the ligand-receptor interaction
    lr_cor_res=paste0(output_path,file_name,"/",file_name,":lr_activity_cor.rds")
    if (!file.exists(lr_cor_res)){
        data=readRDS(file_res)
        data_lr=readRDS(lr_activity_res)
        data_mg=t(as.matrix(data[rownames(data) %in% unique(unlist(MM_list)),rownames(data_lr)]))
        cor=cor(data_lr,data_mg)
        
        get_cor_res=lapply(as.list(interested_MM),get_cor,cor=cor,MM_list=MM_list,class='ligand')
        names(get_cor_res)=interested_MM
        saveRDS(get_cor_res,lr_cor_res)
    }
}

#' tr_target_interaction
#' 
#' @param output_path
#' @param file_name
#' @param lisa_file 
#' @param file_res
#' @param RP_score
#' @param MM_list
#' @param interested_MM
#' 
#' @export

tr_target_interaction=function(output_path,file_name,lisa_file,file_res,RP_score,MM_list,interested_MM){
    cat("Starting tr_target_interaction\n")
    ## get the tr activity
    tr_activity_res=paste0(output_path,file_name,"/",file_name,":tr_activity.rds")
    lisa_res_in_metroscreen <- read.csv(paste0(lisa_file),sep='\t')
    if (!file.exists(tr_activity_res)){
        data=readRDS(file_res)
        tf_gene<-unique(lisa_res_in_metroscreen$factor)

        M_TR=cal_tf_score(tf_gene,data) ## tf activity
    #         cal_target_score_res=lapply(as.list(colnames(M_TR)),
    #                                                 cal_target_score,
    #                                                 lisa_res=lisa_res,
    #                                                 data=data) ## tf target score
        
        cl <- makeCluster(20)
        clusterExport(cl, c("cal_target_score", "lisa_res_in_metroscreen", "data",'RP_score'),envir = environment())
        cal_target_score_res <- parLapply(cl, as.list(colnames(M_TR)), cal_target_score,lisa_res_in_metroscreen=lisa_res_in_metroscreen, data=as.data.frame(data),RP_score=RP_score)
        M_Target=Reduce(function(x,y){cbind(x,y)},cal_target_score_res)
        stopCluster(cl)
        
        lisa_res_in_metroscreen_use=lisa_res_in_metroscreen[!duplicated(lisa_res_in_metroscreen$factor),c('factor','summary_p_value')]
        lisa_res_in_metroscreen_use$score=(-log(lisa_res_in_metroscreen_use[,c('summary_p_value')]))
        rownames(lisa_res_in_metroscreen_use)=lisa_res_in_metroscreen_use[,'factor']
        lisa_res_in_metroscreen_use=lisa_res_in_metroscreen_use[colnames(M_TR),] ## 

        M_actifvity=(M_TR+M_Target) ## sum the tf target score and tf activity score
        M_actifvity_scale=as.data.frame(lapply(M_actifvity, function(x) (x - min(x)) / (max(x) - min(x)))) ##scale
        rownames(M_actifvity_scale)=rownames(M_actifvity)
        S=M_actifvity_scale*lisa_res_in_metroscreen_use$score  ## scripro
        saveRDS(S,tr_activity_res)
    }

    ## get the tr-target interaction
    tr_cor_res=paste0(output_path,file_name,"/",file_name,":tr_activity_cor.rds")
    if (!file.exists(tr_cor_res)){
        data=readRDS(file_res)
        data_tf=readRDS(tr_activity_res)
        data_t=t(as.matrix(data[,rownames(data_tf)]))
        result=data_t %>%
            apply(2,function(x){10*(2**x - 1)}) %>%
            apply(2,function(x){log2(mean(x) + 1)}) 
        expressed_genes <- names(result[order(result, decreasing = TRUE)][1:3000]) 
        data_tf=data_tf[,colnames(data_tf) %in% expressed_genes]
        data_mg=t(as.matrix(data[rownames(data) %in% unique(unlist(MM_list)),rownames(data_tf)]))
        cor=cor(data_tf,data_mg)
        get_cor_res=lapply(as.list(interested_MM),get_cor,cor=cor,MM_list=MM_list,class='tr')

        names(get_cor_res)=interested_MM
        saveRDS(get_cor_res,tr_cor_res)
    }
    cat("Finished tr_target_interaction\n")
}

#' get_causal
#' 
#' @param data
#' @param tr_activity
#' @param lr_activity 
#' @param factor
#' @param MM_list
#' @param pathway_index
#' @param times
#' @param factor_class
#' 
#' @export

get_causal <- function(data, tr_activity, lr_activity, factor, MM_list, pathway_index, times, factor_class) {
    # use switch 
    df_activity = switch(factor_class,
                         'TF' = as.data.frame(tr_activity[, factor, drop = FALSE]),
                         'GG' = as.data.frame(data[, factor, drop = FALSE]),
                         'Ligand' = as.data.frame(lr_activity[, factor, drop = FALSE]))
    rownames(df_activity) = switch(factor_class,
                                   'TF' = rownames(tr_activity),
                                   'GG' = rownames(data),
                                   'Ligand' = rownames(lr_activity))
    colnames(df_activity) = factor

    df_gene = as.data.frame(data[, colnames(data) %in% MM_list[[pathway_index]], drop = FALSE])
    unique_counts = sapply(df_gene, function(col) length(unique(col)))
    df_gene_filtered = df_gene[, unique_counts >= 10]

    # cal cor
    cor = cor(df_activity, df_gene_filtered)
    cor=cor[,order(cor[1,],decreasing = TRUE)]
    gene_max=names(cor[1])
                           
    df_use = data[, colnames(data) %in% c(factor, unique(MM_list[[pathway_index]]))]
    results = vector("list", length = times)

    for (i in 1:times) {
        rows = if (times == 1) 1:nrow(df_use) else sample(1:nrow(df_use), 100)
        df_final=df_use[rows,]
        unique_counts <- sapply(as.data.frame(df_final), function(col) length(unique(col)))
        df_final_filtered <- df_final[, unique_counts >= ifelse(nrow(df_final)*0.1<10,10,nrow(df_final)*0.1)]

        # causal analysis
        kpc <- kpc(suffStat = list(data = df_final_filtered, ic.method = "dcc.gamma"),
                       indepTest = kernelCItest,
                       alpha = 0.1,
                       labels = colnames(df_final_filtered),
                       u2pd = "relaxed",
                       skel.method = "stable.fast",
                       verbose = TRUE)

        g <- graph_from_graphnel(kpc@graph, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
        edges_df <- as.data.frame(get.edgelist(g), stringsAsFactors=FALSE)
        edges_df_use=edges_df[edges_df$V1==factor | edges_df$V2==factor,]
        colnames(edges_df_use)=c('from','to')

        results[[i]] = if (nrow(edges_df_use) > 0 && any(edges_df_use$from == factor & edges_df_use$to == gene_max)) 1 else 0
    }

    result = mean(unlist(results))
    df = data.frame(factor = factor, direction = ifelse(result >= 0.2, 'regulator', 'effector'), stringsAsFactors = FALSE)
    return(df)
}

#' repref_each
#' 
#' @param tr_cor
#' @param tr_activity 
#' @param lr_activity 
#' @param ligand_cor
#' @param gene_cor
#' @param pathway_index
#' @param MM_list
#' @param data
#' @param output_path
#' @param file_name
#' 
#' @export

repref_each=function(tr_cor,tr_activity,lr_activity,ligand_cor,gene_cor,pathway_index,MM_list,data,output_path,file_name){
    tmp_dir <- file.path(output_path,file_name,'MetaRegulon')
    if (!dir.exists(tmp_dir)) {
        dir.create(tmp_dir)
    }
    recom_res=paste0(output_path,file_name,'/MetaRegulon/',file_name,":",pathway_index,'.txt')
    if (!file.exists(recom_res) && pathway_index %in% c(names(tr_cor), names(ligand_cor), names(gene_cor))){
            cor_lists = list(tr_cor[[pathway_index]], ligand_cor[[pathway_index]], gene_cor[[pathway_index]])

            if (all(sapply(cor_lists, nrow) > 0)){
                li_merge = Reduce(function(x, y) merge(x, y, by = 'MR', all = TRUE), cor_lists)
                rownames(li_merge) = li_merge$MR

                score_columns = c("TR_Target_interaction", "Ligand_Receptor_interaction", "Gene_Gene_interaction")
                rank_lists = lapply(score_columns, function(sc) {
                    valid_rows = li_merge[!is.na(li_merge[[sc]] > 0), ]
                    rownames(valid_rows[order(valid_rows[[sc]], decreasing = TRUE), ])
                })
                names(rank_lists) = score_columns
                glist = rank_lists
                                  
                ag=aggregateRanks(glist,exact=TRUE)
                ag=ag[rownames(li_merge),]
                li_merge$ag_score=ag$Score
                li_merge[is.na(li_merge)]<-0
                li_merge<-li_merge[,-1]
                                  
                res_pr=psel(li_merge, high(TR_Target_interaction) *high(Ligand_Receptor_interaction)*high(Gene_Gene_interaction), top_level=30)
                res_pr$Final_score=res_pr$ag_score*res_pr$.level
                res_pr$gene=rownames(res_pr)
                res_pr <- res_pr %>%
                    group_by(.level) %>%
                    arrange(ag_score,.by_group = TRUE) %>%
                    as.data.frame()
                res_pr$rank=1:nrow(res_pr)
                rownames(res_pr)=res_pr$gene

                ## the source of the regulation
                ## identify the regulation source
                res_pr$resource <- ifelse(res_pr$TR_Target_interaction > 0 | res_pr$Gene_Gene_interaction > 0,'intrinsic', 'extrinsic')

                if (length(MM_list[[pathway_index]])<30 & nrow(data)>=10){
                    res_pr$direction='' ## setting the direction of these factor and downstream metabolic reaction
                    times<<-ifelse(nrow(data)>=100,5,1)
                    factor_tf=res_pr[res_pr$TR_Target_interaction>0,'gene']
                    factor_gg=res_pr[res_pr$TR_Target_interaction==0 & res_pr$Ligand_Receptor_interaction==0,'gene']
                    factor_lr_paracrine=res_pr[res_pr$TR_Target_interaction==0 & res_pr$Gene_Gene_interaction==0,'gene']
                    factor_lr_autocrine=res_pr[res_pr$TR_Target_interaction==0 & !res_pr$Gene_Gene_interaction==0 & !res_pr$Ligand_Receptor_interaction==0,'gene']

                    ## for ligand, if the ligand is paracrine, then it is regulator, if the ligand is autocrine, MetroSCREEN will identify it's casauation
                    res_pr$direction=ifelse(rownames(res_pr) %in% factor_lr_paracrine,'regulator',res_pr$direction)

                    metabolic_gene=unique(MM_list[[pathway_index]])

                    if (length(factor_tf)>0){
                        cl <- makeCluster(30)
                        clusterExport(cl, c("get_causal", "data", "tr_activity", "lr_activity", "MM_list", "pathway_index", "times"),envir = environment())
                        clusterEvalQ(cl, {
                            suppressMessages({
                                libraries <- c( "Seurat", "dplyr", "tidyr", "ggplot2",
                                            "nichenetr", "ggpubr", "GSVA", 
                                            "ComplexHeatmap", "limma", "RColorBrewer",
                                            "rPref", "igraph", "RobustRankAggreg", "parallel",
                                            "MetroSCREEN", "reshape2")
                                lapply(libraries, library, character.only = TRUE)
                                # source('/fs/home/tangke/metabolism/tool/MetroSCREEN/R/kpc.r')
                            })
                        })
                        factor_tf_res <- do.call(rbind, parLapply(cl, factor_tf, function(factor){
                            get_causal(data=data, tr_activity=tr_activity, lr_activity=lr_activity, MM_list=MM_list, pathway_index=pathway_index, times=times, factor_class='TF',factor=factor)
                        }))
                        rownames(factor_tf_res) <- factor_tf_res$factor
                        stopCluster(cl)
                        # factor_tf_res=do.call(rbind,
                        #     lapply(as.list(factor_tf),get_causal,data=data,tr_activity=tr_activity,lr_activity=lr_activity,MM_list=MM_list,pathway_index=pathway_index,times=times,factor_class='TF'))
                        # rownames(factor_tf_res)=factor_tf_res$factor
                        res_pr$direction=ifelse(rownames(res_pr) %in% factor_tf,factor_tf_res$direction,res_pr$direction)
                    }

                    if (length(factor_gg)>0){
                        cl <- makeCluster(30)
                        clusterExport(cl, c("get_causal", "data", "tr_activity", "lr_activity", "MM_list", "pathway_index", "times"),envir = environment())
                        clusterEvalQ(cl, {
                            suppressMessages({
                                libraries <- c( "Seurat", "dplyr", "tidyr", "ggplot2",
                                            "nichenetr", "ggpubr", "GSVA", 
                                            "ComplexHeatmap", "limma", "RColorBrewer",
                                            "rPref", "igraph", "RobustRankAggreg", "parallel",
                                            "MetroSCREEN", "reshape2")
                                lapply(libraries, library, character.only = TRUE)
                                # source('/fs/home/tangke/metabolism/tool/MetroSCREEN/R/kpc.r')
                            })
                        })
                        factor_gg_res <- do.call(rbind, parLapply(cl, factor_gg, function(factor){
                            get_causal(data=data, tr_activity=tr_activity, lr_activity=lr_activity, MM_list=MM_list, pathway_index=pathway_index, times=times, factor_class='GG',factor=factor)
                        }))
                        rownames(factor_gg_res) <- factor_gg_res$factor
                        stopCluster(cl)
                        # factor_gg_res=do.call(rbind,
                        #     lapply(as.list(factor_gg),get_causal,data=data,tr_activity=tr_activity,lr_activity=lr_activity,MM_list=MM_list,pathway_index=pathway_index,times=times,factor_class='GG'))
                        # rownames(factor_tf_res)=factor_tf_res$factor
                        res_pr$direction=ifelse(rownames(res_pr) %in% factor_gg,factor_gg_res$direction,res_pr$direction)
                    }

                    if (length(factor_lr_autocrine)>0){
                        cl <- makeCluster(30)
                        clusterExport(cl, c("get_causal", "data", "tr_activity", "lr_activity", "MM_list", "pathway_index", "times"),envir = environment())
                        clusterEvalQ(cl, {
                            suppressMessages({
                                libraries <- c("Seurat", "dplyr", "tidyr", "ggplot2",
                                            "nichenetr", "ggpubr", "GSVA",
                                            "ComplexHeatmap", "limma", "RColorBrewer",
                                            "rPref", "igraph", "RobustRankAggreg", "parallel",
                                            "MetroSCREEN", "reshape2")
                                lapply(libraries, library, character.only = TRUE)
                                # source('/fs/home/tangke/metabolism/tool/MetroSCREEN/R/kpc.r')
                            })
                        })
                        factor_lr_autocrine_res <- do.call(rbind, parLapply(cl, factor_lr_autocrine, function(factor){
                            get_causal(data=data, tr_activity=tr_activity, lr_activity=lr_activity, MM_list=MM_list, pathway_index=pathway_index, times=times, factor_class='Ligand',factor=factor)
                        }))

                        rownames(factor_lr_autocrine_res) <- factor_lr_autocrine_res$factor
                        stopCluster(cl)
                        # factor_tf_res=do.call(rbind,
                        #     lapply(as.list(factor_tf),get_causal,data=data,tr_activity=tr_activity,lr_activity=lr_activity,MM_list=MM_list,pathway_index=pathway_index,times=times,factor_class='TF'))
                        # rownames(factor_tf_res)=factor_tf_res$factor
                        res_pr$direction=ifelse(rownames(res_pr) %in% factor_lr_autocrine,factor_lr_autocrine_res$direction,res_pr$direction)
                    }
                }
                write.csv(res_pr,recom_res)
        }
    }
}

