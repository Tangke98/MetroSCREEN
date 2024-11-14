#' Calculate the MetaModule score for each single cell/Metacell/Sample.
#' 
#' @param object A matrix, the rows are genes and the columns are the samples/Metacells/cells, the value is the expression
#' @param module_list Metabolic reaction list or interested gene sets
#' @param output_path Path to save the MetaModule score result
#' @param file_name File name of the MetaModule score result
#' @param parallel The parallel size to process data
#' 
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import ggpubr
#' @import GSVA
#' @import ComplexHeatmap
#' @import RColorBrewer
#' 
#' @export
#' 

cal_MetaModule=function(object, module_list, output_path, file_name,parallel=20){
    
    if (ncol(object) > 10000){
        tmp_dir <- file.path(output_path, "MM_tmp")
        if (!dir.exists(tmp_dir)) {
            dir.create(tmp_dir)
        }
        if (length(module_list) > 20){
            group_indices <- cut(seq_along(module_list), breaks=seq(1, length(module_list), by=20), labels=FALSE)
            # split the module
            sublists <- split(module_list, group_indices)
            for (i in 1:length(sublists)){
                MM_use = sublists[[i]]
                output_name_use = paste0(file_name, "_", names(sublists[i]))
                gsva <- gsva(object, MM_use, method = "gsva", 
                            min.sz = 1, max.sz = 500, kcdf = "Gaussian", mx.diff = TRUE, 
                            verbose = TRUE, parallel.sz = parallel)
                saveRDS(gsva, paste0(output_path, 'MM_tmp/', output_name_use, ".rds"))
            }
            files = list.files(paste0(output_path, 'MM_tmp/'))
            merge = function(x){
                data <- readRDS(paste0(paste0(output_path, 'MM_tmp/'), x))
                return(data)
            }
            merge_res = do.call(rbind, lapply(as.list(files), merge))
            saveRDS(merge_res, paste0(output_path, file_name, '.rds'))
            unlink(tmp_dir, recursive = TRUE)
        } else {
            gsva <- gsva(object, module_list, method = "gsva", 
                        min.sz = 1, max.sz = 500, kcdf = "Gaussian", mx.diff = TRUE, 
                        verbose = TRUE, parallel.sz = parallel)
            saveRDS(gsva, paste0(output_path, file_name, ".rds"))
        }

    } else if (ncol(object) > 2000) {
        tmp_dir <- file.path(output_path, "MM_tmp")
        if (!dir.exists(tmp_dir)) {
            dir.create(tmp_dir)
        }
        if (length(module_list) > 100){
            group_indices <- cut(seq_along(module_list), breaks=seq(1, length(module_list), by=100), labels=FALSE)
            # split the module
            sublists <- split(module_list, group_indices)
            for (i in 1:length(sublists)){
                MM_use = sublists[[i]]
                output_name_use = paste0(file_name, "_", names(sublists[i]))
                gsva <- gsva(object, MM_use, method = "gsva", 
                            min.sz = 1, max.sz = 500, kcdf = "Gaussian", mx.diff = TRUE, 
                            verbose = TRUE, parallel.sz = parallel)
                saveRDS(gsva, paste0(output_path, 'MM_tmp/', output_name_use, ".rds"))
            }
            files = list.files(paste0(output_path, 'MM_tmp/'))
            merge = function(x){
                data <- readRDS(paste0(paste0(output_path, 'MM_tmp/'), x))
                return(data)
            }
            merge_res = do.call(rbind, lapply(as.list(files), merge))
            saveRDS(merge_res, paste0(output_path, file_name, '.rds'))
            unlink(tmp_dir, recursive = TRUE)
        } else {
            gsva <- gsva(object, module_list, method = "gsva", 
                        min.sz = 1, max.sz = 500, kcdf = "Gaussian", mx.diff = TRUE, 
                        verbose = TRUE, parallel.sz = parallel)
            saveRDS(gsva, paste0(output_path, file_name, ".rds"))
        }
        
    } else {
        gsva <- gsva(object, module_list, method = "gsva", 
                     min.sz = 1, max.sz = 500, kcdf = "Gaussian", mx.diff = TRUE, 
                     verbose = TRUE, parallel.sz = parallel)
        saveRDS(gsva, paste0(output_path, file_name, ".rds"))
    }
}

