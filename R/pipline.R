#' Validate input data for calculate_all_lrlr function.
#'
#' @param Input A matrix or data frame containing gene expression data(Recommend normalized data).
#' @param Coordinate Coordinate matrix.
#' @param st_type A character string indicating the type of spatial transcriptomics data, either "low" or "high". Default is "low".
#' @param deconvolution_result A matrix or data frame containing deconvolution results. Required if st_type is "low".
#' @param scRNA_mean_expression A matrix or data frame containing single-cell RNA-seq mean expression data. Required if st_type is "low".
#' @param lr_db A data frame containing ligand-receptor pairs.
#' @return TRUE if all inputs are valid, otherwise stop with an error message.
#' @export
validate_inputs <- function(Input, Coordinate, st_type, deconvolution_result = NULL, scRNA_mean_expression = NULL, lr_db = NULL) {
  if (!is.matrix(Input)) {
    stop("Input must be a matrix.")
  }
  if (!is.matrix(Coordinate)) {
    stop("Coordinate must be a matrix.")
  }
  if (!identical(colnames(Input), rownames(Coordinate))) {
    stop("The col.names of Input matrix must be identical to the row.names of Coordinate matrix.")
  }
  if (!(st_type %in% c('low', 'high'))) {
    stop("st_type must be either 'low' or 'high'.")
  }
  if (st_type == 'low') {
    if (is.null(deconvolution_result) || !is.matrix(deconvolution_result)) {
      stop("For low resolution data, deconvolution_result must be a matrix.")
    }
    if (is.null(scRNA_mean_expression) || !is.matrix(scRNA_mean_expression)) {
      stop("For low resolution data, scRNA_mean_expression must be a matrix.")
    }
    if (!identical(rownames(deconvolution_result), colnames(Input))) {
      stop("The col.names of Input matrix must be identical to the row.names of deconvolution matrix.")
    }
    if (!identical(colnames(scRNA_mean_expression), colnames(deconvolution_result))) {
      stop("The col.names of deconvolution matrix must be identical to the col.names of scRNA mean expression matrix.")
    }
  }
  if(!is.null(lr_db)){
    if(!all(c('ligand', 'receptor') %in% colnames(lr_db))){
      stop("The lrpair database must contain 'ligand' and 'receptor' columns.")
    }
  }
  return(TRUE)
}



#' Selected cell types in each spot.
#'
#' @param RCTD_result Result from the function 'load_ccRCC_RCTD_result'
#' @param max_n_types The most cell types will be remained in a spot
#' @param cell_ratio_cutoff Cell ratio must higer than a cutoff. Default = 0.1.
#' 
#' @return A matrix containing the cell type name in each spot.
#' @export 
#' @examples
#' spots_celltype_matrix = Construct_spots_celltype(coords)
Construct_spots_celltype = function(RCTD_result, max_n_types = 3, cell_ratio_cutoff = 0.1){
  select_celltype_from_RCTD = function(spot_RCTD){
    spot_remain = spot_RCTD[order(spot_RCTD,decreasing = T)][1:max_n_types]
    spot_remain = spot_remain[spot_remain > cell_ratio_cutoff]
    spot_remain = names(spot_remain)
    if(length(spot_remain) < max_n_types){
      spot_remain = c(spot_remain, rep('none',max_n_types-length(spot_remain)))
    }
    return(spot_remain)
  }
  spots_celltype_matrix = apply(t(RCTD_result),2,select_celltype_from_RCTD) %>% t() %>% as.data.frame()
  colnames(spots_celltype_matrix) = paste0('c',c(1:max_n_types))
  spots_celltype_matrix = as.matrix(spots_celltype_matrix)
  return(spots_celltype_matrix)
}


#' Reconstructed RCTD result based on 'Construct_spots_celltype'.
#'
#' @param RCTD_result Result from the function 'load_ccRCC_RCTD_result'.
#' @param spots_celltype_matrix Result from the function 'Construct_spots_celltype'.
#' @param cell_ratio_cutoff Cell ratio must higer than a cutoff. Default = 0.1.
#' 
#' @return A matrix containing the cell type name in each spot.
#' @export 
#' @examples
#' spots_celltype_matrix = Reconstruct_RCTD(coords)
Reconstruct_RCTD = function(RCTD_result, spots_celltype_matrix){
  all_spots = rownames(RCTD_result)
  all_celltype = colnames(RCTD_result)
  deconvolution_reconsitution = matrix(0, nrow = length(all_spots), ncol = length(all_celltype), dimnames = list(all_spots, all_celltype))
  for (i in 1:length(all_spots)) {
    ct = spots_celltype_matrix[i,,drop = T]
    if ('none' %in% ct) {
      ct = ct[-which(ct == 'none')]
    }
    re = RCTD_result[i,ct]
    re = round(re/sum(re) * (1/min(re)))
    deconvolution_reconsitution[i,ct] = re
  }
  # normarlize
  deconvolution_reconsitution = apply(deconvolution_reconsitution,1,function(x){x/sum(x)}) %>% t() %>% as.matrix()
  return(deconvolution_reconsitution)
}


#' Constucted a expression list of every cell types in each spot.
#'
#' @param st_exp_matrix Expression matrix of ST dataset.
#' @param sc_AverageExp Mean expression from scRNA-seq dataset.
#' @param RCTD_reconsitution Result from the function 'Reconstruct_RCTD' or a matrix like it.
#' 
#' @return A matrix containing the cell type name in each spot.
#' @export 
#' @examples
#' all_celltype_expression_list = Get_reconstruct_ST_expression_ct_list(st_exp_matrix,sc_AverageExp,RCTD_reconsitution)

Get_reconstruct_ST_expression_ct_list = function(st_exp_matrix, sc_AverageExp, deconvolution_reconsitution){
  re = array(0, dim = c(nrow(st_exp_matrix), ncol(st_exp_matrix), ncol(deconvolution_reconsitution)), 
               dimnames = list(rownames(st_exp_matrix), colnames(st_exp_matrix), colnames(deconvolution_reconsitution)))
  for (i in 1:nrow(st_exp_matrix)) {
    for (j in 1:ncol(st_exp_matrix)) {
      exp_all = st_exp_matrix[i,j]
      if(exp_all == 0){next}
      t = sum(deconvolution_reconsitution[j,] * sc_AverageExp[i,])
      if(t > 0){
        a = exp_all/t
      }else{next}
      exp_adj = deconvolution_reconsitution[j,] * sc_AverageExp[i,] * a
      re[i,j,] = exp_adj
    }
  }
  express_recon = list()
  for (i in 1:dim(re)[3]) {
    express_recon[[i]] = as(re[,,i], 'dgCMatrix')
  }
  names(express_recon) = colnames(deconvolution_reconsitution)
  return(express_recon)
}


#' Calculate combined intracelluar Receptor-Liagnd siganal weight
#'
#' @param Input A matrix containing gene expression data. (Recommend scRNA_mean_expression)
#' @param cell_ids The col.names selected to be calculate. If NULL, it will be selected all columns.
#' @param PPI_db A database of protein-protein interaction. 
#' If NULL, it will be inferred from the data which contains a column called 'functional' which means stimulation, inhibition or unknown function by 1,-1,0. 
#' If you want to provide yourself database, it should contain a 'functional' column and set the parameter 1.
#' @param Pro_loc_db A database of protein location. If NULL, it will be inferred from the data.
#' Default is organized Uniprot database which contains a column called 'location' which contains param 'Membrane', 'Intracellular' and 'Nuluar'.
#' If you want to provide yourself database, it should contain a 'location' column and set param.
#' @param lr_db A database of ligand-receptor. If NULL, it will be inferred from the data.
#' Default is the lr_db database from 'Spatalk' which contains two columns called 'ligand' and 'receptor'
#' @param f_edge A function to calculate score of each edge, default is Hill function.
#' @param Kn The parameter of Hill function. Default is 2.
#' @param cutoff1 A quantile to select the probility of R-L pairs over cutoff. Default is 0.5.
#' @param cutoff2 A quantile to build a subgraph containing genes over cutoff. Default is 0.9.
#' @param max_steps The max steps to be find in each R-L signal networks. Default is 5.
#' @param max_num The number of paths to be selected to calculate path scores. Default is 20.
#' @param circuit_scale_factor Adjustment of scaling factor for circuit score. Default is 0.1.
#' @param adjust_alpha A logical parameter that determines whether to adjust alpha. Default is TRUE.
#' @param target_ratio The goal of alpha adjustment. Default is 0.5.
#' @param parallel A logical parameter that determines whether to do parallel. Default is TRUE.
#' @param n_core The number of cores used for parallel computing. If NULL, it will be detectCores() - 2.
#' 
#' @return A dataframe of Receptor-Ligand signal weights of each cell_id
#' @export 
#' @examples
#' utils::data("test_scRNA_mean_expression", package = "stMCCC")
#' results = Calculate_Combined_RL_weight(test_scRNA_mean_expression)


Calculate_Combined_RL_weight = function(Input,
                                        cell_ids = NULL,
                                        PPI_db = NULL, 
                                        Pro_loc_db = NULL, 
                                        lr_db = NULL,
                                        f_edge = NULL, 
                                        Kn = 2,
                                        cutoff1 = 0.5,
                                        cutoff2 = 0.9, 
                                        max_steps = 5, 
                                        max_num = 20, 
                                        circuit_scale_factor = 0.1,
                                        adjust_alpha = TRUE,
                                        target_ratio = 0.5,
                                        parallel = TRUE,
                                        n_core = NULL) {
  
  calculate_path_scores <- function(paths, all_edges, max_num = NULL) {
    if(length(paths) == 0) {return(NULL)}
    if (!is.null(max_num) && max_num > 0 && max_num < length(paths)) {paths = paths[sample(length(paths), max_num)]}
    scores <- sapply(paths, function(path) {
      score <- 0
      for (i in 1:(length(path)-1)) {
        source_gene <- path[i]
        target_gene <- path[i+1]
        weight = all_edges$weight[all_edges$from == source_gene & all_edges$to == target_gene]
        score = score + weight
      }
      return(score/(length(path)-1))
    })
    paths_str <- sapply(paths, paste, collapse = " -> ")
    lengths <- sapply(paths, function(x) length(x)-1)
    return(data.table('Path' = paths_str, 'Score' = scores,'Length' = lengths,stringsAsFactors = FALSE))
  }
  
  
  build_directed_conductance_matrix <- function(graph, expr_vec, f_edge) {
    edge_df <- as_data_frame(graph)
    genes <- unique(c(edge_df$from, edge_df$to))
    A <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
    for (i in seq_len(nrow(edge_df))) {
      from <- edge_df$from[i]
      to <- edge_df$to[i]
      if (from %in% names(expr_vec) && to %in% names(expr_vec)) {
        A[from, to] <- f_edge(expr_vec[from], expr_vec[to])
      }
    }
    return(A)
  }
  # 1. 参数验证与默认数据库加载
  if (is.null(PPI_db)) { utils::data("omnipath_ppi", package = "stMCCC"); PPI_db <- omnipath_ppi[, c("source_genesymbol", "target_genesymbol", "functional")]}
  stopifnot(c("source_genesymbol", "target_genesymbol", "functional") %in% colnames(PPI_db))
  if (is.null(Pro_loc_db)) {utils::data("Uniprot_location", package = "stMCCC"); Pro_loc_db <- Uniprot_location[, c("genesymbol", "location")]}
  stopifnot(c("genesymbol", "location") %in% colnames(Pro_loc_db))
  stopifnot(c("Intracellular", "Membrane", "Nucleus") %in% names(table(Pro_loc_db$location)))
  if (is.null(lr_db)) {utils::data("SpaTalk_lrpairs", package = "stMCCC"); lr_db <- SpaTalk_lrpairs[, c("ligand", "receptor")]}
  stopifnot(c("ligand", "receptor") %in% colnames(lr_db))
  if (is.null(f_edge)) { f_edge <- function(x, y) {return(x*y / (Kn + x * y)) } }
  
  # 2. 构建原始PPI网络图
  ppi_graph <- igraph::graph_from_data_frame(PPI_db[,c(1,2)],directed = TRUE)
  igraph::E(ppi_graph)$interaction <- PPI_db$functional
  
  exp_data = Input
  if (is.null(cell_ids)){cell_ids <- colnames(exp_data)}
  stopifnot(!is.null(cell_ids), length(cell_ids) > 0)
  
  # 4.选择合适的lrpairs和基因
  ligands = unique(lr_db$ligand)
  receptors = unique(lr_db$receptor)
  expressed_genes = rownames(exp_data)[rowSums(exp_data) > 0 & rowSums(exp_data) > quantile(rowSums(exp_data), cutoff1)]
  expressed_genes = expressed_genes[expressed_genes %in% Pro_loc_db$genesymbol]
  
  ligands = ligands[ligands %in% expressed_genes]
  receptors = receptors[receptors %in% expressed_genes]
  rl_pairs = expand.grid('receptor' = receptors, 'ligand' = ligands, stringsAsFactors = F)
  
  if (parallel) {
    if(is.null(n_core)){
      n_core <- detectCores() - 2
    }
    cl <- makeCluster(n_core, type = "FORK")
    results <- do.call(rbind, pbapply::pblapply(cell_ids, cl = cl, function(cid) {
      expr_vec <- exp_data[, cid]
      names(expr_vec) <- rownames(exp_data)
      filter_genes <- names(expr_vec)[expr_vec > 0]
      if(length(filter_genes) > 3000){
        filter_genes = names(expr_vec)[expr_vec > quantile(expr_vec, cutoff2)]
      }
      subgraph <- igraph::induced_subgraph(ppi_graph, which(V(ppi_graph)$name %in% filter_genes))
      # 将subgraph的所有边的权重算一下
      all_edges = as_data_frame(subgraph, what = 'edges')
      all_edges$exp_from = expr_vec[all_edges$from]
      all_edges$exp_to = expr_vec[all_edges$to]
      
      all_edges$weight = 0
      all_edges$weight[all_edges$interaction == 1] = 1+all_edges$exp_from[all_edges$interaction == 1] * all_edges$exp_to[all_edges$interaction == 1]/(Kn + all_edges$exp_from[all_edges$interaction == 1] * all_edges$exp_to[all_edges$interaction == 1])
      all_edges$weight[all_edges$interaction == -1] = Kn/(Kn + all_edges$exp_from[all_edges$interaction == -1] * all_edges$exp_to[all_edges$interaction == -1])
      all_edges$weight[all_edges$interaction == 0] = 1
      
      # 电路分数优化计算
      G <- build_directed_conductance_matrix(subgraph, expr_vec, f_edge)
      nodes <- rownames(G)
      D <- diag(rowSums(G))
      L <- D - G # 构建拉普拉斯矩阵
      mu = 1e-6
      A_reg_full <- t(L) %*% L + mu * diag(length(nodes))
      
      rl_pairs_ct = rl_pairs[(rl_pairs$receptor %in% rownames(L)) & (rl_pairs$ligand %in% rownames(L)),]
      mem = Pro_loc_db$genesymbol[Pro_loc_db$location == 'Membrane']
      rl_pairs_ct = rl_pairs_ct[rl_pairs_ct$receptor %in% mem,]
      
      del_idx = sample(which(!(rownames(L) %in% rl_pairs_ct$receptor) & !(rownames(L) %in% rl_pairs_ct$ligand)), 1)
      if(is.null(del_idx)){del_idx = 1}
      A_reg_reduced = A_reg_full[-del_idx, -del_idx]
      A_reg_reduced_v = solve(A_reg_reduced)
      
      L_trans <- t(L)
      I_full <- rep(0, length(nodes))
      names(I_full) <- nodes
      
      # 对每个细胞计算所有RL对的得分
      inner_results <- pbapply::pblapply(1:nrow(rl_pairs_ct), function(g) {
        receptor <- rl_pairs_ct$receptor[g]
        ligand <- rl_pairs_ct$ligand[g]
        all_paths <- all_simple_paths(subgraph,from = which(V(subgraph)$name == receptor),
                                      to = which(V(subgraph)$name == ligand),
                                      mode = "out", cutoff = max_steps)
        if (length(all_paths) == 0) {
          return(NULL)
        } else {
          # 验证路径并计算分数
          valid_paths <- list()
          for (path in all_paths) {
            gene_path <- V(subgraph)$name[path]
            valid <- TRUE
            # Receptor should on the cell membrane
            receptor_loc <- Pro_loc_db$location[Pro_loc_db$genesymbol == gene_path[1]]
            if (!"Membrane" %in% receptor_loc) valid <- FALSE
            # Ohter molecules should be intracellular molecules
            for (i in 2:(length(gene_path)-2)) {
              gene <- gene_path[i]
              loc <- Pro_loc_db$location[Pro_loc_db$genesymbol == gene]
              if (!"Intracellular" %in% loc) valid <- FALSE
            }
            # The penultimate molecule should be in Nucleus
            tf_loc <- Pro_loc_db$location[Pro_loc_db$genesymbol == gene_path[length(gene_path)-1]]
            if (!"Nucleus" %in% tf_loc) valid <- FALSE
            if (valid) valid_paths <- c(valid_paths, list(gene_path))
          }
          if(length(valid_paths) == 0){
            return(NULL)
          }else{
            path_scores <- calculate_path_scores(paths = valid_paths, all_edges, max_num = max_num)
            path_score = max(path_scores$Score)
          }
          
          # 计算电路分数
          I_full[c(receptor, ligand)] <- c(1, -1)
          b_reg <- L_trans[-del_idx, ] %*% I_full
          V_reduced <- A_reg_reduced_v %*% b_reg
          V_full <- rep(0, length(nodes))
          V_full[-del_idx] <- V_reduced
          names(V_full) <- nodes
          V_full_norm = (V_full-mean(V_full))/sd(V_full)
          circuit_score = (V_full_norm[receptor] - V_full_norm[ligand]) * circuit_scale_factor # 缩放电路分数
          
          # 动态调整alpha
          if (adjust_alpha && path_score > 0 && circuit_score != 0) {
            # 动态调整alpha
            current_ratio <- path_score / abs(circuit_score)
            alpha <- target_ratio / (target_ratio + current_ratio)
          }
          combined_score <- alpha * path_score + (1 - alpha) * circuit_score
          return(data.table('receptor' = receptor,
                            'ligand' = ligand,
                            'cell_id' = cid,
                            'path_score' = path_score,
                            'circuit_score' = circuit_score,
                            'alpha' = alpha,
                            'combined_score' = combined_score))
        }
      })
      return(rbindlist(inner_results, fill = T))
    }
    ))
  }else{
    results <- do.call(rbind, pbapply::pblapply(cell_ids, function(cid) {
      expr_vec <- exp_data[, cid]
      names(expr_vec) <- rownames(exp_data)
      filter_genes <- names(expr_vec)[expr_vec > 0]
      if(length(filter_genes) > 3000){
        filter_genes = names(expr_vec)[expr_vec > quantile(expr_vec, cutoff2)]
      }
      
      subgraph <- igraph::induced_subgraph(ppi_graph, which(V(ppi_graph)$name %in% filter_genes))
      # 将subgraph的所有边的权重算一下
      all_edges = as_data_frame(subgraph, what = 'edges')
      all_edges$exp_from = expr_vec[all_edges$from]
      all_edges$exp_to = expr_vec[all_edges$to]
      
      all_edges$weight = 0
      all_edges$weight[all_edges$interaction == 1] = 1+all_edges$exp_from[all_edges$interaction == 1] * all_edges$exp_to[all_edges$interaction == 1]/(Kn + all_edges$exp_from[all_edges$interaction == 1] * all_edges$exp_to[all_edges$interaction == 1])
      all_edges$weight[all_edges$interaction == -1] = Kn/(Kn + all_edges$exp_from[all_edges$interaction == -1] * all_edges$exp_to[all_edges$interaction == -1])
      all_edges$weight[all_edges$interaction == 0] = 1
      
      # 电路分数优化计算
      G <- build_directed_conductance_matrix(subgraph, expr_vec, f_edge)
      nodes <- rownames(G)
      D <- diag(rowSums(G))
      L <- D - G # 构建拉普拉斯矩阵
      mu = 1e-6
      A_reg_full <- t(L) %*% L + mu * diag(length(nodes))
      
      rl_pairs_ct = rl_pairs[(rl_pairs$receptor %in% rownames(L)) & (rl_pairs$ligand %in% rownames(L)),]
      mem = Pro_loc_db$genesymbol[Pro_loc_db$location == 'Membrane']
      rl_pairs_ct = rl_pairs_ct[rl_pairs_ct$receptor %in% mem,]
      
      del_idx = sample(which(!(rownames(L) %in% rl_pairs_ct$receptor) & !(rownames(L) %in% rl_pairs_ct$ligand)), 1)
      if(is.null(del_idx)){del_idx = 1}
      A_reg_reduced = A_reg_full[-del_idx, -del_idx]
      A_reg_reduced_v = solve(A_reg_reduced)
      
      L_trans <- t(L)
      I_full <- rep(0, length(nodes))
      names(I_full) <- nodes
      
      # 对每个细胞计算所有RL对的得分
      
      inner_results <- pbapply::pblapply(1:nrow(rl_pairs_ct), function(g) {
        receptor <- rl_pairs_ct$receptor[g]
        ligand <- rl_pairs_ct$ligand[g]
        all_paths <- all_simple_paths(subgraph,from = which(V(subgraph)$name == receptor),
                                      to = which(V(subgraph)$name == ligand),
                                      mode = "out", cutoff = max_steps)
        if (length(all_paths) == 0) {
          return(NULL)
        } else {
          # 验证路径并计算分数
          valid_paths <- list()
          for (path in all_paths) {
            gene_path <- V(subgraph)$name[path]
            valid <- TRUE
            # Receptor should on the cell membrane
            receptor_loc <- Pro_loc_db$location[Pro_loc_db$genesymbol == gene_path[1]]
            if (!"Membrane" %in% receptor_loc) valid <- FALSE
            # Ohter molecules should be intracellular molecules
            for (i in 2:(length(gene_path)-2)) {
              gene <- gene_path[i]
              loc <- Pro_loc_db$location[Pro_loc_db$genesymbol == gene]
              if (!"Intracellular" %in% loc) valid <- FALSE
            }
            # The penultimate molecule should be in Nucleus
            tf_loc <- Pro_loc_db$location[Pro_loc_db$genesymbol == gene_path[length(gene_path)-1]]
            if (!"Nucleus" %in% tf_loc) valid <- FALSE
            if (valid) valid_paths <- c(valid_paths, list(gene_path))
          }
          if(length(valid_paths) == 0){
            return(NULL)
          }else{
            path_scores <- calculate_path_scores(paths = valid_paths, all_edges, max_num = max_num)
            path_score = max(path_scores$Score)
          }
          
          # 计算电路分数
          I_full[c(receptor, ligand)] <- c(1, -1)
          b_reg <- L_trans[-del_idx, ] %*% I_full
          V_reduced <- A_reg_reduced_v %*% b_reg
          V_full <- rep(0, length(nodes))
          V_full[-del_idx] <- V_reduced
          names(V_full) <- nodes
          V_full_norm = (V_full-mean(V_full))/sd(V_full)
          circuit_score = (V_full_norm[receptor] - V_full_norm[ligand]) * circuit_scale_factor # 缩放电路分数
          
          # 动态调整alpha
          if (adjust_alpha && path_score > 0 && circuit_score != 0) {
            # 动态调整alpha
            current_ratio <- path_score / abs(circuit_score)
            alpha <- target_ratio / (target_ratio + current_ratio)
          }
          combined_score <- alpha * path_score + (1 - alpha) * circuit_score
          return(data.table('receptor' = receptor,
                            'ligand' = ligand,
                            'cell_id' = cid,
                            'path_score' = path_score,
                            'circuit_score' = circuit_score,
                            'alpha' = alpha,
                            'combined_score' = combined_score))
        }
      })
      return(rbindlist(inner_results, fill = T))
    }
    ))
  }
  return(results)
}


#' Calculate Triple Cell-Cell Interaction (CCI) with Customized Parameters
#'
#' This function calculates the triple cell-cell interaction scores based on the provided input data.
#' It supports both parallel and sequential processing, and can handle different types of spatial transcriptomics data.
#'
#' @param Input A matrix containing gene expression data.(Recommend normalized data)
#' @param Coordinate A data frame containing the spatial coordinates of spots, with columns "x" and "y".
#' @param st_type A character string indicating the type of spatial transcriptomics data, either "low" or "high". Default is "low".
#' @param relay_signal A data frame containing relay signal information.
#' @param n_permutations An integer specifying the number of permutations for background score calculation. Default is 1000.
#' @param future_globals_maxSize A numeric value specifying the maximum size of global variables in future processing. Default is 8e9.
#' @param n_neighbor An integer specifying the number of neighbors to find. Default is 6.
#' @param parallel A logical value indicating whether to use parallel processing. Default is FALSE.
#' @param n_cores An integer specifying the number of cores to use in parallel processing. If NULL, it will use all available cores minus 1.
#' @param Filter A logical value indicating whether to filter the data. Default is TRUE.
#' @param cell_type_list1 A character vector specifying the first list of cell types. If NULL, it will be inferred from the data.
#' @param cell_type_list2 A character vector specifying the second list of cell types. If NULL, it will be inferred from the data.
#' @param cell_type_list3 A character vector specifying the third list of cell types. If NULL, it will be inferred from the data.
#' @param lrpair_list1 A data frame containing the first list of ligand-receptor pairs. If NULL, it will use the built-in data.
#' @param lrpair_list2 A data frame containing the second list of ligand-receptor pairs. If NULL, it will use the built-in data.
#' @param deconvolution_result A matrix or data frame containing deconvolution results. Required if st_type is "low".
#' @param scRNA_mean_expression A matrix or data frame containing single-cell RNA-seq mean expression data. Required if st_type is "low".
#' @param all_celltype_expression_list A list containing the expression data of each cell type. If NULL, it will be calculated.
#' @param double_lrpairs A data frame containing double ligand-receptor pairs. If NULL, it will be generated.
#' @param verbose A logical value indicating whether to print progress messages. Default is TRUE.
#'
#' @return A data.table containing the calculated triple cell-cell interaction scores.
#'
#' @export
#' @examples
#' result <- calculate_tri_cci_customized(Input = expr_data, Coordinate = coord_data, 
#'                                        relay_signal = relay_data, st_type = "low",
#'                                        deconvolution_result = deconv_data,
#'                                        scRNA_mean_expression = scRNA_data)
#'

calculate_tri_cci_customized = function(Input,
                                        Coordinate, 
                                        st_type = c('low', 'high'), 
                                        relay_signal,
                                        n_permutations = 1000,
                                        future_globals_maxSize = 8e9,
                                        n_neighbor = 6,
                                        parallel = FALSE,
                                        n_cores = NULL,
                                        Filter = TRUE,
                                        cell_type_list1 = NULL, 
                                        cell_type_list2 = NULL, 
                                        cell_type_list3 = NULL, 
                                        lrpair_list1 = NULL,
                                        lrpair_list2 = NULL,
                                        deconvolution_result = NULL, 
                                        scRNA_mean_expression = NULL,
                                        all_celltype_expression_list = NULL,
                                        double_lrpairs = NULL,
                                        verbose = TRUE){
  if (st_type == "low") {
    stopifnot(!is.null(deconvolution_result), !is.null(scRNA_mean_expression))
    if (!identical(colnames(deconvolution_result), colnames(scRNA_mean_expression)))
      stop("Error: Colnames of deconvolution_result and scRNA_mean_expression is different!\n")
  }
  if(verbose) cat("=== Start ===\n")
  if(verbose) cat("=== Step0 Data prepare\n")
  # 构建细胞类型组合
  if (is.null(cell_type_list1) || is.null(cell_type_list2) || is.null(cell_type_list3)) {
    cell_types <- colnames(deconvolution_result)
    if (is.null(cell_type_list1)) cell_type_list1 <- cell_types
    if (is.null(cell_type_list2)) cell_type_list2 <- cell_types
    if (is.null(cell_type_list3)) cell_type_list3 <- cell_types
  }
  celltype_pair <- expand.grid(ct1 = cell_type_list1, ct2 = cell_type_list2, ct3 = cell_type_list3, stringsAsFactors = FALSE)
  if(verbose) cat("    There are", nrow(celltype_pair), "triple cell type pairs to calculate.\n")
  
  # 表达预处理
  expressed_genes <- rownames(Input)[rowSums(Input) > 0]
  Input <- Input[expressed_genes, , drop = FALSE]
  
  if (any(is.null(lrpair_list1), is.null(lrpair_list2))) {
    utils::data("SpaTalk_lrpairs", package = "stMCCC")
    if(is.null(lrpair_list1)) lrpair_list1 <- SpaTalk_lrpairs
    if(is.null(lrpair_list2)) lrpair_list2 <- SpaTalk_lrpairs
  }
  lrpair_list1 <- subset(lrpair_list1, ligand %in% expressed_genes & receptor %in% expressed_genes)
  lrpair_list2 <- subset(lrpair_list2, ligand %in% expressed_genes & receptor %in% expressed_genes)
  
  lrpair_list1$LR <- paste0(lrpair_list1$ligand, "-", lrpair_list1$receptor)
  lrpair_list2$LR <- paste0(lrpair_list2$ligand, "-", lrpair_list2$receptor)
  if(verbose) cat("    There are", nrow(lrpair_list1), "lr1 and", nrow(lrpair_list2), "lr2 pairs to calculate.\n")
  
  if(is.null(double_lrpairs)){
    double_lrpairs <- merge(expand.grid(lr1 = lrpair_list1$LR, lr2 = lrpair_list2$LR, stringsAsFactors = FALSE),
                            lrpair_list1, by.x = "lr1", by.y = "LR")
    double_lrpairs <- merge(double_lrpairs, lrpair_list2, by.x = "lr2", by.y = "LR", suffixes = c("1", "2"))
  }
  double_lrpairs$RL <- paste0(double_lrpairs$receptor1, "-", double_lrpairs$ligand2)
  double_lrpairs$LRLR <- paste0(double_lrpairs$lr1, "-", double_lrpairs$lr2)
  relay_signal$RL <- paste0(relay_signal$receptor, "-", relay_signal$ligand)
  relay_signal <- subset(relay_signal, receptor %in% expressed_genes & ligand %in% expressed_genes)
  double_lrpairs <- subset(double_lrpairs, RL %in% relay_signal$RL)
  if(nrow(double_lrpairs) == 0){cat(" There are no lr-lr pairs to calculate.\n"); return(NULL)}
  if(verbose) cat("    There are", nrow(double_lrpairs), "double_lrpairs to calculate.\n")
  
  if(verbose) cat("=== Step1. Start reconsitution the expression of each cell type.\n")
  spots_celltype_matrix <- Construct_spots_celltype(RCTD_result = deconvolution_result)
  deconvolution_reconsitution <- Reconstruct_RCTD(deconvolution_result, spots_celltype_matrix)
  if(is.null(all_celltype_expression_list)){
    all_celltype_expression_list <- Get_reconstruct_ST_expression_ct_list(Input, scRNA_mean_expression, deconvolution_reconsitution)
    names(all_celltype_expression_list) = colnames(deconvolution_result)
  }
  spots_celltype_matrix = as.data.frame(spots_celltype_matrix)
  
  if(verbose) cat("=== Step2. Start find the neighbors. \n")
  Coordinate_matrix <- as.matrix(Coordinate[, c("x", "y")])
  knn_result <- FNN::get.knn(Coordinate_matrix, k = n_neighbor)
  neighbor_indices <- knn_result$nn.index
  neighbor_names <- apply(neighbor_indices, 1, function(idx) rownames(Coordinate)[idx])
  knn_result_df <- as.data.frame(t(neighbor_names))
  rownames(knn_result_df) <- rownames(Coordinate)
  
  if(verbose) cat("=== Step3. Start calculate interaction score.\n")
  if(parallel){
    if(verbose) cat("    Parallel: \n")
    if(is.null(n_cores)){n_cores = availableCores() - 1}
    
    plan(multisession, workers = n_cores)
    future.globals <- c("celltype_pair", "knn_result_df", "spots_celltype_matrix",
                        "all_celltype_expression_list", "relay_signal", "double_lrpairs",
                        "deconvolution_reconsitution", "n_permutations", "process_chunk","Filter")
    
    library(progressr)
    handlers(global = TRUE)
    with_progress({
      res = future_lapply(1:nrow(celltype_pair), function(g){
        tryCatch({
          process_chunk = function(df_sub, all_spot_pair, relay_df){
            L1_mat <- expr_ct1[df_sub$ligand1, all_spot_pair$spot1, drop = FALSE]
            R1_mat <- expr_ct2[df_sub$receptor1, all_spot_pair$spot2, drop = FALSE]
            L2_mat <- expr_ct2[df_sub$ligand2, all_spot_pair$spot2, drop = FALSE]
            R2_mat <- expr_ct3[df_sub$receptor2, all_spot_pair$spot3, drop = FALSE]
            relay_scores <- relay_df$combined_score[match(df_sub$RL, relay_df$RL)]
            score_mat <- sqrt(L1_mat * R1_mat) * sqrt(L2_mat * R2_mat)
            score_mat = apply(score_mat, 2, function(x){x * relay_scores})
            rm(L1_mat, R1_mat, L2_mat,R2_mat)
            L1_mat_permutation <- expr_ct1[df_sub$ligand1, spot_pair_permutation$spot1, drop = FALSE]
            R1_mat_permutation <- expr_ct2[df_sub$receptor1, spot_pair_permutation$spot2, drop = FALSE]
            L2_mat_permutation <- expr_ct2[df_sub$ligand2, spot_pair_permutation$spot2, drop = FALSE]
            R2_mat_permutation <- expr_ct3[df_sub$receptor2, spot_pair_permutation$spot3, drop = FALSE]
            score_mat_permutation <- sqrt(L1_mat_permutation * R1_mat_permutation) * sqrt(L2_mat_permutation * R2_mat_permutation)
            score_mat_permutation = apply(score_mat_permutation, 2, function(x){x * relay_scores})
            result_df = list()
            for (df_r in 1:nrow(df_sub)) {
              se = score_mat[df_r,]
              re = data.frame(all_spot_pair, 'score' = se)
              re = re[re$score > 0,]
              if(nrow(re) == 0){next} else{
                bg_score = score_mat_permutation[df_r,]
                ecdf_permuted_scores <- ecdf(bg_score) 
                re$pvalue = ecdf_permuted_scores(re$score)
                re$LRLR = df_sub$LRLR[df_r]
                result_df[[df_r]] = as.data.table(re)
              }
            }
            result_df = rbindlist(result_df, fill = T)
            result_df$ct_pair = ct_pair
            result_df = as.data.table(result_df)
            result_df = result_df[result_df$pvalue < 0.05,]
            if(nrow(result_df) == 0 ){return(NULL)}
            return(result_df)
          }
          
          
          
          ct1 <- celltype_pair$ct1[g]; ct2 <- celltype_pair$ct2[g]; ct3 <- celltype_pair$ct3[g]
          ct_pair = paste0(ct1,'-',ct2,'-',ct3)
          p <- progressor(along = 1:nrow(celltype_pair))
          p(sprintf("Processing %s", ct_pair))
          
          cat(' ',g,'.',ct_pair, ' \n')
          knn_ct <- knn_result_df[apply(spots_celltype_matrix, 1, function(x) ct2 %in% x), , drop = FALSE]
          knn_ct <- knn_ct[apply(knn_ct, 1, function(x) all(c(ct1, ct3) %in% unlist(spots_celltype_matrix[x, , drop = FALSE]))), , drop = FALSE]
          cat('    Number of ct2 spots: ',nrow(knn_ct),' \n')
          if (nrow(knn_ct) == 0) return(NULL)
          
          if(nrow(knn_ct) > 500 & Filter == TRUE){
            knn_ct = knn_ct[names(deconvolution_reconsitution[order(deconvolution_reconsitution[,ct2],decreasing = T),ct2][1:500]),]
            cat('    Number of filter ct2 spots: 500 \n')
          }
          
          all_spot_pair <- lapply(1:nrow(knn_ct), function(row_i) {
            sp2 <- rownames(knn_ct)[row_i]
            sp_knn <- knn_ct[row_i, ] #%>% t() %>% as.vector()
            sp1 <- sp_knn[which(apply(sp_knn, 2, function(sp) ct1 %in% spots_celltype_matrix[sp,]))] %>% t() %>% as.vector()
            sp3 <- sp_knn[which(apply(sp_knn, 2, function(sp) ct3 %in% spots_celltype_matrix[sp,]))] %>% t() %>% as.vector()
            if (length(sp1) == 0 || length(sp3) == 0) return(NULL)
            expand.grid(spot1 = sp1, spot2 = sp2, spot3 = sp3, KEEP.OUT.ATTRS = F, stringsAsFactors = FALSE) %>% as.data.table()
          })
          all_spot_pair <- data.table::rbindlist(all_spot_pair)
          if (nrow(all_spot_pair) == 0) return(NULL)
          cat('    Number of all spot pairs: ',nrow(all_spot_pair),' \n')
          
          if(nrow(all_spot_pair) > 1000 & Filter == TRUE){
            all_spot_pair = all_spot_pair[sample(1:nrow(all_spot_pair), 1000),]
            cat('    Number of all filter spot pairs: 1000 \n')
          }
          
          # 构造 permutation 对照
          spots_ct1 <- rownames(spots_celltype_matrix)[apply(spots_celltype_matrix, 1, function(x) ct1 %in% x)]
          spots_ct2 <- rownames(spots_celltype_matrix)[apply(spots_celltype_matrix, 1, function(x) ct2 %in% x)]
          spots_ct3 <- rownames(spots_celltype_matrix)[apply(spots_celltype_matrix, 1, function(x) ct3 %in% x)]
          spot_pair_permutation <- data.frame(
            spot1 = sample(spots_ct1, n_permutations, replace = TRUE),
            spot2 = sample(spots_ct2, n_permutations, replace = TRUE),
            spot3 = sample(spots_ct3, n_permutations, replace = TRUE),
            stringsAsFactors = FALSE
          )
          # 表达矩阵提取
          expr_ct1 <- as.matrix(all_celltype_expression_list[[ct1]])
          expr_ct2 <- as.matrix(all_celltype_expression_list[[ct2]])
          expr_ct3 <- as.matrix(all_celltype_expression_list[[ct3]])
          relay_df <- subset(relay_signal, cell_id == ct2)
          df_sub <- double_lrpairs[double_lrpairs$RL %in% relay_df$RL, , drop = FALSE]
          
          if (nrow(df_sub) == 0) return(NULL)
          expr_genes2 <- intersect(rownames(expr_ct2)[rowSums(expr_ct2) > 0], unique(c(df_sub$receptor1, df_sub$ligand2)))
          df_sub <- df_sub[df_sub$receptor1 %in% expr_genes2 & df_sub$ligand2 %in% expr_genes2, , drop = FALSE]
          expr_genes1 <- intersect(rownames(expr_ct1)[rowSums(expr_ct1) > 0], unique(df_sub$ligand1))
          expr_genes3 <- intersect(rownames(expr_ct3)[rowSums(expr_ct3) > 0], unique(df_sub$receptor2))
          df_sub <- df_sub[df_sub$ligand1 %in% expr_genes1 & df_sub$receptor2 %in% expr_genes3, , drop = FALSE]
          if (nrow(df_sub) == 0) return(NULL)
          cat('    Number of all lr-lr pairs: ',nrow(df_sub),' \n')
          
          chunk_size <- 1000
          n_chunks <- ceiling(nrow(df_sub) / chunk_size)
          result_list <- vector("list", n_chunks)
          cat('    n_chunks: ',n_chunks, ' \n')
          for (chunk in seq_len(n_chunks)) {
            idx_start <- (chunk - 1) * chunk_size + 1
            idx_end <- min(chunk * chunk_size, nrow(df_sub))
            df_sub_chunk <- df_sub[idx_start:idx_end, , drop = FALSE]
            chunk_result <- process_chunk(df_sub_chunk, all_spot_pair, relay_df)
            # 保存结果并释放内存
            result_list[[chunk]] <- as.data.table(chunk_result)
            rm(df_sub_chunk, chunk_result)
          }
          result_list = rbindlist(result_list, fill = T)
          rm(knn_ct, all_spot_pair,expr_ct1,expr_ct2,expr_ct3)
          return(as.data.table(result_list))
        }, error = function(e) {
          message(sprintf("Error in celltype_pair %d: %s", g, e$message))
          return(NULL)
        })
      },future.seed = TRUE)
    })
    plan(sequential)
  }else{
    if(verbose) cat("   Sequential: \n")
    res = lapply(1:nrow(celltype_pair), function(g){
      process_chunk = function(df_sub, all_spot_pair, relay_df){
        L1_mat <- expr_ct1[df_sub$ligand1, all_spot_pair$spot1, drop = FALSE]
        R1_mat <- expr_ct2[df_sub$receptor1, all_spot_pair$spot2, drop = FALSE]
        L2_mat <- expr_ct2[df_sub$ligand2, all_spot_pair$spot2, drop = FALSE]
        R2_mat <- expr_ct3[df_sub$receptor2, all_spot_pair$spot3, drop = FALSE]
        relay_scores <- relay_df$combined_score[match(df_sub$RL, relay_df$RL)]
        score_mat <- sqrt(L1_mat * R1_mat) * sqrt(L2_mat * R2_mat)
        score_mat = apply(score_mat, 2, function(x){x * relay_scores})
        rm(L1_mat, R1_mat, L2_mat,R2_mat)
        L1_mat_permutation <- expr_ct1[df_sub$ligand1, spot_pair_permutation$spot1, drop = FALSE]
        R1_mat_permutation <- expr_ct2[df_sub$receptor1, spot_pair_permutation$spot2, drop = FALSE]
        L2_mat_permutation <- expr_ct2[df_sub$ligand2, spot_pair_permutation$spot2, drop = FALSE]
        R2_mat_permutation <- expr_ct3[df_sub$receptor2, spot_pair_permutation$spot3, drop = FALSE]
        score_mat_permutation <- sqrt(L1_mat_permutation * R1_mat_permutation) * sqrt(L2_mat_permutation * R2_mat_permutation)
        score_mat_permutation = apply(score_mat_permutation, 2, function(x){x * relay_scores})
        result_df = list()
        for (df_r in 1:nrow(df_sub)) {
          se = score_mat[df_r,]
          re = data.frame(all_spot_pair, 'score' = se)
          re = re[re$score > 0,]
          if(nrow(re) == 0){next} else{
            bg_score = score_mat_permutation[df_r,]
            ecdf_permuted_scores <- ecdf(bg_score) 
            re$pvalue = ecdf_permuted_scores(re$score)
            re$LRLR = df_sub$LRLR[df_r]
            result_df[[df_r]] = as.data.table(re)
          }
        }
        result_df = rbindlist(result_df, fill = T)
        result_df$ct_pair = ct_pair
        result_df = as.data.table(result_df)
        result_df = result_df[result_df$pvalue < 0.05,]
        if(nrow(result_df) == 0 ){return(NULL)}
        return(result_df)
      }
      
      
      ct1 <- celltype_pair$ct1[g]; ct2 <- celltype_pair$ct2[g]; ct3 <- celltype_pair$ct3[g]
      ct_pair = paste0(ct1,'-',ct2,'-',ct3)
      cat(' ',g,'.',ct_pair, ' \n')
      knn_ct <- knn_result_df[apply(spots_celltype_matrix, 1, function(x) ct2 %in% x), , drop = FALSE]
      knn_ct <- knn_ct[apply(knn_ct, 1, function(x) all(c(ct1, ct3) %in% unlist(spots_celltype_matrix[x, , drop = FALSE]))), , drop = FALSE]
      cat('    Number of ct2 spots: ',nrow(knn_ct),' \n')
      if (nrow(knn_ct) == 0) return(NULL)
      
      if(nrow(knn_ct) > 500 & Filter == TRUE){
        knn_ct = knn_ct[names(deconvolution_reconsitution[order(deconvolution_reconsitution[,ct2],decreasing = T),ct2][1:500]),]
        cat('    Number of filter ct2 spots: 500 \n')
      }
      
      all_spot_pair <- lapply(1:nrow(knn_ct), function(row_i) {
        sp2 <- rownames(knn_ct)[row_i]
        sp_knn <- knn_ct[row_i, ] #%>% t() %>% as.vector()
        sp1 <- sp_knn[which(apply(sp_knn, 2, function(sp) ct1 %in% spots_celltype_matrix[sp,]))] %>% t() %>% as.vector()
        sp3 <- sp_knn[which(apply(sp_knn, 2, function(sp) ct3 %in% spots_celltype_matrix[sp,]))] %>% t() %>% as.vector()
        if (length(sp1) == 0 || length(sp3) == 0) return(NULL)
        expand.grid(spot1 = sp1, spot2 = sp2, spot3 = sp3, KEEP.OUT.ATTRS = F, stringsAsFactors = FALSE) %>% as.data.table()
      })
      all_spot_pair <- data.table::rbindlist(all_spot_pair)
      if (nrow(all_spot_pair) == 0) return(NULL)
      cat('    Number of all spot pairs: ',nrow(all_spot_pair),' \n')
      
      if(nrow(all_spot_pair) > 1000 & Filter == TRUE){
        all_spot_pair = all_spot_pair[sample(1:nrow(all_spot_pair), 1000),]
        cat('    Number of all filter spot pairs: 1000 \n')
      }
      
      # 构造 permutation 对照
      spots_ct1 <- rownames(spots_celltype_matrix)[apply(spots_celltype_matrix, 1, function(x) ct1 %in% x)]
      spots_ct2 <- rownames(spots_celltype_matrix)[apply(spots_celltype_matrix, 1, function(x) ct2 %in% x)]
      spots_ct3 <- rownames(spots_celltype_matrix)[apply(spots_celltype_matrix, 1, function(x) ct3 %in% x)]
      spot_pair_permutation <- data.frame(
        spot1 = sample(spots_ct1, n_permutations, replace = TRUE),
        spot2 = sample(spots_ct2, n_permutations, replace = TRUE),
        spot3 = sample(spots_ct3, n_permutations, replace = TRUE),
        stringsAsFactors = FALSE
      )
      # 表达矩阵提取
      expr_ct1 <- as.matrix(all_celltype_expression_list[[ct1]])
      expr_ct2 <- as.matrix(all_celltype_expression_list[[ct2]])
      expr_ct3 <- as.matrix(all_celltype_expression_list[[ct3]])
      relay_df <- subset(relay_signal, cell_id == ct2)
      df_sub <- double_lrpairs[double_lrpairs$RL %in% relay_df$RL, , drop = FALSE]
      
      if (nrow(df_sub) == 0) return(NULL)
      expr_genes2 <- intersect(rownames(expr_ct2)[rowSums(expr_ct2) > 0], unique(c(df_sub$receptor1, df_sub$ligand2)))
      df_sub <- df_sub[df_sub$receptor1 %in% expr_genes2 & df_sub$ligand2 %in% expr_genes2, , drop = FALSE]
      expr_genes1 <- intersect(rownames(expr_ct1)[rowSums(expr_ct1) > 0], unique(df_sub$ligand1))
      expr_genes3 <- intersect(rownames(expr_ct3)[rowSums(expr_ct3) > 0], unique(df_sub$receptor2))
      df_sub <- df_sub[df_sub$ligand1 %in% expr_genes1 & df_sub$receptor2 %in% expr_genes3, , drop = FALSE]
      if (nrow(df_sub) == 0) return(NULL)
      cat('    Number of all lr-lr pairs: ',nrow(df_sub),' \n')
      chunk_size <- 1000
      n_chunks <- ceiling(nrow(df_sub) / chunk_size)
      result_list <- vector("list", n_chunks)
      cat('    n_chunks: ',n_chunks, ' \n')
      for (chunk in seq_len(n_chunks)) {
        idx_start <- (chunk - 1) * chunk_size + 1
        idx_end <- min(chunk * chunk_size, nrow(df_sub))
        df_sub_chunk <- df_sub[idx_start:idx_end, , drop = FALSE]
        chunk_result <- process_chunk(df_sub_chunk, all_spot_pair, relay_df)
        # 保存结果并释放内存
        result_list[[chunk]] <- as.data.table(chunk_result)
        rm(df_sub_chunk, chunk_result)
      }
      result_list = rbindlist(result_list, fill = T)
      rm(knn_ct, all_spot_pair,expr_ct1,expr_ct2,expr_ct3)
      return(as.data.table(result_list))
    })
    
  }
  
  res = rbindlist(res, use.names = TRUE, fill = TRUE)
  return(res)
}




