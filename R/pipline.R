#' Validate input data for calculate_all_lrlr function.
#'
#' @param Input Normalized expression data matrix.
#' @param Coordinate Coordinate matrix.
#' @param st_type Spatial transcriptome data type ('low' or 'high').
#' @param deconvolution_result Deconvolution result matrix.
#' @param scRNA_mean_expression Mean expression matrix from scRNA-seq.
#' @return TRUE if all inputs are valid, otherwise stop with an error message.
#' @export
validate_inputs <- function(Input, Coordinate, st_type, deconvolution_result = NULL, scRNA_mean_expression = NULL) {
  if (!is.matrix(Input)) {
    stop("Input must be a matrix.")
  }
  if (!is.matrix(Coordinate)) {
    stop("Coordinate must be a matrix.")
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



#' Calculate double cell cell communication
#'
#' @param ct1 Sender cell type name
#' @param ct2 Receiver cell type name
#' @param spots_celltype_matrix Result from the function 'Reconstruct_RCTD'.
#' @param all_celltype_expression_list Result from the function 'Get_reconstruct_ST_expression_ct_list'
#' @param double_spots_pairs Result from the function 'Get_double_spots_pairs'
#' @param lrpairs Database of lrpairs
#' @param n_no_interact Background interaction spots pairs
#' @param p_cutoff  A cutoff for zero-inflated gamma test
#' 
#' @return A.
#' @export 
#' @examples
#' all_celltype_expression_list = Calculate_double_CCC(st_exp_matrix,ct1,ct2, spots_celltype_matrix,all_celltype_expression_list, double_spots_pairs,lrpairs)
Calculate_double_CCC = function(ct1,ct2, spots_celltype_matrix, all_celltype_expression_list, double_spots_pairs,
                     lrpairs, n_no_interact = 10000, p_cutoff = 0.95){
  spots1_remain = rownames(spots_celltype_matrix)[apply(spots_celltype_matrix,1,function(x){(ct1 %in% x)})]
  spots2_remain = rownames(spots_celltype_matrix)[apply(spots_celltype_matrix,1,function(x){(ct2 %in% x)})]
  if((length(spots1_remain) == 0) | (length(spots2_remain) == 0)){stop(paste0(ct1, ' and ',ct2, ' have no interaction spots!'))}
  D_matrix_df_interact = double_spots_pairs$interaction
  D_matrix_df_no_interact = double_spots_pairs$no_interaction

  D_matrix_df_interact_filter = D_matrix_df_interact[(D_matrix_df_interact$Sender %in% spots1_remain) & (D_matrix_df_interact$Receiver %in% spots2_remain),]
  if(nrow(D_matrix_df_interact_filter) == 0){stop(paste0('D_matrix_df_interact_filter have 0 row!'))}
  D_matrix_df_no_interact_filter = D_matrix_df_no_interact[(D_matrix_df_no_interact$Sender %in% spots1_remain) & (D_matrix_df_no_interact$Receiver %in% spots2_remain),]
  if(nrow(D_matrix_df_no_interact_filter) == 0){stop(paste0('D_matrix_df_no_interact_filter have 0 row!'))}
  if (dim(D_matrix_df_no_interact_filter)[1] > n_no_interact) {
    D_matrix_df_no_interact_filter_n = D_matrix_df_no_interact_filter[sample(c(1:dim(D_matrix_df_no_interact_filter)[1]), n_no_interact),]
  }else{
    D_matrix_df_no_interact_filter_n = D_matrix_df_no_interact_filter
  }
  print(paste0('Start calculate ',ct1,'_',ct2,' LR CCC!'))
  all_lr_results = list()
  for (j in 1:dim(lrpairs)[1]) {
    if((j %% 100) == 0){print(j)}
    ligand = lrpairs$ligand[j]
    receptor = lrpairs$receptor[j]
    ligand_exp = all_celltype_expression_list[[ct1]][ligand,D_matrix_df_interact_filter$Sender]
    if (sum(ligand_exp)==0) {next}
    receptor_exp = all_celltype_expression_list[[ct2]][receptor,D_matrix_df_interact_filter$Receiver]
    if (sum(receptor_exp)==0) {next}
    outer_score = ligand_exp * receptor_exp
    if (sum(outer_score) == 0){next}
    
    ligand_exp_null = all_celltype_expression_list[[ct1]][ligand,D_matrix_df_no_interact_filter_n$Sender]
    receptor_exp_null = all_celltype_expression_list[[ct2]][receptor,D_matrix_df_no_interact_filter_n$Receiver]
    null_distribution = ligand_exp_null * receptor_exp_null
    null_distribution = as.vector(null_distribution)
    ## 检验
    if (quantile(null_distribution, p_cutoff) == 0) {
      re = data.frame('ct1' = ct1, 'ct2' = ct2, D_matrix_df_interact_filter, 'ligand' = ligand, 'receptor' = receptor, 'outer_score' = outer_score, 'p_value' = 1)
      re$p_value[re$outer_score > 0] = 0
    }else{
      zero_pro = mean(null_distribution == 0)
      null_distribution_pos = null_distribution[null_distribution > 0]
      null_mean = mean(null_distribution_pos)
      null_var = var(null_distribution_pos)
      p = 1 - pgamma(outer_score, shape = (null_mean * null_mean)/null_var, scale = null_var/null_mean)
      p = p * (1-zero_pro)
      re = data.frame(D_matrix_df_interact_filter, 'ct1' = ct1, 'ct2' = ct2, 'ligand' = ligand, 'receptor' = receptor, 'outer_score' = outer_score, 'p_value' = p)
      re$p_value[re$outer_score == 0] = 1
    }
    all_lr_results[[j]] = as.data.table(re)
  }
  all_lr_results = rbindlist(all_lr_results, fill = TRUE)
  all_lr_results = all_lr_results[all_lr_results$outer_score > 0,]
  return(all_lr_results)
}



#' Calculate triple cell cell communication
#'
#' @param ct1 Sender cell type name
#' @param ct2 Mid cell type name
#' @param ct3 Receiver cell type name
#' @param spots_celltype_matrix Result from the function 'Reconstruct_RCTD'.
#' @param double_spots_pairs Result from the function 'Get_double_spots_pairs'
#' @param mid_spots_cut Database of lrpairs
#' @param bg_mid_spots_cut_min Background interaction spots pairs
#' @param p_cutoff  A cutoff for zero-inflated gamma test
#' 
#' @return A.
#' @export 
#' @examples
#' all_celltype_expression_list = Calculate_double_CCC(st_exp_matrix,ct1,ct2, spots_celltype_matrix,all_celltype_expression_list, double_spots_pairs,lrpairs)
get_triple_CCC_pairs = function(ct1, ct2, ct3, spots_celltype_matrix, 
                                double_spots_pairs, 
                                mid_spots_cut = 5,
                                bg_mid_spots_cut_min = 5, bg_mid_spots_cut_max = 100,
                                bg_spot_pair_cut_max = 10000){
  D_matrix_df_interact = double_spots_pairs$interaction
  D_matrix_df_no_interact = double_spots_pairs$no_interaction
  spots1_remain = rownames(spots_celltype_matrix)[apply(spots_celltype_matrix,1,function(x){(ct1 %in% x)})]
  spots2_remain = rownames(spots_celltype_matrix)[apply(spots_celltype_matrix,1,function(x){(ct2 %in% x)})]
  spots3_remain = rownames(spots_celltype_matrix)[apply(spots_celltype_matrix,1,function(x){(ct3 %in% x)})]
  D_matrix_df_interact_12 = D_matrix_df_interact[(D_matrix_df_interact$Sender %in% spots1_remain) & (D_matrix_df_interact$Receiver %in% spots2_remain),]
  D_matrix_df_interact_23 = D_matrix_df_interact[(D_matrix_df_interact$Sender %in% spots2_remain) & (D_matrix_df_interact$Receiver %in% spots3_remain),]
  spots2_both = intersect(names(table(D_matrix_df_interact_12$Receiver)), names(table(D_matrix_df_interact_23$Sender)))
  if (length(spots2_both) < mid_spots_cut) {stop(paste0(ct1,'_',ct2,'_',ct3,': Trere interacted ct2 was less than ',mid_spots_cut,'.'))}
  
  # 拼接互作spots结果
  D_matrix_df_interact_filter_12 = D_matrix_df_interact_12[D_matrix_df_interact_12$Receiver %in% spots2_both,]
  D_matrix_df_interact_filter_23 = D_matrix_df_interact_23[D_matrix_df_interact_23$Sender %in% spots2_both,]
  D_matrix_df_interact_filter_123 = list()
  for (j in 1:dim(D_matrix_df_interact_filter_12)[1]) {
    mid = D_matrix_df_interact_filter_12$Receiver[j]
    D_23 = D_matrix_df_interact_filter_23[D_matrix_df_interact_filter_23$Sender == mid,]
    re = cbind(D_matrix_df_interact_filter_12[rep(j, dim(D_23)[1]),], D_23)
    colnames(re) = c("Sender","mid","distance1","mid1","Receiver", "distance2")
    re$mid1 = NULL
    D_matrix_df_interact_filter_123[[j]] = as.data.table(re)
  }
  D_matrix_df_interact_filter_123 = rbindlist(D_matrix_df_interact_filter_123, fill = TRUE)
  rm(D_matrix_df_interact_filter_12, D_matrix_df_interact_filter_23)
  print(paste0(ct1,'_',ct2,'_',ct3,' finished get D_matrix_df_interact_filter_123'))
  
  # 构建背景的spots组合
  D_matrix_df_no_interact_filter_12 = D_matrix_df_no_interact[(D_matrix_df_no_interact$Sender %in% spots1_remain) & (D_matrix_df_no_interact$Receiver %in% spots2_remain),]
  D_matrix_df_no_interact_filter_23 = D_matrix_df_no_interact[(D_matrix_df_no_interact$Sender %in% spots2_remain) & (D_matrix_df_no_interact$Receiver %in% spots3_remain),]
  spots2_remain_filter2 = intersect(names(table(D_matrix_df_no_interact_filter_12$Receiver)), names(table(D_matrix_df_no_interact_filter_23$Sender)))
  if (length(spots2_remain_filter2) < bg_mid_spots_cut_min) {stop(paste0(ct1,'_',ct2,'_',ct3,': Trere background interacted ct2 was less than ',bg_mid_spots_cut_min,'.'))}
  
  if (nrow(D_matrix_df_no_interact_filter_12) > bg_spot_pair_cut_max) {
    D_matrix_df_no_interact_filter_12 = D_matrix_df_no_interact_filter_12[sample(1:nrow(D_matrix_df_no_interact_filter_12), bg_spot_pair_cut_max),]
  }
  if (nrow(D_matrix_df_no_interact_filter_23) > bg_spot_pair_cut_max) {
    D_matrix_df_no_interact_filter_23 = D_matrix_df_no_interact_filter_23[sample(1:nrow(D_matrix_df_no_interact_filter_23), bg_spot_pair_cut_max),]
  }
  spots2_remain_filter2 = intersect(names(table(D_matrix_df_no_interact_filter_12$Receiver)), names(table(D_matrix_df_no_interact_filter_23$Sender)))
  if (length(spots2_remain_filter2) > bg_mid_spots_cut_max){
    spots2_remain_filter2 = spots2_remain_filter2[sample(1:length(spots2_remain_filter2), bg_mid_spots_cut_max,replace = F)]
  }
  D_matrix_df_no_interact_filter_12 = D_matrix_df_no_interact_filter_12[D_matrix_df_no_interact_filter_12$Receiver %in% spots2_remain_filter2,]
  D_matrix_df_no_interact_filter_23 = D_matrix_df_no_interact_filter_23[D_matrix_df_no_interact_filter_23$Sender %in% spots2_remain_filter2,]
  
  D_matrix_df_no_interact_filter_123_sereis = list()
  for (j in 1:nrow(D_matrix_df_no_interact_filter_12)) {
    mid = D_matrix_df_no_interact_filter_12$Receiver[j]
    D_23 = D_matrix_df_no_interact_filter_23[D_matrix_df_no_interact_filter_23$Sender == mid,]
    re = cbind(D_matrix_df_no_interact_filter_12[rep(j, dim(D_23)[1]),], D_23)
    colnames(re) = c("Sender","mid","distance1","mid1","Receiver", "distance2")
    re$mid1 = NULL
    D_matrix_df_no_interact_filter_123_sereis[[j]] = as.data.table(re)
  }
  D_matrix_df_no_interact_filter_123_sereis = rbindlist(D_matrix_df_no_interact_filter_123_sereis, fill = TRUE)
  rm(D_matrix_df_no_interact_filter_12,D_matrix_df_no_interact_filter_23)
  
  if (nrow(D_matrix_df_no_interact_filter_123_sereis) > bg_spot_pair_cut_max) {
    D_matrix_df_no_interact_filter_123_sereis = D_matrix_df_no_interact_filter_123_sereis[sample(1:dim(D_matrix_df_no_interact_filter_123_sereis)[1], 1000),]
    spots2_remain_filter2 = unique(D_matrix_df_no_interact_filter_123_sereis$mid)
  }
  D_matrix_df_no_interact_filter_123_sereis = D_matrix_df_no_interact_filter_123_sereis[D_matrix_df_no_interact_filter_123_sereis$mid %in% spots2_remain_filter2,]
  
  print(paste0(ct1,'_',ct2,'_',ct3,' finished get D_matrix_df_no_interact_filter_123_sereis'))
  re = list(D_matrix_df_interact_filter_123,D_matrix_df_no_interact_filter_123_sereis)
  names(re) = c('interaction_spots_pairs', 'background_no_interaction_spots_pairs')
  return(re)
}


#' Prepare intracelluar Receptor Liagnd siganal database
#'
#' @param lrpairs1 Provide the preceding LR pairs. A data frame containing two columns with the name 'Ligand' and 'Receptor'
#' @param lrpairs2 Provide the following LR pairs. A data frame containing two columns with the name 'Ligand' and 'Receptor'
#' @param rl_database Database of rl signal database
#' @param min_network_gene_num Background interaction spots pairs
#' 
#' @return A.
#' @export 
#' @examples
#' all_celltype_expression_list = Calculate_double_CCC(st_exp_matrix,ct1,ct2, spots_celltype_matrix,all_celltype_expression_list, double_spots_pairs,lrpairs)

Prepare_lr_rl_lr <- function(lrpairs1, lrpairs2, rl_database, min_network_gene_num = 35, min_cutoff = 5){
  rl_database_fil = rl_database[rl_database$receptor %in% lrpairs1$receptor & rl_database$ligand %in% lrpairs2$ligand,]
  rl_database_fil = rl_database_fil[rl_database_fil$network_gene_num2 > min_network_gene_num,]
  rl_database_fil = as.data.frame(rl_database_fil)
  if(nrow(rl_database_fil) < min_cutoff){
    stop(paste0('Receptor-Liagnd may have no intracelluar signal in rl_database(less than',min_cutoff,')'))
  }
  RL_mid_genes = separate_rows(rl_database_fil, network_genes2, sep = "_")
  RL_mid_genes = RL_mid_genes[,c(1,2,6,7)]
  colnames(RL_mid_genes)[3] = 'network_genes'
  
  merge_lrpairs = merge(lrpairs1, lrpairs2, by = NULL)
  colnames(merge_lrpairs) = c('ligand1','receptor1','ligand2','receptor2')
  merge_lrpairs$rl = paste0(merge_lrpairs$receptor1,"_",merge_lrpairs$ligand2)
  merge_lrpairs = merge_lrpairs[merge_lrpairs$rl %in% rl_database_fil$rl,]
  results = list('Merged_lipairs' = merge_lrpairs, 'Intracelluar_siganal_genes'=RL_mid_genes)
  return(results)
}








#' Prepare a object to calculate triple cell cell communication
#'
#' @param ct1 Sender cell type name
#' @param ct2 Mid cell type name
#' @param ct3 Receiver cell type name
#' @param all_celltype_expression_list Result from the function 'Get_reconstruct_ST_expression_ct_list'.
#' @param triple_interaction_pairs Result from the function 'Get_double_spots_pairs'
#' @param lrpairs1 Provide the preceding LR pairs. A data frame containing two columns with the name 'Ligand' and 'Receptor'
#' @param lrpairs2 Provide the following LR pairs. A data frame containing two columns with the name 'Ligand' and 'Receptor'
#' @param rl_database Database of rl signal database
#' @param bg_mid_spots_cut_min Background interaction spots pairs
#' @param p_cutoff  A cutoff for zero-inflated gamma test
#' 
#' @return A.
#' @export 
#' @examples
#' all_celltype_expression_list = Calculate_double_CCC(st_exp_matrix,ct1,ct2, spots_celltype_matrix,all_celltype_expression_list, double_spots_pairs,lrpairs)

Prepare_calculate_triple_CCC <- function(ct1, ct2, ct3, all_celltype_expression_list, 
                                 triple_interaction_pairs,
                                 p_cutoff = 0.95){
  interaction_spots_pairs = triple_interaction_pairs$interaction_spots_pairs
  bg_spots_pairs = triple_interaction_pairs$background_no_interaction_spots_pairs
  
  # 从上一步互作的pairs中得到中间细胞的像个基因表达
  ct2_spots = unique(interaction_spots_pairs$mid)
  inner_exp = cbind(RL_mid_genes, as.matrix(all_celltype_expression_list[[ct2]][RL_mid_genes$network_gene, ct2_spots]))
  
  rl_remain = unique(inner_exp$rl)
  mid_spot_inner_score = list()
  for (i in 1:length(rl_remain)) {
    if ((i %% 100) == 0) {print(i)}
    x = inner_exp[inner_exp$rl == rl_remain[i],]
    x = x[,-c(1:2)] %>% as.matrix()
    score = apply(x, 2, function(x){
      sum(x)})
    re = data.frame('rl' = rl_remain[i], t(as.data.frame(score)))
    mid_spot_inner_score[[i]] = as.data.table(re)
  }
  mid_spot_inner_score = rbindlist(mid_spot_inner_score, fill = T)
  mid_spot_inner_score = as.data.frame(mid_spot_inner_score)
  mid_spot_inner_score = mid_spot_inner_score[,-1]
  rownames(mid_spot_inner_score) = rl_remain
  colnames(mid_spot_inner_score) = ct2_spots
  
  
  bg_spots = unique(bg_spots_pairs$mid)
  bg_inner_exp = cbind(RL_mid_genes, as.matrix(all_celltype_expression_list[[ct2]][RL_mid_genes$network_gene, bg_spots]))
  bg_mid_spot_inner_score = list()
  for (i in 1:length(rl_remain)) {
    if ((i %% 100) == 0) {print(i)}
    x = bg_inner_exp[bg_inner_exp$rl == rl_remain[i],]
    x = x[,-c(1:2)] %>% as.matrix()
    score = apply(x, 2, function(x){
      sum(x)
    })
    re = data.frame('rl' = rl_remain[i], t(as.data.frame(score)))
    bg_mid_spot_inner_score[[i]] = as.data.table(re)
  }
  bg_mid_spot_inner_score = rbindlist(bg_mid_spot_inner_score, fill = T)
  bg_mid_spot_inner_score = as.data.frame(bg_mid_spot_inner_score)
  bg_mid_spot_inner_score = bg_mid_spot_inner_score[,-1]
  rownames(bg_mid_spot_inner_score) = rl_remain
  colnames(bg_mid_spot_inner_score) = bg_spots
  
  
  # 需要进行过滤，大部分lr和lr没有串联证据
  merge_lrpairs_filter = merge_lrpairs[merge_lrpairs$RL %in% rl_remain,]
  x1 = all_celltype_expression_list[[ct1]][merge_lrpairs_filter$ligand1, interaction_spots_pairs$Sender]
  merge_lrpairs_filter = merge_lrpairs_filter[which(rowSums(x1) != 0),]
  x2 = all_celltype_expression_list[[ct2]][merge_lrpairs_filter$receptor1, interaction_spots_pairs$mid]
  merge_lrpairs_filter = merge_lrpairs_filter[which(rowSums(x2) != 0),]
  x3 = all_celltype_expression_list[[ct2]][merge_lrpairs_filter$ligand2, interaction_spots_pairs$mid]
  merge_lrpairs_filter = merge_lrpairs_filter[which(rowSums(x3) != 0),]
  x4 = all_celltype_expression_list[[ct3]][merge_lrpairs_filter$receptor2, interaction_spots_pairs$Receiver]
  merge_lrpairs_filter = merge_lrpairs_filter[which(rowSums(x4) != 0),]
  
  all_re = list()
  for (i in 1:nrow(merge_lrpairs_filter)) {
    if((i %% 100) == 0){print(i)}
    exp_ligand1 = all_celltype_expression_list[[ct1]][merge_lrpairs_filter$ligand1[i],interaction_spots_pairs$Sender]
    exp_receptor1 = all_celltype_expression_list[[ct2]][merge_lrpairs_filter$receptor1[i],interaction_spots_pairs$mid]
    exp_ligand2 = all_celltype_expression_list[[ct2]][merge_lrpairs_filter$ligand2[i],interaction_spots_pairs$mid]
    exp_receptor2 = all_celltype_expression_list[[ct3]][merge_lrpairs_filter$receptor2[i],interaction_spots_pairs$Receiver]
    RL_score = mid_spot_inner_score[merge_lrpairs_filter$RL[i],interaction_spots_pairs$mid] %>% t() %>% as.vector()
    final_score = exp_ligand1 * exp_receptor1 * exp_ligand2 * exp_receptor2 * RL_score * 10000
    final_score2 = which(final_score > 0)
    if (length(final_score2) > 0) {
      re = data.frame(merge_lrpairs_filter[rep(i,length(final_score2)),], interaction_spots_pairs[final_score2,], 'score' = final_score[final_score2], 'p'=1)
    }else{
      next
    }
    exp_ligand1_null = all_celltype_expression_list[[ct1]][merge_lrpairs_filter$ligand1[i], bg_spots_pairs$Sender]
    exp_receptor1_null = all_celltype_expression_list[[ct2]][merge_lrpairs_filter$receptor1[i], bg_spots_pairs$mid]
    exp_ligand2_null = all_celltype_expression_list[[ct2]][merge_lrpairs_filter$ligand2[i], bg_spots_pairs$mid]
    exp_receptor2_null = all_celltype_expression_list[[ct3]][merge_lrpairs_filter$receptor2[i], bg_spots_pairs$Receiver]
    RL_score_null = bg_mid_spot_inner_score[merge_lrpairs_filter$RL[i],bg_spots_pairs$mid] %>% t() %>% as.vector()
    final_score_null = exp_ligand1_null * exp_receptor1_null * exp_ligand2_null * exp_receptor2_null * RL_score_null * 10000
    if (quantile(final_score_null,0.95) == 0) {
      re$p = 0
      all_re[[i]] = as.data.table(re)
    }else{
      zero_pro = mean(final_score_null == 0)
      null_distribution_pos = final_score_null[final_score_null > 0]
      null_mean = mean(null_distribution_pos)
      null_var = var(null_distribution_pos)
      p = 1 - pgamma(re$p, shape = (null_mean * null_mean)/null_var, scale = null_var/null_mean)
      p = p * (1-zero_pro)
      re$p = p
      all_re[[i]] = as.data.table(re)
    }
  }
  return(all_re)
}






#' Calculate triple cell cell communication
#'
#' @param ct1 Sender cell type name
#' @param ct2 Mid cell type name
#' @param ct3 Receiver cell type name
#' @param spots_celltype_matrix Result from the function 'Reconstruct_RCTD'.
#' @param double_spots_pairs Result from the function 'Get_double_spots_pairs'
#' @param mid_spots_cut Database of lrpairs
#' @param bg_mid_spots_cut_min Background interaction spots pairs
#' @param p_cutoff  A cutoff for zero-inflated gamma test
#' 
#' @return A.
#' @export 
#' @examples
#' all_celltype_expression_list = Calculate_double_CCC(st_exp_matrix,ct1,ct2, spots_celltype_matrix,all_celltype_expression_list, double_spots_pairs,lrpairs)

Calculate_triple_CCC <- function(ct1, ct2, ct3, all_celltype_expression_list, 
         triple_interaction_pairs,
         lr_rl_lr_list){
  interaction_spots_pairs = triple_interaction_pairs$interaction_spots_pairs
  bg_spots_pairs = triple_interaction_pairs$background_no_interaction_spots_pairs
  Merged_lrpairs = lr_rl_lr_list$Merged_lipairs
  Intracelluar_siganal_genes = lr_rl_lr_list$Intracelluar_siganal_genes
  
  # 根据表达中存在的gene，进一步过滤merged lrpairs interacelluar signal gene
  all_genes = rownames(all_celltype_expression_list[[ct1]])
  Intracelluar_siganal_genes = Intracelluar_siganal_genes[Intracelluar_siganal_genes$network_genes %in% all_genes,]
  rl_remain = unique(Intracelluar_siganal_genes$rl)
  Merged_lrpairs = Merged_lrpairs[Merged_lrpairs$rl %in% rl_remain,]
  
  ### 根据lrpairs1 和 lrpairs2 的任意基因表达情况再进行过滤
  x1 = all_celltype_expression_list[[ct1]][unique(Merged_lrpairs$ligand1), interaction_spots_pairs$Sender] %>% as.matrix()
  x1 = rownames(x1)[rowSums(x1) == 0]
  
  x2 = all_celltype_expression_list[[ct2]][unique(Merged_lrpairs$receptor1), interaction_spots_pairs$mid] %>% as.matrix()
  x2 = rownames(x2)[rowSums(x2) == 0]
  
  x3 = all_celltype_expression_list[[ct2]][unique(Merged_lrpairs$ligand2), interaction_spots_pairs$mid] %>% as.matrix()
  x3 = rownames(x3)[rowSums(x3) == 0]
  
  x4 = all_celltype_expression_list[[ct3]][unique(Merged_lrpairs$receptor2), interaction_spots_pairs$Receiver] %>% as.matrix()
  x4 = rownames(x4)[rowSums(x4) == 0]
  
  Merged_lrpairs = Merged_lrpairs[!(Merged_lrpairs$ligand1 %in% x1) & !(Merged_lrpairs$receptor1 %in% x2) & 
                                    !(Merged_lrpairs$ligand2 %in% x3) & !(Merged_lrpairs$receptor2 %in% x4),]
  Intracelluar_siganal_genes = Intracelluar_siganal_genes[Intracelluar_siganal_genes$rl %in% Merged_lrpairs$rl,]

  rl_remain = unique(Intracelluar_siganal_genes$rl)
  
  # 中介细胞的细胞内信号gene表达
  ct2_spots = unique(interaction_spots_pairs$mid)
  inner_exp = cbind(Intracelluar_siganal_genes, as.matrix(all_celltype_expression_list[[ct2]][Intracelluar_siganal_genes$network_genes, ct2_spots]))
  mid_spots_intracelluar_score = inner_exp[,-c(1:3)]
  mid_spots_intracelluar_score = mid_spots_intracelluar_score %>%
    group_by(rl) %>%
    summarise(across(c(colnames(mid_spots_intracelluar_score)[-1]), sum)) %>% 
    as.data.frame()
  mid_spots_intracelluar_score = mid_spots_intracelluar_score[,-1]
  rownames(mid_spots_intracelluar_score) = rl_remain
  colnames(mid_spots_intracelluar_score) = ct2_spots
  
  ### 根据lrpairs1 和 lrpairs2 的任意基因表达情况再进行过滤
  bg_spots = unique(bg_spots_pairs$mid)
  bg_inner_exp = cbind(Intracelluar_siganal_genes, as.matrix(all_celltype_expression_list[[ct2]][Intracelluar_siganal_genes$network_genes, bg_spots]))
  bg_mid_spots_intracelluar_score = bg_inner_exp[,-c(1:3)]
  bg_mid_spots_intracelluar_score = bg_mid_spots_intracelluar_score %>%
    group_by(rl) %>%
    summarise(across(c(colnames(bg_mid_spots_intracelluar_score)[-1]), sum)) %>% 
    as.data.frame()
  bg_mid_spots_intracelluar_score = bg_mid_spots_intracelluar_score[,-1]
  rownames(bg_mid_spots_intracelluar_score) = rl_remain
  colnames(bg_mid_spots_intracelluar_score) = bg_spots
  
  all_re = list()
  print(paste0('There are ',nrow(Merged_lrpairs), ' lr1-lr2 series interactions!'))
  for (n_M_lr in 1:nrow(Merged_lrpairs)) {
    if((n_M_lr %% 100) == 0){print(n_M_lr)}
    exp_ligand1 = all_celltype_expression_list[[ct1]][Merged_lrpairs$ligand1[n_M_lr],interaction_spots_pairs$Sender]
    exp_receptor1 = all_celltype_expression_list[[ct2]][Merged_lrpairs$receptor1[n_M_lr],interaction_spots_pairs$mid]
    exp_ligand2 = all_celltype_expression_list[[ct2]][Merged_lrpairs$ligand2[n_M_lr],interaction_spots_pairs$mid]
    exp_receptor2 = all_celltype_expression_list[[ct3]][Merged_lrpairs$receptor2[n_M_lr],interaction_spots_pairs$Receiver]
    RL_score = mid_spots_intracelluar_score[Merged_lrpairs$rl[n_M_lr],interaction_spots_pairs$mid] %>% t() %>% as.vector()
    final_score = exp_ligand1 * exp_receptor1 * exp_ligand2 * exp_receptor2 * RL_score * 10000
    final_score2 = which(final_score > 0)
    if (length(final_score2) > 0) {
      re = data.frame(Merged_lrpairs[rep(n_M_lr,length(final_score2)),], interaction_spots_pairs[final_score2,], 'score' = final_score[final_score2], 'p'=1)
    }else{
      next
    }
    exp_ligand1_null = all_celltype_expression_list[[ct1]][Merged_lrpairs$ligand1[n_M_lr], bg_spots_pairs$Sender]
    exp_receptor1_null = all_celltype_expression_list[[ct2]][Merged_lrpairs$receptor1[n_M_lr], bg_spots_pairs$mid]
    exp_ligand2_null = all_celltype_expression_list[[ct2]][Merged_lrpairs$ligand2[n_M_lr], bg_spots_pairs$mid]
    exp_receptor2_null = all_celltype_expression_list[[ct3]][Merged_lrpairs$receptor2[n_M_lr], bg_spots_pairs$Receiver]
    RL_score_null = bg_mid_spots_intracelluar_score[Merged_lrpairs$rl[n_M_lr],bg_spots_pairs$mid] %>% t() %>% as.vector()
    final_score_null = exp_ligand1_null * exp_receptor1_null * exp_ligand2_null * exp_receptor2_null * RL_score_null * 10000
    if (quantile(final_score_null,0.95) == 0) {
      re$p = 0
      all_re[[n_M_lr]] = as.data.table(re)
    }else{
      zero_pro = mean(final_score_null == 0)
      null_distribution_pos = final_score_null[final_score_null > 0]
      null_mean = mean(null_distribution_pos)
      null_var = var(null_distribution_pos)
      p = 1 - pgamma(re$p, shape = (null_mean * null_mean)/null_var, scale = null_var/null_mean)
      p = p * (1-zero_pro)
      re$p = p
      all_re[[n_M_lr]] = as.data.table(re)
    }
  }
  all_re = rbindlist(all_re, fill = T) %>% as.data.frame()
  all_re = all_re[all_re$p < 0.05,]
  
  Merged_lrpairs = Merged_lrpairs[Merged_lrpairs$rl %in% all_re$rl,]
  Intracelluar_siganal_genes = Intracelluar_siganal_genes[Intracelluar_siganal_genes$rl %in% all_re$rl,]
  interaction_spots_pairs = unique(all_re[,c(6:10)])
  
  results = list('All_score' = all_re, 'Merged_lrpairs' = Merged_lrpairs, 
                 'Intracelluar_siganal_genes' = Intracelluar_siganal_genes, 
                 'interaction_spots_pairs' = interaction_spots_pairs)
  return(results)
}


#' Calculate combined intracelluar Receptor-Liagnd siganal weight
#'
#' @param Input Normalized expression data matrix
#' @param Anno If mode is 'CellType', you need to provide the annotation of cells/spots
#' @param PPI_database A database of protein-protein interaction. 
#' Default is omnipath ppi database which contains a column called 'functional' which means stimulation, inhibition or unknown function by 1,-1,0.
#' If you want to provide yourself database, it should contain a 'functional' column and set the parameter 1.
#' @param Pro_loc_db A database of protein location.
#' Default is organized Uniprot database which contains a column called 'father_location' which contains three param 'Membrane', 'Intracellular' and 'Nuluar'.
#' If you want to provide yourself database, it should contain a 'father_location' column and set param.
#' @param lr_db A database of ligand-receptor.
#' Default is the lr_db database from 'Spatalk' which contains two columns called 'ligand' and 'receptor'
#' @param f_edge A function to calculate score of each edge, default is Hill function.
#' @param expression_cutoff A quantile to build a subgraph containing genes over cutoff. Default is 0.9.
#' @param cutoff1 A quantile to select the probility of R-L pairs over cutoff. Default is 0.5.
#' @param max_steps find_simple_path max steps in each R-L pairs
#' @param max_num number of paths to calculate path scores
#' @param alpha Path score and circuit score coefficients 
#' 
#' 
#' @return A dataframe of Receptor-Ligand signal weights of each cell or celltypes
#' @export 
#' @examples
#' utils::data("test_exp_data", package = "stMCCC")
#' utils::data("test_anno", package = "stMCCC")
#' results = Calculate_Combined_RL_weight(test_exp_data[,1,drop = F], 'SingleCell')
#' results = Calculate_Combined_RL_weight(test_exp_data, test_anno$Label, 'CellType')

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
    results <- do.call(rbind, pbapply::pblapply(cell_ids, cl = cl,function(cid) {
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





#' Calculate muti-cellular interactions in celltype resolution
#'
#' @param Input Normalized expression data matrix
#' @param Anno If mode is 'CellType', you need to provide the annotation of cells/spots
#' @param PPI_database A database of protein-protein interaction. 
#' Default is omnipath ppi database which contains a column called 'functional' which means stimulation, inhibition or unknown function by 1,-1,0.
#' If you want to provide yourself database, it should contain a 'functional' column and set the parameter 1.
#' @param Location_database A database of protein location.
#' Default is organized Uniprot database which contains a column called 'father_location' which contains three param 'Membrane', 'Intracellular' and 'Nuluar'.
#' If you want to provide yourself database, it should contain a 'father_location' column and set param.
#' @param lrpairs A database of ligand-receptor.
#' Default is the lrpairs database from 'Spatalk' which contains two columns called 'ligand' and 'receptor'
#' @param f_edge A function to calculate score of each edge, default is Hill function.
#' @param expression_cutoff A quantile to build a subgraph containing genes over cutoff. Default is 0.9.
#' @param lr_expression_cutoff A quantile to select the probility of R-L pairs over cutoff. Default is 0.5.
#' @param max_steps find_simple_path max steps in each R-L pairs
#' @param max_num number of paths to calculate path scores
#' @param alpha Path score and circuit score coefficients 
#' 
#' 
#' @return A dataframe of Receptor-Ligand signal weights of each cell or celltypes
#' @export 
#' @examples
#' utils::data("test_exp_data", package = "stMCCC")
#' utils::data("test_anno", package = "stMCCC")
#' results = Calculate_Combined_RL_weight(test_exp_data[,1,drop = F], 'SingleCell')
#' results = Calculate_Combined_RL_weight(test_exp_data, test_anno$Label, 'CellType')

calculate_multi_cci = function(Input, deconvolution_result, scRNA_mean_expression, relay_signal = NULL ,lrpairs = NULL) {
  cat("=== Start calculate muticellular cell cell interaction ===\n")
  start_time <- Sys.time()
  
  # 1. 加载配体-受体对数据
  if (is.null(lrpairs)) {
    cat("··· Load lrpairs database \n")
    utils::data("SpaTalk_lrpairs", package = "stMCCC")
    lrpairs <- SpaTalk_lrpairs[, c("ligand", "receptor")]
    cat("  Finished! There are ", nrow(lrpairs), " lr pairs. \n")
  } else {
    cat("  Finished! There are ", nrow(lrpairs), " lr pairs. \n")
  }
  
  # 2. 准备细胞类型信息
  all_celltype = colnames(deconvolution_result)
  cat("··· Detected these cell types: \n", paste(all_celltype, collapse = ", "), "\n")
  
  # 3. 表达重构
  cat("··· Start reconstruct expression...\n")
  spots_celltype_matrix = do.call("Construct_spots_celltype", list(deconvolution_result))
  deconvolution_reconsitution = do.call("Reconstruct_RCTD", list(deconvolution_result, spots_celltype_matrix))
  all_celltype_expression_list = do.call("Get_reconstruct_ST_expression_ct_list", 
                                         list(Input, scRNA_mean_expression, deconvolution_reconsitution))
  cat("  Finished! Dim: ", dim(deconvolution_reconsitution)[1], "spots ×", dim(deconvolution_reconsitution)[2], "cell types \n")
  
  # 4. 计算中继信号
  if (is.null(relay_signal)) {
    cat("··· Start calculate relay signal \n")
    relay_signal = do.call("Calculate_Combined_RL_weight", list(scRNA_mean_expression, 'SingleCell'))
    relay_signal = relay_signal[relay_signal$combined_score > 0, ]
    relay_signal$RL = paste0(relay_signal$receptor, '-', relay_signal$ligand)
    cat("   Finished! There are ", nrow(relay_signal), "RL pairs. \n")
  } else {
    relay_signal$RL = paste0(relay_signal$receptor, '-', relay_signal$ligand)
    cat("··· There are ", nrow(relay_signal), "RL pairs of input relay_signal. \n")
  }
  
  # 5. 构建relay矩阵
  cat("··· Construct relay matrix \n")
  relay_matrix = relay_signal[,-c(1,2,4,5,6)] %>% 
    pivot_wider(names_from = cell_id, values_from = combined_score, values_fill = 0) %>% 
    as.data.frame()
  rownames(relay_matrix) = relay_matrix$RL
  relay_matrix = relay_matrix[,-1] %>% as.matrix()
  cat("  Finished! Dim: ", dim(relay_matrix)[1], "RL ×", dim(relay_matrix)[2], "cell types \n")
  
  # 6. 计算空间基因表达均值
  cat("··· Start calculate spatial gene mean expression \n")
  spatial_mean_expression = matrix(0, nrow = nrow(Input), ncol = ncol(deconvolution_reconsitution), 
                                   dimnames = list(rownames(Input), colnames(deconvolution_reconsitution)))
  for (i in 1:ncol(spatial_mean_expression)) { 
    spatial_mean_expression[,i] = apply(Input,1,function(x){x * deconvolution_reconsitution[,i]}) %>% t() %>% rowMeans()
  }
  cat("  Finished! \n")
  
  # 7. 过滤表达的基因
  cat("··· Start filter genes \n")
  all_expressed_gene = rownames(spatial_mean_expression)[rowSums(spatial_mean_expression) != 0]
  lrpairs = lrpairs[lrpairs$ligand %in% all_expressed_gene & lrpairs$receptor %in% all_expressed_gene,]
  lrpairs$LR = paste0(lrpairs$ligand,'-', lrpairs$receptor)
  rl_pairs = relay_signal$RL[relay_signal$receptor %in% all_expressed_gene & relay_signal$ligand %in% all_expressed_gene] %>% unique()
  rl_pairs = rl_pairs[rl_pairs %in% relay_signal$RL]
  cat("  Finished! There are ", nrow(lrpairs), "LR pairs \n")
  
  # 8. 计算两两细胞互作
  cat("··· Calculate double cell type interaction score \n")
  double_celltype = expand.grid('ct1' = all_celltype, 'ct2' = all_celltype, stringsAsFactors = F)
  double_celltype$d_ctp = paste0(double_celltype$ct1,'-',double_celltype$ct2)
  double_celltype_interaction = matrix(0, nrow = nrow(lrpairs), ncol = nrow(double_celltype), 
                                       dimnames = list(lrpairs$LR, double_celltype$d_ctp))
  for (i in 1:ncol(double_celltype_interaction)) {
    double_celltype_interaction[,i] = spatial_mean_expression[lrpairs$ligand, double_celltype$ct1[i]] * 
      spatial_mean_expression[lrpairs$receptor, double_celltype$ct2[i]] * 10^4
  }
  cat("  Finished! Dim: ", dim(double_celltype_interaction)[1], "LR ×", dim(double_celltype_interaction)[2], "cell pairs\n")
  
  # 9. 准备三细胞组合
  cat("··· Prepair triple cell type pairs \n")
  tri_cell_pairs = expand.grid('ct1' = all_celltype, 'ct2' = all_celltype, 'ct3' = all_celltype, stringsAsFactors = F)
  tri_cell_pairs$ct1_ct2 = paste0(tri_cell_pairs$ct1,'-',tri_cell_pairs$ct2)
  tri_cell_pairs$ct2_ct3 = paste0(tri_cell_pairs$ct2,'-',tri_cell_pairs$ct3)
  tri_cell_pairs$ct1_ct2_ct3 = paste0(tri_cell_pairs$ct1,'-',tri_cell_pairs$ct2,'-',tri_cell_pairs$ct3)
  cat(" Finished! There are ", nrow(tri_cell_pairs), "triple cell type pairs. \n")
  
  # 10. 准备LR对组合
  cat("··· Prepair double lr pairs \n")
  double_lrpairs = expand.grid('lr1' = lrpairs$LR, 'lr2' = lrpairs$LR, stringsAsFactors = F)
  double_lrpairs$ligand1 = lrpairs$ligand[match(double_lrpairs$lr1, lrpairs$LR)]
  double_lrpairs$receptor1 = lrpairs$receptor[match(double_lrpairs$lr1, lrpairs$LR)]
  double_lrpairs$ligand2 = lrpairs$ligand[match(double_lrpairs$lr2, lrpairs$LR)]
  double_lrpairs$receptor2 = lrpairs$receptor[match(double_lrpairs$lr2, lrpairs$LR)]
  double_lrpairs$RL = paste0(double_lrpairs$receptor1,'-',double_lrpairs$ligand2)
  double_lrpairs = double_lrpairs[double_lrpairs$RL %in% rl_pairs,]
  double_lrpairs$LRLR = paste0(double_lrpairs$lr1,'-',double_lrpairs$lr2)
  cat("  Finished! There are ", nrow(double_lrpairs), "lr-lr pairs. \n")
  
  # 11. 计算三细胞互作
  cat("··· Calculate triple cell type interaction score \n")
  tri_cell_interaction = matrix(0, nrow = nrow(double_lrpairs), ncol = nrow(tri_cell_pairs))
  colnames(tri_cell_interaction) = tri_cell_pairs$ct1_ct2_ct3
  rownames(tri_cell_interaction) = double_lrpairs$LRLR
  
  progress_interval = max(round(nrow(tri_cell_pairs)/10), 100)  # 进度汇报间隔
  for (i in 1:nrow(tri_cell_pairs)) {
    if(i %% progress_interval == 0) {
      cat("  进度:", i, "/", nrow(tri_cell_pairs), 
          paste0("(", round(i/nrow(tri_cell_pairs)*100, 1), "%)"), "\n")
    }
    ct1_ct2 = tri_cell_pairs$ct1_ct2[i]
    ct2_ct3 = tri_cell_pairs$ct2_ct3[i]
    mid = tri_cell_pairs$ct2[i]
    lr1_score = double_celltype_interaction[double_lrpairs$lr1, ct1_ct2]
    lr2_score = double_celltype_interaction[double_lrpairs$lr2, ct2_ct3]
    mid_score = relay_matrix[double_lrpairs$RL, mid]
    tri_cell_interaction[,i] = lr1_score * mid_score * lr2_score
  }
  
  # 计算完成
  end_time <- Sys.time()
  cat("=== Finished! Time: ", round(as.numeric(end_time - start_time), 1), "seconds ===\n")
  
  return(list(
    'triple_cell_pairs' = tri_cell_pairs,
    'double_lrpais' = double_lrpairs,
    'interaction_scores' = tri_cell_interaction
  ))
}



#' Calculate multi cell cell single cell interactions
#'
#' @param Input Normalized expression data matrix
#' 
#' 
#' @return A dataframe of Receptor-Ligand signal weights of each cell or celltypes
#' @export 
#' @examples
#' utils::data("test_exp_data", package = "stMCCC")
#' utils::data("test_anno", package = "stMCCC")
#' results = Calculate_Combined_RL_weight(test_exp_data[,1,drop = F], 'SingleCell')
#' results = Calculate_Combined_RL_weight(test_exp_data, test_anno$Label, 'CellType')
calculate_single_ct_single_lrlr = function(Input, Coordinate, deconvolution_result, 
                                           ct1, ct2, ct3, L1, R1, L2, R2, 
                                           scRNA_mean_expression, relay_signal){
  ct_pair = paste0(ct1,'-',ct2,'-',ct3)
  lr_lr_pair = paste0(L1,'-',R1,'-',L2,'-',R2)
  
  relay_signal_score = relay_signal$combined_score[relay_signal$RL == paste0(R1,'-',L2) & relay_signal$cell_id == ct2]
  ### 计算三种细胞类型下的所有可能的空间中spots组合
  spots_celltype_matrix = Construct_spots_celltype(RCTD_result = deconvolution_result)
  
  Coordinate_matrix <- as.matrix(Coordinate[, c("x", "y")])
  knn_result <- get.knn(Coordinate_matrix, k = 6)
  neighbor_indices <- knn_result$nn.index
  neighbor_names <- apply(neighbor_indices, 1, function(idx) {rownames(Coordinate)[idx]  })
  knn_result_df <- data.frame(neighbors = t(neighbor_names))
  rownames(knn_result_df) = rownames(Coordinate)
  
  knn_result_df = knn_result_df[apply(spots_celltype_matrix,1,function(x){ct2 %in% x}),]
  knn_result_df = knn_result_df[apply(knn_result_df,1,function(x){ct1 %in% sapply(x, function(y){as.vector(spots_celltype_matrix[y,])})}),]
  knn_result_df = knn_result_df[apply(knn_result_df,1,function(x){ct3 %in% sapply(x, function(y){as.vector(spots_celltype_matrix[y,])})}),]
  knn_result_df = as.matrix(knn_result_df)
  if (nrow(knn_result_df) == 0) {
    return(NULL)
  }
  
  senders = matrix(0, nrow = nrow(knn_result_df), ncol = ncol(knn_result_df), dimnames = list(rownames(knn_result_df), colnames(knn_result_df)))
  receivers = matrix(0, nrow = nrow(knn_result_df), ncol = ncol(knn_result_df), dimnames = list(rownames(knn_result_df), colnames(knn_result_df)))
  for (i in 1:nrow(senders)) {
    for (k in 1:ncol(senders)) {
      sp = knn_result_df[i,k]
      if (ct1 %in% spots_celltype_matrix[sp,]) {
        senders[i,k] = 1
      }
      if (ct3 %in% spots_celltype_matrix[sp,]) {
        receivers[i,k] = 1
      }
    }
  }
  
  all_spot_pair = list()
  for (i in 1:nrow(knn_result_df)) {
    spot1 = knn_result_df[i,][senders[i,] == 1]
    spot2 = rownames(knn_result_df)[i]
    spot3 = knn_result_df[i,][receivers[i,] == 1]
    spot_pair = expand.grid(spot1, spot2, spot3, stringsAsFactors = F)
    all_spot_pair[[i]] = as.data.table(spot_pair)
  }
  all_spot_pair = rbindlist(all_spot_pair)
  colnames(all_spot_pair) = c('spot1','spot2','spot3')
  
  if (nrow(all_spot_pair) == 0) {
    spot_interaction = data.frame()
    spot_interaction$ct_pair = ct_pair
    spot_interaction$lr_lr_pair = lr_lr_pair
    return(spot_interaction)
  }
  
  all_spot_pair$x1 = Coordinate$x[match(all_spot_pair$spot1,rownames(Coordinate))]
  all_spot_pair$y1 = Coordinate$y[match(all_spot_pair$spot1,rownames(Coordinate))]
  all_spot_pair$x2 = Coordinate$x[match(all_spot_pair$spot2,rownames(Coordinate))]
  all_spot_pair$y2 = Coordinate$y[match(all_spot_pair$spot2,rownames(Coordinate))]
  all_spot_pair$x3 = Coordinate$x[match(all_spot_pair$spot3,rownames(Coordinate))]
  all_spot_pair$y3 = Coordinate$y[match(all_spot_pair$spot3,rownames(Coordinate))]
  
  deconvolution_reconsitution = Reconstruct_RCTD(deconvolution_result, spots_celltype_matrix)
  all_celltype_expression_list = Get_reconstruct_ST_expression_ct_list(Input, scRNA_mean_expression, deconvolution_reconsitution)
  
  all_spot_pair$exp_L1 = all_celltype_expression_list[[ct1]][L1, all_spot_pair$spot1]
  all_spot_pair$exp_R1 = all_celltype_expression_list[[ct2]][R1, all_spot_pair$spot2]
  all_spot_pair$exp_L2 = all_celltype_expression_list[[ct2]][L2, all_spot_pair$spot2]
  all_spot_pair$exp_R2 = all_celltype_expression_list[[ct3]][R2, all_spot_pair$spot3]
  all_spot_pair$interaction_score = all_spot_pair$exp_L1 * all_spot_pair$exp_R1 * 10^4 * relay_signal_score *
    all_spot_pair$exp_L2 * all_spot_pair$exp_R2 * 10^4
  all_spot_pair$ct_pair = ct_pair
  all_spot_pair$lr_lr_pair = lr_lr_pair
  spot_interaction = all_spot_pair[all_spot_pair$interaction_score > 0,,drop = F]
  return(spot_interaction)
}




#' Calculate multi cell cell interactions of all lr-lr pairs
#'
#' @param Input Normalized expression data matrix
#' @param Coordinate Normalized expression data matrix
#' @param st_type Normalized expression data matrix
#' 
#' 
#' @return A dataframe of Receptor-Ligand signal weights of each cell or celltypes
#' @export 
#' @examples
#' utils::data("test_exp_data", package = "stMCCC")
#' utils::data("test_anno", package = "stMCCC")
#' results = Calculate_Combined_RL_weight(test_exp_data[,1,drop = F], 'SingleCell')
#' results = Calculate_Combined_RL_weight(test_exp_data, test_anno$Label, 'CellType')

calculate_all_lrlr = function(Input, Coordinate, st_type = c('low', 'high'), deconvolution_result = NULL, scRNA_mean_expression = NULL,
                              ct1, ct2, ct3, relay_signal, lrpairs = NULL, n_permutations = 1000) {
  if(st_type == 'low' & (is.null(deconvolution_result) | is.null(scRNA_mean_expression))){
    cat('Low resolution spatial transcriptome data requires deconvolution result and scRNA mean expression!')
    break
  }
  
  ct_pair = paste0(ct1,'-',ct2,'-',ct3)
  
  # 根据提供的lrpairs预先构建lr-lr组合
  expressed_genes = rownames(Input)[rowSums(Input) > 0]
  if(is.null(lrpairs)){
    utils::data("SpaTalk_lrpairs", package = "stMCCC")
    lrpairs <- SpaTalk_lrpairs[, c("ligand", "receptor")]
    lrpairs$LR = paste0(lrpairs$ligand,'-',lrpairs$receptor)
    lrpairs = lrpairs[lrpairs$ligand %in% expressed_genes & lrpairs$receptor %in% expressed_genes,]
  }
  
  # double_lrpairs 构建
  lrpairs_idx <- seq_len(nrow(lrpairs))
  idx_grid <- expand.grid(lr1 = lrpairs_idx, lr2 = lrpairs_idx, stringsAsFactors = FALSE)
  double_lrpairs <- data.frame(
    lr1 = lrpairs$LR[idx_grid$lr1],
    lr2 = lrpairs$LR[idx_grid$lr2],
    ligand1 = lrpairs$ligand[idx_grid$lr1],
    receptor1 = lrpairs$receptor[idx_grid$lr1],
    ligand2 = lrpairs$ligand[idx_grid$lr2],
    receptor2 = lrpairs$receptor[idx_grid$lr2]
  )
  double_lrpairs$RL = paste0(double_lrpairs$receptor1,'-',double_lrpairs$ligand2)
  double_lrpairs = double_lrpairs[double_lrpairs$RL %in% relay_signal[relay_signal$cell_id == ct2,]$RL,]
  double_lrpairs$LRLR = paste0(double_lrpairs$lr1,'-',double_lrpairs$lr2)
  
  
  
  if (st_type == 'low') {
    spots_celltype_matrix = Construct_spots_celltype(RCTD_result = deconvolution_result)
    
    # 0. 计算每一个细胞类型在spots中的表达（仅低分辨率需要，同时也是低分辨率数据置换检验的必要步骤）
    deconvolution_reconsitution = Reconstruct_RCTD(deconvolution_result, spots_celltype_matrix)
    all_celltype_expression_list = Get_reconstruct_ST_expression_ct_list(Input, scRNA_mean_expression, deconvolution_reconsitution)
    cat('=== Finished calculate the spatial expression of each cell type! \n')
    
    # 1. 计算所有spots的临近关系
    Coordinate_matrix <- as.matrix(Coordinate[, c("x", "y")])
    knn_result <- get.knn(Coordinate_matrix, k = 6)
    neighbor_indices <- knn_result$nn.index
    neighbor_names <- apply(neighbor_indices, 1, function(idx) {rownames(Coordinate)[idx]  })
    knn_result_df <- data.frame(neighbors = t(neighbor_names))
    rownames(knn_result_df) = rownames(Coordinate)
    
    # 2. 选出符合ct1-ct2-ct3要求的spots组合
    knn_result_df = knn_result_df[apply(spots_celltype_matrix,1,function(x){ct2 %in% x}),]
    knn_result_df = knn_result_df[apply(knn_result_df,1,function(x){ct1 %in% sapply(x, function(y){as.vector(spots_celltype_matrix[y,])})}),]
    knn_result_df = knn_result_df[apply(knn_result_df,1,function(x){ct3 %in% sapply(x, function(y){as.vector(spots_celltype_matrix[y,])})}),]
    knn_result_df = as.matrix(knn_result_df)
    if (nrow(knn_result_df) == 0) {
      return(NULL)
    }
    senders = matrix(0, nrow = nrow(knn_result_df), ncol = ncol(knn_result_df), dimnames = list(rownames(knn_result_df), colnames(knn_result_df)))
    receivers = matrix(0, nrow = nrow(knn_result_df), ncol = ncol(knn_result_df), dimnames = list(rownames(knn_result_df), colnames(knn_result_df)))
    for (i in 1:nrow(senders)) {
      for (k in 1:ncol(senders)) {
        sp = knn_result_df[i,k]
        if (ct1 %in% spots_celltype_matrix[sp,]) {
          senders[i,k] = 1
        }
        if (ct3 %in% spots_celltype_matrix[sp,]) {
          receivers[i,k] = 1
        }
      }
    }
    all_spot_pair = list()
    for (i in 1:nrow(knn_result_df)) {
      spot1 = knn_result_df[i,][senders[i,] == 1]
      spot2 = rownames(knn_result_df)[i]
      spot3 = knn_result_df[i,][receivers[i,] == 1]
      spot_pair = expand.grid(spot1, spot2, spot3, stringsAsFactors = F)
      all_spot_pair[[i]] = as.data.table(spot_pair)
    }
    all_spot_pair = rbindlist(all_spot_pair)
    colnames(all_spot_pair) = c('spot1','spot2','spot3')
    if (nrow(all_spot_pair) == 0) {
      spot_interaction = data.frame()
      spot_interaction$ct_pair = ct_pair
      spot_interaction$lr_lr_pair = 'null'
      return(spot_interaction)
    }
    
    # 3. 记录每个spot的坐标值
    all_spot_pair$x1 = Coordinate$x[match(all_spot_pair$spot1,rownames(Coordinate))]
    all_spot_pair$y1 = Coordinate$y[match(all_spot_pair$spot1,rownames(Coordinate))]
    all_spot_pair$x2 = Coordinate$x[match(all_spot_pair$spot2,rownames(Coordinate))]
    all_spot_pair$y2 = Coordinate$y[match(all_spot_pair$spot2,rownames(Coordinate))]
    all_spot_pair$x3 = Coordinate$x[match(all_spot_pair$spot3,rownames(Coordinate))]
    all_spot_pair$y3 = Coordinate$y[match(all_spot_pair$spot3,rownames(Coordinate))]
    
    cat('=== Start building permutation test 0 distribution data \n')
    # 4. 构建置换检验的数据
    spots_ct1 = rownames(spots_celltype_matrix)[apply(spots_celltype_matrix,1,function(x){ct1 %in% x})]
    spots_ct2 = rownames(spots_celltype_matrix)[apply(spots_celltype_matrix,1,function(x){ct2 %in% x})]
    spots_ct3 = rownames(spots_celltype_matrix)[apply(spots_celltype_matrix,1,function(x){ct3 %in% x})]
    
    spot_pair_permutation = data.frame('spot1' = sample(spots_ct1, n_permutations, replace = T),
                                       'spot2' = sample(spots_ct2, n_permutations, replace = T),
                                       'spot3' = sample(spots_ct3, n_permutations, replace = T))
    
    # 5. 计算所有lr-lr组合分数
    cat('=== Start calculate all lr-lr score \n')
    progress_interval = max(round(nrow(double_lrpairs)/10), 100)  # 进度汇报间隔
    re = list()
    for (g in 1:nrow(double_lrpairs)) {
      if(g %% progress_interval == 0) {
        cat("  Process:", g, "/", nrow(double_lrpairs), paste0("(", round(g/nrow(double_lrpairs)*100, 1), "%)"), "\n")
      }
      double_lrpairs_select = double_lrpairs[g,,drop = F]
      L1 = double_lrpairs_select$ligand1
      R1 = double_lrpairs_select$receptor1
      L2 = double_lrpairs_select$ligand2
      R2 = double_lrpairs_select$receptor2
      lr_lr_pair = double_lrpairs_select$LRLR
      exp_L1 = all_celltype_expression_list[[ct1]][L1, all_spot_pair$spot1]
      exp_R1 = all_celltype_expression_list[[ct2]][R1, all_spot_pair$spot2]
      exp_L2 = all_celltype_expression_list[[ct2]][L2, all_spot_pair$spot2]
      exp_R2 = all_celltype_expression_list[[ct3]][R2, all_spot_pair$spot3]
      if(sum(exp_L1) == 0 | sum(exp_R1) == 0 | sum(exp_L2) == 0 | sum(exp_R2) == 0){
        re[[g]] = NULL
        next
      }
      all_spot_pair_lrlr = all_spot_pair
      all_spot_pair_lrlr$exp_L1 = exp_L1
      all_spot_pair_lrlr$exp_R1 = exp_R1
      all_spot_pair_lrlr$exp_L2 = exp_L2
      all_spot_pair_lrlr$exp_R2 = exp_R2
      
      R1_L2_st_exp_mean = mean(all_celltype_expression_list[[ct2]][R1,]) * mean(all_celltype_expression_list[[ct2]][L2,])
      relay_signal_score = relay_signal$combined_score[relay_signal$RL == paste0(R1,'-',L2) & relay_signal$cell_id == ct2]
      if(R1_L2_st_exp_mean == 0){
        re[[g]] = NULL
        next
      }
      all_spot_pair_lrlr$relay_score = all_spot_pair_lrlr$exp_R1 * all_spot_pair_lrlr$exp_L2/R1_L2_st_exp_mean * relay_signal_score
      all_spot_pair_lrlr$interaction_score = all_spot_pair_lrlr$exp_L1 * all_spot_pair_lrlr$exp_R1 * all_spot_pair_lrlr$relay_score *
        all_spot_pair_lrlr$exp_L2 * all_spot_pair_lrlr$exp_R2
      all_spot_pair_lrlr$ct_pair = ct_pair
      all_spot_pair_lrlr$lr_lr_pair = lr_lr_pair
      spot_interaction = all_spot_pair_lrlr[all_spot_pair_lrlr$interaction_score > 0,,drop = F]
      
      # 置换检验
      permuted_scores <- numeric(n_permutations)
      # 重新计算表达值
      permuted_exp_L1 = all_celltype_expression_list[[ct1]][L1, spot_pair_permutation$spot1]
      permuted_exp_R1 = all_celltype_expression_list[[ct2]][R1, spot_pair_permutation$spot2]
      permuted_exp_L2 = all_celltype_expression_list[[ct2]][L2, spot_pair_permutation$spot2]
      permuted_exp_R2 = all_celltype_expression_list[[ct3]][R2, spot_pair_permutation$spot3]
      
      if(sum(permuted_exp_L1) == 0 | sum(permuted_exp_R1) == 0 | sum(permuted_exp_L2) == 0 | sum(permuted_exp_R2) == 0){
        permuted_scores <- numeric(n_permutations)
      }
      permuted_R1_L2_st_exp_mean = mean(all_celltype_expression_list[[ct2]][R1,]) * mean(all_celltype_expression_list[[ct2]][L2,])
      if(permuted_R1_L2_st_exp_mean == 0){
        permuted_scores <- numeric(n_permutations)
        next
      }
      permuted_relay_score = permuted_exp_R1 * permuted_exp_L2 / permuted_R1_L2_st_exp_mean * relay_signal_score
      permuted_interaction_score = permuted_exp_L1 * permuted_exp_R1 * permuted_relay_score *permuted_exp_L2 * permuted_exp_R2
      permuted_scores <- permuted_interaction_score
      
      # 计算 p 值
      ecdf_permuted_scores <- ecdf(permuted_scores) 
      quantile_position <- ecdf_permuted_scores(spot_interaction$interaction_score)
      p_value <- 1-quantile_position
      spot_interaction$p_value <- p_value
      re[[g]] = as.data.table(spot_interaction)
    }
    re = rbindlist(re)
    return(re)
    
  }else if(st_type == 'high') {
    
  }
}

#' Calculate Triple Cell-Cell Interaction (CCI) with Customized Parameters
#'
#' This function calculates the triple cell-cell interaction scores based on the provided input data.
#' It supports both parallel and sequential processing, and can handle different types of spatial transcriptomics data.
#'
#' @param Input A matrix or data frame containing gene expression data.
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