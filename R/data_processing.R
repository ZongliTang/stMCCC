#' Load ccRCC Visium test dataset 
#'
#' @param pathway Pathway of the dataset. 
#' In this pathway, there is a filtered_feature_bc_matrix.h5 dgCmatrix file, a TLS annotation file
#'  and a spatial fold containing tissue_positions_list.csv and tissue_lowres_image.png
#' @return A seurat object of the test dataset
#' @export 
#' @examples
#' st_seu_obj = load_ccRCC_Visium_data()  # Return a seurat Object
load_ccRCC_Visium_data = function(pathway = '/Users/tangzongli_macmini/Desktop/2025 work/My_packages/stMCCC/test/ccRCC'){
  exp_mat = Read10X_h5(file.path(pathway,"filtered_feature_bc_matrix.h5"))
  coords = read.csv(file.path(pathway,"spatial/tissue_positions_list.csv"),header = F)
  coords = coords[coords$V2 == "1",c(1,5,6)]
  rownames(coords) = coords$V1
  coords$V1 = NULL
  colnames(coords) = c("x","y")
  spatial_img <- Read10X_Image(image.dir = file.path(pathway,"spatial"),image.name = "tissue_lowres_image.png",filter.matrix = F)
  
  st_seu_obj <- CreateSeuratObject(counts = exp_mat, assay = "Spatial", project = 'ccRCC', min.cells = 3, min.features = 100)
  coords <- coords[Cells(st_seu_obj), ]
  spatial_img <- spatial_img[Cells(st_seu_obj)]
  st_seu_obj[["slice1"]] <- spatial_img
  DefaultAssay(st_seu_obj) <- "Spatial"
  TLS_anno = read.csv(file.path(pathway,'TLS_annotation.csv'))
  TLS_anno = TLS_anno[TLS_anno$Barcode %in% rownames(coords),]
  st_seu_obj$TLS = TLS_anno$TLS_2_cat
  return(st_seu_obj)
}


#' Load ccRCC reference scRNA-seq dataset 
#'
#' @param pathway Pathway of the dataset. 
#' In this pathway, there is a annotation file and a 'GSE159115_RAW' fold containing
#' filtered_gene_bc_matrices_h5 files
#' @return A scRNA-seq seurat object of ccRCC dataset
#' @export 
#' @examples
#' sc_seu_obj = load_ccRCC_scRNA_seq_data()
load_ccRCC_scRNA_seq_data = function(pathway = '/Users/tangzongli_macmini/Desktop/2025 work/My_packages/stMCCC/test/ccRCC_scRNA'){
  sc_anno = read.csv(file.path(pathway,'GSE159115_ccRCC_anno.csv'), header = T)
  T1 = Read10X_h5(file.path(pathway,'GSE159115_RAW/GSM4819725_SI_18854_filtered_gene_bc_matrices_h5.h5'))
  T2 = Read10X_h5(file.path(pathway,'GSE159115_RAW/GSM4819727_SI_18855_filtered_gene_bc_matrices_h5.h5'))
  T3 = Read10X_h5(file.path(pathway,'GSE159115_RAW/GSM4819729_SI_19703_filtered_gene_bc_matrices_h5.h5'))
  T4 = Read10X_h5(file.path(pathway,'GSE159115_RAW/GSM4819734_SI_22368_filtered_gene_bc_matrices_h5.h5'))
  T5 = Read10X_h5(file.path(pathway,'GSE159115_RAW/GSM4819736_SI_22604_filtered_gene_bc_matrices_h5.h5'))
  T6 = Read10X_h5(file.path(pathway,'GSE159115_RAW/GSM4819737_SI_23459_filtered_gene_bc_matrices_h5.h5'))
  T7 = Read10X_h5(file.path(pathway,'GSE159115_RAW/GSM4819738_SI_23843_filtered_gene_bc_matrices_h5.h5'))
  colnames(T1) = paste0("SI_18854_", colnames(T1))
  colnames(T2) = paste0("SI_18855_", colnames(T2))
  colnames(T3) = paste0("SI_19703_", colnames(T3))
  colnames(T4) = paste0("SI_22368_", colnames(T4))
  colnames(T5) = paste0("SI_22604_", colnames(T5))
  colnames(T6) = paste0("SI_23459_", colnames(T6))
  colnames(T7) = paste0("SI_23843_", colnames(T7))
  sc_data_all = cbind(T1,T2,T3,T4,T5,T6,T7)
  sc_data_all = sc_data_all[,colnames(sc_data_all) %in% sc_anno$cell]
  sc_data_all = sc_data_all[,sc_anno$cell]
  sc_seu_obj = CreateSeuratObject(sc_data_all, assay = "RNA", meta.data = sc_anno,min.cells = 3, min.features = 200)
  sc_seu_obj = subset(sc_seu_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & pct_MT < 0.1)
  return(sc_seu_obj)
}

#' Load RCTD deconvelution result 
#'
#' @param pathway Pathway of the result. 
#' 
#' @return A matrix with n_spots * c_cell_types
#' @export 
#' @examples
#' RCTD_result = load_ccRCC_RCTD_result() 
load_ccRCC_RCTD_result = function(pathway = '/Users/tangzongli_macmini/Desktop/2025 work/My_packages/stMCCC/test/ccRCC/RCTD/RCTD_matrix.csv'){
  RCTD_re = read.csv(pathway)
  rownames(RCTD_re) = RCTD_re$X
  RCTD_re$X = NULL
  RCTD_re = as.matrix(RCTD_re)
  RCTD_re = RCTD_re[all_spots,]
  return(RCTD_re)
}


#' Get interaction double spots pairs
#'
#' @param coords Coordinate of the Visium data 
#' @param n_hop Selected the spots around a single spot. Default n_hop = 1. n_hop can be 1,2,3.
#' 
#' @return A list
#' @export 
#' @examples
#' double_spots_pairs = Get_double_spots_pairs(coords)
Get_double_spots_pairs = function(coords, n_hop = 1){
  D_matrix_df = stats::dist(coords) %>% as.matrix()
  D_matrix_df = reshape2::melt(D_matrix_df) %>% as.data.frame()
  colnames(D_matrix_df) = c('Sender', 'Receiver','distance')
  D_matrix_df$Sender = as.character(D_matrix_df$Sender)
  D_matrix_df$Receiver = as.character(D_matrix_df$Receiver)
  distance_cutoff = Get_interaction_cutoff(coords,n_hop)
  D_matrix_df_interact = D_matrix_df[D_matrix_df$distance < distance_cutoff,]
  D_matrix_df_no_interact = D_matrix_df[D_matrix_df$distance >= distance_cutoff,]
  double_spots_pairs = list('interaction' = D_matrix_df_interact, 'no_interaction'= D_matrix_df_no_interact)
  return(double_spots_pairs)
}




