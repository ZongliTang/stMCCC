#' Get interaction distance cutoff
#'
#' @param coords Coordinate of the Visium data 
#' @param n_hop Selected the spots around a single spot. Default n_hop = 1. n_hop can be 1,2,3.
#' 
#' @return A distance length for the interaction area and print the mean of interacted spots of a single spot
#' @export 
#' @examples
#' RCTD_result = load_ccRCC_RCTD_result() 
Get_interaction_cutoff = function(coords, n_hop = 1){
  all_spots = rownames(coords)
  D_matrix_df = dist(coords) %>% as.matrix()
  D_matrix_df = reshape2::melt(D_matrix_df) %>% as.data.frame()
  colnames(D_matrix_df) = c('Sender', 'Receiver','distance')
  D_matrix_df$Sender = as.character(D_matrix_df$Sender)
  D_matrix_df$Receiver = as.character(D_matrix_df$Receiver)
  d1 = max(coords$x)/sqrt(length(all_spots))
  d2 = max(coords$y)/sqrt(length(all_spots))
  d = round(max(d1,d2))
  if (n_hop == 1) {
    cut = d
    print(sum(0 < D_matrix_df$distance & D_matrix_df$distance < cut)/length(all_spots))
  }else if(n_hop == 2){
    cut = d*2
    print(sum(0 < D_matrix_df$distance & D_matrix_df$distance < cut)/length(all_spots))
  }else{
    cut = d*3
    print(sum(0 < D_matrix_df$distance & D_matrix_df$distance < cut)/length(all_spots))
  }
  return(cut)
}