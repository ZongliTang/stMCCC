#' Plot hierarchical interaction of all celltypes
#'
#' @param celltype_multi_cci Result from the function 'calculate_multi_cci'. We need the interaction scores and triple cell pairs.
#' @param Mode 'Score' or 'Num'. Represent the mean score or the number of selected lr-lr pairs in all cell type pairs
#' 
#' @return A ggplot object
#' @export 
#' @examples
#' plot = plot_hierarchical_all_celltypes_interaction(celltype_multi_cci, Mode = 'Score')
#' plot = plot_hierarchical_all_celltypes_interaction(celltype_multi_cci, Mode = 'Num')

plot_hierarchical_all_celltypes_interaction <- function(celltype_multi_cci, Mode = NULL){
  interaction_scores = celltype_multi_cci$interaction_scores
  triple_cell_pairs = celltype_multi_cci$triple_cell_pairs
  if (Mode == 'Score') {
    if (nrow(plot1) > 1) {
      plot1 = colMeans(interaction_scores)
    }
    plot1 = data.frame(triple_cell_pairs[,c(1:3)], 'score' = plot1)
    plot1 = plot1[plot1$score > 0,]
    all_cell_types <- unique(c(plot1$ct1, plot1$ct2, plot1$ct3)) %>% sort()
    nodes <- data.frame(cell_type = all_cell_types,y = seq_along(all_cell_types))
    nodes_full <- expand.grid(cell_type = all_cell_types, x = 1:3,stringsAsFactors = FALSE) %>% 
      left_join(nodes, by = "cell_type")
    edges1 <- plot1 %>% select(from = ct1, to = ct2, score) %>% mutate(x_start = 1, x_end = 2) %>% 
      mutate(log_score = log10(score + 1e-10)) %>% 
      left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
      left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
    edges2 <- plot1 %>% select(from = ct2, to = ct3, score) %>% mutate(x_start = 2, x_end = 3) %>% 
      mutate(log_score = log10(score + 1e-10)) %>% 
      left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
      left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
    plot = ggplot() +
      geom_point(data = nodes_full, aes(x = x, y = y), size = 5, fill = "white", shape = 21, stroke = 1.2)+ # 绘制三列节点
      geom_segment(data = edges1, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                   linewidth = 1.2,
                   lineend = "round")+ # 绘制第一列和第二列的连边和分数
      geom_segment(data = edges2, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                   linewidth = 1.2,
                   lineend = "round")+ # 绘制第二列和第三列的连边和分数
      geom_text_repel(
        data = nodes_full,
        aes(x = x, y = y, label = cell_type),
        direction = "y",  # 仅垂直方向调整
        nudge_y = 0.3,    # 垂直偏移量
        segment.color = NA,
        size = 3.5,
        fontface = "bold"
      ) +
      scale_color_gradientn(
        colors = c("darkblue","#F0F0F0", "#CB3425"),
        name = "Log10(Score)",
        breaks = scales::pretty_breaks(n = 5)
      ) +
      scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
      scale_x_continuous(
        breaks = 1:3,
        labels = c("CT1 (Sender)", "CT2 (Relay)", "CT3 (Receiver)"),
        limits = c(0.8, 3.2)
      ) +
      labs(title = "Cell-Cell Interaction Network across Three Hierarchical Levels", x = "", y = "") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  }else if(Mode == 'Num'){
    plot1 = colSums(interaction_scores > 0)
    plot1 = data.frame(triple_cell_pairs[,c(1:3)], 'score' = plot1)
    plot1 = plot1[plot1$score > 0,]
    all_cell_types <- unique(c(plot1$ct1, plot1$ct2, plot1$ct3)) %>% sort()
    nodes <- data.frame(cell_type = all_cell_types,y = seq_along(all_cell_types))
    nodes_full <- expand.grid(cell_type = all_cell_types, x = 1:3,stringsAsFactors = FALSE) %>% 
      left_join(nodes, by = "cell_type")
    edges1 <- plot1 %>% select(from = ct1, to = ct2, score) %>% mutate(x_start = 1, x_end = 2) %>% 
      mutate(log_score = log10(score + 1e-10)) %>% 
      left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
      left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
    edges2 <- plot1 %>% select(from = ct2, to = ct3, score) %>% mutate(x_start = 2, x_end = 3) %>% 
      mutate(log_score = log10(score + 1e-10)) %>% 
      left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
      left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
    
    plot = ggplot() +
      geom_point(data = nodes_full, aes(x = x, y = y), size = 5, fill = "white", shape = 21, stroke = 1.2)+ # 绘制三列节点
      geom_segment(data = edges1, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                   linewidth = 1.2,
                   lineend = "round")+ # 绘制第一列和第二列的连边和分数
      geom_segment(data = edges2, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                   linewidth = 1.2,
                   lineend = "round")+ # 绘制第二列和第三列的连边和分数
      geom_text_repel(
        data = nodes_full,
        aes(x = x, y = y, label = cell_type),
        direction = "y",  # 仅垂直方向调整
        nudge_y = 0.3,    # 垂直偏移量
        segment.color = NA,
        size = 3.5,
        fontface = "bold"
      ) +
      scale_color_gradientn(
        colors = c("darkblue","#F0F0F0", "#CB3425"),
        name = "Log10(Num)",
        breaks = scales::pretty_breaks(n = 5)
      ) +
      scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
      scale_x_continuous(
        breaks = 1:3,
        labels = c("CT1 (Sender)", "CT2 (Relay)", "CT3 (Receiver)"),
        limits = c(0.8, 3.2)
      ) +
      labs(title = "Cell-Cell Interaction Network across Three Hierarchical Levels", x = "", y = "") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  }
  return(plot)
}

#' Plot hierarchical interaction of selected celltypes
#'
#' @param ct1s A list of selected sender cell types.
#' @param ct2s A list of selected relay cell types.
#' @param ct3s A list of selected receiver cell types.
#' @param celltype_multi_cci Result from the function 'calculate_multi_cci'. We need the interaction scores and triple cell pairs.
#' @param Mode 'Score' or 'Num'. Represent the mean score or the number of selected lr-lr pairs in all cell type pairs
#' 
#' @return A ggplot object
#' @export 
#' @examples
#' ct1s = c('Tumor', 'Bcell','Tcell')
#' ct2s = c('Macro','Tcell')
#' ct3s = c('Endo_PLVAP','Plasma')
#' plot_hierarchical_all_celltypes_interaction(celltype_multi_cci, Mode = 'Score')
#' plot_hierarchical_all_celltypes_interaction(celltype_multi_cci, Mode = 'Num')

plot_hierarchical_selected_celltypes_interaction <- function(ct1s, ct2s, ct3s, celltype_multi_cci, Mode = NULL){
  interaction_scores = celltype_multi_cci$interaction_scores
  triple_cell_pairs = celltype_multi_cci$triple_cell_pairs
  triple_cell_pairs = triple_cell_pairs[triple_cell_pairs$ct1 %in% ct1s & 
                                          triple_cell_pairs$ct2 %in% ct2s & 
                                          triple_cell_pairs$ct3 %in% ct3s,]
  interaction_scores = interaction_scores[,triple_cell_pairs$ct1_ct2_ct3,drop = F]
  
  
  if (Mode == 'Score') {
    if(nrow(interaction_scores) > 1){
      plot1 = colMeans(interaction_scores)
    }else{
      plot1 = interaction_scores[,1]
    }
    plot1 = data.frame(triple_cell_pairs[,c(1:3)], 'score' = plot1)
    plot1 = plot1[plot1$score > 0,]
    all_cell_types <- unique(c(plot1$ct1, plot1$ct2, plot1$ct3)) %>% sort()
    nodes <- data.frame(cell_type = all_cell_types,y = seq_along(all_cell_types))
    nodes_full <- expand.grid(cell_type = all_cell_types, x = 1:3,stringsAsFactors = FALSE) %>% 
      left_join(nodes, by = "cell_type")
    nodes_full = nodes_full[(nodes_full$cell_type %in% ct1s & nodes_full$x == 1) |
                              (nodes_full$cell_type %in% ct2s & nodes_full$x == 2) |
                              (nodes_full$cell_type %in% ct3s & nodes_full$x == 3),]
    edges1 <- plot1 %>% select(from = ct1, to = ct2, score) %>% mutate(x_start = 1, x_end = 2) %>% 
      mutate(log_score = log10(score + 1e-10)) %>% 
      left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
      left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
    edges2 <- plot1 %>% select(from = ct2, to = ct3, score) %>% mutate(x_start = 2, x_end = 3) %>% 
      mutate(log_score = log10(score + 1e-10)) %>% 
      left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
      left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
    plot = ggplot() +
      geom_point(data = nodes_full, aes(x = x, y = y), size = 5, fill = "white", shape = 21, stroke = 1.2)+ # 绘制三列节点
      geom_segment(data = edges1, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                   linewidth = 1.2,
                   lineend = "round")+ # 绘制第一列和第二列的连边和分数
      geom_segment(data = edges2, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                   linewidth = 1.2,
                   lineend = "round")+ # 绘制第二列和第三列的连边和分数
      geom_text_repel(
        data = nodes_full,
        aes(x = x, y = y, label = cell_type),
        direction = "y",  # 仅垂直方向调整
        nudge_y = 0.3,    # 垂直偏移量
        segment.color = NA,
        size = 3.5,
        fontface = "bold"
      ) +
      scale_color_gradientn(
        colors = c("darkblue","#F0F0F0", "#CB3425"),
        name = "Log10(Score)",
        breaks = scales::pretty_breaks(n = 5)
      ) +
      scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
      scale_x_continuous(
        breaks = 1:3,
        labels = c("CT1 (Sender)", "CT2 (Relay)", "CT3 (Receiver)"),
        limits = c(0.8, 3.2)
      ) +
      labs(title = "Cell-Cell Interaction Network across Three Hierarchical Levels", x = "", y = "") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  }else if(Mode == 'Num'){
    plot1 = colSums(interaction_scores > 0)
    plot1 = data.frame(triple_cell_pairs[,c(1:3)], 'score' = plot1)
    plot1 = plot1[plot1$score > 0,]
    all_cell_types <- unique(c(plot1$ct1, plot1$ct2, plot1$ct3)) %>% sort()
    nodes <- data.frame(cell_type = all_cell_types,y = seq_along(all_cell_types))
    nodes_full <- expand.grid(cell_type = all_cell_types, x = 1:3,stringsAsFactors = FALSE) %>% 
      left_join(nodes, by = "cell_type")
    nodes_full = nodes_full[(nodes_full$cell_type %in% ct1s & nodes_full$x == 1) |
                              (nodes_full$cell_type %in% ct2s & nodes_full$x == 2) |
                              (nodes_full$cell_type %in% ct3s & nodes_full$x == 3),]
    edges1 <- plot1 %>% select(from = ct1, to = ct2, score) %>% mutate(x_start = 1, x_end = 2) %>% 
      mutate(log_score = log10(score + 1e-10)) %>% 
      left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
      left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
    edges2 <- plot1 %>% select(from = ct2, to = ct3, score) %>% mutate(x_start = 2, x_end = 3) %>% 
      mutate(log_score = log10(score + 1e-10)) %>% 
      left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
      left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
    plot = ggplot() +
      geom_point(data = nodes_full, aes(x = x, y = y), size = 5, fill = "white", shape = 21, stroke = 1.2)+ # 绘制三列节点
      geom_segment(data = edges1, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                   linewidth = 1.2,
                   lineend = "round")+ # 绘制第一列和第二列的连边和分数
      geom_segment(data = edges2, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                   linewidth = 1.2,
                   lineend = "round")+ # 绘制第二列和第三列的连边和分数
      geom_text_repel(
        data = nodes_full,
        aes(x = x, y = y, label = cell_type),
        direction = "y",  # 仅垂直方向调整
        nudge_y = 0.3,    # 垂直偏移量
        segment.color = NA,
        size = 3.5,
        fontface = "bold"
      ) +
      scale_color_gradientn(
        colors = c("darkblue","#F0F0F0", "#CB3425"),
        name = "Log10(Num)",
        breaks = scales::pretty_breaks(n = 5)
      ) +
      scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
      scale_x_continuous(
        breaks = 1:3,
        labels = c("CT1 (Sender)", "CT2 (Relay)", "CT3 (Receiver)"),
        limits = c(0.8, 3.2)
      ) +
      labs(title = "Cell-Cell Interaction Network across Three Hierarchical Levels", x = "", y = "") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  }
  return(plot)
}








#' Plot hierarchical selected interaction of all celltypes
#'
#' @param celltype_multi_cci Result from the function 'calculate_multi_cci'. We need the interaction scores, double_lrpais and triple cell pairs.
#' @param plot_row The charactor name of lr-lr or the number of lr-lr in double_lrpais.
#' 
#' @return A ggplot object
#' @export 
#' @examples 
#' plot_hierarchical_all_celltypes_lr_lr(celltype_multi_cci, plot_row = 1)
#' plot_hierarchical_all_celltypes_lr_lr(celltype_multi_cci, plot_row = 'SEMA3F-NRP1-SEMA3F-PLXNA3')

plot_hierarchical_all_celltypes_lr_lr <- function(celltype_multi_cci, plot_row = NULL){
  interaction_scores = celltype_multi_cci$interaction_scores
  double_lrpais = celltype_multi_cci$double_lrpais
  triple_cell_pairs = celltype_multi_cci$triple_cell_pairs
  interaction_scores = interaction_scores[plot_row,,drop = F]
  if(is.numeric(plot_row)){
    lr_lr_name = double_lrpais$LRLR[plot_row]
  }else{
    lr_lr_name = plot_row
  }
  
  plot1 = interaction_scores[1,]
  plot1 = data.frame(triple_cell_pairs[,c(1:3)], 'score' = plot1)
  plot1 = plot1[plot1$score > 0,]
  all_cell_types <- unique(c(plot1$ct1, plot1$ct2, plot1$ct3)) %>% sort()
  nodes <- data.frame(cell_type = all_cell_types,y = seq_along(all_cell_types))
  nodes_full <- expand.grid(cell_type = all_cell_types, x = 1:3,stringsAsFactors = FALSE) %>% 
    left_join(nodes, by = "cell_type")
  ct1s = all_cell_types[all_cell_types %in% plot1$ct1]
  ct2s = all_cell_types[all_cell_types %in% plot1$ct2]
  ct3s = all_cell_types[all_cell_types %in% plot1$ct3]
  
  nodes_full = nodes_full[(nodes_full$cell_type %in% ct1s & nodes_full$x == 1) |
                            (nodes_full$cell_type %in% ct2s & nodes_full$x == 2) |
                            (nodes_full$cell_type %in% ct3s & nodes_full$x == 3),]
  edges1 <- plot1 %>% select(from = ct1, to = ct2, score) %>% mutate(x_start = 1, x_end = 2) %>% 
    mutate(log_score = log10(score + 1e-10)) %>% 
    left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
    left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
  edges2 <- plot1 %>% select(from = ct2, to = ct3, score) %>% mutate(x_start = 2, x_end = 3) %>% 
    mutate(log_score = log10(score + 1e-10)) %>% 
    left_join(nodes_full, by = c("from" = "cell_type", 'x_start' = 'x')) %>% rename(y_start = y) %>% 
    left_join(nodes_full, by = c("to" = "cell_type", "x_end" = "x")) %>% rename(y_end = y)
  
  plot = ggplot() +
    geom_point(data = nodes_full, aes(x = x, y = y), size = 5, fill = "white", shape = 21, stroke = 1.2)+ # 绘制三列节点
    geom_segment(data = edges1, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                 linewidth = 1.2,
                 lineend = "round")+ # 绘制第一列和第二列的连边和分数
    geom_segment(data = edges2, aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = log_score, alpha = log_score),
                 linewidth = 1.2,
                 lineend = "round")+ # 绘制第二列和第三列的连边和分数
    geom_text_repel(
      data = nodes_full,
      aes(x = x, y = y, label = cell_type),
      direction = "y",  # 仅垂直方向调整
      nudge_y = 0.3,    # 垂直偏移量
      segment.color = NA,
      size = 3.5,
      fontface = "bold"
    ) +
    scale_color_gradientn(
      colors = c("darkblue","#F0F0F0", "#CB3425"),
      name = "Log10(Score)",
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
    scale_x_continuous(
      breaks = 1:3,
      labels = c("CT1 (Sender)", "CT2 (Relay)", "CT3 (Receiver)"),
      limits = c(0.8, 3.2)
    ) +
    labs(title = paste0(lr_lr_name,' interactions of cell types'), x = "", y = "") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  return(plot)
}

