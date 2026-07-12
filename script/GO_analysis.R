# Load required libraries
library(genekitr)
library(geneset)
library(ggplot2)
library(Cairo)

plotEnrichAdv2 <- function(up_enrich_df,
                          down_enrich_df,
                          plot_type = c("one", "two"),
                          term_metric = c("FoldEnrich", "GeneRatio", "Count", "RichFactor"),
                          stats_metric = c("p.adjust", "pvalue", "qvalue"),
                          wrap_length = NULL,
                          xlim_left = NULL,
                          xlim_right = NULL,
                          color,
                          ...) {
  #--- args ---#
  lst <- list(...) # store outside arguments in list
  stopifnot(
    "The input enrichment analysis is not data frame!" =
      is.data.frame(up_enrich_df) | is.data.frame(down_enrich_df)
  )
  plot_type <- match.arg(plot_type)
  if(any(grepl("nes",colnames(up_enrich_df),ignore.case = T))) term_metric <- "Count"
  if(any(grepl("nes",colnames(down_enrich_df),ignore.case = T))) term_metric <- "Count"
  
  #--- codes ---#
  tryCatch(
    {
      up_enrich_df$Description <- stringr::str_replace(up_enrich_df$Description, "^\\w{1}", toupper)
      down_enrich_df$Description <- stringr::str_replace(down_enrich_df$Description, "^\\w{1}", toupper)
    },
    error = function(e) {
      message(paste0(
        "We need the 'Description' column which means pathway detailed description", "\n",
        "Maybe you should rename the column name..."
      ))
    }
  )
  
  x_label <- ifelse(stats_metric == "pvalue", "-log10(Pvalue)",
                    ifelse(stats_metric == "p.adjust", "-log10(P.adjust)", "-log10(FDR)")
  )
  
  
  #--- plot ---#
  if (plot_type == "two") {
    if (missing(color)) color <- c("#3665a6", "#d5e4ef", "#a32a31", "#f7dcca")
    
    left <- suppressMessages(plotEnrich(down_enrich_df,
                                        plot_type = "bar",
                                        term_metric = term_metric,
                                        stats_metric = stats_metric,
                                        up_color = color[1], down_color = color[2], ...
    ) +
      #scale_y_discrete(limits = rev) +
      scale_x_reverse() +
      theme(
        axis.title.y = element_blank(),
        legend.position = c(0.2, 0.2)
      ))
    
    right <- suppressMessages(plotEnrich(up_enrich_df,
                                         plot_type = "bar",
                                         term_metric = term_metric,
                                         stats_metric = stats_metric,
                                         up_color = color[3], down_color = color[4], ...
    ) +
      scale_y_discrete(position = "right") +
      theme(
        axis.title.y = element_blank(),
        legend.position = c(0.8, 0.2)
      ))
    
    p <- cowplot::plot_grid(left, right, ncol = 2)
  } else {
    if (missing(color)) color <- c("#3665a6", "#a32a31")
    if (!"main_text_size" %in% names(lst)) lst$main_text_size <- 8
    
    up_go <- dplyr::mutate(up_enrich_df, change = "up")
    down_go <- dplyr::mutate(down_enrich_df, change = "down")
    df <- rbind(up_go, down_go) %>%
      dplyr::mutate(new_x = ifelse(change == "up", -log10(eval(parse(text = stats_metric))), log10(eval(parse(text = stats_metric))))) %>%
      dplyr::arrange(change, new_x) %>%
      dplyr::mutate(Description = factor(Description,
                                         levels = unique(Description),
                                         ordered = TRUE
      ))
    
    if (is.null(xlim_left) & is.null(xlim_right)) {
      tmp <- with(df, labeling::extended(range(new_x)[1], range(new_x)[2], m = 5))
      lm <- tmp[c(1, length(tmp))]
      lm <- c(floor(min(df$new_x)), ceiling(max(df$new_x)))
    } else {
      if (is.null(xlim_left)) xlim_left <- abs(floor(min(df$new_x)) + 1)
      if (is.null(xlim_right)) xlim_right <- 15
      tmp <- seq(-abs(xlim_left), xlim_right, 10)
      lm <- c(-abs(xlim_left), xlim_right)
    }
    
    p <- suppressMessages(ggplot(df, aes(x = Description, y = new_x, fill = change)) +
                            geom_bar(stat = "identity", width = 0.8) +
                            scale_fill_manual(
                              values = color,
                              name = "change",
                              labels = c("Down-regulated pathways", "Up-regulated pathways")
                            ) +
                            guides(fill = guide_legend(reverse = TRUE)) +
                            scale_x_discrete(expand = expansion(add = .5)) +
                            coord_flip() +
                            scale_y_continuous(
                              breaks = tmp, labels = abs(tmp),
                              limits = lm
                            ) +
                            geom_text(
                              data = subset(df, change == "up"),
                              aes(x = Description, y = 0, label = paste0(Description, "  "), color = change),
                              size = lst$main_text_size / 3.6,
                              hjust = "inward", show.legend = FALSE
                            ) +
                            geom_text(
                              data = subset(df, change == "down"),
                              aes(x = Description, y = 0, label = paste0("  ", Description), color = change),
                              size = lst$main_text_size / 3.6, hjust = "outward", show.legend = FALSE
                            ) +
                            scale_colour_manual(values = c("black", "black")) +
                            labs(x = "", y = x_label) +
                            plot_theme(remove_grid = T, remove_legend = F, ...) +
                            theme(
                              axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              legend.title = element_blank()
                            ))
  }
  
  # wrap long text
  if (!is.null(wrap_length) & is.numeric(wrap_length)) {
    p <- p + scale_y_discrete(labels = text_wraper(wrap_length))
  }
  
  return(p)
}

genders <- c('knockdown_flow')

conditions <- c('LPHN2')


# Prepare gene sets (same for all combinations)
gene_sets <- list(
  GO_BP = geneset::getGO(org = "human", ont = "bp"),
  GO_MF = geneset::getGO(org = "human", ont = "mf"),
  GO_CC = geneset::getGO(org = "human", ont = "cc"),
  KEGG_PATHWAY = geneset::getKEGG(org = "mouse",category = "pathway")
  DisGeNET = geneset::getEnrichrdb(org = "human", library = "DisGeNET")
)

# Create a single output directory if it doesn't exist
output_dir <- "HERE IS PATH FOR OUTPUT FOLDER"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through all conditions and genders
for (condition in conditions) {
  for (gender in genders) {
    # Construct file path for the input gene list
    file_path <- sprintf("HERE IS PATH FOR INPUT CSV FILE", condition, gender)
    
    # Check if file exists before proceeding
    if (!file.exists(file_path)) {
      cat(sprintf("File not found: %s\n", file_path))
      next # Skip to the next iteration
    }
    
    cat(sprintf("Processing: %s - %s\n", condition, gender))
    
    # 1st step: Read gene lists from CSV
    gene_df <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Ensure columns are named 'Up' and 'Down'
    up_genes <- na.omit(gene_df$Up)
    up_genes <- up_genes[up_genes != ""]
    down_genes <- na.omit(gene_df$Down)
    down_genes <- down_genes[down_genes != ""]
    
    # Combine the two vectors and get only the unique genes
    both_genes <- unique(c(up_genes, down_genes))
    
    # Check if we have genes to analyze
    if (length(up_genes) == 0 && length(down_genes) == 0) {
      cat(sprintf("No genes found for %s - %s\n", condition, gender))
      next
    }
    
    # 3rd and 4th steps: ORA analysis and plotting for each gene set
    for (gs_name in names(gene_sets)) {
      gs <- gene_sets[[gs_name]]
      
      # ORA analysis
      up_go <- genORA(up_genes, geneset = gs, p_cutoff = 0.4, q_cutoff = 0.4)
      down_go <- genORA(down_genes, geneset = gs, p_cutoff = 0.4, q_cutoff = 0.4)
      both_go <- genORA(both_genes, geneset = gs, p_cutoff = 0.4, q_cutoff = 0.4)
      
      # Print dimensions for checking
      cat(sprintf("Gene set: %s\n", gs_name))
      cat("Up-regulated ORA result dim: ", dim(up_go), "\n")
      cat("Down-regulated ORA result dim: ", dim(down_go), "\n")
      cat("Up&Down-regulated ORA result dim: ", dim(both_go), "\n")
      
      # Check if we have results
      if (nrow(up_go) == 0 && nrow(down_go) == 0) {
        cat(sprintf("No significant results for %s - %s - %s\n", condition, gender, gs_name))
        next
      }
      
      # Take top 10 terms
      up_go_10 <- head(up_go, 10)
      down_go_10 <- head(down_go, 10)
      both_go_10 <- head(both_go, 10)
        
      combined_go <- rbind(up_go_10, down_go_10) 
      sorted_combined_go <- combined_go[order(combined_go$FoldEnrich, decreasing = TRUE), ]  
      combined_go_10 <- head(sorted_combined_go, 10)
      
      # Plot 1: "one" type
      p1 <- plotEnrichAdv(
        up_go_10, down_go_10,
        plot_type = "one",
        term_metric = "FoldEnrich",
        stats_metric = "p.adjust",
        xlim_left = 14, xlim_right = 8,
        color = c("#1f77b4", "#d62728")
      ) +
        theme(legend.position = "top")
      
      ggsave(
        filename = file.path(output_dir, sprintf("%s_%s_plot_one_%s.png", condition, gender, gs_name)),
        plot = p1,
        width = 6, height = 6, dpi = 300
      )
      
      # Plot 2: "two" type
      p2 <- plotEnrichAdv2(up_go_10, down_go_10,
                           plot_type = "two",
                           term_metric = "FoldEnrich",
                           stats_metric = "qvalue",
                           main_text_size = 15) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.margin = unit(c(0.75, 2, 0.75, 2), "cm"),
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA))
      
      ggsave(
        filename = file.path(output_dir, sprintf("%s_%s_plot_two_%s.PNG", condition, gender, gs_name)), 
        plot = p2, 
        width = 24, height = 6, dpi = 600
      )
      
      # Plot 3: "two" type
      p3 <- plotEnrich(both_go_10,
                       plot_type = "bar",
                       term_metric = "FoldEnrich",
                       stats_metric = "qvalue",
                       up_color = "#a32a31",  down_color = "#f7dcca", 
                       main_text_size = 12) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.margin = unit(c(0.75, 2, 0.75, 2), "cm"),
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA))
      
      ggsave(
        filename = file.path(output_dir, sprintf("%s_%s_plot_both_%s.png", condition, gender, gs_name)), 
        plot = p3, 
        width = 12, height = 6, dpi = 600
      )

      # Plot 4: "two" type
      p4 <- plotEnrich(combined_go_10,
                       plot_type = "bar",
                       term_metric = "FoldEnrich",
                       stats_metric = "qvalue",
                       up_color = "#a32a31",  down_color = "#f7dcca", 
                       main_text_size = 12) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.margin = unit(c(0.75, 2, 0.75, 2), "cm"),
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA))
      
      ggsave(
        filename = file.path(output_dir, sprintf("%s_%s_plot_combined_%s.png", condition, gender, gs_name)), 
        plot = p4, 
        width = 12, height = 6, dpi = 600
      )  
      # Save the results as CSV
      if (nrow(up_go) > 0) {
        write.csv(up_go, file.path(output_dir, sprintf("%s_%s_upregulated_%s.csv", condition, gender, gs_name)), row.names = FALSE)
      }
      if (nrow(down_go) > 0) {
        write.csv(down_go, file.path(output_dir, sprintf("%s_%s_downregulated_%s.csv", condition, gender, gs_name)), row.names = FALSE)
      }
      if (nrow(both_go) > 0) {
        write.csv(both_go, file.path(output_dir, sprintf("%s_%s_up&downregulated_%s.csv", condition, gender, gs_name)), row.names = FALSE)
      }
    }
  }
}

cat("Analysis complete!\n")