#' get amplicon sequence from json file or log file of run parameters
#'
#' @param path path to the sample directory which generated by CRISPResso2
#' @return The sequence of \code{amplicon}, that is reference sequence.
#' @export

get_ampseq <- function(path){
  json_file <- paste0(path, "/CRISPResso2_info.json")
  log_file <- paste0(path, "/CRISPResso_RUNNING_LOG.txt")
  if (!file.exists(json_file)) {
    message(paste0("WANING: file can not found: ", json_file))
    message(paste0("WANING: trying use: ", log_file))
    if (!file.exists(log_file)) {
      stop(paste0("ERROR: can not found file '", json_file, "' or '", log_file, "'"))
    }else {
      run_log <- read.table(log_file, sep = "\t", header = F)
      amplicon_seq <- run_log$V1[3] %>% stringr::str_split(" --") %>% unlist %>% grep(pattern = "amplicon_seq", value = T) %>% stringr::str_remove("amplicon_seq ")
    }
  }else {
    run_info <- jsonlite::read_json(path = json_file)
    amplicon_seq <- run_info$running_info$args$value$amplicon_seq
  }
  return(amplicon_seq)
}

#' Plotting the modification count vectors
#'
#' @param path path to the sample directory which generated by CRISPResso2
#' @param start_pos start position of the highlited mutation window
#' @param end_pos end position of the highlited mutation window
#' @export

plot_MCV <- function(path, start_pos, end_pos, show_gRNA = T){
  require(ggplot2)
  require(magrittr)
  mcv_file <-  paste0(path, "/Modification_count_vectors.txt")
  if (!file.exists(mcv_file)) {
    stop(paste0("can not find file: ", mcv_file))
  }
  MCV <- read.table(mcv_file, sep = "\t", row.names = 1, header = T, check.names = F)
  seq_char <- colnames(MCV)
  MCV <- MCV %>% t %>% as.data.frame
  MCV$Seq <- seq_char

  max_y <- max((MCV$Insertions + MCV$Insertions_Left) / MCV$Total, MCV$Deletions / MCV$Total, MCV$Substitutions / MCV$Total) * 100

  p <- ggplot(MCV)+
    geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = 0, ymax = max_y + 5), fill = "#F8F8FF", alpha = 0.1)

  if (show_gRNA) {
    json_file <- paste0(path, "/CRISPResso2_info.json")
    if (!file.exists(json_file)) {
      stop(paste0("file can not found: ", json_file))
    }else {
      run_info <- jsonlite::read_json(path = json_file)
      gRNA_start_pos <- run_info$results$refs$Reference$sgRNA_intervals[[1]][[1]]
      gRNA_end_pos <- run_info$results$refs$Reference$sgRNA_intervals[[1]][[2]]

      p <- p + geom_rect(aes(xmin = gRNA_start_pos, xmax = gRNA_end_pos, ymin = 0, ymax = max_y + 5), fill = "#FFF8DC")
    }
  }

  p <- p + geom_vline(xintercept = start_pos, lty = 2, col = "grey")+
    geom_vline(xintercept = end_pos, lty = 2, col = "grey")+
    geom_line(aes(x = 1:length(Seq), y = (Insertions + Insertions_Left) / Total * 100, color = "Insertions"), group = 1)+
    geom_line(aes(x = 1:length(Seq), y = Deletions / Total * 100, color = "Deletions"), group = 1)+
    geom_line(aes(x = 1:length(Seq), y = Substitutions / Total * 100, color = "Substitutions"), group = 1)+
    labs(color = NULL, x = "Reference amplicon position (bp)", y = "Sequences: % Total (no.)")+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    theme_test()+
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.margin = margin(5, 5, 5, 5),
          text = element_text(size = 15))
  return(p)
}
