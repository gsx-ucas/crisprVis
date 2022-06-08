#' get amplicon sequence from json file or log file of run parameters
#'
#' @param path path to the output directory of CRISPResso2
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
}
