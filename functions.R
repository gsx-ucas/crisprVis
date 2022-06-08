# library(BSgenome)
# library(Biostrings)
# library(GenomicRanges)
# library(stringi)
# library(stringr)
# library(stats)
# library(utils)
# library(dplyr)
# library(readr)
# library(tidyr)

gather <- function (file) 
{
  dfx <- data.frame(read.delim(file, header = T, sep = "\t"))
  dfx <- dfx[dfx$Read_Status == "MODIFIED" & dfx$n_deleted > 0 & 
               dfx$n_inserted < 1 & dfx$n_mutated < 1 & grepl("-", 
                                                              dfx$Aligned_Sequence, fixed = T), ]
  buffer <- as.character("-----")
  dfx$Aligned_Sequence <- paste0(buffer, dfx$Aligned_Sequence, 
                                 buffer)
  dfx$Reference_Sequence <- paste0(buffer, dfx$Reference_Sequence, 
                                   buffer)
  return(dfx)
}

index <- function (aligned, reference) 
{
  Del_index1 <- (data.frame(stri_locate_all(pattern = "-", 
                                            aligned, fixed = TRUE))[, 1])
  Del_index2 <- split(Del_index1, cumsum(seq_along(Del_index1) %in% 
                                           (which(diff(Del_index1) > 1) + 1)))
  I <- max(data.frame(Del_index2[1]))
  II <- (max(data.frame(Del_index2[3])) - min(data.frame(Del_index2[3])))
  DNA_1 <- stri_sub(aligned, I + 1)
  DNA_2 <- substr(DNA_1, 1, nchar(DNA_1) - (II + 1))
  ref1 <- stri_sub(reference, I + 1)
  Ref <- substr(ref1, 1, nchar(ref1) - (II + 1))
  Del_index3 <- data.frame(stri_locate_all(pattern = "-", 
                                           DNA_2, fixed = TRUE))
  Del_index4 <- Del_index3[, 1]
  index_out <- list(delSeq = DNA_2, refSeq = Ref, delInd = Del_index4)
  return(index_out)
}

is.simple <- function (index_out) 
{
  runLen <- rle(diff(index_out$delInd))
  delInd <- index_out$delInd
  simple_out <- any(runLen$lengths >= (length(delInd) - 1) & 
                      runLen$values == 1) | length(delInd) == 1
  return(simple_out)
}

is.mut.prep <- function (index_out) 
{
  delInd <- index_out$delInd
  Del_index_L <- min(delInd)
  Del_index_R <- max(delInd)
  Seq <- DNAString(index_out$delSeq)
  Ref <- DNAString(index_out$refSeq)
  Del_test <- c(subseq(Seq, start = Del_index_L - 10, end = Del_index_L - 
                         1), subseq(Seq, start = Del_index_R + 1, end = Del_index_R + 
                                      10))
  Ref_test <- c(subseq(DNAString(Ref), start = Del_index_L - 
                         10, end = Del_index_L - 1), subseq(DNAString(Ref), start = Del_index_R + 
                                                              1, end = Del_index_R + 10))
  mut_prep_out <- list(Del_test = Del_test, Ref_test = Ref_test)
  return(mut_prep_out)
}

is.mutated <- function (a, b, exclude = "-") 
{
  a <- as.character(a)
  b <- as.character(b)
  if (nchar(a) != nchar(b)) 
    stop("Lengths of input strings differ. Please check your input.")
  split_seqs <- strsplit(c(a, b), split = "")
  only.diff <- (split_seqs[[1]] != split_seqs[[2]])
  only.diff[(split_seqs[[1]] %in% exclude) | (split_seqs[[2]] %in% 
                                                exclude)] <- NA
  diff.info <- data.frame(which(is.na(only.diff) | only.diff), 
                          split_seqs[[1]][only.diff], split_seqs[[2]][only.diff])
  names(diff.info) <- c("position", "poly.seq.a", "poly.seq.b")
  is.mut.seq <- paste(diff.info$poly.seq.a, sep = "", collapse = "")
  is.mut <- is.mut.seq != ""
  return(is.mut)
}

MHQuant.prep1 <- function (index_out) 
{
  Seq <- index_out$delSeq
  Ref <- index_out$refSeq
  Del_index_L <- min(index_out$delInd)
  Del_index_R <- max(index_out$delInd)
  sizeOfDel <- length(index_out$delInd)
  L <- subseq(Seq, start = 1, end = Del_index_L - 1)
  R <- subseq(Seq, start = Del_index_R + 1, end = nchar(Seq))
  longDel <- (sizeOfDel >= nchar(L) | sizeOfDel >= nchar(R))
  shortSizeSeq <- (min(nchar(L), nchar(R)) - 1)
  fullSizeSeq <- sizeOfDel
  if (longDel) {
    MHSeqSize <- shortSizeSeq
  }
  else {
    MHSeqSize <- fullSizeSeq
  }
  return(MHSeqSize)
}

MHQuant.prep2 <- function (index_out, MHSeqSize) 
{
  Seq <- index_out$delSeq
  Ref <- index_out$refSeq
  Del_index_L <- min(index_out$delInd)
  Del_index_R <- max(index_out$delInd)
  L <- subseq(Seq, start = 1, end = Del_index_L - 1)
  R <- subseq(Seq, start = Del_index_R + 1, end = nchar(Seq))
  seq1A <- DNAString(subseq(Seq, start = Del_index_L - MHSeqSize, 
                            end = Del_index_L - 1))
  seq1B <- DNAString(subseq(Ref, start = Del_index_L, end = Del_index_L + 
                              MHSeqSize - 1))
  seq2A <- DNAString(subseq(Ref, start = Del_index_R - MHSeqSize + 
                              1, end = Del_index_R))
  seq2B <- DNAString(subseq(Seq, start = Del_index_R + 1, 
                            end = Del_index_R + MHSeqSize))
  MHQuant_seqs <- list(seq1A = seq1A, seq1B = seq1B, seq2A = seq2A, 
                       seq2B = seq2B, MHSeqSize = MHSeqSize)
  return(MHQuant_seqs)
}

MHQuant_sub <- function (MHQuant_seqs) 
{
  MH_A <- MHQuant_L(MHQuant_seqs)
  MH_B <- MHQuant_R(MHQuant_seqs)
  max_MH <- max(MH_A, MH_B)
  MHQuant_out <- list(max_MH = max_MH, MH_A = MH_A, MH_B = MH_B)
  return(MHQuant_out)
}

MHSeq.return <- function (MHQuant_out, index_out) 
{
  delInd <- index_out$delInd
  Del_index_L <- min(delInd)
  Del_index_R <- max(delInd)
  Seq <- DNAString(index_out$delSeq)
  max_MH <- MHQuant_out$max_MH
  MH_A <- MHQuant_out$MH_A
  MH_B <- MHQuant_out$MH_B
  if (max_MH == 0) {
    MHSeq <- "No_MH"
  }
  if (max_MH > 0 & max_MH == MH_A) {
    MHSeq <- as.character(subseq(Seq, start = Del_index_L - 
                                   max_MH, end = Del_index_L - 1))
  }
  if (max_MH > 0 & max_MH == MH_B) {
    MHSeq <- as.character(subseq(Seq, start = Del_index_R + 
                                   1, end = (Del_index_R + max_MH)))
  }
  return(MHSeq)
}

altMH <- function (MHSeq_out, index_out) 
{
  NULL
  testSeq <- index_out$delSeq
  refSeq <- index_out$refSeq
  MH <- MHSeq_out
  Del_index <- index_out$delInd
  Del_index_L <- min(index_out$delInd)
  Del_index_R <- max(index_out$delInd)
  refDel <- subseq(DNAString(refSeq), start = Del_index_L, 
                   end = Del_index_R)
  altMH_out <- countPattern(MH, refDel) - 1
  return(altMH_out)
}

MHQuant_L <- function (MHQuant_seqs) 
{
  MH_A <- 0
  for (i in 1:MHQuant_seqs$MHSeqSize) {
    if (subseq(MHQuant_seqs$seq1A, start = i, end = i) == 
        subseq(MHQuant_seqs$seq2A, start = i, end = i)) {
      MH_A <- MH_A + 1
    }
    else {
      MH_A <- 0
    }
  }
  return(MH_A)
}

MHQuant_R <- function (MHQuant_seqs) 
{
  MH_B <- 0
  for (i in 1:MHQuant_seqs$MHSeqSize) {
    if (subseq(MHQuant_seqs$seq1B, start = i, end = i) == 
        subseq(MHQuant_seqs$seq2B, start = i, end = i)) {
      MH_B <- MH_B + 1
    }
    else {
      MH_B <- MH_B
      break
    }
  }
  return(MH_B)
}

MHQuant_CRISPResso2 <- function(directory){
  file <- paste0(directory, "/Alleles_frequency_table.zip")
  if (!file.exists(file)) {
    stop("File Alleles_frequency_table.zip does not exist, please check your output of CRISPResso2 ...")
  }else {
    data <- gather(unzip(file))
    
    strings <- strsplit(data$Aligned_Sequence, "-")
    logical_list <- lapply(strings, function(x){
      str_tab <- table(nchar(x) > 0)
      if (unname(str_tab[2]) > 2) {
        return(FALSE)
      }else {
        return(TRUE)
      }
    }) %>% unlist
    data <- data[logical_list, ]
    
    results_int <- NULL
    
    for (i in 1:nrow(data)) {
      df <- data[i, ]
      index_out <- index(df$Aligned_Sequence, df$Reference_Sequence)
      simple_out <- is.simple(index_out)
      if (!simple_out) {
        next
      }else {
        mut_prep_out <- is.mut.prep(index_out)
        is_mutated <- is.mutated(a = mut_prep_out$Del_test, 
                                 b = mut_prep_out$Ref_test, exclude = "-")
        if (is_mutated) {
          next
        }else {
          MHSeqSize <- MHQuant.prep1(index_out)
          MHQuant_seqs <- MHQuant.prep2(index_out, MHSeqSize)
          MHQuant_out <- MHQuant_sub(MHQuant_seqs)
          MHSeq_out <- MHSeq.return(MHQuant_out, index_out)
        }
        if (MHSeq_out == "No_MH") {
          altMH_out <- "No_MH"
        }else {
          altMH_out <- altMH(MHSeq_out, index_out)
        }
        
        results_int <- rbind(results_int, data.frame(MutantSequence = index_out$delSeq, 
                                                     ReferenceSequence = index_out$refSeq, SizeOfDeletion = df$n_deleted, 
                                                     NumberOfReads = df$X.Reads, MH_amount = MHQuant_out$max_MH, 
                                                     MH_sequence = MHSeq_out, altMH_count = altMH_out))
        results_int <- results_int[results_int$MH_sequence != "No_MH", ]
      }
    }
    dir.create(paste0(directory, "/MHQuant_results"), recursive = T)
    output_file <- paste0(directory, "/MHQuant_results/MHQuant_results.csv")
    write.csv(results_int, output_file)
    return(results_int) 
  }
}