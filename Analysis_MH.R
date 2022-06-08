library(tidyverse)
library(stringi)

df <- read.csv(unzip("./Alleles_frequency_table.zip"), sep = "\t", header = T)

## 选择indel序列
df <- subset(df, n_deleted > 0)

aligned <- df$Aligned_Sequence
reference <- df$Reference_Sequence

Del_index1 <- data.frame(stri_locate_all(pattern = "-", aligned, fixed = TRUE))[, 1] # 调取‘--’的位置信息
Del_index2 <- split(Del_index1, cumsum(seq_along(Del_index1) %in% (which(diff(Del_index1) > 1) + 1))) # 根据位置切分‘---’

I <- max(data.frame(Del_index2[1])) #序列截取起始位点
II <- (max(data.frame(Del_index2[3])) - min(data.frame(Del_index2[3]))) # 序列截取终止位点
