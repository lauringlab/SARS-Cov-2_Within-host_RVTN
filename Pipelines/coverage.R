
library (dplyr)
filename <-snakemake@input[[1]]
filename2 <- sub(".coverage.csv", "", filename)

filename_vec <- strsplit(filename2, split = "/")[[1]]

if (file.info(snakemake@input[[1]])$size == 0) {
	summary.df <-data.frame (sample = filename_vec[4], mean = "0")
        coverage.df <- data.frame ("Seg", "Pos", "Coverage", "sample")
}else {coverage.df <- read.table(snakemake@input[[1]],stringsAsFactors=F,comment.char = '#') 
	coverage.df <- rename (coverage.df, Seg = V1, Pos =V2, Coverage= V3)
	coverage.df$sample <- filename_vec[4]
	Avg <- mean (coverage.df$Coverage)
	summary.df <-  data.frame ("sample"= filename_vec[4], "mean"= Avg)
}
write.table (coverage.df, snakemake@output[[1]], quote=F, row.names=F)

write.table (summary.df, snakemake@output[[2]], quote=F, row.names=F)

