
hpo.matrix <- data.frame()
input.path <- "~/Desktop/depict2/hpo_predictions/"
input.files <- list.files(input.path, pattern="height*")

for (file in input.files) {
  cat(file, "\n")
  data <- read.table(paste0(input.path, file), header=T, row.names =1)
  if (ncol(hpo.matrix) == 0) {
    hpo.matrix <- data
  } else {
    hpo.matrix <- cbind(hpo.matrix, data[rownames(hpo.matrix),])
  }
}

colnames(hpo.matrix) <- gsub("\\_hpoAUC.txt\\_hpo.txt", "", gsub("2018\\_.*\\_hg19", "", input.files))


hpo.cors <- cor(hpo.matrix, use="pairwise.complete.obs")


pairs(hpo.matrix)

hist(apply(hpo.matrix, 1, sd), breaks=100, xlim=c(0,0.1))