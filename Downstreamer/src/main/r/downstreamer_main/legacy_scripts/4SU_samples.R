table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/SRA_Metadata/metadata_sra1.txt", delim = "\t", quote = "", guess_max = 20000)
sraMeta1 <- as.data.frame(table_tmp[,-1])
rm(table_tmp)

table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/SRA_Metadata/metadata_sra2.txt", delim = "\t", quote = "", guess_max = 20000)
sraMeta2 <- as.data.frame(table_tmp[,-1])
rm(table_tmp)


sraSharedCol <- intersect(colnames(sraMeta2), colnames(sraMeta1))

sraMeta <- rbind(sraMeta1[,sraSharedCol], sraMeta2[,sraSharedCol])


dim(sraMeta)



columnsToCheckForSu4 <- c("sra.library_construction_protocol", 
                          "sra.study_abstract", 
                          "sra.experiment_title", 
                          "sra.design_description",
                          "sra.sample_description",
                          "sra.library_construction_protocol",
                          "sra.sample_attributes",
                          "sra.sample_title")


match4su <- Reduce("|", lapply(columnsToCheckForSu4, function(col, pattern, ...){grepl(pattern,sraMeta[,col],...)}, pattern = "4su|thiouridine", ignore.case = T))
sum(match4su)

write.table(sraMeta[match4su,"external_id",drop =F], file = "4su_samples.txt", sep = "\t", quote  = FALSE, row.names = F)



