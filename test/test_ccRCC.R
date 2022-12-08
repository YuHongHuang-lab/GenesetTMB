#library(devtools)
#devtools::load_all()

library(GSVA)


# data prepare
mut_exp = readxl::read_xlsx("E:/source_code/data/NIHMS1611472_ccRCC_MAF.xlsx",
                            col_names = TRUE)
colnames(mut_exp) = mut_exp[1, ]
mut_exp = mut_exp[2:nrow(mut_exp), ]

GEP = readxl::read_xlsx("E:/source_code/data/NIHMS1611472_ccRCC_GEP.xlsx")
colnames(GEP) = GEP[1, ]
GEP = GEP[2:nrow(GEP), ]
rownames(GEP) = NULL
GEP = GEP%>%distinct(gene_name, .keep_all = TRUE)%>%
  column_to_rownames(var = "gene_name")


sample_info = read.csv("E:/source_code/data/NIHMS1611472_ccRCC_phenotype.csv")
colnames(sample_info) = sample_info[1, ]
sample_info = sample_info[2:nrow(sample_info), ]

# ID transformer
sample_info = sample_info[!is.na(sample_info$MAF_Tumor_ID) & !is.na(sample_info$RNA_ID), ]
rownames(sample_info) = NULL

sample_info_1 = column_to_rownames(sample_info, var = 'MAF_Tumor_ID')
sample_info_2 = column_to_rownames(sample_info, var = 'RNA_ID')
sample_info_3 = column_to_rownames(sample_info, var = "SUBJID")

# MAF
mut_exp$Tumor_Sample_Barcode = sample_info_1[mut_exp$Tumor_Sample_Barcode, "SUBJID"]
mut_exp = mut_exp[!is.na(mut_exp$Tumor_Sample_Barcode), ]

# GEP
GEP_id = colnames(GEP)
GEP_id_new = sample_info_2[GEP_id, "SUBJID"]
GEP = GEP[, -which(is.na(GEP_id_new))]
GEP_id_new = GEP_id_new[!is.na(GEP_id_new)]
colnames(GEP) = GEP_id_new

# com id
com_id = intersect(unique(mut_exp$Tumor_Sample_Barcode), colnames(GEP))
mut_exp = mut_exp[mut_exp$Tumor_Sample_Barcode%in%com_id, ]
mut_exp = mut_exp[mut_exp$Variant_Classification!="Silent", ]
mut_exp = mut_exp[, 1:2]
GEP = GEP[, com_id]
for (name in colnames(GEP)){
  GEP[, name] = as.numeric(GEP[, name])
}
GEP = rownames_to_column(GEP, var = 'geneid')


sample_info = sample_info_3[com_id, c("ORR", "OS", "OS_CNSR")]
rm(sample_info_1, sample_info_2, sample_info_3)

ccRCC_PathObj = CreatePathObject(gmt_path="E:/source_code/data/h.all.v7.4.symbols.gmt",
                                gene_exp_counts=GEP,
                                gene_map='E:/source_code/data/Homo_sapiens.GRCh38.99.chr.gtf',
                                gene_exp_TPM=GEP,
                                ppi_data_path="E:/source_code/data/ppi_dataset_all.txt",
                                cor_method=1,
                                mut_exp=mut_exp,
                                count_log=FALSE,
                                TPM_log=FALSE,
                                GeneSetCollect=FALSE,
                                sep = "\t")
ccRCC_PathObj = DegreePathmutscore(ccRCC_PathObj)


ccRCC_IWHMB = ccRCC_PathObj@module1$pathway_mutscore_stdrow
ccRCC_IWHMB = t(ccRCC_IWHMB)
ccRCC_IWHMB = as.data.frame(ccRCC_IWHMB)

sample_info = cbind(sample_info, ccRCC_IWHMB[rownames(sample_info), c("HALLMARK_HEDGEHOG_SIGNALING",
                                                                      "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                                                                      "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                                                      "HALLMARK_P53_PATHWAY")])
sample_info = sample_info[sample_info$ORR!="NE", ]
boxplot(HALLMARK_P53_PATHWAY~ORR, data = sample_info)
sample_info_2 = sample_info[sample_info$ORR%in%c("CR", "CRPR", "SD"), ]
sample_info_2$ORR_2 = if_else(sample_info_2$ORR%in%c("CR", "CRPR", "PR"), 1, 2)
boxplot(HALLMARK_HEDGEHOG_SIGNALING~ORR_2, data = sample_info_2)
t.test(HALLMARK_HEDGEHOG_SIGNALING~ORR_2, data = sample_info_2)

sample_info_2$shh_2 = if_else(sample_info_2$HALLMARK_HEDGEHOG_SIGNALING>0, 1, 2)
table(sample_info_2$shh_2, sample_info_2$ORR_2)


community_1 = readRDS("E:/source_code/step4_RWR/cluster_df.rds")
community_1 = community_1[community_1$cluster==1, "geneid"]
community_1 = c(community_1, "CD276")



#GSVA
GEP = column_to_rownames(GEP, var = "geneid")
GEP = log(GEP+1)
GEP = as.matrix(GEP)

Comm_1 = gsva(GEP, list(ECM=community_1), method="gsva", kcdf="Gaussian",
              abs.ranking	= TRUE)
Comm_1 = as.data.frame(t(Comm_1))

sample_info$ECM = Comm_1[rownames(sample_info), ]
boxplot(ECM~ORR, data = sample_info)

plot(sample_info$ECM, sample_info$HALLMARK_HEDGEHOG_SIGNALING)

# all Hallmark Gene Set




