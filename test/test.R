#usethis::use_gpl_license()
#usethis::use_tidy_eval()
#usethis::use_pipe()

#library(devtools)
#devtools::load_all()
#install.packages("rmarkdown")

devtools::install_github("YuHongHuang-lab/IWMHB")


library(IWHMB)


library(ggsci)
my_pal = c(pal_npg()(12), pal_aaas()(12), pal_futurama()(12))
my_pal = my_pal[!is.na(my_pal)]
my_pal = my_pal[!duplicated(my_pal)]


TCGA_PathObj = CreatePathObject(gmt_path="test/data/h.all.v7.4.symbols.gmt",
                                gene_exp_counts='test/data/TCGA_HNSC_gene_exp_count.csv',
                                gene_map='test/data/Homo_sapiens.GRCh38.99.chr.gtf',
                                gene_exp_TPM="test/data/TCGA_HNSC_gene_exp_TPM_tumor.csv",
                                ppi_data_path="test/data/ppi_dataset_all.txt",
                                cor_method=1,
                                mut_exp="test/data/TCGA-HNSC_UCSC.varscan2_snv.tsv.gz",
                                count_log=FALSE,
                                TPM_log=FALSE,
                                GeneSetCollect=FALSE,
                                sep = "\t")



TCGA_PathObj = NetFeatureCalculate(TCGA_PathObj)

#small world
plot_point_samllworld(TCGA_PathObj, my_pal=my_pal)

#Significance
plot_point_fetsig(TCGA_PathObj, my_pal=my_pal)

TCGA_PathObj = DegreePathmutscore(TCGA_PathObj)

TCGA_PathObj = ConsensusCluster(obj = TCGA_PathObj, maxK = 9, repcount = 50, pSample = 0.8, pRow = 1, distance = 'euclidean', clusterAlg = 'km',
                                seed = 123456)
TCGA_PathObj = calculate_ICI(TCGA_PathObj)
TCGA_PathObj = Hclust(TCGA_PathObj, nboot = 10)

#plot cluster
TCGA_PathObj = SelectBestClusterNum(obj = TCGA_PathObj,
                                    method="hclust",
                                    k=12)

#somatic mutation
mut_exp_TCGA = fread("test/data/TCGA.HNSC.varscan.5296cf00-4d8c-4db3-80d7-930a4b44f90d.DR-10.0.somatic.maf",
                     header = TRUE)
mut_exp_TCGA = mut_exp_TCGA[mut_exp_TCGA$Variant_Classification!="Silent", ]
classification = as.data.frame(table(mut_exp_TCGA$Variant_Classification))
mut_exp_TCGA$mut_classifcation = if_else(mut_exp_TCGA$Variant_Classification%in%c("3'Flank",
                                                                                  "3'UTR",
                                                                                  "5'Flank",
                                                                                  "5'UTR",
                                                                                  "Intron",
                                                                                  "Translation_Start_Site"),
                                         "non_coding", "coding")
mut_exp_TCGA = mut_exp_TCGA[, c("Hugo_Symbol", "mut_classifcation")]

mut_classifcation = mut_exp_TCGA

#plot coding and non coding
plot_mut_classification(obj = TCGA_PathObj, my_pal, mut_classifcation)


#plot mut in network
plot_mut_distribution(obj = TCGA_PathObj, my_pal, mut_classifcation)

#plot mut enrich
mut_enrich(obj = TCGA_PathObj, my_pal=my_pal, mut_classifcation=mut_classifcation,
           degree=2)


#mean sample geneset mut
plot_bar = function(df, my_pal){

  colnames(df) = c("pathid", "FDR", "cluster")
  #df$FDR = -log10(df$FDR)
  df = arrange(df, desc(FDR))
  df$pathid = factor(df$pathid,
                     levels = df$pathid)

  theme_main = function(){
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8, colour = "black"),
      axis.text.y = element_text(size = 8, colour = "black"),
      panel.border = element_rect(size = 0.7, linetype = 'solid', colour = "black"),
      #legend.position = 'right'
    )
  }

  theme_legend = function(){
    theme(
      legend.position = 'top',
      legend.direction = "horizontal"
    )
  }

  ggplot(df, aes(pathid, FDR))+geom_col(aes(fill=cluster), width = 0.4)+
    scale_fill_manual(values = c(my_pal[1], "gray60", my_pal[2]))+
    ylab("FDR")+xlab("")+scale_y_continuous(expand = c(0, 0))+
    theme_bw()+
    ggpubr::rotate_x_text(angle = 45)+
    geom_hline(yintercept = c(1.5, 2),
               linetype="dotdash", color='gray30')+
    theme_main()+theme_legend()
}





library(GSEABase)
gmt = getGmt("test/data/h.all.v7.4.symbols.gmt")
name_id = names(gmt)
gmt = lapply(gmt, function(x)x@geneIds)
names(gmt) = name_id
mut_exp = fread("test/data/TCGA-HNSC_UCSC.varscan2_snv.tsv.gz", header = TRUE,
                data.table = FALSE)
mut_exp = mut_exp[mut_exp$effect!="synonymous_variant", ]
mut_exp = mut_exp[, c(1, 2)]
mut_list = split(mut_exp$gene, mut_exp$Sample_ID)

outTab_list = list()
for (geneset in names(gmt)){
  outTab_tmp = data.frame()
  i = 0
  for (id in names(mut_list)){
    i = i+1
    geneset_tmp = gmt[[geneset]]
    mut_genes = mut_list[[id]]
    com_id = intersect(geneset_tmp, mut_genes)
    outTab_tmp = rbind(outTab_tmp, cbind(sample_id=id, genesetid=geneset, len=length(com_id)))
    tmp_id = str_c(geneset, i, sep = "_")
    print(tmp_id)
  }
  outTab_list[[geneset]] = outTab_tmp
}

outTab_list = lapply(outTab_list, function(x){
  x$len = as.numeric(x$len)
  x = x[x$len>0, ]
  x = x[, c(2,3)]
  return(x)
})

outTab_df = do.call("rbind", outTab_list)
outTab_df_mean = outTab_df%>%group_by(genesetid)%>%summarise(mean=mean(len))


outTab_df_mean$direction = case_when(outTab_df_mean$mean>=2~"high",
                                     outTab_df_mean$mean>=1.5&
                                       outTab_df_mean$mean<2~"middle",
                                     outTab_df_mean$mean<1.5~"low")
outTab_df_mean_TCGA = outTab_df_mean
outTab_df_mean_TCGA$genesetid = str_remove_all(outTab_df_mean_TCGA$genesetid,
                                               "HALLMARK_")

plot_bar(df = outTab_df_mean_TCGA, my_pal = my_pal)


