#' Calculate Individual Weight Hallmark Gene Set mutation burden for every sample
#'
#' @param  obj PathObj object
#' @param ... custom parameters
#'
#' @return PathObj
#' @export
#'
#' @import igraph
#'
setMethod("DegreePathmutscore", "PathObject", function(obj, ...){

  geneset_names = obj@module1$geneset_names
  ppi_exp_cor = obj@module1[["ppi_exp_cor"]]
  sample_mut_list = obj@module1[["sample_mut_list"]]
  ppi_dataset = obj@module1[["ppi_dataset"]]
  geneset = obj@module1[["geneset"]]

  pathway_mutscore = sapply(geneset, .pathway_degreemut, ppi_exp_cor=ppi_exp_cor,
                            sample_mut_list=sample_mut_list,
                            geneset_names=geneset_names,
                            ppi_dataset=ppi_dataset,
                            USE.NAMES = TRUE)
  obj@module1[["pathway_mutscore"]] = pathway_mutscore

  pathway_mutscore_stdrow = apply(pathway_mutscore, 1, function(x){(x-mean(x))/(sd(x))})
  pathway_mutscore_stdrow[is.na(pathway_mutscore_stdrow)] = 0
  obj@module1[["pathway_mutscore_stdrow"]] = pathway_mutscore_stdrow

  return(obj)
})


#' Calculate Network Features of Hallmark Gene Set or Pathways for current cohort
#'
#' @param  obj PathObj object
#' @param ... custom parameters
#'
#' @return new PathObj
#' @export
#'
#'
setMethod("NetFeatureCalculate", "PathObject", function(obj, ...){

  geneset_names = obj@module1$geneset_names
  ppi_exp_cor = obj@module1[["ppi_exp_cor"]]
  sample_mut_list = obj@module1[["sample_mut_list"]]
  ppi_dataset = obj@module1[["ppi_dataset"]]
  geneset = obj@module1[["geneset"]]

  pathway_tmp = geneset[[1]]

  pathway_obj_list = lapply(geneset, .graph_obj_construct,
                            geneset_names=geneset_names,
                            ppi_dataset=ppi_dataset,
                            ppi_exp_cor=ppi_exp_cor)
  pathway_feature_list = lapply(pathway_obj_list, .network_summary)

  obj@module1$pathway_feature = pathway_feature_list

  return(obj)
})




#' ConsinusCluster encapsulates the functions in the ConsinusClusterPlus package
#'
#' @param obj PathObject
#' @param all All parameters are the same as those of the R package ConsensusClusterPlus function
#'
#' @return PathObject
#' @export
#'
setMethod("ConsensusCluster", "PathObject", function(obj, maxK=9, repcount=50, pSample=0.8, weightItems=NULL,
                                                     pRow=1, weightFeatures=NULL, distance='euclidean', clusterAlg='km',
                                                     corUse=NULL, innerLinkage="average", seed=123456,
                                                     out_dir = "test/ConsensusCluster"){
  if (!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  mat = obj@module1$pathway_mutscore_stdrow

  Consensus_obj = ConsensusClusterPlus::ConsensusClusterPlus(mat, maxK = maxK,
                                       reps = repcount,
                                       pItem = pSample,
                                       pFeature = pRow,
                                       clusterAlg = clusterAlg,
                                       innerLinkage = innerLinkage,
                                       seed = seed,
                                       plot = "pdf",
                                       distance = distance,
                                       weightsFeature = weightFeatures,
                                       title = out_dir)
  obj@module1$Consensus_obj = Consensus_obj

  return(obj)
})

#' Calculate ICI index encapsulates the functions in the ConsinusClusterPlus package
#'
#' @param obj PathObject
#' @param ... custom parameters
#'
#' @return PathObject
#' @export
#'
setMethod("calculate_ICI", "PathObject", function(obj, out_dir="test/ConsensusCluster",...){
  Consensus_cluster = obj@module1$Consensus_obj
  Consensus_ICL = ConsensusClusterPlus::calcICL(Consensus_cluster, out_dir,
                              plot = "pdf")

  obj@module1$ICI = Consensus_ICL
  return(obj)
})


#' This function is hierarchical clustering Mainly depends on pvlcust
#'
#' @param obj PathObject
#' @param method.hclust same in pvlcust function
#' @param method.dist same in pvlcust function
#' @param nboot same in pvlcust function
#'
#' @return PathObject
#' @export
#'
setMethod("Hclust", "PathObject", function(obj, method.hclust="complete",
                                           method.dist="euclidean",
                                           nboot = 100){
  mat = obj@module1$pathway_mutscore_stdrow
  fit = pvclust::pvclust(mat, method.hclust = method.hclust,
                method.dist = method.dist, nboot = nboot)

  plot(fit, label=NULL, cex=0.1)
  return(obj)
})


#' SelectBestClusterNum for cluster
#'
#' @param obj PathObject
#' @param method "hclust" or "Consensus"
#' @param k the number of cluster you want
#'
#' @return PathObject
#' @export
#'
setMethod("SelectBestClusterNum", "PathObject", function(obj, method="hclust",
                                                         k=NULL, ...){
  cluster <- NULL
  mat = obj@module1$pathway_mutscore_stdrow
  if (method=="hclust"){

    cluster_row = stats::hclust(stats::dist(t(mat), method = "euclidean"), method = "complete")
    cluster_id = stats::cutree(cluster_row, k = k)

    annotation_col = data.frame(cluster=as.character(cluster_id), row.names = names(cluster_id))
    annotation_col = annotation_col%>%dplyr::arrange(cluster)
    p = pheatmap::pheatmap(mat[, rownames(annotation_col)], show_rownames = TRUE,
                       show_colnames = FALSE,
                       cluster_rows = TRUE,
                       cluster_cols = FALSE,
                       annotation_col = annotation_col)
    print(p)
    obj@module1$sample_info = annotation_col
    return(obj)
  }else if (method=="Consensus"){
    Consensus_cluster = obj@module1$Consensus_obj
    annotation_col = Consensus_cluster[[k]][["consensusClass"]]
    annotation_col = data.frame(cluster=as.character(annotation_col),
                                row.names = names(annotation_col))
    annotation_col = annotation_col%>%dplyr::arrange(cluster)
    p = pheatmap::pheatmap(mat[, rownames(annotation_col)], show_rownames = TRUE,
                       show_colnames = FALSE,
                       cluster_rows = TRUE,
                       cluster_cols = FALSE,
                       annotation_col = annotation_col)
    print(p)
    obj@module1$sample_info = annotation_col
    return(obj)
  }
})


#' This function is plot the distributation of coding and non-coding muation in Hallmark Gene Set or parhway network
#'
#' @param obj PathObject
#' @param my_pal Color table
#' @param mut_classifcation The somatic mutation classification data table of the cohort
#' is composed of two columns, the first column is the gene name, and the second column
#' is the mutation division of the gene
#'
#' @export
#'
#'
#' @import ggplot2
setMethod("plot_mut_classification", "PathObject", function(obj, my_pal, mut_classifcation){
  pathid <- mut_classification <- non_coding_count_freq <- NULL
  stopifnot(class(obj)=="PathObject")
  stopifnot(ncol(mut_classifcation)==2)
  if (is.null(obj@module1$pathway_feature)){
    stop("please compute network features first", call. = FALSE)
  }
  colnames(mut_classifcation) = c('geneid', 'mut_classification')
  pathway_feature = obj@module1$pathway_feature

  path_mut_classifcation = lapply(names(pathway_feature), function(x){
    path_feature_tmp = pathway_feature[[x]]
    feature_tmp = path_feature_tmp$node_df
    path_mut_classifcation_tmp = mut_classifcation[mut_classifcation$geneid%in%feature_tmp$geneid, ]
    path_mut_classifcation_tmp$pathid = x
    return(path_mut_classifcation_tmp)
  })

  path_mut_classifcation = do.call('rbind.data.frame', path_mut_classifcation)
  path_mut_classifcation$pathid = stringr::str_remove_all(path_mut_classifcation$pathid,
                                                 "HALLMARK_")

  path_mut_classifcation = path_mut_classifcation[, 2:3]

  path_mut_classifcation = path_mut_classifcation%>%dplyr::group_by(pathid, mut_classification)%>%
    dplyr::summarise(count=n())

  path_factor_1 = path_mut_classifcation[path_mut_classifcation$mut_classification=="coding", ]%>%
    tibble::column_to_rownames(var = 'pathid')%>%dplyr::select(coding_count=count)
  path_factor_2 = path_mut_classifcation[path_mut_classifcation$mut_classification=="non_coding", ]%>%
    tibble::column_to_rownames(var = 'pathid')%>%dplyr::select(non_coding_count=count)


  if (!purrr::is_empty(setdiff(rownames(path_factor_1), rownames(path_factor_2)))){
    path_mut_classifcation_add = data.frame(pathid=setdiff(rownames(path_factor_1), rownames(path_factor_2)),
                                            mut_classification='non_coding',
                                            count=as.integer(0))
    path_mut_classifcation = as.data.frame(path_mut_classifcation)
    path_mut_classifcation = rbind.data.frame(path_mut_classifcation,
                                              path_mut_classifcation_add)
    path_mut_classifcation = as.data.frame(path_mut_classifcation)
    path_factor_2[setdiff(rownames(path_factor_1), rownames(path_factor_2)), "non_coding_count"] = as.integer(0)

  }

  path_factor_1 = path_factor_1[rownames(path_factor_2), , drop=FALSE]
  path_factor = cbind.data.frame(path_factor_1, path_factor_2)
  path_factor$non_coding_count_freq = path_factor$non_coding_count/(path_factor$coding_count+path_factor$non_coding_count)
  path_factor = tibble::rownames_to_column(path_factor, var = 'pathid')

  path_factor = path_factor%>%arrange(non_coding_count_freq)
  path_factor = factor(path_factor$pathid, levels = path_factor$pathid)

  path_mut_classifcation$pathid = factor(path_mut_classifcation$pathid,
                                         levels = levels(path_factor))

  path_mut_classifcation = arrange(path_mut_classifcation, pathid)

  theme_main = function(){
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #axis.text.x = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      panel.border = element_rect(size = 0.7, linetype = 'solid', colour = "black"),
      legend.position = 'right'
    )
  }
  ggplot(path_mut_classifcation, aes(x = pathid,
                                     y = count,
                                     fill = mut_classification))+
    geom_bar(stat = "identity", width = 0.8, position = "fill")+
    scale_fill_manual(values = my_pal[1:length(unique(path_mut_classifcation$mut_classification))])+
    labs(x="", y="Ration")+theme_bw()+
    scale_y_continuous(position = "right")+
    theme_main()+ggpubr::rotate_x_text(angle = 45)


})


#' This function is plot the distributation of muation in or out of Hallmark Gene Set or parhway network
#'
#' @param obj PathObject
#' @param my_pal Color table in coding and non-coding
#' @param mut_classifcation The somatic mutation classification data table of the cohort
#' is composed of two columns, the first column is the gene name, and the second column
#' is the mutation division of the gene
#'
#' @export
#'
#'
#' @import ggplot2
setMethod("plot_mut_distribution", "PathObject", function(obj, my_pal, mut_classifcation){
  pathid <- mut_classification <- innet <- NULL
  stopifnot(class(obj)=="PathObject")
  stopifnot(ncol(mut_classifcation)==2)
  if (is.null(obj@module1$pathway_feature)){
    stop("please compute network features first", call. = FALSE)
  }
  colnames(mut_classifcation) = c('geneid', 'mut_classification')
  pathway_feature = obj@module1$pathway_feature
  geneset = obj@module1$geneset



  path_mut_classifcation = lapply(names(pathway_feature), function(x){
    path_feature_tmp = pathway_feature[[x]]
    feature_tmp_innet = path_feature_tmp$node_df
    feature_tmp_innet = mut_classifcation[mut_classifcation$geneid%in%feature_tmp_innet$geneid, ]
    feature_tmp_innet = feature_tmp_innet%>%
      dplyr::distinct(geneid, .keep_all = TRUE)

    feature_tmp_outnet = geneset[[x]]
    feature_tmp_outnet = setdiff(feature_tmp_outnet, path_feature_tmp$node_df$geneid)
    feature_tmp_outnet = mut_classifcation[mut_classifcation$geneid%in%feature_tmp_outnet, ]
    feature_tmp_outnet = feature_tmp_outnet%>%
      dplyr::distinct(geneid, .keep_all = TRUE)

    feature_tmp_innet$pathid = x
    feature_tmp_outnet$pathid = x

    feature_tmp_innet$mut_classification = 'innet'
    feature_tmp_outnet$mut_classification = 'outnet'

    feature_net = rbind.data.frame(feature_tmp_innet, feature_tmp_outnet)

    return(feature_net)
  })

  path_mut_classifcation = do.call('rbind.data.frame', path_mut_classifcation)
  path_mut_classifcation$pathid = stringr::str_remove_all(path_mut_classifcation$pathid,
                                                 "HALLMARK_")

  path_mut_classifcation = path_mut_classifcation[, 2:3]

  path_mut_classifcation = path_mut_classifcation%>%dplyr::group_by(pathid, mut_classification)%>%
    dplyr::summarise(count=n())

  path_factor_1 = path_mut_classifcation[path_mut_classifcation$mut_classification=="innet", ]%>%
    tibble::column_to_rownames(var = 'pathid')%>%dplyr::select(innet=count)
  path_factor_2 = path_mut_classifcation[path_mut_classifcation$mut_classification=="outnet", ]%>%
    tibble::column_to_rownames(var = 'pathid')%>%dplyr::select(outnet=count)

  path_factor = cbind.data.frame(path_factor_1, path_factor_2)
  path_factor$innet = path_factor$innet/(path_factor$innet+path_factor$outnet)
  path_factor = tibble::rownames_to_column(path_factor, var = 'pathid')

  path_factor = path_factor%>%dplyr::arrange(desc(innet))
  path_factor = factor(path_factor$pathid, levels = path_factor$pathid)

  path_mut_classifcation$pathid = factor(path_mut_classifcation$pathid,
                                         levels = levels(path_factor))

  path_mut_classifcation = dplyr::arrange(path_mut_classifcation, pathid)

  theme_main = function(){
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #axis.text.x = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      panel.border = element_rect(size = 0.7, linetype = 'solid', colour = "black"),
      legend.position = 'right'
    )
  }
  ggplot(path_mut_classifcation, aes(x = pathid,
                                     y = count,
                                     fill = mut_classification))+
    geom_bar(stat = "identity", width = 0.8, position = "fill")+
    scale_fill_manual(values = my_pal[1:length(unique(path_mut_classifcation$mut_classification))])+
    labs(x="", y="Ration")+theme_bw()+
    scale_y_continuous(position = "right")+
    theme_main()+ggpubr::rotate_x_text(angle = 45)
})


#' Enrichment of somatic mutation in each Hallmark Gene Set or parhway network
#'
#' @param obj PathObject
#' @param my_pal Color table in coding and non-coding
#' @param mut_classifcation The somatic mutation classification data table of the cohort
#' is composed of two columns, the first column is the gene name, and the second column
#' is the mutation division of the gene
#' @param degree Number of genes connected to somatic mutant genes in Hallmark Gene Set Network
#'
#' @export
#'
#'
#' @import ggplot2
setMethod("mut_enrich", "PathObject", function(obj, my_pal, mut_classifcation,
                                               degree=NULL){
  FDR <- pathid <- cluster <- NULL
  stopifnot(class(obj)=="PathObject")
  stopifnot(ncol(mut_classifcation)==2)
  if (is.null(obj@module1$pathway_feature)){
    stop("please compute network features first", call. = FALSE)
  }
  colnames(mut_classifcation) = c('geneid', 'mut_classification')
  pathway_feature = obj@module1$pathway_feature
  num_samples = dim(obj@gene_exp_TPM_obj)[2]

  mut_classifcation_freq = as.data.frame(table(mut_classifcation$geneid))
  mut_classifcation_freq = mut_classifcation_freq[mut_classifcation_freq$Freq>num_samples*0.01, ]
  high_mut_genes = mut_classifcation_freq$Var1
  high_mut_genes = as.character(high_mut_genes)

  mut_enrich_pvalue = lapply(names(pathway_feature), function(x){
    path_feature_tmp = pathway_feature[[x]]
    gene_tmp_innet = path_feature_tmp$node_df
    if (!is.null(degree)){
      gene_tmp_innet = gene_tmp_innet[gene_tmp_innet$degree>=degree, ]
    }

    num_innet_gene = nrow(gene_tmp_innet)

    gene_tmp_innet_mut = high_mut_genes[high_mut_genes%in%gene_tmp_innet$geneid]

    pvalue = stats::phyper(q = length(gene_tmp_innet_mut)-1,
                    m = length(high_mut_genes),
                    n = 25000-length(high_mut_genes),
                    k = num_innet_gene,
                    lower.tail = FALSE)
    return(pvalue)
  })

  names(mut_enrich_pvalue) = names(pathway_feature)
  mut_enrich_pvalue = data.frame(pvalue = unlist(mut_enrich_pvalue),
                                 pathid = names(mut_enrich_pvalue))
  mut_enrich_pvalue$pathid = stringr::str_remove_all(mut_enrich_pvalue$pathid,
                                            "HALLMARK_")
  mut_enrich_pvalue$FDR = stats::p.adjust(as.numeric(mut_enrich_pvalue$pvalue),
                                   method = "fdr")

  mut_enrich_pvalue$cluster = if_else(mut_enrich_pvalue$FDR<=0.1,
                                      "sign", if_else(mut_enrich_pvalue$FDR<=0.3,
                                                      "middle", "non"))

  mut_enrich_pvalue$FDR = -log10(mut_enrich_pvalue$FDR)
  mut_enrich_pvalue = dplyr::arrange(mut_enrich_pvalue, desc(FDR))
  mut_enrich_pvalue$pathid = factor(mut_enrich_pvalue$pathid,
                                    levels = mut_enrich_pvalue$pathid)

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

  ggplot(mut_enrich_pvalue, aes(pathid, FDR))+geom_col(aes(fill=cluster), width = 0.4)+
    scale_fill_manual(values = c(my_pal[2], "gray60", my_pal[1]))+
    ylab("FDR")+xlab("")+scale_y_continuous(expand = c(0, 0))+
    theme_bw()+
    ggpubr::rotate_x_text(angle = 45)+
    geom_hline(yintercept = -log10(0.1),
               linetype="dotdash", color='gray30')+
    theme_main()+theme_legend()
})


#' Describe the small world characteristics of Hallmark Gene Set or parhway network
#'
#' @param obj PathObject
#' @param my_pal color table
#'
#' @export
#'
#'
#' @import ggplot2

setMethod("plot_point_samllworld", "PathObject", function(obj, my_pal){
  alphen <- cluster_coef <- NULL
  nettop_list = obj@module1$pathway_feature
  nettop_feature1 = sapply(nettop_list, function(x){
    x = x$small_world_sig
    return(x)
  }, simplify = TRUE)

  nettop_feature1 = t(nettop_feature1)

  nettop_feature1[, 'parameter1'] = nettop_feature1[, 'parameter1']
  nettop_feature1[nettop_feature1>1] = 1
  nettop_feature1 = as.data.frame(nettop_feature1)
  colnames(nettop_feature1) = c('alphen', 'cluster_coef', 'plaw')
  nettop_feature1 = dplyr::arrange(nettop_feature1, desc(alphen))
  nettop_feature1_TCGA = nettop_feature1
  rownames(nettop_feature1_TCGA) = stringr::str_remove_all(rownames(nettop_feature1_TCGA), "HALLMARK_")


  #my_palatte = colorRampPalette(brewer.pal(9,"YlOrRd"))(nrow(nettop_feature1_TCGA))
  nettop_feature1_TCGA = as.data.frame(nettop_feature1_TCGA)
  colnames(nettop_feature1_TCGA) = c('alphen', 'cluster_coef', 'plaw')
  nettop_feature1_TCGA$pathid = stringr::str_remove(rownames(nettop_feature1_TCGA), 'HALLMARK_')
  ggplot(nettop_feature1_TCGA, aes(alphen, cluster_coef))+geom_point(alpha=0.8, color=my_pal[1])+
    cowplot::theme_cowplot()+theme(legend.title = element_blank(),
                          legend.position = 'none')+
    geom_vline(xintercept = 0.5, linetype='dotdash', color='black')+
    geom_hline(yintercept = 1, linetype='dotdash', color='black')
})

#' Describe the signaficance of  characteristics of Hallmark Gene Set or parhway network
#'
#' @param obj PathObject
#' @param my_pal color table
#'
#' @export
#'
#'
#' @import ggplot2

setMethod("plot_point_fetsig", "PathObject", function(obj, my_pal){
  aplen <- alphen <- cluster_coef <- NULL
  nettop_list = obj@module1$pathway_feature
  nettop_feature2 = sapply(nettop_list, function(x){
    x = x$feature_sig
    return(x)
  }, simplify = TRUE)

  nettop_feature2 = as.data.frame(t(nettop_feature2))
  nettop_feature2[nettop_feature2$aplen>0.5, 'aplen'] = 1-nettop_feature2[nettop_feature2$aplen>0.5, 'aplen']
  nettop_feature2 = dplyr::arrange(nettop_feature2, desc(aplen))
  nettop_feature2_TCGA = nettop_feature2
  rownames(nettop_feature2_TCGA) = stringr::str_remove_all(rownames(nettop_feature2_TCGA), "HALLMARK_")

  nettop_feature2_TCGA = as.data.frame(nettop_feature2_TCGA)
  colnames(nettop_feature2_TCGA) = c('alphen', 'cluster_coef')
  nettop_feature2_TCGA$pathid = stringr::str_remove(rownames(nettop_feature2_TCGA), 'HALLMARK_')
  ggplot(nettop_feature2_TCGA, aes(alphen, cluster_coef))+geom_point(alpha=0.8, color=my_pal[1])+
    cowplot::theme_cowplot()+theme(legend.title = element_blank(),
                          legend.position = 'none')+
    geom_vline(xintercept = 0.05, linetype='dotdash', color='black')+
    geom_hline(yintercept = 0.05, linetype='dotdash', color='black')
})
