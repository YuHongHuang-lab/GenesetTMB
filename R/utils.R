.count_to_TPM = function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}


.pathway_degreemut = function(pathway_tmp, ppi_exp_cor,
                              sample_mut_list, geneset_names=geneset_names,
                              ppi_dataset=ppi_dataset){
  pathway_tmp_name = geneset_names[Reduce('paste', pathway_tmp)]
  pathway_tmp = pathway_tmp[pathway_tmp%in%colnames(ppi_exp_cor)]
  if (is_empty(pathway_tmp)){
    cat("The gene set does not contain any intersection genes")
    gene_netfeature_tmp = data.frame(net_all=rep(1, length(pathway_tmp)), row.names = pathway_tmp)
  }else {
    ppi_dataset_tmp = ppi_dataset[(ppi_dataset$gene_A%in%pathway_tmp & ppi_dataset$gene_B%in%pathway_tmp), ]
    if (nrow(ppi_dataset_tmp)==0){
      sprintf("The protein interaction of %s was not recorded in the protein database",
              pathway_tmp_name)
      gene_netfeature_tmp = data.frame(net_all=rep(1, length(pathway_tmp)), row.names = pathway_tmp)
    }else {
      ppi_dataset_tmp$cor = apply(ppi_dataset_tmp, 1, function(x){ppi_exp_cor[x[1], x[2]]})
      ppi_dataset_tmp = ppi_dataset_tmp[ppi_dataset_tmp$cor!=1, ]
      ppi_dataset_tmp = distinct(ppi_dataset_tmp, cor, .keep_all = TRUE)
      ppi_dataset_tmp = ppi_dataset_tmp[abs(ppi_dataset_tmp$cor)>=0.4, ]

      if (nrow(ppi_dataset_tmp)!=0){
        ppi_dataset_tmp = dplyr::select(ppi_dataset_tmp, from='gene_A', to='gene_B')

        ppi_dataset_tmp_relation = igraph::graph_from_data_frame(ppi_dataset_tmp,
                                                         directed = FALSE)
        gene_degree = igraph::degree(ppi_dataset_tmp_relation)
        gene_betweenness = igraph::betweenness(ppi_dataset_tmp_relation)
        gene_closeness = igraph::closeness(ppi_dataset_tmp_relation)

        gene_netfeature = data.frame(degree=gene_degree,
                                     betweenness=gene_betweenness,
                                     closeness=gene_closeness,
                                     row.names = names(gene_degree))
        gene_netfeature_tmp = apply(gene_netfeature, 2, function(x){tmp=scale(x);tmp-min(tmp)})
        gene_netfeature_tmp[is.nan(gene_netfeature_tmp)] = 0
        rownames(gene_netfeature_tmp) = rownames(gene_netfeature)
        gene_netfeature_tmp = as.data.frame(gene_netfeature_tmp)
        gene_netfeature_tmp$net_all = gene_netfeature_tmp[, 1] + gene_netfeature_tmp[, 2] +
          gene_netfeature_tmp[, 3]
        gene_netfeature_tmp[pathway_tmp[!pathway_tmp%in%rownames(gene_netfeature_tmp)], 'net_all'] = min(gene_netfeature_tmp$net_all)
        gene_netfeature_tmp = gene_netfeature_tmp[, 'net_all', drop=FALSE]
        gene_netfeature_tmp$net_all = gene_netfeature_tmp$net_all+1
      }else{
        gene_netfeature_tmp = data.frame(net_all=rep(1, length(pathway_tmp)), row.names = pathway_tmp)
      }
    }
  }

  sapply(sample_mut_list, function(x){
    if (!is_empty(intersect(x, rownames(gene_netfeature_tmp)))){
      com_genes = intersect(x, rownames(gene_netfeature_tmp))
      result = sum(gene_netfeature_tmp[com_genes, 'net_all']) / length(pathway_tmp)
      return(result)
    }else return(0)}, simplify = TRUE)
}



.sample_Cols = function(d, pSample=NULL, weightItems=NULL,
                        pRow=NULL, weightFeatures=NULL){
  space = ncol(d)
  sampleN = floor(space*pSample)
  sampleCols = sort(sample(space, sampleN, replace = FALSE,
                           prob = weightItems))
  this_samples <- sampleRows <- NA
  if ((!is.null(pRow)) && (pRow<1) || (!is.null(weightFeatures))){
    space = nrow(d)
    sampleN = floor(space*pRow)
    sampleRows = sort(sample(space, sampleN, replace = FALSE,
                             prob = weightFeatures))
    this_samples = d[sampleRows, sampleCols]
    dimnames(this_samples) = NULL
  }


  return(list(submat=this_samples, subrows=sampleRows, subcols=sampleCols))
}



.connect_mat = function(clusterAssignments, mCount, sampleKey){
  names(clusterAssignments) = sampleKey
  cls = lapply(unique(clusterAssignments), function(i)as.numeric(names(clusterAssignments[clusterAssignments%in%
                                                                                            i])))
  for (i in 1:length(cls)){
    nelts = 1:ncol(mCount)
    cl = as.numeric(nelts%in%cls[[i]])
    updt = outer(cl, cl)
    mCount = mCount + updt
  }
  return(mCount)
}

.triangle = function(m, mode){
  n = dim(m)[1]
  m_3 = m
  m_2 = m
  m_2[upper.tri(m_2)] = NA
  diag(m_2) = NA
  m_1 = m[lower.tri(m)]

  if (mode==3){
    return(m_3)
  }else if (mode==2){
    return(m_2)
  }else if (mode==1){
    return(m_1)
  }
}

.graph_obj_construct = function(pathway_tmp, geneset_names,
                                ppi_dataset, ppi_exp_cor){
  pathway_tmp_name = geneset_names[Reduce('paste', pathway_tmp)]
  pathway_tmp = pathway_tmp[pathway_tmp%in%colnames(ppi_exp_cor)]
  if (is_empty(pathway_tmp)){
    cat("The gene set does not contain any intersection genes")
    pathway_graph_tmp = NULL
  }else {
    ppi_dataset_tmp = ppi_dataset[(ppi_dataset$gene_A%in%pathway_tmp & ppi_dataset$gene_B%in%pathway_tmp), ]
    if (nrow(ppi_dataset_tmp)==0){
      sprintf("The protein interaction of %s was not recorded in the protein database",
              pathway_tmp_name)
      pathway_graph_tmp = NULL
    }else {
      ppi_dataset_tmp$cor = apply(ppi_dataset_tmp, 1, function(x){ppi_exp_cor[x[1], x[2]]})
      ppi_dataset_tmp = ppi_dataset_tmp[ppi_dataset_tmp$cor!=1, ]
      ppi_dataset_tmp = distinct(ppi_dataset_tmp, cor, .keep_all = TRUE)
      ppi_dataset_tmp = ppi_dataset_tmp[abs(ppi_dataset_tmp$cor)>=0.4, ]
      ppi_dataset_tmp = ppi_dataset_tmp%>%na.omit()

      if (nrow(ppi_dataset_tmp)!=0){
        ppi_dataset_tmp = dplyr::select(ppi_dataset_tmp, from='gene_A', to='gene_B')

        pathway_graph_tmp = graph_from_data_frame(ppi_dataset_tmp,
                                                  directed = FALSE)
      }else{
        pathway_graph_tmp = NULL
      }
    }
  }
  cat("Net summary of ", pathway_tmp_name, "is done", sep = " ")
  return(pathway_graph_tmp)
}

#object = pathway_obj_list[[21]]
#' Calculate the network properties of each Hallmark Gene Set network only in windows
#'
#' @param object PathObj
#' @param nthreads Number of threads
#'
#' @export
#'
#' @import igraph
.network_summary = function(object, nthreads=4){
  if (attr(object, "class")!="igraph"){
    stop("invalid net object", call. = FALSE)
  }

  object = igraph::delete.vertices(object, names(igraph::degree(object)[igraph::degree(object)==0]))
  object = igraph::simplify(object)

  igraph::V(object)$degree = igraph::degree(object)
  #hist(V(object)$degree)
  pvalue_plaw = .Plaw_pvalue(x = igraph::V(object)$degree, nthreads = nthreads)

  igraph::V(object)$ndegree = graph.knn(object, igraph::V(object))$knn
  igraph::V(object)$betweenness_centrality = igraph::betweenness(object)
  igraph::V(object)$eigenvector_centrality = igraph::evcent(object)$vector
  node_df = data.frame(geneid = names(igraph::V(object)),
                       degree = igraph::V(object)$degree,
                       ndegree = igraph::V(object)$ndegree,
                       betweenness_centrality = igraph::V(object)$betweenness_centrality,
                       eigenvector_centrality = igraph::V(object)$eigenvector_centrality)
  node_df[is.nan(node_df$ndegree), "ndegree"] = 0

  num_node = length(igraph::V(object))
  num_edge = length(igraph::E(object))
  aplen = igraph::average.path.length(object)
  diameter = igraph::diameter(object)
  density = igraph::graph.density(object)
  cluster_coef = igraph::transitivity(object)

  net_list = list(num_node = num_node,
                  num_edge = num_edge,
                  aplen = aplen,
                  diameter = diameter,
                  density = density,
                  cluster_coef = cluster_coef)
  deg = node_df$degree
  aplen_boot = vector()
  cluster_coef_boot = vector()

  for (i in 1:1000){
    g_boot = igraph::degree.sequence.game(deg, method="simple")
    aplen_boot = c(aplen_boot, average.path.length(g_boot))
    cluster_coef_boot = c(cluster_coef_boot, transitivity(g_boot))
  }

  aplen_small_world = vector()
  cluster_coef_small_world = vector()

  for (i in 1:1000){
    g_small_world = erdos.renyi.game(net_list$num_node, net_list$num_edge, type="gnm",
                                     directed = FALSE)
    aplen_small_world = c(aplen_small_world, average.path.length(g_small_world))
    cluster_coef_small_world = c(cluster_coef_small_world, transitivity(g_small_world))
  }

  aplen_boot[is.na(aplen_boot)] = 0
  cluster_coef_boot[is.na(cluster_coef_boot)] = 0
  aplen_small_world[is.na(aplen_small_world)] = 0
  cluster_coef_small_world[is.na(cluster_coef_small_world)] = 0

  #feature significance
  pvalue_aplen = sum(aplen_boot>net_list$aplen)/length(aplen_boot)
  pvalue_cluster_coef = sum(cluster_coef_boot>net_list$cluster_coef)/length(cluster_coef_boot)

  #small world
  pvalue_aplen = pvalue_aplen
  smallworld_parameter1 = sum(aplen_small_world>=net_list$aplen)/sum(aplen_small_world<=net_list$aplen)
  smallworld_parameter1 = .sigmoid(smallworld_parameter1-1)


  if (net_list$cluster_coef>mean(cluster_coef_small_world)){
    smallworld_parameter2 = sum(cluster_coef_small_world<=net_list$cluster_coef)/sum(cluster_coef_small_world>=net_list$cluster_coef)
    smallworld_parameter2 = .sigmoid(smallworld_parameter2-1)
  }else {
    smallworld_parameter2 = sum(cluster_coef_small_world>=net_list$cluster_coef)/sum(cluster_coef_small_world<=net_list$cluster_coef)
    smallworld_parameter2 = .sigmoid(smallworld_parameter2-1)
  }
  feature_sig = c(pvalue_aplen, pvalue_cluster_coef)
  names(feature_sig) = c("aplen", "cluster_coef")

  small_world_sig = c(smallworld_parameter1, smallworld_parameter2, pvalue_plaw)
  names(small_world_sig) = c("parameter1", "parameter2", "plaw")

  return(list(node_df=node_df,
              net_list=net_list,
              feature_sig=feature_sig,
              small_world_sig=small_world_sig))
}

#' Functions for calculating power law distribution are not open to the public
#'
#' @param x Degree
#' @param nthreads number of nthreads
#'
#' @import poweRlaw
.Plaw_pvalue = function(x, nthreads){
  x = x[x!=0]
  m_pl = poweRlaw::displ$new(x)
  est = poweRlaw::estimate_xmin(m_pl)
  m_pl$setXmin(est)
  bs_p = poweRlaw::bootstrap_p(m_pl, threads = nthreads)$p
  return(bs_p)
}

.sigmoid = function(x){1/(1+exp(-x))}



