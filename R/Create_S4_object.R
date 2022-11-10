#' Create a PathObj S4 Object
#'
#' @param gmt_path a gmt type file path It is recommended to use the 50Hallmark gene set
#' @param gene_exp_counts a data.frame or matrix gene exprssion or path the first column should be geneid
#' @param gene_map a gtf type file or data.frame with 5 columns included gene length start end chrid and geneid
#' @param gene_exp_TPM same as gene_exp_counts
#' @param ppi_data_path two columns node A and node B
#' @param cor_method 1 pearson 2 spearman 3 kendall
#' @param mut_exp MAF or binary matrix
#' @param count_log if log2 Transform
#' @param TPM_log if log2 Transform
#' @param GeneSetCollect gmt transform to GeneSetCollect
#' @param sep Tabs
#'
#' @return PathObj S4 object
#' @import dplyr methods data.table SummarizedExperiment GenomicRanges tidyfst stringr purrr tibble IRanges S4Vectors
#' @importClassesFrom GenomicRanges GRanges
#' @export

CreatePathObject = function(gmt_path=NULL, gene_exp_counts=NULL,
                            gene_map=NULL,
                            gene_exp_TPM=NULL,
                            ppi_data_path=NULL,
                            cor_method=1,
                            mut_exp=NULL,
                            count_log=TRUE,
                            TPM_log=FALSE,
                            GeneSetCollect=FALSE,
                            sep = "\t"){


  cor_method = switch(cor_method, "pearson", "kendall", "spearman")
  if (is.null(cor_method)){
    stop("cor method must be one of pearson, kendall, spearman")
  }
  if (is.null(gmt_path)){
    stop('Gene set or pathway file is empty', call. = FALSE)
  }
  lines = strsplit(readLines(gmt_path), sep)
  if (any(sapply(lines, length)<2)){
    txt = paste("all records in the GMT file must have >= 2 fields",
                "\n first invalid line: %s\n", collapse = "")
    .stopf(txt, lines[sapply(lines, length) < 2][[1]])
  }
  describle_lines = sapply(lines, function(x)as.character(x[2]), simplify = TRUE)
  if (any(str_detect(describle_lines, 'http://www.')==FALSE)){
    cat('No gene set interpretation row detected or empty interpretation row')
  }
  dups = new.env(parent = emptyenv())
  lines = lapply(lines, function(elt, dups){
    if (any(d<-duplicated(elt[-c(1:2)]))){
      dups[[elt[1]]] = unique(elt[-c(1:2)][d])
      elt = c(eltp[1:2], unique(elt[-c(1:2)]))
    }
    elt
  }, dups)

  if (length(dups)>0){
    cat('A single gene set or signal pathway has duplicate genes')
  }

  if (isTRUE(GeneSetCollect)){
    template = GeneSet(geneIdType=NullIdentifier(),
                       collectionType=NullCollection())
    GeneSetCollection(lapply(lines, function(x){
      initialize(template, geneIds=unlist(x[-c(1:2)]),
                 setName=x[1], shortDescription=x[2])
    }))
  }else {
    geneset_names = sapply(lines, function(x)as.character(x[1]), simplify = TRUE)
    geneset = lapply(lines, function(x)as.character(x[-c(1:2)]))
    names(geneset) = geneset_names
  }
  names(geneset_names) = as.vector(unlist(lapply(geneset, function(x){Reduce("paste", x)})))
  if (length(unique(names(geneset_names)))<length(names(geneset_names))){
    warning("There is a consistent gene set in which all genes are the same", call. = FALSE)
  }

  if (is.null(gene_exp_TPM)){
    cat("The TPM expression matrix is not provided \n  which will be calculated according to the count matrix")
  }else {
    if (is.character(gene_exp_TPM)){
      if (file.exists(gene_exp_TPM)){
        gene_exp_TPM = fread(gene_exp_TPM, quote = FALSE,
                             header = TRUE, data.table = FALSE)
        colnames(gene_exp_TPM)[1] = 'geneid'
      }else {stop('Invalid expression matrix input path', call. = FALSE)}
    }else if (attr(gene_exp_TPM, 'class')%in%c('tbl_df', 'tbl', 'data.frame', 'matrix', 'array')){
      gene_exp_TPM = as.data.frame(gene_exp_TPM)
      colnames(gene_exp_TPM)[1] = 'geneid'
    }else {stop("invalid structure")}
  }

  if (!is.null(gene_exp_counts)){
    if (is.character(gene_exp_counts)){
      if (file.exists(gene_exp_counts)){
        gene_exp_counts = fread(gene_exp_counts, quote = FALSE,
                                header = TRUE, data.table = FALSE)
        colnames(gene_exp_counts)[1] = 'geneid'
      }else {stop('Invalid expression matrix input path', call. = FALSE)}
    }else if (attr(gene_exp_counts, 'class')%in%c('tbl_df', 'tbl', 'data.frame', 'matrix', 'array')){
      gene_exp_counts = as.data.frame(gene_exp_counts)
      colnames(gene_exp_counts)[1] = 'geneid'
    }else {stop('invalid expression matrix', call. = FALSE)}
  }else {
    cat('no gene exp count supply')
    if (is.null(gene_exp_TPM)){
      stop('gene exp count and  gene exp TPM must supply one', call. = FALSE)
    }
    gene_exp_counts = gene_exp_TPM
  }


  if (is.character(gene_map)){
    if (file.exists(gene_map)){
      if (str_detect(gene_map, "gtf")){
        gene_map = rtracklayer::import(gene_map)
        gene_map = as.data.frame(gene_map)%>%distinct(gene_name, .keep_all = TRUE)%>%
          dplyr::select(id=gene_id, gene=gene_name, chr=seqnames, start=start, end=end)
      }else {
        gene_map = fread(gene_map, sep = '\t', quote = FALSE,
                         header = TRUE, data.table = FALSE)
        gene_map = gene_map[, 1:5]
      }
    }else {stop('Invalid gene map input path', call. = FALSE)}
  }else if (attr(gene_map, 'class')%in%c('tbl_df', 'tbl', 'data.frame', 'matrix', 'array') &&
            ncol(gene_map)==6){
    gene_map = as.data.frame(gene_map)
    gene_map = gene_map[, 1:5]
  }else if (is.null(gene_map)){
    cat('Confirm that the expression matrix ID has been converted to symbol form')
  }else {stop('invalid gene map structure', call. = FALSE)}

  if (is.character(mut_exp)){
    if (file.exists(mut_exp)){
      mut_exp = fread(mut_exp, sep = '\t', quote = FALSE,
                      header = TRUE, data.table = FALSE)
      mut_exp = mut_exp[, c(1, 2, 9)]
      colnames(mut_exp) = c('Sample_ID', 'gene', 'effect')
      mut_exp = mut_exp[!mut_exp$effect%in%'synonymous_variant', ]
      mut_exp = mut_exp[, 1:2]
      mut_exp = as.data.frame(table(mut_exp$Sample_ID, mut_exp$gene))%>%
        dplyr::filter(Freq>0)%>%dplyr::select(Sample_ID=Var1, gene=Var2)
      mut_exp$Sample_ID = as.character(mut_exp$Sample_ID)
      mut_exp$gene = as.character(mut_exp$gene)
    }else {stop('Invalid gene map input path', call. = FALSE)}
  }else if (attr(mut_exp, 'class')%in%c('tbl_df', 'tbl', 'data.frame', 'matrix', 'array') &&
            ncol(mut_exp)==2){
    mut_exp = as.data.frame(mut_exp)
    #mut_exp = mut_exp[, 1:3]
    colnames(mut_exp) = c('Sample_ID', 'gene')
    mut_exp = as.data.frame(table(mut_exp$Sample_ID, mut_exp$gene))%>%
      dplyr::filter(Freq>0)%>%dplyr::select(Sample_ID=Var1, gene=Var2)
    mut_exp$Sample_ID = as.character(mut_exp$Sample_ID)
    mut_exp$gene = as.character(mut_exp$gene)

  }else {stop('invalid mut exp structure', call. = FALSE)}


  colnames(gene_map) = c('id', 'gene', 'chr', 'start', 'end')
  gene_map$start = as.numeric(gene_map$start)
  gene_map$end = as.numeric(gene_map$end)

  if (!str_detect(gene_map$chr, "(C|c)hr")){
    gene_map$chr = str_c("chr", gene_map$chr, sep = '')
  }

  if (!is.null(gene_map))
    com_ids = intersect(as.character(gene_exp_counts[, 1]), gene_map$id)
  if (!is_empty(com_ids)){
    gene_exp_counts = gene_exp_counts[gene_exp_counts$geneid%in%com_ids, ]
    gene_map = gene_map[gene_map$id%in%com_ids, ]

    gene_id = gene_map$gene
    names(gene_id) = gene_map$id
    samples = colnames(gene_exp_counts)[2:ncol(gene_exp_counts)]

    gene_exp_counts = gene_exp_counts%>%left_join(gene_map, by = c('geneid'='id'))%>%distinct(gene, .keep_all = TRUE)%>%
      dplyr::select(gene, samples)%>%filter(!is.na(gene))%>%column_to_rownames(var = 'gene')

    gene_map = gene_map%>%distinct(gene, .keep_all = TRUE)%>%
      column_to_rownames(var = 'gene')
    com_gene = intersect(rownames(gene_exp_counts), rownames(gene_map))

    gene_exp_counts = gene_exp_counts[com_gene, ]
    gene_map = gene_map[com_gene, ]

    gene_map$gene_len = gene_map$end - gene_map$start
  }else {
    com_ids = intersect(as.character(gene_exp_counts[, 1]), gene_map$gene)
    if (!is_empty(com_ids)){
      gene_exp_counts = gene_exp_counts[gene_exp_counts$geneid%in%com_ids, ]
      gene_map = gene_map[gene_map$gene%in%com_ids, ]

      gene_id = gene_map$gene
      names(gene_id) = gene_map$id
      samples = colnames(gene_exp_counts)[2:ncol(gene_exp_counts)]

      gene_exp_counts = gene_exp_counts%>%left_join(gene_map, by = c('geneid'='gene'))%>%distinct(geneid, .keep_all = TRUE)%>%
        dplyr::select(geneid, samples)%>%filter(!is.na(geneid))%>%column_to_rownames(var = 'geneid')%>%na.omit()

      gene_map = gene_map%>%distinct(gene, .keep_all = TRUE)%>%
        column_to_rownames(var = 'gene')
      com_gene = intersect(rownames(gene_exp_counts), rownames(gene_map))

      gene_exp_counts = gene_exp_counts[com_gene, ]
      gene_map = gene_map[com_gene, ]

      gene_map$gene_len = gene_map$end - gene_map$start
    }else {stop('Gene transformation ID mismatch', call. = FALSE)}
  }

  if (count_log){
    gene_exp_counts = 2^gene_exp_counts-1
  }

  if (is.null(gene_exp_TPM)){
    gene_exp_TPM = apply(gene_exp_counts, 2, .count_to_TPM, effLen=gene_map$gene_len)
  }else {
    gene_exp_TPM = column_to_rownames(gene_exp_TPM, var = 'geneid')
  }

  com_gene = Reduce(intersect, list(rownames(gene_exp_counts),
                                    rownames(gene_exp_TPM),
                                    rownames(gene_map)))
  gene_exp_counts = gene_exp_counts[com_gene, ]
  gene_exp_TPM = gene_exp_TPM[com_gene, ]
  gene_map = gene_map[com_gene, ]

  if (isTRUE(TPM_log)){
    gene_exp_TPM = log2(gene_exp_TPM+1)
  }


  rowRanges = GRanges(seqnames = gene_map$chr,
                      ranges = IRanges(gene_map$start, width = gene_map$gene_len),
                      strand = "*",
                      SYMBOL=rownames(gene_map),
                      ENSEMBL=gene_map$id)
  gene_exp_counts_obj = SummarizedExperiment(assays = SimpleList(counts=as.matrix(gene_exp_counts)),
                                             rowRanges = rowRanges)
  gene_exp_TPM_obj = SummarizedExperiment(assays = SimpleList(counts=as.matrix(gene_exp_TPM)),
                                          rowRanges = rowRanges)

  unique_gene = unique(as.vector(unlist(geneset)))

  com_gene = intersect(rownames(gene_exp_TPM), unique_gene)
  gene_exp_TPM_cor = gene_exp_TPM[com_gene, ]

  ppi_dataset = fread(ppi_data_path,header = TRUE, quote = FALSE,
                      sep = '\t', data.table = FALSE)
  ppi_dataset = ppi_dataset[ppi_dataset$gene_A%in%com_gene & ppi_dataset$gene_B%in%com_gene, ]

  cor_test = rlang::call2("cor", x=expr(t(gene_exp_TPM_cor)), method=expr(cor_method))

  ppi_exp_cor = eval(cor_test)

  sample_mut_list = split(mut_exp$gene, mut_exp$Sample_ID)

  module1 = list(ppi_exp_cor=ppi_exp_cor,
                 sample_mut_list=sample_mut_list,
                 geneset=geneset,
                 ppi_dataset=ppi_dataset,
                 geneset_names=geneset_names)

  PathObjectnew = new("PathObject", gene_exp_counts_obj=gene_exp_counts_obj,
                      gene_exp_TPM_obj=gene_exp_TPM_obj,
                      module1=module1)
  return(PathObjectnew)
}

utils::globalVariables(c("eltp", "GeneSet", "NullIdentifier", "NullCollection",
                         "GeneSetCollection", ".stopf", "gene_name",
                         "Freq", "Var1", "Var2", "gene", "geneid"))



