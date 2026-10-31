#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig2_E5_mCHIRA_clusters/plot_mCHIRA_PRISM.R

library(hdf5r)
library(ggseqlogo)
library(ggplot2)
library(egg)
library(gridExtra)

PRISM_H5 <- "input_files/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_PRISM.h5"
CLUSTER_H5 <- "output/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_clusters.h5"
OUT_DIR <- "output/Fig2_E5_mCHIRA_clusters/"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

prism_h5 <- H5File$new(PRISM_H5,"r")
cluster_h5 <- H5File$new(CLUSTER_H5,"r")
read_ids <- cluster_h5[["pred/CRE_6005/read_ids"]]$read()
labels <- cluster_h5[["pred/CRE_6005/cluster_labels"]]$read()
new_cluster <- c("2"="A","1"="B","3"="C","0"="D","4"="E","5"="F","6"="G","7"="H")
bases <- c("A","C","G","T")
clusters <- sort(unique(labels))

get_delta <- function(i) prism_h5[[paste0("read_",i,"/deltas")]]$read()

ref <- apply(abs(get_delta(0)),2,which.min)

logo_list <- list(); summary_list <- list()
n_reads <- setNames(integer(length(clusters)),clusters)

for(cl in clusters){
  ids <- read_ids[labels==cl]
  avg <- -apply(simplify2array(lapply(ids,get_delta)),c(1,2),mean)

  pos_effect <- sapply(seq_len(ncol(avg)), \(j) mean(avg[-ref[j],j]))
  logo <- matrix(0,4,ncol(avg)); logo[cbind(ref,seq_len(ncol(avg)))] <- pos_effect
  rownames(logo) <- bases
  name <- unname(new_cluster[as.character(cl)])

  logo_list[[name]] <- logo
  n_reads[as.character(cl)] <- length(ids)
  summary_list[[name]] <- data.frame(cluster=name, position=seq_len(ncol(logo)), base=bases[ref], signal=pos_effect)
}

summary_df <- do.call(rbind,summary_list)
rownames(summary_df) <- NULL

write.table(summary_df,file.path(OUT_DIR,"CRE_6005_cluster_PRISM_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

ylims <- range(unlist(logo_list),na.rm=TRUE)

for(cl in clusters){
  name <- unname(new_cluster[as.character(cl)])
  p <- ggseqlogo(logo_list[[name]],method="custom",seq_type="dna") +
    labs(title=paste0("CRE_6005 cluster ", name," (n=",n_reads[as.character(cl)],")"), x="Position region", y="ISM score") +
    scale_x_continuous(breaks=c(1,125,250)) +
    coord_cartesian(ylim=ylims) +
    theme_classic() +
    theme(legend.text=element_text(size=6), panel.border=element_rect(color="black",fill=NA,linewidth=.5),
                      axis.line=element_blank(), plot.title=element_text(size=6), axis.text=element_text(size=6),
                      title=element_text(size=6), axis.title=element_text(size=6), strip.text=element_text(size=6),
                      strip.background=element_blank(), plot.margin=grid::unit(c(0,0,0,0),"cm"))

  pdf(file.path(OUT_DIR,paste0("CRE_6005_cluster_",name,"_WT_PRISM_logo.pdf")),
      width=20/2.54,height=20/2.54)

  grid.arrange(set_panel_size(p,height=grid::unit(2,"cm"),width=grid::unit(4,"cm")))
  dev.off()
}

prism_h5$close_all()
cluster_h5$close_all()
