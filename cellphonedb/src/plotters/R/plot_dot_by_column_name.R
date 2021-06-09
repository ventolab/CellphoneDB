library(ggplot2)
dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
                    width = 8,
                    height = 10,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){

  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)

  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }

  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

  ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  if (output_extension == '.pdf') {
      ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
      ggsave(filename, width = width, height = height, limitsize=F)
  }
}

