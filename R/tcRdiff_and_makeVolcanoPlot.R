#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title tcrDiff: Performs differential motif analysis
#' @description tcRdiff helps with an in-depth characterization of T-cell repertoires. It uses the output from GLIPH2 to distinguishing clonally expanded motifs. Clonally expanded and enriched motifs are identified by comparing the summed contribution scores of samples from each cluster (summed template frequency) using exact poisson statistics. Enriched motifs' log2-fold change and p-values can be visualized as volcano plots to identify expanded TCRs. The function outputs a data.table that can be visualized using the makeVolcanoPlot function.
#' @export dplyr @export ggplot2 @export ggpubr @export data.table @export stats @export ggrepel
#' @param data data.frame. GLIPH output .csv file
#' @param ident.1 string. Identity class with a sample name of interest.
#' @param ident.2 string. A second identity class for comparison.
#' @param save set save=TRUE to save results to a .csv file
#' @param expansion.score.threshold numeric. To threshold enrichment of clonal expansion within the cluster.
#' @param fisher.score.threshold numeric. A threshold for the p-value obtained by performing the Fisher's exact test with a contingency table containing unique_cdr3_sample, unique_cdr3_ref, the number of remaining sample sequences and the number of remaining reference sequences.The score reports the probability that a random sample of the same size as the sample set but into the reference set (i.e. naive repertoire) would generate an enrichment of the given motif at least as high as has been observed in the sample set.
#' @examples

#'
#' # Read GLIPH output with 'sample_name' column (with only sample information)
#' data <- read.csv("36604540_GLIPH2_Output.csv")
#'
#' # Use function
#' df <- tcRDiff(data, ident.1 = "MtbLys",ident.2 = "D360", save = F, expansion.score.threshold = 0.05, fisher.score.threshold = 0.05)
#'


tcRDiff <- function(
    data,
    ident.1 = NULL,
    ident.2 = NULL,
    save=FALSE,
    expansion.score.threshold = Inf,
    fisher.score.threshold = Inf,
    ...
    ) {

  if(ident.1 %in% unique(data$sample_name)==FALSE|ident.2 %in% unique(data$sample_name)==FALSE){
    tryCatch(stop(cat('ident.1 and/or ident.2 not present in the \"sample_name\" column. Enter the right identity name(s) to proceed.')),
             error=function(cond) {})
  }


  df_filtered <- data[data$sample_name %in% c(ident.1, ident.2), ]

  tryCatch(stat_summary_df <- df_filtered %>% group_by(pattern) %>% mutate(freq_summed=sum(Freq), contribution=(Freq/freq_summed)*100,
                                                                  "{ident.1}":= sum(.data$contribution[grep(ident.1,.data$Sample)])/length(grep(ident.1,.data$Sample)),
                                                                  "{ident.2}":= sum(.data$contribution[grep(ident.2,.data$Sample)])/length(grep(ident.2,.data$Sample)),
                                                                  "{ident.1}_by_{ident.2}":= .data[[ident.1]]/.data[[ident.2]],
                                                                  "log2FC_{ident.1}_by_{ident.2}" := log2(.data[[ident.1]]/.data[[ident.2]]),
                                                                  "summed_{ident.1}_contribution" := round(sum(.data$contribution[grep(ident.1,.data$Sample)])),
                                                                  "summed_{ident.2}_contribution":= round(sum(.data$contribution[grep(ident.2,.data$Sample)]))),
           error=function(cond) {print("Error! Check file format.")})

  stat_summary_df_lfc <- stat_summary_df %>% filter(!(.data[[ident.1]]/.data[[ident.2]]== "NaN"))
  gliph_stats <- NULL

  for(i in unique(stat_summary_df_lfc$pattern)){
    print(paste0("Calculating statistics for: ",i))
    df <- stat_summary_df %>% filter(pattern==i)
    gliph_stats[[i]] <- tryCatch(poisson.test(c(unique(round(sum(df$contribution[grep(ident.1,df$Sample)]))),
                                                unique(round(sum(df$contribution[grep(ident.2,df$Sample)])))),
                                              c(length(grep(ident.1,df$Sample)),length(grep(ident.2,df$Sample)))),
                                 error=function(cond) {print("Error! Check for NA values in columns")})

  }

  tryCatch(gliph_stats_merged <- rbindlist(gliph_stats, idcol = T) %>% distinct(.id, p.value,estimate, parameter) %>% dplyr::rename(pattern=.id),
           error=function(cond) {print("Error! Check file format")})


  gliph_stats_final <- merge(stat_summary_df_lfc, gliph_stats_merged, by="pattern") %>%
    filter(Fisher_score < fisher.score.threshold & expansion_score < expansion.score.threshold) %>%  data.table()

  if(save==TRUE){
    write.csv(gliph_stats_final,file = paste0(ident.1,"_vs_",ident.2,".csv"))
  }

  return(gliph_stats_final)


}



#' @title makeVolcanoPlot: Makes a volcano plot using tcRdiff's output
#' @description Uses the opuput data from the tcRdiff function to make a volcano plot to distinguish between clonally expanded motifs.
#' @param data data.table. Output from tcRdiff function
#' @param ident.1 string. Identity class with a sample name of interest.
#' @param ident.2 string. A second identity class for comparison.
#' @param expansion.score.threshold numeric. To threshold enrichment of clonal expansion within the cluster
#' @param log2FC.threshold numeric. A threshold for the log2FC values calculated by the tcRdiff function
#' @param fisher.score.threshold numeric. A threshold for the p-value obtained by performing the Fisher's exact test with a contingency table containing unique_cdr3_sample, unique_cdr3_ref, the number of remaining sample sequences and the number of remaining reference sequences.The score reports the probability that a random sample of the same size as the sample set but into the reference set (i.e. naive repertoire) would generate an enrichment of the given motif at least as high as has been observed in the sample set.
#' @param limits.x x-axis limits set to c(NA, NA) by default
#' @param colors set to c("red", "blue", alpha("grey85",0.8)) by default
#' @param point.size.range set to c(5,2) by default
#' @examples
#' # example code
#'
#' # Read GLIPH output with 'sample_name' column (with only sample information)
#' data <- read.csv("36604540_GLIPH2_Output.csv")
#'
#' # Calculate statistics
#' df <- tcRDiff(data, ident.1 = "MtbLys",ident.2 = "D360", save = F, expansion.score.threshold = 0.05, fisher.score.threshold = 0.05)
#'
#' # Make volcano plot
#' makeVolcanoPlot(df, ident.1 = "MtbLys",ident.2 = "D360", log2FC.threshold = 0, expansion.score.threshold = 0.05)
#' dev.off()



makeVolcanoPlot <- function(
    data,
    ident.1 = NULL,
    ident.2 = NULL,
    expansion.score.threshold = Inf,
    log2FC.threshold=0,
    fisher.score.threshold = Inf,
    limits.x = c(NA, NA),
    colors = c("red", "blue", alpha("grey85",0.8)),
    point.size.range = c(5,2),
    ...
    ){

  data %>% filter(!(pattern=="single")) %>% select(pattern, grep("log2FC", names(data), value = TRUE), p.value, expansion_score) %>%  distinct() -> data

  names(colors) <- c(paste0("Up in ", ident.1), paste0("Up in ", ident.2), "NA")
  data$diffexpressed <- "NA"
  data$diffexpressed[data[[grep("log2FC", names(data), value = TRUE)]]> log2FC.threshold & data$p.value < 0.05] <- paste0("Up in ", ident.1)
  data$diffexpressed[data[[grep("log2FC", names(data), value = TRUE)]]< -log2FC.threshold & data$p.value < 0.05] <- paste0("Up in ", ident.2)

  data[[grep("log2FC", names(data), value = TRUE)]] <- as.numeric(data[[grep("log2FC", names(data), value = TRUE)]])
  data[[grep("p.value", names(data), value = TRUE)]] <- as.numeric(data[[grep("p.value", names(data), value = TRUE)]])
  data[["diffexpressed"]] <- as.factor(data[["diffexpressed"]])


  data %>% filter(expansion_score<0.05) %>% ggplot(aes(x=!!sym(paste0("log2FC_",ident.1,"_by_",ident.2)), y=-log10(p.value), color=diffexpressed), label=pattern)  + geom_point(position="jitter", aes(size = as.numeric(expansion_score))) +
    theme_minimal()  + geom_text_repel(aes(label=pattern)) + geom_vline(xintercept=c(-log2FC.threshold,log2FC.threshold), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") + scale_color_manual(values = colors) + scale_size(trans = 'reverse') +
    labs(x = "log2 fold change", y = "-log10(p-value)") + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_x_continuous(breaks = scales::pretty_breaks(n =8), limits = limits.x) +
    guides(colour = guide_legend(override.aes = list(size=4))) + scale_size(name="expansion_score", breaks = c(0.001,0.01,0.05),
                                                                            range=point.size.range, limits = c(0.001,0.05), labels = c("p \u2264 0.001", "p \u2264 0.01", "p \u2264 0.05"))


}



