library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("caller", "truvari", "mode", "metric", "value"))
res$caller = factor(res$caller, levels=c('dipcall', 'svim'), labels=c('Dipcall', 'SVIM'))
res$mode = factor(res$mode, levels=c('calls', 'gtcomp'), labels=c('call', 'genotype'))

res %>%
    filter(truvari == "truvari") %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0) %>%
    mutate(precision = 100*precision, recall = 100*recall) %>%
    mutate(f1 = 2*precision*recall/(precision+recall)) %>%
    pivot_longer(c("precision", "recall", "f1"), names_to="metric", values_to="value") %>%
    mutate(metric = factor(metric, levels=c('precision', 'recall', 'f1'), labels=c('Precision', 'Recall', 'F1 Score'))) %>%
    ggplot(aes(mode, value, fill=caller, linetype=mode)) +
      geom_col(color="black", position=position_dodge()) +
      scale_fill_manual(values=c("deepskyblue3", "firebrick2")) +
      labs(y = "Value", x = "Evaluation mode", fill = "Tool", linetype = "Evaluation mode") +
      facet_wrap(~metric) +
      scale_y_continuous(minor_breaks = seq(0 , 100, 5), breaks = seq(0, 100, 10)) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9))

ggsave(args[2], width=8, height=4)
