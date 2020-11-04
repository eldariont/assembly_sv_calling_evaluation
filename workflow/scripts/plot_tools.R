library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("assembly", "caller", "tp_matching_gt", "tp_differing_gt", "fp", "fn"))
res$assembly = factor(res$assembly, levels=c('wenger', 'garg'), labels=c('Asm A', 'Asm B'))
res$caller = factor(res$caller, levels=c('dipcall', 'svim'), labels=c('DipCall', 'SVIM-asm'))


metrics_calls <- res %>%
    mutate(precision = (tp_matching_gt + tp_differing_gt) / (tp_matching_gt + tp_differing_gt + fp)) %>%
    mutate(recall = (tp_matching_gt + tp_differing_gt) / (tp_matching_gt + tp_differing_gt + fn)) %>%
    mutate(f1 = 2*precision*recall/(precision+recall)) %>%
    mutate(mode = 'calls') %>%
    pivot_longer(c("precision", "recall", "f1"), names_to="metric", values_to="value")

metrics_gt <- res %>%
    mutate(precision = (tp_matching_gt) / (tp_matching_gt + tp_differing_gt + fp)) %>%
    mutate(recall = (tp_matching_gt) / (tp_matching_gt + tp_differing_gt + fn)) %>%
    mutate(f1 = 2*precision*recall/(precision+recall)) %>%
    mutate(mode = 'gtcomp') %>%
    pivot_longer(c("precision", "recall", "f1"), names_to="metric", values_to="value")

final <- rbind(metrics_calls, metrics_gt) %>%
    mutate(metric = factor(metric, levels=c('precision', 'recall', 'f1'), labels=c('Precision', 'Recall', 'F1 Score'))) %>%
    mutate(mode = factor(mode, levels=c('calls', 'gtcomp'), labels=c('Calls', 'Calls + Genotypes')))

final%>%
    ggplot(aes(assembly, value, fill=caller)) +
      geom_col(color="black", position=position_dodge()) +
      scale_fill_manual(values=c("gray", "black")) +
      labs(y = "Value", x = "Assembly", fill = "Tool") +
      facet_grid(mode~metric) +
      scale_y_continuous(minor_breaks = seq(0 , 1, 0.1), breaks = seq(0, 1, 0.2)) +
      theme_bw() +
      theme(text = element_text(size=18))
ggsave(args[2], width=8, height=8)

final %>% write_tsv(args[3])
