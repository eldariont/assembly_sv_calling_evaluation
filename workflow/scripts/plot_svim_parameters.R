library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("parameters", "assembly", "caller", "tp_matching_gt", "tp_differing_gt", "fp", "fn"))
res$assembly = factor(res$assembly, levels=c('wenger', 'garg'), labels=c('Asm A', 'Asm B'))
res$caller = factor(res$caller, levels=c('dipcall', 'svim'), labels=c('Dipcall', 'SVIM-asm'))

metrics_calls <- res %>%
    mutate(precision = (tp_matching_gt + tp_differing_gt) / (tp_matching_gt + tp_differing_gt + fp)) %>%
    mutate(recall = (tp_matching_gt + tp_differing_gt) / (tp_matching_gt + tp_differing_gt + fn)) %>%
    mutate(f1 = 2*precision*recall/(precision+recall)) %>%
    mutate(mode = 'calls')

metrics_gt <- res %>%
    mutate(precision = (tp_matching_gt) / (tp_matching_gt + tp_differing_gt + fp)) %>%
    mutate(recall = (tp_matching_gt) / (tp_matching_gt + tp_differing_gt + fn)) %>%
    mutate(f1 = 2*precision*recall/(precision+recall)) %>%
    mutate(mode = 'gtcomp')

final <- rbind(metrics_calls, metrics_gt) %>%
    mutate(mode = factor(mode, levels=c('calls', 'gtcomp'), labels=c('Calls', 'Calls + Genotypes')))

final %>%
    ggplot(aes(recall, precision, color=parameters)) +
      geom_point(size=1.0) +
      scale_color_manual(values=c("deepskyblue3", "goldenrod2", "firebrick2")) +
      labs(y = "Precision", x = "Recall", fill = "Parameters") +
lims(x=c(0.5,1), y=c(0.5,1)) +
      facet_grid(assembly~mode) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9))

ggsave(args[2], width=8, height=4)
