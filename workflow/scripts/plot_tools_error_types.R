library(tidyverse)
library(scales)
library(viridis)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("assembly", "caller", "tp_matching_gt", "tp_differing_gt", "fp", "fn"))
res$assembly = factor(res$assembly, levels=c('wenger', 'garg'), labels=c('Asm A', 'Asm B'))
res$caller = factor(res$caller, levels=c('dipcall', 'svim'), labels=c('DipCall', 'SVIM-asm'))

final <- res %>%
    pivot_longer(cols=c("tp_matching_gt", "tp_differing_gt", "fp", "fn")) %>%
    mutate(name = factor(name, levels=c('fp', 'fn', 'tp_differing_gt', 'tp_matching_gt'), labels=c('FP', 'FN', 'Wrong GT', 'TP')))

final%>%
    ggplot(aes(caller, value, fill=name)) +
      geom_bar(position="stack", stat="identity") +
      scale_fill_manual(values=c("#482677FF", "#33638DFF", "#FDE725FF", "#55C667FF")) +
      labs(y = "SV calls", x = "Tool", fill = "Group") +
      facet_grid(~assembly) +
      theme_bw() +
      theme(text = element_text(size=16))

ggsave(args[2], width=6, height=4)
    
final %>% write_tsv(args[3])
