library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("size", "type", "count", "cumlength", "tool"))
res[res$type == "DUP", "type"] <- "DUP:TANDEM"
res[res$type == "DUP_INT", "type"] <- "DUP:INT"
res[res$type == "DEL/INV", "type"] <- "Complex"
res[res$type == "DUP/INS", "type"] <- "Complex"
res[res$type == "INVDUP", "type"] <- "Complex"
res[res$type == "INV/INVDUP", "type"] <- "Complex"
res[res$type == "cnv", "type"] <- "Complex"
res$size = factor(res$size, levels=c('tiny', 'small', 'medium', 'large', 'huge', 'all'), labels=c('tiny (50-199bp)', 'small (200-999bp)', 'medium (1-9kb)', 'large (10-99kb)', 'huge (100kb-)', 'all'))
res$type = factor(res$type, levels=c('DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP:INT', 'Complex'), labels=c('DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP:INT', 'Complex'))
res$tool = factor(res$tool, levels=c('DipCall', 'SVIM'), labels=c('DipCall', 'SVIM'))

totals <- res %>%
            group_by(tool, size) %>%
            summarize(total = sum(count))
res%>%
    ggplot(aes(tool, count, fill=type)) +
      geom_col() +
      scale_fill_brewer(palette="RdYlBu") +
      labs(y = "Count", x = "Tool", fill = "SV class") +
      facet_wrap(~size, ncol=3) +
      theme_bw() +
      geom_text(aes(tool, total+3000, label = total, fill = NULL), size=2, data = totals) +
      theme(text = element_text(size=8), axis.text.x = element_text(size=7), axis.text.y = element_text(size=7), legend.position="bottom") +
      guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(args[2], width=4, height=4)
res %>% write_tsv(args[3])

