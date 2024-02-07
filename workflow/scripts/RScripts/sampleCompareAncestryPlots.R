library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
user.input.1 <- args[1]
user.input.2 <- args[2]
user.input.3 <- args[3]
user.input.4 <- args[4]

t <- read.table(user.input.1, sep = "\t", header = T, comment.char = "")
pairs <- read.table(user.input.2, sep = "\t", header = T)
pairs <- pairs %>% mutate(key = paste0(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = ""))
pairs <- pairs[duplicated(pairs[, "key"]), ]
samples <- unique(c(pairs$Sample1, pairs$Sample2))

mapping <- list("EUR" = "European", "EAS" = "East Asian", "AMR" = "American", "SAS" = "South Asian", "AFR" = "African")
t$predAncestry <- unlist(mapping[t$predicted_ancestry])
t$color <- t$predAncestry
t$color[t$X.sample_id %in% samples] <- "UserSamples"
p <- plot_ly(t,
  x = ~PC1, y = ~PC2, color = as.factor(t$color), colors = c("#0C4B8E", "#FF0000", "#f1a340", "#43a2ca", "#8856a7", "grey"),
  hoverinfo = "text", text = ~ paste("</br> Sample Id:", t$X.sample_id, "</br> Ancestry:", t$predAncestry), type = "scatter", mode = "markers"
) %>%
  # add_trace(marker = list(size = 12)) %>%
  layout(scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2")))
htmlwidgets::saveWidget(p, user.input.3)

samplesAncestry <- t[t$X.sample_id %in% samples, c(1, 4:8)]
mapping <- list("EUR_prob" = "European", "EAS_prob" = "East Asian", "AMR_prob" = "American", "SAS_prob" = "South Asian", "AFR_prob" = "African")
colnames(samplesAncestry) <- c(colnames(samplesAncestry)[1], unlist(mapping[colnames(samplesAncestry)[2:6]]))
d <- data.frame(rbind(
  cbind(pairs$Sample1, c(rep("Sample1", length(pairs$Sample1))), paste0(pairs$Sample1, "\nvs\n", pairs$Sample2)),
  cbind(pairs$Sample2, c(rep("Sample2", length(pairs$Sample2))), paste0(pairs$Sample1, "\nvs\n", pairs$Sample2))
))
mDF <- merge(x = d, y = samplesAncestry, by.x = "X1", by.y = "X.sample_id", all.x = TRUE, all.y = FALSE)
gData <- mDF %>% tidyr::pivot_longer(c(4:8), names_to = "Ancestry", values_to = "Somalier.Score")
g <- ggplot(gData, aes(fill = Ancestry, y = Somalier.Score, x = X2)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = 5)) +
  facet_wrap(~X3)
ggsave(user.input.4, g)
