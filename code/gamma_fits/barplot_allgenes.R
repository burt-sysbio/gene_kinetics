require(dplyr)
require(tidyr)
require(ggplot2)
require(readr)
require(RColorBrewer)

df <- read.csv("output/gamma_fits/feb2021/delays_allstudies.csv")

g <- ggplot(df, aes(x = ynorm, y = reorder(gene, ynorm)))

g1 <- g + geom_bar(aes(fill = category), stat = "identity") +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c("grey", "blue", "red")) +
  xlab("delayed (% of studies)") +
  ylab("gene")

print(g1)

