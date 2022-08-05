# plot FACS kinetics from Ahmed
require(readxl)

require(ggplot2)
require(tidyr)
require(dplyr)
require(stringr)
df <- read_xlsx("data/tbet_gata_kinetics_AH.xlsx", sheet = 2)

colnames(df) <- c("id", "cell", "stain", "T-bet", "Gata-3")


df2 <- df %>% pivot_longer(cols = c("T-bet", "Gata-3"))

# separate IDS
df2 <- df2 %>% separate(id, into = c(NA, NA, "time", NA, NA), sep = "_|h")

# clean names

df2$cell[df2$cell == "TH1"] <- "Th1"

df2$cell[df2$cell == "ThMix"] <- "Th1/2"


df3 <- df2 %>% pivot_wider(names_from = stain, values_from = value)

df3$geo_mean <- df3$A / df3$B

df3$time <- as.numeric(df3$time)


# kick out naive cells
df4 <- df3 %>% filter(cell != "naive")
df5 <- df3 %>% filter(cell == "naive")

df3$geo_mean_norm = df3$geo_mean
df3$geo_mean_norm[df3$name == "T-bet"] = df3$geo_mean_norm[df3$name == "T-bet"] / 1.38
df3$geo_mean_norm[df3$name == "Gata-3"] = df3$geo_mean_norm[df3$name == "Gata-3"] / 7.78


g <- ggplot(data = df3, aes(x = time, y = geo_mean_norm, colour = cell))

g2 <- g + geom_point() + geom_line() + facet_wrap(~name) + theme_bw(base_size = 12) + 
  ylab("geo mean index") + xlab("time (h)")
print(g2)
