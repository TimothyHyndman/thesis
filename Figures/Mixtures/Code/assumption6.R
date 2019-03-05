# Plot graphs of cauchy and normal derivatives showing that they satisfy assumption 6 in mixture chapter, 2d result.
#setwd("C:/Users/timot/Git/thesis/Figures/Mixtures/Code")
setwd("C:/Users/Timothy/Git/thesis/Figures/Mixtures/Code")
library(tidyverse)

# Normal
x <- seq(-1, 4, length.out = 1000)
y1 <- -x * exp(-x^2/2)
y2 <- (x - 2) * exp(-(x - 2)^2/2)

df <- tibble(x = x, y1 = y1, y2 = y2) %>%
    gather(y1, y2, key = "line", value = "y") %>%
    mutate(
        line = str_replace(line, "y1", "f_norm'(x)"),
        line = str_replace(line, "y2", "-f_norm'(x - 2)"))


df %>% ggplot(aes(x, y, colour = line)) + geom_line()
ggsave("../assumption6_normal.png")

# Cauchy
y1 <- -2*x / (1 + x^2)^2
y2 <- 2*(x - 2/sqrt(3)) / (1 + (x - 2/sqrt(3))^2)^2

df <- tibble(x = x, y1 = y1, y2 = y2) %>%
    gather(y1, y2, key = "line", value = "y") %>%
    mutate(
        line = str_replace(line, "y1", "f_Cauchy'(x)"),
        line = str_replace(line, "y2", "-f_Cauchy'(x - 2)"))

df %>% ggplot(aes(x, y, colour = line)) + geom_line()
ggsave("../assumption6_cauchy.png")
