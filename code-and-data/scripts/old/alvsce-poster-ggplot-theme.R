# title: alvsce poster ggplot theme
# author: madeleine wallace
# date: 2024-03-28


# PACKAGES ---------------------------------------------------------------------

library(ggplot2)
data(iris)

# THEME ------------------------------------------------------------------------

ggplot(iris, aes(x = Petal.Length)) +
  geom_histogram() 

?theme

alvsce_poster_theme <- function() {
  theme(
    # add border 1)
    panel.border = element_rect(colour = "black", fill = NA, linetype = 1),
    # color background 2)
    panel.background = element_rect(fill = "white"),
    # modify grid 3)
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "Verdana", size = 24),
    axis.title = element_text(colour = "black", face = "plain", family = "Verdana", size = 32),
    axis.ticks = element_line(colour = "black"),
    # aspect ratio
    # aspect.ratio = 3/5,
    # legend at the bottom 6)
    legend.position = "bottom"
  )
}

ggplot(iris, aes(x = Petal.Length)) +
  geom_histogram(fill = "#6A814B", col = "#6A814B") +
  alvsce_poster_theme()
