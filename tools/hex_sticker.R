# make a hex sticker for the cd4backcalc package
# pak::pak("emilioxavier/hexSticker")
library(hexSticker)
library(ggplot2)

# incidence curve
x <- seq(0, 13.5, length.out = 100)
y <- dnorm(x, mean = 5, sd = 1.3)

# delayed diagnosis curve
x2 <- seq(0, 13.6, length.out = 100)
y2 <- dnorm(x2, mean = 7, sd = 2)

# loop: starts at (anchor_x, 0), curves up, around, and returns to zero
t <- seq(0, 2 * pi, length.out = 600)
loop_df <- data.frame(
  x = 16.5 + 0.9 * (1 - cos(t)) * sin(t),
  y = 0.06 + 0.05 * (1 - cos(t))
)

tails <- rbind(
  data.frame(
    x = seq(13.5, 16.5, length.out = 50),
    y = exp(seq(-1, 3, length.out = 50)) * 0.003
  ),
  data.frame(
    x = seq(16.5, 18.65, length.out = 50),
    y = rev(exp(seq(0.85, 3, length.out = 50)) * 0.003)
  )
)

p <- ggplot(data = data.frame(x, y), aes(x = x, y = y)) +
  geom_line(linewidth = 2, color = "#444") +
  xlim(0, 20) +
  ylim(0, max(y) * 1.1) +
  geom_line(
    data = data.frame(x = x2, y = y2),
    aes(x = x, y = y),
    linewidth = 2,
    color = "red"
  ) +
  geom_path(
    data = loop_df,
    aes(x = x, y = y),
    color = "red",
    linewidth = 2
  ) +
  geom_path(
    data = tails,
    aes(x = x, y = y),
    color = "red",
    linewidth = 2
  ) +
  theme_void() +
  theme_transparent()

p

# create hex sticker
s <- sticker(
  p,
  package = "",
  url = "cd4backcalc",
  u_x = 0.23,
  u_y = 1.35,
  u_size = 13.5,
  s_x = 1.05,
  s_y = 1,
  s_width = 1.6,
  s_height = 1,
  h_fill = "#eee",
  h_color = "#444",
  filename = "figures/hex_sticker.png"
)

plot(s)
