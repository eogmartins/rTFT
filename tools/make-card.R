# make-card.R - the rTFT social card
#
# The image other applications show when the package URL is pasted into them --
# Slack, Teams, Bluesky, an email preview. It is NOT displayed on the pkgdown
# site; it is named in the page's og:image metadata. 1200 x 630 is the standard
# frame, so the hexagon is composed onto it rather than letterboxed into it.
#
# The hexagon is placed as the PLOT returned by tools/brand.R, not as the saved
# logo.png. The card therefore cannot go stale against a logo it never reads, and
# no image-reading package is needed.
#
# The tagline is DESCRIPTION's Title field, verbatim, so the card cannot
# contradict the package.
#
# Inputs:
#   tools/brand.R, DESCRIPTION
#   tools/fonts/JetBrainsMono-Regular.ttf, tools/fonts/JetBrainsMono-Bold.ttf
#     (read by brand.R)
#
# Outputs:
#   man/figures/card.png   (1200 x 630)
#
# Run from the package root:  Rscript tools/make-card.R

# Setup ------------------------------------------------------------------------

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root, not from tools/", call. = FALSE)
}
source(file.path("tools", "brand.R"))
library(grid)

# Turning the hexagon into a drawing that can be placed on the card needs a
# graphics device to measure it. With none open, R opens its default, which writes
# a stray Rplots.pdf in the working directory. A device that writes to no file is
# opened here in its place and closed once the card is saved.
grDevices::pdf(NULL)
measuring_device <- grDevices::dev.cur()

# Layout -----------------------------------------------------------------------

card_w <- 1200; card_h <- 630
hex_frac <- 0.230                     # the fraction of the card width it occupies
# No name or name band inside the hexagon: "rTFT" is set beside it, and twice is
# once too many. The border is cream because deep blue would vanish into the ground.
hex      <- ggplotGrob(rtft_hex_plot(wordmark = FALSE, scale = hex_frac,
                                     border_colour = rtft_cream))

# The rule under the wordmark is the thermal ramp itself, on the same scale the
# palette is defined over, so it is the palette rather than a decoration.
rule <- tibble(
  x = seq(0.335, 0.545, length.out = 300),
  y = 0.475
) |>
  mutate(v = seq(-rtft_x_edge, rtft_x_edge, length.out = 300))

tagline <- paste(strwrap(rtft_title(), width = 34), collapse = "\n")

p <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = rtft_deep_blue, colour = NA) +
  annotation_custom(hex, xmin = 0.075, xmax = 0.305, ymin = 0.155, ymax = 0.845) +
  annotate("text", x = 0.335, y = 0.655, label = "rTFT", hjust = 0, vjust = 0.5,
           colour = rtft_cream, family = rtft_font, fontface = "bold", size = 11.5) +
  geom_line(data = rule, aes(x, y, colour = v),
            linewidth = 0.9, lineend = "round") +
  annotate("text", x = 0.335, y = 0.330, label = tagline, hjust = 0, vjust = 0.5,
           colour = rtft_cream, family = "sans", size = 3.5,
           alpha = 0.82, lineheight = 1.3) +
  scale_colour_gradientn(colours = rtft_thermal,
                         limits = c(-rtft_x_edge, rtft_x_edge), guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0), legend.position = "none")

# Write ------------------------------------------------------------------------

dir.create(file.path("man", "figures"), recursive = TRUE, showWarnings = FALSE)
out_png <- file.path("man", "figures", "card.png")

ggsave(out_png, p, width = card_w / 300, height = card_h / 300, dpi = 300,
       device = ragg::agg_png)

grDevices::dev.off(measuring_device)

if (!file.exists(out_png)) stop("the card was not written", call. = FALSE)
if (file.exists("Rplots.pdf")) stop("a stray Rplots.pdf was written", call. = FALSE)
message("wrote ", out_png, "  (", card_w, " x ", card_h, ")")
message("tagline taken from DESCRIPTION: ", rtft_title())
