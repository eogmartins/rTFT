# make-logo.R - the rTFT hex sticker
#
# Saves the hexagon defined in tools/brand.R. The geometry, the palette and the
# assertions live there; this file only chooses the output size and writes it.
#
# Inputs:
#   tools/brand.R
#   tools/fonts/JetBrainsMono-Regular.ttf, tools/fonts/JetBrainsMono-Bold.ttf
#     (read by brand.R)
#
# Outputs:
#   man/figures/logo.png   (1200 x 1384, transparent, R sticker proportion)

# Setup ------------------------------------------------------------------------

if (!file.exists("DESCRIPTION")) {
  stop("run this from the package root, not from tools/", call. = FALSE)
}
source(file.path("tools", "brand.R"))

# Write ------------------------------------------------------------------------

dir.create(file.path("man", "figures"), recursive = TRUE, showWarnings = FALSE)
out_png <- file.path("man", "figures", "logo.png")

ggsave(out_png, rtft_hex_plot(), width = 1200 / 300, height = 1385 / 300,
       dpi = 300, bg = "transparent", device = ragg::agg_png)

if (!file.exists(out_png)) stop("the sticker was not written", call. = FALSE)
message("wrote ", out_png)
