# brand.R - the rTFT visual identity
#
# The palette, the hexagon geometry and the plot that draws it. Sourced by
# tools/make-logo.R and tools/make-card.R so that the colours and the shape have
# one definition and cannot drift apart.
#
# The sticker is a hexagon filled with the thermal gradient, cold on the left and
# hot on the right. Across it run two curves that tell the story of the package:
# the failure rate, which rises exponentially with temperature, and the survival
# probability, which falls as the exposure accumulates. They cross once, and the
# dot marks where. A deep blue band along the bottom carries the name.
#
# Every mark is built as a polygon in the hexagon's own coordinates (circumradius
# 1, so the sticker is 2 units tall) and cut to the hexagon, rather than as a line
# of a given width on a device. Widths are therefore fixed fractions of the
# sticker height, 1/300 of a unit per 1/600 of the height, and nothing can spill
# past the border however large or small the plot is drawn. The gradient is a few
# hundred thin strips for the same reason: it does not depend on the graphics
# device supporting gradient fills.
#
# The name is set in JetBrains Mono Bold, read from tools/fonts/ and registered for
# the PNG device. The script stops if the font files are missing; it never falls
# back to another face. systemfonts, ragg, sf and scales are used by these scripts
# only and are not dependencies of the package.
#
# Nothing here writes a file. The callers do that.
#
# Inputs:
#   DESCRIPTION (the Title field, used verbatim on the card)
#   tools/fonts/JetBrainsMono-Regular.ttf
#   tools/fonts/JetBrainsMono-Bold.ttf
#
# Outputs:
#   none

# Setup ------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tibble)
library(sf)
library(scales)

# The palette ------------------------------------------------------------------

# Cold to hot. These four stops are the package's colour identity: the ground of
# the sticker, the rule on the card and any temperature axis in the manuscript are
# drawn from this vector rather than from a fresh guess.

rtft_stops <- tribble(
  ~position, ~hex,       ~role,
  0.00,      "#3D74B8",  "cold end of the assay range",
  0.35,      "#2E8B8B",  "lower-middle",
  0.70,      "#E3B23C",  "upper-middle",
  1.00,      "#D1462F",  "hot end, at or past the threshold"
)

rtft_thermal   <- rtft_stops$hex
rtft_deep_blue <- "#1D3F6E"   # border, name band, rate curve, outlines; card ground
rtft_cream     <- "#F2E8DC"   # survival curve, dot, name

# The font ---------------------------------------------------------------------

# JetBrains Mono (SIL Open Font License, tools/fonts/OFL.txt) is not a system font.
# The Regular and Bold files, from the official 2.304 release, are kept in
# tools/fonts/ and registered with systemfonts under the family name
# "jetbrains_mono", which is where the PNG device the callers save with (ragg)
# looks fonts up. Nothing is fetched at run time. Registration is checked by asking
# systemfonts which family it resolves the name to.

rtft_font_files <- c(regular = file.path("tools", "fonts", "JetBrainsMono-Regular.ttf"),
                     bold    = file.path("tools", "fonts", "JetBrainsMono-Bold.ttf"))

if (!all(file.exists(rtft_font_files))) {
  stop("JetBrains Mono is missing from tools/fonts/: ",
       paste(rtft_font_files[!file.exists(rtft_font_files)], collapse = ", "),
       ". The sticker is drawn in that face only; run from the package root with the ",
       "font files in place.", call. = FALSE)
}

systemfonts::register_font("jetbrains_mono",
                           plain = normalizePath(rtft_font_files[["regular"]]),
                           bold  = normalizePath(rtft_font_files[["bold"]]))
rtft_font <- "jetbrains_mono"

local({
  resolved <- systemfonts::font_info(rtft_font, weight = "bold")
  stopifnot(
    "the font family resolved is not JetBrains Mono" =
      grepl("JetBrains Mono", resolved$family, fixed = TRUE) &&
      grepl("JetBrainsMono-Bold", resolved$path, fixed = TRUE)
  )
})

# The hexagon: pointy-top, circumradius 1, the R sticker proportion -------------

# Seven points, not six: the last repeats the first so that the outline closes.

rtft_hex_angles <- seq(90, 450, by = 60) * pi / 180

rtft_hexagon <- tibble(
  x = cos(rtft_hex_angles) * sqrt(3) / 2 / cos(pi / 6),
  y = sin(rtft_hex_angles)
) |>
  (\(d) { d[7, ] <- d[1, ]; d })()      # the last point is the first, exactly

rtft_x_edge <- max(rtft_hexagon$x)

# Half of the vertical extent of the hexagon at x; the boundary is |y| = this.
rtft_hex_half_height <- function(x) 1 - abs(x) / sqrt(3)

stopifnot(
  "the hexagon is not the R sticker proportion" =
    isTRUE(all.equal(rtft_x_edge, sqrt(3) / 2, tolerance = 1e-6)),
  "the border path is not closed" =
    nrow(rtft_hexagon) == 7L &&
    isTRUE(all.equal(unlist(rtft_hexagon[1, ]), unlist(rtft_hexagon[7, ]))),
  "the hexagon boundary is not |y| = 1 - |x| / sqrt(3)" =
    isTRUE(all.equal(rtft_hexagon$y[2], rtft_hex_half_height(rtft_hexagon$x[2]),
                     tolerance = 1e-6))
)

rtft_hex_sf <- st_sfc(st_polygon(list(as.matrix(rtft_hexagon))))

# Widths, as fractions of the sticker height (600 units; 1 unit = 1/300) --------

rtft_unit <- 1 / 300

rtft_w <- list(
  border         = 12 * rtft_unit,
  rate           = 12 * rtft_unit,
  survival_cream = 12 * rtft_unit,
  survival_blue  = 19 * rtft_unit,
  dot_radius     = 17 * rtft_unit,
  dot_outline    =  6 * rtft_unit,
  name_height    = 72 * rtft_unit
)

# The name band and the name ---------------------------------------------------

rtft_band_top   <- -0.464
rtft_name_width <- 0.72
rtft_name_label <- "rTFT"

# JetBrains Mono is monospaced: every glyph advances 0.6 em and a capital stands
# 0.73 em tall. The em that makes four glyphs 0.72 wide, and the cap height that
# follows, are read off those two numbers.
rtft_name_em  <- rtft_name_width / (nchar(rtft_name_label) * 0.6)
rtft_name_cap <- 0.73 * rtft_name_em

# The name is centred at y = -0.69 when it can be. The hexagon narrows toward the
# bottom point, and at that height the baseline corners of the outer letters would
# touch the border at the stated width. The centre is therefore raised, if need be,
# until the baseline clears the inner edge of the border by 0.03.
rtft_name_clearance <- 0.03
rtft_name_centre <- local({
  half_w <- rtft_name_width / 2 + rtft_name_clearance + rtft_w$border / 2
  y_base_min <- -(1 - half_w / sqrt(3))        # lowest baseline that still clears
  max(-0.69, y_base_min + rtft_name_cap / 2)
})

stopifnot(
  "the name does not clear the hexagon border at its baseline" =
    rtft_name_width / 2 + rtft_name_clearance + rtft_w$border / 2 <=
    sqrt(3) * (1 - abs(rtft_name_centre - rtft_name_cap / 2)) + 1e-9,
  "the name is not inside the band" =
    rtft_name_centre + rtft_name_cap / 2 < rtft_band_top
)

# The two curves ---------------------------------------------------------------

# Rate of failure: M(x) = c + A exp(k x), an exponential in temperature. It sits
# just above the name band on the cold side (floor c), passes through
# (-sqrt(3)/2, -0.36) and leaves through the upper right at (0.79, 0.79).

rtft_rate_floor <- -0.385
rtft_rate_k <- log((0.79 - rtft_rate_floor) / (-0.36 - rtft_rate_floor)) /
  (0.79 + rtft_x_edge)
rtft_rate_A <- (-0.36 - rtft_rate_floor) / exp(-rtft_rate_k * rtft_x_edge)

rtft_rate_fn <- function(x) rtft_rate_floor + rtft_rate_A * exp(rtft_rate_k * x)

# Survival: a logistic fall from a high shoulder (hi) on the cold side to a low
# one (lo) on the hot side, with midpoint m and steepness s. The four values were
# fitted so that the curve passes through (-sqrt(3)/2, 0.44), below the top of the
# vertical side at 0.5 so that the whole stroke meets that side with clearance, and
# through (sqrt(3)/2, -0.39), and so that it cuts the rate curve near x = 0.35.

rtft_surv_hi <- 0.446
rtft_surv_lo <- -0.52
rtft_surv_s  <- 4.0
rtft_surv_m  <- 0.397

rtft_survival_fn <- function(x) {
  rtft_surv_lo + (rtft_surv_hi - rtft_surv_lo) / (1 + exp(rtft_surv_s * (x - rtft_surv_m)))
}

# The crossing, found rather than typed.
rtft_cross_x <- uniroot(function(x) rtft_rate_fn(x) - rtft_survival_fn(x),
                        interval = c(-rtft_x_edge, rtft_x_edge), tol = 1e-10)$root
rtft_cross_y <- rtft_rate_fn(rtft_cross_x)

# Guardrails: one crossing, inside the hexagon, rate up, survival down.
local({
  x <- seq(-1, 1, length.out = 20001)
  d <- rtft_rate_fn(x) - rtft_survival_fn(x)
  xi <- seq(-rtft_x_edge, rtft_x_edge, length.out = 20001)
  stopifnot(
    "the rate curve is not strictly increasing" = all(diff(rtft_rate_fn(xi)) > 0),
    "the survival curve is not strictly decreasing" =
      all(diff(rtft_survival_fn(xi)) < 0),
    "the curves do not cross exactly once in the hexagon" =
      sum(diff(sign(rtft_rate_fn(xi) - rtft_survival_fn(xi))) != 0) == 1L,
    "the crossing is not inside the hexagon" =
      abs(rtft_cross_y) < rtft_hex_half_height(rtft_cross_x),
    "the dot is not on both curves" =
      abs(rtft_cross_y - rtft_survival_fn(rtft_cross_x)) < 1e-8
  )
})

# Geometry helpers -------------------------------------------------------------

# A curve of the given width, cut flat at the ends and to the hexagon.
rtft_ribbon <- function(fn, width) {
  x <- seq(-1, 1, length.out = 4001)
  line <- st_sfc(st_linestring(cbind(x, fn(x))))
  line |>
    st_buffer(dist = width / 2, endCapStyle = "FLAT", joinStyle = "ROUND",
              nQuadSegs = 16) |>
    st_intersection(rtft_hex_sf)
}

# Region under a curve, down to below the hexagon, cut to the hexagon.
rtft_under <- function(fn) {
  x <- seq(-1, 1, length.out = 4001)
  st_sfc(st_polygon(list(rbind(cbind(x, fn(x)), c(1, -1.2), c(-1, -1.2), c(-1, fn(-1)))))) |>
    st_intersection(rtft_hex_sf)
}

# Polygon coordinates as a data frame for geom_polygon (holes through subgroup).
rtft_sf_df <- function(geom) {
  co <- st_coordinates(st_cast(st_sfc(st_cast(geom, "MULTIPOLYGON")), "MULTIPOLYGON"))
  tibble(x = co[, "X"], y = co[, "Y"], poly = co[, "L2"], ring = co[, "L1"])
}

rtft_poly_layer <- function(geom, fill, alpha = 1) {
  d <- rtft_sf_df(geom)
  geom_polygon(data = d, aes(x, y, group = poly, subgroup = ring),
               fill = fill, colour = NA, alpha = alpha)
}

# The thermal ground: thin vertical strips, each cut to the hexagon at its x.
rtft_ground_strips <- function(n_strips = 400) {
  edges <- seq(-rtft_x_edge, rtft_x_edge, length.out = n_strips + 1)
  ramp <- gradient_n_pal(rtft_thermal, values = rtft_stops$position)
  tibble(x0 = head(edges, -1), x1 = tail(edges, -1)) |>
    mutate(strip = row_number(),
           fill  = ramp(((x0 + x1) / 2 + rtft_x_edge) / (2 * rtft_x_edge)),
           x1    = pmin(x1 + 3 * (x1 - x0), rtft_x_edge)) |>   # each strip runs under
                                                               # its next three, so no
                                                               # seam can open
    rowwise() |>
    reframe(strip = strip, fill = fill,
            x = c(x0, x0, x1, x1),
            y = c(-rtft_hex_half_height(x0), rtft_hex_half_height(x0),
                  rtft_hex_half_height(x1), -rtft_hex_half_height(x1)))
}

# The hexagon plot -------------------------------------------------------------

# Returns the sticker as a ggplot. Both callers use this: make-logo.R saves it,
# make-card.R places it on the card. Because the card composes the plot and not a
# saved PNG, the card cannot go stale against a logo it never reads.
#
# Only the name has a size in absolute units. `scale` is the fraction of the
# reference width (4 inches, the logo) the plot will occupy, so the name keeps its
# proportion when the plot is placed small. `border_colour` is deep blue on the
# logo and cream on the card, where deep blue would vanish into the ground.

rtft_ref_width_in <- 4
rtft_ref_height_px <- 1385
rtft_ref_width_px  <- 1200

rtft_hex_plot <- function(wordmark = TRUE, scale = 1, border_colour = rtft_deep_blue) {
  pad  <- rtft_w$border / 2
  x_lim <- rtft_x_edge + pad
  y_lim <- x_lim * rtft_ref_height_px / rtft_ref_width_px
  stopifnot("the plot window clips the border" = y_lim >= 1 + pad)

  hex_ring <- st_difference(st_buffer(rtft_hex_sf, pad, joinStyle = "ROUND"),
                            st_buffer(rtft_hex_sf, -pad))
  band <- st_intersection(rtft_hex_sf,
                          st_sfc(st_polygon(list(rbind(c(-2, rtft_band_top), c(2, rtft_band_top),
                                                       c(2, -2), c(-2, -2), c(-2, rtft_band_top))))))
  cross_pt <- st_sfc(st_point(c(rtft_cross_x, rtft_cross_y)))
  dot_blue  <- st_buffer(cross_pt, rtft_w$dot_radius, nQuadSegs = 32)
  dot_cream <- st_buffer(cross_pt, rtft_w$dot_radius - rtft_w$dot_outline, nQuadSegs = 32)

  strips <- rtft_ground_strips()

  p <- ggplot() +
    geom_polygon(data = strips, aes(x, y, group = strip, fill = I(fill)), colour = NA) +
    rtft_poly_layer(rtft_under(rtft_survival_fn), rtft_cream, alpha = 0.38) +
    rtft_poly_layer(rtft_ribbon(rtft_rate_fn, rtft_w$rate), rtft_deep_blue) +
    rtft_poly_layer(rtft_ribbon(rtft_survival_fn, rtft_w$survival_blue), rtft_deep_blue) +
    rtft_poly_layer(rtft_ribbon(rtft_survival_fn, rtft_w$survival_cream), rtft_cream) +
    rtft_poly_layer(dot_blue, rtft_deep_blue) +
    rtft_poly_layer(dot_cream, rtft_cream)

  if (wordmark) {
    # em in data units -> points: the plot is 2 * x_lim data units across.
    size_pt <- rtft_name_em * (rtft_ref_width_in * scale / (2 * x_lim)) * 72
    p <- p +
      rtft_poly_layer(band, rtft_deep_blue) +
      annotate("text", x = 0, y = rtft_name_centre, label = rtft_name_label,
               colour = rtft_cream, family = rtft_font, fontface = "bold",
               size = size_pt / .pt, hjust = 0.5, vjust = 0.5)
  }

  p +
    rtft_poly_layer(hex_ring, border_colour) +
    coord_fixed(xlim = c(-x_lim, x_lim), ylim = c(-y_lim, y_lim), expand = FALSE) +
    theme_void() +
    theme(plot.background  = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.margin      = margin(0, 0, 0, 0),
          legend.position  = "none")
}

# The package title, read from DESCRIPTION so the card cannot contradict it -----

rtft_title <- function(root = ".") {
  dcf <- file.path(root, "DESCRIPTION")
  if (!file.exists(dcf)) stop("DESCRIPTION not found; run from the package root",
                              call. = FALSE)
  unname(read.dcf(dcf, fields = "Title")[1, 1])
}
