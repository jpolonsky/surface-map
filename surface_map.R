## surface_map(): smoothed geospatial surface of prevalence, incidence or a rate
## from survey data, clustered or not.
##
## Method (two steps)
## 1. Cluster level (clustered = TRUE). Cluster is the independent sampling unit, so individual rows
##    are first reduced to one value per cluster. Raw cluster values are unstable
##    when clusters are small (one positive of two people gives "50%"), so each
##    cluster's value is the empirical-Bayes (partial-pooling) prediction of a
##    random-intercept model, mgcv::gam(y ~ s(cluster, bs = "re")), which shrinks
##    small clusters toward the survey-wide value. The likelihood follows the type:
##      prevalence : binomial,  y = 0/1 per row
##      incidence  : Poisson with offset(log(exposure)), y = event count per row,
##                   exposure = follow-up time (person-years, person-days, ...)
##      rate       : same model as incidence, for any events-per-denominator rate
##                   (e.g. deaths per 10,000 person-months); only labels differ
##      mean       : Gaussian, y = continuous value per row (e.g. antibody level)
## 2. Surface. Cluster values, at one coordinate per cluster, are kernel-smoothed
##    over `polygon` with spatstat::Smooth(), each cluster weighted by its sample
##    size (or exposure). Smoothing runs in a projected CRS (`work_crs`); the
##    result is returned in the CRS of `polygon`.
##    Unclustered data (clustered = FALSE): every row is its own point, with no
##    aggregation or shrinkage, and the kernel smoother does all the pooling. For
##    incidence/rate the rows are weighted by exposure, which makes the surface
##    exactly smoothed events / smoothed exposure.
##
## Arguments
##   dat        data frame with one row per person (or household), including
##              coordinates, the outcome and, if clustered, a cluster id.
##   clustered  TRUE (default) if rows are grouped in sampling clusters; FALSE if
##              each row is an independent point.
##   cluster    unquoted cluster-id column; only read when clustered = TRUE.
##   polygon    sf polygon(s) bounding the surface (e.g. the study area).
##   type       "prevalence", "incidence", "rate" or "mean".
##   outcome    unquoted column. 0/1 or TRUE/FALSE for prevalence; non-negative
##              event count for incidence/rate; numeric for mean.
##   exposure   unquoted column of follow-up time or denominator per row; required
##              for incidence and rate, ignored otherwise. Rows with exposure <= 0
##              are dropped from the model.
##   lon, lat   unquoted coordinate columns, in `data_crs` (default WGS84). If `dat`
##              is an sf object they are not needed: coordinates come from its
##              geometry (centroids if not points), in the object's own CRS, and
##              `data_crs` is ignored.
##   cluster_xy "first" takes each cluster's first non-missing coordinate; "mean"
##              takes the mean of its coordinates. Ignored if clustered = FALSE.
##   per        multiplier for the legend scale. Defaults: prevalence 100 (percent),
##              incidence 1000, rate 10000, mean 1.
##   unit       legend wording after the quantity name, e.g. "per 1,000 person-years".
##              Defaults to "%" for prevalence, otherwise "per <per>".
##   legend_name  full legend title, overriding the automatic one.
##   shrink     FALSE uses raw cluster values (ignored if clustered = FALSE, where
##              there is nothing to shrink). Shrinkage also falls back to raw
##              values, with a message, if the model cannot be fitted (fewer than
##              two clusters, or a constant outcome).
##   gam_method smoothing-parameter method passed to mgcv::gam() for the shrinkage
##              model. "REML" (default) is recommended for random effects; "GCV.Cp"
##              is mgcv's own default and reproduces older scripts that did not set one.
##   smooth_weights  "n" (rows per cluster), "exposure" (summed exposure; incidence
##              and rate only) or "none". Default: "n" if clustered; if unclustered,
##              "exposure" for incidence/rate and "none" otherwise.
##   sigma      kernel bandwidth in units of `work_crs` (metres for 3857). NULL
##              lets spatstat choose by cross-validation, separately per call; pass
##              one shared value to smooth several panels at the same scale.
##   palette    "viridis" (default), any other viridis option ("magma", "inferno",
##              "plasma", "cividis", "mako", "rocket", "turbo"), or "rog"
##              (red-orange-green). reverse_palette flips it (e.g. for a measure
##              where high is good under "rog").
##   zero_floor TRUE floors the colour scale at 0 (absolute level); FALSE stretches
##              it to the surface's own range. `fill_limits` = c(min, max) overrides
##              both, e.g. for a shared scale across panels.
##   contours, contour_labels, contour_binwidth   dashed contour lines, their
##              numeric labels, and optional bin width.
##
## Returns a list: `plot` (ggplot), `surface` (sf grid of smoothed values, column
## `value`, already multiplied by `per`), `points` (sf points that were smoothed:
## cluster [row id if unclustered], n, exposure, raw, shrunk, value), `settings`.
##
## Example
##   res <- surface_map(df, study_area, type = "prevalence", outcome = rdt_pos,
##                      cluster = cluster, lon = lon, lat = lat)
##   res$plot
##   res <- surface_map(df, study_area, type = "prevalence", outcome = rdt_pos,
##                      clustered = FALSE, palette = "rog", reverse_palette = TRUE)
##   res <- surface_map(df, study_area, type = "incidence", outcome = cases,
##                      exposure = person_years, per = 1000,
##                      unit = "per 1,000 person-years")

pacman::p_load(sf, mgcv, spatstat, stars, metR, tidyverse)

surface_map <- function(
 dat, polygon,
 type = c("prevalence", "incidence", "rate", "mean"),
 outcome, exposure = NULL,
 clustered = TRUE, cluster = cluster, lon = lon, lat = lat,
 cluster_xy = c("first", "mean"),
 per = NULL, unit = NULL, legend_name = NULL,
 shrink = TRUE, gam_method = "REML", smooth_weights = NULL, sigma = NULL,
 data_crs = 4326, work_crs = 3857,
 palette = "viridis", reverse_palette = FALSE,
 zero_floor = TRUE, fill_limits = NULL,
 contours = TRUE, contour_labels = TRUE, contour_binwidth = NULL,
 boundary_linewidth = 0.5
) {
 type <- match.arg(type)
 cluster_xy <- match.arg(cluster_xy)
 is_rate <- type %in% c("incidence", "rate")
 is_sf <- inherits(dat, "sf")
 if (is_sf) {
  geom <- st_geometry(dat)
  if (is.na(st_crs(geom))) stop("`dat` is an sf object with no CRS.", call. = FALSE)
  if (!all(st_geometry_type(geom) %in% "POINT")) geom <- suppressWarnings(st_centroid(geom))
  sf_xy <- st_coordinates(geom)[, 1:2, drop = FALSE]
  data_crs <- st_crs(geom)
  dat <- st_drop_geometry(dat) |> mutate(.sf_x = sf_xy[, 1], .sf_y = sf_xy[, 2])
 }
 smooth_weights <-
  if (is.null(smooth_weights)) {
   if (clustered) "n" else if (is_rate) "exposure" else "none"
  } else {
   match.arg(smooth_weights, c("n", "exposure", "none"))
  }
 shrink <- clustered && shrink
 viridis_options <- c("viridis", "magma", "inferno", "plasma", "cividis", "mako", "rocket", "turbo")
 if (!palette %in% c("rog", viridis_options)) {
  stop("Invalid `palette`: use 'rog' or a viridis option (", str_c(viridis_options, collapse = ", "), ").", call. = FALSE)
 }

 exposure_quo <- enquo(exposure)
 if (is_rate && rlang::quo_is_null(exposure_quo)) {
  stop("`exposure` is required when type is 'incidence' or 'rate'.", call. = FALSE)
 }
 if (smooth_weights == "exposure" && !is_rate) {
  stop("smooth_weights = 'exposure' needs type 'incidence' or 'rate'.", call. = FALSE)
 }
 if (!inherits(polygon, c("sf", "sfc"))) stop("`polygon` must be an sf object.", call. = FALSE)

 if (is.null(per)) per <- c(prevalence = 100, incidence = 1000, rate = 10000, mean = 1)[[type]]
 if (is.null(unit)) unit <- if (type == "prevalence") "%" else if (type == "mean") "" else str_c("per ", format(per, big.mark = ","))
 if (is.null(legend_name)) {
  legend_name <- c(prevalence = "Prevalence", incidence = "Incidence", rate = "Rate", mean = "Mean")[[type]]
  if (unit != "") legend_name <- str_c(legend_name, " (", unit, ")")
 }

 ## ---- individual-level frame -------------------------------------------------
 df <-
  dat |>
  transmute(
   .cl = if (clustered) factor({{ cluster }}) else factor(row_number()),
   .y = as.numeric({{ outcome }}),
   .e = if (is_rate) as.numeric(!!exposure_quo) else NA_real_,
   .lon = if (is_sf) .sf_x else {{ lon }},
   .lat = if (is_sf) .sf_y else {{ lat }}
  )

 df_model <- df |> filter(!is.na(.cl), !is.na(.y), if (is_rate) !is.na(.e) & .e > 0 else TRUE) |> mutate(.cl = droplevels(.cl))
 if (nrow(df_model) == 0) stop("No rows left with a non-missing cluster, outcome and exposure.", call. = FALSE)
 if (type == "prevalence" && !all(df_model$.y %in% c(0, 1))) {
  stop("For type 'prevalence', `outcome` must be 0/1 or TRUE/FALSE (convert factors first).", call. = FALSE)
 }
 if (is_rate && any(df_model$.y < 0)) stop("For incidence/rate, `outcome` must be a non-negative count.", call. = FALSE)

 ## ---- cluster level ------------------------------------------------------------
 df_cluster <-
  df_model |>
  group_by(.cl) |>
  summarise(
   n = n(),
   exposure = sum(.e),
   raw = if (is_rate) sum(.y) / sum(.e) else mean(.y),
   .groups = "drop"
  )

 shrunk <- NA_real_
 can_fit <- shrink && n_distinct(df_model$.cl) >= 2 && n_distinct(df_model$.y) >= 2
 if (shrink && !can_fit) {
  message("surface_map(): shrinkage skipped (fewer than 2 clusters or constant outcome); using raw cluster values.")
 }
 if (can_fit) {
  fit <- tryCatch(
   if (is_rate) {
    gam(.y ~ s(.cl, bs = "re") + offset(log(.e)), data = df_model, family = poisson(), method = gam_method)
   } else {
    gam(.y ~ s(.cl, bs = "re"), data = df_model, family = if (type == "prevalence") binomial() else gaussian(), method = gam_method)
   },
   error = \(e) {
    message("surface_map(): shrinkage model failed (", conditionMessage(e), "); using raw cluster values.")
    NULL
   }
  )
  if (!is.null(fit)) {
   ## exposure = 1 makes the Poisson prediction a rate per one unit of exposure
   shrunk <- as.numeric(predict(fit, newdata = tibble(.cl = df_cluster$.cl, .e = 1), type = "response"))
  }
 }
 df_cluster <- df_cluster |> mutate(shrunk = shrunk, value = coalesce(shrunk, raw) * per)

 ## ---- one coordinate per cluster ------------------------------------------------
 df_xy <-
  df |>
  filter(!is.na(.cl), !is.na(.lon), !is.na(.lat)) |>
  group_by(.cl) |>
  summarise(
   .lon = if (cluster_xy == "first") first(.lon) else mean(.lon),
   .lat = if (cluster_xy == "first") first(.lat) else mean(.lat),
   .groups = "drop"
  )

 clusters <- df_cluster |> inner_join(df_xy, by = ".cl")
 n_dropped <- nrow(df_cluster) - nrow(clusters)
 if (n_dropped > 0) message("surface_map(): ", n_dropped, " cluster(s)/row(s) dropped for missing coordinates.")
 if (nrow(clusters) < 3) stop("Fewer than 3 clusters/rows with coordinates; cannot smooth a surface.", call. = FALSE)

 clusters <-
  clusters |>
  rename(cluster = .cl) |>
  st_as_sf(coords = c(".lon", ".lat"), crs = data_crs)

 ## ---- smoothing -----------------------------------------------------------------
 out_crs <- st_crs(polygon)
 polygon_proj <- polygon |> st_union() |> st_make_valid() |> suppressWarnings(st_collection_extract("POLYGON")) |> st_transform(work_crs) |> st_sf()
 clusters_proj <- st_transform(clusters, work_crs)
 xy <- st_coordinates(clusters_proj)

 ## as.owin() on the sf object handles multi-part and holed polygons
 win <- as.owin(polygon_proj)
 keep <- inside.owin(xy[, 1], xy[, 2], win)
 if (sum(keep) < 3) stop("Fewer than 3 clusters fall inside `polygon`; check `data_crs` and `polygon`.", call. = FALSE)
 if (any(!keep)) message("surface_map(): ", sum(!keep), " cluster(s)/row(s) fall outside `polygon` and are ignored.")
 w <- switch(smooth_weights, n = clusters_proj$n, exposure = clusters_proj$exposure, none = rep(1, nrow(clusters_proj)))
 ppp_obj <- ppp(x = xy[keep, 1], y = xy[keep, 2], marks = clusters_proj$value[keep], window = win)
 surface_proj <-
  ppp_obj |>
  Smooth(sigma = sigma, weights = w[keep]) |>
  stars::st_as_stars() |>
  st_as_sf() |>
  st_set_crs(work_crs) |>
  rename(value = v)

 ## Contour lines need an exactly regular x/y grid. Transforming each cell
 ## centroid separately leaves floating-point jitter (stat_contour then fails),
 ## so the grid axes are transformed once each and re-expanded. This assumes the
 ## plot CRS axes run parallel to `work_crs` axes (true for 3857 -> 4326).
 cen <- surface_proj |> mutate(geometry = st_centroid(geometry))
 cen_xy <- st_coordinates(cen)
 xs <- sort(unique(round(cen_xy[, 1], 6)))
 ys <- sort(unique(round(cen_xy[, 2], 6)))
 axis_lon <- st_sfc(map2(xs, ys[1], \(x, y) st_point(c(x, y))), crs = work_crs) |> st_transform(out_crs) |> st_coordinates()
 axis_lat <- st_sfc(map2(xs[1], ys, \(x, y) st_point(c(x, y))), crs = work_crs) |> st_transform(out_crs) |> st_coordinates()
 grid_points <-
  tibble(
   x = axis_lon[match(round(cen_xy[, 1], 6), xs), 1],
   y = axis_lat[match(round(cen_xy[, 2], 6), ys), 2],
   value = cen$value
  )

 surface <- st_transform(surface_proj, out_crs)

 ## ---- plot ------------------------------------------------------------------------
 vmin <- min(grid_points$value, na.rm = TRUE)
 vmax <- max(grid_points$value, na.rm = TRUE)
 limits <- if (is.null(fill_limits)) c(if (zero_floor) min(0, vmin) else vmin, vmax) else fill_limits

 fill_scale <-
  if (palette == "rog") {
   cols <- c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850")
   scale_fill_gradientn(colours = if (reverse_palette) rev(cols) else cols, name = legend_name, limits = limits)
  } else {
   scale_fill_viridis_c(option = palette, direction = if (reverse_palette) -1 else 1, name = legend_name, limits = limits)
  }

 contour_layers <- list()
 if (contours) {
  contour_layers <- list(
   geom_contour(data = grid_points, aes(x = x, y = y, z = value), colour = "grey20", linewidth = 0.3, linetype = "dashed", binwidth = contour_binwidth)
  )
  if (contour_labels) {
   contour_layers <- c(
    contour_layers,
    list(geom_text_contour(data = grid_points, aes(x = x, y = y, z = value), colour = "black", size = 3, skip = 0, stroke = 0.1, binwidth = contour_binwidth))
   )
  }
 }

 p <-
  ggplot() +
  geom_sf(data = surface, aes(fill = value), colour = NA) +
  geom_sf(data = st_boundary(st_union(polygon)), linewidth = boundary_linewidth) +
  contour_layers +
  fill_scale +
  coord_sf(crs = out_crs) +
  theme_void()

 list(
  plot = p,
  surface = surface,
  points = clusters |> select(cluster, n, exposure, raw, shrunk, value),
  settings = list(type = type, clustered = clustered, per = per, unit = unit, legend_name = legend_name, sigma = sigma, shrink = shrink && !is.na(shrunk[1]), work_crs = work_crs)
 )
}
