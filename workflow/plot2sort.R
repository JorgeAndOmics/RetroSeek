# ----------------------------
# PLOT GENERATION SCRIPT
# ----------------------------


# ----------------------------
# LOAD LIBRARIES & OPTIONS
# ----------------------------
options(warn = -1)  # Suppress warnings

suppressMessages({
  library(argparse)     # Command-line argument parser
  library(arrow)        # For reading Parquet files
  library(tidyverse)    # Data manipulation and visualization
  library(yaml)         # For reading YAML configuration
  library(ggsci)        # Scientific color palettes (e.g., Futurama)
  library(ggalluvial)   # Sankey/Alluvial plots
  library(ggdist)       # Raincloud plots and advanced distributions
})


# ----------------------------
# 2. COLOR PALETTE FUNCTION
# ----------------------------
futurama_unlimited_palette <- function(input_colour_number = 12, output_colour_number) {
  planet_express <- pal_futurama("planetexpress")(input_colour_number)
  output_colour  <- colorRampPalette(planet_express)(output_colour_number)
  return(output_colour)
}


# ----------------------------
# 3. PARSE COMMAND-LINE ARGS
# ----------------------------
parser <- ArgumentParser(description = "Generate plots from enERVate Parquet results")

parser$add_argument("--input",   required = TRUE, help = "Input directory with *_main.parquet and *_accessory.parquet files")
parser$add_argument("--output",  required = TRUE, help = "Directory to save output plots")
parser$add_argument("--config",  required = TRUE, help = "YAML config file with plot parameters")

args <- parser$parse_args()

args.input.dir       <- args$input
args.output.dir      <- args$output
args.config.file     <- args$config

# Read threshold from YAML
plot_filter_threshold <- yaml::read_yaml(args.config.file)$plots$omit_lower_percent %||% 0.99
plot_dpi              <- yaml::read_yaml(args.config.file)$plots$dpi %||% 300
plot_height           <- yaml::read_yaml(args.config.file)$plots$height %||% 12
plot_width            <- yaml::read_yaml(args.config.file)$plots$width %||% 15

cat("Generating RetroSeek plots...\n")


# ----------------------------
# 4. READ PARQUET INPUT FILES
# ----------------------------
main_files      <- list.files(args.input.dir, pattern = "_main\\.parquet$", full.names = TRUE)
accessory_files <- list.files(args.input.dir, pattern = "_accessory\\.parquet$", full.names = TRUE)

if (length(main_files) == 0)      stop("No main parquet files found.")
if (length(accessory_files) == 0) stop("No accessory parquet files found.")

all.main      <- purrr::map_dfr(main_files, arrow::read_parquet)
all.accessory <- purrr::map_dfr(accessory_files, arrow::read_parquet)
all.full      <- bind_rows(all.main, all.accessory)


# ----------------------------
# 4B. MAP ACCESSION TO SPECIES NAME
# ----------------------------

# Load species mapping from YAML config
species_map <- yaml::read_yaml(args.config.file)$species

# Convert named list to tibble for join
species_df <- tibble(
  species      = names(species_map),
  species_name = unname(unlist(species_map))
)

# Join species_name into all datasets
all.main      <- all.main      %>% left_join(species_df, by = "species")
all.accessory <- all.accessory %>% left_join(species_df, by = "species")
all.full      <- all.full      %>% left_join(species_df, by = "species")

all.main <- all.main %>% select(-species) %>% rename(species = species_name)
all.accessory <- all.accessory %>% select(-species) %>% rename(species = species_name)
all.full <- all.full %>% select(-species) %>% rename(species = species_name)

# ----------------------------
# 5. QUANTILE CALCULATIONS
# ----------------------------
q_stats <- function(df) {
  list(
    mean   = mean(df$mean_bitscore),
    q1     = quantile(df$mean_bitscore, 0.25),
    median = quantile(df$mean_bitscore, 0.50),
    q3     = quantile(df$mean_bitscore, 0.75)
  )
}

bit.main      <- q_stats(all.main)
bit.accessory <- q_stats(all.accessory)


# ----------------------------
# 6. COUNT PROBES PER GROUP
# ----------------------------
group_count <- function(df) {
  df %>%
    group_by(species, virus, probe, family, abbreviation) %>%
    summarise(count = n(), .groups = "keep") %>%
    ungroup()
}

all.counted_probe      <- group_count(all.full)
main.counted_probe     <- group_count(all.main)
accessory.counted_probe <- group_count(all.accessory)

# ----------------------------
# 7. DENSITY PLOTS
# ----------------------------
density_bitscore_plot <- function(data, q1, median, q3) {
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))
  
  ggplot(data, aes(x = mean_bitscore, fill = probe, colour = probe, weight = mean_identity)) +
    geom_density(alpha = 0.4, adjust = 3) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    geom_vline(xintercept = c(q1, median, q3), linetype = "dashed", linewidth = 0.2) +
    annotate("text", x = q1 + 20,     y = 1e-4, label = "Q1") +
    annotate("text", x = median + 20, y = 1e-4, label = "Q2") +
    annotate("text", x = q3 + 20,     y = 1e-4, label = "Q3") +
    scale_x_continuous(breaks = seq(0, max(data$mean_bitscore, na.rm = TRUE), by = 100)) +
    labs(x = "HSP Bitscore", y = "Density") +
    theme_minimal() +
    theme(
      text            = element_text(face = "bold"),
      axis.title.x    = element_text(size = 15, margin = margin(b = 15)),
      axis.title.y    = element_text(size = 15, margin = margin(l = 10)),
      axis.text       = element_text(size = 12)
    )
}

all.main.density.plot      <- density_bitscore_plot(all.main, bit.main$q1, bit.main$median, bit.main$q3)
all.accessory.density.plot <- density_bitscore_plot(all.accessory, bit.accessory$q1, bit.accessory$median, bit.accessory$q3)


# ----------------------------
# 8. BAR PLOTS
# ----------------------------
bar_plot <- function(data) {
  data <- data %>%
    group_by(species, family) %>%
    summarise(count = sum(count), .groups = "drop")
  
  ggplot(data) +
    aes(x = species, y = count, fill = family, alpha = count) +
    geom_col(color = "black", linewidth = 0.2) +
    scale_fill_futurama() +
    scale_alpha_continuous(range = c(0.5, 1)) +
    theme_void() +
    labs(fill = "Group", y = "Count") +
    theme(
      axis.title     = element_text(size = 12, face = "bold"),
      axis.text.x    = element_text(size = 11, angle = 45, hjust = 1),
      axis.text.y    = element_text(size = 11)
    ) +
    guides(alpha = "none")
}

all.full.bar.plot      <- bar_plot(all.counted_probe)
all.main.bar.plot      <- bar_plot(main.counted_probe)
all.accessory.bar.plot <- bar_plot(accessory.counted_probe)


# ----------------------------
# 9. RAINCLOUD PLOTS
# ----------------------------
raincloud_bitscore_plot <- function(data) {
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))
  
  ggplot(data, aes(y = probe, x = mean_bitscore, fill = probe, color = probe)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = 0, alpha = 0.5, .width = c(0.5, 0.95)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    ggdist::stat_dots(side = "right", dotsize = 0.1, alpha = 0.01, binwidth = 0.2) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    theme_minimal() +
    labs(x = "HSP Bitscore", y = "Probe") +
    theme(
      text        = element_text(face = "bold"),
      axis.title  = element_text(size = 12),
      axis.text   = element_text(size = 11)
    )
}

main.full.raincloud.bitscore_probe.plot      <- raincloud_bitscore_plot(all.main)
accessory.full.raincloud.bitscore_probe.plot <- raincloud_bitscore_plot(all.accessory)


# ----------------------------
# 10. BALLOON PLOTS
# ----------------------------
balloon_virus_species_plot <- function(data) {
  manual_colours <- futurama_unlimited_palette(5, length(unique(data$probe)))
  
  ggplot(data, aes(x = abbreviation, y = factor(species, levels = rev(unique(species))))) +
    geom_point(aes(color = probe, fill = probe, alpha = count, size = count^3), shape = 15) +
    facet_wrap(~ probe) +
    scale_color_manual(values = manual_colours) +
    scale_fill_manual(values = manual_colours) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    theme_minimal() +
    labs(x = "Virus", y = "Species") +
    theme(
      text         = element_text(face = "bold"),
      axis.title   = element_text(size = 12),
      axis.text.x  = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 9)
    ) +
    guides(size = "none")
}

main.full.balloon.virus_species.plot      <- balloon_virus_species_plot(main.counted_probe)
accessory.full.balloon.virus_species.plot <- balloon_virus_species_plot(accessory.counted_probe)


# ----------------------------
# 11. SANKEY PLOTS
# ----------------------------

# Prepare data for each sankey
data.main.sankey.species.probe <-main.counted_probe %>% group_by(species, probe) %>% summarise(count = sum(count), .groups = "drop")
data.accessory.sankey.species.probe <- accessory.counted_probe %>% group_by(species, probe) %>% summarise(count = sum(count), .groups = "drop")
data.main.sankey.family.probe <- main.counted_probe %>% group_by(family, probe) %>% summarise(count = sum(count), .groups = "drop")
data.accessory.sankey.family.probe <- accessory.counted_probe %>% group_by(family, probe) %>% summarise(count = sum(count), .groups = "drop")
data.main.sankey.species.family <- main.counted_probe %>% group_by(species, family) %>% summarise(count = sum(count), .groups = "drop")
data.accessory.sankey.species.family <- accessory.counted_probe %>% group_by(species, family) %>% summarise(count = sum(count), .groups = "drop")

sankey_species_probe_plot <- function(data, filter_threshold = plot_filter_threshold) {
  grand_total <- sum(data$count)
  
  ordered_species <- sort(unique(data$species), decreasing = TRUE)
  ordered_probe <- sort(unique(data$probe), decreasing = TRUE)
  
  data <- data %>%
  mutate(
    species = factor(species, levels = ordered_species),
    probe   = factor(probe,   levels = ordered_probe)
  )

  num_colours     <- max(length(unique(data$species)), length(unique(data$probe)))
  manual_colours  <- futurama_unlimited_palette(12, num_colours)
  
  ggplot(data, aes(axis1 = species, axis2 = probe, y = count)) +
    geom_alluvium(aes(fill = species, alpha = count)) +
    geom_stratum(alpha = 0, color = "black", linewidth = 0.2) +
    geom_text(
      aes(
        size  = after_stat(count),
        label = sprintf("%s [%d]", after_stat(stratum), after_stat(count))
      ),
      stat      = "stratum",
      nudge_x   = -0.15,
      fontface  = "bold",
      hjust     = 0,
      direction = "y"
    ) +
    scale_size_continuous(range = c(2, 5)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    scale_fill_manual(values = manual_colours) +
    theme_void() +
    theme(legend.position = "none")
}


sankey_family_probe_plot <- function(data, filter_threshold = plot_filter_threshold) {
  grand_total <- sum(data$count)
  
  ordered_family <- sort(unique(data$family), decreasing = TRUE)
  ordered_probe <- sort(unique(data$probe), decreasing = TRUE)

  data <- data  %>%
  mutate(
    family = factor(family, levels = ordered_family),
    probe  = factor(probe,  levels = ordered_probe)
  )
  
  num_colours     <- max(length(unique(data$family)), length(unique(data$probe)))
  manual_colours  <- futurama_unlimited_palette(12, num_colours)
  
  ggplot(data, aes(axis1 = family, axis2 = probe, y = count)) +
    geom_alluvium(aes(fill = family, alpha = count)) +
    geom_stratum(alpha = 0, color = "black", linewidth = 0.2) +
    geom_text(
      aes(
        size  = after_stat(count),
        label = sprintf("%s [%d]", after_stat(stratum), after_stat(count))
      ),
      stat      = "stratum",
      nudge_x   = -0.15,
      fontface  = "bold",
      hjust     = 0,
      direction = "y"
    ) +
    scale_size_continuous(range = c(2, 5)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    scale_fill_manual(values = manual_colours) +
    theme_void() +
    theme(legend.position = "none")
}


sankey_species_family_plot <- function(data, filter_threshold = plot_filter_threshold) {
  grand_total <- sum(data$count)
  
  ordered_species <- sort(unique(data$species), decreasing = TRUE)
  ordered_family <- sort(unique(data$family), decreasing = TRUE)
  
    data <- data %>%
    mutate(
      species = factor(species, levels = ordered_species),
      family  = factor(family,  levels = ordered_family)
    )

  num_colours     <- max(length(unique(data$species)), length(unique(data$family)))
  manual_colours  <- futurama_unlimited_palette(12, num_colours)
  
  ggplot(data, aes(axis1 = species, axis2 = family, y = count)) +
    geom_alluvium(aes(fill = species, alpha = count)) +
    geom_stratum(alpha = 0, color = "black", linewidth = 0.2) +
    geom_text(
      aes(
        size  = after_stat(count),
        label = sprintf("%s [%d]", after_stat(stratum), after_stat(count))
      ),
      stat      = "stratum",
      nudge_x   = -0.15,
      fontface  = "bold",
      hjust     = 0,
      direction = "y"
    ) +
    scale_size_continuous(range = c(2, 5)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    scale_fill_manual(values = manual_colours) +
    theme_void() +
    theme(legend.position = "none")
}

main.full.sankey.species_probe_plot   <- sankey_species_probe_plot(data.main.sankey.species.probe)
accessory.full.sankey.species_probe_plot <- sankey_species_probe_plot(data.accessory.sankey.species.probe)

main.full.sankey.family_probe_plot    <- sankey_family_probe_plot(data.main.sankey.family.probe)
accessory.full.sankey.family_probe_plot <- sankey_family_probe_plot(data.accessory.sankey.family.probe)

main.full.sankey.species_family_plot  <- sankey_species_family_plot(data.main.sankey.species.family)
accessory.full.sankey.species_family_plot <- sankey_species_family_plot(data.accessory.sankey.species.family)

# ----------------------------
# 12. EXPORT PLOTS
# ----------------------------
save_plot <- function(name, plot, w = plot_width, h = plot_height, r = plot_dpi) {
  ggsave(filename = file.path(args.output.dir, name), plot = plot, width = w, height = h, dpi = r)
}

# Density
save_plot("main_density.png",      all.main.density.plot)
save_plot("accessory_density.png", all.accessory.density.plot)

# Raincloud
save_plot("main_raincloud.png",      main.full.raincloud.bitscore_probe.plot)
save_plot("accessory_raincloud.png", accessory.full.raincloud.bitscore_probe.plot)

# Bar
save_plot("full_bar.png",      all.full.bar.plot)
save_plot("main_bar.png",      all.main.bar.plot)
save_plot("accessory_bar.png", all.accessory.bar.plot)

# Balloon
save_plot("main_balloon.png",      main.full.balloon.virus_species.plot)
save_plot("accessory_balloon.png", accessory.full.balloon.virus_species.plot)

# Sankey
save_plot("main_sankey_a.png", main.full.sankey.species_probe_plot)
save_plot("main_sankey_b.png", main.full.sankey.family_probe_plot)
save_plot("main_sankey_c.png", main.full.sankey.species_family_plot)
save_plot("accessory_sankey_a.png", accessory.full.sankey.species_probe_plot)
save_plot("accessory_sankey_b.png", accessory.full.sankey.family_probe_plot)
save_plot("accessory_sankey_c.png", accessory.full.sankey.species_family_plot)
