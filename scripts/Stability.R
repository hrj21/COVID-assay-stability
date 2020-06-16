
# Setup and restore R environment -----------------------------------------

#renv::init() # initialize r enironment file in current project
renv::restore() # install all packages listed in the environment lock file

# Load packages -----------------------------------------------------------

library(tidyverse)
library(gghighlight)

renv::snapshot()

# Read data ---------------------------------------------------------------

elisa <- read_csv("data/5 ELISA data (updated).csv")

# Plot development --------------------------------------------------------

elisa %>%
  ggplot(aes(Time, Abs, col = Type, linetype = Well_replicate, 
             group = interaction(Well, Well_replicate))) +
  facet_grid(Plate ~ Plate_replicate) +
  geom_line(size = 0.2) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(panel.spacing.x = unit(0, "lines"),
        panel.grid = element_blank())

# Intra-assay stability --------------------------------------------------

intra_assay <- elisa %>%
  mutate(
    Row = case_when(
      Type == "Blank" ~ 1L,
      Type == "Negative" ~ 6L,
      Type == "Positive" ~ 7L,
      TRUE ~ match(str_sub(Well, 1, 1), LETTERS[1:8])), 
    
    Column = case_when(
      Type == "Blank" ~ 6L,
      Type == "Negative" ~ 6L,
      Type == "Positive" ~ 6L,
      TRUE ~ as.integer(str_sub(Well, 2, 2))),
    
    Well = case_when(
      Type == "Blank" ~ "Blank wells",
      Type == "Positive" ~ "Positive wells",
      TRUE ~ Well)
    ) %>%
  
  group_by(Plate, Time, Plate_replicate, Type, Well, Column, Row) %>%
  summarize(Max = max(Abs), 
            Min = min(Abs),
            Range = Max - Min,
            Percent_diff = Range / ((Max + Min) / 2) * 100,
            CV = sd(Abs) / mean(Abs) * 100)

# Plot intra-assay stability metrics as lines -----------------------------

plot_lines <- function(data, variable) {
  data %>%
    ggplot(aes(Time, {{variable}}, col = Type, group = Well)) +
    facet_grid(Plate ~ Plate_replicate) +
    geom_line(size = 0.2) +
    geom_point(size = 1.5) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(panel.spacing.x = unit(0, "lines"),
          panel.grid = element_blank())
}

plot_lines(intra_assay, Range)
plot_lines(intra_assay, Percent_diff)
plot_lines(intra_assay, CV)

# Plot intra-assay stability metrics as heatmaps --------------------------

plot_heat <- function(data, variable, threshold) {
  data %>% 
    mutate(fill = min(c({{variable}}, {{threshold}}))) %>%
    ggplot(aes(Column, Row, fill = fill)) +
    facet_grid(Time ~ interaction(Plate, Plate_replicate), switch = "x") +
    geom_tile() +
    scale_x_continuous(breaks = 1:6,
                       labels = as.character(1:6),
                       expand = c(0, 0),
                       position = "top") +
    scale_y_reverse(breaks = 1:8, labels = LETTERS[1:8],
                    expand = c(0, 0)) +
    scale_fill_viridis_c() +
    coord_equal() +
    theme_bw() +
    labs(fill = deparse(substitute(variable))) +
    theme(axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.ontop = TRUE,
          panel.background = element_rect(fill = "transparent"))
}

plot_heat(intra_assay, Range, 0.5)
plot_heat(intra_assay, Percent_diff, 50)
plot_heat(intra_assay, CV, 25)

# Plot intra-assay stability metrics as points ----------------------------

plot_point <- function(data, variable) {
  data %>%
    group_by(Plate, Time, Type, Well, Plate_replicate) %>%
    summarize(across(.cols = c(Max, Min, Range, Percent_diff, CV),
                     .fun = mean)) %>%
    ggplot(aes(Type, {{variable}}, col = Plate_replicate)) +
    facet_grid(Time ~ Plate) +
    scale_color_brewer(type = "qual", palette = "Set1") +
    geom_point(position = position_dodge(0.3), size = 2) +
    theme_bw()
}

plot_point(intra_assay, Range)
plot_point(intra_assay, Percent_diff)
plot_point(intra_assay, CV)


# Calculate intra-assay CV ------------------------------------------------

intra_assay %>%
  group_by(Type) %>%
  summarize(across(Range:CV, mean))

# Inter-assay stability ---------------------------------------------------

inter_assay <- elisa %>%
  mutate(
    Row = match(str_sub(Well, 1, 1), LETTERS[1:8]), 
    
    Column = case_when(
      Well_replicate == "A" ~ str_sub(Well, 2, 2),
      Well_replicate == "B" ~ str_sub(Well, 5, -1)),
    
    Column = as.integer(Column)
  ) %>%
  
  group_by(Plate, Time, Plate_replicate, Type, Column, Row) %>%
  summarize(Mean = mean(Abs)) %>%
  group_by(Plate, Time, Type, Column, Row) %>%
  summarize(Range = range(Mean),
            Percent_diff = Range / ((max(Mean) + min(Mean)) / 2) * 100,
            CV = sd(Mean) / mean(Mean) * 100)

inter_assay

# Plot stability metrics as heatmaps --------------------------------------

plot_heat_inter <- function(data, variable, threshold) {
  data %>% 
    mutate(fill = min(c({{variable}}, {{threshold}}))) %>%
    ggplot(aes(Column, Row, fill = fill)) +
    facet_grid(Time ~ Plate, switch = "x") +
    geom_tile() +
    scale_x_continuous(breaks = 1:12,
                       labels = as.character(1:12),
                       expand = c(0, 0),
                       position = "top") +
    scale_y_reverse(breaks = 1:8, labels = LETTERS[1:8],
                    expand = c(0, 0)) +
    scale_fill_viridis_c() +
    coord_equal() +
    labs(fill = deparse(substitute(variable))) +
    theme_bw() +
    theme(panel.spacing.x = unit(0.5, "lines"),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.ontop = TRUE,
          panel.background = element_rect(fill = "transparent"))
}

plot_heat_inter(inter_assay, Range, 4)
plot_heat_inter(inter_assay, Percent_diff, 500)
plot_heat_inter(inter_assay, CV, 500)


# Plot stability metrics as points ----------------------------------------

plot_point_inter <- function(data, variable) {
  data %>% 
    #group_by(Plate, Time, Type) %>%
    summarize(across(Range:CV, mean)) %>%
    ggplot(aes(Type, {{variable}})) +
    facet_grid(Time ~ Plate) +
    scale_color_brewer(type = "qual", palette = "Set1") +
    geom_point(position = position_jitter(0.2), size = 2) +
    theme_bw()
}

plot_point_inter(inter_assay, Range)
plot_point_inter(inter_assay, CV)

# Calculate inter-assay CV ------------------------------------------------

inter_assay %>%
  group_by(Type, Time) %>%
  summarize(across(Range:CV, mean))
