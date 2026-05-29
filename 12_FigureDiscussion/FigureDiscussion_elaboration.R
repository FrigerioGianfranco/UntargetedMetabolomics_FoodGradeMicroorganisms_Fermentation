
### Developing a new summary figure


library(tidyverse)
library(forcats)
library(ggtext)



polar_results <- read_csv("results_combined.csv")

lipid_results <- read_csv("results_combined_known_even_non_sign_lipids.csv")



discussion_tbl <- read_tsv("molecules_discussion.txt") %>%
  mutate(discussion_order = row_number(),
         section = str_replace(section,
                               "Lipid metabolism and membrane-associated signalling",
                               "Lipid metabolism and membrane-associated signaling"))


clean_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("γ", "gamma") %>%
    str_replace_all("–|—", "-") %>%
    iconv(from = "", to = "ASCII//TRANSLIT") %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish()
}

paste_existing_cols <- function(df, cols) {
  cols <- intersect(cols, names(df))
  
  if (length(cols) == 0) {
    return(rep("", nrow(df)))
  }
  
  df %>%
    select(all_of(cols)) %>%
    mutate(across(everything(), ~replace_na(as.character(.x), ""))) %>%
    unite("search_text", everything(), sep = " | ", remove = FALSE) %>%
    pull(search_text)
}




polar_prepared <- polar_results %>%
  mutate(source_table = "polar_metabolomics",
         result_row_id = row_number(),
         primary_name_text = paste_existing_cols(., c("Name", "Name.x", "Name.y", "Name_bis.x", "Name_bis.y")),
         extended_name_text = paste_existing_cols(., c("Name", "Name.x", "Name.y", "Name_bis.x", "Name_bis.y", "MSDialName")))

lipid_prepared <- lipid_results %>%
  mutate(source_table = "lipidomics",
         result_row_id = row_number(),
         primary_name_text = paste_existing_cols(., c("Name", "Name_bis")),
         extended_name_text = paste_existing_cols(., c("Name", "Name_bis", "MSDialName")))

all_results <- bind_rows(polar_prepared, lipid_prepared) %>%
  mutate(primary_name_clean = clean_name(primary_name_text),
         extended_name_clean = clean_name(extended_name_text))



molecule_aliases <- tribble(~molecule_or_class, ~alias,
                            
                            "Serotonin / 5-HT", "serotonin",
                            "Serotonin / 5-HT", "5-HT",
                            "Serotonin / 5-HT", "5-hydroxytryptamine",
                            
                            "Epinephrine", "epinephrine",
                            "Epinephrine", "adrenaline",
                            
                            "Octopamine", "octopamine",
                            
                            "Noradrenaline", "noradrenaline",
                            "Noradrenaline", "norepinephrine",
                            
                            "Tyramine", "tyramine",
                            
                            "Azelaic acid", "azelaic acid",
                            "Ricinoleic acid", "ricinoleic acid",
                            
                            "trans-Ferulic acid", "trans-ferulic acid",
                            "trans-Ferulic acid", "ferulic acid",
                            
                            "12-Hydroxydodecanoic acid", "12-hydroxydodecanoic acid",
                            
                            "13-HODE", "13-HODE",
                            "13-HODE", "13-hydroxyoctadecadienoic acid",
                            
                            "Linoleic acid", "linoleic acid",
                            "Palmitoleic acid", "palmitoleic acid",
                            
                            "Corticosterone", "corticosterone",
                            
                            "Orotic acid", "orotic acid",
                            
                            "6-(γ,γ-dimethylallylamino) purine", "6-(gamma,gamma-dimethylallylamino)purine",
                            "6-(γ,γ-dimethylallylamino) purine", "dimethylallylamino purine",
                            
                            "Biliverdin", "biliverdin",
                            
                            "Ophthalmic acid", "ophthalmic acid",
                            "Ophthalmic acid", "ophthalmate",
                            
                            "Histamine", "histamine",
                            
                            "Glutathione", "glutathione") %>%
  mutate(alias_clean = clean_name(alias))


pairwise_direction_cols <- names(all_results) %>%
  str_subset("^MIC_.+_vs_MIC_.+$") %>%
  str_subset("Pvalue|PvalueFDR|others", negate = TRUE)

pairwise_fdr_cols <- paste0(pairwise_direction_cols, "_PvalueFDR")
pairwise_fdr_cols <- intersect(pairwise_fdr_cols, names(all_results))

pairwise_pvalue_cols <- paste0(pairwise_direction_cols, "_Pvalue")
pairwise_pvalue_cols <- intersect(pairwise_pvalue_cols, names(all_results))



matches <- merge(all_results, molecule_aliases, by = NULL) %>%
  as_tibble() %>%
  mutate(match_primary = str_detect(primary_name_clean,
                                    regex(paste0("\\b", alias_clean, "\\b"))),
         match_extended = str_detect(extended_name_clean,
                                     regex(paste0("\\b", alias_clean, "\\b"))),
         match_quality = case_when(match_primary ~ "primary_name_match",
                                   match_extended ~ "extended_annotation_match",
                                   TRUE ~ NA_character_)) %>%
  filter(match_primary | match_extended)


summary_with_pairwise_results_long <- discussion_tbl %>%
  left_join(matches %>%
              select(molecule_or_class,
                     matched_alias = alias,
                     match_quality,
                     source_table,
                     result_row_id,
                     chromatographic_run,
                     Feature_ID,
                     Level,
                     MSDialName,
                     Name,
                     primary_name_text,
                     all_of(pairwise_direction_cols),
                     all_of(pairwise_fdr_cols)),
            by = "molecule_or_class")

write_tsv(summary_with_pairwise_results_long,
          "discussion_molecules_with_pairwise_results_LONG.tsv")


collapse_unique <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- x[x != ""]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(x, collapse = "; ")
}

keep_best_level <- function(x) {
  x <- na.omit(as.character(x))
  x <- x[x != ""]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  x_split <- str_split(x, pattern = ";") %>%
    unlist() %>%
    str_squish()
  
  numeric_levels <- suppressWarnings(as.numeric(x_split))
  numeric_levels <- numeric_levels[!is.na(numeric_levels)]
  
  if (length(numeric_levels) == 0) {
    return(paste(unique(x_split), collapse = "; "))
  }
  
  as.character(min(numeric_levels))
}

level_to_balls <- function(level) {
  dplyr::case_when(is.na(level) ~ "",
                   level == "1" ~ "●",
                   level == "2" ~ "●●",
                   level == "3" ~ "●●●",
                   TRUE ~ "")
}

summary_with_pairwise_results_compact <- summary_with_pairwise_results_long %>%
  group_by(section,
           molecule_or_class,
           associated_strain_or_context,
           category_level) %>%
  summarise(across(c(matched_alias,
                     match_quality,
                     source_table,
                     chromatographic_run,
                     Feature_ID,
                     MSDialName,
                     Name,
                     all_of(pairwise_direction_cols),
                     all_of(pairwise_fdr_cols)),
                   collapse_unique),
            .groups = "drop")

write_tsv(summary_with_pairwise_results_compact,
          "discussion_molecules_with_pairwise_results_COMPACT.tsv")

unmatched_molecules <- discussion_tbl %>%
  filter(category_level == "Specific molecule") %>%
  anti_join(matches %>% distinct(molecule_or_class),
            by = "molecule_or_class")


section_order <- c("Aromatic amino acid metabolism and neuroactive compounds",
                   "Organic acids and antimicrobial or barrier-associated functions",
                   "Lipid metabolism and membrane-associated signaling",
                   "Nucleotide metabolism and cellular processes",
                   "Redox-related metabolism and oxidative stress")

summary_with_pairwise_results_long <- summary_with_pairwise_results_long %>%
  mutate(section = str_replace(section,
                               "Lipid metabolism and membrane-associated signalling",
                               "Lipid metabolism and membrane-associated signaling"))

pairwise_direction_cols <- names(summary_with_pairwise_results_long) %>%
  str_subset("^MIC_.+_vs_MIC_.+$") %>%
  str_subset("Pvalue|PvalueFDR|others", negate = TRUE)

plot_tbl_raw <- summary_with_pairwise_results_long %>%
  filter(category_level == "Specific molecule") %>%
  filter(!is.na(Feature_ID)) %>%
  filter(section %in% section_order) %>%
  pivot_longer(cols = all_of(pairwise_direction_cols),
               names_to = "comparison",
               values_to = "pairwise_result")

plot_tbl <- plot_tbl_raw %>%
  group_by(section,
           discussion_order,
           molecule_or_class,
           associated_strain_or_context,
           comparison) %>%
  summarise(n_matched_features = n_distinct(Feature_ID),
            matched_features = paste(sort(unique(Feature_ID)), collapse = "; "),
            annotation_level = keep_best_level(Level),
            pairwise_result_collapsed = collapse_unique(pairwise_result),
            n_distinct_pairwise_results = n_distinct(na.omit(pairwise_result)),
            .groups = "drop") %>%
  mutate(annotation_level = replace_na(annotation_level, "not available"),
         annotation_balls = level_to_balls(annotation_level),
         section = factor(section, levels = section_order),
         
         first_group = str_extract(comparison, "^MIC_.+(?=_vs_)"),
         second_group = str_extract(comparison, "(?<=_vs_)MIC_.+$"),
         
         direction = case_when(is.na(pairwise_result_collapsed) ~ "Not significant",
                               n_distinct_pairwise_results > 1 ~ "Multiple directions",
                               str_detect(pairwise_result_collapsed, paste0(first_group, " > ", second_group)) ~ "Higher in first group",
                               str_detect(pairwise_result_collapsed, paste0(second_group, " > ", first_group)) ~ "Higher in second group",
                               TRUE ~ "Check manually"))


group_abbrev <- c(MIC_BIF_BREV = "BBRE",
                  MIC_LAC_AMYL = "LAMY",
                  MIC_LAC_KEFI = "LKEF",
                  MIC_LAC_PLAN = "LPLA",
                  MIC_LAC_SALI = "LSAL")

group_full_names <- c(MIC_BIF_BREV = "BBRE = B. breve",
                      MIC_LAC_AMYL = "LAMY = L. amylovorus",
                      MIC_LAC_KEFI = "LKEF = L. kefiri",
                      MIC_LAC_PLAN = "LPLA = L. plantarum",
                      MIC_LAC_SALI = "LSAL = L. salivarius")

group_colours <- c(BBRE = "#1B9E77",
                   LAMY = "#D95F02",
                   LKEF = "#2F6FDB",
                   LPLA = "#E7298A",
                   LSAL = "#7570B3",
                   `Not significant` = "grey90",
                   `Mixed` = "grey50",
                   `Check manually` = "black")

level_legend_text <- paste("<b>Annotation levels:</b>",
                           "● = Level 1: confirmed by retention time, accurate mass, and MS/MS match to an analysed reference standard.",
                           "●● = Level 2: supported by accurate mass and MS/MS similarity, without retention-time confirmation.",
                           "●●● = Level 3: tentative annotation based on accurate mass, with or without retention-time agreement.",
                           sep = "<br>")



plot_tbl <- plot_tbl %>%
  mutate(first_group_abbrev = unname(group_abbrev[first_group]),
         second_group_abbrev = unname(group_abbrev[second_group]),
         
         comparison_label = paste(first_group_abbrev, "vs", second_group_abbrev),
         
         higher_group = case_when(direction == "Higher in first group" ~ first_group_abbrev,
                                  direction == "Higher in second group" ~ second_group_abbrev,
                                  direction == "Not significant" ~ "Not significant",
                                  direction == "Multiple directions" ~ "Mixed",
                                  TRUE ~ "Check manually"),
         
         tile_label = case_when(direction == "Higher in first group" ~ paste(first_group_abbrev, ">", second_group_abbrev),
                                direction == "Higher in second group" ~ paste(second_group_abbrev, ">", first_group_abbrev),
                                direction == "Not significant" ~ "n.s.",
                                direction == "Multiple directions" ~ "mixed",
                                TRUE ~ "check"))



comparison_order <- tibble(comparison = pairwise_direction_cols) %>%
  mutate(first_group = str_extract(comparison, "^MIC_.+(?=_vs_)"),
         second_group = str_extract(comparison, "(?<=_vs_)MIC_.+$"),
         first_group_abbrev = unname(group_abbrev[first_group]),
         second_group_abbrev = unname(group_abbrev[second_group]),
         comparison_label = paste(first_group_abbrev, "vs", second_group_abbrev)) %>%
  pull(comparison_label)

plot_tbl <- plot_tbl %>%
  mutate(comparison_label = factor(comparison_label, levels = comparison_order),
         higher_group = factor(higher_group,
                               levels = c("BBRE", "LAMY", "LKEF", "LPLA", "LSAL",
                                          "Not significant", "Mixed", "Check manually")))

molecule_order <- plot_tbl %>%
  distinct(section, discussion_order, molecule_or_class, annotation_level, annotation_balls) %>%
  arrange(section, discussion_order) %>%
  mutate(molecule_label = paste0(molecule_or_class,
    "<br><span style='font-size:6pt;color:#555555'>",
    annotation_balls,
    "</span>")) %>%
  pull(molecule_label)

molecule_order <- rev(unique(molecule_order))

plot_tbl <- plot_tbl %>%
  mutate(molecule_label = paste0(molecule_or_class,
                                 "<br><span style='font-size:6pt;color:#555555'>",
                                 annotation_balls,
                                 "</span>"),
         molecule_label = factor(molecule_label, levels = molecule_order))


discussion_molecule_plot <- ggplot(plot_tbl,
                                   aes(x = comparison_label, y = molecule_label)) +
  geom_tile(aes(fill = higher_group),
            colour = "white",
            linewidth = 0.6) +
  geom_text(aes(label = tile_label),
            size = 2.4,
            colour = "black",
            fontface = "bold") +
  facet_grid(section ~ .,
             scales = "free_y",
             space = "free_y") +
  scale_fill_manual(values = group_colours,
                    labels = c(
                      BBRE = "BBRE = B. breve higher",
                      LAMY = "LAMY = L. amylovorus higher",
                      LKEF = "LKEF = L. kefiri higher",
                      LPLA = "LPLA = L. plantarum higher",
                      LSAL = "LSAL = L. salivarius higher",
                      `Not significant` = "Not significant",
                      Mixed = "Mixed direction across matched features",
                      `Check manually` = "Check manually"),
                    name = "Direction of significant difference") +
  labs(title = "Pairwise differences for metabolites discussed in the manuscript",
       subtitle = "Each tile reports the group with higher abundance; grey indicates non-significant comparisons",
       x = "Pairwise comparison",
       y = NULL,
       caption = paste0("<b>Abbreviations:</b> ",
                        paste(group_full_names, collapse = "; "),
                        "<br><br>",
                        level_legend_text)) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggtext::element_markdown(size = 9),
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.text.y = element_text(face = "bold", angle = 0, hjust = 0),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        plot.caption = ggtext::element_markdown(hjust = 0,
                                                size = 8,
                                                lineheight = 1.15),
        plot.caption.position = "plot")


ggsave(filename = "discussion_molecules_pairwise_summary.png",
       plot = discussion_molecule_plot,
       width = 13,
       height = 7,
       dpi = 300)
