library(tidyverse)
library(patchwork)


# Functions ---------------------------------------------------------------
theme_mine <- function(){
        theme(axis.title = element_text(size = 18),
              axis.text = ggtext::element_markdown(color = "black", size = 16),
              
              strip.text = ggtext::element_markdown(size = 16, face = "bold", color = "white"),
              strip.background = element_rect(fill = "grey30", color = "black"),
              
              panel.spacing.x = unit(1.5, units = "cm"),
              
              legend.text = ggtext::element_markdown(size = 16),
              legend.title = element_text(size = 18),
              legend.key.spacing.y = unit(0.3, "cm"),
              
              
              text = element_text(color = "black")
        )
}


# Import data -------------------------------------------------------------

metadata_semipolars <- readxl::read_xlsx(path = "cuesta2025_probiotic_appetite/data/metadata.xlsx") %>% 
  mutate(bacteria = ifelse(is.na(bacteria), yes = "Control", no = bacteria),
         bacteria = gsub(x = bacteria, pattern = "[.] ", replacement = "_"),
         bacteria = factor(x = bacteria, 
                           levels = c("Control", 
                                      "B_longum", 
                                      "L_reuteri"),
                           labels = c("Control", 
                                      "*B. longum*<br>APC1472", 
                                      "*L. reuteri*<br>ATCC PTA 6475")),
         
         media = factor(x = media, levels = c("CFSs", "CSSs")),
         time = ifelse(is.na(time), "control", time),
         time = factor(x = time, levels = c("control", "early", "late")),
         full_id = gsub(x = full_id, pattern = "[.] | ", replacement = "_"))



lod_df <- readxl::read_xlsx(path = "cuesta2025_probiotic_appetite/data/semipolar_mbx.xlsx",
                            sheet = "LOD") %>% 
  rename_with(tolower) %>% 
  relocate(name) %>% 
  mutate(name = gsub(x = name, pattern = " ", replacement = "_"))


# level 1 of annotation
semipolar_mbx_1 <- readxl::read_xlsx(path = "cuesta2025_probiotic_appetite/data/semipolar_mbx.xlsx", 
                                   sheet = "Annotation_level_1_R") %>% 
  rename_with(.fn = ~gsub(x = ., pattern = " ", replacement = "_")) %>% 
  mutate(full_id = gsub(x = full_id, pattern = "[.] | ", replacement = "_")) %>% 
  column_to_rownames("full_id")

# level 2a of annotation
semipolar_mbx_2a <- readxl::read_xlsx(path = "cuesta2025_probiotic_appetite/data/semipolar_mbx.xlsx", 
                  sheet = "Annotation_level_2a_R") %>% 
  rename_with(.fn = ~gsub(x = ., pattern = " ", replacement = "_")) %>% 
  mutate(full_id = gsub(x = full_id, pattern = "[.] | ", replacement = "_")) %>% 
  column_to_rownames("full_id")

# merging both levels
semipolar_mbx <- inner_join(x = rownames_to_column(.data = semipolar_mbx_1, var = "dummy"),
                            y = rownames_to_column(.data = semipolar_mbx_2a, var = "dummy"),
                            by = "dummy") %>% 
  column_to_rownames(var = "dummy")



# Cleaning mbx data based on LOD (limit of detection)
clean_semipolar_mbx <- semipolar_mbx %>% 
  rownames_to_column("full_id") %>% 
  as_tibble() %>% 
  pivot_longer(cols = !full_id,
               names_to = "name",
               values_to = "raw_abundance") %>% 
  full_join(x = .,
            y = lod_df,
            by = "name") %>% 
  mutate(consider = ifelse(raw_abundance > lod, yes = TRUE, no = FALSE)) %>% 
  filter(consider) %>% 
  select(full_id, name, raw_abundance) %>% 
  pivot_wider(id_cols = full_id,
              names_from = "name",
              values_from = "raw_abundance",
              values_fill = 0)



clean_semipolar_mbx_clr <- clean_semipolar_mbx %>% 
        column_to_rownames("full_id") %>% 
        t() %>% 
        vegan::decostand(x = .,
                         MARGIN = 2,
                         method = "clr",
                         pseudocount = 2/3)



scfa_mbx <- readxl::read_xlsx(path = "cuesta2025_probiotic_appetite/data/scfa_mbx.xlsx",
                  sheet = "absolute_abundances_R", 
                  na = "<LOD") %>% 
        janitor::clean_names() %>% 
        mutate(full_id = gsub(x = full_id, pattern = "[.] ", replacement = "_"),
               full_id = gsub(x = full_id, pattern = " ", replacement = "_")) %>% 
        mutate(across(where(is.numeric), 
                      .fn = ~ifelse(is.na(.x), yes = 0, no = .x))) %>% 
        column_to_rownames("full_id") %>% 
        t()

scfa_mbx <- scfa_mbx[rowSums(scfa_mbx, na.rm = TRUE) > 0,]

rowSums(scfa_mbx > 0)


scfa_mbx

# remove those metabolites where every sample was below the LoD
clean_scfa_mbx <- scfa_mbx[!is.na(rowSums(scfa_mbx)),]


clean_scfa_mbx_clr <- clean_scfa_mbx %>% 
        vegan::decostand(x = .,
                         method = "clr",
                         MARGIN = 2,
                         pseudocount = 2/3) %>% 
        as.data.frame()

# Fig 2 : PCA mbx of strains & growth media -------------------------------

# PCA
clean_semipolar_prcomp <- clean_semipolar_mbx_clr %>% 
  t() %>% 
  prcomp()

semipolar_pc1 <- round(clean_semipolar_prcomp$sdev[1]^2/sum(clean_semipolar_prcomp$sdev^2), 4) * 100
semipolar_pc2 <- round(clean_semipolar_prcomp$sdev[2]^2/sum(clean_semipolar_prcomp$sdev^2), 4) * 100


# TODO: CHANGE THE SHAPES SO THEY CAN BE FILLED
# MAKE THE LEGEND TO SHOW COLORS OF MEDIA

# THREE WAYS
semipolar_pca_plot_df <- clean_semipolar_prcomp$x[,1:2] %>% 
  as.data.frame() %>% 
  rownames_to_column("full_id") %>% 
  inner_join(x = .,
             y = metadata_semipolars,
             by = "full_id") %>% 
  filter(bacteria != "Control")

semipolar_arrow_pca_df <- semipolar_pca_plot_df %>% 
        filter(!grepl(x = full_id, 
                     pattern = "L_reuteriB_mMRS_EarlyStationary|L_reuteriB_mMRS_LateStationary|B_longumA_mMRS_EarlyStationary|B_longumA_mMRS_LateStationary"))
  

fig2a <- ggplot(data = semipolar_pca_plot_df,
         aes(x = PC1, y = PC2, 
             fill = bacteria, 
             shape = media, 
             group = interaction(media, bacteria, bio_rep))) +
  
  geom_point(alpha = 0.75,
             size = 6,
             color = "black"
             #aes(size = time)
             ) +
        
  geom_path(data = semipolar_arrow_pca_df,
            aes(x = PC1, y = PC2, 
                group = interaction(media, bacteria, bio_rep)),
            arrow = arrow(length = unit(0.5, "cm"), angle = 15, type = "closed")) +      
        
  labs(x = paste0("PC 1: ", semipolar_pc1, "%"),
       y = paste0("PC 2: ", semipolar_pc2, "%"),
       fill = "Bacteria",
       shape = "Supernatant\ntype") +
  
  scale_fill_manual(values = c("#DA1902",
                                "#0278DA")) +
  scale_shape_manual(values=c(21, 23)) +
        
  scale_x_continuous(limits = c(-45, 45),
                     breaks = seq(from = -40, to = 40, by = 20)) +
        
  scale_y_continuous(limits = c(-45, 45),
                     breaks = seq(from = -40, to = 40, by = 20)) +
  
  guides(fill = guide_legend(override.aes = list(shape = 21, 
                                                 size = 5)),
         
         shape = guide_legend(override.aes = list(size = 5,
                                                  fill = "grey30"))) +
  
  theme_bw(base_size = 18) + 
  # theme(axis.text = element_text(color = "black"),
  #       
  #       panel.grid.minor = element_blank(),
  #       
  #       legend.text = ggtext::element_markdown(),
  #       legend.key.spacing.y = unit(0.3, "cm"),
  #       
  #       )
        theme_mine(); fig2a

ggsave(filename = "cuesta2025_probiotic_appetite/outputs/figs/semipolar_pca.svg", 
       plot = fig2a, 
       width = 9, height = 6.5, 
       units = "in")




# PCA
clean_scfa_prcomp <- clean_scfa_mbx_clr %>% 
        t() %>% 
        prcomp()

scfa_pc1 <- round(clean_scfa_prcomp$sdev[1]^2/sum(clean_scfa_prcomp$sdev^2), 4) * 100
scfa_pc2 <- round(clean_scfa_prcomp$sdev[2]^2/sum(clean_scfa_prcomp$sdev^2), 4) * 100


# TODO: CHANGE THE SHAPES SO THEY CAN BE FILLED
# MAKE THE LEGEND TO SHOW COLORS OF MEDIA

# THREE WAYS
scfa_pca_plot_df <- clean_scfa_prcomp$x[,1:2] %>% 
        as.data.frame() %>% 
        rownames_to_column("full_id") %>% 
        inner_join(x = .,
                   y = metadata_semipolars,
                   by = "full_id") %>% 
        filter(bacteria != "Control")


scfa_arrow_pca_df <- scfa_pca_plot_df %>% 
        pivot_wider(id_cols = c(bacteria, media, bio_rep),
                    names_from = "time",
                    values_from = c(PC1, PC2))


fig2b <- ggplot(data = scfa_pca_plot_df,
                aes(x = PC1, y = PC2, 
                    fill = bacteria, 
                    shape = media, 
                    group = interaction(media, bacteria, bio_rep)
                    )) +
        
        geom_point(alpha = 0.75,
                   size = 6
        ) +
        
        geom_segment(data = scfa_arrow_pca_df,
                     aes(x = PC1_early, y = PC2_early,
                         xend = PC1_late, yend = PC2_late,
                         group = interaction(media, bacteria, bio_rep)),
                     arrow = arrow(length = unit(0.5, "cm"), angle = 15, type = "closed")) +
        
        labs(x = paste0("PC 1: ", scfa_pc1, "%"),
             y = paste0("PC 2: ", scfa_pc2, "%"),
             fill = "Bacteria",
             shape = "Supernatant\ntype") +
        
        scale_fill_manual(values = c("#DA1902",
                                              "#0278DA")) +
                                                      scale_shape_manual(values=c(21, 23)) +
        
        
        scale_y_continuous(limits = c(-.6, .6),
                           breaks = seq(from = -.6, to = .6, by = .3)) +
        
        guides(fill = guide_legend(override.aes = list(shape = 21, 
                                                       size = 5)),
               
               shape = guide_legend(override.aes = list(size = 5,
                                                        fill = "grey30"))) +
        
        theme_bw(base_size = 18) + 
        # theme(axis.text = element_text(color = "black"),
        #       
        #       panel.grid.minor = element_blank(),
        #       
        #       legend.text = ggtext::element_markdown(),
        #       legend.key.spacing.y = unit(0.3, "cm"),
        #       
        # )
        theme_mine(); fig2b

ggsave(filename = "cuesta2025_probiotic_appetite/outputs/figs/scfa_pca.svg", 
       plot = fig2b, 
       width = 9, height = 6.5, 
       units = "in")



# PERMANOVA ---------------------------------------------------------------

clean_semipolar_mbx_clr

# reorder metadata
metadata_permanova <- column_to_rownames(metadata_semipolars, var = "full_id")
rownames(metadata_permanova) <- colnames(clean_semipolar_mbx_clr)

#class(clean_semipolar_mbx_clr)
dist(t(clean_semipolar_mbx_clr))

semipolar_permanova <- vegan::adonis2(formula = t(clean_semipolar_mbx_clr) ~ media * bacteria * time,
               data = metadata_permanova,
               permutations = 1e4,
               method = "euclidian")


semipolar_permanova %>% 
        broom::tidy() %>% 
        write_csv(x = .,
                  file = "cuesta2025_probiotic_appetite/outputs/permanova_semipolar.csv")


# SCFA
clean_scfa_mbx_clr

scfa_permanova <- vegan::adonis2(formula = t(clean_scfa_mbx_clr) ~ media * bacteria * time,
                                 data = metadata_permanova,
                                 permutations = 1e4,
                                 method = "euclidian")

scfa_permanova %>% 
        broom::tidy() %>% 
        write_csv(x = .,
                  file = "cuesta2025_probiotic_appetite/outputs/permanova_scfa.csv")






# Fig 3A : GBMs in the strains of interest --------------------------------

# custom functions to read and transform the gbm dataframes
read_gbms_df <- function(file){
        
        read_tsv(file = file) %>% 
                mutate(strain = gsub(x = file,
                                     pattern = "_gbm.*",
                                     replacement = ""),
                       strain = word(string = strain,
                                     sep = "/",
                                     start = 3,
                                     end = 3))
}
reshape_gbm_df <- function(df){
        
        df %>% 
                select(!Value) %>% 
                rename_with(tolower)       
}

gbm_files <- c("cuesta2025_probiotic_appetite/data/b_longum_gbm_coverage.modules", 
               "cuesta2025_probiotic_appetite/data/l_plantarum_gbm_coverage.modules")


gbms_db <- omixerRpm::loadDB("GBMs.v1.0")
gbms_info <- as.data.frame(gbms_db@module.names) %>% 
        rownames_to_column("module") %>% 
        rename(name = "V2")

gbms_info


# Plot pathway coverage if there is SCFA production...
gbms_coverage <- gbm_files %>% 
        purrr::map(.x = .,
                   .f = read_gbms_df) %>% 
        purrr::map(.x = .,
                   .f = reshape_gbm_df) %>% 
        purrr::reduce(rbind) %>% 
        inner_join(x = .,
                   y = gbms_info,
                   by = "module") %>% 
        # rename gbms to remove the parentheses and the text within them
        mutate(name = gsub(x = name, pattern = " \\(.*\\)", replacement = "")) %>% 
        
        # change the names of the strains for plotting
        mutate(strain = ifelse(grepl(x = strain,
                                     pattern = "b_longum"),
                               yes = "*B. longum*<br>APC1472",
                               no = "*L. reuteri*<br>ATCC PTA 6475")) %>% 
        # make hierarchical classification of the gbms and transform them to factors
        mutate(GBMs_group = case_when(grepl(x = name, 
                                            pattern = "Tryptophan|Quinolinic") ~ "Tryptophan metabolism",
                                      grepl(x = name, 
                                            pattern = "Glutamate|GABA") ~ "GABA metabolism",
                                      grepl(x = name, 
                                            pattern = "Acetate|Isovaleric") ~ "SCFA",
                                      T ~ "Others"),
        GBMs_group = factor(x = GBMs_group,
                            levels = c("Tryptophan metabolism",
                                       "GABA metabolism",
                                       "SCFA",
                                       "Others")))


gbms_order <- gbms_coverage %>% 
        arrange(GBMs_group, name) %>% 
        pull(name) %>% 
        unique()        



fig3a <- gbms_coverage %>% 
        mutate(GBMs = factor(x = name, 
                             levels = rev(gbms_order)),
               coverage_cat = case_when(coverage <= .1 ~ 0.1,
                                        coverage > .1 & coverage <= .3 ~ 0.3,
                                        coverage > .3 & coverage <= .5 ~ 0.5,
                                        coverage > .5 & coverage <= .8 ~ 0.8,
                                        coverage > .8 ~ 1),
               coverage_cat = as.character(coverage_cat)) %>% 
        ggplot(aes(x = strain, y = GBMs, fill = coverage_cat)) + 
        geom_tile(color = "black") +
        
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        
        scale_fill_manual(values = c("#FFFFAF", "#C2E699", "#78C679", "#31A354", "#006837")) +

        # scale_fill_manual(values = c("0" = "#FFFFAF", "1" = "white", "2" = "#78C679", "3" = "#238443"),
        #                   #"#FFFFCC" "#C2E699" "#78C679" "#31A354" "#006837"
        #                   limits = as.character(0:3)) +
                
        labs(x = NULL, y = NULL, fill = "Pathway\ncoverage") +
        
        theme_bw(base_size = 18) +
        # theme(axis.text = element_text(color = "black"),
        #       axis.text.x = ggtext::element_markdown(color = "black"),
        #       axis.text.y = element_text(color = "black"),
        #       
        #       axis.ticks = element_blank(),
        #       
        #       panel.grid = element_blank())+
        theme_mine() + 
        theme(panel.grid = element_blank()); fig3a


ggsave(filename = "cuesta2025_probiotic_appetite/outputs/figs/fig3a.svg", 
       plot = fig3a, 
       width = 8, height = 6, 
       units = "in")



# Fig 3B : Semipolar MBX based on GBMs ------------------------------------

metabolites_of_interest <- c("Glutamic_acid",
                             "??-Aminobutyric_acid",
                             "Tryptophan")

abs_vs_clr_semipolar_abundances_df <- clean_semipolar_mbx %>% 
        pivot_longer(cols = !full_id,
                     names_to = "metabolite",
                     values_to = "abs_abundance") %>% 
        mutate(clr_abundance = vegan::decostand(x = abs_abundance, 
                                                method = "clr", 
                                                MARGIN = 2, 
                                                pseudocount = 2/3), 
               .by = full_id) %>% 
        filter(metabolite %in% metabolites_of_interest) %>% 
        inner_join(x = .,
                   y = metadata_semipolars,
                   by = "full_id") %>% 
        #pull(bacteria) %>% unique()
        mutate(x_axis = paste0(bacteria, "<br>", "(", media, ")"),#) %>% pull(x_axis)
               x_axis = factor(x_axis,
                               levels = c("Control<br>(CSSs)",
                                          "*B. longum*<br>APC1472<br>(CSSs)",
                                          "*L. reuteri*<br>ATCC PTA 6475<br>(CSSs)",
                                          "Control<br>(CFSs)",
                                          "*B. longum*<br>APC1472<br>(CFSs)",
                                          "*L. reuteri*<br>ATCC PTA 6475<br>(CFSs)")),
               
               metabolite = case_when(grepl(x = metabolite, pattern = "Glutamic") ~ "Glutamate",
                                      grepl(x = metabolite, pattern = "Tryptophan") ~ metabolite,
                                      grepl(x = metabolite, pattern = "Aminobutyric") ~ "GABA")
        ) %>% 
        
        #pull(media) %>% unique()
        pivot_longer(cols = matches(match = "abundance"),
                     names_to = "abundance_type",
                     values_to = "abundance_value") %>% 
        mutate(abundance_type = factor(x = abundance_type,
                                       levels = c("abs_abundance",
                                                  "clr_abundance"),
                                       labels = c("Absolute abundance",
                                                  "CLR-transformed abundance"))
               )





fig3b <- abs_vs_clr_semipolar_abundances_df %>% 
        # only the CLR abundance of the groups grown in CFSs media will be shown
        filter(abundance_type == "CLR-transformed abundance" & media == "CFSs") %>% 
        ggplot(aes(x = x_axis, y = abundance_value, 
                   #color = bacteria, 
                   fill = bacteria,
                   shape = time)) + 
        
        geom_point(size = 4, alpha = 0.7) +
        
        # geom_vline(xintercept = 3.5, 
        #            color = "black",
        #            linewidth = 1.2,
        #            linetype = "dashed") + 
        
        facet_grid(.~metabolite, scales = "free_y") +
        
        scale_fill_manual(values = c("grey40", "#DA1902", "#0278DA")) +
        scale_shape_manual(values = c(4, 21, 22)) +
        
        guides(fill = guide_legend(override.aes = list(shape = 21, 
                                                       size = 5)),
               
               shape = guide_legend(override.aes = list(size = 5,
                                                        fill = "grey30"))) +
        
        labs(y = "CLR-transformed abundance", x = NULL,
             shape = "Collection\ntime",
             fill = "Bacterial\nstrain") +
        
        theme_bw(base_size = 14) + 
        # theme(axis.text = element_text(size = 14),
        #       axis.title = element_text(size = 18),
        #       
        #       strip.text = element_text(size = 18, face = "bold", color = "white"),
        #       strip.background = element_rect(fill = "grey30", color = "black"),
        #       
        #       legend.text = ggtext::element_markdown(),
        #       legend.title = element_text(size = 14),
        #       legend.key.spacing.y = unit(0.3, "cm"),
        #       
        #       axis.text.x = ggtext::element_markdown(color = "black"),
        #       axis.text.y = element_text(color = "black"),
        #       text = element_text(color = "black")
        # ) + 
        theme_mine(); fig3b

ggsave(filename = "cuesta2025_probiotic_appetite/outputs/figs/semipolar_abundances_fig_3b.svg", 
       plot = fig3b, 
       width = 16, height = 8, 
       units = "in")






figsup1 <- abs_vs_clr_semipolar_abundances_df %>% 
        ggplot(aes(x = x_axis, y = abundance_value, 
                   #color = bacteria, 
                   fill = bacteria,
                   shape = time)) + 
        
        geom_rect(xmin = 0, xmax = 3.48, 
                  ymin = -Inf, ymax = Inf,
                  color = "grey75",
                  fill = "grey75",
                  alpha = .35) +
        
        geom_rect(xmin = 3.48, xmax = Inf, 
                  ymin = -Inf, ymax = Inf,
                  color = "grey90",
                  fill = "grey90",
                  alpha = .15) +
        
        geom_point(size = 4, alpha = 0.7) +
        
        geom_vline(xintercept = 3.5, 
                   color = "black",
                   linewidth = 1.2,
                   linetype = "dashed") + 
        
        facet_grid(abundance_type~metabolite, scales = "free_y") +
        
        scale_fill_manual(values = c("grey40", "#DA1902", "#0278DA")) +
        scale_shape_manual(values = c(4, 21, 22)) +
        
        guides(fill = guide_legend(override.aes = list(shape = 21, 
                                                       size = 5)),
               
               shape = guide_legend(override.aes = list(size = 5,
                                                        fill = "grey30"))) +
        
        
        labs(x = NULL, y = "Abundance",
             shape = "Collection\ntime",
             fill = "Bacterial\nstrain") +
        
        theme_bw(base_size = 11) + 
        
        theme(axis.title = element_text(size = 18),
              axis.text.x = ggtext::element_markdown(color = "black", size = 14, angle = 45, hjust = 1),
              axis.text.y = element_text(color = "black", size = 14),
                      
              strip.text = ggtext::element_markdown(size = 16, face = "bold", color = "white"),
              strip.background = element_rect(fill = "grey30", color = "black"),
              
              panel.spacing.x = unit(1.5, units = "cm"),
              
              legend.text = ggtext::element_markdown(size = 16),
              legend.title = element_text(size = 18),
              legend.key.spacing.y = unit(0.3, "cm"),
              
              
              text = element_text(color = "black")
        ); figsup1

ggsave(filename = "cuesta2025_probiotic_appetite/outputs/figs/supplementary_fig_1.svg", 
       plot = figsup1, 
       width = 18, height = 6, 
       units = "in")


# DAA on metabolites
clean_semipolar_mbx_glms <- clean_semipolar_mbx %>% 
        pivot_longer(cols = !full_id,
                     names_to = "metabolite",
                     values_to = "abs_abundance") %>% 
        mutate(clr_abundance = vegan::decostand(x = abs_abundance, 
                                                method = "clr", 
                                                MARGIN = 2, 
                                                pseudocount = 2/3), 
               .by = full_id) %>% 
        inner_join(x = .,
                   y = metadata_semipolars,
                   by = "full_id") %>% 
        mutate(x_axis = paste0(bacteria, "<br>", "(", media, ")"),
               x_axis = factor(x_axis,
                               levels = c("Control<br>(CSSs)",
                                          "*B. longum*<br>APC1472<br>(CSSs)",
                                          "*L. reuteri*<br>ATCC PTA 6475<br>(CSSs)",
                                          "Control<br>(CFSs)",
                                          "*B. longum*<br>APC1472<br>(CFSs)",
                                          "*L. reuteri*<br>ATCC PTA 6475<br>(CFSs)"))) %>% 
        group_by(metabolite) %>% 
        nest() %>% 
        mutate(lm = purrr::map(.x = data,
                               .f = ~lm(data = .x,
                                        formula = clr_abundance~bacteria * media) %>% 
                                       broom::tidy())) %>% 
        ungroup() %>% 
        unnest(lm) %>% 
        select(!data) %>% 
        filter(!is.na(estimate)) %>% 
        filter(term != "(Intercept)") %>% 
        mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
        arrange(metabolite, p.adj) %>% 
        mutate(term = gsub(x = term, pattern = "\\*|\\. ", ""),
               term = gsub(x = term, pattern = "<br>| ", ""))
        
clean_semipolar_mbx_glms %>% 
        xlsx::write.xlsx(x = .,
                         file = "cuesta2025_probiotic_appetite/outputs/semipolar_glms.xlsx")
        
clean_semipolar_mbx_glms %>%         
        filter(metabolite %in% metabolites_of_interest) %>% 
        #filter(grepl(x = metabolite, pattern = "Amino")) %>% 
        #filter(grepl(x = metabolite, pattern = "Tryptophan"))
        filter(p.adj < 0.05)
        
        
        

# Fig 3C : SCFA MBX based on GBMs -----------------------------------------

abs_vs_clr_scfa_abundances_df <- clean_scfa_mbx %>% 
        t() %>% 
        as.data.frame() %>% 
        rownames_to_column("full_id") %>% 
        pivot_longer(cols = !full_id,
                     names_to = "metabolite",
                     values_to = "abs_abundance") %>% 
        mutate(clr_abundance = vegan::decostand(x = abs_abundance, 
                                                method = "clr", 
                                                MARGIN = 2, 
                                                pseudocount = 2/3),
               .by = full_id) %>% 
        inner_join(x = .,
                   y = metadata_semipolars,
                   by = "full_id") %>% 
        #pull(media) %>% unique()
        pivot_longer(cols = matches(match = "abundance"),
                     names_to = "abundance_type",
                     values_to = "abundance_value") %>% 
        mutate(abundance_type = factor(x = abundance_type,
                                       levels = c("abs_abundance",
                                                  "clr_abundance"),
                                       labels = c("Absolute abundance",
                                                  "CLR-transformed abundance")),
               
               x_axis = paste0(bacteria, "<br>", "(", media, ")"),
               x_axis = factor(x_axis,
                               levels = c("Control<br>(CSSs)",
                                          "*B. longum*<br>APC1472<br>(CSSs)",
                                          "*L. reuteri*<br>ATCC PTA 6475<br>(CSSs)",
                                          "Control<br>(CFSs)",
                                          "*B. longum*<br>APC1472<br>(CFSs)",
                                          "*L. reuteri*<br>ATCC PTA 6475<br>(CFSs)")),
               
               metabolite_clean = case_when(grepl(x = metabolite, pattern = "acetic") ~ "Acetate",
                                      grepl(x = metabolite, pattern = "butanoic") ~ "Butyrate",
                                      grepl(x = metabolite, pattern = "propanoic") ~ "Propionate"
                                      )
               )


abs_vs_clr_scfa_abundances_df %>% 
        # only the CLR abundance of the groups grown in CFSs media will be shown
        filter(abundance_type == "CLR-transformed abundance" & media == "CFSs") %>% 
        pull(metabolite) %>% unique()


fig3c <- abs_vs_clr_scfa_abundances_df %>% 
        # only the CLR abundance of the groups grown in CFSs media will be shown
        filter(abundance_type == "CLR-transformed abundance" & media == "CFSs") %>% 
        filter(grepl(x = metabolite, pattern = "acetic_acid|^propanoic_acid|x3_methylbutanoic_acid")) %>% 
        mutate(metabolite = factor(x = metabolite,
                                   levels = c("acetic_acid",
                                              "propanoic_acid",
                                              "x3_methylbutanoic_acid"),
                                   labels = c("Acetate",
                                              "Propionate",
                                              "Isovaleric acid"))) %>% 
        
        ggplot(aes(x = x_axis, y = abundance_value, 
                   #color = bacteria, 
                   fill = bacteria,
                   shape = time)) + 
        
        geom_point(size = 4, alpha = 0.7) +
        
        
        facet_grid(.~metabolite, scales = "free_y") +
        
        scale_fill_manual(values = c("grey40", "#DA1902", "#0278DA")) +
        scale_shape_manual(values = c(4, 21, 22)) +
        
        guides(fill = guide_legend(override.aes = list(shape = 21,
                                                       size = 5)),
               
               shape = guide_legend(override.aes = list(size = 5,
                                                        fill = "grey30"))) +
        
        labs(y = "CLR-transformed abundance", x = NULL,
             shape = "Collection\ntime",
             fill = "Bacterial\nstrain") +
        
        theme_bw(base_size = 14) + 
        theme_mine(); fig3c

ggsave(filename = "cuesta2025_probiotic_appetite/outputs/figs/scfa_abundances_fig_3c.svg", 
       plot = fig3c, 
       width = 16, height = 8, 
       units = "in")






figsup2 <- abs_vs_clr_scfa_abundances_df %>% 
        mutate(metabolite_clean = gsub(x = metabolite, pattern = "_", replacement = " "),
               metabolite_clean = gsub(x = metabolite_clean, pattern = "x[0-9] ", replacement = "")) %>% 
        ggplot(aes(x = x_axis, y = abundance_value, 
                   #color = bacteria, 
                   fill = bacteria,
                   shape = time)) + 
        
        geom_rect(xmin = 0, xmax = 3.48, 
                  ymin = -Inf, ymax = Inf,
                  color = "grey75",
                  fill = "grey75",
                  alpha = .35) +
        
        geom_rect(xmin = 3.48, xmax = Inf, 
                  ymin = -Inf, ymax = Inf,
                  color = "grey90",
                  fill = "grey90",
                  alpha = .15) +
        
        geom_point(size = 4, alpha = 0.7) +
        
        geom_vline(xintercept = 3.5, 
                   color = "black",
                   linewidth = 1.2,
                   linetype = "dashed") + 
        
        facet_grid(abundance_type~metabolite_clean, scales = "free_y") +
        
        scale_fill_manual(values = c("grey40", "#DA1902", "#0278DA")) +
        scale_shape_manual(values = c(4, 21, 22)) +
        
        guides(fill = guide_legend(override.aes = list(shape = 21, 
                                                       size = 5)),
               
               shape = guide_legend(override.aes = list(size = 5,
                                                        fill = "grey30"))) +
        
        
        labs(y = "Abundance", x = NULL,
             shape = "Collection\ntime",
             fill = "Bacterial\nstrain") +
        
        theme_bw(base_size = 11) + 
        
        theme(axis.title = element_text(size = 18),
              axis.text.x = ggtext::element_markdown(color = "black", size = 11, angle = 45, hjust = 1),
              axis.text.y = element_text(color = "black", size = 14),
              
              strip.text = ggtext::element_markdown(size = 11, face = "bold", color = "white"),
              strip.background = element_rect(fill = "grey30", color = "black"),
              
              panel.spacing.x = unit(1.5, units = "cm"),
              
              legend.text = ggtext::element_markdown(size = 16),
              legend.title = element_text(size = 18),
              legend.key.spacing.y = unit(0.3, "cm"),
              
              
              text = element_text(color = "black")
              ); figsup2


ggsave(filename = "cuesta2025_probiotic_appetite/outputs/figs/supplementary_fig_2.svg", 
       plot = figsup2, 
       width = 38, height = 6, 
       units = "in")



# Fig 3 : Put everything together -----------------------------------------

# (fig3a | (fig3b / fig3c)) + 
#         plot_layout(width = c(1, 15))
# 
# ggsave(filename = "cuesta2025_probiotic_appetite/outputs/figs/fig3.svg",
#        plot = last_plot(),
#        width = 20, height = 10, units = "in")


# Fig 4 : Metabolite changes with annotation ------------------------------

metabolites_annotation <- readxl::read_xlsx(path = "cuesta2025_probiotic_appetite/data/dictionary_metabolite_classification.xlsx", sheet = "Sheet1--reduced") %>% 
        rename(metabolite = feature)

metabolomics_df <- clean_semipolar_mbx %>% 
        pivot_longer(!full_id,
                     names_to = "metabolite",
                     values_to = "abundances") %>% 
        mutate(clr_abundance = vegan::decostand(x = abundances, method = "clr", MARGIN = 2, pseudocount = 2/3),
               .by = full_id) %>% 
        inner_join(x = .,
                   y = metadata_semipolars,
                   by = "full_id") %>% 
        inner_join(x = .,
                   y = metabolites_annotation,
                   by = "metabolite") %>% 
        mutate(x_axis = paste0(media, bacteria, time, bio_rep),
               x_axis = factor(x = x_axis,
                               levels = c("CSSsControlcontrolA",
                                          "CSSsControlcontrolB",
                                          
                                          "CSSs*B. longum*<br>APC1472earlyA",
                                          "CSSs*B. longum*<br>APC1472earlyB",
                                          
                                          "CSSs*B. longum*<br>APC1472lateA",
                                          "CSSs*B. longum*<br>APC1472lateB",
                                          
                                          
                                          "CSSs*L. reuteri*<br>ATCC PTA 6475earlyA",
                                          "CSSs*L. reuteri*<br>ATCC PTA 6475earlyB",
                                          
                                          "CSSs*L. reuteri*<br>ATCC PTA 6475lateA",
                                          "CSSs*L. reuteri*<br>ATCC PTA 6475lateB",
                                          
                                          "CFSsControlcontrolA",
                                          "CFSsControlcontrolB",
                                          
                                          "CFSs*B. longum*<br>APC1472earlyA",
                                          "CFSs*B. longum*<br>APC1472earlyB",
                                          
                                          "CFSs*B. longum*<br>APC1472lateA",
                                          "CFSs*B. longum*<br>APC1472lateB",
                                          
                                          
                                          "CFSs*L. reuteri*<br>ATCC PTA 6475earlyA",
                                          "CFSs*L. reuteri*<br>ATCC PTA 6475earlyB",
                                          
                                          "CFSs*L. reuteri*<br>ATCC PTA 6475lateA",
                                          "CFSs*L. reuteri*<br>ATCC PTA 6475lateB")),
               group_reduced = factor(x = group_reduced, 
                                      levels = c("fatty acids and conjugates",
                                                 "amino acid and derivatives",
                                                 "purines, pyrimidines and derivatives",
                                                 "lipids and lipid-like molecules",
                                                 "biotin and derivatives (vitamins)",
                                                 "carbohydrates and carbohydrate conjugates",
                                                 "others")))


metabolite_order <- metabolomics_df %>% 
        arrange(group_reduced, metabolite) %>% 
        pull(metabolite) %>% 
        unique()


mbx_media_bar <- metabolomics_df %>% 
        mutate(media = factor(x = media, levels = c("CSSs", "CFSs"))) %>% 
        ggplot(aes(x = x_axis, y = 1, fill = media)) + 
        geom_tile() + 
        scale_fill_manual(values = c("CSSs" = "grey70", "CFSs" = "#FECD34")) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        labs(x = NULL, y = NULL, fill = "Growth media") +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              #axis.text = element_text(angle = 45),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              legend.direction = "horizontal",
              legend.byrow = TRUE); mbx_media_bar


mbx_bacteria_bar <- metabolomics_df %>% 
        #mutate(media = factor(x = media, levels = c("CSSs", "CFSs"))) %>% 
        ggplot(aes(x = x_axis, y = 1, fill = bacteria)) + 
        geom_tile() + 
        scale_fill_manual(values = c("Control" = "grey40",
                                     "*B. longum*<br>APC1472" = "#DA1902",
                                     "*L. reuteri*<br>ATCC PTA 6475" = "#0278DA")) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        labs(x = NULL, y = NULL, fill = "Bacteria") +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.title = element_text(size = 16),
              legend.text = ggtext::element_markdown(size = 16),
              legend.direction = "horizontal",
              legend.byrow = TRUE); mbx_bacteria_bar

group_reduced_colors <- c(RColorBrewer::brewer.pal(n = 6, name = "Dark2"), "grey50")

mbx_metabolite_annotation_bar <- metabolomics_df %>% 
        mutate(metabolite = factor(x = metabolite, levels = rev(metabolite_order))) %>% 
        ggplot(aes(x = 1, y = metabolite, fill = group_reduced)) + 
        geom_tile() + 
        # scale_fill_manual(values = c("Control" = "grey40",
        #                              "*B. longum*<br>APC1472" = "#DA1902",
        #                              "*L. reuteri*<br>ATCC PTA 6475" = "#0278DA")) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        #scale_fill_brewer(palette = "Dark2") +
        scale_fill_manual(values = group_reduced_colors) +
        labs(x = NULL, y = NULL, fill = "Metabolites group") +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20),
              #legend.position = "lef",
              #legend.direction = "horizontal"
              ); mbx_metabolite_annotation_bar


# z scale per media.
mbx_metabolite_abundances <- metabolomics_df %>% 
        mutate(metabolite = factor(x = metabolite, levels = rev(metabolite_order))) %>% 
        #filter(media == "CFSs") %>% 
        ggplot(aes(x = x_axis, y = metabolite, 
                   fill = clr_abundance)) + 
        geom_tile() + 
        
        scale_fill_gradientn(colours = c("#004961",
                                                  "#009BCE",
                                                  "#FFFFBF",
                                                  "#F56400",
                                                  "#9B4101"),
                                                  # high = "darkred",
                             # low = "darkblue",
                             limits = c(-18, 18),
                             breaks = c(-10, 10),
                             labels = c("Low abundance", "High abundance")
                             ) +
        # scale_fill_gradient2(high = "#9B4101", low = "#004961", mid = "#FFFFBF") +
        
        labs(fill = "CLR-abundance",
             y = NULL,
             X = NULL) +
        
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        
        theme_bw() + 
        #theme_mine() +
        theme(
              
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.title = element_blank(),
              legend.text = element_text(size = 16),
              legend.title.position = "top",
              # legend.justification.bottom = 0.5,
              legend.key.width = unit(3, units = "cm")); mbx_metabolite_abundances


plot_blank <- ggplot() +
        theme_void()


(((mbx_media_bar / mbx_bacteria_bar) | plot_blank)
        + plot_layout(heights = c(1, 1, 2),
                      width = c(28, 1),
                      guides = "collect")) /
        
        ((mbx_metabolite_abundances | mbx_metabolite_annotation_bar) + 
                 plot_layout(widths = c(25, 0.5))) + 
        plot_layout(heights = c(1, 10),
                    widths = c(20, 1))

ggsave(filename = "cuesta2025_probiotic_appetite/outputs/figs/fig4.svg",
       plot = last_plot(),
       height = 10, width = 20, units = "in")


# Supplementary fig 1 : Growth curves -------------------------------------

plates_metadata <- readxl::read_xlsx(path = "cuesta2025_probiotic_appetite/data/growth_curves_microplate.xlsx",
                  skip = 4, 
                  n_max = 1
                  ) %>% 
        janitor::clean_names() %>% 
        select(!c(scheduled_end_timestamp, x2, ovf, x4)) %>% 
        t() %>% 
        as.data.frame() %>% 
        rownames_to_column("condition") %>% 
        rename(well = V1) %>% 
        mutate(condition = ifelse(grepl(x = condition, pattern = "^x"),
                                  yes = "media",
                                  no = condition)) %>% 
        relocate(well) %>% 
        filter(condition != "media") %>% 
        mutate(bacteria = case_when(grepl(x = condition, pattern = "l_reuteri") ~ "*L. reuteri* ATCC PTA 6475",
                                    grepl(x = condition, pattern = "b_long") ~ "*B. longum* APC1472",
                                    T ~ "Control"),
               bacteria = factor(x = bacteria, 
                                 levels = c("Control", 
                                            "*B. longum* APC1472",
                                            "*L. reuteri* ATCC PTA 6475"
                                            )),
               seeding_od = as.numeric(paste0("0.", word(string = condition, sep = "_", start = 4, end = 4))),
               replicate = word(string = condition, sep = "_", start = 5, end = 5))



growth_curves_df <- readxl::read_xlsx(path = "cuesta2025_probiotic_appetite/data/growth_curves_microplate.xlsx",
                  skip = 5, 
                  na = "*", 
                  n_max = 587) %>% 
        janitor::clean_names() %>% 
        select(!c(unix_timestamp, air_temp, x3)) %>% 
        pivot_longer(cols = !time_min,
                     values_to = "OD",
                     names_to = "well") %>% 
        mutate(well = toupper(well),
               time_hours = time_min/60) %>% 
        inner_join(x = .,
                   y = plates_metadata,
                   by = "well") %>% 
        filter(seeding_od == 0.01)

supp_fig_growth_curves <- growth_curves_df %>% 
        ggplot(aes(x = time_hours, 
                   y = OD,
                   color = bacteria,
                   #fill = bacteria,
                   group = replicate)) +
        
        geom_line(linewidth = 1.5,
                  show.legend = TRUE) +
        geom_point(show.legend = FALSE,
                   size = 1,
                   shape = 21) + 
        
        facet_wrap(. ~ bacteria) + 
        
        scale_y_continuous(limits = c(0, 1.85),
                           breaks = seq(from = 0, to = 1.8, by = .3),
                           expand = c(0,0)) +
        
        scale_x_continuous(limits = c(0, 30),
                           breaks = seq(from = 0, to = 30, by = 6),
                           expand = c(0,0)) +
        
        scale_color_manual(values = c("*B. longum* APC1472" = "#DA1902",
                                      "*L. reuteri* ATCC PTA 6475" = "#0278DA")) +
        
        labs(x = "Time (hours)", 
             y = expression(OD["600"]),
             color = "Bacterial strain") + 
        
        guides(#color = guide_legend(override.aes = list(shape = 21, size = 5)),
               color = guide_legend(override.aes = list(linewidth = 5))
               ) +
        
        theme_bw(base_size = 13) + 
        # theme(axis.title = element_text(size = 18),
        #       axis.text = ggtext::element_markdown(color = "black", size = 16),
        #       
        #       strip.text = ggtext::element_markdown(size = 16, face = "bold", color = "white"),
        #       strip.background = element_rect(fill = "grey30", color = "black"),
        #       
        #       panel.spacing.x = unit(1.5, units = "cm"),
        #       
        #       legend.text = ggtext::element_markdown(size = 16),
        #       legend.position = "bottom",
        #       legend.title = element_text(size = 18),
        #       legend.key.spacing.y = unit(0.3, "cm"),
        #       
        #       
        #       text = element_text(color = "black")
        #       ) +
        theme_mine() + 
        theme(legend.position = "bottom"); supp_fig_growth_curves

ggsave(file = "cuesta2025_probiotic_appetite/outputs/figs/supplementary_fig3_growth_curves.svg",
       plot = supp_fig_growth_curves,
       width = 12, height = 6, units = "in")



# supplementary tables : Stats on MBX -------------------------------------

# SCFA dataframe & metadata integration
scfa_df <- clean_scfa_mbx_clr %>% 
        rownames_to_column("metabolite") %>% 
        pivot_longer(cols = !metabolite, 
                     names_to = "full_id",
                     values_to = "clr_abundance") %>% 
        inner_join(x = .,
                   y = metadata_semipolars,
                   by = "full_id")


# effect of bacteria X media -- SCFA
scfa_anova_tukeys_results <- scfa_df %>% 
        filter(bacteria != "Control") %>% 
        nest(.by = "metabolite") %>% 
        mutate(
                lm = purrr::map(.x = data,
                               .f = ~lm(formula = clr_abundance ~ bacteria * media,
                                        data = .x) %>%
                                       car::Anova(type = "II") %>%
                                       broom::tidy()),
                
                tukey = purrr::map(.x = data,
                                   .f = ~aov(formula = clr_abundance ~ bacteria * media,
                                        data = .x) %>% 
                                           TukeyHSD() %>% 
                                           broom::tidy()
                                       )
               ) %>% 
        select(!data)

scfa_anova_tukeys_results %>% 
        select(c(metabolite, lm)) %>% 
        unnest(lm) %>% 
        mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>% 
        xlsx::write.xlsx(x = .,
                         file = "cuesta2025_probiotic_appetite/outputs/mbx_stats_bacteria_X_media.xlsx",
                         sheetName = "scfa_anova", append = FALSE)

scfa_anova_tukeys_results %>% 
        select(c(metabolite, tukey)) %>% 
        unnest(tukey) %>% 
        #mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>% 
        xlsx::write.xlsx(x = .,
                         file = "cuesta2025_probiotic_appetite/outputs/mbx_stats_bacteria_X_media.xlsx",
                         sheetName = "scfa_tukeys", append = TRUE)




# effect of bacteria X media -- semipolar metabolites
semipolars_anova_tukeys_results <- metabolomics_df %>% 
  filter(bacteria != "Control") %>% 
  nest(.by = "metabolite") %>% 
  mutate(
    lm = purrr::map(.x = data,
                    .f = ~lm(formula = clr_abundance ~ bacteria * media,
                             data = .x) %>%
                      car::Anova(type = "II") %>%
                      broom::tidy()),
    
    tukey = purrr::map(.x = data,
                       .f = ~aov(formula = clr_abundance ~ bacteria * media,
                                 data = .x) %>% 
                         TukeyHSD() %>% 
                         broom::tidy()
    )
  ) %>% 
  select(!data)

semipolars_anova_tukeys_results %>% 
  select(c(metabolite, lm)) %>% 
  unnest(lm) %>% 
  mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>% 
  #pull(metabolite) %>% unique() %>% sort()
  xlsx::write.xlsx(x = .,
                   file = "cuesta2025_probiotic_appetite/outputs/mbx_stats_bacteria_X_media.xlsx",
                   sheetName = "semipolar_anova", append = TRUE)

semipolars_anova_tukeys_results %>% 
  select(c(metabolite, tukey)) %>% 
  unnest(tukey) %>% 
  #mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>% 
  xlsx::write.xlsx(x = .,
                   file = "cuesta2025_probiotic_appetite/outputs/mbx_stats_bacteria_X_media.xlsx",
                   sheetName = "semipolar_tukeys", append = TRUE)






# effect of bacteria X time -- SCFA
scfa_anova_tukeys_results <- scfa_df %>% 
  filter(bacteria != "Control") %>% 
  nest(.by = "metabolite") %>% 
  mutate(
    lm = purrr::map(.x = data,
                    .f = ~lm(formula = clr_abundance ~ bacteria * time,
                             data = .x) %>%
                      car::Anova(type = "II") %>%
                      broom::tidy()),
    
    tukey = purrr::map(.x = data,
                       .f = ~aov(formula = clr_abundance ~ bacteria * time,
                                 data = .x) %>% 
                         TukeyHSD() %>% 
                         broom::tidy()
    )
  ) %>% 
  select(!data)

scfa_anova_tukeys_results %>% 
  select(c(metabolite, lm)) %>% 
  unnest(lm) %>% 
  mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>% 
  xlsx::write.xlsx(x = .,
                   file = "cuesta2025_probiotic_appetite/outputs/mbx_stats_bacteria_X_time.xlsx",
                   sheetName = "scfa_anova", append = FALSE)

scfa_anova_tukeys_results %>% 
  select(c(metabolite, tukey)) %>% 
  unnest(tukey) %>% 
  #mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>% 
  xlsx::write.xlsx(x = .,
                   file = "cuesta2025_probiotic_appetite/outputs/mbx_stats_bacteria_X_time.xlsx",
                   sheetName = "scfa_tukeys", append = TRUE)



# effect of bacteria X time -- semipolar metabolites
semipolars_anova_tukeys_results <- metabolomics_df %>% 
  filter(bacteria != "Control") %>% 
  nest(.by = "metabolite") %>% 
  mutate(
    lm = purrr::map(.x = data,
                    .f = ~lm(formula = clr_abundance ~ bacteria * time,
                             data = .x) %>%
                      car::Anova(type = "II") %>%
                      broom::tidy()),
    
    tukey = purrr::map(.x = data,
                       .f = ~aov(formula = clr_abundance ~ bacteria * time,
                                 data = .x) %>% 
                         TukeyHSD() %>% 
                         broom::tidy()
    )
  ) %>% 
  select(!data)

semipolars_anova_tukeys_results %>% 
  select(c(metabolite, lm)) %>% 
  unnest(lm) %>% 
  mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>% 
  #pull(metabolite) %>% unique() %>% sort()
  xlsx::write.xlsx(x = .,
                   file = "cuesta2025_probiotic_appetite/outputs/mbx_stats_bacteria_X_time.xlsx",
                   sheetName = "semipolar_anova", append = TRUE)

semipolars_anova_tukeys_results %>% 
  select(c(metabolite, tukey)) %>% 
  unnest(tukey) %>% 
  #mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>% 
  xlsx::write.xlsx(x = .,
                   file = "cuesta2025_probiotic_appetite/outputs/mbx_stats_bacteria_X_time.xlsx",
                   sheetName = "semipolar_tukeys", append = TRUE)
