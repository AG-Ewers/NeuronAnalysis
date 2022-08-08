
library(tidyverse)
library(glue)
library(cowplot)

### Load manual data -----
dirs <- list.dirs(recursive = F)
dirs <- dirs[-1] # remove .Rproj



manual_data <- read_delim(file.path("All_images_unique_name", "Random_Manual_tracings.txt"), delim = "\t")


 ### Load software data ----
 
software_data <- read_delim(file.path("Results_software", "Results_all2661.txt"), delim = " ", trim_ws = T) %>% 
  select(!starts_with("Ratio") & !last_col()) %>% 
  extract(NeuronIndex, into = "NeuronIndex", regex = "SegmentedImage_(.+)") 

 # adjust column names to be comparable 
 varnames <- read_delim("Variable_names.txt", delim = "\t")
 
 name_match = match(names(software_data), varnames$Software)
 names(software_data)[na.omit(name_match)] <- varnames$Manual[!is.na(name_match)]


 ## calculate missing values in software
 software_data <- software_data %>% 
   rowwise() %>% 
   mutate(
     total_branch_length = sum(total_primary_length, total_secondary_length, TertBranchLen),
     total_axon = sum(primary_axon, total_branch_length),
     Total_branches = sum(Number_primary, Number_secondary, Number_tertiary),
     Total_branches_by_total_length = Total_branches / total_axon,
     total_dendrite = total_dendrite + DendBranchLen,
     mean_dendrite = total_dendrite / Number_dendrites,
     Branches_per_dendrite = Number_dendrite_branches / Number_dendrites,
     axon_by_dendrite = total_axon / (total_dendrite + 1),
     Branchpoints_by_total_axon = BranchPointsNum / total_axon,
     Tool = "Software",
     unique_name = NeuronIndex
     ) 
  
 # merge software and automatic
 total <- full_join(manual_data, software_data )

 total_wide <- pivot_wider(total, names_from = "Tool", 
                           values_from = c(total_axon:Total_branches_by_total_length, total_dendrite:axon_by_dendrite), 
                           id_cols = c("unique_name")) %>% 
   filter(total_axon_Software != is.na(total_axon_Software)) 
 


 ### Correlation  of important measures ----
 
 ## Calculate correlation, return only R²
 
 calculate_r2 <- function(manual, software) {

   model <- lm(total_wide[[software]] ~ total_wide[[manual]], data = total_wide)
   broom::glance(model)$adj.r.squared

   
 }

 calculate_pval <- function(readout) {
   fulldata <- total %>% add_count(unique_name) %>% filter(n > 1)
   
   model <- lm(fulldata[[readout]] ~ Tool, data = fulldata)
   broom::glance(model)$p.value
 }
 
 
 ## plot correlation
 
 my_theme <- function() {
   theme_minimal() + theme(text = element_text(size = 7), legend.position = "none", 
                           plot.title = element_text(hjust = 0.5, face = "bold"),
                           plot.subtitle = element_text(hjust = 0.5))
 }
 theme_set(my_theme())
 
 

 
 
plot_correlations <- function(readout) {
  manual <- glue("{readout}_Manual")
  software <- glue("{readout}_Software")
  max_val <- total %>% filter(!str_detect(unique_name, "DIV7")) %>% select(readout) %>% max(na.rm = T)
  min_val <- total %>% filter(!str_detect(unique_name, "DIV7")) %>% select(readout) %>% min(na.rm = T) 
  #print(min_val)
  R2 <- calculate_r2(software, manual) 
  pval <- calculate_pval(readout)
  
  total_wide %>% 
   ggplot(aes_(x = as.name(manual), y = as.name(software))) + 
   geom_point() +
  geom_abline(intercept = 0, slope = 1) +
    geom_smooth(method = "lm") +
    
    labs(x = "Manual", y = "Software", title = readout,
         subtitle = glue("Correlation to manual data R² = {format(R2, digits = 2)}\nDifference to manual p = {format(pval, digits = 2 )}")) +
    coord_fixed(xlim = c(min_val,max_val),  ylim = c(min_val,max_val))
}
 
# Lengths 
 A <- plot_correlations("total_neurite_length") 
 B <- plot_correlations("total_axon") 
 C <- plot_correlations("total_dendrite") 
 D <- plot_correlations("axon_by_dendrite") + scale_x_log10() + scale_y_log10()
 
# Axon Branches
 E <- plot_correlations("Total_branches")
 EF <-plot_correlations("Total_branches_by_total_length") 
 
 # Dendrite Branches
  G <- plot_correlations("Number_dendrite_branches")
  H <- plot_correlations("Branches_per_dendrite")
  

cor_plot <- plot_grid(A, B, C, D, E, EF, G, H,
          labels = "AUTO", ncol = 4, align = "hv", label_size = 10)

cor_plot 

ggsave("Correlations_plot_all2661.png", cor_plot, device = "png", scale = 1, width = 210, height = 140, units = "mm" )
# ggsave("Correlations_plot.pdf", cor_plot, device = "pdf", scale = 1, width = 210, height = 140, units = "mm" )

