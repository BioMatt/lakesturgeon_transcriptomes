library(tidyverse)
library(fishualize) # Colour palettes based on fish species

# Read in the Excel file of Busco scores
buscos <- readxl::read_excel("../manuscript/transcriptome_summary.xlsx", sheet = "BUSCO_scores")

# Melt the dataframe for ggplot, by tissue
melted_buscos <- reshape2::melt(buscos, id = 'Tissue') %>% 
  filter(variable != "Complete") %>% 
  mutate(Tissue = fct_rev(factor(Tissue, levels = c("Overall", "Brain", "Gill", "Head_Kidney", "Heart", "Liver", "White_Muscle", "Esophagus", "Stomach_Proximal", "Stomach_Distal", "Pyloric_Ceca", "Anterior_Intestine", "Spiral_Valve", "Rectum")))) %>% 
  mutate(BUSCO = (factor(variable, levels = c("Missing", "Fragmented", "Duplicated Complete", "Single Complete", "Complete"))))

# Plot BUSCO scores in a sideways stacked barplot
busco_plot <- ggplot(melted_buscos, aes(x = Tissue, y = value, fill = BUSCO, group = BUSCO)) +
  geom_bar(position = 'stack', stat = 'identity') +
  theme_classic() +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"), expand = expansion(mult = 0.009)) +
  ylab("% BUSCOs") +
  scale_x_discrete(labels = rev(c("Overall", "Brain", "Gill", "Head Kidney", "Heart", "Liver", "White Muscle", "Esophagus", "Glandular Stomach", "Muscular Stomach", "Pyloric Caecum", "Anterior Intestine", "Spiral Valve", "Rectum"))) +
  theme(text=element_text(size = 20)) +
  guides(fill = guide_legend(reverse=T)) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1) + 
  coord_flip()
busco_plot

# Pull the lake sturgeon icon from the fishualize package for addition in Inkscape, since I cannot add it above the legend here
lkst_icon <- ggplot() + 
  add_fishape(family = "Acipenseridae", option = "Acipenser_fulvescens") +
  theme_void()

# Save the icon and plot
ggsave(filename = "busco_plot.pdf", plot = busco_plot, dpi = 4000, width = 10, height = 7)
ggsave(filename = "lkst_icon.pdf", plot = lkst_icon, dpi = 4000)
