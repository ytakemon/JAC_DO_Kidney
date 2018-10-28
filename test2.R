mrna <- read.csv("~/Dropbox/JAX/TheAgingKidneyData/ANOVA/mrna.kidney_anova_table.csv")
protein <- read.csv("~/Dropbox/JAX/TheAgingKidneyData/ANOVA/protein.kidney_anova_table.csv")

load("~/Dropbox/JAX/TheAgingKidneyData/RData/DO188b_kidney_noprobs.RData")
samples <- annot.samples


gene <- c("Psmd2")
mrna %>% filter(symbol %in% gene)
#protein %>% filter(symbol == gene)

#plot
samples %>% mutate(
  expr = expr.mrna[,annot.mrna[annot.mrna$symbol == gene,]$id]) %>%
  ggplot(., aes(x = as.factor(Age), y = expr, colour = Sex, by = Sex))+
    geom_line()+
    facet_wrap(~ Age)+
    scale_fill_aaas()

samples %>% mutate(
  expr = expr.protein[,annot.protein[annot.protein$symbol == gene,]$id]) %>%
  ggplot(., aes(x = Sex, y = expr, fill = Sex))+
    geom_boxplot()+
    facet_wrap(~ Age)+
    scale_fill_aaas()

Upheno <- read.delim("~/Dropbox/JAX/TheAgingKidneyData/Phenotype/JAC_CS_urine_chem_v1.txt", sep = "\t") %>%
  filter(Age.Urine.Chem %in% c(6,12,18)) %>%
  arrange(mouse.id) %>%
  mutate(Age.Urine.Chem = as.factor(Age.Urine.Chem))

fit <- Upheno %>% na.omit() %>% aov(log(phs.cr.u) ~ Age.Urine.Chem, data = .)
summary(fit)
TukeyHSD(fit, "Age.Urine.Chem", conf.level = 0.95)

Upheno %>% na.omit() %>%
  ggplot(., aes(x = as.factor(as.character(Age.Urine.Chem)), y = log(phs.cr.u), fill = Age.Urine.Chem))+
    geom_boxplot()



Upheno <- read.delim("~/Dropbox/JAX/TheAgingKidneyData/Phenotype/JAC_CS_urine_chem_v1.txt", sep = "\t") %>%
mutate(mouse.id = as.character(mouse.id))

annot <- read.csv("~/Dropbox/JAX/TheAgingKidneyData/Phenotype/JAC_DO_all_sampleinfo.csv") %>%
filter(Mouse.ID %in% Upheno$mouse.id) %>%
rename(mouse.id = Mouse.ID) %>%
mutate(mouse.id = as.character(mouse.id))
Upheno <- filter(Upheno, mouse.id %in% annot$mouse.id)


upheno <- inner_join(Upheno, annot, by="mouse.id") %>%
  mutate(Age.Urine.Chem = as.factor(Age.Urine.Chem))


fit <- upheno %>%
  filter(phs.cr.u != "NA" & Sex == "F") %>%
  aov(log(cr.u) ~ Age.Urine.Chem, data = .)
summary(fit)
TukeyHSD(fit, "Age.Urine.Chem", conf.level = 0.95)

pdf("~/Desktop/Phos_data_redo.pdf", width = 6, height =4)
upheno %>%
  na.omit() %>%
  ggplot(., aes(x = Age.Urine.Chem, y = log(phs.cr.u), fill = Age.Urine.Chem))+
  geom_boxplot()+
  facet_wrap(~Sex)
dev.off()

# Age by sex interaction examples
gene <- "Rps28"
protein %>% filter(symbol %in% gene)

plot <- samples %>% mutate(
  expr = expr.protein[,annot.protein[annot.protein$symbol == gene,]$id]) %>%
  ggplot(., aes(x = Age, y = expr, colour = Sex))+
    geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
    geom_point(position = position_jitterdodge(0.4)) +
    theme_bw() +
    labs(y = paste(gene, "protein"),
         x = "Age") +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black")) +
    scale_colour_aaas()

pdf(paste0("~/Desktop/temp/",gene,"_example.pdf"), width = 5, height =4)
plot
dev.off()
