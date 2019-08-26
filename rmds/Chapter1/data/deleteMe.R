# setwd("Desktop/teaching/BIOL_3272_2019/bookSlides/Ch1/lecture_slides/data/")

div.col <-  '#beaed4'
mic.col <-  '#fdc086'
library(readr)
library(readxl)
library(tidyverse)
library(plotly)
library(ggthemes)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
yb_theme <-   theme(axis.text.x=element_text(size=12),
                    axis.text.y=element_text(size=12),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=16),
                    legend.title = element_text(size=12),
                    legend.text = element_text(size=10),
                    plot.title = element_text(size=18),
                    strip.text =  element_text(size=18)) 


admixture <- read_csv("admix_div_mic_named.csv") 
admixture %>%
  filter(trustworthy) %>%
  select(-trustworthy) %>%
  rename(nominal_species = sp)%>% 
  mutate(short_id = factor(short_id, levels = c(paste("d",1:22,sep =""), paste("m",1:5,sep=""))))%>%
  gather(key = sp, value = q, -id,-nominal_species ,-short_id) %>%
  ggplot(aes(x = short_id, y = q, fill = sp))+
  geom_bar(stat = "identity",show.legend = FALSE)+ 
  theme_tufte() + 
  theme(axis.ticks = element_line(linetype = "blank"),
        axis.text.x = element_text(face = "bold",colour = rep(c(div.col, mic.col), times = c(20,5))))+
  scale_fill_manual(values = c(div.col,mic.col))+
  scale_y_continuous(expand = c(0,0))+
  xlab("")+
  annotate(geom = "text",x = 9, y = -.1, label ="H. divaricatus", size = 4.5, color = div.col, fontface="italic")+
  annotate(geom = "text",x = 23, y = -.1, label ="H. microcephalus", size = 4.5, color = mic.col, fontface="italic")+
  coord_cartesian(ylim = c(0, 1), # This focuses the x-axis on the range of interest
                  clip = 'off') +   # This keeps the labels from disappearing
  theme(plot.margin = unit(c(1,1,1,1), "lines")) + # This widens the right margin
  ggtitle("ADMIXTURE: H. divaricatus H. microcephalus", 
          subtitle = "Whole transcriptomes of GRIN samples")+
  ylab("Q (admixture proportion)")



putative.hybrids <- admixture %>% filter(mic > 0.001 & div > 0.001, trustworthy)


transcriptome.pi <- read_excel(path = "sunflowertranscriptome.xlsx") %>%
  mutate(sp1 = case_when(str_detect(ind_1,"664603|547197")  ~ "gross",
                         str_detect(ind_1,"503211")  ~ "micro",
                         str_detect(ind_1 ,"HA"  )       ~ "annuus",
                         str_detect(ind_1 ,"649973|664604|664605|664645|503214|503215|435675|503209")~"div"))%>%
  mutate(sp2 = case_when(str_detect(ind_2,"664603|547197")  ~ "gross",
                         str_detect(ind_2,"503211")  ~ "micro",
                         str_detect(ind_2 ,"HA"  )       ~ "annuus",
                         str_detect(ind_2 ,"649973|664604|664605|664645|503214|503215|435675|503209")~"div"))%>%
  mutate(comparison = factor(case_when( (sp1 == "annuus" & sp2 == "annuus")~ "an.",
                                 (sp1 == "div" & sp2 == "div")~ "div.",
                                 (sp1 == "gross" & sp2 == "gross")~ "gross.",
                                 (sp1 == "micro" & sp2 == "micro")~ "micro.",
                                 ((sp1 == "annuus" & sp2 == "div") | (sp2 == "annuus" & sp1 == "div"))~ "an.\ndiv",
                                 ((sp1 == "annuus" & sp2 == "gross") | (sp2 == "annuus" & sp1 == "gross"))~ "an.\ngross",
                                 ((sp1 == "annuus" & sp2 == "micro") | (sp2 == "annuus" & sp1 == "micro"))~ "an.s\nmicro",
                                 ((sp1 == "gross" & sp2 == "div") | (sp2 == "gross" & sp1 == "div"))~ "div\ngross",
                                 ((sp1 == "micro" & sp2 == "div") | (sp2 == "micro" & sp1 == "div"))~ "div\nmicro",
                                 ((sp1 == "micro" & sp2 == "gross") | (sp2 == "micro" & sp1 == "gross"))~ "micro\ngross"
                                 )))%>%
  mutate(comparison = fct_reorder(.f = comparison, .x = dxy, .fun  = mean))%>%
  mutate(accession_1 = str_remove(ind_1,pattern = "[:alpha:].*"),
         accession_2 = str_remove(ind_2,pattern = "[:alpha:].*"),
         same_accession = accession_1 == accession_2)%>%
  filter(sp1 != "annuus" & sp2 != "annuus") %>%
  filter(ind_1 != ind_2) %>%
  filter(  !(str_detect(string = ind_1,"664604")  | str_detect(string = ind_1,"664605") |str_detect(string = ind_2,"664604")  | str_detect(string = ind_2,"664605")  )   )%>%
  mutate(admx =str_detect(ind_1,"503215| 503213|503214A")| str_detect(ind_2,"503215| 503213|503214A")) %>%
  mutate(putative_hybrid = case_when(
    ind_1 == "503214A"  | ind_2 ==  "503214A" ~ "d22 (q = .604)",
    ind_1 == 503215  | ind_2 ==  503215 ~ "d21 (q = .916)",
    ind_1 == "547197CGRO" |ind_2 ==  "547197CGRO" ~ "gross1"
    )) %>%
  mutate(putative_hybrid =ifelse(is.na(putative_hybrid ),"none",putative_hybrid))%>%
  mutate(putative_hybrid = factor(putative_hybrid, levels = c("d22 (q = .604)","d21 (q = .916)","gross1","none")))%>%  
  mutate(div_same_accession = same_accession & sp1=="div" & sp2 =="div")%>%
  mutate(`With putative\nhybrid` = putative_hybrid)

0.00482

ggplot( transcriptome.pi ,aes(x = comparison, y = dxy, color =`With putative\nhybrid`, shape = `With putative\nhybrid`,alpha = div_same_accession, label = paste(ind_1,ind_2)))+
  geom_jitter(height = 0,width = .2)+
  theme_light()+
  xlab("")+
  ylim(c(0,0.015))+
  scale_alpha_manual(values = c(1,.3))+
  scale_color_manual(values = c(gg_color_hue(3),"#a7adb2" ))+
  labs(x = NULL, y = expression(D[xy]),alpha = "Same accession\n& divaricatus") +
  ggtitle("Pairwise sequence divegence within &\nbetween perrenial sunflower species")+
  theme_tufte()+
  theme(axis.ticks.x = element_blank())+
  yb_theme


ggplotly(a)

#weirdo  664604   664605
mutate(admx =str_detect(ind_1,"503215| 503213|503214A")| str_detect(ind_2,"503215| 503213|503214A"))

grep(503215
     503213
     503214A)

547197



transcriptome.pi %>% 
  filter(!div_same_accession)%>%
  group_by(putative_hybrid, comparison)%>%
  summarise(mean(dxy))
