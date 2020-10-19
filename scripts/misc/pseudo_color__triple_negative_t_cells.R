source("/home/byrne/halo/dev/melanoma/source_all.R")
stDat <- loadStudyData(read_yaml("input/config/study_config.yaml"))

allDat <- stDat$annCells %>% 
          filter(FOV_ID %in% c("1_1_07", "1_2_25", "4_1_04", "4_4_03")) %>%
          mutate(X = (XMin + XMax)/2, Y = -(YMin + YMax)/2)

tcells <- allDat %>%
          #filter(Cell_type == "T cell") %>%
          filter(Subtype == "CD8+ T cell") %>%
          mutate(TripleNeg = ifelse(grepl("PD1|TIM3|LAG3", PositiveMarkers), 0, 1)) %>%
          group_by(FOV_ID, Lesion_response) %>%
          mutate(TCD8cellCount = n(), TripleNegCount = sum(TripleNeg), Fraction = TripleNegCount/TCD8cellCount)

print(tcells)

counts <- tcells %>%
          select(FOV_ID, Lesion_response, TcellCount, TripleNegCount, Fraction) %>%
          unique() %>%
          arrange(FOV_ID)

exh <- tcells %>% filter(TripleNeg == 0)
tn <- tcells %>% filter(TripleNeg == 1)       


p <- ggplot() +
     #geom_point(allDat, mapping = aes(x = X, y = Y), color = "#383838", size = 0.5) +
     geom_point(exh, mapping = aes(x = X, y = Y), color = "#CC0066", size = 0.1) +
     geom_point(tn, mapping = aes(x = X, y = Y), color = "white", size = 0.1) +
     facet_grid(~FOV_ID) +
     theme(panel.background = element_rect(fill = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks.x = element_blank(), 
           axis.ticks.y = element_blank())

pdf("cd8_tcells_white_tripleNeg__non-tripNeg_purple.pdf", height = 2.4, width = 11)
print(p)
dev.off()


