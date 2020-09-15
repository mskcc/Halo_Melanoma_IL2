source("/home/byrne/halo/dev/melanoma/source_all.R")

### TEMP
source("/home/byrne/halo/dev/haloAnalysis/R/condition_detail_figure.R")

library(argparse)

##### GET USER INPUT #####
if(!interactive()){
    ap <- ArgumentParser()
    ap$add_argument("-m", "--manifest", type="character",
                    help = "study manifest")
    ap$add_argument("-sc", "--statistics_config_file", type="character",
                    help = "YAML file containing config for this single analysis")
    ap$add_argument("-fc", "--figure_config_file", type="character",
                    help = "YAML file containing config for figure generation; must include at minimum, 1) list 'comparisons' describing comparisons to be visualized, 2) list 'facet_order' describing how to sort groups of rows, 3) list 'calculations' describing calculations to be visualized (fractions and/or densities), and 4) list 'bio_filter' describing how to choose cell population fractions or densities to be plotted")
    ap$add_argument("-q", "--allQuestions", type="character",
                    help = "XLSX file describing how to subset data in order answer allQuestions statistically")
    ap$add_argument("-c", "--conditions", type="character",
                    help = "XLSX file describing which populations/subpopulations to compare between sample groups")
    args <- ap$parse_args()
} else {
    args <- list(manifest = "input/config/study_config.yaml", 
                 statistics_config_file = "input/config/stats_config.yaml", 
                 figure_config_file = "input/config/stats_detail_figure_config.yaml", 
                 data_out_file = "Fig3.xlsx")
}
###########################
numConds <- function(dat){
    length(unique(dat$`Cell State ID`))
}

log_filtering <- function(indiv, signif, bio, hmp, min.subpop, min.pop, max.fdr, min.diff){

    log_debug(paste("        [FILTER 0] TOTAL NUMBER OF CONDITIONS:", numConds(indiv)))

    notSignif <- numConds(indiv) - numConds(signif)
    log_debug(paste("        [FILTER 1] removed overall conditions with no significance (zero lesion level FDR < 0.05):", notSignif))

    failedBio <- numConds(signif) - numConds(bio)
    log_debug(paste("        [FILTER 2] biological filters:"))
    log_debug(paste("                     min overall subpopulation count = ", min.subpop))
    log_debug(paste("                     min overall population count    = ", min.pop))
    log_debug(paste("                     min overall median difference   = ", min.diff))
    log_debug(paste("                     FDR less than                   = ", max.fdr))
    log_debug(paste("                   conditions removed:", failedBio))

    log_debug(paste("        HMP calculated for [",numConds(hmp), "] conditions"))

    filtHMP <- hmp %>% filter(FDR < 0.05)
    failedHMP <- numConds(bio) - numConds(filtHMP)
    log_debug(paste("        [FILTER 3] removed conditions with adjusted HMP either NA or >= 0.05:", failedHMP)) 

    totalRem <- numConds(indiv) - numConds(filtHMP)
    log_debug(paste("        [TOTAL CONDITIONS REMOVED FROM THIS SET OF COMPARISONS]", totalRem))
    #finalNum <- numConds(filtHMP)
    #log_debug(paste("        [CONDITIONS REMAINING]", finalNum))

}

###### CONFIGURATION & LOAD INPUT DATA ######
cfg  <- resolveConfig(read_yaml(args$manifest), 
                      read_yaml(args$statistics_config_file), 
                      read_yaml(args$figure_config_file), 
                      args)

calcUnit      <- "FOV_ID"
outDir        <- cfg$statistics_detail_dir
statsFiles    <- getFiles(path = cfg$statistics_tables_dir, pattern = ".xlsx")
sheet         <- "all_fractions_and_densities"

facetOrder    <- cfg$facet_order
popOrder      <- cfg$facet_order$Tag
comps         <- cfg$comparisons
calcInfo      <- cfg$calculations

stDat <- loadStudyData(cfg, questions = T, conditions = T)

idMap       <- getIDmap(stDat$allQuestions,
                        qNames = c(comps[[1]]$indiv_question_numbers, comps[[2]]$indiv_question_numbers), 
                        qPre = c(comps[[1]]$indiv_question_prefix, comps[[2]]$indiv_question_prefix))
statUnit    <- "Question"
###############################################

### PREP DATA
###### - determine statuses of individual comparisons vs overall comparisons (same or opposite direction)
###### - calculate harmonic pvalues
###### - filter for conditions where at least one harmonic pval < 0.05 
###### - format for final table
finalDat <- list("fractions" = tibble(`Cell State ID` = seq(1:509)),
                 "densities" = tibble(`Cell State ID` = seq(1:509)))
allIndivDat <- list("fractions" = tibble(), "densities" = tibble())

for(ct in names(calcInfo)){
    log_debug(ct)

    for(x in names(comps)){
        cmp <- comps[[x]]
        qu  <- cmp$question_number
        log_debug(paste0("  ", cmp$name))

        quest     <- stDat$allQuestions[[qu]]
        sGrps     <- getQuestionGroups(quest, stDat$sampAnn)

        statsFile <- statsFiles[grepl(paste0(qu, ".xlsx"), statsFiles)]
        log_debug(paste("      Cohort level stats file:", statsFile))
        log_debug(paste("      Reading cohort level stats from sheet:", sheet))

        c2p       <- stDat$conds %>% filter(AnalysisType == "general") %>% pull(`Cell State ID`) %>% as.numeric()

        if(is.null(cmp$indiv_question_numbers)){
            indivFiles <- getFiles(path = cfg$statistics_tables_dir, pattern = cmp$indiv_question_pattern)
        } else {
            indivFiles <- file.path(cfg$statistics_tables_dir, paste0(cmp$indiv_question_numbers, ".xlsx"))
        }
        log_debug("      Individual level stats files:")
        for(ind in indivFiles){
            log_debug(paste("        ", ind))            
        }       


        ## LESION COMPARISON TO SELECT CONDITIONS TO BE PLOTTED
        ##    create a table containing overall data AND individual data, and how they compare
        log_debug("      Formatting cohort vs individual comparisons...")
        indivDat <- formatForOverallvsIndiv(statsFile, indivFiles, c2p, 
                                            dataCol = calcInfo[[ct]]$calc_column,
                                            effectCol = calcInfo[[ct]]$effect_column, 
                                            statUnit = statUnit)
        allIndivDat[[ct]] <- allIndivDat[[ct]] %>% 
                             bind_rows(indivDat %>% mutate(CompGroup = cmp$name))

        numSame <- indivDat %>% 
                   group_by(`Cell State ID`, Status) %>% 
                   summarize(Num = n()) %>% 
                   spread(Status, Num, fill = 0) %>% 
                   gather(2:ncol(.), key = "Status", value = "NumSame") %>% 
                   filter(grepl("same", Status)) %>% 
                   group_by(`Cell State ID`) %>% 
                   summarize(!!as.name(paste(cmp$name, "NumSame")) := sum(NumSame))

        log_debug("      Filtering data for condition selection...")
        ## remove conditions with zero significant lesion level results
        noSignif <- condsNoSignif(indivDat, statUnit, "signif")
        signif <- indivDat %>% filter(!`Cell State ID` %in% noSignif)

        ## apply "biological" filters
        bio <- signif %>% 
               bioFilter(calc = ct,
                         subpop.col       = 'OverallSubpopulationCount',
                         pop.col          = 'OverallPopulationCount',
                         fdr.col          = 'OverallFDR', 
                         median.cols      = c('OverallGroup1Median', 'OverallGroup2Median'), 
                         max.fdr          = cfg$bio_filter$fdr_cutoff,
                         min.subpop.count = cfg$bio_filter$minimum_subpopulation_count, 
                         min.pop.count    = cfg$bio_filter$minimum_population_count,
                         min.median.diff  = cfg$bio_filter$minimum_median_difference)
     
        ## calculate HMP for all remaining conditions 
        hmp <- bio %>% getHMP("IndivPval", na.rm = FALSE)

        ## log filtering, including adjusted HMP which will actually happen later
        log_filtering(indivDat, signif, bio, hmp, 
                      cfg$bio_filter$minimum_subpopulation_count, 
                      cfg$bio_filter$minimum_population_count, 
                      cfg$bio_filter$fdr_cutoff, 
                      cfg$bio_filter$minimum_median_difference)

        ## compile all data for this comparison
        iSpread <- indivDat %>% 
                   left_join(numSame, by = "Cell State ID") %>%
                   left_join(hmp, by = intersect(names(hmp), names(.))) %>%
                   left_join(read.xlsx(statsFile, sheetName = sheet, check.names = F) %>%
                             as_tibble() %>%
                             mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
                             select_at(c("Cell State ID", 
                                         paste(calcInfo[[ct]]$calc_column, "CI.low"), 
                                         paste(calcInfo[[ct]]$calc_column, "CI.high"))),
                             by = "Cell State ID") %>% 
                   spreadIndivDat(cmp$name,  
                                  cmp$group_1_label, 
                                  cmp$group_2_label, 
                                  calcInfo[[ct]]$calc_column, 
                                  calcInfo[[ct]]$effect_abbreviation, 
                                  idMap, statUnit, "Sample_ID")

        finalDat[[ct]] <- finalDat[[ct]] %>% 
                          left_join(iSpread, by = intersect(names(.), names(iSpread))) %>%
                          select(-dplyr::matches("vs Overall"))
    }

    filt <- paste0("filtered_",ct)
    ### FILTER ON HMP
    finalDat[[filt]] <- finalDat[[ct]] %>%
                        filter_at(vars(names(.)[grepl("adjusted HMP", names(.))]), any_vars(. < 0.05)) %>%
                        select(`Cell State ID`, `Cell State`, Population,  
                        dplyr::matches("Median"),
                        dplyr::matches(comps$comp1$name), 
                        dplyr::matches(comps$comp2$name), 
                        everything())
    hmpPassed <- numConds(finalDat[[filt]])

    ### REMOVE CONDITIONS THAT HAVE LESS THAN X% INDIVIDUAL COMPARISONS GOING IN SAME 
    ### DIRECTION AS CORRESPONDING OVERALL COMPARISON
    numIndiv <- grep(calcInfo[[ct]]$effect_abbreviation, names(finalDat[[ct]])) %>% length() - 2
    mnSmFrac <- 0.6 ### to do: MAKE THIS DYNAMIC
    mnSmCount <- round(numIndiv * mnSmFrac)
    nsCols <- c(paste(comps$comp1$name, "NumSame"), paste(comps$comp2$name, "NumSame"))

    finalDat[[filt]] <- 
        finalDat[[filt]] %>%
        mutate(!!as.name(nsCols[1]) := ifelse(is.na(!!as.name(nsCols[1])), 0, !!as.name(nsCols[1])),
               !!as.name(nsCols[2]) := ifelse(is.na(!!as.name(nsCols[2])), 0, !!as.name(nsCols[2])),
               TotalSame = !!as.name(nsCols[1]) + !!as.name(nsCols[2])) %>%
        filter(TotalSame >= mnSmCount)                  
    final <- numConds(finalDat[[filt]])
    removed <- hmpPassed - final
    log_debug(paste("    [FILTER 4] removed conditions with < ", mnSmCount, " indiv comparisons in same direction:", removed))

    log_debug(paste("    [TOTAL] conditions considered for plotting:", final))

}

fileName <- cfg$data_out_file
openxlsx::write.xlsx(finalDat, fileName)

######################################
######### PLOTTING
######################################
cfg$plot_colors <- getPlotColors(plotColors = cfg$plot_colors, 
                                 plotConfig = cfg$plot_config_file)

xOrder <- c("0_3","1_3","4_3","4_4","4_5","6_2","6_3",
            "0_2","2_2","2_4","1_2","3_2","5_2")

for(ct in names(calcInfo)){
    log_debug(paste("Calculation:",ct))
    suppressMessages(
        plotSeparateCohortVsIndividual(finalDat[[paste0("filtered_",ct)]], 
                                       allIndivDat[[ct]],
                                       cfg, statsFiles, sheet, comps, ct, stDat$cellTypes, stDat$conds, 
                                       idMap %>% filter(Lesion_response != "UT"),
                                       facetY     = calcInfo[[ct]]$y_facets, 
                                       effectCol  = calcInfo[[ct]]$effect_column, 
                                       effectAbb  = calcInfo[[ct]]$effect_abbreviation,
                                       calcCol    = calcInfo[[ct]]$calc_column,
                                       facetOrder = facetOrder, 
                                       popOrder   = popOrder, 
                                       orderBy    = calcInfo[[ct]]$effect_column,
                                       xVar       = "Sample_ID",
                                       xOrder     = xOrder
                                      )
    )
}




