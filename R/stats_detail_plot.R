#' Based on statistics results, determine which conditions should be visualized
#' 
#' Determine which conditions/cell states to visualize by filtering statistics by specified criteria
#' including minimum subpopulation and population counts, minimum percentage difference between 
#' median values, and FDR
#'
#' @param  statsRes    tibble of statistics results as read in from XLSX statistics reports
#' @param  calcColumn  calculation type as it appears in the column names of statsRes; default = "Fraction"
#' @param  statsFilter list of filtering criteria including keys:
#'                         minimum_subpopulation_count
#'                         minimum_population_count
#'                         minimum_median_difference
#'                         fdr_cutoff
#' @return  vector of Cell State IDs for conditions/cell states that pass all filters
conditionsToPlot <- function(statsRes, calcColumn = "Fraction", statsFilter = NULL){
    medCols <- names(statsRes)[grepl(paste0("median ", tolower(calcColumn)), names(statsRes))]
    fdrCol  <- paste(calcColumn, "adjusted p.value")

    if(is.null(statsFilter)){
        ## run filter with defaults
        statsRes %>%
        bioFilter(median.cols = medCols,
                  fdr.col     = fdrCol) %>%
        pull(`Cell State ID`)
    } else {
        statsRes %>%
        bioFilter(min.subpop.count = statsFilter$minimum_subpopulation_count,
                  min.pop.count    = statsFilter$minimum_population_count,
                  min.median.diff  = statsFilter$minimum_median_difference,
                  max.fdr          = statsFilter$fdr_cutoff,
                  median.cols      = medCols,
                  fdr.col          = fdrCol) %>%
        pull(`Cell State ID`)
    }
}


#' Extract and format all data needed to generate stats detail figure for a single question
#' 
#' Filter and format condition labels, FOV values, and statistics for a specified set of 
#' cell states/conditions included in a previously-run statistical analysis
#' 
#' @param statsFile        XLSX file containing all results from a statistical analysis
#' @param sampAnn          tibble containing all sample & FOV annotation
#' @param annCells         tibble of all annotated cells
#' @param question         list describing the question/comparison made in the statistical analysis
#' @param conds            tibble describing all conditions/cell states compared in the analysis
#' @param cellTypes        tibble of cell type definitions
#' @param analyses         list with one element that is a tibble of all of the cell states run
#'                         ONLY for the calculation of interest (e.g., "fractions" or "densities")
#' @param markers          vector of marker names; used to distinguish markers from cell classes
#' @param metricsDir       root directory containing all pre-calculated metrics
#' @param sheetName        name of XLSX workbook sheet containing stats results of interest
#' @param calcUnit         calculation unit used in analysis (default: FOV_ID)
#' @param calculation      calculation of interest; must match name of element in analyses list;
#'                         default = "fractions"
#' @param calcColumn       character string; the substring in column names of XLSX results that
#'                         identifies the results for this calculation; default = "Fraction"
#' @param cellStateIDs     vector of Cell State IDs to include in data; if not provided, results
#'                         will automatically chosen by filtering as specified in statsFilters; 
#'                         default = NULL
#' @param statsFilters     list of filters indicating how to filter results when cellStateIDs are
#'                         not provided; Should include values for keys:
#'                             minimum_subpopulation_count
#'                             minimum_population_count
#'                             minimum_median_difference
#'                             fdr_cutoff
#' @param facetOrder       list with one element for each value in facets, where each element is a
#'                         vector of facet values in the order they should appear on the figure
#' @param popOrder         vector of Cell_type/Tag indicating the order the populations should
#'                         appear in the figure within each facet
#' @param orderBy          when provided, conditions of the same population will be ordered by
#'                         this column
#' @param facets           vector of column names from cellTypes tibble, the values of which will
#'                         be grouped in the visualization
#' @param nbhdCounts       neighborhood counts table; to be specified if comparison involved 
#'                         neighborhood counts; default = NULL
#' @param tumorNbhdCells   vector of UUIDs of cells found in a tumor neighborhood
#' 
#' @return list with three elements: labels, fovVals, and stats
getStatsDetailData <- function(statsFile, sampAnn, annCells, question, conds, 
                               cellTypes, analyses, markers, metricsDir, 
                               sheetName = "all_fractions_and_densities",
                               calcUnit = "FOV_ID", calculation = "fractions", 
                               calcColumn = "Fraction", cellStateIDs = NULL,
                               statsFilters = NULL, facetOrder = NULL, popOrder = NULL, 
                               orderBy = "Odds Ratio", facets = NULL, 
                               nbhdCounts = NULL, tumorNbhdCells = NULL){

    ## get all stats results
    res <- read.xlsx(statsFile, sheetName = sheetName, check.names = F) %>% as_tibble()

    ## determine which cell states to plot
    c2p <- cellStateIDs
    if(is.null(c2p)){
        c2p <- conditionsToPlot(res, calcColumn = calcColumn, statsFilter = statsFilters)
    }
    if(is.null(c2p) || length(c2p) == 0){
        log_warn(paste0("Could not determine which conditions to include in the figure."))
        return(NULL)
    }
    
    sGrps <- getQuestionGroups(question, sampAnn, tumorNbhdCells = tumorNbhdCells)

    labs <- formatConditionsForPlotting(conds, c2p, calculation, cellTypes, 
                                        statsFile = statsFile, sheet = sheetName, 
                                        facetY = facets, orderBy = orderBy,
                                        facetOrder = facetOrder, popOrder = popOrder)

    stats <- formatStatsForPlottingEffects(statsFile, calculation, calcColumn, conds,
                                           c2p, sheet = sheetName, facetY = facets,
                                           cellTypes = cellTypes, facetOrder = facetOrder,
                                           popOrder = popOrder, orderBy = orderBy)

    fovDat <- formatFOVdataForPlottingDetail(question, sGrps, annCells, c2p, conds, 
                                             analyses, calculation, calcColumn, markers,
                                             metricsDir, calcUnit = calcUnit, cellTypes = cellTypes,
                                             nbhdCounts = nbhdCounts, tumorNbhdCells = tumorNbhdCells,
                                             facetY = facets, facetOrder = facetOrder,
                                             popOrder = popOrder, orderBy = orderBy,
                                             statsFile = statsFile, sheet = sheetName)
print("finished")
    list(labels = labs, stats = stats, fovVals = fovDat) 
}


#' ggplot theme to use as a base for stats detail plot
#' 
#' Get a standard ggplot theme for a stats detail plot
#' 
#' @param fontsize     default = 12
#' 
#' @return ggplot theme object
statsDetailTheme <- function(fontsize = 12){
    theme(text                = element_text(size = fontsize, color = "black"),#, family = fontFam),
          axis.line.y         = element_blank(),
          axis.ticks.y        = element_blank(),
          axis.text.y         = element_blank(),
          axis.title.y        = element_blank(),
          axis.line.x         = element_line(size = 0.25),
          axis.text.x         = element_text(size = fontsize, color = "black"),
          axis.title.x.top    = element_text(size = fontsize, color = "black"),
          axis.title.x.bottom = element_blank(),
          panel.grid.major    = element_blank(),
          panel.grid.minor    = element_blank(),
          plot.background     = element_rect(fill = "transparent", color = NA),
          legend.key.height   = unit(5,"mm"),
          legend.key.width    = unit(3,"mm"),
          legend.position     = "top",
          legend.direction    = "vertical",
          legend.text         = element_text(size = fontsize, color = "black"),
          plot.margin         = margin(0,0,0,0))
}


#' Format cell type label for figure
#'
#' Given the plain text version of a cell type, generate a bquote
#' object to be used when plotting. Use the cellTypes definition tibble
#' to determine the correct abbreviation and/or subscript formatting.
#'
#' @param ct         plain text version of a single cell type
#' @param cellTypes  tibble mapping all cell type labels including Category, 
#'                   Cell_type, Subtype, Tag, Abbrev and Subscript
#' @param subtype    value from Subtype column of cellTypes that applies to ct
#' 
#' @return a bquote object
getCellTypeLabel <- function(ct, cellTypes, subtype = NULL){

    idCols <- c("Category", "Cell_type", "Subtype", "Tag")
    ctInfo <- cellTypes %>%
              select_at(c(idCols, "Abbrev", "Subscript")) %>%
              filter_at(vars(idCols), any_vars(. %in% c(ct, subtype))) %>%
              unique()

    if(nrow(ctInfo) == 0){ 
        return(ct)
    }  
    if(nrow(ctInfo) > 1){
        if(any(ctInfo$Subtype == "All")){
            ctInfo <- ctInfo %>% filter(Subtype == "All")
            if(nrow(ctInfo) > 1){ return(bquote(.(ct))) }
        } else {
            warning(paste0("Could not resolve label for condition label [", ct, "]"))
            return(bquote(.(ct)))
        }
    }
    if(!is.na(ctInfo$Abbrev) && ctInfo$Abbrev != ""){
        bquote(.(ctInfo$Abbrev)[.(ctInfo$Subscript)])
    } else {
        bquote(.(ct))
    }

}


#' Get a label as an expression for a single cell state
#' 
#' Given a string representing a single cell state or condition, format 
#' by getting appropriate cell type label according to any abbreviations,
#' in cell type definition table, superscripting positive/negative signs for markers,
#' and converting result to an expression object for correct visualization in 
#' a figure
#' 
#' @param labelStr   plain text character string describing cell state
#' @param cellTypes  table of cell type definitions including columns for 
#'                   Category, Cell_type, Subtype, Tag, Abbrev and Subscript
#' @param delim      for fraction conditions, the delimiter used to distinguish
#'                   subpopulation from population
#' 
#' @return expression object to be parsed for labeling in a ggplot object
conditionLabel <- function(labelStr, cellTypes, delim = " | "){

    #if(is.na(labelStr) || labelStr == "" || !grepl(" \\| ", labelStr)){ return(labelStr) }
    if(is.na(labelStr) || labelStr == ""){ return(labelStr) }

    num_den <- unlist(strsplit(labelStr, " \\| "))

    fin <- list()
    for(x in seq(num_den)){
        lstr <- num_den[x]
        parts <- unlist(strsplit(lstr, ","))
        label <- c()
        for(pt in parts){
            if(grepl("\\+$", pt)){
                mrkr <- gsub("\\+$", "", pt)
                fin[[length(fin) + 1]] <- bquote(.(mrkr)^.("+"))
            } else if(grepl("\\-$", pt)){
                mrkr <- gsub("\\-$", "", pt)
                fin[[length(fin) + 1]] <- bquote(.(mrkr)^.("-"))
            } else {
                fin[[length(fin) + 1]] <- getCellTypeLabel(pt, cellTypes)
            }
        }
        if(x == 1 & length(num_den) > 1){
            fin[[length(fin) + 1]] <- delim
        }
    }

    final <- fin[[1]]
    if(length(fin) > 1){
        for(x in 2:length(fin)){
            sep <- ifelse(bquote(" | ") %in% c(fin[[x]], fin[[x - 1]]), "", ",")
            final <- bquote(.(final)*.(sep)*.(fin[[x]]))
        }
    }

    as.expression(final)

}


#' Format data table to include a striped background
#' 
#' Add to a data table columns for alternating colors and 
#' y coordinates of tops and bottoms of horizontal stripes
#' 
#' @param dat     tibble containing data to be visualized, including at 
#'                minimum y-coordinates
#' @param colors  vector of two colors to be alterated in background of plot
#' @param yVar    column name of y-variable; this column does not necessarily
#'                need to be declared numeric, but must contain values that
#'                can be converted to numbers
#'
#' @return dat tibble, with additional columns: 
#'            y2: numeric version of yVar column
#'            ystart: y-coordinates to use as bottoms of stripe rectangles
#'            yend: y-coordinates to use as tops of stripe rectangles
#'            col: colors to use as fill for stripe rectangles
formatForStripedBackground <- function(dat, colors = c("white", "#f0f0f0"), yVar = "y"){
    dat %>%
    mutate(y2 = as.numeric(as.vector(!!as.name(yVar))),
           ystart = y2 - 0.5, yend = y2 + 0.5,
           col = ifelse(as.numeric(labelY) %% 2 == 0, colors[1], colors[2]))
}


#' Add rows to table that will serve as spacer rows between outermost
#' facets
#'
#' Given a data tibble to be used in a facetted figure, add a row that
#' will be located at the bottom of the facet (as long as there is another
#' facet below it) that can be used to draw a blank rectangle serving
#' as a spacer between facets. Data for these rows will be empty except for
#' the column representing the outermost facet, x/y coordinates, and 'col' 
#'
#' @param dat     tibble containing one or more columns to be used as 
#'                row facets and x and y coordinates
#' @param facetY  vector of column names from dat that will be used to
#'                facet figure rows
#' @param xCol    column name for x values
#' @param yCol    column name for y values; must contain values that
#'                can be converted to numbers
#' @param bgColor main background color
#'
#' @return  table formatted to include facet spacer rows
addFacetSpacerRows <- function(dat, facetY, xCol = "Column", yCol = "y", bgColor = "white"){

    if(is.null(facetY)){ return(dat) }

    dat[[yCol]] <- as.numeric(dat[[yCol]])

    usedLvls <- levels(dat[[facetY[1]]])[levels(dat[[facetY[1]]]) %in% dat[[facetY[1]]]]
    last <- usedLvls[length(usedLvls)]
    fcts <- dat %>%
            select_at(c(facetY[1], xCol)) %>%
            unique() %>%
            filter_at(vars(facetY[1]), all_vars(!. == last)) %>%
            mutate(col = bgColor, !!as.name(yCol) := 1, labelY = 0) %>%
            formatForStripedBackground(colors = rep(bgColor, 2))

    if(length(facetY) > 1){
        for(f in facetY[-1]){
            fcts <- fcts %>% mutate(!!as.name(f) := "")
            fcts[[f]] <- factor(fcts[[f]], levels = levels(dat[[f]]))
        }
    } else {
        yCols <- names(dat)[grep("^y", names(dat))]
        for(yc in yCols){ dat[[yc]] <- dat[[yc]] + 1 }
    }

    fcts <- bind_rows(fcts, dat) 

    fcts[[yCol]] <- factor(fcts[[yCol]], levels = sort(unique(fcts[[yCol]])))
    fcts
}



#' Plot results of a study question 
#'
#' Generate a three- or four-panel figure visualizing results of a single 
#' study sample comparison. Panel 1 is a table describing each
#' cell state or condition compared between groups. Panel 2 shows
#' the distribution of values within each sample group used to
#' statistically compare groups. Generally these are values per
#' FOV. This panel includes medians & IQRs for each group. Panel 3
#' shows effect sizes (OR or FC) and confidence intervals. Optionally,
#' Panel 4 can show a heatmap of additional comparisons where each
#' box indicates whether the effect size was in the same direction
#' as the main comparison and if so, which direction that is. Generally,
#' Panel 4 is used to show consistency between patient level and cohort
#' level results.
#'
#' @param dat             tibble of all densities or fractions used in 
#'                        sample group comparisons, with one row per
#'                        cell state/condition and one column per calc
#'                        unit (generally per FOV); this data is what
#'                        is returned from the formatFOVdataForPlottingDetail()
#' @param stats           tibble of formatted statistics; value returned
#'                        from formatStatsForPlottingEffects()
#' @param conditions      table of formatted and indexed conditions
#' @param clrs            list of colors, with one color for each sample group
#' @param groupVar        name of column containing group labels that match 
#'                        names of clrs list
#' @param calcType        type of calculation being compared; 'fractions' or 'densities'
#' @param cellTypes       tibble mapping all cell type labels including Category, 
#'                        Cell_type, Subtype, Tag, Abbrev and Subscript
#' @param xVar            column name of x values
#' @param yVar            column name of y values
#' @param yNudge          a fraction indicating the vertical adjustment of data group
#'                        summaries (median + IQR plot elements)
#' @param facetY          vector of data column name(s) to be used for row facetting; 
#'                        recommended: Cell_type and Tag
#' @param fontsize        numeric; font size
#' @param spacerColor     character string; R color name or hex color value, including '#';
#'                        fill color of rectangle drawn as a spacer between outer facets
#' @param stripeBG        logical indicating whether to add a horizontal striped background,
#'                        with two alternating colors 
#' @param bgColor         main background color (default = white)
#' @param stripeColor     when stripeBG is TRUE, this color will be one of the two alternating
#'                        row colors; default = "#f0f0f0" (very light gray)
#' @param comparisonData  this optional tibble contains data to be displayed as a heatmap
#'                        in panel 4, where each column in the heatmap is a separate comparison
#'                        of the same conditions, generally with subsets of sample groups (e.g.,
#'                        when main comparison is cohort level, this data might be patient-level
#'                        data). This table must be pre-formatted and must contain columns
#'                        for yVar and compCol, where compCol is to be used as the x variable.
#'                        default = NULL 
#' @param compCol         column name in comparisonData tibble that distinguishes comparison 
#'                        sets, to be used as the x variable in the heatmap; default = NULL
#' @param keepLegend      logical indicating whether to display figure legend(s)
#' @param compFill        column name in comparisonData tibble to be used as the fill variable
#' @param rel_widths      vector of relative widths of all panels 
#' 
#' @return a plot_grid object including all figure panels
plotQuestionResults <- function(dat, stats, conditions, clrs, groupVar, calcType = "fractions",
                                xVar = "Fraction", yVar = "Condition", cellTypes = NULL,
                                yNudge=0.1, facetY = "Cell_type", fontsize = 8, spacerColor = "#D8D8D8",
                                stripeBG = TRUE, bgColor = "white", stripeColor = "#f0f0f0",
                                comparisonData = NULL, compCol = NULL, keepLegend = T, compFill = "Value",
                                rel_widths = c(1, 0.7, 0.6)){

    condTbl   <- conditionTable(conditions, cellTypes, facetY = facetY, fontsize = fontsize,
                                stripeBG = stripeBG, bgColor = bgColor, stripeColor = stripeColor,
                                spacerColor = spacerColor, calcType = calcType)

    detPlot   <- plotQuestionDetail(dat, groupVar, clrs, xVar = xVar, yVar = yVar, yNudge = yNudge,
                                    facetY = facetY, fontsize = fontsize,
                                    stripeBG = stripeBG, bgColor = bgColor, stripeColor = stripeColor,
                                    spacerColor = spacerColor, keepLegend = keepLegend)

    statPlot  <- plotEffectSize(stats, 
                                effectCol = ifelse(xVar %in% c("Fraction","MeanFraction"), 
                                                   "Odds Ratio", "Fold Change"),
                                yVar = yVar, 
                                stripeBG = stripeBG, bgColor = bgColor, stripeColor = stripeColor,
                                facetY = facetY, fontsize = fontsize, spacerColor = spacerColor,
                                keepLegend = keepLegend)

    condAndDetList <- matchPanelHeights(condTbl, detPlot)
    detAndStatList <- matchPanelHeights(detPlot, statPlot)

    plotList <- list(condAndDetList[[1]], condAndDetList[[2]], detAndStatList[[2]])

    if(!is.null(comparisonData)){
        ptPlot        <- plotComparisonHeatmap(comparisonData, compCol, yVar = yVar, facetY = facetY, 
                                               spacerColor = spacerColor, fill = compFill, 
                                               calcType = calcType)
        statAndPtList <- matchPanelHeights(statPlot, ptPlot)
        plotList      <- c(plotList, list(statAndPtList[[2]]))
        if(!length(rel_widths) == 4){
            rel_widths    <- c(2.3, 1, 0.8, 1.5) ## this can be made smaller once legend labels change
        }
    }
    plot_grid(plotlist = plotList, align = 'hv', axis = 'bt', nrow = 1, rel_widths = rel_widths)

}


#' Build ggplot table object displaying formatted and facetted cell state/condition
#' labels of statistics being visualized
#' 
#' Format, facet and align labels of cell states/conditions to be plotted in a
#' two-column table of Cell States and Cell State IDs
#' 
#' @param dat             formatted condition labels returned from formatConditionsForPlotting()
#' @param cellTypes       tibble mapping all cell type labels including Category, 
#'                        Cell_type, Subtype, Tag, Abbrev and Subscript 
#' @param facetY          vector of column names from dat that will be used to
#'                        facet figure rows
#' @param stripeBG        logical indicating whether to add a horizontal striped background,
#'                        with two alternating colors 
#' @param bgColor         main background color (default = white)
#' @param stripeColor     when stripeBG is TRUE, this color will be one of the two alternating
#'                        row colors; default = "#f0f0f0" (very light gray)
#' @param fontsize        numeric; font size; default = 8
#' @param spacerColor     character string; R color name or hex color value, including '#'; 
#'                        fill color of rectangle drawn as a spacer between outer facets
#' @param calcType        type of calculation being compared; 'fractions' or 'densities' 
#' @param hideXaxisTop    logical indicating whether to remove X axis from the top of the figure
#' @param hideXaxisBottom logical indicating whether to remove X axis from the bottom of the figure
#' @param stripWidth      unit(); the fixed width value of left facet strips
#' 
#' @return a ggplot table object
conditionTable <- function(dat, cellTypes, facetY = NULL, stripeBG = TRUE, bgColor = "white",
                           stripeColor = "#f0f0f0", fontsize=8, spacerColor = "#D8D8D8",
                           calcType = "fractions", hideXaxisTop = FALSE, hideXaxisBottom = FALSE, 
                           stripWidth = unit(2,"cm")){

    tmp <- dat

    clrs <- c(bgColor, stripeColor)
    tmp <- tmp %>%
           formatForStripedBackground(colors = clrs, yVar = "y") %>%
           addFacetSpacerRows(facetY, xCol = c("Column"), yCol = "y", bgColor = spacerColor) %>%
           mutate(hjust = ifelse(Column == "Cond", 0, 1),
                  x = ifelse(Column == "Cond", 0, 1),
                  Abbrev = ifelse(is.na(Abbrev)|Abbrev == "", as.character(Tag), Abbrev),
                  Subscript = ifelse(is.na(Subscript), "", Subscript)) 
    if(!is.null(facetY)){
        tmp <- tmp %>%
               filter_at(vars(facetY[1]), all_vars(!is.na(.)))
    } 

    clrs <- c(bgColor, stripeColor, spacerColor)
    clrs <- unique(clrs)
    names(clrs) <- clrs

    lbls <- ggplot(tmp) +
            geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ystart, ymax=yend, fill = col), size = 0.25) +
            geom_hline(yintercept = Inf, size = 0.5, color = "darkgray") +
            geom_hline(yintercept = -Inf, size = 0.5, color = "darkgray") +
            geom_vline(xintercept = -Inf, size = .5, color = "darkgray") +
            scale_y_discrete(expand = c(0, 0)) +
            scale_fill_manual(values = clrs, guide=FALSE) +
            scale_color_manual(values = clrs, guide = FALSE) +
            statsDetailTheme(fontsize = fontsize) +
            theme(axis.text.x = element_text(size = fontsize, hjust = 0.5),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  plot.margin = margin(r = -0.4, t = 3, b = 0.5, l = 0.5, "cm"))

   if(!is.null(facetY)){
       if(length(facetY) == 1 && facetY == "Center"){
           lbls <- lbls + 
                   facet_grid(as.formula(paste(facetY, " ~ .")),
                              scales = "free", space = "free", switch = "y")
       } else {
           lbls <- lbls +
                   facet_grid(as.formula(paste(paste(c(facetY,"Abbrev","Subscript"), collapse=" + "), "~ .")),
                              scales = "free", space = "free", switch = "y",
                              labeller = label_bquote(.(Abbrev)[.(Subscript)])) 
       }
       lbls <- lbls +
               theme(panel.spacing = unit(0, "lines"),
                     strip.placement = "outside",
                     strip.text.y.left = element_text(size = fontsize, angle = 0, margin = margin(0,0.75,0,0,"mm")),
                     strip.background = element_rect(fill = "white", color = "white"))
    }

    if(any(is.na(tmp$Value))){
        lbls <- lbls +
                geom_rect(data = tmp %>% filter(is.na(Value)),
                          aes(xmin = -Inf, xmax = Inf, ymin = ystart, ymax = yend),
                          fill = spacerColor)
    }

    ## add labels one facet at a time to prevent them from being added to all facets
    facets <- tmp %>% select(all_of(facetY)) %>% unique() %>% filter(Tag != "")
    for(i in 1:nrow(facets)){
        fct <- facets[i,] %>% left_join(tmp, by = facetY)
        for(rw in 1:nrow(fct)){
            row <- lapply(as.list(fct[rw,]), function(x) as.vector(x))
            lbls <- lbls +
                    geom_text(dat = fct[rw,],
                              aes(x = x, y = y2, hjust = hjust),
                              vjust = 0.5,
                              label = conditionLabel(row$Value, cellTypes))
        }
    }
    title <- ifelse(calcType == "fractions",
                    paste0(paste(rep(" ", 42), collapse = ""), "Subpopulation | Population"),
                    paste0(paste(rep(" ", 16), collapse = ""), "Population"))
    lbls <- lbls + 
            scale_x_continuous(breaks = c(0,1), limits = c(0,1.05), expand = c(0.01,0.01),
                               labels = c(title, "ID    "),
                               sec.axis = dup_axis()) +
            theme(axis.title.x = element_blank(),
                  axis.title.x.top = element_blank(),
                  axis.text.x.bottom = element_blank())

    if(hideXaxisTop){
        lbls <- lbls +
                theme(axis.text.x.top = element_blank(),
                      axis.ticks.x.top = element_blank(),
                      plot.margin = margin(r = -0.4, t = 0.5, b = 0.5, l = 0.5, "cm"))
    }
    if(hideXaxisBottom){ 
        lbls <- lbls + 
                theme(axis.text.x.bottom = element_blank(), 
                      axis.ticks.x.bottom = element_blank()) }
 
    condTbl <- ggplot_gtable(ggplot_build(lbls))
    for(i in which(grepl("strip-l", condTbl$layout$name))){
       condTbl$grobs[[i]]$layout$clip <- "off"
    }
    condTbl$layout$clip[grep("panel",condTbl$layout$name)] <- "off"
    condTbl$widths[4] <- stripWidth      ## fix strip width for vertical plot alignment
    condTbl
}


#' Generate plot showing distribution of values used to compare cell states/conditions
#' between two sample groups.
#' 
#' Plot real fractions or densities, one per calculation unit (generally FOV), and
#' show median and IQR for both groups. Points and summary elements are colored
#' by sample group.
#' 
#' @param points      data tibble formatted by formatFOVdataForPlotting()
#' @param groupVar    name of column identifying sample groups
#' @param clrs        list of colors named by values in groupVar column
#' @param xVar        column name of x values
#' @param yVar        column name of y values
#' @param facetY      vector of data column name(s) to be used for row facetting; 
#'                    recommended: Cell_type and Tag
#' @param stripeBG    logical indicating whether to add a horizontal striped background,
#'                        with two alternating colors 
#' @param bgColor     main background color (default = white)
#' @param stripeColor when stripeBG is TRUE, this color will be one of the two alternating
#'                    row colors; default = "#f0f0f0" (very light gray)
#' @param fontsize    numeric; font size; default = 8
#' @param yNudge      a fraction indicating the vertical adjustment of data group
#'                    summaries (median + IQR plot elements); default = 0.1
#' @param spacerColor character string; R color name or hex color value, including '#';
#'                    fill color of rectangle drawn as a spacer between outer facets
#' @param keepLegend  logical indicating whether to display figure legend(s)
#' 
#' @return ggplot table object
plotQuestionDetail <- function(points, groupVar, clrs, xVar = "Fraction", yVar = "Condition",
                               facetY = "Subtype", yNudge = 0.1, fontsize = 8,
                               stripeBG=T, bgColor = "white", stripeColor = "#f0f0f0",
                               spacerColor = "#D8D8D8", keepLegend = T){

    if(!stripeBG){ stripeColor <- bgColor }

    points <- points %>%
              formatForStripedBackground(colors = c(bgColor, stripeColor), yVar = "y") %>%
              addFacetSpacerRows(facetY, xCol = NULL, yCol = "y", bgColor = spacerColor) %>%
              mutate(smryY = ifelse(!!as.name(groupVar) == levels(points[[groupVar]])[1],
                                    y2 + yNudge, y2 - yNudge))
    if(!is.null(facetY)){
        points <-  points %>%
                   filter_at(vars(facetY[1]), all_vars(!is.na(.)))
    }

    bgColors <- c(bgColor, stripeColor, spacerColor)
    names(bgColors) <- bgColors

    gridColor <- "white"
    gridSize  <- 1

    nr <- length(unique(points$`Cell State ID`))
    legendY <- 1 + (0.1 * (60/nr))

    breaks <- list("Fraction" = c(0, 0.25, 0.5, 0.75, 1),
                   "MeanFraction" = c(0, 0.25, 0.5, 0.75, 1),
                   "Density" = c(0, 1, 10, 100, 1000, 10000))

    labels <- list("Fraction" = c("0","0.25","0.5","0.75","1"),
                   "MeanFraction" = c("0","0.25","0.5","0.75","1"),
                   "Density" = c("0", "1","10","100",parse(text = "10^3"), parse(text="10^4")))

    xTrans     <- ifelse(xVar == "Density", "log1p", "identity")
    xAxisMin   <- ifelse(xVar == "Density", breaks[[xVar]][1], 0)
    xStripeMin <- ifelse(xVar == "Density", breaks[[xVar]][1], -Inf)
    xMax <- breaks[[xVar]][length(breaks[[xVar]])*1.03]

    ## init
    p <- ggplot(points, 
                aes_string(x = xVar, y = "y", color = groupVar, fill = groupVar)) 

    ## colors & axes
    xLbl <- paste(xVar,"Marker Positive Cells")
    if(xVar == "Density"){ xLbl <- bquote(atop(.(xLbl)*phantom(0), "("*cells/mm^2*")")) }
    p <- p +
         scale_y_discrete(expand = c(0,0)) +
         scale_x_continuous(breaks = breaks[[xVar]], labels = labels[[xVar]], limits = c(xAxisMin, xMax),
                            sec.axis = dup_axis(), trans = xTrans, expand = c(0.03, 0.03)) +
         scale_color_manual(name = groupVar, values = c(bgColors, clrs),
                            breaks = names(clrs), labels = names(clrs)) +
         scale_fill_manual(name = groupVar, values = c(bgColors, clrs),
                           breaks = names(clrs), labels = names(clrs)) +
         ylab("") +
         xlab(xLbl)

    ## add background (likely striped)
    p <- p + 
         geom_rect(aes(xmin = xStripeMin, xmax = Inf, ymin = ystart, 
                       ymax = yend, fill = col, color = col),
                   size = 0.25, show.legend = F)

    ## add light border on left side 
    p <- p +  
         geom_vline(xintercept = -Inf, size = 0.5, color = "lightgray") 

    ## add white vertical guide lines
    for(x in seq(breaks[[xVar]])){
        p <- p +
             geom_vline(xintercept = breaks[[xVar]][x], size = gridSize, color = gridColor) 
    }

    ## actual data
    p <- p +
         geom_jitter(height = 0.25, size = 1.25, alpha = 0.30, shape = 19) +
         geom_segment(aes_string(x = "min", xend = "max", y = "smryY", yend = "smryY", color = groupVar),
                      size = 0.6) +
         geom_point(aes_string(x = "Median", y = "smryY", fill = groupVar),
                    color = "black", shape = 23, size = 3, stroke = 0.65) 

    ## top & bottom border lines
    p <- p +
         geom_hline(yintercept = Inf, size = 0.5, color = "darkgray") +
         geom_hline(yintercept = -Inf, size = 0.5, color = "darkgray") 

    ## theme
    p <- p +
         statsDetailTheme(fontsize = fontsize) +
         theme(legend.title = element_blank(),
               legend.position = c(0.5, legendY),
               legend.justification = c(0.5,0),
               panel.background = element_blank(),
               panel.spacing = unit(0,"lines"),
               plot.margin = margin(l = 0, r = -0.8, t=3, b=0.5, unit = "cm")) +
         guides(shape = guide_legend(groupVar), 
                linetype = guide_legend(groupVar), 
                color = guide_legend(groupVar))

    ## add blank spacer rows between facets
    if(any(!is.na(points$labelY) && points$labelY == 0)){
        p <- p +
             geom_rect(data = points %>% filter(labelY == 0),
                       aes(xmin = xStripeMin, xmax = Inf, ymin = ystart, ymax = yend),
                       color = "transparent", fill = spacerColor)
    }
    if(!is.null(facetY)){
        p <- p +
             facet_grid(as.formula(paste(paste0(facetY,collapse="+"), " ~ .")), 
                        scales="free", space="free") 
    }

    if(!keepLegend){ p <- p + theme(legend.position = "none") }

    g1 <- ggplot_gtable(ggplot_build(p))
    g1$layout$clip[grep("panel",g1$layout$name)] <- "off"

    return(g1)

}


#' Plot effect size for a set of conditions/cell states
#'
#' Plot effect size (odds ratios or fold changes) and confidence 
#' intervals for a specific set of cell states or conditions. 
#'
#' @param dat             formatted statistics from formatStatsForPlottingEffects()
#' @param effectCol       column name of effect sizes (e.g.,  "Odds Ratio" or "Fold Change")
#' @param yVar            column of y values
#' @param clrs            named list of colors to use for effect and confidence intervals
#' @param colorBy         column containing names of clrs list
#' @param facetY          vector of data column name(s) to be used for row facetting; 
#'                        recommended: Cell_type and Tag
#' @param stripeBG        logical indicating whether to add a horizontal striped background,
#'                        with two alternating colors  
#' @param bgColor         main background color (default = white)
#' @param stripeColor     when stripeBG is TRUE, this color will be one of the two alternating
#'                        row colors; default = "#f0f0f0" (very light gray) 
#' @param fontsize        numeric; font size; default = 8
#' @param spacerCol       character string; R color name or hex color value, including '#';
#'                        fill color of rectangle drawn as a spacer between outer facets
#' @param keepLegend      logical indicating whether to display figure legend(s)
#' @param hideXaxisTop    logical indicating whether to remove x axis ticks and labels from top
#'                        of panel; default = FALSE
#' @param hideXaxisBottom logical indicating whether to remove x axis ticks and labels from bottom
#'                        of panel; default = FALSE
#' @param hideXtitle      logical indicating whether to remove x title from top of panel 
#'                        (default = FALSE)
#' 
#' @return ggplot table object
plotEffectSize <- function(dat, effectCol="Odds Ratio", yVar="Condition",
                           clrs = NULL, colorBy = "signif", facetY = "Subtype", 
                           stripeBG = TRUE, bgColor = "white", stripeColor = "#f0f0f0",
                           fontsize = 8, spacerColor = "#D8D8D8", keepLegend = T,
                           hideXaxisTop = FALSE, hideXaxisBottom = FALSE, hideXtitle = FALSE){

    if(!stripeBG){ stripeColor <- bgColor }

    resTbl <- dat %>%
              formatForStripedBackground(colors = c(bgColor, stripeColor), yVar = "y") %>%
              addFacetSpacerRows(facetY, xCol = NULL, yCol = "y", bgColor = spacerColor) 
    if(!is.null(facetY)){
        resTbl <- resTbl %>%
                  filter_at(vars(facetY[1]), all_vars(!is.na(.)))
    }

    bgClrs <- c(bgColor, stripeColor, spacerColor)
    names(bgClrs) <- bgClrs

    breaks <- c(0.01, 0.1, 1, 10, 100)
    labels <- c(parse(text = "''<=10^-2"), 
                parse(text = "10^-1"), "1", 
                parse(text = "10^1"), 
                parse(text = "''>=10^2"))

    mn <- breaks[1]
    mx <- breaks[length(breaks)]

    ## set up colors just so we can put significance in the legend; 
    ## they must be three different colors but will be transparent so it 
    ## doesn't matter what they are
    if(is.null(clrs)){
        clrs <- c("***" = "black", "**" = "#808080", "*" = "#808080")
        clrBreaks <- c("***", "**")
        clrLabels <- c("  < 0.001", "  < 0.05")
        scaleName <- "FDR"
    } else {
        clrBreaks <- clrLabels <- names(clrs)
        scaleName <- ""
        first <- which(!is.na(resTbl[[colorBy]]) & resTbl[[colorBy]] == clrLabels[1])
        second <- which(!is.na(resTbl[[colorBy]]) & resTbl[[colorBy]] == clrLabels[2])

        resTbl$y2[first] <- resTbl$y2[first] - 0.01
        resTbl$y2[second] <- resTbl$y2[second] + 0.01
    }

    allClrs <- c(clrs, bgClrs)

    resTbl$CI.low[resTbl$CI.low < breaks[1]] <- breaks[1]
    resTbl$CI.high[resTbl$CI.high > breaks[length(breaks)]] <- breaks[length(breaks)]
    resTbl$signif[resTbl$signif == ""] <- NA

    gridSize <- 1
    gridColor <- "white"

    nr <- length(unique(resTbl$`Cell State ID`))
    legendY <- 1 + (0.035 * (100/nr))

    p <- ggplot(resTbl, aes_string(x = add_ticks(effectCol), y = "y")) +
         geom_rect(aes(xmin = 0, xmax = Inf, ymin = ystart, ymax=yend, 
                       fill = col, color = col),
                       size = 0.25, show.legend = F) +
         geom_vline(xintercept = 0, size = 0.5, color = "lightgray") +
         geom_vline(xintercept = breaks[1], size = gridSize, color = gridColor) +
         geom_vline(xintercept = breaks[2], size = gridSize, color = gridColor) +
         geom_vline(xintercept = breaks[3], size = gridSize, color = gridColor) +
         geom_vline(xintercept = breaks[4], size = gridSize, color = gridColor) +
         geom_vline(xintercept = breaks[5], size = gridSize, color = gridColor) +
         geom_line(aes_string(color = colorBy), alpha=0, size = 0.01) +
         geom_hline(yintercept = Inf, size = .5, color = "darkgray") +
         geom_hline(yintercept = -Inf, size = .5, color = "darkgray") +
         geom_vline(xintercept = Inf, size = .5, color = "darkgray") +
         geom_vline(xintercept = 1, size=0.25) +
         geom_point(aes_string(add_ticks(effectCol), y = "y2", color = colorBy), size=2, show.legend = T) +
         geom_segment(aes_string(x = "CI.low", xend = "CI.high", y = "y2", yend = "y2", color = colorBy), 
                      size = 0.6, show.legend = T) +
         scale_y_discrete(expand = c(0,0)) +
         scale_shape_manual(values = 19) +
         scale_linetype_manual(values = "1") +
         scale_color_manual(name = scaleName, values = allClrs, breaks = clrBreaks, labels = clrLabels) +
         scale_fill_manual(values = allClrs, guide = FALSE) +
         scale_x_continuous(trans = "log", oob = scales::squish, breaks=breaks, labels=labels, 
                            limits=c(mn,mx), expand = c(0.1, 0.1),
                            sec.axis = dup_axis()) +
         ylab("") +
         xlab(paste0(effectCol, " & \n95% CI")) +
         statsDetailTheme(fontsize = fontsize) +
         theme(legend.direction = "vertical",
                 legend.position = c(1, legendY),
                 legend.justification = c(1,0),
                 legend.text = element_text(size = fontsize, hjust = 1),
                 legend.title = element_text(size = fontsize, hjust = 1),
                 legend.key = element_rect(color = "transparent", fill = "white"),
                 legend.key.size = unit(0.5, "cm"),
                 legend.key.width = unit(0.75,"cm"),
                 legend.background = element_rect(color = "transparent", fill = "transparent"),
                 panel.background = element_blank(),
                 panel.spacing = unit(0,"lines"),
                 plot.margin = margin(l = -0.5, r = -1.1, t=3, b=0.5, unit = "cm")) 

    if(colorBy == "signif"){
         p <- p +
              guides(shape = guide_legend("FDR"), 
              linetype = guide_legend("FDR"), 
              color = guide_legend("FDR"))#, fill = NULL)
    }

    ## shade space between facets
    if(any(!is.na(resTbl$labelY) && resTbl$labelY == 0)){
        p <- p +
             geom_rect(data = resTbl %>% filter(labelY == 0),
                       aes(xmin = 0, xmax = Inf, ymin = ystart, ymax = yend),
                       fill = spacerColor)
    }

    if(!is.null(facetY)){
        p <- p +
             facet_grid(as.formula(paste(paste0(facetY,collapse="+"), " ~ .")), scales = "free", space="free") +
             theme(strip.background = element_blank(),
                   strip.text = element_text(color = "transparent"))
    }

    if(!keepLegend){ p <- p + theme(legend.position = "none") }
    if(hideXaxisTop){ 
        p <- p + 
             theme(axis.text.x.top = element_blank(), 
                   axis.ticks.x.top = element_blank(),
                   plot.margin = margin(l = -0.5, r = -1.1, t=0.5, b=0.5, unit = "cm")) 
    }
    if(hideXaxisBottom){ p <- p + theme(axis.text.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank()) }
    if(hideXtitle){ p <- p + theme(axis.title.x = element_blank(), axis.title.x.top = element_blank()) }

    gt <- ggplot_gtable(ggplot_build(p))

    # Code to override clipping
    gt$layout$clip[grep("panel",gt$layout$name)] <- "off"

    gt
}


#' Generate heatmap highlighting the similarities or dissimilarities
#' between multiple comparisons
#'
#' Create a heatmap to accompany a stats detail figure that depicts the level
#' of consistency or lack thereof between multiple sample sets. Generally this
#' plot is used to compare patient level sample comparisons to cohort comparisons.
#' The heatmap is colored according to the direction of the effect size when
#' the patient level comparison is consistent with the cohort level results and
#' is colored gray when the effect size is in the opposite direction of the cohort
#' level result.
#' 
#' @param compDat         tibble containing stats results for multiple sets of
#'                        samples; formatted by [INSERT FUNCTION NAME HERE]
#' @param compCol         column distinguishing comparisons, to be used as x variable
#' @param yVar            column of y values
#' @param facetY          vector of data column name(s) to be used for row facetting; 
#'                        recommended: Cell_type and Tag
#' @param fontsize        numeric; font size; default = 12
#' @param stripeColor     when stripeBG is TRUE, this color will be one of the two alternating
#'                        row colors; default = "#f0f0f0" (very light gray)
#' @param bgColor         main background color (default = white)
#' @param calcType        type of calculation being compared; 'fractions' or 'densities'
#' @param spacerColor     character string; R color name or hex color value, including '#';
#'                        fill color of rectangle drawn as a spacer between outer facets
#' @param keepLegend      logical indicating whether to display figure legend(s)
#' @param fill            name of column containing "up/down/opposite" values to be used 
#'                        for coloring heatmap boxes; colors are red/blue/gray, respectively
#' @param hideXaxisTop    logical indicating whether to remove x axis ticks and labels from top of plot
#' @param hideXaxisBottom logical indicating whether to remove x axis ticks and labels from bottom of plot
#' @param hideXtitle      logical indicating whether to remove x title from plot
#' 
#' @return ggplot table object
plotComparisonHeatmap <- function(compDat, compCol, yVar = "y", facetY = NULL, fontsize = 12, 
                                  stripeColor = "#f0f0f0", bgColor = "white", calcType = "fractions",
                                  spacerColor = "#d8d8d8", keepLegend = T, fill = Value,
                                  hideXaxisTop = FALSE, hideXaxisBottom = FALSE, hideXtitle = FALSE){

    breaks <- as.numeric(levels(compDat$x))
    maxx   <- as.numeric(max(breaks))
    minx   <- as.numeric(min(breaks))

    dat <- compDat %>%
           formatForStripedBackground(colors = c(bgColor, stripeColor), yVar = yVar) %>%
           addFacetSpacerRows(facetY, xCol = "x", yCol = yVar, bgColor = spacerColor)
    dat$x <- as.numeric(dat$x)
    xLabels <- dat %>% 
               select_at(c("x", compCol)) %>% 
               filter(!is.na(x), !is.na(!!as.name(compCol))) %>% 
               unique %>% arrange(x)
    xLbls <- xLabels$x
    names(xLbls) <- xLabels[[compCol]]

    gridColor <- "darkgray"
    tileH <- 0.95
    tileW <- 0.95

    nr <- length(unique(dat$`Cell State ID`))
    legendY <- 1 + (0.1 * (35/nr))
    legendLabels <- c("OR > 1", "OR < 1", "OR opposite direction")
    if(calcType == "densities"){
        legendLabels <- gsub("OR", "FC", legendLabels)
    }

    clrs <- c('up' = "#b2182b", 'down' = '#2166ac', 'opposite' = '#D8DDDF')
    p <- ggplot(dat %>% filter(!is.na(!!as.name(compCol))), 
                aes_string(x = "x", y = "y2")) + 
         geom_tile(data = dat %>% filter(is.na(!!as.name(fill))), fill = "white", color = gridColor, 
                   height = tileH, width = tileW) 

    if(is.numeric(dat[[fill]])){
        dat[[fill]] <- cut(dat[[fill]], breaks = c(-2, 0, 0.01, 0.5, 2, 100, max(dat[[fill]], na.rm=T)))
        dat[[fill]] <- factor(dat[[fill]], levels = rev(levels(dat[[fill]])))

        clrs <- c("#b2182b","#f4a582","#f7f7f7","#92c5de","#2166ac","#B8B8B8")
        names(clrs) <- levels(dat[[fill]])

        p <- p +
             geom_tile(data = dat, 
                       aes(fill = !!as.name(fill)), stat = "identity", size = 0.25, color = gridColor,
                       height = tileH, width = tileW) +
             geom_point(data = dat %>% filter(!!as.name(fill) != "(-2,0]"),
                        aes(shape = signif), color = "black", size = 1) +
             scale_fill_manual(values = clrs,
                               labels = c("OR > 100", "OR > 2", "0.5 < OR < 2", "OR < 0.5", "0R < 0.1", "OR opposite"), drop = TRUE) 
    } else {
        if(!all(levels(dat[[fill]]) %in% c(NA,names(clrs)))){
            stop(paste0("Fill variable can not have more than four categories."))
        }
        dat[[fill]] <- factor(dat[[fill]], levels = names(clrs))
        p <- p +
             geom_tile(data = dat %>% filter(!is.na(!!as.name(fill))),
                       aes(fill = !!as.name(fill)), stat = "identity", size = 0.25, color = gridColor,
                       height = tileH, width = tileW) +
             geom_point(data = dat %>% filter(!!as.name(fill) != "opposite"),
                       aes(shape = signif), stat = "identity", color = "black", size = 1) +
             scale_fill_manual(values = clrs, labels = legendLabels) 
    }

    p <- p +
         geom_segment(data = dat %>% filter(is.na(!!as.name(fill))),
                      aes(x = x - 0.45, xend = x + 0.45, 
                          y = y2 + 0.45, yend = y2 - 0.45), color = "lightgray") +
         geom_segment(data = dat %>% filter(is.na(!!as.name(fill))),
                      aes(x = x + 0.45, xend = x - 0.45, 
                          y = y2 + 0.45, yend = y2 - 0.45), color = "lightgray") +
         xlab(compCol) +
         scale_y_discrete(expand = c(0,0)) +
         scale_x_continuous(breaks = breaks, limits = c(minx - 0.5, maxx + 0.50),
                            expand = c(0,0), sec.axis = dup_axis(), labels = names(xLbls)) +
         scale_shape_manual(name = "signif", values = 8, breaks = "*", labels = "FDR < 0.05") +
         statsDetailTheme(fontsize = fontsize) +
         theme(strip.background = element_blank(),
               legend.position = c(0.5,legendY),
               legend.justification = c(0.5,0),
               legend.direction = "vertical",
               legend.key.width = unit(0.15, "inches"),
               legend.text = element_text(size = fontsize),
               legend.margin = margin(t = 0, b = 0, unit = "lines"),
               legend.title = element_blank(),
               strip.text = element_text(color = "transparent"),
               panel.spacing = unit(0, "lines"),
               plot.margin = margin(l = -0.2, r = 1, unit = "cm"), 
               axis.text.x.top = element_text(hjust = 0, vjust = 0.5, angle = 90),
               axis.text.x.bottom = element_text(hjust = 1, vjust = 0.5, angle = 90)) 

    if(any(dat$labelY == 0)){
        p <- p +
             geom_rect(data = dat %>% filter(labelY == 0),
                       aes(xmin = -Inf, xmax = Inf, ymin = ystart, ymax = yend),
                       fill = spacerColor, color = spacerColor) +
             geom_hline(yintercept = Inf, size = 0.5, color = "darkgray") +
             geom_hline(yintercept = -Inf, size = 0.5, color = "darkgray") +
             geom_vline(xintercept = Inf, size = 0.5, color = "darkgray")
    }
    if(!is.null(facetY)){
        p <- p +
             facet_grid(as.formula(paste(paste0(facetY,collapse="+"), " ~ .")),
                        scales = "free_y",
                        space="free_y")
    }
 
    if(!keepLegend){ p <- p + theme(legend.position = "none") }
    if(hideXaxisTop){ 
        p <- p + 
             theme(axis.text.x.top = element_blank(), 
                   axis.ticks.x.top = element_blank(),
                   plot.margin = margin(l = -0.2, r = 1, t = 0.5, b = 0.5, unit = "cm")) 
    }
    if(hideXaxisBottom){ p <- p + theme(axis.text.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank()) }
    if(hideXtitle){ p <- p + theme(axis.title.x = element_blank(), axis.title.x.top = element_blank()) }


    gt <- ggplot_gtable(ggplot_build(p))
    gt$layout$clip[grep("panel",gt$layout$name)] <- "off"

    gt

}


#' Plot summary bar charts for each cell state condition in a stats detail plot
#'
#' Create a figure panel in which each condition row contains its own bar chart, with one
#' bar for each sample group in the comparison(s). The values reflected in the bar chart
#' are arbitrary but generally are a summary value like mean or median.
#'
#' @param dat               tibble containing all data to be shown on plot including 
#'                          columns for summary value, y axis, x axis (sample groups)
#'                          and background color(s)
#' @param clrs              named list of R or hex color values, named by a values
#'                          on the x axis
#' @param fillCol           column containing names of clrs list
#' @param xVar              usually the same column as fillCol, the one containing
#'                          values to be shown on x axis
#' @param yVar              column of y values
#' @param title             panel title
#' @param facetY            vector of data column name(s) to be used for row facetting; 
#'                          recommended: Cell_type and Tag
#' @param fontsize          numeric; font size
#' @param spacerColor       character string; R color name or hex color value, including '#';
#'                          fill color of rectangle drawn as a spacer between outer facets
#' @param stripeBG          logical indicating whether to add a horizontal striped background,
#'                          with two alternating colors 
#' @param bgColor           main background color (default = white)
#' @param stripeColor       when stripeBG is TRUE, this color will be one of the two alternating
#' @param keepLegend        logical indicating whether to display figure legend(s)
#' @param hideXaxisTop      remove x axis ticks and labels from top of panel
#' @param hideXaxisBottom   remove x axis ticks and labels from bottom of panel
#' @param hideXtitle        remove x title from panel
#'
#' @param a ggplot table object
plotSummaryBars <- function(dat, clrs, fillCol = "Group", xVar = "Group", yVar = "y", title = "Median Fraction",
                           facetY = NULL, stripeBG = TRUE, bgColor = "white", stripeColor = "#f0f0f0", 
                           fontsize = 12, spacerColor = "#d8d8d8", keepLegend = TRUE,
                           hideXaxisTop = FALSE, hideXaxisBottom = FALSE, hideXtitle = FALSE){

    if(!stripeBG){ stripeColor <- bgColor }

    resTbl <- dat %>%
              formatForStripedBackground(colors = c(bgColor, stripeColor), yVar = yVar) %>%
              addFacetSpacerRows(facetY, xCol = NULL, yCol = yVar, bgColor = spacerColor)
    if(!is.null(facetY)){
        resTbl <- resTbl %>%
                  filter_at(vars(facetY[1]), all_vars(!is.na(.)))
    }

    bgClrs <- c(bgColor, stripeColor, spacerColor)
    names(bgClrs) <- bgClrs

    allClrs <- c(clrs, bgClrs)

    nr <- length(unique(resTbl$`Cell State ID`))
    legendY <- 1 + (0.03 * (100/nr))

    xMax <- max(as.numeric(resTbl[[xVar]]), na.rm=T) + 0.75

    p <- ggplot(data = resTbl, aes_string(x = xVar, y = yVar)) +
         geom_rect(aes(xmin = 0.25, xmax = xMax, ymin = ystart, ymax=yend,
                       fill = col, color = col), size = 0.25, show.legend = F) +
         geom_rect(aes(xmin = as.numeric(x) - 0.45, xmax = as.numeric(x) + 0.45, 
                       ymin = bottom, ymax = top, fill = !!as.name(fillCol))) +
         geom_hline(yintercept = 0.5, size = 0.5, color = "darkgray") +
         geom_hline(yintercept = Inf, size = 0.5, color = "darkgray") +
         geom_vline(xintercept = 0.25, size = 0.5, color = "darkgray") + 
         geom_vline(xintercept = xMax, size = 0.5, color = "darkgray") + 
         xlab(gsub(" ", "\n", title)) +
         scale_color_manual(values = allClrs, guide = FALSE) +
         scale_fill_manual(name = "", values = allClrs, breaks = names(clrs)) +
         scale_y_discrete(expand = c(0,0)) + 
         scale_x_continuous(expand = c(0, 0), limits = c(0.25, xMax), sec.axis = dup_axis()) +
         theme_minimal() +
         statsDetailTheme(fontsize = fontsize) +
         theme(legend.direction = "vertical", 
               legend.position = c(0,legendY),
               legend.justification = c(0,0),
               legend.text = element_text(size = fontsize, hjust = 1),
               legend.title = element_text(size = fontsize, hjust = 1),
               legend.key = element_rect(color = "transparent", fill = "white"),
               legend.key.size = unit(0.5, "cm"),
               legend.key.width = unit(0.5,"cm"),
               legend.background = element_rect(color = "transparent", fill = "transparent"),
               panel.background = element_blank(),
               panel.spacing = unit(0,"lines"),
               axis.text.y.left = element_blank(),
               axis.title.y = element_blank(),
               axis.text.x = element_blank(),
               axis.text.x.top = element_blank(),
               axis.title.x.bottom = element_blank(),
               axis.title.x.top = element_text(size = fontsize, hjust = 0.5),
               plot.margin = margin(l = -0.2, r = -0.83, t=3, b=0.5, unit = "cm"))

    
    if(!is.null(facetY)){
        p <- p +
             facet_grid(as.formula(paste(paste0(facetY,collapse="+"), " ~ .")), space = "free_y", scales = "free_y") +
             theme(strip.background = element_blank(),
                   strip.text = element_text(color = "transparent")) 
    }

    ## shade space between facets
    if(any(!is.na(resTbl$labelY) && resTbl$labelY == 0)){
        p <- p +
             geom_rect(data = resTbl %>% filter(labelY == 0),
                       aes(xmin = 0.25, xmax = xMax, ymin = ystart, ymax = yend),
                       fill = spacerColor)
    }

    if(!keepLegend){ p <- p + theme(legend.position = "none") }
    if(hideXaxisTop){ 
        p <- p + 
             theme(axis.text.x.top = element_blank(), 
                   axis.ticks.x.top = element_blank(),
                   plot.margin = margin(l = -0.2, r = -0.83, t=0.5, b=0.5, unit = "cm")) 
    }
    if(hideXaxisBottom){ p <- p + theme(axis.text.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank()) }
    if(hideXtitle){ p <- p + theme(axis.title.x = element_blank(), axis.title.x.top = element_blank()) }

    gt <- ggplot_gtable(ggplot_build(p))
    gt$layout$clip[grep("panel",gt$layout$name)] <- "off"

    gt

}


#' Check to ensure conditions are properly placed across all figure panels
#' 
#' Ensure that every condition/cell state is assigned the same y value and 
#' facets in every plot data table. 
#'
#' @param labs    tibble of formatted condition labels including columns 'Column',
#'                'labelY', 'y' and 'Value'
#' @param stats   tibble of formatted statistics including columns 'Cell State ID',
#'                'labelY' and 'y'
#' @param indiv   tibble of multiple comparison results, to be used for heatmap
#'                panel in stats detail figure including columns 'Cell State ID',
#'                'labelY' and 'y'
#' @param meds    tibble of formatted summary (median) values to be plotted
#'                in summary bar charts including columns 'Cell State ID', 'labelY' 
#'                and 'y'
#' @param points  tibble of formatted calculation unit data (e.g., FOV) to be
#'                plotted in panel 2 or 3, depending on whether summary bar chart
#'                is to be included; columns must include 'Cell State ID', 'labelY'
#'                and 'y'
#' 
#' @return logical; TRUE when all conditions are properly aligned
conditionsAligned <- function(labs, stats, indiv, meds = NULL, points = NULL){
    ## manual check to make sure all rows are in the correct order
    labLayout <- labs %>%
                 filter(Column == "Cell State ID") %>%
                 select(labelY, y, Value) %>%
                 mutate(`Cell State ID` = as.numeric(as.character(Value)), y = as.numeric(as.character(y))) %>%
                 select(`Cell State ID`, labelY, y)

    statsLayout <- stats %>%
                   select(`Cell State ID`, labelY, y) %>%
                   mutate(y = as.numeric(as.character(y))) %>%
                   arrange(labelY) %>%
                   unique()

    if(!is.null(meds)){
        medsLayout <- meds %>%
                      select(`Cell State ID`, labelY, y) %>%
                      ungroup() %>%
                      mutate(`Cell State ID` = as.numeric(as.character(`Cell State ID`)),
                             y = as.numeric(as.character(y))) %>%
                      arrange(labelY) %>%
                      unique()
    }
    if(!is.null(points)){
        pointsLayout <- points %>%
                        ungroup() %>%
                        select(`Cell State ID`, labelY, y) %>%
                        mutate(`Cell State ID` = as.numeric(as.character(`Cell State ID`)),
                             y = as.numeric(as.character(y))) %>%
                      arrange(labelY) %>%
                      unique()
    }

    indivLayout <- indiv %>%
                   select(`Cell State ID`, labelY, y) %>%
                   mutate(y = as.numeric(as.character(y))) %>%
                   arrange(labelY) %>%
                   unique()

    good <- all(c(identical(labLayout, statsLayout),
                 (!is.null(meds) && identical(labLayout, medsLayout) | is.null(meds)),
                 (!is.null(points) && identical(labLayout, pointsLayout) | is.null(points)),
                 identical(labLayout, indivLayout)))

    if(!good){
        log_error("ROW ORDER DOES NOT MATCH BETWEEN PANELS")
        return(FALSE)
    }
    TRUE
}


plotSeparateCohortVsIndividual <- function(dat, indivDat, cfg, statsFiles, sheet, comps, calcType, cellTypes, conds,
                                           idMap, effectCol = "Odds Ratio", effectAbb = "OR", calcCol = "Fraction",
                                           facetY = NULL, facetOrder = NULL, popOrder = NULL, orderBy = NULL,
                                           xVar = "Sample_ID", xOrder = NULL){

    idOrder <- NULL
    figs <- list()
    pdfHeight <- 0
    rel_heights <- c()
    widths <- list(fractions = c(2, 0.4, 0.8, 1.2),
                   densities = c(1.6, 0.4, 0.8, 1.2))

    comp1 <- comps$comp1$name
    comp2 <- comps$comp2$name

    ## create two separate plots for conditions that are up/up or down/down and opposite directions
    for(type in c("same", "opposite")){
        log_debug(paste("Comparison Type: comparison 1 and comparison 2 in", toupper(type), "direction"))

        excl <- dat %>%
                filter(!!as.name(paste(comp1, "Overall", effectAbb)) == 1 |
                       !!as.name(paste(comp2, "Overall", effectAbb)) == 1)
        log_debug(paste("        [FILTER] removed conditions where at least one effect == 1:", numConds(excl)))

        if(type == "same"){
            c2p <- dat %>%
                   filter((!!as.name(paste(comp1, "Overall", effectAbb)) > 1 &
                           !!as.name(paste(comp2, "Overall", effectAbb)) > 1) |
                          (!!as.name(paste(comp1, "Overall", effectAbb)) < 1 &
                             !!as.name(paste(comp2, "Overall", effectAbb)) < 1))  %>%
                   pull(`Cell State ID`)
            log_debug(paste("        [PLOT] conditions with", effectAbb, ">1 in both",
                                     comp1, "and", comp2, ":", length(c2p)))

        } else {
            c2p <- dat %>%
                   filter((!!as.name(paste(comp1, "Overall", effectAbb)) > 1 &
                            !!as.name(paste(comp2, "Overall", effectAbb)) < 1) |
                          (!!as.name(paste(comp1, "Overall", effectAbb)) < 1 &
                           !!as.name(paste(comp2, "Overall", effectAbb)) > 1)) %>%
                   pull(`Cell State ID`)
            log_debug(paste("        [PLOT] conditions with", effectAbb,
                                     "<1 in one comparison and >1 in the other:", length(c2p)))
        }

        ### CONDITION LABELLING
        ## for labeling, just use any single stats file (they are all the same)
        idOrder <- NULL
        statsFile <- statsFiles[grepl(paste0(comps$comp1$question_number, ".xlsx"), statsFiles)]

        labs   <- formatConditionsForPlotting(conds, c2p, calcType, cellTypes, facetY = facetY, idOrder = idOrder,
                                              facetOrder = facetOrder, popOrder = popOrder, orderBy = orderBy,
                                              statsFile = statsFile, sheet = sheet)
        idOrder <- labs %>%
                   filter(Column == "Cell State ID") %>%
                   select(Value, labelY) %>%
                   pull(Value) %>%
                   as.numeric %>% rev %>% unique

        ## get effect sizes & CIs for both overall comparsions
        allStats <- tibble()
        for(cmp in names(comps)){
            qu <- comps[[cmp]]$question_number
            statsFile <- statsFiles[grepl(paste0(qu,".xlsx"), statsFiles)]
            stats <- formatStatsForPlottingEffects(statsFile, sheet, calcType, calcCol, conds, c2p,
                                                   cellTypes = cellTypes, idOrder = idOrder,
                                                   facetY = facetY, facetOrder = facetOrder,
                                                   popOrder = popOrder, orderBy = orderBy) %>%
                     mutate(CompGroup = comps[[cmp]]$name)
            allStats <- allStats %>% bind_rows(stats)
        }
        insig <- which((allStats[[effectCol]] >= 1 & allStats$CI.low <= 1) |
                       (allStats[[effectCol]] <= 1 & allStats$CI.high >= 1))
        allStats$CompGroup[insig] <- "not significant"
        layout <- allStats %>% select_at(c("Cell State ID", facetY, "labelY", "y")) %>% unique()

        grpIdxs <- grep("_label", names(unlist(comps)))
        grps    <- unlist(comps)[grpIdxs] %>% unique
        allMeds <- formatMedians(dat, c2p, grps, calcCol, layout)

        ## format heatmap of individual comparison "same" or "opposite" effect direction when 
        ## compared to the overall comparison
        indivDatT <- indivDat %>%
                     filter(`Cell State ID` %in% c2p) %>%
                     formatIndivDatForPlotting(idMap, layout, xVar = xVar, xOrder = xOrder)

        if(!conditionsAligned(labs, allStats, indivDatT, meds = allMeds)){
            stop()
        }

        ## make figure
        clrs <- c('#32cd32', '#f5d000', 'darkgray')
        names(clrs) <- c(comps$comp1$name, comps$comp2$name, 'not significant')

        medClrs <- c("#32cd32", "#525253", "#f5d000")
        names(medClrs) <- c(comps$comp1$group_1_label,
                            comps$comp1$group_2_label,
                            comps$comp2$group_1_label)
       log_info("Generating figure...")

        strpBG  <- TRUE
        bg <- "white"
        stripe <- "#f0f0f0"
        spacer <- "#e0e0e0"
        fs <- 12

        plotParams <- list("same" = list(hideXtitle = FALSE, hideXaxisTop = FALSE,
                                         hideXaxisBottom = TRUE, keepLegend = TRUE),
                           "opposite" = list(hideXtitle = TRUE, hideXaxisTop = TRUE,
                                             hideXaxisBottom = FALSE, keepLegend = FALSE))
        pp <- plotParams[[type]]

        options(warn = -1)
        ## panel 1
        condTbl   <- conditionTable(labs, cellTypes, facetY = facetY, fontsize = fs,
                                    stripeBG = strpBG, bgColor = bg, stripeColor = stripe,
                                    spacerColor = spacer, calcType = calcType, hideXaxisTop = pp$hideXaxisTop,
                                    hideXaxisBottom = pp$hideXaxisBottom)
        ## panel 2
        title = ifelse(calcType == "fractions", "Median Fraction Scaled", "Median Density Scaled")
        medians <- plotMedianBars(allMeds, medClrs, fillCol = "Group", xVar = "x", yVar = "y",
                                  title = title, facetY = facetY, fontsize = fs, stripeBG = strpBG,
                                  bgColor = bg, stripeColor = stripe, spacerColor = spacer,
                                  keepLegend = pp$keepLegend, hideXtitle = pp$hideXtitle,
                                  hideXaxisTop = pp$hideXaxisTop, hideXaxisBottom = pp$hideXaxisBottom)
        ## panel 3
        allStats$CompGroup <- factor(allStats$CompGroup, levels = names(clrs))
        statPlot  <- plotEffectSize(allStats, clrs = clrs, colorBy = "CompGroup", fontsize = fs,
                                    effectCol = effectCol, facetY = facetY, yVar = "y",
                                    stripeBG = strpBG, bgColor = bg,
                                    stripeColor = stripe, spacerColor = spacer, hideXtitle = pp$hideXtitle,
                                    hideXaxisTop = pp$hideXaxisTop, hideXaxisBottom = pp$hideXaxisBottom,
                                    keepLegend = pp$keepLegend)
        ## panel 4
        ptPlot    <- plotComparisonHeatmap(indivDatT, "Sample_ID", yVar = "y", facetY = facetY,
                                           spacerColor = spacer, fill = "status", calcType = calcType,
                                           hideXtitle = pp$hideXtitle, hideXaxisTop = pp$hideXaxisTop,
                                           hideXaxisBottom = pp$hideXaxisBottom, keepLegend = pp$keepLegend)

        ## match all heights for proper alignment
        condAndStatList <- matchPanelHeights(condTbl, statPlot)
        condAndMedsList <- matchPanelHeights(condTbl, medians)
        statAndPtList   <- matchPanelHeights(condAndStatList[[1]], ptPlot)

        plotList <- list(condAndStatList[[1]], condAndMedsList[[2]], condAndStatList[[2]], statAndPtList[[2]])
        rel_widths <- widths[[calcType]]

        pg <- plot_grid(plotlist = plotList, align = 'hv', axis = 'bt', nrow = 1, rel_widths = rel_widths)
        figs[[type]] <- pg
        exHeight <- ifelse(type == "same", 2.5, 0.9)
        pdfHeight <- pdfHeight + (0.11 * nrow(labs)) + exHeight
        rel_heights <- c(rel_heights, (0.11 * nrow(labs)) + exHeight)
        log_debug(paste0("Relative heights: ", paste(rel_heights, collapse = ", ")))
        log_debug(paste0("PDF height: ", pdfHeight))
    }
    ### Plotting all passing data
    pdfFile <- file.path(outDir, paste(comps$comp1$question_number, comps$comp2$question_number,
                                       ct, "combined_detail_figure.pdf", sep="_"))
    pdfWidth <- 11.5

    pdf(pdfFile, height = pdfHeight, width = pdfWidth, onefile=FALSE)
    pg2 <- plot_grid(plotlist = figs, align = 'hv', axis = 'lr', ncol = 1, rel_heights = rel_heights)
    grid.draw(pg2)
    dev.off()
    #options(warn = 0)
}
 

