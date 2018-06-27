################################################################################
## Helper functions for mindusa_rct.Rmd
##  Stored separately to help with debugging, clarity
################################################################################

## -- General helper functions (summary stats, simple calculations) ------------
factor_tf <-
  function(x){ factor(as.numeric(x), levels = 1:0, labels = c("Yes", "No")) }
format_comma <- partial(format, big.mark = ",")
sum_na <- partial(sum, na.rm = TRUE)
q25 <- partial(quantile, probs = 0.25, na.rm = TRUE)
q50 <- partial(quantile, probs = 0.50, na.rm = TRUE)
q75 <- partial(quantile, probs = 0.75, na.rm = TRUE)
mean_na <- partial(mean, na.rm = TRUE)
sd_na <- partial(sd, na.rm = TRUE)
rndformat <- function(x, digits = 2){ format(round(x, digits), nsmall = digits) }
get_npct <- function(num, denom){
  sprintf("%s (%s%%)", num, round((num / denom) * 100))
}
describe_cont <- function(
  v_q50, v_q25, v_q75, v_mean, v_sd, dig_iqr = 0, dig_msd = 2
){
  sprintf(
    "%s [%s, %s]\\\n%s +/- %s",
    rndformat(v_q50, digits = dig_iqr),
    rndformat(v_q25, digits = dig_iqr),
    rndformat(v_q75, digits = dig_iqr),
    rndformat(v_mean, digits = dig_msd),
    rndformat(v_sd, digits = dig_msd)
  )
}

## -- Printing helper functions ------------------------------------------------
## Wrapper with my preferred options for printing summaryM objects
my_html <- function(sMobj, caption){
  html(
    sMobj,
    exclude1 = FALSE, long = TRUE, digits = 2, what = "%", npct = "both",
    prmsd = TRUE, brmsd = TRUE, middle.bold = TRUE,
    ## These options don't seem to be working.
    msdsize = mu$smaller2, outer.size = mu$smaller2, rowsep = TRUE,
    caption = caption
  )
}

## kable_styling wrapper to ensure all tables are consistently styled
mykablestyle <- function(obj, stripes = FALSE, ...){
  boptions <- c("hover", "responsive", "condensed", "bordered")
  if(stripes){ boptions <- c(boptions, "striped") }
  
  kable_styling(
    obj,
    bootstrap_options = boptions,
    full_width = FALSE,
    ...
  ) %>%
    row_spec(0, bold = TRUE, background = palette_colors[["lgray"]])
  
}

## Format p-values per NEJM style
formatp_nejm <- function(p){
  ifelse(p < 0.0001, '<0.0001',
  ifelse(p < 0.001, '<0.001',
  ifelse(p < 0.01, rndformat(p, digits = 3),
         rndformat(p, digits = 2))))
}

## -- Continuous outcomes: Modeling/results helper functions -------------------
## Create string with results of Kruskal-Wallis test
kw_results <- function(kw_obj){
  
  if(!inherits(kw_obj, "htest")){
    stop("kw_obj must be the result of kruskal.test()", call. = FALSE)
  }
  
  ## Couldn't get multiple lines to work
  ##  (doesn't take \n; tried list of bquote objects)
  bquote(
    "Kruskal-Wallis test:" ~
      X^2 * "," ~ .(rndformat(kw_obj$statistic, 2)) * "; df," ~
      .(kw_obj$parameter) * "; P," ~ .(formatp(kw_obj$p.value))
  )
}

## Plot distributions of mental status variables by treatment, including
## results of Kruskal-Wallis test
mental_plot <- function(xvar, xtitle, kw_obj){
  ## No patient should have 14 DCFDs (must have delirium to be randomized),
  ## but could have 14 delirium/coma days
  if(xvar == "dcfd_int"){
    xmax <- 13
  } else{
    xmax <- 14
  }
  
  ggplot(data = model_df, aes_string(x = xvar)) +
    facet_wrap(~ trt, nrow = 1) +
    geom_histogram(aes(fill = trt), binwidth = 1) +
    scale_fill_manual(values = mindusa_trtcols_long) +
    scale_x_continuous(name = xtitle, breaks = seq(0, xmax, 2)) +
    scale_y_continuous(name = "Patient Count") +
    theme(
      legend.position = "none",
      plot.caption = element_text(size = 10)
    ) +
    labs(
      title = glue("{xtitle} by Treatment"),
      subtitle = kw_results(kw_obj)
    )
  
}

## Function to plot adjusted medians or odds ratios for mental status variables
## df = object created by quantile_orm_df()
mental_medians_plot <- function(df, mod, text_results = TRUE){
  if(!inherits(mod, "rms")){
    stop("mod must be an rms model object", call. = FALSE)
  }
  
  ## Reorder treatment variable
  df$trt <- fct_relevel(df$trt, c("Placebo", "Haloperidol", "Ziprasidone"))
  
  ## Max should be 13 for DCFDs, 14 for delirium/coma
  ## Font needs to be smaller for delirium duration to fit entire title
  if(mod$sformula[[2]] == "dcfd_int"){
    xmax <- 13
    base_pct <- 1.2
    outcome <- "DCFDs"
  } else{
    xmax <- 14
    if(mod$sformula[[2]] == "del_int_all"){
      base_pct <- 1.1
      outcome <- "Delirium Duration"
    } else if(mod$sformula[[2]] == "hyperdel_int_all"){
      base_pct <- 1.0
      outcome <- "Duration of Hyperactive Delirium"
    } else if(mod$sformula[[2]] == "hypodel_int_all"){
      base_pct <- 1.0
      outcome <- "Duration of Hypoactive Delirium"
    } else{
      base_pct <- 1.2
      outcome <- "Coma Duration"
    }
  }
  
  ## Add exact results as part of Y axis labels, if desired
  if(text_results){
    ## Create named vector of new factor labels
    med_labels <- df %>%
      mutate(
        med_label = glue(
          "\n\n{trt}\n",
          "{rndformat(quantile)} ({rndformat(lb)}, {rndformat(ub)})"
        )
      ) %>%
      pull(med_label) %>%
      set_names(df$trt)
    
    df$trt <- forcats::fct_relabel(df$trt, ~ med_labels[.])
  }
  
  p <- ggplot(data = df, aes(y = quantile, x = trt)) +
    geom_pointrange(
      aes(ymin = lb, ymax = ub),
      shape = 16, size = 0.8, colour = as.character(palette_colors["dred"])
    ) +
    scale_y_continuous(limits = c(0, xmax), breaks = seq(0, xmax, 2)) +
    labs(
      title = glue("Adjusted Median {outcome} by Treatment"),
      subtitle = glue("P: {formatp(anova(mod)['trt', 'P'])}"),
      y = glue("Adjusted Median (95% Confidence Interval)"),
      caption =
        "\nAdjusted analysis using proportional odds logistic regression."
    ) +
    theme(
      legend.position = "none",
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = basetext_size * base_pct),
      plot.caption = element_text(size = basetext_size * 0.7)
    ) +
    coord_flip()
  
  return(p)
}

## Calculate medians, bounds for model with interaction term (effect modifier)
calc_medians_em <- function(mod, df, em_vals, em_var){
  quantile_orm_df(
    mod = mod,
    new.data = set_names(df, str_subset(names(coef(mod)), "^[^y>=]")),
    trt_levels = rep(levels(model_df$trt), each = length(em_vals))
  ) %>%
    bind_cols(dplyr::select(df, one_of(em_var))) %>%
    mutate(
      outcome = as.character(mod$sformula[[2]]),
      p_int =
        anova(mod)[paste(em_var, "* trt  (Factor+Higher Order Factors)"), "P"]
    )
}

## Plot medians, bounds for mental status outcomes with interaction terms
## median_df assumed to have cols: quantile, lb, ub, trt, em_text, outcome_text
mental_medians_plot_em <- function(median_df, em_string){
  ggplot(data = median_df, aes(y = quantile, x = trt_short)) +
    facet_grid(em_text ~ outcome_text) +
    geom_pointrange(
      aes(ymin = lb, ymax = ub),
      colour = as.character(palette_colors["dred"]), size = 0.5
    ) +
    scale_y_continuous(
      breaks = seq(0, 14, 2), name = "Adjusted Median (95% Confidence Interval)"
    ) +
    coord_flip() +
    labs(
      caption = glue(
        "Adjusted analysis using proportional odds logistic regression.\n",
        "P-values for overall {em_string} * treatment interaction."
      )
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text = element_text(size = basetext_size * 0.7),
      panel.spacing.y = unit(0.5, "cm"),
      strip.text.x = element_text(vjust = 0),
      plot.caption = element_text(face = "italic")
    )
}

## -- Helper functions for any rms model object --------------------------------
## Wrapper for rms_calc_comparisons; adds outcome variable as a column
## (Tested on lrm, cph fits)
rms_comparisons_addoutcome <- function(rmsObj, ...){
  outcome <- gsub("()", "", rmsObj$sformula[2], fixed = TRUE)
  rms_calc_comparisons(rmsObj, ...) %>%
    mutate(outcome = outcome)
}

## Print rms_model_results using kableExtra
rms_results_kable <- function(rmsObj, ...){
  ## Label for effect size column
  output_type <- case_when(
    inherits(rmsObj, "lrm") ~ "Odds Ratio",
    inherits(rmsObj, "orm") ~ "Odds Ratio",
    inherits(rmsObj, "cph") ~ "Hazard Ratio",
    inherits(rmsObj, "ols") ~ "Difference",
    TRUE                    ~ "Effect"
  )
  
  rms_model_results(rmsObj, ...) %>%
    ## Remove rows that give number of observations at every outcome level;
    ##  we already describe the distribution of ordinal outcomes
    filter(is.na(as.numeric(label))) %>%
    kable(
      format = "html",
      align = c("l", rep("r", 6)),
      col.names = c(
        "Variable", "Reference", "Comparison",
        sprintf("%s (95%% CI)", output_type),
        "X^2^", "df", "P"
      )
    ) %>%
    mykablestyle()  
}

## Calculate ratios for both treatments vs placebo;
##   result = tidy data frame we can plot! Hooray!
trt_ratios <- function(mod){
  map_df(
    list(mod),
    rms_comparisons_addoutcome,
    df = model_df, getRatios = TRUE, vname = "trt", refVal = "Placebo"
  ) %>%
    mutate(comp.c = factor(comp.c, levels = trt_levels$trt_actual)) %>%
    ## Keep only one row for reference group
    distinct(ref, comp, effect, is.ref, .keep_all = TRUE)
}

## Plot results of trt_ratios
## - Treatment on Y axis
## - Point estimates + CIs for ziprasidone, haloperidol vs placebo on X axis
## - Black square for placebo
plot_trt_ratios <- function(
  ratio_df,     ## data.frame w/ one row per treatment, likely from trt_ratios()
  ## columns include effect, lcl, ucl, comp.c, ref.c
  outcome_text, ## Text describing outcome (eg, "Delirium/Coma-Free Days")
  extra_text = "", ## Additional text to include on the second line of caption
  text_results = TRUE, ## whether to include actual results in Y axis labels
  mod           ## rms model object from which P for trt will be extracted
){
  ## Extract X axis label from model type: orm/lrm = odds, cph = hazard
  if(inherits(mod, "cph")){
    ratio_type <- "Hazard"
    ## Hacky way to check whether cph() fit uses Fine-Gray approach or not:
    ## if it does, there will be weights. This will not extend well if we ever
    ## want to use weights for a different scenario, but works for now.
    if(is.null(mod$weights)){
      model_type <- "Cox proportional hazards"
    } else{
      model_type <- "competing risks"
    }
  } else{
    ratio_type <- "Odds"
    model_type <- "proportional odds logistic"
  }
  
  ## Reorder treatment variable
  ratio_df$comp.c <- fct_relevel(
    ratio_df$comp.c, c("Placebo", "Haloperidol", "Ziprasidone")
  )
  
  ## Add exact results as part of Y axis labels, if desired
  if(text_results){
    ## Create named vector of new factor labels
    comp_labels <- ratio_df %>%
      mutate(
        comp_label = ifelse(
          is.ref, glue("\n\n{comp.c}\n\n"),
          glue("\n\n{comp.c}\n",
               "{rndformat(effect)} ({rndformat(lcl)}, {rndformat(ucl)})")
        )
      ) %>%
      pull(comp_label) %>%
      set_names(ratio_df$comp.c)
    
    ratio_df$comp.c <- forcats::fct_relabel(
      ratio_df$comp.c,
      ~ comp_labels[.]
    )
  }
  
  p <- ggplot(data = ratio_df, aes(y = effect, x = comp.c)) +
    ## Fake row to set up order properly
    geom_point(shape = NA) +
    ## Reference line at 1 (no effect)
    geom_hline(yintercept = 1, linetype = "solid",
               colour = palette_colors["lgray"], size = 2) +
    ## Plot a point for control group
    geom_point(
      shape = 15, colour = "black", size = 4,
      data = ratio_df %>% filter(is.ref)
    ) +
    ## Add ratios, CIs for treatment groups vs control
    geom_pointrange(
      aes(ymin = lcl, ymax = ucl),
      position = position_dodge(width = 0.5),
      shape = 16, size = 1, colour = as.character(palette_colors["dred"]),
      data = ratio_df %>% filter(!is.ref)
    ) +
    labs(
      title = glue("Treatment vs {outcome_text}"),
      subtitle = glue("P: {formatp(anova(mod)['trt', 'P'])}"),
      y = glue("{ratio_type} Ratio (95% Confidence Interval)"),
      caption = glue(
        "\nAdjusted analysis using {model_type} regression.\n{extra_text}"
      )
    ) +
    theme(
      legend.position = "none",
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      # axis.text.y = element_text(vjust = 1),
      plot.caption = element_text(size = basetext_size * 0.7)
    ) +
    coord_flip()

  return(p)  
}

plot_trt_ratios_em <- function(
  ratio_df,     ## data.frame w/ one row per treatment per outcome
                ## columns include effect, lcl, ucl, comp.c, ref.c
  em_string,    ## Text describing effect modifier
  ratio_type = c("Hazard", "Odds"),
  facet_formula = "em_text ~ outcome_text" ## string specifying how to facet
){
  
  if(ratio_type == "Hazard"){
    caption_text <- glue(
      "Adjusted mortality outcomes use standard Cox proportional hazards regression.\n",
      "Other outcomes use Fine-Gray competing risks regression, with competing risk of death\n",
      "(and ICU discharge without the event for MV liberation and ICU readmission).\n\n",
      "P-values for overall {em_string} * treatment interaction."
    )
  } else{
    caption_text <- glue(
      "Adjusted analysis using proportional odds logistic regression.\n",
      "P-values for overall {em_string} * treatment interaction."
    )
  }

  p <- ggplot(data = ratio_df, aes(y = effect, x = comp.c_short)) +
    ## Facet for each outcome, interacting value
    facet_grid(facet_formula) +
    ## Fake row to set up order properly
    geom_point(shape = NA) +
    ## Reference line at 1 (no effect)
    geom_hline(yintercept = 1, linetype = "solid",
               colour = palette_colors["lgray"], size = 1) +
    ## Plot a point for control group
    geom_point(
      shape = 15, colour = "black", size = 2,
      data = ratio_df %>% filter(comp.c == ref.c)
    ) +
    ## Add ratios, CIs for treatment groups vs control
    geom_pointrange(
      aes(ymin = lcl, ymax = ucl),
      position = position_dodge(width = 0.5),
      shape = 16, size = 0.5, colour = as.character(palette_colors["dred"]),
      data = ratio_df %>% filter(comp.c != ref.c)
    ) +
    labs(
      title = glue("Treatment vs {capitalize(em_string)}"),
      y = glue("{ratio_type} Ratio (95% Confidence Interval)"),
      caption = caption_text
    ) +
    theme(
      legend.position = "none",
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = basetext_size * 0.7),
      strip.text.x = element_text(vjust = 0),
      panel.spacing.y = unit(0.5, "cm"),
      plot.caption = element_text(size = basetext_size * 0.7)
    ) +
    coord_flip()
  
  return(p)  
}

## Tabular results of POLR, Cox models: Ratios for each treatment group
table_trt_ratios <- function(
  ratio_df, ## data.frame w/ one row per treatment, likely from trt_ratios()
  ## columns include effect, lcl, ucl, comp.c, ref.c
  ratio_type = c("Hazard", "Odds") ## Type of ratio (hazard or odds)
){
  ## Set ratio_type to default if not specified
  ##  (Primary outcomes include more Cox models than POLR models)
  ratio_type <- match.arg(ratio_type)
  
  ratio_df %>%
    ## Create character strings for prettier tables
    mutate(
      comparison = glue("{comp.c} vs {ref.c}"),
      effect_ci = glue(
        "{rndformat(effect)} ",
        "({rndformat(lcl)}, ",
        "{rndformat(ucl)})"
      )
    ) %>%
    ## Arrange in consistent order
    arrange(desc(ref), comp) %>%
    ## Only want newly created character columns
    dplyr::select(comparison, effect_ci) %>%
    ## Printing options: HTML format, set column names/alignment
    kable(
      format = "html",
      col.names = c("Comparison", glue("{ratio_type} Ratio (95% CI)")),
      align = c("l", "r")
    ) %>%
    mykablestyle(position = "left") %>%
    ## Include "placebo vs placebo" row for clarity, but gray out
    row_spec(3, color = palette_colors["lgray"])
}

## -- Time to event outcomes: Modeling/results helper functions ----------------

## Unadjusted analyses: Kaplan-Meier curves for time to death ------------------
km_plot_death <- function(
  sf,              ## survfit object
  lr,              ## log-rank test object
  time = c(30, 90) ## time frame of interest
){
  
  if(!inherits(sf, "survfit")){
    stop("sf must be a survfit object", call. = FALSE)
  }
  
  if(!inherits(lr, "survdiff")){
    stop("lr must be a survdiff object", call. = FALSE)
  }
  
  ## -- Get chi-square statistic, p-value from log-rank test -------------------
  ## df = number of treatment groups - 1
  df <- length(names(lr$n)) - 1
  chis <- lr$chisq
  pval <- 1 - pchisq(chis, df = df)
  
  ## -- For 30 days, break every 5 days; for 90 days, break every 15 -----------
  time_breaks <- ifelse(time == 30, 5, 15)
  
  ## -- Create base plot using survminer package -------------------------------
  km_death <- survminer::ggsurvplot(
    sf, ## survfit() object
    
    ## CIs
    conf.int = TRUE, conf.int.alpha = 0.20,
    
    ## Show N, cumulative events every time_breaks days; color text by trt
    xlim = c(0, time + 1),
    break.time.by = time_breaks, risk.table = "nrisk_cumevents",
    risk.table.fontsize = 3, risk.table.col = "strata", tables.height = 0.35,
    legend.labs = gsub("trt=", "", names(sf$strata)),
    
    ## Use specified themes, colors
    palette = mindusa_trtcols_long,
    ggtheme = mindusa_theme(), tables.theme = mindusa_theme(),
    font.family = basetext_family
  )
  
  ## -- More finely control plot options ---------------------------------------
  km_death$plot <- km_death$plot +
    ## Change Y axis to % labels, not proportion
    scale_y_continuous(
      name = "Percent Alive",
      breaks = seq(0, 1, 0.25),
      labels = paste0(seq(0, 1, 0.25) * 100, "%")
    ) +
    labs(
      title = glue("Kaplan-Meier Curve, {time}-Day All-Cause Death"),
      ## Add log-rank test results to subtitle, using bquote for superscript/
      ## inclusion of variables
      subtitle = bquote(
        "Log-rank test for difference between treatments:" ~ X^2 ~ ", " ~
          .(rndformat(chis, 1)) * "; df, " ~ .(df) * "; P, " ~ .(formatp(pval))
      )
    ) +
    ## Place legend in bottom lefthand corner, vertically so full trt names fit
    theme(
      ## X axis label will be at the bottom of the table; no need to duplicate
      axis.title.x = element_blank(),
      legend.position = c(0.01, 0.025),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.justification = c(0, 0),
      legend.direction = "vertical"
    )
  
  ## Add labels to and remove unnecessary reference lines from nrisk table
  km_death$table <- km_death$table +
    labs(title = "Number at risk (cumulative number of deaths)") +
    xlab("Days after Randomization") +
    theme(
      plot.title = element_text(size = basetext_size * 0.8),
      axis.title.x = element_text(vjust = 0),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  return(km_death)
  
}

## Graphically check proportional hazards assumption; mod_zph is class cox.zph
plot_ph_assume <- function(mod_zph){
  par(mfrow = c(3, 3))
  plot(mod_zph, col = 'red', lwd = 2)
  abline(v = 0, lty = 2, col = 'gray')
  par(mfrow = c(1, 1))
}

## Extract p-value for global PH assumption; create reproducible statement
text_ph_assume <- function(mod_zph){
  global_phtest <- mod_zph$table["GLOBAL", "p"]
  global_phtext <- ifelse(
    global_phtest < 0.05,
    glue(
      "The global test for proportional hazards indicates the assumption is ",
      "not fully met (p: {formatp(global_phtest)}); please see the figures ",
      "below for further detail."
    ),
    glue(
      "The global test for proportional hazards indicates no major concern ",
      "(p = {formatp(global_phtest)}); please see the figures below for ",
      "further detail."
    )
  )
  
  return(global_phtext)
}

## Create text examples of HR interpretation, given result of trt_ratios()
hr_example <- function(hr_df, outcome){
  ## Death is modeled using regular Cox regression, thus provides a hazard
  ## Other outcomes use competing risks/Fine-Gray approach, thus provide a
  ##  subdistribution HR for the relative incidence
  if(outcome == "death"){
    hr_type <- "hazard"
  } else{
    hr_type <- "incidence"
  }
  
  hr_df %>%
    dplyr::select(-outcome) %>%
    mutate(
      direction = case_when(
        is.ref     ~ "",
        effect > 1 ~ "higher",
        TRUE       ~ "lower"
      ),
      pct_change = case_when(
        is.ref     ~ 0,
        effect > 1 ~ (effect - 1) * 100,
        TRUE       ~ (1 - effect) * 100
      )
    ) %>%
    glue_data(
      "For example, an HR of {rndformat(effect)} for {comp.c} vs placebo ",
      "indicates that, on average, patients in the {comp.c} group have a ",
      "{round(pct_change, 0)}% {direction} relative {hr_type} of {outcome} ",
      "vs patients in the {ref.c} group, given that they have not yet ",
      "experienced the outcome of interest."
    ) %>%
    set_names(hr_df %>% pull(comp.c))
}

## -- Competing risks regression -----------------------------------------------

## String describing test results from cmprsk::cuminc()
cuminc_test <- function(cuminc_obj){
  df <- cuminc_obj$Tests %>%
    ## as_tibble makes it easy to keep rownames as new variable
    as_tibble(rownames = "subdist")
  
  ## Could not figure out how to automate/map this into a bquote() object,
  ##  and I'm sad about it. I know my competing risks will always have two
  ##  events so I'll live.
  ## These^ are famous last words, since we added the additional competing risk
  ## of hospital discharge. Hilarious. **At this point** we always have 2 or 3.
  if(nrow(df) == 2){
    bquote(
      .(df$subdist[1]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[1], 2)) * "; df," ~ .(df$df[1]) * "; P," ~ .(formatp(df$pv[1])) * "." ~ .(df$subdist[2]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[2], 2)) * "; df," ~ .(df$df[2]) * "; P," ~ .(formatp(df$pv[2])) * "."
    )
  } else if(nrow(df) == 3){
    bquote(
      .(df$subdist[1]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[1], 2)) * "; df," ~ .(df$df[1]) * "; P," ~ .(formatp(df$pv[1])) * "." ~ .(df$subdist[2]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[2], 2)) * "; df," ~ .(df$df[2]) * "; P," ~ .(formatp(df$pv[2])) * "." ~ .(df$subdist[3]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[3], 2)) * "; df," ~ .(df$df[3]) * "; P," ~ .(formatp(df$pv[3])) * "."
    )
  }
  
}

## -- Functions to extract and plot the N at risk, cumulative events given -----
## -- a summary(survfit(...)) object -------------------------------------------
cr_risktable <- function(
  sf_sum,    ## summary(survfit(...))
  main_event ## character string; one of the event types in sf_sum
){
  ## Get character strings of event types
  event_types <- colnames(sf_sum$p0)[1:(ncol(sf_sum$n.event) - 1)]

  ## For our purposes, we assume competing risks are "Death" and possibly
  ## "Discharge". Issue a warning if this is not the case.
  ## There is probably a better way to do this. Refactoring opportunity!
  
  if(!(length(unique(c(event_types, "Discharge", "Death"))) == 3)){
    warning(
      "Competing risks are assumed to be `Death` and possibly `Discharge`, but events in `sf_sum` do not match this assumption",
      call. = FALSE
    )
  }
     
  df <- data.frame(
    trt = str_replace(sf_sum$strata, "trt=", ""),
    time = sf_sum$time,
    n_risk = sf_sum$n.risk[, ncol(sf_sum$n.risk)],
    n_censored = sf_sum$n.censor
  ) %>%
    bind_cols(
      as.data.frame(sf_sum$n.event[, 1:(ncol(sf_sum$n.event) - 1)]) %>%
        set_names(event_types)
    ) %>%
    group_by(trt) %>%
    mutate_at(vars(event_types, "n_censored"), cumsum) %>%
    ungroup() %>%
    set_names(
      str_replace(names(.), main_event, "primary")
    )
  
  ## -- Create string for N at risk, etc ---------------------------------------
  ## If time = 0, only include N at risk
  ## If death is only competing risk, include (events; death) on line 2
  ## If death&discharge are competing risks: N -> (events) -> (death; discharge)
  ## (tried glue_data, couldn't get it to work)
  if("Discharge" %in% names(df)){
    df$risk_string <- with(df, {
      ifelse(time == 0, as.character(n_risk),
             sprintf("%s\n(%s)\n(%s; %s)", n_risk, primary, Death, Discharge))
    })
  } else{
    df$risk_string <- with(df, {
      ifelse(time == 0, as.character(n_risk),
             sprintf("%s\n(%s; %s)", n_risk, primary, Death))
    })
  }

  return(df)
  
}

cr_risktable_plot <- function(sf_sum, main_event, event_string){
  df <- cr_risktable(sf_sum, main_event = main_event)
  
  # ## Put plot in correct order - BUT can't quickly figure this out for
  # ##  ggcompetingrisks, so leave alphabetical for now :(
  # df$trt <- factor(df$trt, levels = c("Placebo", "Haloperidol", "Ziprasidone"))
  
  ## Set plot title
  plot_title <- ifelse(
    "Discharge" %in% names(df),
    glue(
    "Number at risk (cumulative number of {event_string}) (deaths; discharges)"
    ),
    glue(
      "Number at risk (cumulative number of {event_string}; deaths)"
    )
  )
  
  ggplot(data = df, aes(x = time, y = 1, colour = trt)) +
    facet_wrap(~ trt, nrow = 1) +
    geom_text(
      aes(label = risk_string),
      size = basetext_size * 0.15, family = basetext_family,
      color = palette_colors[["dgray"]]
    ) +
    # scale_colour_manual(values = mindusa_trtcols_long) +
    scale_x_continuous(
      name = "Days after Randomization",
      breaks = unique(sf_sum$time),
      labels = unique(sf_sum$time)
    ) +
    scale_y_continuous(breaks = 1, labels = "100%") +
    labs(
      title = plot_title,
      ylab = " " ## faking it for spacing purposes
    ) +
    theme(plot.title = element_text(size = basetext_size * 0.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          ## Set Y axis text to background color; need placeholder text for
          ## spacing, so that facets line up
          axis.title.y = element_text(colour = "white"),
          axis.text.y = element_text(colour = "white"),
          axis.ticks.y = element_blank(),
          legend.position = "none")
  
}

## -- Function to create plot for cumulative incidence functions ---------------
## Relies on survminer::ggcompetingrisks
## test_string should be result of cuminc_test()
cif_plot <- function(
  cuminc_obj,    ## class cuminc
  legend_string, ## event of interest, for plot legend
  event_string,  ## event of interest, for plot title
  test_string,   ## to include in subtitle; often result of cuminc_test()
  caption_string = NULL ## optional plot caption
){
  ## What are our event types?
  event_types <- rownames(cuminc_obj$Tests)
  
  ## Which colors to use? Red always main event. Grays = competing risks.
  risk_colors <- c(
    as.character(palette_colors["lgray"]),
    as.character(palette_colors["dgray"]),
    as.character(palette_colors["dred"])
  ) %>%
    set_names(
      c("Discharge", "Death", setdiff(event_types, c("Death", "Discharge")))
    )

  ## -- Gave up on this piece for now; caused inconsistent color matching. -----
  ## -- Stick with what's in the original object, revisit later. ---------------
  ## Label death, discharge with "competing risks"
  # which_competing <- grep("^(Death|Discharge)$", event_types)

  ## When it was just death, used italics for competing risk, but couldn't get
  ## it to work with a flexible # competing risks:
  # labels = c(
  #   expression(paste("Death ", italic("(competing risk)"))), legend_string
  # )
  
  # event_labels <- map_at(
  #   event_types,
  #   which_competing,
  #   ~ paste(., "(competing risk)")
  # )
  # ## Reverse order for proper legends
  # event_labels <- event_labels[length(event_labels):1]
  
  ggcompetingrisks(
    cuminc_obj,
    conf.int = TRUE,
    gnames = sub(" ", "thankscole!", setdiff(names(cuminc_obj), "Tests")),
    gsep = "thankscole!",
    ggtheme = mindusa_theme()
  ) +
    scale_color_manual(name = "", values = risk_colors) +
    scale_fill_manual(name = "", values = risk_colors) +
    scale_y_continuous(
      limits = -0.01:1,
      breaks = seq(0, 1, 0.2),
      labels = paste0(seq(0, 1, 0.2) * 100, "%")
    ) +
    labs(
      title = glue("Cumulative Incidence of {event_string}"),
      subtitle = test_string,
      x = "Days after Randomization",
      y = "Cumulative Probability",
      caption = caption_string
    ) +
    theme(legend.position = c(0.01, 1),
          legend.justification = c(0, 1),
          legend.text.align = 0,
          legend.background = element_blank(),
          legend.direction = "vertical",
          plot.subtitle = element_text(size = basetext_size * 0.75))
}

## Function to combine results of cif_plot, cr_risktable_plot
cr_combine_plots <- function(cif, rt, xmax, time_breaks, caption_string){
  multiplot(
    plotlist = list(
      ## Remove X axis info, caption from CIF plot; only include one at bottom
      cif +
        scale_x_continuous(limits = c(0, xmax), breaks = time_breaks) +
        theme(plot.caption = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank()),
      
      ## Add caption to risk table plot
      rt +
        scale_x_continuous(
          limits = c(0, xmax),
          breaks = time_breaks,
          name = "Days after Randomization"
        ) +
        labs(caption = paste0("\n", caption_string))
    ),
    
    ## Layout: single column, CIF takes up 2/3 space
    layout = matrix(c(1, 1, 2), ncol = 1)
  )
}
