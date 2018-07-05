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

## -- Wrapper functions for facet axis breaks (thanks, Jim Hester!) ------------
cif_breaks <- function(x) if (max(x) <= 40) seq(0, 30, 5) else seq(0, 90, 15)
tteem_breaks <- function(x){
  sort(unique(c(1, round(seq(min(x), max(x), length.out = 3), 0.5))))
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
      .(kw_obj$parameter) * "; P," ~ .(formatp_nejm(kw_obj$p.value))
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
      subtitle = glue("P: {formatp_nejm(anova(mod)['trt', 'P'])}"),
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
mental_medians_plot_em <- function(median_df, em_string, title_string = NULL){
  p <- ggplot(data = median_df, aes(y = quantile, x = trt_short)) +
    facet_grid(em_text ~ outcome_text) +
    geom_pointrange(
      aes(ymin = lb, ymax = ub),
      colour = as.character(palette_colors["dred"]), size = 0.5
    ) +
    scale_y_continuous(
      limits = c(0, 14), breaks = c(0, 3, 7, 10, 14),
      name = "Adjusted Median (95% Confidence Interval)"
    ) +
    coord_flip() +
    labs(
      subtitle =
        glue("P-values test total {em_string} * treatment interaction."),
      caption = glue(
        "Adjusted analysis using proportional odds logistic regression.\n",
        "P-values test total {em_string} * treatment interaction."
      )
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text = element_text(size = basetext_size * 0.7),
      panel.spacing.y = unit(0.5, "cm"),
      strip.text.x = element_text(vjust = 0),
      plot.subtitle = element_text(face = "italic"),
      plot.caption = element_text(face = "italic")
    )
  
  if(!is.null(title_string)){
    p <- p + labs(title = title_string)
  }
  
  return(p)
  
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
      subtitle = glue("P: {formatp_nejm(anova(mod)['trt', 'P'])}"),
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
  title_string = NULL, ## Text for plot title
  ratio_type = c("Hazard", "Odds"),
  facet_formula = "em_text ~ outcome_text" ## string specifying how to facet
){
  
  if(ratio_type == "Hazard"){
    caption_text <- glue(
      "Adjusted mortality outcomes use standard Cox proportional hazards regression.\n",
      "Other outcomes use Fine-Gray competing risks regression, with competing risk of death\n",
      "(and ICU discharge without the event for MV liberation and ICU readmission).\n\n",
      "P-values test total {em_string} * treatment interaction.\n\n",
      "NOTE: MV liberation includes only the {sum_na(ptsummary_df$on_mv_rand24)}",
      " patients on any type of MV at or within 24h of randomization.\n",
      "ICU readmission includes only the {sum_na(ptsummary_df$elig_readm)} ",
      "patients who survived and were discharged from their first ICU stay.\n",
      "All other analyses include all {n_rand} randomized patients."
    )
  } else{
    caption_text <- glue(
      "Adjusted analysis using proportional odds logistic regression.\n",
      "P-values test total {em_string} * treatment interaction."
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
    # scale_y_continuous(breaks = tteem_breaks) + ## didn't actually go so well
    labs(
      subtitle =
        glue("P-values test total {em_string} * treatment interaction."),
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
      plot.subtitle = element_text(face = "italic"),
      plot.caption = element_text(size = basetext_size * 0.7)
    ) +
    coord_flip()
  
  if(!is.null(title_string)){
    p <- p + labs(title = title_string)
  }
  
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
    
    # linetype = "strata",
    
    ## CIs
    conf.int = TRUE, conf.int.alpha = 0.15,
    
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
      name = "Probability of Overall Survival (%)",
      breaks = seq(0, 1, 0.25),
      labels = paste0(seq(0, 1, 0.25) * 100, "%")
    ) +
    labs(
      title = glue("Kaplan-Meier Curve, {time}-Day All-Cause Death"),
      ## Add log-rank test results to subtitle, using bquote for superscript/
      ## inclusion of variables
      subtitle = bquote(
        "Log-rank test for difference between treatments:" ~ X^2 ~ ", " ~
          .(rndformat(chis, 1)) * "; df, " ~ .(df) * "; P, " ~ .(formatp_nejm(pval))
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
      "not fully met (p: {formatp_nejm(global_phtest)}); please see the figures ",
      "below for further detail."
    ),
    glue(
      "The global test for proportional hazards indicates no major concern ",
      "(p = {formatp_nejm(global_phtest)}); please see the figures below for ",
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
      .(df$subdist[1]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[1], 2)) * "; df," ~ .(df$df[1]) * "; P," ~ .(formatp_nejm(df$pv[1])) * "." ~ .(df$subdist[2]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[2], 2)) * "; df," ~ .(df$df[2]) * "; P," ~ .(formatp_nejm(df$pv[2])) * "."
    )
  } else if(nrow(df) == 3){
    bquote(
      .(df$subdist[1]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[1], 2)) * "; df," ~ .(df$df[1]) * "; P," ~ .(formatp_nejm(df$pv[1])) * "." ~ .(df$subdist[2]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[2], 2)) * "; df," ~ .(df$df[2]) * "; P," ~ .(formatp_nejm(df$pv[2])) * "." ~ .(df$subdist[3]) * ":" ~ X^2 * "," ~ .(rndformat(df$stat[3], 2)) * "; df," ~ .(df$df[3]) * "; P," ~ .(formatp_nejm(df$pv[3])) * "."
    )
  }
  
}

## -- Functions to extract and plot the N at risk, cumulative events given -----
## -- a summary(survfit(...)) object -------------------------------------------
cr_risktable <- function(
  sf_sum,       ## summary(survfit(...))
  main_event ## character string; one of the event types in sf_sum
){
  ## Get character strings of event types
  event_types <- colnames(sf_sum$p0)[1:(ncol(sf_sum$n.event) - 1)]
  
  ## -- Create dummy dataset for merging treatments/times ----------------------
  ## Get all time points specified in any group in sf_sum
  sf_times <- unique(sf_sum$time)
  ## Want N, etc at risk at all specified time points; LOCF time points in each
  ## treatment group if all pts are out of risk set
  dummy_times <- cross_df(
    list(
      trt = unique(str_replace(sf_sum$strata, "trt=", "")),
      time = sort(unique(sf_sum$time))
    )
  )

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
    ## Add Ns for each event type
    bind_cols(
      as.data.frame(sf_sum$n.event[, 1:(ncol(sf_sum$n.event) - 1)]) %>%
        set_names(event_types)
    ) %>%
    ## Fill in values for times after last time available in sf_sum
    right_join(dummy_times, by = c("trt", "time")) %>%
    arrange(trt, time) %>%
    group_by(trt) %>%
    fill(-trt, -time) %>%
    ## Calculate *cumulative* Ns at each time point
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
             sprintf("%s\n(%s)\n(%s)", n_risk, primary, Death))
    })
  }

  return(df)
  
}

cr_risktable_plot <- function(
  sf_sum,              ## survfit() object
  main_event,          ## main event of interest as represented in ftype variable
  event_string,        ## Text for main event
  text_size = 0.2,     ## Multiplier for plot text size (default = in report)
  order_trts = TRUE    ## Whether to order treatments the way we really want
    ## If using cif_plot()/ggcompetingrisks, probably want to use FALSE
){
  df <- cr_risktable(sf_sum, main_event = main_event)

  if(order_trts){
    ## Put plot in correct order - BUT can't quickly figure this out for
    ##  ggcompetingrisks, so leave it as an option
    df$trt <- factor(df$trt, levels = c("Placebo", "Haloperidol", "Ziprasidone"))
  }  
  
  ## Set plot title
  plot_title <- ifelse(
    "Discharge" %in% names(df),
    glue(
    "Number at risk (cumulative number {event_string}) (deaths; discharges)"
    ),
    glue(
      "Number at risk (cumulative number {event_string}) (deaths)"
    )
  )
  
  ggplot(data = df, aes(x = time, y = 1, colour = trt)) +
    facet_wrap(~ trt, nrow = 1) +
    geom_text(
      aes(label = risk_string),
      size = basetext_size * text_size, family = basetext_family,
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

## -- Create plot for cumulative incidence functions with ggcompetingrisks -----
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

################################################################################
## ggcompetingrisks() (and cuminc()) cut off graph at the final observed time,
## even if it's well before the end of the time period. We need some data
## management to extend our data to the very end, then plot manually.
################################################################################

## -- Goal: Data set w/ one record per time + trt, one column per outcome ------
## -- Function to do dataset prep for one outcome ------------------------------
## Arguments can take *lists* - intended to be used with purrr::pmap() functions
prep_cuminc_df <- function(
  cuminc_obj,    ## cmprsk::cuminc() object
  censor_time,   ## administrative censoring time
  cph_mod = NULL ## cph() model object, if p-value desired
){
  
  if(!is.null(cph_mod)){
    ## Extract p-value for treatment, if cph() model object provided
    p_trt <- anova(cph_mod)["trt", "P"]
  } else{
    p_trt <- NA
  }
  
  ## We want to extract data from all elements of the cuminc() object except
  ##  the test string
  cif_names <- setdiff(names(cuminc_obj), "Tests")
  
  ## 1. Create data.frame out of time, est, var components
  df <- map(
    cuminc_obj[names(cuminc_obj) != "Tests"],
    bind_cols
  ) %>%
    ## For each df, add record for time of administrative censoring if not present
    map(
      ~ bind_rows(., tail(., n = 1) %>% mutate(time = censor_time)) %>% unique()
    ) %>%
    ## Add primary outcome, treatment, overall outcome to each df;
    ##  combine into one big df
    map2_df(
      .y = cif_names,
      ~ .x %>%
        mutate(
          trt_curve = .y,
          tte_outcome = rownames(cuminc_obj$Tests)[1],
          p_trt = p_trt
        )
    ) %>%
    ## Separate outcome, treatment variables
    ##  extra = "merge": keep everything before first space separate; merge
    ##   everything else into "subdist" column
    separate(
      trt_curve, into = c("trt", "subdist"), sep = " ", extra = "merge"
    ) %>%
    ## Add CIs
    mutate(
      lb = est - qnorm(0.975) * sqrt(var),
      lb = ifelse(lb < 0, 0, lb),
      ub = est + qnorm(0.975) * sqrt(var),
      ## Order factors
      trt = factor(trt, levels = c("Placebo", "Haloperidol", "Ziprasidone")),
      outcome_order = case_when(
        tte_outcome == "ICU Discharge" ~ 1,
        tte_outcome == "Hospital Discharge" ~ 2,
        tte_outcome == "Readmission" ~ 3,
        TRUE ~ 4
      ),
      outcome_label = LETTERS[outcome_order],
      tte_outcome2 = case_when(
        tte_outcome == "MV Liberation" ~ "Liberation from MV",
        tte_outcome == "Readmission" ~ "ICU Readmission",
        TRUE ~ tte_outcome
      ),
      ## Two-line, one-line versions of outcome + p-values
      p_string = ifelse(
        !is.na(p_trt),
        glue("P{ifelse(p_trt < 0.001, , '=')}{formatp_nejm(p_trt)}"),
        ""
      ),
      tte_outcome_2l = glue(
        "{tte_outcome2}{ifelse(!is.na(p_trt), '\n\n', '')}{p_string}"
      ),
      tte_outcome_2l = fct_reorder(tte_outcome_2l, .x = outcome_order),
      tte_outcome_1l = glue(
        "{outcome_label}: ",
        "{tte_outcome2}{ifelse(!is.na(p_trt), '; ', '')}{p_string}"
      ),
      tte_outcome_1l = fct_reorder(tte_outcome_1l, .x = outcome_order),
      ## Indicator for whether this is the subdistribution of primary interest
      is_subdist = tte_outcome == subdist
    )
  
  return(df)
  
}

## -- Function to "manually" create CIF curves with all events -----------------
## ggcompetingrisks() doesn't extend all the way to the end of the time period,
##  if all patients are out of the risk set beforehand; also can't figure out
##  how to order facets appropriately
create_cif_manual <- function(
  cuminc_obj,            ## cmprsk::cuminc() object
  main_event,            ## string indicating primary event of interest
  event_string,          ## text to describe primary event of interest
  # color_bytrt = FALSE, ## Change colors for main curve by treatment?
  event_color = palette_colors[["dred"]],
  caption_string = NULL, ## text for caption
  ## color for curve of main event, if same for all treatments
  test_string = NULL,    ## to include in subtitle; often result of cuminc_test()
  cph_mod = NULL         ## optional: Cox model fit; will extract P-value
){
  ## Not currently an option
  # if(color_bytrt){
  #   cif_colors <- mindusa_trtcols_long
  # } else{
  ## Graphics prep: named character vector with color scheme
  ## ggplot2 will get confused if we have two things named ICU Discharge
  # if(main_event == "ICU Discharge"){
  #   cif_colors <- c(
  #     event_color, palette_colors[["dgray"]]
  #   ) %>%
  #     set_names(c(main_event, "Death"))
  # } else{
  cif_colors <- c(
    event_color, palette_colors[["dgray"]], palette_colors[["lgray"]]
  ) %>%
    set_names(c(main_event, "Death", "Discharge"))
  # }
  # }
  
  ## Data preparation: Need to
  ## - Put all data in a data.frame for plotting
  ## - Create a final record for the end of the time frame, if it doesn't exist
  
  ## All TTE outcomes censored at 90 days except MV liberation
  cens_time <- case_when(
    main_event == "MV Liberation" ~ 30.01,
    TRUE ~ 90.01
  )
  
  cuminc_data <- prep_cuminc_df(
    cuminc_obj = cuminc_obj,
    censor_time = cens_time,
    cph_mod = cph_mod
  ) %>%
    ## Factor-ize subdist so legend displays in desired order
    mutate(
      subdist_order = case_when(
        subdist %in% c(
          "ICU Discharge", "Hospital Discharge", "MV Liberation", "Readmission"
        ) ~ 1,
        subdist == "Death" ~ 2,
        TRUE ~ 3
      ),
      subdist = fct_reorder(subdist, subdist_order)
    )
  
  ## Create string for main event in plot title
  title_string <-
    paste(capitalize(str_split(event_string, " ")[[1]]), collapse = " ")
  if(is.null(test_string)){
    test_string <- sprintf(
      "P for treatment, %s: %s",
      event_string,
      formatp_nejm(unique(cuminc_data$p_trt))
    )
  }
  
  ## CIF curves for all events
  ## Competing risks will be grayscale; event of interest will be a color
  gg_cif <-
    ggplot(data = cuminc_data, aes(x = time, y = est, group = subdist)) +
    facet_wrap(~ trt, nrow = 1) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = subdist), alpha = 0.15) +
    geom_step(aes(color = subdist)) +
    scale_color_manual(values = cif_colors) +
    scale_fill_manual(values = cif_colors) +
    scale_x_continuous(breaks = cif_breaks) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = paste0(seq(0, 1, 0.2) * 100, "%"),
      name = "Probability of Event (%)"
    ) +
    labs(
      title = glue(
        "Cumulative Incidence of {title_string} and Competing Risks"
      ),
      subtitle = test_string
    ) +
    theme(
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.position = c(0.01, 0.99),
      legend.direction = "vertical",
      legend.justification = c(0, 1),
      plot.title = element_text(size = basetext_size * 1.1),
      axis.title.x = element_blank()
    )
  
  return(gg_cif)
  
}

## Function to combine results of cif_plot/create_cif_manual, cr_risktable_plot
cr_combine_plots <- function(cif, rt, xmax, time_breaks, caption_string = NULL){
  ## Add caption to risk table, if specified
  if(!is.null(caption_string)){
    rt <- rt +
      labs(caption = paste0("\n", caption_string))
  }
  
  xlims <- c(max(time_breaks) - xmax, xmax)
  
  multiplot(
    plotlist = list(
      ## Remove X axis info, caption from CIF plot; only include one at bottom
      cif +
        scale_x_continuous(limits = xlims, breaks = time_breaks) +
        theme(plot.caption = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              legend.position = c(0.01, 0.98)),
      
      ## Add caption to risk table plot
      rt +
        scale_x_continuous(
          limits = xlims,
          breaks = time_breaks,
          name = "Days after Randomization"
        )
    ),
    
    ## Layout: single column, CIF takes up 2/3 space
    layout = matrix(c(1, 1, 2), ncol = 1)
  )
}
