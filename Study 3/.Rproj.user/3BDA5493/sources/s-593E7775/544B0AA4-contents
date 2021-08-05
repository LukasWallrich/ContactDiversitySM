# ------------
# Introduction
# ------------

NAME <- '3_models' ## Name of the R file goes here (without the file extension!)

# ------------
# Sources
# ------------
if (!require("pacman")) install.packages("pacman")
pacman::p_install_version_gh(c("lukaswallrich/rNuggets"),
                             c("0.1.8"))
pacman::p_load(tidyverse, magrittr, here, foreign, survey, srvyr, lavaan, lavaan.survey, rNuggets)


source(here("1_tools", "managementFunctions.R"))


#Set up pipeline folder if missing
pipeline <- createPipelineDir(NAME)
datadir <- "0_data"
pipelinedir <- "3_pipeline"

# ----------------------
# Load data
# ----------------------

imp_long_list <- read_rds(here(pipelinedir, "1_data_prep", "out", "imp_long_list.RDS"))
imp_long_mids <- read_rds(here(pipelinedir, "1_data_prep", "out", "imp_long_mids.RDS"))
ac_svy <- read_rds(here(pipelinedir, "1_data_prep", "out", "ac_svy.RDS"))

imp_long_list_with_dummy <- imp_long_list %>% map(function(x) x$foreign_neighbours %>% forcats::fct_recode(none = "(Almost) no foreigners", some = "Some foreigners", many = "Many foreigners", mostly = "Mostly foreigners") %>% psych::dummy.code() %>% data.frame() %>% select(-none) %>% cbind(x, .))

ac_svy_dummy <- ac_svy$variables$foreign_neighbours %>% forcats::fct_recode(none = "(Almost) no foreigners", some = "Some foreigners", many = "Many foreigners", mostly = "Mostly foreigners") %>% psych::dummy.code() %>% data.frame() %>% select(-none) %>% cbind(ac_svy$variables, .)

my_scale <- function(x) c(scale(x))

names_col <- names(imp_long_list_with_dummy[[1]])

imp_long_list_with_dummy_sd <- imp_long_list_with_dummy %>% map(function(x) x %>% mutate_all(as.numeric) %>% mutate_at(vars(-.id, -.imp, -wt), my_scale))

ac_svy_dummy_sd <- ac_svy_dummy %>% haven::zap_labels() %>% mutate_all(as.numeric) %>% mutate_at(vars(-wt), my_scale)

model <- (' # direct effect
             nb_scoreY ~ some + many + mostly + c1*posCont + c2*negCont + b1*divprefinstr + eastwest +  age  + sex + educN + leftright + b2*for_att

           # mediator
             divprefinstr ~ some + many + mostly + a1*posCont + a2*negCont + eastwest +  age  + sex + educN + leftright
             for_att ~ some + many + mostly + d1*posCont + d2*negCont + eastwest +  age  + sex + educN + leftright

           # direct effect
             direct_pos := c1
             direct_neg := c2

           # indirect effect div (a*b1)
             ind_pos_div := a1*b1
             ind_neg_div := a2*b1

           # indirect effect att (d*b2)
             ind_pos_att := d1*b2
             ind_neg_att := d2*b2

           # total effect
             total_pos := a1*b1 + d1*b2 + c1
             total_neg := a2*b1 + d2*b2 + c2

            total_diff := total_pos + total_neg

            #Pairwise comp indirect effects
            pair_ind_pos := ind_pos_att - ind_pos_div
            pair_ind_neg := ind_neg_att - ind_neg_div
            pair_ind_div := ind_pos_div - ind_neg_div
            pair_ind_att := ind_pos_att - ind_neg_att

         ')

mod_complete <-lavaan::sem(model, ac_svy_dummy_sd)

imp_svy <- svydesign(ids = ~1, weights = imp_long_list_with_dummy_sd[[1]]$wt, data = mitools::imputationList(imp_long_list_with_dummy_sd))

mod_weighted <- lavaan.survey::lavaan.survey(mod_complete, imp_svy)


ind_CIs <- semTools::monteCarloCI(mod_weighted, nRep = 2e4) %>% tibble::rownames_to_column("name") %>%
  rename(est.std = est) %>% filter(str_detect(name, "^ind_|^total|^pair_"))

ind_p <- semTools::monteCarloCI(mod_weighted, return.samples = TRUE, nRep = 2e4) %>% select(matches("^ind_|^total|^pair_")) %>%
  map2_dfr(., names(.), function(x, name) {
    est <- ind_CIs$est.std[ind_CIs$name == name]
    tibble(name = name, pvalue = mean(abs(x-mean(x))>abs(est)), values = list(x))
    })

ind_CIs <- left_join(ind_CIs, ind_p) %>% tibble()


res <- parameterestimates(mod_weighted, ci = TRUE) %>%
  filter(str_detect(lhs, ("^direct"))) %>%
  select(name = lhs, est.std = est, ci.lower, ci.upper, pvalue) %>%
  rbind(ind_CIs %>% select(-values) %>% filter(!name == "total_diff"&!str_detect(name, "^pair_")))

res_tbl <- res %>%
  tidyr::separate(name, c("type", "pred", "mod"), fill = "right") %>%
  mutate(type = coalesce(mod, type), fmt = paste(sprintf("%.2f", round(est.std, 2)), rNuggets::sigstars(pvalue), rNuggets:::.fmt_ci(ci.lower, ci.upper, 2))) %>%
  select(type, pred, fmt) %>%
  spread(type, fmt) %>%
  select(pred, direct, everything(), total) %>%
  arrange(desc(pred)) %>%
  gt::gt() %>%
  gt::cols_label(pred = "Measure") %>%
  gt::tab_spanner(gt::md("**Paths** (std. coefficients)"), 2:ncol(.[["_data"]])) %>%
  gt::fmt_markdown(everything()) %>%
  gt::tab_source_note(gt::md(rNuggets:::.make_stars_note()))

res_tbl %>% gt_apa_style() %>% gt::gtsave(filename = here(pipeline, "out", "mediation.html"))

ind_CIs %>% filter(name == "total_diff"|str_detect(name, "^pair_")) %>%
  mutate(Difference = glue::glue("{round_(est.std)} {rNuggets:::.fmt_ci(ci.lower, ci.upper)}"), Significance = glue::glue("p {fmt_p(pvalue)}")) %>%
  select(name, Difference, Significance) %>%
  gt::gt() %>%
  gt::tab_header("Significance tests for differences between indirect and total effects") %>%
  gt::gtsave(filename = here(pipeline, "out", "differences_ind_effects.html"))

graph_parameters <- parameterestimates(mod_weighted, ci = TRUE) %>% filter(nchar(label)>0)

write_rds(graph_parameters, here(pipeline, "out", "graph_parameters.RDS"))
