# ------------
# Introduction
# ------------

## Create scales and descriptives

NAME <- '1_scales'

# ------------
# Sources
# ------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, magrittr, here, readr, lavaan)
pacman::p_install_version(c("haven", "modelsummary", "stringr", "purrr", "gt"),
                          c("2.3.1", "0.5.1", "1.4.0", "0.3.4", "0.2.1"))
pacman::p_install_version_gh(c("lukaswallrich/rNuggets", "easystats/report"),
                             c("0.1.8", "0.1.0"))

source(here("1_tools/managementFunctions.R"))

#Set up pipeline folder if missing
pipeline <- createPipelineDir(NAME)
datadir <- "0_data"
pipelinedir <- "3_pipeline"

AllCases <- haven::read_sav(here(datadir, "CombinedData.sav"))

notes <- character()

# ------------
# Create scales
# ------------

AllCases %<>% mutate_at(c("Ethn"), haven::as_factor)

AllCases %<>% dplyr::filter(Ethn == 'White British' | Ethn == 'Other White')

scalesList <- list("BlackPOS_1" = c("BlackComplimented", "BlackBefriended", "BlackWelcome",
                                  "BlackSupported", "BlackHelped"),
                  "BlackPOS_2" =  c("BlackComplimented_2", "BlackBefriended_2",
                                    "BlackWelcome_2", "BlackSupported_2", "BlackHelped_2"),
                  "BlackNEG_1" = c("BlackRidiculed", "BlackUnwanted", "BlackVerbAbused",
                                   "BlackIntimidated", "BlackThreatened"),
                  "BlackNEG_2" = c("BlackRidiculed_2", "BlackUnwanted_2", "BlackVerbAbused_2",
                                   "BlackIntimidated_2", "BlackThreatened_2"),
                  "ValDiv_1" = c("ValDiv1", "ValDiv2", "ValDiv3"),
                  "ValDiv_2" = c("ValDiv1_2", "ValDiv2_2", "ValDiv3_2"))


scales <- rNuggets::make_mult_scales(AllCases, scalesList, print_hist = F, print_desc = F)

AllCases <- cbind(AllCases, scales$scores)


scales$descriptives %>% rNuggets::round_df(3) %>% gt::gt() %>% gt::tab_header("Scale descriptives") %>%
  gt::tab_source_note("Reliabilities are Cronbach's alpha for scales with more than two items, Spearman-Brown for scales with two items") %>%
  gt::gtsave(here(pipeline, "out/scale_descriptives.html"))


write_rds(AllCases, here(pipeline, "out/WhiteCasesWScales.RDS"))

# ------------
# Assess & desribe missing data
# ------------

# Testing whether drop-out was systematic, i.e. related to any of the other variables
AllCases$drop_out <- F
AllCases$drop_out[!AllCases$BothWaves] <- T

notes %<>% c(paste("Summary of participant drop-out: ", report::report(AllCases$drop_out)))

mod <- glm(drop_out ~ Age + Gender + src + BlackPOS_1 + BlackNEG_1 + ValDiv_1, data = AllCases, family = "binomial")

modelsummary::modelsummary(mod, output="gt", statistic = "p.value", stars = rNuggets:::std_stars, title = "Odds ratios and p-values for logistic regression predicting drop_out", exponentiate = TRUE) %>% gt::gtsave(here(pipeline, "out/drop_out_logreg.html"))



# ------------
# Calculate descriptives (with max likelihood)
# ------------


notes %<>% c(AllCases %>% mutate(Gender = ifelse(Gender == 1, "M", "F") %>% factor()) %>% report::report_participants(age = "Age", sex = "Gender"))

desc_mod <- ("
            #Variances
            BlackPOS_1 ~~ BlackPOS_1
            BlackPOS_2 ~~ BlackPOS_2
            BlackNEG_1 ~~ BlackNEG_1
            BlackNEG_2 ~~ BlackNEG_2
            ValDiv_1 ~~ ValDiv_1
            ValDiv_2 ~~ ValDiv_2

             #Co-cariances/correlations
             BlackPOS_1 ~~ BlackPOS_2 + BlackNEG_1 + BlackNEG_2 + p1d1*ValDiv_1 + p1d2*ValDiv_2
             BlackPOS_2 ~~ BlackNEG_1 + BlackNEG_2 + p2d1*ValDiv_1 + p2d2*ValDiv_2
             BlackNEG_1 ~~ BlackNEG_2 + n1d1*ValDiv_1 + n1d2*ValDiv_2
             BlackNEG_2 ~~ n2d1*ValDiv_1 + n2d2*ValDiv_2
             ValDiv_1 ~~ ValDiv_2

             #Pos-neg abs differences
             d11 := p1d1 + n1d1
             d22 := p2d2 + n2d2
             d21 := p1d2 + n1d2
             d12 := p2d1 + n2d1
             ")



mod <- sem(desc_mod, AllCases, missing = "fiml", meanstructure = TRUE)

Ms <- parameterestimates(mod) %>%
  filter(op == "~1") %>%
  select(var = lhs, M = est)

desc <- parameterestimates(mod) %>%
  filter(op == "~~" & lhs == rhs) %>%
  transmute(var = lhs, SD = sqrt(est)) %>%
  left_join(Ms, ., by = "var")


cors <- standardizedsolution(mod, ci = TRUE) %>%
  filter(op == "~~" & lhs != rhs) %>%
  transmute(lhs, rhs, cors = est.std, p.values = pvalue, ci.low = ci.lower, ci.high = ci.upper)

write.csv(cors, here(pipeline, "out/cors.csv"))


named_matrix <- matrix(rep(0, 36), nrow = 6) %>% set_rownames(union(cors$lhs, cors$rhs) %>% unique()) %>% set_colnames(union(cors$lhs, cors$rhs) %>% unique())

cor_matrix <- purrr::map(3:6, function(x) {
  named_matrix[lower.tri(named_matrix)] <- cors[[x]]
  named_matrix
}) %>% set_names(names(cors)[3:6])

cor_matrix[["desc"]] <- desc

rNuggets::apa_cor_table(cor_matrix, ci = "given", filename = here(pipeline, "out/correlations.html"))

 parameterestimates(mod) %>%
  filter(op == ":=") %>% select(comparison = lhs, everything(), -c(op, rhs, label, se, z)) %>% rNuggets::round_df(3) %>% gt::gt() %>% gt::tab_source_note("Significance tests for asymetry in correlation coefficients between pos/neg contact and diversity. The numbers indicate timepoints, e.g., d21 stands for the correlation of diversity at T2 with contact measures at T1") %>% gt::gtsave(here(pipeline, "out/sig_tests_asymetries.html"))

 # ------------
 # Save notes
 # ------------

 writeLines(notes, here(pipeline, "out/notes.txt"))
