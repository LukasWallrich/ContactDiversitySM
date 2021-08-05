if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, purrr)

set.seed(123456)

files <- list.files(here("2_code"), pattern = "\\.R$")[-1]

map(files, function(x) {
  message(crayon::white(crayon::bgGreen("Now running ", x)))
  source(here("2_code", x))})

notes <- character()

notes <- c(notes, "Last complete run:", timestamp(), "\n\n\n", capture.output(sessionInfo()))

writeLines(notes, here("2_code", "last_complete_run.txt"))

