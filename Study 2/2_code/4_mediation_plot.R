# ------------
# Introduction
# ------------

NAME <- '4_mediation_plot' ## Name of the R file goes here (without the file extension!)

# ------------
# Sources
# ------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, DiagrammeR, DiagrammeRsvg, here)
pacman::p_install_version_gh(c("lukaswallrich/rNuggets"),
                             c("0.1.9.1000"))

source(here("1_tools", "managementFunctions.R"))


#Set up pipeline folder if missing
pipeline <- createPipelineDir(NAME)
datadir <- "0_data"
pipelinedir <- "3_pipeline"

# ----------------------
# Load data
# ----------------------

graph_parameters <- read_rds(here(pipelinedir, "3_models", "out", "graph_parameters.RDS"))

fmt_param <- function(param) {
  paste0(
  sprintf("%.2f",  graph_parameters$est[graph_parameters$label == param]),
  rNuggets::sigstars(graph_parameters$pvalue[graph_parameters$label == param]))
}

plot_code_edited <- glue::glue("


digraph  {{

graph [layout = 'neato',
outputorder = 'edgesfirst',
bgcolor = 'white', rankdir=LR,]

node [fontname = 'Helvetica',
fontsize = '10',
shape = 'circle',
fixedsize = 'true',
width = '0.5',
style = 'filled',
fillcolor = 'white',
color = 'black',
fontcolor = 'black']



'x1' [label = <Valuing <br /> diversity>, color = 'SaddleBrown', fontcolor = 'SaddleBrown', shape = 'rectangle', height = '0.5', width = '1.5', pos = '2.5,2.4!']
'x2' [label = <Positive attitude <br /> re foreigners>, color = 'SteelBlue', fontcolor = 'SteelBlue', shape = 'rectangle', height = '0.5', width = '1.5', pos = '2.5,0!']


'y' [label = <Approach <br /> intentions>, color = 'black', shape = 'rectangle', height = '0.5', width = '1.5', pos = '5,1.2!']

'M1' [label = <Negative <br /> contact<br />>,
     color = 'black', shape = 'rectangle', height = '0.5', width = '1.5',
     pos = '0,0.5!']
'M2' [label = <Positive <br /> contact<br />>,
     color = 'black', shape = 'rectangle', height = '0.5', width = '1.5',
     pos = '0,1.9!']


'c1' [label = <{fmt_param('c1')}>,
     color = 'black', shape = 'plaintext', fillcolor='transparent',
     pos = '2.4,1.75!']
'c2' [label = <{fmt_param('c2')}>,
     color = 'black', shape = 'plaintext', fillcolor='transparent',
     pos = '2.4,0.6!']

'a1' [label = <{fmt_param('a1')}>,
     color = 'black', shape = 'plaintext', fillcolor='transparent',
     pos = '1.3,2.15!', fontcolor = 'SaddleBrown']
'd2' [label = <{fmt_param('d2')}>,
     color = 'black', shape = 'plaintext', fillcolor='transparent',
     pos = '1.2,0.2!', fontcolor = 'SteelBlue']

'b1' [label = <{fmt_param('b1')}>,
     color = 'black', shape = 'plaintext', fillcolor='transparent',
     pos = '3.9, 1.9!', fontcolor = 'SaddleBrown' ]
'b2' [label = <{fmt_param('b2')}>,
     color = 'black', shape = 'plaintext', fillcolor='transparent',
     pos = '3.91, 0.44!', fontcolor = 'SteelBlue']

'd1' [label = <{fmt_param('d1')}>,
     color = 'black', shape = 'plaintext', fillcolor='transparent',
     pos = '0.92, 1.4!', fontcolor = 'SteelBlue']
'a2' [label = <{fmt_param('a2')}>,
     color = 'black', shape = 'plaintext', fillcolor='transparent',
     pos = '0.92, 1.0!', fontcolor = 'SaddleBrown']





'CV' [label = <<B>Covariates:</B> <BR ALIGN='LEFT' />&#8226; Gender<BR ALIGN='LEFT' />&#8226; Age<BR ALIGN='LEFT' />&#8226; Education<BR ALIGN='LEFT' />&#8226; Region (East/West)<BR ALIGN='LEFT' />&#8226; Political orientation<BR ALIGN='LEFT' />&#8226; Neighborhood diversity<BR ALIGN='LEFT' />>, color = 'black', shape = 'rectangle', height = '1.1', width = '2.2', pos = '3.8,-1!']

'CV1' [style = invis,height = '0', width = '0', pos = '3, -0.5!']
'CV2' [style = invis,height = '0', width = '0', pos = '3.75, -0.5!']
'CV3' [style = invis,height = '0', width = '0', pos = '4.5, -0.5!']

edge [fontname = 'Helvetica',
fontsize = '10',
len = '1.5',
color = 'black',
arrowsize = '0.5']

M1->x1 [tailport = 'e', color = 'SaddleBrown']
M2->x1 [tailport = 'e', color = 'SaddleBrown']
M1->x2 [tailport = 'e', color = 'SteelBlue']
M2->x2 [tailport = 'e', color = 'SteelBlue']
M1->y [tailport = 'e']
M2->y [headport = 'w', tailport = 'e']

x1->y  [style = solid, color = 'SaddleBrown']
x2->y  [style = solid, color = 'SteelBlue']
CV3->y  [style=dashed]
CV2->x1  [style=dashed]
CV1->x2  [style=dashed]

 }}
")

rNuggets:::.grViz_and_save(plot_code_edited, here(pipeline, "out", "mediation_model.svg"))

