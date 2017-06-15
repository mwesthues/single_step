if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "pophelper")

# STRUCTURE files (do not use this command to read local files)
sfiles <- list.files(
  path = "./data/processed/maizego",
  pattern = "popstruc[0-9]",
  full.names = TRUE
  )
# basic usage
slist <- readQ(
  files = sfiles,
  indlabfromfile = TRUE,
  filetype = "structure"
  )

# Bar plots
slist %>% 
  plotQ(
    qlist = .,
    imgoutput = "join"
    )

# Multiline plots
slist %>% 
  .[1] %>% 
  plotQMultiline()

# Evanno
slist %>% 
  tabulateQ() %>% 
  summariseQ() %>% 
  evannoMethodStructure(
    data = .,
    exportplot = TRUE,
    writetable = TRUE,
    na.rm = TRUE
  )
