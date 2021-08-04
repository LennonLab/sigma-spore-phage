library(here)
library(tidyverse)
library(lubridate)

# read Synergy time series and convert to table
read_synergy_txt <- function (txt_file){
  # number of times plate was read
  xpt.reads = grep("Reads", readLines(txt_file), value = T)%>%
    str_remove(pattern = ".*,") %>% parse_number()
  
  # start line of OD table
  start.line = grep("^Time.*A1", readLines(txt_file))
  
  # end line of OD table
  aprox.end.line.n= floor(start.line+xpt.reads)-start.line
  
  # read data and organize it
  input <- 
    read_tsv(txt_file,
             skip = start.line-1,
             n_max = aprox.end.line.n-1,
             col_types = cols(Time = "c"))%>%
    # rename temp column
    rename("temp"=2)%>%
    # Convert time to hms
    mutate(Time=lubridate::hms(Time))%>%
    mutate(Time=as.duration(Time))%>%
    # remove empty lines
    filter(Time != 0)
  
  return(input)
}


#---------------------------------------------------

# average time series OD data

smooth_synergy_txt <- function(raw_input, window_size){
  #number of points after averaging
  n.group <-floor( nrow(raw_input)/window_size)
  
  # remove last rows (if needed) so it divides by window_size
  raw_input <- raw_input[1:(window_size*n.group),]
  
  input.avg <- 
    raw_input%>%
    #add grouping variable 
    mutate(g=ntile(1:(window_size*n.group),n.group))%>%
    group_by(g)%>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))%>%
    #remove grouping var
    select(-g) %>%
    # convert time to hours
    mutate(Time = Time/3600)
  
  return(input.avg)
}

