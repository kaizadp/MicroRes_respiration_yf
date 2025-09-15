
# load packages 
library(tidyverse)

# import files and layout ----
import_plate_results <- function(FILEPATH){
  
  filepaths = list.files(path=FILEPATH, pattern = c(".xlsx"), full.names = TRUE, recursive = TRUE)
  do.call(bind_rows, lapply(filepaths, function(path){
    df <- readxl::read_excel(path, skip = 23) %>% 
      mutate(file = basename(path),
             source = path
      )
    df
  }))
}
import_layout <- function(FILEPATH){
  
  filepaths = list.files(path=FILEPATH, pattern = c(".xlsx"), full.names = TRUE, recursive = TRUE)
  
  inoculation = do.call(bind_rows, lapply(filepaths, function(path){
    df <- readxl::read_excel(path, sheet = "inoculation") %>% 
      mutate(file = basename(path)
      )
    df
  }))
  sources = do.call(bind_rows, lapply(filepaths, function(path){
    df <- readxl::read_excel(path, sheet = "sources") %>% 
      mutate(file = basename(path)
      )
    df
  }))
  
  
  x = 
    bind_rows(inoculation %>% mutate(details = "inoculation"),
                sources %>% mutate(details = "carbon_source")) %>% 
    rename(letter = `...1`) %>% 
    pivot_longer(cols = c(-letter, -file, -details),
                 names_to = "number") %>% 
    mutate(number = str_pad(number, 2, pad = "0"),
           well = paste0(letter, number)) %>% 
    dplyr::select(-letter, -number) %>% 
    pivot_wider(names_from = "details")
  
  x
  
}

df = import_plate_results(FILEPATH = "1-data/plate_results")
layout = import_layout(FILEPATH = "1-data/layout")


# process the layout
layout_processed = 
  layout %>% 
  mutate(file = str_remove(file, ".xlsx")) %>% 
  separate(file, sep = "_", into = c("experiment", "date")) %>% 
  mutate(date = ymd(date),
         strain = case_when(grepl("Control", inoculation, ignore.case = T) ~ "Control",
                            grepl("Str", inoculation, ignore.case = T) ~ str_extract(inoculation, "Str_[0-9]{3}"))) %>% 
  drop_na()

# process the results
df_processed = 
  df %>% 
  mutate(date = str_extract(source, "[0-8]{8}"),
         date = ymd(date)) %>% 
  dplyr::select(-`...14`, -source) %>% 
  rename(letter = `...1`) %>% 
  pivot_longer(cols = -c(letter, file, date),
               names_to = "number",
               values_to = "abs_570") %>% 
  mutate(number = str_pad(number, 2, pad = "0"),
         well = paste0(letter, number)) %>% 
  dplyr::select(-letter, -number) %>% 
  mutate(file = str_remove(file, ".xlsx")) %>% 
  separate(file, sep = " ", into = c("experiment", "plate")) %>% 
  separate(plate, sep = "_", into = c("timepoint_hr", "plate")) %>% 
  mutate(timepoint_hr = parse_number(timepoint_hr)) %>% 
  # merge the layout
  # left_join(layout_processed) %>% 
  force()

# do calculations ----

## 1. normalize the data
## normalized data (Ai) = (At6/At0) * mean (At0)

mean_t0 = 
  df_processed %>% 
  filter(timepoint_hr == 0) %>% 
  group_by(experiment, plate, date) %>% 
  dplyr::summarise(mean_t0 = mean(abs_570))

t0 = 
  df_processed %>% 
  filter(timepoint_hr == 0) %>% 
  rename(t0 = abs_570) %>% 
  dplyr::select(-timepoint_hr)

abs_calculated = 
  df_processed %>% 
  filter(!timepoint_hr == 0)  %>% 
  left_join(t0) %>% 
  left_join(mean_t0) %>% 
  mutate(abs_relative = (abs_570/t0) * mean_t0)

abs_summary = 
  abs_calculated %>% 
  left_join(layout_processed) %>% 
  group_by(experiment, timepoint_hr, plate, date, inoculation, strain, carbon_source) %>% 
  dplyr::summarise(abs_mean = mean(abs_relative)) %>% 
  ungroup() %>% 
  mutate(CO2_percent = -0.2265 +(-1.606/(1-6.771*abs_mean)),
         CO2_percent2 = A + B/(1 + (D*abs_mean)),
         # calculate CO2 ug in headspace
         # 44 = molar mass of CO2. 1 mole -- 44 g, or 1 umole -- 44 ug
         # 22.4 = volume occupied by one mole of CO2/ideal gas. 1 mole -- 22.4 L, or 1 umole -- 22.4 uL
         CO2_ugC = (CO2_percent/100) * headspace_uL * (44/22.4) * (12/44) * (273/(273 + temp)),
         CO2_ugC_uL = CO2_ugC / 500,
         CO2_ugC_uL_hr = CO2_ugC_uL/6,
         CO2_ngC_uL_hr = CO2_ugC_uL_hr * 1000
         ) %>% 
  filter(!is.na(carbon_source))

A = -0.2265
B = -1.606
D = -6.771
headspace_uL = 945
temp = 20

## graph
abs_summary %>% 
  ggplot(aes(x = strain, y = CO2_ngC_uL_hr, color = strain, shape = plate))+
  geom_point(size = 3,
             position = position_dodge(width = 0.4))+
  facet_wrap(~carbon_source)


l = lm(CO2_ngC_uL_hr ~ strain + carbon_source + plate, data = abs_summary)
a = car::Anova(l)
a
