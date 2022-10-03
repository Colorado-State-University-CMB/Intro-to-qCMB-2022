# Section 1.6
total_cases <- coronavirus %>%
  filter(type != "recovery") %>% # "!=" (not equal to)
  group_by(type) %>%
  summarise(cases = sum(cases)) %>%
  mutate(type = factor(type, levels = c("confirmed", "death"))) 
# Factors are used to represent categorical data which cannot be done with a character vector

       # Simplified
total_cases <- coronavirus %>%
  filter(type != "recovery") %>% # "!=" (not equal to)
  group_by(type) %>%
  summarise(cases = sum(cases)) 




# Section 2
coronavirus %>%
  filter(type != "recovery") %>%
  group_by(type, continent_name) %>%
  summarise(cases = sum(cases), .groups = "drop") %>%
  mutate(type = factor(type, levels = c("confirmed", "death"))) %>%
  pivot_wider(names_from = type, values_from = cases) %>%
  mutate(death_rate = death / confirmed) %>%
  filter(!is.na(continent_name)) %>%
  arrange(-confirmed) %>%
  datatable(rownames = FALSE,
            colnames = c("Continent", "Confrimed Cases", "Death Cases","Death Rate %")) %>%
  formatPercentage("death_rate", 2)

       # Simplified
coronavirus %>%
  filter(type != "recovery") %>%
  group_by(type, country) %>%
  summarise(cases = sum(cases)) %>%
  mutate(type = factor(type, levels = c("confirmed", "death"))) %>%
  pivot_wider(names_from = type, values_from = cases) %>%
  mutate(death_rate = death / confirmed) %>%
  datatable(rownames = FALSE,
            colnames = c("Country", "Confrimed Cases", "Death Cases","Death Rate %")) %>%
  formatPercentage("death_rate", 0)









