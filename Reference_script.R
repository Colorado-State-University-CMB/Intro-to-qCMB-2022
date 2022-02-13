# Section 2

coronavirus %>%
  filter(type != "recovery") %>%
  group_by(type, continent_name) %>%
  summarise(cases = sum(cases), .groups = "drop") %>%
  mutate(type = factor(type, levels = c("confirmed", "death"))) %>%
  pivot_wider(names_from = type, values_from = cases) %>%
  mutate(death_rate = death / confirmed) %>%
  filter(!is.na(continent_name)) %>%
  arrange(-death_rate) %>%
  datatable(rownames = FALSE,
            colnames = c("Continent", "Confrimed Cases", "Death Cases","Death Rate %")) %>%
  formatPercentage("death_rate", 2)

# Simplified
coronavirus %>%
  filter(type != "recovery") %>%
  group_by(type, country) %>%
  summarise(cases = sum(cases)) %>%
  pivot_wider(names_from = type, values_from = cases) %>%
  mutate(death_rate = death / confirmed) %>%
  datatable(rownames = FALSE,
            colnames = c("Country", "Confrimed Cases", "Death Cases","Death Rate %")) %>%
  formatPercentage("death_rate", 0)