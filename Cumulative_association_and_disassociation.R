# libraries #########################
library("RMySQL")
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(qvalue)
library(signnet)
library(graphlayouts) # for multilevel layout

# Functions #######################
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}

### DATA WRANGLING #####################
# GET DATA FROM DB
# Remote DB with password
con <- dbConnect(MySQL(), 
                 user  = "guest",
                 password    = "guest",
                 dbname="meadows",
                 port = 3306,
                 host   = "sxouse.ddns.net")

q <- sprintf('select quadrat_id, survey_id, assembly_name, species.species_name from surveys
      join quadrats on quadrats.survey_id = surveys_id
      join visit_dates on quadrats.vd_id = visit_dates.vds_id
      join records on records.quadrat_id = quadrats_id
      join species on species.species_id = records.species_id where species.species_id != 4 and
      major_nvc_community = "MG5" and quadrat_size = "2x2";')
####### major_nvc_community like "MG%%" and quadrat_size = "2x2";')

# NOTE the double %% to escape the % formatting character

plot_title <- "Meadows DB extract MG5" # Change title when changing extract

rs1 = dbSendQuery(con, q)
the_data <- as_tibble(fetch(rs1, n=-1))
dbDisconnectAll()
rm(con, q, rs1)

# Count (n) hits for each species (columns) in each quadrat (rows)
d <- (the_data %>% select(quadrat_id, species_name)
      %>% group_by(quadrat_id, species_name)
      %>% summarise(n=n())
      %>% ungroup()
      %>% pivot_wider(names_from = species_name, values_from = n))

# At this point, d has either a 1 (for a hit) or NA (for a miss) for each quadrat (rows) and species (columns).
# Replace anything numeric with 1, and any NA with 0
d <- (d %>% select(-quadrat_id) %>% replace(., is.na(.), 0)) # Replace NAs with 0)

n_quadrats <- length(unlist((the_data |> distinct(quadrat_id))))

species_data <- the_data |> distinct(species_name)
species_data <- species_data |> 
  add_column(hits = colSums(d)) |>
  mutate(freq = hits/n_quadrats) |>
  rename(name = species_name) |>
  mutate(jit = jitter(hits))

# What follows thanks to Brian Shalloway: 
# https://www.bryanshalloway.com/2020/06/03/tidy-2-way-column-combinations/#fn4
# Manipulating the data into contingency tables for species pairs (dyads)
# and getting the fisher.test statistics using map not for-loops.

# dfLists: 2 variables, var and vector. Var is species name. Column "vector" is a list of binary vectors for each species 
df_lists <- d %>%
  summarise_all(list) %>% # this is what turns the columns of d into lists.
  pivot_longer(cols = everything(), 
               names_to = "var", 
               values_to = "vector")

# Combine to get all possible combinations of two species. Includes A-B and B-A (upper and lower triangles)
df_lists_comb <- expand(df_lists, # Expand data frame to include all possible combinations of values
                        nesting(var, vector),
                        nesting(var2 = var, vector2 = vector))

# Get rid of lower triangle here?
# Eliminate lower triangle duplicates
# Thank you https://stackoverflow.com/questions/22756392/deleting-reversed-duplicates-with-r
df_lists_comb <- df_lists_comb %>%
  rowwise() %>%
  mutate(grp = paste(sort(c(var, var2)), collapse = "_")) %>%
  group_by(grp) %>%
  slice(1) %>%
  ungroup() %>%
  select(-grp) %>%
  filter(var != var2)                # 0 on main diagonal

is_event <- function(...){ # See map2_int below. If this sum is zero, there are no co-occurrences
  sum(..1 & ..2)
}

df_lists_comb_as <- df_lists_comb %>%
  mutate(events = map2_int(.x = vector, .y = vector2, .f = is_event)) %>%
  filter((events > 0)) # & events < length(df_lists$vector[[1]]))  %>%
# select(-events)
rm(is_event)

# Arrange so that A is the commonest species of the pair
# This is not strictly necessary as odds ratio is symmetric, but could be convenient in future
# if needed the graphs to be directed from most to least frequent species.
df_lists_comb_as <- df_lists_comb_as |>
  rowwise() |>
  mutate(sum1 = sum(vector), sum2 = sum(vector2)) |>
  ungroup()

id <- seq(1:length(df_lists_comb_as$var))
d_list <- map(.x = id, ~ {ifelse(df_lists_comb_as$sum1[.x ] < df_lists_comb_as$sum2[.x], 
                                 temp <- c(df_lists_comb_as$var2[.x], df_lists_comb_as$var[.x], df_lists_comb_as$vector2[.x], df_lists_comb_as$vector[.x]),
                                 temp <- c(df_lists_comb_as$var[.x], df_lists_comb_as$var2[.x], df_lists_comb_as$vector[.x], df_lists_comb_as$vector2[.x]))
  return(temp)})

d2 <- as_tibble(do.call(rbind, d_list))
dyads <- tibble(from = unlist(d2$V1), to = unlist(d2$V2), vector = d2$V3, vector2 = d2$V4)
## rm(d_list, d2, df_lists, df_lists_comb, df_lists_comb_as, id)

# function to make the contingency tables.
xtab <- function(...){ # See map2 below
  t <- table(..1, ..2) 
}
# And map them into dyads.
dyads <- dyads %>% 
  mutate(xtable = map2(vector, vector2, xtab))

# Important function to SAFELY apply fisher.test, ensuring we get
# some return value for dyad.
sft <- safely(.f = fisher.test, otherwise = NULL, quiet = FALSE)

# Map fisher.test over dyads
dyads <- dyads |> mutate(f_test = map(xtable, sft))

# Temporary data structures needed to extract fisher.test statistics
# from <hresult> class and transfer them to edges (dyads) data where they belong.
x <- tibble(pval = rep(0, length(dyads$from)))
x <- x %>% mutate(odds_ratio = pval)

# Extract p.val and lor from fisher.test <hresult>
# Would prefer to do using map.
for(i in 1:length(dyads$from)) { 
  x$pval[i] <- (ifelse(!is.null(dyads$f_test[[i]][["error"]]), NA, dyads$f_test[[i]][["result"]][["p.value"]]))
  x$odds_ratio[i] <- (ifelse(!is.null(dyads$f_test[[i]][["error"]]), NA, dyads$f_test[[i]][["result"]][["estimate"]]))
}
dyads <- dyads |> add_column(pval = x$pval, lor = log10(x$odds_ratio))
dyads <- dyads |> select(-3, -4, -6) # Might need the xtables (5) again
dyads <- dyads |> filter(pval < 1)   # needed to make qval work???
# Clean up...
rm(i, x, sft, xtab)

# FDR adjustment
q_obj <- qvalue(dyads$pval) #, pi0 = 0.5)
plot(hist(q_obj))
dyads <- dyads |>
  add_column(lfdr = q_obj$lfdr) |>
  filter(lfdr < 0.01)


# Make the signed graph
g <- as_tbl_graph(graph.data.frame(d = dyads, directed = FALSE))

# Add species names to edges from - to
g <- g %>% activate(edges) %>% mutate(A = .N()$name[from], B = .N()$name[to])

# Make g signed
g <- g |> activate(edges) |> mutate(sign = ifelse(lor > 0, 1, -1))

#### Negative edge count as species added ################
### Add species, starting with the least common #####
g <- g |>
  activate(nodes) |>
  left_join(species_data) # Add hits and jittered hits

thresh <- g |> 
  activate(nodes) |>
  as_tibble()

# Increase the threshold so more species added on each iteration.
# collect cumulative disassociative edge count
balance <- map_df(thresh$jit, ~{
  g1 <- g |> activate(nodes) |> filter(jit <= .x) 
  df <- g1 |> activate(edges) |> as_tibble() |> group_by(sign) |> count(sign)
  return(data.frame(.x, negative = df$n[1], positive = df$n[2], res = nrow(df)))
})

balance <- balance %>%
  filter(res == 2) |>
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))

p0 <- balance |> ggplot(aes(.x, negative)) + 
  geom_point() + 
  ggtitle(paste(plot_title, "Cumulative disassociative dyads", sep = ", ")) +
  labs(x = "Hits", y = "Disassociative edge count")
  
plot(p0)

p1 <- balance |> ggplot(aes(.x, positive)) + 
  geom_point() + 
  ggtitle(paste(plot_title, "Cumulative associative dyads", sep = ", ")) +
  labs(x = "Hits", y = "Associative edge count")

plot(p1)

