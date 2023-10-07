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

### EXAMINE SIGNED TRIANGLES ################
triangles <- count_signed_triangles(g)|>as_tibble() |> mutate(type=c('+++', '++-', '+--', '---'))
triangles <- triangles |> mutate(balance=ifelse(type %in% c('+++', '+--'), 'balanced', 'unbalanced'))

# EXAMINE REDUCED SPECIES COUNTS ################
# Use jittered species hit counts as thresholds
### Add species, starting with the least common #####
g <- g |>
  activate(nodes) |>
  left_join(species_data)

thresh <- g |> 
  activate(nodes) |>
  as_tibble() |>
  filter(hits >= 50)

# Increase the threshold so more species added on each iteration.
balance <- map(thresh$jit, ~{
  g1 <- g |> activate(nodes) |> filter(jit <= .x) 
  return(c(.x, balance_score(g1, method = "triangles")))
})

b <- as.data.frame(do.call(rbind, balance)) |>
  as_tibble() |>
  rename(threshold = 1, balance = 2) |> # jittered threshold
  mutate(name = thresh$name)

b <- b |> filter(!is.na(balance))

p1 <- b |> ggplot(aes(threshold, balance)) +
  ggtitle(paste(plot_title, "Balance with increasing species recruitment", sep = ", ")) +
  geom_line(colour = "blue") +
  geom_point()

plot(p1 + geom_text(label = b$name, vjust = 0, nudge_y = 0.001, angle = 20))


##### Triangle analysis

g <- g |>  activate(nodes) |>
 mutate(neighborhood = local_members(mindist = 1))

g1 <- g |> filter(hits > 219) # All beyond T. Pratense

# Now I want a list of all the main nodes and their neighbours.
tbl <- g1 |> activate(nodes) |> as_tibble() # tbl required for the for-loop

g3 <- g1
for(i in 1:length(g1)) {
  g2 <- g |> filter(name %in% species_data$name[tbl$neighborhood[[i]]])
  g3 <- g3 |> graph_join(g2)
}
rm(g2)

# Need to reconstruct the graph at this point, otherwise it will be directed 
# and stuff like count_signed_triangles won't work
tbl <- g3 |> activate(edges) |> as_tibble() |> select(A, B, lor, sign)

g4 <- as_tbl_graph(graph.data.frame(d = tbl, directed = FALSE))
g4 <- g4 |> activate(nodes) |> left_join(species_data)

# Add column to distinguish influential species
g4 <- g4 |> activate(nodes) |> mutate(influential = hits > 219)

g4 %>% 
  ggraph(layout = 'kk') + 
  geom_edge_link(aes(colour = factor(sign))) + 
  geom_node_point(aes(colour = influential), size = 8) +
  ggtitle('Impact species with neighbours') + 
  theme_graph()

# Edge type distribution
# tbl <- g4 |> filter(!influential) |> activate(edges) |> as_tibble()
g4 <- g4 |> activate(edges) |> 
  mutate(A = .N()$name[from]) |>
  mutate(B = .N()$name[to])
g4 <- g4 |> 
  mutate(influencer_A = .N()$influential[from]) |>
  mutate(influencer_B = .N()$influential[to]) |>
  mutate(intra_I = influencer_A & influencer_B) |>
  mutate(intra_not_I = !(influencer_A & influencer_B)) |>
  mutate(inter_I_not_I = xor(influencer_A, influencer_B))

g4_edges <- g4 |> activate(edges) |> as_tibble()

g4 %>% 
  ggraph(layout = 'kk') + 
  geom_edge_link(aes(colour = inter_I_not_I, lty = factor(sign) )) + 
  geom_node_point(aes(colour = influential), size = 8) +
  ggtitle('Impact species with neighbours') + 
  theme_graph()

# Multilevel
g4 <- g4 |> activate(nodes) |> mutate(lvl = ifelse(influential, 1, 2))

xy <- layout_as_multilevel(g4, type = "all", alpha = 25, beta = 45)

ggraph(g4, "manual", x = xy[, 1], y = xy[, 2]) +
  geom_edge_link0(
    aes(filter = (node1.lvl == 1 & node2.lvl == 1)),
    edge_colour = "firebrick3",
    alpha = 0.5,
    edge_linewidth = 0.3
  ) +
  geom_edge_link0(
    aes(filter = (node1.lvl != node2.lvl)),
    alpha = 0.3,
    edge_linewidth = 0.1,
    edge_colour = "black"
  ) +
  geom_edge_link0(
    aes(filter = (node1.lvl == 2 &
                    node2.lvl == 2)),
    edge_colour = "goldenrod3",
    edge_linewidth = 0.3,
    alpha = 0.5
  ) +
  geom_node_point(aes(shape = as.factor(lvl)), fill = "grey25", size = 3) +
  scale_shape_manual(values = c(21, 22)) +
  theme_graph() +
  coord_cartesian(clip = "off", expand = TRUE) +
  theme(legend.position = "none")


# bm <- signed_blockmodel(g4, 2, alpha = 0.5, annealing = T)
# g4 <- g4 |> activate(nodes) |> mutate(grp = bm$membership)
# 
# ggblock(g4, bm$membership,  show_blocks = TRUE, show_labels = TRUE)



