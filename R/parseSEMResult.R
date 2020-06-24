parseSEMResult <- function(resultDf = NULL) {
paths <- resultDf %>%
  select(lhs, op, rhs, Unstand_Est, Unstand_SE, STD_Genotype, STD_Genotype_SE, p_value)

#fix variable types
paths$Unstand_SE<-as.numeric(paths$Unstand_SE)
paths$STD_Genotype_SE<-as.numeric(paths$STD_Genotype_SE)
paths$p_value<-as.numeric(paths$p_value)

# Latent variables: left-hand side of "=~" paths
latent <- paths %>%
  filter(op == "=~") %>%
  select(nodes = lhs) %>%
  distinct %>%
  mutate(type="latent",
         shape = "oval",
         label = nodes,
         fontname = 'Garamond',
         fontsize = '10',
         fixedsize = 'true',
         width = '1',
         fillcolor = 'moccasin',
         rank = 'min'
  )

# Manifest variables: not latent variables
`%not_in%` <- Negate(`%in%`)
manifest <- paths %>%
  filter(op != "~1", lhs %not_in% latent$nodes) %>%
  select(nodes = lhs) %>%
  distinct %>%
  mutate(type="manifest",
         shape = "circle", #since already not original measurements
         label = nodes,
         fontname = 'Garamond',
         fontsize = '10',
         fixedsize = 'true',
         width = '0.8',
         fillcolor = 'GhostWhite',
         rank = '1'
  )

# Residuals: "~~" paths with 
residual <- paths %>%
  filter(op == "~~", lhs==rhs, lhs %in% manifest$nodes) %>%
  select(nodes = lhs) %>%
  #distinct %>%
  mutate(type="residual",
         shape = "circle",
         label = paste0("U(",nodes,")"),
         fontname = 'Garamond',
         fontsize = '8',
         fixedsize = 'true',
         width = '0.6',
         fillcolor = 'aliceblue',
         rank = 'max'
  )

# Nodes are prepared
node_set <- combine_edfs(latent, manifest, residual)

# Edges, labeled by the factor loading estimates
edges <- paths %>%
  filter(op != "~1") %>%
  mutate(
    values = round(STD_Genotype, 2),
    values.se = round(STD_Genotype_SE, 2),
    label = paste0(values," (",values.se,")")) %>%
  left_join(node_set[,c("nodes","id")],by = c("lhs" = "nodes")) %>%
  mutate(from=id) %>%
  select(-id) %>%
  left_join(node_set[,c("nodes","id")],by = c("rhs" = "nodes")) %>%
  mutate(to=id) %>%
  select(-id) %>%
  filter(from %not_in% node_set$id[which(node_set$type=="residual")]) %>%
  filter(to %not_in% node_set$id[which(node_set$type=="residual")]) %>%
  left_join(node_set[which(node_set$type=="residual"),c("nodes","id")],by = c("lhs" = "nodes")) %>%
  mutate(tofrom.residual=id, style = if_else(p_value<0.05, "solid","dashed")) %>%
  select(-id) %>%
  select(-Unstand_Est)

#might need more filters here

loadings <- edges %>%
  filter(op == "=~") %>%
  mutate(type="loading", dir="forward", rel = "leading_to", minlen="5", headport="n", weight="2")

# Regressions: "~" lines
regressions <- edges %>%
  filter(op == "~") %>%
  mutate(type="regression", dir="forward", rel = "leading_to", headport="n", weight="2")

# Covariances: ~~ for non manifest variable residuals
covariances <- edges %>%
  filter(op == "~~") %>%
  filter(is.na(tofrom.residual)) %>%
  mutate(type="covariance", rel = "related", dir="both")

# Set longer length between latent factors
covariances <- covariances %>%
  mutate(minlen = if_else( lhs %in% node_set$nodes[which(node_set$type=="latent")],"5","0.5"))

# Set variance port
covariances[which(covariances$lhs==covariances$rhs),c("headport")]<-"n"
covariances[which(covariances$lhs==covariances$rhs),c("tailport")]<-"n"


# Residual loadings: ~~ for non manifest variable residuals
residual_loadings <- edges %>%
  filter(op == "~~") %>%
  filter(!is.na(tofrom.residual)) %>%
  mutate(type="residual_loading", rel = "related", dir="forward", from=tofrom.residual, label="1", headport="s", tailport="n")

# Residual self loadings: ~~ for non manifest variable residuals
residual_self_loadings <- edges %>%
  filter(op == "~~") %>%
  filter(!is.na(tofrom.residual)) %>%
  mutate(type="residual_loading", rel = "related", dir="both", from=tofrom.residual, to=tofrom.residual, headport="s", tailport="s")


edge_set <- combine_edfs(loadings, regressions, covariances, residual_loadings, residual_self_loadings)

toReturn<-c(nodeDf=node_set, edgeDf=edge_set)
return(toReturn)

}