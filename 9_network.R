#Network treatments-seasons
reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()
#first filter by gr > than x 
#filter dark treatments

pres_abs_long <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(treatment != 'CD' &
           treatment != 'PD') %>%
  mutate(presence_absense = case_when(slope_chosen_days > 0 ~ '1',
                                      slope_chosen_days < 0  ~ '0')) %>%
  filter(presence_absense != is.na(presence_absense))

##SanketNetwork-----
library(networkD3)
##entre seasons i treatments
pres_abs_long_filt <- pres_abs_long %>%
  #mutate(season_treatment = paste(season,'_',treatment)) %>%
  mutate(season_treatment = paste(season,'',treatment)) %>%
  select(season_treatment, asv_num, presence_absense) #phylum_f, class_f, order_f, family_f

pres_abs_long_filt %>%
  head()

colnames(pres_abs_long_filt) <- c("source", "target", "value")
pres_abs_long_filt$target <- paste(pres_abs_long_filt$target, " ", sep="")
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(pres_abs_long_filt$source), as.character(pres_abs_long_filt$target)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
pres_abs_long_filt$IDsource=match(pres_abs_long_filt$source, nodes$name)-1 
pres_abs_long_filt$IDtarget=match(pres_abs_long_filt$target, nodes$name)-1

# prepare colour scale
ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

# Make the Network
sankeyNetwork(Links = pres_abs_long_filt, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20)


###provo sankey
pres_abs_long_filt <- pres_abs_long %>%
  #mutate(season_treatment = paste(season,'_',treatment)) %>%
  mutate(season_treatment = paste(season,'',treatment)) %>%
  select(season_treatment, asv_num, presence_absense) #phylum_f, class_f, order_f, family_f

pres_abs_long_filt %>%
  head()
subset_ed_f_ed2 <- subset_ed_f %>%
  mutate(seas_treat_2 = seas_treat) %>%
  select(seas_treat, seas_treat_2, asv_num)

colnames(subset_ed_f_ed2) <- c("source", "target", "value")
subset_ed_f_ed2$target <- paste(subset_ed_f_ed2$target, " ", sep="")
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(subset_ed_f_ed2$source), as.character(subset_ed_f_ed2$target)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
subset_ed_f_ed2$IDsource=match(subset_ed_f_ed2$source, nodes$name)-1 
subset_ed_f_ed2$IDtarget=match(subset_ed_f_ed2$target, nodes$name)-1

ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

sankeyNetwork(Links = subset_ed_f_ed2, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20)

# Charge the circlize library------
library(circlize)
pres_abs_long %>%
  colnames()

pres_abs_long_filt %>%
  unique()
nm <- unique(unlist(pres_abs_long_filt$season))

group = structure(gsub('\\d', '', nm), names = nm)
group

grid.col = structure(c(rep(2, 5), rep(3, 5), rep(4, 5)),
                     names = c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5)))

chordDiagram(pres_abs_long_filt, group = group, grid.col = grid.col)

pres_abs_long_filt %>%
  glimpse()

pres_abs_long_filt <- pres_abs_long %>%
  mutate(asv_num = as_factor(asv_num),
         presence_absense = as.numeric(presence_absense),
         treatment_season = as_factor(paste(treatment,'', season))) %>%
  select(treatment_season, presence_absense)

circos.clear()

circos.par(start.degree = 90, gap.degree = 1, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# Base plot
chordDiagram(
  x = pres_abs_long_filt, 
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

# Make the circular plot
data_long <- pres_abs_long %>%
  mutate(season_treatment_asv = paste(season,'_',treatment,'_',asv_num)) %>%
  select(season_treatment_asv, presence_absense)
library(magrittr)
data_long %$%
  season_treatment %>%
  unique()

pres_abs_long %>% 
  select(season, treatment, asv_num, presence_absense) %>%
  chordDiagram( transparency = 0.5) 
#, order = c("Winter _ CL", "Winter _ PL", "Winter _ DL", "Winter _ VL", 
# "Spring _ CL", "Spring _ PL", "Spring _ DL", "Spring _ VL",
# "Summer _ CL", "Summer _ PL", "Summer _ DL", "Summer _ VL", 
# "Fall _ CL", "Fall _ PL", "Fall _ DL", "Fall _ VL")

#Another try
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag) 

data_long <- pres_abs_long_filt %>%
  mutate(season_treatment = paste(season,'_',treatment)) %>%
  select(season_treatment, presence_absense) #phylum_f, class_f, order_f, family_f

data_long %>%
  glimpse()

data_long %>%
  head()

data_long %$%
  season_treatment %>%
  unique()

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 0, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor <- viridis(24, alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- mycolor[sample(1:24)]

chordDiagram(
  x = data_long, 
  grid.col = 24,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)


chordDiagram(
  x = data_long, 
  grid.col = 24,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE, 
  order = )

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
      minor.ticks = 1, 
      major.tick.percentage = 0.5,
      labels.niceFacing = FALSE)
  }
)

##computacionalment lent de crear el gràfic-----
circos.clear()
circos.info() # chordDiagram() creates two tracks, one track for labels and one track for grids with axes

#These two tracks can be controlled by annotationTrack argument. Available values for this argument are grid, 
#name and axis. The height of annotation tracks can be set by annotationTrackHeight which is the percentage 
#to the radius of unit circle and can be set by mm_h() function with an absolute unit. Axes are only added 
#if grid is set in annotationTrack
#From verion 0.4.10 of the circlize package, there is a new group argument in chordDiagram() function which 
#is very convenient for making multiple-group Chord diagrams.

#dcast converts a dataframe to a matrix (we need a matrix for a chordplot)
data_long <- pres_abs_long %>%
  mutate(season_treatment = paste(season,'_',treatment)) %>%
  select(season_treatment, presence_absense, asv_num, season, treatment) #phylum_f, class_f, order_f, family_f

data_long %>%
  head()
library(reshape2)

data_long %>%
  pivot_wider(id_cols = c('asv_num', 'season', 'treatment'), values_from = presence_absense, names_from = c('asv_num', 'season', 'treatment'))

nm = unique(unlist(dimnames(data_long)))
group = structure(gsub("\\d", "", nm), names = nm)
group


### Define ranges of circos sectors and their colors (both of the sectors and the links)-----
data_long$xmin <- 0
data_long$xmax <- rowSums(data_long) + colSums(data_long)
n <- nrow(df1)
df1$rcol<-rgb(df1$r, df1$g, df1$b, max = 255)
df1$lcol<-rgb(df1$r, df1$g, df1$b, alpha=200, max = 255)

### Plot sectors (outer part)
par(mar=rep(0,4))
circos.clear()
### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)

### Sector details
circos.initialize(factors = data_long$season, x = data_long$presence_absense) #, xlim = cbind(data_long$slope_chosen, data_long$slope_chosen)

### Plot sectors
circos.trackPlotRegion(ylim = c(0, 1), factors = df1$country, track.height=0.1,
                       #panel.fun for each sector
                       panel.fun = function(x, y) {
                         #select details of current sector
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         
                         #text direction (dd) and adjusmtents (aa)
                         theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                         dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                         aa = c(1, 0.5)
                         if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                         
                         #plot country labels
                         circos.text(x=mean(xlim), y=1.7, labels=name, facing = dd, cex=0.6,  adj = aa)
                         
                         #plot main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                     col = df1$rcol[i], border=df1$rcol[i])
                         
                         #blank in part of main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2]-rowSums(m)[i], ytop=ylim[1]+0.3, 
                                     col = "white", border = "white")
                         
                         #white line all the way around
                         circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")
                         
                         #plot axis
                         circos.axis(labels.cex=0.6, direction = "outside", major.at=seq(from=0,to=floor(df1$xmax)[i],by=5), 
                                     minor.ticks=1, labels.away.percentage = 0.15)
                       })

### Plot links (inner part)
### Add sum values to df1, marking the x-position of the first links
### out (sum1) and in (sum2). Updated for further links in loop below.
df1$sum1 <- colSums(m)
df1$sum2 <- numeric(n)

### Create a data.frame of the flow matrix sorted by flow size, to allow largest flow plotted first
df2 <- cbind(as.data.frame(m),orig=rownames(m),  stringsAsFactors=FALSE)
df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
               timevar="dest", time=rownames(m),  v.names = "m")
df2 <- arrange(df2,desc(m))

### Keep only the largest flows to avoid clutter
df2 <- subset(df2, m > quantile(m,0.6))

### Plot links
for(k in 1:nrow(df2)){
  #i,j reference of flow matrix
  i<-match(df2$orig[k],df1$country)
  j<-match(df2$dest[k],df1$country)
  
  #plot link
  circos.link(sector.index1=df1$country[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(m[i, j])),
              sector.index2=df1$country[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(m[i, j])),
              col = df1$lcol[i])
  
  #update sum1 and sum2 for use when plotting the next link
  df1$sum1[i] = df1$sum1[i] + abs(m[i, j])
  df1$sum2[j] = df1$sum2[j] + abs(m[i, j])
}

### nested proportions parallel sets plot

##another try----
Draw1() <-function {
  
  /*First disable click event on clicker button*/
    stopClicker();
  
  /*Show and run the progressBar*/
    runProgressBar(time=700*11);
  
  changeTopText(newText = "These days most people switch phones every few years. " + 
                  "Some people stay loyal, but many also switch to a different phone brand...",
                loc = 4/2, delayDisappear = 0, delayAppear = 1);
  
  changeTopText(newText = "In the next few steps I would like to introduce you to the flows of people between the phone brands ",
                loc = 8/2, delayDisappear = 9, delayAppear = 10, finalText = true);
  
  changeBottomText(newText = "Let's start by drawing out the division of the 1846 respondents, that have had at least 2 phones, among the biggest 7 brands ",
                   loc = 1/2, delayDisappear = 0, delayAppear = 10);
  
  //Remove arcs again
  d3.selectAll(".arc")
  .transition().delay(9*700).duration(2100)
  .style("opacity", 0)
  .each("end", function() {d3.selectAll(".arc").remove();});
  
}

###circular network
library(igraph)

data <- matrix(sample(0:1, 400, replace=TRUE, prob=c(0.8,0.2)), nrow=20)
network <- graph_from_adj_list(pres_abs_long_filt, mode='all', duplicate = T )

# When ploting, we can use different layouts:
par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(network, layout=layout.sphere, main="sphere")
plot(network, layout=layout.circle, main="circle")
plot(network, layout=layout.random, main="random")
plot(network, layout=layout.fruchterman.reingold, main="fruchterman.reingold")


#multiple-group chord digram


###hierarchical edge bundling------
# Libraries
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(ggraph)
library(igraph)


reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()

reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(slope_chosen_days > 0)

reg_all_slopes_chosen_silva_tax_1perc %>%
  glimpse()

connections <-reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(slope_chosen_days > 0 &
           pvalue_slope_chosen < 0.05) %>%
  select(asv_num, treatment)



d1 <- data.frame(from="origin", to=paste("group", seq(1,10), sep=""))
d2 <- data.frame(from=rep(d1$to, each=10), to=paste("subgroup", seq(1,100), sep="_"))
hierarchy <- rbind(d1, d2)

mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )

# create a vertices data.frame. One line per object of our hierarchy, giving features of nodes.
vertices <- data.frame(name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))) ) 


##test amb subset
subset <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(slope_chosen_days > 0 &
           pvalue_slope_chosen < 0.05) %>%
  dplyr::select(treatment, season, asv_num, slope_chosen, family_f) %>%
  filter(treatment == 'CD')

check<- subset %>%
  group_by(treatment, asv_num, family_f) %>%
  filter(n() >3)

subset_w <- subset %>%
  group_by(treatment, asv_num, family_f, season) %>% 
  mutate(slope_chosen = round(slope_chosen, digits = 2)) %>%
  distinct(treatment, asv_num, family_f, season, slope_chosen) %>%
  #summarize(slope_chosen = mean(slope_chosen)) %>%
  pivot_wider(id_cols = c('treatment', 'asv_num', 'family_f'), names_from = season, values_from = slope_chosen)

subset_w %>%
  create_layout()

ggraph(subset) + 
  geom_edge_link(aes(colour = factor('asv_num'))) #+ 
 # geom_node_point(data = subset)

subset_w %>%
  head()

library(tidygraph)

flaregraph <- tbl_graph(subset_w$treatment, subset_w$asv_num)
from <- match(subset_w$Winter, subset_w$Summer)
to <- match(subset_w$Summer, subset_w$Winter)
ggraph(flaregraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha = 0.1) + 
  coord_fixed()


subset_w %>%
  head()
  
  
library(ggnetwork)
  # Create a vertex attribute
subset_w %>%
  igraph::set.vertex.attribute('group', index = V(subset_w), value = Winter)


plot(subset_w, vertex.color="green",  # Changes node color
     edge.color = 'black',     # Changes edge color
     vertex.size = 10,         # Changes vertex size
     vertex.shape = 'circle',  # Changes vertex shape
     asp = 0,                  # Spread out nodes
     layout = layout_in_circle)# Format nodes in a circle


# Plot igraph object with ggplot
ggplot(subset_w, aes(x = seasa, y = Summer, xend = Spring, yend = Fall)) + 
  geom_edges() +
  geom_nodes()

subset %>%
  head()

ggplot(subset, aes(x = Winter, y = Summer, xend = Spring, yend = Fall)) + 
  geom_edges() +
  geom_nodes()

##testing circular plots-----
library(circlize)
circos.initialize(factor(subset_w$treatment), x = subset_w$Winter)
circos.track(treataments)

##igraph nº asv compartides entre treatments
subset <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(slope_chosen_days > 0 &
           pvalue_slope_chosen < 0.05) %>%
  dplyr::select(treatment, season, asv_num, slope_chosen_days, family_f) %>%
  filter(treatment %in% c('CD', 'CL', 'PD', 'PL')) %>%
  dplyr::select(-slope_chosen_days, -family_f)


plot(subset,
             layout.circle())


  vertex.connectivity(test, source = asv_num)


g <-make_ring(4)
g %>%
  class()
plot(g, layout = layout.circle, vertex.color="green")


g <- make_ring(10)
plot(g, layout=layout_with_kk, vertex.color="green")


##another option----
library("phyloseqGraphTest")
library("igraph")
library("ggnetwork")
net <- make_network(subset, max.dist=0.35)
sampledata <- data.frame(sample_data(subset))
V(subset)$id <- sampledata[names(V(subset)), "host_subject_id"]
V(subset)$litter <- sampledata[names(V(subset)), "family_relationship"]


##---
plot_net(subset, color = 'Treatment')


plot_net(prevotella, color="Species", type="taxa")


##---
subset<- subset %>%
  mutate(seas_tret = paste(treatment, season)) %>%
  select(seas_tret, asv_num)


t <- graph_from_data_frame(subset_w, vertices = subset_w$treatment)#, vertices = subset$seas_tret

plot(t, layout=layout_in_circle,
     vertex.label.dist=0.5, vertex.color="red")

##nº d'ASV per treatment i season
number_asv_treatment_season <- reg_all_slopes_chosen_silva_tax_filt %>%
  group_by(treatment, season, asv_num) %>%
  dplyr::summarize(n = n()) %>%
  mutate(treat_seas = paste(treatment, season)) %>%
  ungroup() %>%
  dplyr::select(treat_seas, asv_num)

t <- graph_from_data_frame(number_asv_treatment_season, edges = number_asv_treatment_season$asv_num)#, vertices = subset$seas_tret

plot(t, layout=layout_in_circle,
     vertex.label.dist=0.5, vertex.color="red")
  
  
##other tryyy
# Libraries
library(ggraph)
library(igraph)

subset <- subset %>%
mutate(treat_seas = paste(treatment, season))
# create a data frame giving the hierarchical structure of your individuals. 

subset %>%
  head()
# Origin on top, then groups, then subgroups
d1 <- data.frame(from=unique(subset$treat_seas), to=paste("group", seq(1,16), sep=""))
d2 <- data.frame(from=rep(unique(subset$treat_seas), each=16), to=paste("subgroup", seq(1,256), sep="_"))
hierarchy <- rbind(d1, d2)

# create a vertices data.frame. One line per object of our hierarchy, giving features of nodes.
vertices <- data.frame(name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))) ) 

# Create a graph object with the igraph library
mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )
# This is a network object, you visualize it as a network like shown in the network section!

# With igraph: 
plot(mygraph, vertex.label="", edge.arrow.size=0, vertex.size=2)

ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(alpha=0.1) +
  geom_conn_bundle(data = get_con(from = c(18,20,30), to = c(19, 50, 70)), alpha=1, width=1, colour="skyblue", tension = 1) +
  theme_void()


##tutorial
g=make_empty_graph(16)
plot(g)
Adj_m=g[]

g=make_ring_graph(16)
plot(g)

g=sample_gnm(10,40) #defining nodes and m edges placing
##the number of edges should not exceed the total possible edge number which is nC2
#for 10 nodes it is (10*9/2) = 45 edges
subset <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(slope_chosen_days > 0 &
           pvalue_slope_chosen < 0.05) %>%
  dplyr::select(treatment, season, asv_num, slope_chosen_days, family_f) %>%
  filter(treatment %in% c('CD', 'CL', 'PD', 'PL')) %>%
  dplyr::select(-family_f)

subset %>%
  colnames()

subset_ed <- subset %>%
  mutate(seas_treat = paste(treatment, season)) %>%
  select(seas_treat, asv_num, slope_chosen_days) %>%
  as_tibble()

subset %>%
  head()

# subset_f <- subset %>%
#   dplyr::select(treat_seas, asv_num)
g <- graph_from_data_frame(subset_f, vertices = subset_f$treat_seas)
plot(g)

subset_ed %>%
  glimpse()

Use the following dplyr code to identify duplicates.
subset_ed %>%
  dplyr::group_by(seas_treat, asv_num) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 

subset_ed_w <- subset_ed %>% 
  dplyr::group_by(seas_treat, asv_num) %>%
  dplyr::summarise(slope_chosen_days = mean(slope_chosen_days)) %>%
  #mutate(n = 1:nrow(subset_f)) %>%
  pivot_wider(id_cols = seas_treat, names_from = asv_num, values_from = slope_chosen_days)

subset_ed_w %$%
  seas_treat %>%
  unique()

subset_ed_w %>%
  dim()

links <- data.frame(from = c('CD Fall', 'CD Spring', 'CD Summer', 'CD Winter', 
                                'CL Fall', 'CL Spring', 'CL Summer', 'CL Winter',
                                'PD Fall', 'PD Spring', 'PD Summer', 'PD Winter',
                                'PL Fall', 'PL Spring', 'PL Summer', 'PL Winter', 
                             'CD Fall', 'CD Spring', 'CD Summer', 'CD Winter', 
                             'CL Fall', 'CL Spring', 'CL Summer', 'CL Winter',
                             'PD Fall', 'PD Spring', 'PD Summer', 'PD Winter',
                             'PL Fall', 'PL Spring', 'PL Summer', 'PL Winter'),
                       to = c('PL Winter', 'CD Fall', 'CD Spring', 'CD Summer', 
                              'CD Winter', 'CL Fall', 'CL Spring', 'CL Summer', 
                              'CL Winter','PD Fall', 'PD Spring', 'PD Summer', 
                              'PD Winter','PL Fall', 'PL Spring', 'PL Summer', 
                              'PL Summer','PL Winter', 'CD Fall', 'CD Spring', 
                              'CD Summer', 'CD Winter', 'CL Fall', 'CL Spring', 
                              'CL Summer',  'CL Winter','PD Fall', 'PD Spring', 
                              'PD Summer', 'PD Winter', 'PL Fall', 'PL Spring', 
                              'PL Spring','PL Summer', 'PL Summer','PL Winter', 
                              'CD Fall', 'CD Spring', 'CD Summer', 'CD Winter', 
                              'CL Fall', 'CL Spring', 'CL Summer', 'CL Winter',
                              'PD Fall', 'PD Spring', 'PD Summer', 'PD Winter', 
                              'PL Fall', 'PL Spring', 'PD Winter', 'PL Fall'),
                       weight = subset_ed_w)
links$weight
nodes <- data.frame(subset_ed_w$seas_treat)
)

subset_ed_w %>%
  dim()
network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
g<-tbl_graph(nodes, links, directed = FALSE)

g <- graph_from_data_frame(links, vertices = subset_ed_w)

plot(graph, layout = layout.circle, label.vertex = TRUE)

##another try
subset_ed_w_f <- 
  subset_ed_w %>%
  filter(seas_treat %in% c('CD Fall', 'PL Winter'))


links <- data.frame(from = c('CD Fall', 'PL Winter' ),
                    to = c('PL Winter', 'CD Fall'),
                    weight = subset_ed_w_f[,2:137])
links$weight
nodes <- data.frame(subset_ed_w_f[,2:137])
)

subset_ed_w %>%
  dim()
network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
g<-tbl_graph(nodes, links, directed = FALSE)

g <- graph_from_data_frame(links, vertices = subset_ed_w_f)

plot(g, layout = layout.circle, label.vertex = TRUE)

#ggraph
g %>% ggraph(circular = TRUE) +
  geom_node_point() +
  geom_edge_link()+#data = subset_ed_w
  geom_node_point()+ 
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph()

E(g)$weigth <-  seq(ecount(g))
strength(g)
strength(g, mode="out")

## The typical case is that these tables are read in from files....-------
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
print(g, e=TRUE, v=TRUE)
plot(graph)



##other-----

ggraph(subset_ed_w, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_link(size=0.4, alpha=0.1) +
  #geom_conn_bundle(data = get_con(from = from_head, to = to_head), alpha = 1, colour="#69b3a2", width=2, tension=0) + 
  geom_node_text(aes(x = x*1.01, y=y*1.01, filter = leaf, label=shortName, angle = angle, hjust=hjust), size=1.5, alpha=1) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))


##other----
subset_ed_f <- subset_ed %>%
  select(-slope_chosen_days)
chordDiagram(
  x = subset_ed_f, 
  grid.col = palette_class_assigned,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

    
##other
    flaregraph <- tbl_graph(common_asv_sum)
    from <- match(common_asv_sum$from, common_asv_sum$relative_connectivity)
    to <- match(common_asv_sum$to, common_asv_sum$relative_connectivity)
    ggraph(flaregraph, layout = 'dendrogram', circular = TRUE) + 
      geom_conn_bundle(data = get_con(from = from, to = to), alpha = 0.1) + 
      coord_fixed()

    subset_ed_w
    flaregraph <- tbl_graph(nodes = common_asv_sum$sample, edges = common_asv_sum$relative_connectivity)
    from <- match(flare$imports$from, flare$vertices$name)
    to <- match(flare$imports$to, flare$vertices$name)
    ggraph(flaregraph, layout = 'dendrogram', circular = TRUE) + 
      geom_conn_bundle(data = get_con(from = from, to = to), alpha = 0.1) + 
      coord_fixed()  
     
    ##other
    subset_ed_f_f  <-  subset_ed_f %>%
      group_by(asv_num) %>%
      dplyr::summarize(n = n()) %>%
      filter (n>1) %>%
      left_join(subset_ed_f) %>%
      select(seas_treat, asv_num)
    
    subset_ed_f_rev <- subset_ed_f %>%
      select(asv_num, seas_treat)
    # Create graph of highschool friendships
    graph <- as_tbl_graph(subset_ed_f_f) %>% 
      mutate(Popularity = centrality_degree())# %>%
      filter(Popularity > 0 )
    ##degree sentralty is the simples centraality measaure to compute. Recaall that a node's degree is simply a count
    ##of how many social connections (i.e., edges) it has. The degree centraaliry for a node is simply its degree. A node with
    ## 10 social connections would have a degree centrality of 10. A node with 1 edge would have a degree centrality of 1.
   
    # plot using ggraph
    ggraph(graph) + #,, layout = 'kk' 
      geom_edge_fan(aes(alpha = after_stat(index)), show.legend = TRUE) + 
      geom_node_point(aes()) + #size = Popularity
      #facet_edges(~) + 
      theme_graph(foreground = 'steelblue', fg_text_colour = 'white')
 
    
    graph <- as_tbl_graph(common_asv_sum) %>%
      mutate(Popularity = centrality_degree())# %>%
    filter(Popularity > 0 )
    
    ##ultim intent abans del cap de setmana-----
    
    
      connect <- subset_ed_f %>% 
      gather(key="to", value="value", -1) %>%
      mutate(to = gsub("\\.", " ",to)) %>%
      na.omit() 
    
    # Number of connection per person
    c( as.character(subset_ed_f$from), as.character(subset_ed_f$to)) %>%
      as.tibble() %>%
      group_by(value) %>%
      summarize(n=n()) -> coauth
    colnames(coauth) <- c("name", "n")
    
    # Create a graph object with igraph
    mygraph <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
    
    # Find community
    com <- walktrap.community(mygraph)
    
    #Reorder dataset and make the graph
    coauth <- coauth %>% 
      mutate( grp = com$membership) %>%
      arrange(grp) %>%
      mutate(name=factor(name, name))
    
    # keep only 10 first communities
    coauth <- coauth %>% 
      filter(grp<16)
    
    # keep only this people in edges
    connect <- connect %>%
      filter(from %in% coauth$name) %>%
      filter(to %in% coauth$name)
    
    # Add label angle
    number_of_bar=nrow(coauth)
    coauth$id = seq(1, nrow(coauth))
    angle= 360 * (coauth$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    coauth$hjust <- ifelse(angle > 90 & angle<270, 1, 0)
    coauth$angle <- ifelse(angle > 90 & angle<270, angle+180, angle)
    
    # Create a graph object with igraph
    mygraph <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
    
    # prepare a vector of n color in the viridis scale
    mycolor <- colormap(colormap=colormaps$viridis, nshades=max(coauth$grp))
    mycolor <- sample(mycolor, length(mycolor))
    
    # Make the graph
    ggraph(mygraph, layout="circle") + 
      geom_edge_link(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=FALSE) +
      geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.9) +
      scale_size_continuous(range=c(0.5,8)) +
      scale_color_manual(values=mycolor) +
      geom_node_text(aes(label=paste("    ",name,"    "), angle=angle, hjust=hjust), size=2.3, color="black") +
      theme_void() +
      theme(
        legend.position="none",
        plot.margin=unit(c(0,0,0,0), "null"),
        panel.spacing=unit(c(0,0,0,0), "null")
      ) +
      expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) 
    
    
    subset_ed_w_f %>%
      ungroup() %>%
      dplyr::select(-seas_treat) %>%
      mutate(everything(case_when(is.numeric() ~ 1,
                           !is.numeric() ~ 0)))
 
    
####Calculating how many ASV were common between treatments and seasons-------
    #trec els duplicats de WINTER PD (comprovar que no estigui malament a tot el dataset)
    subset_ed %$%
      seas_treat %>%
      unique()
    
     subset_ed2 <- subset_ed %>%
      group_by(seas_treat, asv_num) %>%
      distinct(seas_treat, asv_num)
    
    num_asv_seas_treat <- subset_ed2 %>%
      group_by(seas_treat) %>%
      dplyr::summarize(n = n()) %>%
      as_tibble()

##count number of common ASVs
    library(glue)
    common.asvs <- function(data, sample1, sample2){
      
      { 
        data1 <- data %>%
        subset(seas_treat == sample1)
      data2 <- data %>%
        subset(seas_treat == sample2)
      }
      {
        sample12 <- paste0(sample1, '_', sample2) %>%
          as.character() 
        result <-  data1 %>%
          left_join(data2, by = 'asv_num') %>%
          filter(seas_treat.y != is.na(seas_treat.y)) %>%
          group_by(seas_treat.x) %>%
          dplyr::summarise(num_common_asvs = n()) %>%
          dplyr::mutate(seas_treat.x = str_replace(seas_treat.x, !!{sample1}, !!{sample12}))
        }
    return(result)
    }
    
    subset_ed2 %>%
      colnames()

  ##calculating nº of common ASVs between treatments and seasons
CD_CL_W <- common.asvs(data = subset_ed2, sample1 = 'CD Winter', sample2 = 'CL Winter')      
CD_CL_Sp <- common.asvs(data = subset_ed2, sample1 = 'CD Spring', sample2 = 'CL Spring')      
CD_CL_Su <- common.asvs(data = subset_ed2, sample1 = 'CD Summer', sample2 = 'CL Summer')  
CD_CL_F <- common.asvs(data = subset_ed2, sample1 = 'CD Fall', sample2 = 'CL Fall') 
    
PD_PL_W <- common.asvs(data = subset_ed2, sample1 = 'PD Winter', sample2 = 'PL Winter')      
PD_PL_Sp <- common.asvs(data = subset_ed2, sample1 = 'PD Spring', sample2 = 'PL Spring')      
PD_PL_Su <- common.asvs(data = subset_ed2, sample1 = 'PD Summer', sample2 = 'PL Summer')  
PD_PL_F <- common.asvs(data = subset_ed2, sample1 = 'PD Fall', sample2 = 'PL Fall') 

CD_DL_W <- common.asvs(data = subset_ed2, sample1 = 'CD Winter', sample2 = 'DL Winter')      
CD_DL_Sp <- common.asvs(data = subset_ed2, sample1 = 'CD Spring', sample2 = 'DL Spring')      
CD_DL_Su <- common.asvs(data = subset_ed2, sample1 = 'CD Summer', sample2 = 'DL Summer')  
CD_DL_F <- common.asvs(data = subset_ed2, sample1 = 'CD Fall', sample2 = 'DL Fall') 

CD_CD_Sp_F <- common.asvs(data = subset_ed2, sample1 = 'CD Spring', sample2 = 'CD Fall')  
CD_CD_Sp_Su <- common.asvs(data = subset_ed2, sample1 = 'CD Spring', sample2 = 'CD Summer') 

common_asv_sum <- rbind( CD_CL_W, CD_CL_Sp, CD_CL_Su, CD_CL_F,
  PD_PL_W, PD_PL_Sp, PD_PL_Su, PD_PL_F,
  CD_DL_W, CD_DL_Sp, CD_DL_Su, CD_DL_F, CD_CD_Sp_F, CD_CD_Sp_Su) %>%
  separate(seas_treat.x, '_', into = c('from', 'to')) %>%
  left_join(num_asv_seas_treat, by = c('from' = 'seas_treat')) %>%
  left_join(num_asv_seas_treat, by = c('to' = 'seas_treat')) %>%
  mutate(division = case_when( n.x < n.y ~ n.x,
                               n.y < n.x ~ n.y)) %>%
  mutate(relative_connectivity = num_common_asvs/division) %>%
  mutate(sample = paste0(from, to)) %>%
  select(from, to, relative_connectivity) 
  

common_asv_sum 

## The typical case is that these tables are read in from files....-------
# actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
#                             "Esmeralda"),
#                      age=c(48,33,45,34,21),
#                      gender=c("F","M","F","M","F"))
# relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
#                                "David", "Esmeralda"),
#                         to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
#                         same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
#                         friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))


# vertices <- data.frame(rbind(as_tibble(common_asv_sum$from), 
#                              as_tibble(common_asv_sum$to))) 
# 
# g <- graph_from_data_frame(common_asv_sum, vertices = vertices)
# 
# print(g, e=TRUE, v=TRUE, layout = circular)
# plot(g, e=TRUE, v=TRUE, layout = circular)


# library(tidyr)
# 
# sample1 <- 'PD Winter'
# sample2 <-  'PL Winter'
# sample12 <- paste(sample1, sample2) %>%
#   as.character()
# test <- subset_ed2 %>%
# mutate(seas_treat = str_replace(seas_treat, !!{sample1}, !!{sample12}))
# 
# str(sample12)
# 
# data1 <- subset_ed2 %>%
#   subset(seas_treat == 'PD Winter')
# 
# data2 <- subset_ed2 %>%
#   subset(seas_treat == 'PL Winter')
# 
# 
#   result <-  data1 %>%
#     left_join(data2) %>%
#     filter(asv_num != is.na(asv_num)) %>%
#     group_by(seas_treat) %>%
#     dplyr::summarize(n = n()) %>%
#     mutate(seas_treat = str_replace(seas_treat, !!{sample1}, !!{sample12})) 


###----
# The flare dataset is provided in ggraph
edges <- common_asv_sum %>%
  select(from, to)
vertices <- common_asv_sum %>%
  select(sample) %>%
  mutate(size = 10) %>%
  as_tibble()

vertices <- data.frame(rbind(as_tibble(common_asv_sum$from), 
                             as_tibble(common_asv_sum$to))) %>%
  group_by(value) %>%
  unique() %>%
  mutate(size = 10) %>%
  as_tibble()

connections <- common_asv_sum %>%
  select(from, to, relative_connectivity)

# # Preparation to draw labels properly:
# vertices$id=NA
# myleaves=which(is.na( match(vertices$name, edges$from) ))
# nleaves=length(myleaves)
# vertices$id[ myleaves ] = seq(1:nleaves)
# vertices$angle= 90 - 360 * vertices$id / nleaves
# vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)
# vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

# Build a network object from this dataset:
mygraph <- graph_from_data_frame(edges, vertices = vertices)

# The connection object must refer to the ids of the leaves:
from = match( connections$from, vertices$value)
to = match( connections$to, vertices$value)

# Basic dendrogram
# ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
#   geom_conn_bundle(data = get_con(from = from, to = to, paths = connections$relative_connectivity), alpha = 0.1, colour="#69b3a2") +
#   #geom_edge_link(size=0.4, alpha=0.1) +
#   #geom_node_text(aes(x = x*1.01, y=y*1.01, filter = leaf, label=vertices$value, angle = angle, hjust=hjust), size=1.5, alpha=1) +
#   coord_fixed() +
#   theme_void() +
#   theme(
#     legend.position="none",
#     plot.margin=unit(c(0,0,0,0),"cm"),
#   ) +
#   expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

##CIRCULAR CONNECTIONS 
ggraph(mygraph, layout="linear", circular = TRUE) + 
  geom_edge_arc(aes(width = connections$relative_connectivity), edge_colour="black", edge_alpha=0.2, fold=TRUE) +
  #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  #geom_ar_link(aes(width = connections$relative_connectivity))
  #scale_color_manual(values=mycolor) +
  coord_fixed()+
  geom_node_text(aes(label=name)) + #, angle=65, hjust=1, nudge_y = -1.1, size=2.3
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0), "null"),
    panel.spacing=unit(c(0,0,3.4,0), "null")) ##expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)

#####AMB TOTES LES DADES------
##importing functions
source("src/calculate_common_asv_conditions.R")
source("src/calculate_exclusive_asv_condition.R")

data_for_common_asvs <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0 &
           pvalue_slope_chosen < 0.05) %>%
  dplyr::select(treatment, season, asv_num) %>% #, slope_chosen_days, family_f
  #filter(treatment %in% c('CD', 'CL', 'PD', 'PL')) %>%
  #dplyr::select(-slope_chosen_days, -family_f) %>%
  mutate(seas_treat = paste(treatment, season)) %>%
  group_by(seas_treat, asv_num) %>% ##sembla que PD Winter está repetit al dataset original (trobar a on es duplica, té efecte aals gràfics?)
  distinct(seas_treat, asv_num) 
  
num_asv_seas_treat <- data_for_common_asvs %>%
  group_by(seas_treat) %>%
  dplyr::summarize(n = n()) %>%
  as_tibble()

##calculating nº of common ASVs between treatments and seasons------
##entre estacions mateix treatment
CD_CD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CD Spring')      
CD_CD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CD Summer')      
CD_CD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CD Fall')  
CD_CD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CD Summer') 
CD_CD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CD Fall') 
CD_CD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'CD Fall') 

CL_CL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CL Spring')      
CL_CL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CL Summer')      
CL_CL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CL Fall')  
CL_CL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CL Summer') 
CL_CL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CL Fall') 
CL_CL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'CL Fall') 

PL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PL Spring')      
PL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PL Summer')      
PL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PL Fall')  
PL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PL Summer') 
PL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PL Fall') 
PL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'PL Fall') 

PD_PD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PD Spring')      
PD_PD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PD Summer')      
PD_PD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PD Fall')  
PD_PD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PD Summer') 
PD_PD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PD Fall') 
PD_PD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'PD Fall')

DL_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'DL Spring')      
DL_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'DL Summer')      
DL_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'DL Fall')  
DL_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'DL Summer') 
DL_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'DL Fall') 
DL_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'DL Fall')

VL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'VL Spring')      
VL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'VL Summer')      
VL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'VL Fall')  
VL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'VL Summer') 
VL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'VL Fall') 
VL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'VL Fall')

##entre treatments
###CD_CL
CD_CL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Winter')      
CD_CL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CL Spring')      
CD_CL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'CL Summer')  
CD_CL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Fall', sample2 = 'CL Fall') 
CD_CL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Spring') 
CD_CL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CD Spring') 
CD_CL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Summer') 
CD_CL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CD Summer') 
CD_CL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Fall') 
CD_CL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CD Fall') 
CD_CL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CL Summer') 
CD_CL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CD Summer') 
CD_CL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CL Fall') 
CD_CL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CD Fall') 
CD_CL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'CL Fall') 
CD_CL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'CD Fall') 

##CD_PL
CD_PD_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Winter')      
CD_PD_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Spring')      
CD_PD_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'CL Summer')  
CD_PD_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'CL Fall') 
CD_PD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Spring') 
CD_PD_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Spring') 
CD_PD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Summer') 
CD_PD_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Summer') 
CD_PD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Fall') 
CD_PD_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Fall') 
CD_PD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Summer') 
CD_PD_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Summer') 
CD_PD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Fall') 
CD_PD_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Fall') 
CD_PD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'CL Fall') 
CD_PD_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PD Fall')

##CD_PL
CD_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Winter')      
CD_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Spring')      
CD_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Summer')  
CD_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Fall', sample2 = 'CL Fall') 
CD_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Spring') 
CD_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Spring') 
CD_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Summer') 
CD_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Summer') 
CD_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Fall') 
CD_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Fall') 
CD_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Summer') 
CD_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Summer') 
CD_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Fall') 
CD_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Fall') 
CD_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Fall') 
CD_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Fall') 

##CD_DL
CD_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Winter')      
CD_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'DL Spring')      
CD_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'DL Summer')  
CD_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Fall', sample2 = 'DL Fall') 
CD_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Spring') 
CD_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CD Spring') 
CD_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Summer') 
CD_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CD Summer') 
CD_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Fall') 
CD_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CD Fall') 
CD_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'DL Summer') 
CD_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CD Summer') 
CD_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'DL Fall') 
CD_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CD Fall') 
CD_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'DL Fall') 
CD_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'CD Fall') 

##CD_VL (16)
CD_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Winter')      
CD_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'VL Spring')      
CD_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'VL Summer')  
CD_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Fall', sample2 = 'VL Fall') 
CD_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Spring') 
CD_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CD Spring') 
CD_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Summer') 
CD_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CD Summer') 
CD_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Fall') 
CD_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CD Fall') 
CD_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'VL Summer') 
CD_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CD Summer') 
CD_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'VL Fall') 
CD_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CD Fall') 
CD_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'VL Fall') 
CD_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'CD Fall') 


##CL_PL (16)
CL_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Winter')      
CL_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Spring')      
CL_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Summer')  
CL_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'PL Fall') 
CL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Spring') 
CL_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Spring') 
CL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Summer') 
CL_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Summer') 
CL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Fall') 
CL_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Fall') 
CL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Summer') 
CL_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Summer') 
CL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Fall') 
CL_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Fall') 
CL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Fall') 
CL_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Fall') 

##CL_PD (16)
CL_PD_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Winter')      
CL_PD_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Spring')      
CL_PD_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PD Summer')  
CL_PD_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'PD Fall') 
CL_PD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Spring') 
CL_PD_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Spring') 
CL_PD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Summer') 
CL_PD_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Summer') 
CL_PD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Fall') 
CL_PD_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Fall') 
CL_PD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Summer') 
CL_PD_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Summer') 
CL_PD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Fall') 
CL_PD_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Fall') 
CL_PD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PD Fall') 
CL_PD_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'CL Fall') 

##CL_PL (16)
CL_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Winter')      
CL_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Spring')      
CL_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Summer')  
CL_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'PL Fall') 
CL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Spring') 
CL_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Spring') 
CL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Summer') 
CL_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Summer') 
CL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Fall') 
CL_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Fall') 
CL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Summer') 
CL_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Summer') 
CL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Fall') 
CL_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Fall') 
CL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Fall') 
CL_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Fall') 

##CL_DL (16)
CL_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Winter')      
CL_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Spring')      
CL_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Summer')  
CL_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'DL Fall') 
CL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Spring') 
CL_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Spring') 
CL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Summer') 
CL_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Summer') 
CL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Fall') 
CL_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Fall') 
CL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Summer') 
CL_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Summer') 
CL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Fall') 
CL_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Fall') 
CL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Fall') 
CL_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'CL Fall') 

##CL_DL (16)
CL_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Winter')      
CL_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Spring')      
CL_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Summer')  
CL_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'DL Fall') 
CL_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Spring') 
CL_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Spring') 
CL_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Summer') 
CL_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Summer') 
CL_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Fall') 
CL_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Fall') 
CL_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Summer') 
CL_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Summer') 
CL_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Fall') 
CL_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Fall') 
CL_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Fall') 
CL_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'CL Fall') 

##CL_VL (16)
CL_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Winter')      
CL_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'VL Spring')      
CL_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'VL Summer')  
CL_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'VL Fall') 
CL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Spring') 
CL_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CL Spring') 
CL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Summer') 
CL_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CL Summer') 
CL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Fall') 
CL_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CL Fall') 
CL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'VL Summer') 
CL_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CL Summer') 
CL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'VL Fall') 
CL_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CL Fall') 
CL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'VL Fall') 
CL_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'CL Fall')

##PD_PL (16)
PD_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Winter')      
PD_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PL Spring')      
PD_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'PL Summer')  
PD_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'PL Fall') 
PD_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Spring') 
PD_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PD Spring') 
PD_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Summer') 
PD_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PD Summer') 
PD_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Fall') 
PD_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PD Fall') 
PD_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PL Summer') 
PD_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PD Summer') 
PD_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PL Fall') 
PD_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PD Fall') 
PD_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'PL Fall') 
PD_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'PD Fall') 

##PD_DL (16)
PD_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Winter')      
PD_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'DL Spring')      
PD_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'DL Summer')  
PD_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'DL Fall') 
PD_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Spring') 
PD_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PD Spring') 
PD_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Summer') 
PD_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PD Summer') 
PD_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Fall') 
PD_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PD Fall') 
PD_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'DL Summer') 
PD_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PD Summer') 
PD_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'DL Fall') 
PD_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PD Fall') 
PD_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'DL Fall') 
PD_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'PD Fall') 

##PD_VL (16)
PD_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Winter')      
PD_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'VL Spring')      
PD_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'VL Summer')  
PD_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'VL Fall') 
PD_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Spring') 
PD_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PD Spring') 
PD_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Summer') 
PD_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PD Summer') 
PD_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Fall') 
PD_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PD Fall') 
PD_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'VL Summer') 
PD_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PD Summer') 
PD_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'VL Fall') 
PD_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PD Fall') 
PD_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'VL Fall') 
PD_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'PD Fall') 

##PL_DL (16)
PL_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Winter')      
PL_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'DL Spring')      
PL_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'DL Summer')  
PL_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Fall', sample2 = 'DL Fall') 
PL_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Spring') 
PL_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PL Spring') 
PL_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Summer') 
PL_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PL Summer') 
PL_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Fall') 
PL_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PL Fall') 
PL_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'DL Summer') 
PL_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PL Summer') 
PL_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'DL Fall') 
PL_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PL Fall') 
PL_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'DL Fall') 
PL_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'PL Fall') 

##PL_VL (16)
PL_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Winter')      
PL_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'VL Spring')      
PL_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'VL Summer')  
PL_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Fall', sample2 = 'VL Fall') 
PL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Spring') 
PL_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PL Spring') 
PL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Summer') 
PL_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PL Summer') 
PL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Fall') 
PL_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PL Fall') 
PL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'VL Summer') 
PL_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PL Summer') 
PL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'VL Fall') 
PL_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PL Fall') 
PL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'VL Fall') 
PL_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'PL Fall') 

##DL_VL
DL_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Winter')      
DL_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'VL Spring')      
DL_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'VL Summer')  
DL_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Fall', sample2 = 'VL Fall') 
DL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Spring') 
DL_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'DL Spring') 
DL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Summer') 
DL_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'DL Summer') 
DL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Fall') 
DL_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'DL Fall') 
DL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'VL Summer') 
DL_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'DL Summer') 
DL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'VL Fall') 
DL_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'DL Fall') 
DL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'VL Fall') 
DL_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'DL Fall') 

common_asv_sum <- 
  bind_rows(CD_CD_W_Sp, CD_CD_W_Su, CD_CD_W_F, CD_CD_Sp_Su, CD_CD_Sp_F, CD_CD_Su_F,
         CL_CL_W_Sp, CL_CL_W_Su, CL_CL_W_F, CL_CL_Sp_Su, CL_CL_Sp_F, CL_CL_Su_F,
         PL_PL_W_Sp, PL_PL_W_Su, PL_PL_W_F, PL_PL_Sp_Su, PL_PL_Sp_F, PL_PL_Su_F,
         PD_PD_W_Sp, PD_PD_W_Su, PD_PD_W_F, PD_PD_Sp_Su, PD_PD_Sp_F, PD_PD_Su_F,
         DL_DL_W_Sp, DL_DL_W_Su, DL_DL_W_F, DL_DL_Sp_Su, DL_DL_Sp_F, DL_DL_Su_F,
         VL_VL_W_Sp, VL_VL_W_Su, VL_VL_W_F, VL_VL_Sp_Su, VL_VL_Sp_F, VL_VL_Su_F,
         CD_CL_W, CD_CL_Sp, CD_CL_Su, CD_CL_F, 
         CD_CL_W_Sp, CD_CL_Sp_W, CD_CL_W_Su, CD_CL_Su_W, CD_CL_W_F, CD_CL_F_W, CD_CL_Sp_Su , CD_CL_Su_Sp, CD_CL_Sp_F, CD_CL_F_Sp, CD_CL_Su_F, CD_CL_F_Su, 
         CD_PD_W, CD_PD_Sp, CD_PD_Su, CD_PD_F,
         CD_PD_W_Sp, CD_PD_Sp_W, CD_PD_W_Su, CD_PD_Su_W, CD_PD_W_F, CD_PD_F_W, CD_PD_Sp_Su , CD_PD_Su_Sp, CD_PD_Sp_F, CD_PD_F_Sp, CD_PD_Su_F, CD_PD_F_Su,
         CD_PL_W, CD_PL_Sp, CD_PL_Su, CD_PL_F,
         CD_PL_W_Sp, CD_PL_Sp_W, CD_PL_W_Su, CD_PL_Su_W, CD_PL_W_F, CD_PL_F_W, CD_PL_Sp_Su , CD_PL_Su_Sp, CD_PL_Sp_F, CD_PL_F_Sp, CD_PL_Su_F, CD_PL_F_Su,
         CD_DL_W, CD_DL_Sp, CD_DL_Su, CD_DL_F,
         CD_DL_W_Sp, CD_DL_Sp_W, CD_DL_W_Su, CD_DL_Su_W, CD_DL_W_F, CD_DL_F_W, CD_DL_Sp_Su , CD_DL_Su_Sp, CD_DL_Sp_F, CD_DL_F_Sp, CD_DL_Su_F, CD_DL_F_Su, 
         CD_VL_Sp, CD_VL_Su, CD_VL_F,
         CD_VL_W_Sp, CD_VL_Sp_W, CD_VL_W_Su, CD_VL_Su_W, CD_VL_W_F, CD_VL_F_W, CD_VL_Sp_Su , CD_VL_Su_Sp, CD_VL_Sp_F, CD_VL_F_Sp, CD_VL_Su_F, CD_VL_F_Su, 
         PD_PL_W, PD_PL_Sp, PD_PL_Su, PD_PL_F,
         PD_PL_W_Sp, PD_PL_Sp_W, PD_PL_W_Su, PD_PL_Su_W, PD_PL_W_F, PD_PL_F_W, PD_PL_Sp_Su , PD_PL_Su_Sp, PD_PL_Sp_F, PD_PL_F_Sp, PD_PL_Su_F, PD_PL_F_Su,
         PD_DL_W, PD_DL_Sp, PD_DL_Su, PD_DL_F,
         PD_DL_W_Sp, PD_DL_Sp_W, PD_DL_W_Su, PD_DL_Su_W, PD_DL_W_F, PD_DL_F_W, PD_DL_Sp_Su , PD_DL_Su_Sp, PD_DL_Sp_F, PD_DL_F_Sp, PD_DL_Su_F, PD_DL_F_Su, 
         PD_VL_Sp, PD_VL_Su, PD_VL_F,
         PD_VL_W_Sp, PD_VL_Sp_W, PD_VL_W_Su, PD_VL_Su_W, PD_VL_W_F, PD_VL_F_W, PD_VL_Sp_Su , PD_VL_Su_Sp, PD_VL_Sp_F, PD_VL_F_Sp, PD_VL_Su_F, PD_VL_F_Su,                 
         PL_DL_W, PL_DL_Sp, PL_DL_Su, PL_DL_F,
         PL_DL_W_Sp, PL_DL_Sp_W, PL_DL_W_Su, PL_DL_Su_W, PL_DL_W_F, PL_DL_F_W, PL_DL_Sp_Su , PL_DL_Su_Sp, PL_DL_Sp_F, PL_DL_F_Sp, PL_DL_Su_F, PL_DL_F_Su, 
         PL_VL_Sp, PL_VL_Su, PL_VL_F,
         PL_VL_W_Sp, PL_VL_Sp_W, PL_VL_W_Su, PL_VL_Su_W, PL_VL_W_F, PL_VL_F_W, PL_VL_Sp_Su , PL_VL_Su_Sp, PL_VL_Sp_F, PL_VL_F_Sp, PL_VL_Su_F, PL_VL_F_Su,                 
         DL_VL_Sp, DL_VL_Su, DL_VL_F,
         DL_VL_W_Sp, DL_VL_Sp_W, DL_VL_W_Su, DL_VL_Su_W, DL_VL_W_F, DL_VL_F_W, DL_VL_Sp_Su , DL_VL_Su_Sp, DL_VL_Sp_F, DL_VL_F_Sp, DL_VL_Su_F, DL_VL_F_Su) %>%
  separate(seas_treat.x, '_', into = c('from', 'to')) %>%
  left_join(num_asv_seas_treat, by = c('from' = 'seas_treat')) %>%
  left_join(num_asv_seas_treat, by = c('to' = 'seas_treat')) %>%
  mutate(division = case_when( n.x < n.y ~ n.x,
                               n.y < n.x ~ n.y,
                               n.x == n.y ~ n.x)) %>%
  mutate(relative_connectivity = num_common_asvs/division) %>%
  mutate(sample = paste0(from, to))  #%>%
 # select(from, to, relative_connectivity) 

common_asv_sum 

##prepare graph
edges <- common_asv_sum %>%
  select(from, to) %>%
  filter(from != is.na(from))
vertices <- common_asv_sum %>%
  select(sample) %>%
  mutate(size = 10) %>%
  as_tibble()

vertices <- data.frame(rbind(as_tibble(common_asv_sum$from), 
                             as_tibble(common_asv_sum$to))) %>%
  group_by(value) %>%
  unique() %>%
  mutate(size = 10) %>%
  separate(value, ' ', into = c('Treatment', 'Season'), remove = F) %>%
  mutate(Treatment = as.factor(Treatment),
         Season = as.factor(Season),
         value = as.factor(value)) %>%
  select(value, Season, Treatment) %>%
  filter(value != is.na(value)) %>%
  as_tibble()
  
vertices$Season <- vertices$Season %>% 
  factor(levels = c('Winter', 'Spring', 'Summer', 'Fall'))

vertices$Treatment <- vertices$Treatment %>% 
  factor(levels = c('CD', 'CL', 'PD', 'PL', 'DL', 'VL'))

vertices <- vertices %>%
  mutate(value = fct_relevel(value, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                      "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                      "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                      "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
  
  arrange(value) %>%
  as_tibble()

connections <- common_asv_sum %>%
  separate(from, ' ', into = c('Treatment', 'Season'), remove = F) %>%
  filter(relative_connectivity != is.na(relative_connectivity)) %>%
  select(from, to, relative_connectivity, Treatment, Season) 
  

num_connections1 <- common_asv_sum %>%
  group_by(from) %>%
  dplyr::summarize(n = n()) %>%
  right_join(vertices, by = c('from' = 'value'))


num_connections2 <- common_asv_sum %>%
  group_by(to) %>%
  dplyr::summarize(n = n()) %>%
  right_join(vertices, by = c('to' = 'value'))
  
num_connections12 <- num_connections1 %>%
  left_join(num_connections2, by=c('Season' = 'Season', 'Treatment' = 'Treatment')) %>%
  mutate(n.x = case_when(#!is.na(n.x) ~ n.x,
    is.na(n.x) ~ '0',
    TRUE ~ as.character(n.x)),
    n.y = case_when(#!is.na(n.x) ~ n.x,
      is.na(n.y) ~ '0',
      TRUE ~ as.character(n.y))
    ) %>%
  mutate(num_connections = (as.numeric(n.x)+as.numeric(n.y)))

common_asv_sum$from

common_asv_sum$to

# # Preparation to draw labels properly:
# vertices$id=NA
# myleaves=which(is.na( match(vertices$name, edges$from) ))
# nleaves=length(myleaves)
# vertices$id[ myleaves ] = seq(1:nleaves)
# vertices$angle= 90 - 360 * vertices$id / nleaves
# vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)
# vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

# Build a network object from this dataset:
mygraph <- graph_from_data_frame(edges, vertices = vertices)
#graph_test <-  tbl_graph(edges, vertices)
##reorder nodes
# s <- arrange(names(V(mygraph)), vertices$Season)
# mygraph <- permute(mygraph, match(V(mygraph)$name, s))

# V(mygraph)$node_label <- names(V(mygraph)) 

# The connection object must refer to the ids of the leaves:
#from = match( connections$Treatment$from, vertices$Treatment$value) (no funciona)
from = match( connections$from, vertices$value)
to = match( connections$to, vertices$value)

# Basic dendrogram
# ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
#   geom_conn_bundle(data = get_con(from = from, to = to, paths = connections$relative_connectivity), alpha = 0.1, colour="#69b3a2") +
#   #geom_edge_link(size=0.4, alpha=0.1) +
#   #geom_node_text(aes(x = x*1.01, y=y*1.01, filter = leaf, label=vertices$value, angle = angle, hjust=hjust), size=1.5, alpha=1) +
#   coord_fixed() +
#   theme_void() +
#   theme(
#     legend.position="none",
#     plot.margin=unit(c(0,0,0,0),"cm"),
#   ) +
#   expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

##CIRCULAR CONNECTIONS 
specific_range <- common_asv_sum %$%
  relative_connectivity %>%
  range()

palete_gradient <- c(#"#005a47",
  #"#678078",
  #"#aaaaaa",
  "#e6e6e6",
           #"#7b7b7b",
           "#003318") 
                              

common_asvs <- ggraph(mygraph, layout="linear", circular = TRUE) + 
  # geom_conn_bundle(data = get_con(from = from, to = to, alpha = 0.1, colour="#e6e6e6") #alpha=0.2, , tension = 0.9,
  # )+
  geom_edge_arc(aes(width = as.numeric(connections$relative_connectivity), 
                    color = as.numeric(connections$relative_connectivity)
  ),
  fold=FALSE,
  edge_alpha=0.8,
  lineend = 'square',
  linejoin = 'mitre', check_overlap = TRUE, #strength = 0.5
  # title = 'Relative connectivity'
  # edge_width = 1
  ) + #, position = position_jitter(0.5)
  #geom_edge_density(aes(fill = as.numeric(connections$relative_connectivity)))+
  # geom_edge_fan0(aes(#alpha = after_stat(connections$relative_connectivity,
  #                                      color = connections$Treatment))+
  #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
  # scale_color_manual(values = palette_treatments_remei)+
  scale_edge_colour_gradientn(colours = palete_gradient, name = 'Relative connectivity')+
  scale_edge_width_continuous(limits = c(0.3, 0.89), guide = 'none')+
  #geom_edge_diagonal(aes(width = connections$relative_connectivity, color = connections$relative_connectivity))+
  #geom_ar_link(aes(width = connections$relative_connectivity))
  #scale_color_manual(values=mycolor) +
  coord_fixed()+
  geom_node_point(aes(color = Season, size = num_connections12$num_connections), alpha = 0.9)+#, size = 0.75
  scale_size(range=c(14,19), limits = c(14,19)) +
  geom_node_text(aes(label=Treatment), vjust = -0.5, show.legend = TRUE) + #, angle=65, hjust=1, nudge_y = -1.1, size=2.3
  scale_color_manual(values = palette_seasons_4)+
  guides(color = guide_legend(ncol = 1, size = 10,
                              override.aes = aes(label = '')),
         shape = guide_legend(ncol = 1, size = 1),
         # edge = guide_legend(title = 'Relative connectivity',
         #                           override.aes = aes(label = ''),
         #                     name = 'Relative connectivity'),
         size = 'none',
         #width = guide_legend(title = 'Relative connectivity')
  )+
  theme_void() +
  theme(legend.text = element_text(size = 5), 
        legend.title = element_text(size =7))

ggsave(filename = 'common_asvs_all.pdf', plot = common_asvs, 
       path = 'results/figures/corrected_asv_num/',
       width = 188, height = 188, units = 'mm')

###same graph but only with highest relations----
edges <- common_asv_sum %>%
  filter(relative_connectivity > 0.3) %>%
  select(from, to)
vertices <- common_asv_sum %>%
  filter(relative_connectivity > 0.3) %>%
  select(sample) %>%
  mutate(size = 10) %>%
  as_tibble()


vertices <- data.frame(rbind(as_tibble(common_asv_sum$from), 
                             as_tibble(common_asv_sum$to))) %>%
  group_by(value) %>%
  unique() %>%
  mutate(size = 10) %>%
  separate(value, ' ', into = c('Treatment', 'Season'), remove = F) %>%
  mutate(Season = as.factor(Season),
    Treatment = as.factor(Treatment),
     value = as.factor(value)) %>%
  as_tibble() %>%
  dplyr::select(value, Season, Treatment)

vertices$Season <- vertices$Season %>% 
  factor(levels = c('Winter', 'Spring', 'Summer', 'Fall'))


vertices$Treatment <- vertices$Treatment %>% 
  factor(levels = c('CD', 'CL', 'PD', 'PL', 'DL', 'VL'))

connections <- common_asv_sum %>%
  select(from, to, relative_connectivity) %>%
  filter(relative_connectivity > 0.3)

vertices <- vertices %>%
  mutate(value = fct_relevel(value, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                       "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                       "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                       "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%

  arrange(value) %>%
  mutate(name = factor(Season)) %>%
  as_tibble()

  vertices$value
# # Preparation to draw labels properly:
# vertices$id=NA
# myleaves=which(is.na( match(vertices$name, edges$from) ))
# nleaves=length(myleaves)
# vertices$id[ myleaves ] = seq(1:nleaves)
# vertices$angle= 90 - 360 * vertices$id / nleaves
# vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)
# vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

# Build a network object from this dataset:
mygraph_filt <- graph_from_data_frame(edges, vertices = vertices)

##reorder nodes
# s <- arrange(names(V(mygraph)), vertices$Season)
# mygraph <- permute(mygraph, match(V(mygraph)$name, s))

# V(mygraph)$node_label <- names(V(mygraph)) 

# The connection object must refer to the ids of the leaves:
from = match( connections$from, vertices$name)
to = match( connections$to, vertices$name)

# Basic dendrogram
# ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
#   geom_conn_bundle(data = get_con(from = from, to = to, paths = connections$relative_connectivity), alpha = 0.1, colour="#69b3a2") +
#   #geom_edge_link(size=0.4, alpha=0.1) +
#   #geom_node_text(aes(x = x*1.01, y=y*1.01, filter = leaf, label=vertices$value, angle = angle, hjust=hjust), size=1.5, alpha=1) +
#   coord_fixed() +
#   theme_void() +
#   theme(
#     legend.position="none",
#     plot.margin=unit(c(0,0,0,0),"cm"),
#   ) +
#   expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

##CIRCULAR CONNECTIONS 
specific_range <- common_asv_sum %>%
filter(relative_connectivity > 0.3) %$%
  relative_connectivity %>%
  range()

palete_gradient <- c(#"#005a47",
  #"#678078",
  #"#aaaaaa",
  "#e6e6e6",
           "#7b7b7b",
           "#003318") 
           
common_asvs_filt <- 
  ggraph(mygraph_filt, layout="linear", circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to, alpha = 0.1, colour="#e6e6e6") #alpha=0.2, , tension = 0.9,
                   )+
  geom_edge_arc(aes(width = as.numeric(connections$relative_connectivity), 
                    color = as.numeric(connections$relative_connectivity)
                    ),
                fold=FALSE,
                edge_alpha=0.8,
                lineend = 'square',
                linejoin = 'mitre', check_overlap = TRUE, #strength = 0.5
               # title = 'Relative connectivity'
                # edge_width = 1
  ) + #, position = position_jitter(0.5)
  #geom_edge_density(aes(fill = as.numeric(connections$relative_connectivity)))+
  # geom_edge_fan0(aes(#alpha = after_stat(connections$relative_connectivity,
  #                                      color = connections$Treatment))+
  #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
  # scale_color_manual(values = palette_treatments_remei)+
  scale_edge_colour_gradientn(colours = palete_gradient, name = 'Relative connectivity')+
  scale_edge_width_continuous(limits = c(0.3, 0.89), guide = 'none')+
  #geom_edge_diagonal(aes(width = connections$relative_connectivity, color = connections$relative_connectivity))+
  #geom_ar_link(aes(width = connections$relative_connectivity))
  #scale_color_manual(values=mycolor) +
  coord_fixed()+
  geom_node_point(aes(color = Season, size = num_connections12$num_connections), alpha = 0.9)+#, size = 0.75
  scale_size(range=c(14,19), limits = c(14,19)) +
  geom_node_text(aes(label=Treatment), vjust = -0.5, show.legend = TRUE) + #, angle=65, hjust=1, nudge_y = -1.1, size=2.3
  scale_color_manual(values = palette_seasons_4)+
    guides(color = guide_legend(ncol = 1, size = 10,
                              override.aes = aes(label = '')),
           shape = guide_legend(ncol = 1, size = 1),
           # edge = guide_legend(title = 'Relative connectivity',
           #                           override.aes = aes(label = ''),
           #                     name = 'Relative connectivity'),
        size = 'none',
  #width = guide_legend(title = 'Relative connectivity')
  )+
    theme_void() +
    theme(legend.text = element_text(size = 5), 
          legend.title = element_text(size =7))
   
  
# theme(
   #   legend.position="right")
  #   plot.margin=unit(c(0,0,0,0), "null"),
  #   panel.spacing=unit(c(0,0,0,0), "null"))


  ggsave(filename = 'common_asvs_filt_0.3.pdf', plot = common_asvs_filt, 
       path = 'results/figures/corrected_asv_num/',
       width = 188, height = 188, units = 'mm')


###same graph but only with GR>1-------------
  data_for_common_asvs <- reg_all_slopes_chosen_silva_tax %>%
    filter(slope_chosen_days > 1 &
             pvalue_slope_chosen < 0.05) %>%
    dplyr::select(treatment, season, asv_num) %>% #, slope_chosen_days, family_f
    #filter(treatment %in% c('CD', 'CL', 'PD', 'PL')) %>%
    #dplyr::select(-slope_chosen_days, -family_f) %>%
    mutate(seas_treat = paste(treatment, season)) %>%
    group_by(seas_treat, asv_num) %>% ##sembla que PD Winter está repetit al dataset original (trobar a on es duplica, té efecte aals gràfics?)
    distinct(seas_treat, asv_num) 
  
  num_asv_seas_treat <- data_for_common_asvs %>%
    group_by(seas_treat) %>%
    dplyr::summarize(n = n()) %>%
    as_tibble()
  
  ##calculating nº of common ASVs between treatments and seasons------
  ##entre estacions mateix treatment
  CD_CD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CD Spring')      
  CD_CD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CD Summer')      
  CD_CD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CD Fall')  
  CD_CD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CD Summer') 
  CD_CD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CD Fall') 
  CD_CD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'CD Fall') 
  
  CL_CL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CL Spring')      
  CL_CL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CL Summer')      
  CL_CL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CL Fall')  
  CL_CL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CL Summer') 
  CL_CL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CL Fall') 
  CL_CL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'CL Fall') 
  
  PL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PL Spring')      
  PL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PL Summer')      
  PL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PL Fall')  
  PL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PL Summer') 
  PL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PL Fall') 
  PL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'PL Fall') 
  
  PD_PD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PD Spring')      
  PD_PD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PD Summer')      
  PD_PD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PD Fall')  
  PD_PD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PD Summer') 
  PD_PD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PD Fall') 
  PD_PD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'PD Fall')
  
  DL_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'DL Spring')      
  DL_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'DL Summer')      
  DL_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'DL Fall')  
  DL_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'DL Summer') 
  DL_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'DL Fall') 
  DL_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'DL Fall')
  
  VL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'VL Spring')      
  VL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'VL Summer')      
  VL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'VL Fall')  
  VL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'VL Summer') 
  VL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'VL Fall') 
  VL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'VL Fall')
  
  ##entre treatments
  ###CD_CL
  CD_CL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Winter')      
  CD_CL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CL Spring')      
  CD_CL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'CL Summer')  
  CD_CL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Fall', sample2 = 'CL Fall') 
  CD_CL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Spring') 
  CD_CL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CD Spring') 
  CD_CL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Summer') 
  CD_CL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CD Summer') 
  CD_CL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Fall') 
  CD_CL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CD Fall') 
  CD_CL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CL Summer') 
  CD_CL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CD Summer') 
  CD_CL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CL Fall') 
  CD_CL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CD Fall') 
  CD_CL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'CL Fall') 
  CD_CL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'CD Fall') 
  
  ##CD_PL
  CD_PD_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Winter')      
  CD_PD_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Spring')      
  CD_PD_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'CL Summer')  
  CD_PD_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'CL Fall') 
  CD_PD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Spring') 
  CD_PD_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Spring') 
  CD_PD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Summer') 
  CD_PD_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Summer') 
  CD_PD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Fall') 
  CD_PD_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Fall') 
  CD_PD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Summer') 
  CD_PD_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Summer') 
  CD_PD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Fall') 
  CD_PD_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Fall') 
  CD_PD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'CL Fall') 
  CD_PD_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PD Fall')
  
  ##CD_PL
  CD_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Winter')      
  CD_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Spring')      
  CD_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Summer')  
  CD_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Fall', sample2 = 'CL Fall') 
  CD_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Spring') 
  CD_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Spring') 
  CD_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Summer') 
  CD_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Summer') 
  CD_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Fall') 
  CD_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Fall') 
  CD_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Summer') 
  CD_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Summer') 
  CD_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Fall') 
  CD_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Fall') 
  CD_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Fall') 
  CD_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Fall') 
  
  ##CD_DL
  CD_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Winter')      
  CD_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'DL Spring')      
  CD_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'DL Summer')  
  CD_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Fall', sample2 = 'DL Fall') 
  CD_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Spring') 
  CD_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CD Spring') 
  CD_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Summer') 
  CD_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CD Summer') 
  CD_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Fall') 
  CD_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CD Fall') 
  CD_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'DL Summer') 
  CD_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CD Summer') 
  CD_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'DL Fall') 
  CD_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CD Fall') 
  CD_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'DL Fall') 
  CD_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'CD Fall') 
  
  ##CD_VL (16)
  CD_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Winter')      
  CD_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'VL Spring')      
  CD_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'VL Summer')  
  CD_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Fall', sample2 = 'VL Fall') 
  CD_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Spring') 
  CD_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CD Spring') 
  CD_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Summer') 
  CD_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CD Summer') 
  CD_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Fall') 
  CD_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CD Fall') 
  CD_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'VL Summer') 
  CD_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CD Summer') 
  CD_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'VL Fall') 
  CD_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CD Fall') 
  CD_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'VL Fall') 
  CD_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'CD Fall') 
  
  
  ##CL_PL (16)
  CL_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Winter')      
  CL_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Spring')      
  CL_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Summer')  
  CL_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'PL Fall') 
  CL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Spring') 
  CL_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Spring') 
  CL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Summer') 
  CL_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Summer') 
  CL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Fall') 
  CL_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Fall') 
  CL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Summer') 
  CL_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Summer') 
  CL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Fall') 
  CL_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Fall') 
  CL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Fall') 
  CL_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Fall') 
  
  ##CL_PD (16)
  CL_PD_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Winter')      
  CL_PD_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Spring')      
  CL_PD_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PD Summer')  
  CL_PD_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'PD Fall') 
  CL_PD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Spring') 
  CL_PD_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Spring') 
  CL_PD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Summer') 
  CL_PD_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Summer') 
  CL_PD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Fall') 
  CL_PD_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Fall') 
  CL_PD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Summer') 
  CL_PD_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Summer') 
  CL_PD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Fall') 
  CL_PD_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Fall') 
  CL_PD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PD Fall') 
  CL_PD_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'CL Fall') 
  
  ##CL_PL (16)
  CL_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Winter')      
  CL_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Spring')      
  CL_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Summer')  
  CL_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'PL Fall') 
  CL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Spring') 
  CL_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Spring') 
  CL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Summer') 
  CL_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Summer') 
  CL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Fall') 
  CL_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Fall') 
  CL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Summer') 
  CL_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Summer') 
  CL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Fall') 
  CL_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Fall') 
  CL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Fall') 
  CL_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Fall') 
  
  ##CL_DL (16)
  CL_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Winter')      
  CL_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Spring')      
  CL_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Summer')  
  CL_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'DL Fall') 
  CL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Spring') 
  CL_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Spring') 
  CL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Summer') 
  CL_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Summer') 
  CL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Fall') 
  CL_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Fall') 
  CL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Summer') 
  CL_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Summer') 
  CL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Fall') 
  CL_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Fall') 
  CL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Fall') 
  CL_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'CL Fall') 
  
  ##CL_DL (16)
  CL_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Winter')      
  CL_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Spring')      
  CL_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Summer')  
  CL_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'DL Fall') 
  CL_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Spring') 
  CL_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Spring') 
  CL_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Summer') 
  CL_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Summer') 
  CL_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Fall') 
  CL_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Fall') 
  CL_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Summer') 
  CL_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Summer') 
  CL_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Fall') 
  CL_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Fall') 
  CL_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Fall') 
  CL_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'CL Fall') 
  
  ##CL_VL (16)
  CL_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Winter')      
  CL_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'VL Spring')      
  CL_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'VL Summer')  
  CL_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'VL Fall') 
  CL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Spring') 
  CL_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CL Spring') 
  CL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Summer') 
  CL_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CL Summer') 
  CL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Fall') 
  CL_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CL Fall') 
  CL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'VL Summer') 
  CL_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CL Summer') 
  CL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'VL Fall') 
  CL_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CL Fall') 
  CL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'VL Fall') 
  CL_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'CL Fall')
  
  ##PD_PL (16)
  PD_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Winter')      
  PD_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PL Spring')      
  PD_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'PL Summer')  
  PD_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'PL Fall') 
  PD_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Spring') 
  PD_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PD Spring') 
  PD_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Summer') 
  PD_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PD Summer') 
  PD_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Fall') 
  PD_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PD Fall') 
  PD_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PL Summer') 
  PD_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PD Summer') 
  PD_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PL Fall') 
  PD_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PD Fall') 
  PD_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'PL Fall') 
  PD_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'PD Fall') 
  
  ##PD_DL (16)
  PD_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Winter')      
  PD_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'DL Spring')      
  PD_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'DL Summer')  
  PD_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'DL Fall') 
  PD_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Spring') 
  PD_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PD Spring') 
  PD_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Summer') 
  PD_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PD Summer') 
  PD_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Fall') 
  PD_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PD Fall') 
  PD_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'DL Summer') 
  PD_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PD Summer') 
  PD_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'DL Fall') 
  PD_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PD Fall') 
  PD_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'DL Fall') 
  PD_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'PD Fall') 
  
  ##PD_VL (16)
  PD_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Winter')      
  PD_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'VL Spring')      
  PD_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'VL Summer')  
  PD_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'VL Fall') 
  PD_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Spring') 
  PD_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PD Spring') 
  PD_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Summer') 
  PD_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PD Summer') 
  PD_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Fall') 
  PD_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PD Fall') 
  PD_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'VL Summer') 
  PD_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PD Summer') 
  PD_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'VL Fall') 
  PD_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PD Fall') 
  PD_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'VL Fall') 
  PD_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'PD Fall') 
  
  ##PL_DL (16)
  PL_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Winter')      
  PL_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'DL Spring')      
  PL_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'DL Summer')  
  PL_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Fall', sample2 = 'DL Fall') 
  PL_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Spring') 
  PL_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PL Spring') 
  PL_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Summer') 
  PL_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PL Summer') 
  PL_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Fall') 
  PL_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PL Fall') 
  PL_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'DL Summer') 
  PL_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PL Summer') 
  PL_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'DL Fall') 
  PL_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PL Fall') 
  PL_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'DL Fall') 
  PL_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'PL Fall') 
  
  ##PL_VL (16)
  PL_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Winter')      
  PL_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'VL Spring')      
  PL_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'VL Summer')  
  PL_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Fall', sample2 = 'VL Fall') 
  PL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Spring') 
  PL_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PL Spring') 
  PL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Summer') 
  PL_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PL Summer') 
  PL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Fall') 
  PL_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PL Fall') 
  PL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'VL Summer') 
  PL_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PL Summer') 
  PL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'VL Fall') 
  PL_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PL Fall') 
  PL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'VL Fall') 
  PL_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'PL Fall') 
  
  ##DL_VL
  DL_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Winter')      
  DL_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'VL Spring')      
  DL_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'VL Summer')  
  DL_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Fall', sample2 = 'VL Fall') 
  DL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Spring') 
  DL_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'DL Spring') 
  DL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Summer') 
  DL_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'DL Summer') 
  DL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Fall') 
  DL_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'DL Fall') 
  DL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'VL Summer') 
  DL_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'DL Summer') 
  DL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'VL Fall') 
  DL_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'DL Fall') 
  DL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'VL Fall') 
  DL_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'DL Fall') 
  
  common_asv_sum <- 
    bind_rows(CD_CD_W_Sp, CD_CD_W_Su, CD_CD_W_F, CD_CD_Sp_Su, CD_CD_Sp_F, CD_CD_Su_F,
              CL_CL_W_Sp, CL_CL_W_Su, CL_CL_W_F, CL_CL_Sp_Su, CL_CL_Sp_F, CL_CL_Su_F,
              PL_PL_W_Sp, PL_PL_W_Su, PL_PL_W_F, PL_PL_Sp_Su, PL_PL_Sp_F, PL_PL_Su_F,
              PD_PD_W_Sp, PD_PD_W_Su, PD_PD_W_F, PD_PD_Sp_Su, PD_PD_Sp_F, PD_PD_Su_F,
              DL_DL_W_Sp, DL_DL_W_Su, DL_DL_W_F, DL_DL_Sp_Su, DL_DL_Sp_F, DL_DL_Su_F,
              VL_VL_W_Sp, VL_VL_W_Su, VL_VL_W_F, VL_VL_Sp_Su, VL_VL_Sp_F, VL_VL_Su_F,
              CD_CL_W, CD_CL_Sp, CD_CL_Su, CD_CL_F, 
              CD_CL_W_Sp, CD_CL_Sp_W, CD_CL_W_Su, CD_CL_Su_W, CD_CL_W_F, CD_CL_F_W, CD_CL_Sp_Su , CD_CL_Su_Sp, CD_CL_Sp_F, CD_CL_F_Sp, CD_CL_Su_F, CD_CL_F_Su, 
              CD_PD_W, CD_PD_Sp, CD_PD_Su, CD_PD_F,
              CD_PD_W_Sp, CD_PD_Sp_W, CD_PD_W_Su, CD_PD_Su_W, CD_PD_W_F, CD_PD_F_W, CD_PD_Sp_Su , CD_PD_Su_Sp, CD_PD_Sp_F, CD_PD_F_Sp, CD_PD_Su_F, CD_PD_F_Su,
              CD_PL_W, CD_PL_Sp, CD_PL_Su, CD_PL_F,
              CD_PL_W_Sp, CD_PL_Sp_W, CD_PL_W_Su, CD_PL_Su_W, CD_PL_W_F, CD_PL_F_W, CD_PL_Sp_Su , CD_PL_Su_Sp, CD_PL_Sp_F, CD_PL_F_Sp, CD_PL_Su_F, CD_PL_F_Su,
              CD_DL_W, CD_DL_Sp, CD_DL_Su, CD_DL_F,
              CD_DL_W_Sp, CD_DL_Sp_W, CD_DL_W_Su, CD_DL_Su_W, CD_DL_W_F, CD_DL_F_W, CD_DL_Sp_Su , CD_DL_Su_Sp, CD_DL_Sp_F, CD_DL_F_Sp, CD_DL_Su_F, CD_DL_F_Su, 
              CD_VL_Sp, CD_VL_Su, CD_VL_F,
              CD_VL_W_Sp, CD_VL_Sp_W, CD_VL_W_Su, CD_VL_Su_W, CD_VL_W_F, CD_VL_F_W, CD_VL_Sp_Su , CD_VL_Su_Sp, CD_VL_Sp_F, CD_VL_F_Sp, CD_VL_Su_F, CD_VL_F_Su, 
              PD_PL_W, PD_PL_Sp, PD_PL_Su, PD_PL_F,
              PD_PL_W_Sp, PD_PL_Sp_W, PD_PL_W_Su, PD_PL_Su_W, PD_PL_W_F, PD_PL_F_W, PD_PL_Sp_Su , PD_PL_Su_Sp, PD_PL_Sp_F, PD_PL_F_Sp, PD_PL_Su_F, PD_PL_F_Su,
              PD_DL_W, PD_DL_Sp, PD_DL_Su, PD_DL_F,
              PD_DL_W_Sp, PD_DL_Sp_W, PD_DL_W_Su, PD_DL_Su_W, PD_DL_W_F, PD_DL_F_W, PD_DL_Sp_Su , PD_DL_Su_Sp, PD_DL_Sp_F, PD_DL_F_Sp, PD_DL_Su_F, PD_DL_F_Su, 
              PD_VL_Sp, PD_VL_Su, PD_VL_F,
              PD_VL_W_Sp, PD_VL_Sp_W, PD_VL_W_Su, PD_VL_Su_W, PD_VL_W_F, PD_VL_F_W, PD_VL_Sp_Su , PD_VL_Su_Sp, PD_VL_Sp_F, PD_VL_F_Sp, PD_VL_Su_F, PD_VL_F_Su,                 
              PL_DL_W, PL_DL_Sp, PL_DL_Su, PL_DL_F,
              PL_DL_W_Sp, PL_DL_Sp_W, PL_DL_W_Su, PL_DL_Su_W, PL_DL_W_F, PL_DL_F_W, PL_DL_Sp_Su , PL_DL_Su_Sp, PL_DL_Sp_F, PL_DL_F_Sp, PL_DL_Su_F, PL_DL_F_Su, 
              PL_VL_Sp, PL_VL_Su, PL_VL_F,
              PL_VL_W_Sp, PL_VL_Sp_W, PL_VL_W_Su, PL_VL_Su_W, PL_VL_W_F, PL_VL_F_W, PL_VL_Sp_Su , PL_VL_Su_Sp, PL_VL_Sp_F, PL_VL_F_Sp, PL_VL_Su_F, PL_VL_F_Su,                 
              DL_VL_Sp, DL_VL_Su, DL_VL_F,
              DL_VL_W_Sp, DL_VL_Sp_W, DL_VL_W_Su, DL_VL_Su_W, DL_VL_W_F, DL_VL_F_W, DL_VL_Sp_Su , DL_VL_Su_Sp, DL_VL_Sp_F, DL_VL_F_Sp, DL_VL_Su_F, DL_VL_F_Su) %>%
    separate(seas_treat.x, '_', into = c('from', 'to')) %>%
    left_join(num_asv_seas_treat, by = c('from' = 'seas_treat')) %>%
    left_join(num_asv_seas_treat, by = c('to' = 'seas_treat')) %>%
    mutate(division = case_when( n.x < n.y ~ n.x,
                                 n.y < n.x ~ n.y,
                                 n.x == n.y ~ n.x)) %>%
    mutate(relative_connectivity = num_common_asvs/division) %>%
    mutate(sample = paste0(from, to))  #%>%
  # select(from, to, relative_connectivity) 
  
  common_asv_sum 
  
  ##prepare graph
  edges <- common_asv_sum %>%
    select(from, to) %>%
    filter(from != is.na(from))
  vertices <- common_asv_sum %>%
    select(sample) %>%
    mutate(size = 10) %>%
    as_tibble()
  
  vertices <- data.frame(rbind(as_tibble(common_asv_sum$from), 
                               as_tibble(common_asv_sum$to))) %>%
    group_by(value) %>%
    unique() %>%
    mutate(size = 10) %>%
    separate(value, ' ', into = c('Treatment', 'Season'), remove = F) %>%
    mutate(Treatment = as.factor(Treatment),
           Season = as.factor(Season),
           value = as.factor(value)) %>%
    select(value, Season, Treatment) %>%
    filter(value != is.na(value)) %>%
    as_tibble()
  
  vertices$Season <- vertices$Season %>% 
    factor(levels = c('Winter', 'Spring', 'Summer', 'Fall'))
  
  vertices$Treatment <- vertices$Treatment %>% 
    factor(levels = c('CD', 'CL', 'PD', 'PL', 'DL', 'VL'))
  
  vertices <- vertices %>%
    mutate(value = fct_relevel(value, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                        "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                        "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                        "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
    
    arrange(value) %>%
    as_tibble()
  
  connections <- common_asv_sum %>%
    separate(from, ' ', into = c('Treatment', 'Season'), remove = F) %>%
    filter(relative_connectivity != is.na(relative_connectivity)) %>%
    select(from, to, relative_connectivity, Treatment, Season) 
  
  
  num_connections1 <- common_asv_sum %>%
    group_by(from) %>%
    dplyr::summarize(n = n()) %>%
    right_join(vertices, by = c('from' = 'value'))
  
  
  num_connections2 <- common_asv_sum %>%
    group_by(to) %>%
    dplyr::summarize(n = n()) %>%
    right_join(vertices, by = c('to' = 'value'))
  
  num_connections12 <- num_connections1 %>%
    left_join(num_connections2, by=c('Season' = 'Season', 'Treatment' = 'Treatment')) %>%
    mutate(n.x = case_when(#!is.na(n.x) ~ n.x,
      is.na(n.x) ~ '0',
      TRUE ~ as.character(n.x)),
      n.y = case_when(#!is.na(n.x) ~ n.x,
        is.na(n.y) ~ '0',
        TRUE ~ as.character(n.y))
    ) %>%
    mutate(num_connections = (as.numeric(n.x)+as.numeric(n.y)))
  
  common_asv_sum$from
  common_asv_sum$to
  

  # Build a network object from this dataset:
  mygraph_rel <- graph_from_data_frame(edges, vertices = vertices)
  
  # The connection object must refer to the ids of the leaves:
  from = match( connections$from, vertices$value)
  to = match( connections$to, vertices$value)
  
  ##CIRCULAR CONNECTIONS 
  specific_range <- common_asv_sum %$%
    relative_connectivity %>%
    range()
  
  palete_gradient <- c(#"#005a47",
    #"#678078",
    #"#aaaaaa",
    "#e6e6e6",
             #"#7b7b7b",
             "#003318") 
             
##dibuixem només les connexions de més del 30%  
  common_asvs_1gr_relative <- ggraph(mygraph_rel, layout="linear", circular = TRUE) + 
    # geom_conn_bundle(data = get_con(from = from, to = to, alpha = 0.1, colour="#e6e6e6") #alpha=0.2, , tension = 0.9,
    # )+
    geom_edge_arc(aes(width = as.numeric(connections$relative_connectivity), 
                      color = as.numeric(connections$relative_connectivity)
    ),
    fold=FALSE,
    edge_alpha=0.8,
    lineend = 'square',
    linejoin = 'mitre', check_overlap = TRUE, #strength = 0.5
    # title = 'Relative connectivity'
    # edge_width = 1
    ) + #, position = position_jitter(0.5)
    #geom_edge_density(aes(fill = as.numeric(connections$relative_connectivity)))+
    # geom_edge_fan0(aes(#alpha = after_stat(connections$relative_connectivity,
    #                                      color = connections$Treatment))+
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    # scale_color_manual(values = palette_treatments_remei)+
    scale_edge_colour_gradientn(colours = palete_gradient, name = 'Relative\nconnectivity')+
    scale_edge_width_continuous(limits = c(0.3, 0.89), guide = 'none')+
    #geom_edge_diagonal(aes(width = connections$relative_connectivity, color = connections$relative_connectivity))+
    #geom_ar_link(aes(width = connections$relative_connectivity))
    #scale_color_manual(values=mycolor) +
    coord_fixed()+
    geom_node_point(aes(color = Season, size = num_connections12$num_connections/2), alpha = 0.8)+#, size = 0.75
    scale_size(range=c(11/2,19/2), limits = c(11/2,19/2)) +
    geom_node_text(aes(label=Treatment), vjust = -0.1, show.legend = TRUE) + #, angle=65, hjust=1, nudge_y = -1.1, size=2.3
    scale_color_manual(values = palette_seasons_4)+
    guides(color = guide_legend(ncol = 1, size = 10,
                                override.aes = aes(label = '')),
           shape = guide_legend(ncol = 1, size = 1),
           # edge = guide_legend(title = 'Relative connectivity',
           #                           override.aes = aes(label = ''),
           #                     name = 'Relative connectivity'),
           size = 'none',
           #width = guide_legend(title = 'Relative connectivity')
    )+
    theme_void() +
    theme(legend.text = element_text(size = 5), 
          legend.title = element_text(size =7),
          legend.box.spacing = unit(1, 'mm'))
  
  ggsave(filename = 'common_asvs_all_1gr_relative_02.pdf', plot = common_asvs_1gr_relative, 
         path = 'results/figures/corrected_asv_num/',
         width = 188, height = 188, units = 'mm')
  
  ##common ASVs absolut-----
  ##prepare graph
  edges <- common_asv_sum %>%
    select(from, to) %>%
    filter(from != is.na(from))
  vertices <- common_asv_sum %>%
    select(sample) %>%
    mutate(size = 10) %>%
    as_tibble()
  
  vertices <- data.frame(rbind(as_tibble(common_asv_sum$from), 
                               as_tibble(common_asv_sum$to))) %>%
    group_by(value) %>%
    unique() %>%
    mutate(size = 10) %>%
    separate(value, ' ', into = c('Treatment', 'Season'), remove = F) %>%
    mutate(Treatment = as.factor(Treatment),
           Season = as.factor(Season),
           value = as.factor(value)) %>%
    select(value, Season, Treatment) %>%
    filter(value != is.na(value)) %>%
    as_tibble()
  
  vertices$Season <- vertices$Season %>% 
    factor(levels = c('Winter', 'Spring', 'Summer', 'Fall'))
  
  vertices$Treatment <- vertices$Treatment %>% 
    factor(levels = c('CD', 'CL', 'PD', 'PL', 'DL', 'VL'))
  
  vertices <- vertices %>%
    mutate(value = fct_relevel(value, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                        "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                        "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                        "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
    
    arrange(value) %>%
    as_tibble()
  
  connections <- common_asv_sum %>%
    separate(from, ' ', into = c('Treatment', 'Season'), remove = F) %>%
    filter(num_common_asvs != is.na(num_common_asvs)) %>%
    select(from, to, num_common_asvs, Treatment, Season) 
  
  
  num_connections1 <- common_asv_sum %>%
    group_by(from) %>%
    dplyr::summarize(n = n()) %>%
    right_join(vertices, by = c('from' = 'value'))
  
  
  num_connections2 <- common_asv_sum %>%
    group_by(to) %>%
    dplyr::summarize(n = n()) %>%
    right_join(vertices, by = c('to' = 'value'))
  
  num_connections12 <- num_connections1 %>%
    left_join(num_connections2, by=c('Season' = 'Season', 'Treatment' = 'Treatment')) %>%
    mutate(n.x = case_when(#!is.na(n.x) ~ n.x,
      is.na(n.x) ~ '0',
      TRUE ~ as.character(n.x)),
      n.y = case_when(#!is.na(n.x) ~ n.x,
        is.na(n.y) ~ '0',
        TRUE ~ as.character(n.y))
    ) %>%
    mutate(num_connections = (as.numeric(n.x)+as.numeric(n.y)))
  
  common_asv_sum$from
  common_asv_sum$to
  
  # Build a network object from this dataset:
  mygraph <- graph_from_data_frame(edges, vertices = vertices)
  
  # The connection object must refer to the ids of the leaves:
  from = match( connections$from, vertices$value)
  to = match( connections$to, vertices$value)
  
##CIRCULAR CONNECTIONS 
  specific_range <- common_asv_sum %$%
  num_common_asvs %>%
    range() ##el 30 percent de 132 és 39,6 dibuixo a partir d'aquí per poder comparar amb el de relative
  ##i el 20 percent es 26,4
  
  palete_gradient <- c(#"#005a47",
    #"#678078",
    #"#aaaaaa",
    "#e6e6e6",
             #"#7b7b7b",
             "#003318") 
             
  
  common_asvs_1gr_absolut <- ggraph(mygraph, layout="linear", circular = TRUE) + 
    # geom_conn_bundle(data = get_con(from = from, to = to, alpha = 0.1, colour="#e6e6e6") #alpha=0.2, , tension = 0.9,
    # )+
    geom_edge_arc(aes(width = as.numeric(connections$num_common_asvs), 
                      color = as.numeric(connections$num_common_asvs)
    ),
    fold=FALSE,
    edge_alpha=0.8,
    lineend = 'square',
    linejoin = 'mitre', check_overlap = TRUE, #strength = 0.5
    # title = 'Relative connectivity'
    # edge_width = 1
    ) + #, position = position_jitter(0.5)
    #geom_edge_density(aes(fill = as.numeric(connections$relative_connectivity)))+
    # geom_edge_fan0(aes(#alpha = after_stat(connections$relative_connectivity,
    #                                      color = connections$Treatment))+
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    # scale_color_manual(values = palette_treatments_remei)+
    scale_edge_colour_gradientn(colours = palete_gradient, name = 'Nº of\ncommon\nASVs')+
    scale_edge_width_continuous(limits = c(39.6, 132), guide = 'none')+
    #geom_edge_diagonal(aes(width = connections$relative_connectivity, color = connections$relative_connectivity))+
    #geom_ar_link(aes(width = connections$relative_connectivity))
    #scale_color_manual(values=mycolor) +
    coord_fixed()+
    geom_node_point(aes(color = Season, size = num_connections12$num_connections/2), alpha = 0.8)+#, size = 0.75
    scale_size(range=c(11/2,19/2), limits = c(11/2,19/2)) +
    geom_node_text(aes(label=Treatment), vjust = 0.1, show.legend = TRUE) + #, angle=65, hjust=1, nudge_y = -1.1, size=2.3
    scale_color_manual(values = palette_seasons_4)+
    guides(color = guide_legend(ncol = 1, size = 10,
                                override.aes = aes(label = '')),
           shape = guide_legend(ncol = 1, size = 1),
           # edge = guide_legend(title = 'Relative connectivity',
           #                           override.aes = aes(label = ''),
           #                     name = 'Relative connectivity'),
           size = 'none',
           #width = guide_legend(title = 'Relative connectivity')
    )+
    theme_void() +
    theme(legend.text = element_text(size = 5), 
          legend.title = element_text(size =7),
          legend.box.spacing = unit(1, 'mm'))
  
  ggsave(filename = 'common_asvs_all_1gr_absolut_30perc.pdf', plot = common_asvs_1gr_absolut, 
         path = 'results/figures/corrected_asv_num/',
         width = 188, height = 188, units = 'mm')
  
  
   ##same graph REAL BLOOMERS EXPERIMENTS -----
  ##dades de growth rates
  reg_all_slopes_chosen_silva_tax <- read.csv("data/intermediate_files/reg_all_slopes_chosen_silva_tax_corrected_asv_num.csv", sep=",") %>%
    filter(season != "Early_fall")
  reg_all_slopes_chosen_silva_tax_filt <- reg_all_slopes_chosen_silva_tax %>%
    filter(slope_chosen_days > 0,
           pvalue_slope_chosen < 0.05)
  
  ##filter by asv present at least >1% relative abundance at some treatment or station or replicate
  #rem_fc_relabun <- read_rds("data/rem_fc_relabun.rds")
  rem_relabun_melt <-  read.table("data/rem_relabun_melt.txt", sep="\t", header = TRUE) %>%
    filter(season != "Early_fall")
  
  rem_relabun_melt %$%
    asv_num %>%
    unique() #4594 asv al meu dataset
  
  abundant_asv_tf <- rem_relabun_melt %>% 
    filter(time %in% c('t3', 't4') &
             Abundance > 0.01) %>% #more than 1% of the community at time 3 or t4
    dplyr::select(asv_num, Abundance, treatment, season, time) %>%
    group_by(asv_num, treatment, season, time) %>%
    dplyr::summarize(mean_abundance = mean(Abundance)) %>%
    group_by(asv_num, treatment, season, time, mean_abundance) %>%
    distinct() %>%
    as_tibble() %>%
    pivot_wider(names_from = time, values_from = mean_abundance, id_cols = c(asv_num, treatment, season)) %>%
    mutate(high_mean = case_when (t3 > t4 ~ t3,
                                  t4 > t3 ~ t4, 
                                  is.na(t3) ~ t4,
                                  is.na(t4) ~ t3))
  
  abundant_asv_t0 <- 
    rem_relabun_melt %>% 
    filter(time == 't0') %>%
    dplyr::select(asv_num, Abundance, treatment, season, time) %>%
    group_by(asv_num, treatment, season, time) %>%
    dplyr::summarize(mean_abundance = mean(Abundance)) %>%
    filter(mean_abundance  < 0.01)
  
  abundant_asv_tf_t0 <- abundant_asv_tf %>%
    left_join(abundant_asv_t0, by = c('season', 'treatment', 'asv_num')) #%>%
  # filter(mean_abundance != is.na(mean_abundance)) %>%
  # distinct()
  
  # rem_relabun_melt_1perc <- rem_relabun_melt %>%
  #   right_join(abundant_asv, by = "asv_num", copy = TRUE) #per comprovar amb quin % de relabund estic treballant
  
  # reg_all_slopes_chosen_silva_tax_1perc %>%
  #   colnames()
  
  abundant_asv_tf_t0 %>%
    dim()
  reg_all_slopes_chosen_silva_tax_filt %>%
    dim()
  
  reg_all_slopes_chosen_silva_tax_tf <- reg_all_slopes_chosen_silva_tax_filt %>% 
    right_join(abundant_asv_tf_t0, by = c('treatment' = 'treatment', 'season' = 'season', 'asv_num' = 'asv_num'), copy = TRUE) %>% #by = c"asv_num"
    mutate(#rank_abund = as.numeric(order(rel_abund_ed, decreasing = TRUE)), ##general per tot el dataset
      rank_abund = rank(-high_mean)) %>%
    distinct()
  
  
  data_for_common_asvs <- reg_all_slopes_chosen_silva_tax_tf %>% ##no està filtada per gr >2 
    filter(slope_chosen_days > 2 &
             pvalue_slope_chosen < 0.05) %>%
    dplyr::select(treatment, season, asv_num) %>% #, slope_chosen_days, family_f
    #filter(treatment %in% c('CD', 'CL', 'PD', 'PL')) %>%
    #dplyr::select(-slope_chosen_days, -family_f) %>%
    mutate(seas_treat = paste(treatment, season)) %>%
    group_by(seas_treat, asv_num) %>% ##sembla que PD Winter está repetit al dataset original (trobar a on es duplica, té efecte aals gràfics?)
    distinct(seas_treat, asv_num) 
  
  num_asv_seas_treat <- data_for_common_asvs %>%
    group_by(seas_treat) %>%
    dplyr::summarize(n = n()) %>%
    as_tibble()
  
  ##calculating nº of common ASVs between treatments and seasons------
  ##entre estacions mateix treatment
  CD_CD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CD Spring')      
  CD_CD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CD Summer')      
  CD_CD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CD Fall')  
  CD_CD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CD Summer') 
  CD_CD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CD Fall') 
  CD_CD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'CD Fall') 
  
  CL_CL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CL Spring')      
  CL_CL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CL Summer')      
  CL_CL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CL Fall')  
  CL_CL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CL Summer') 
  CL_CL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CL Fall') 
  CL_CL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'CL Fall') 
  
  PL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PL Spring')      
  PL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PL Summer')      
  PL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PL Fall')  
  PL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PL Summer') 
  PL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PL Fall') 
  PL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'PL Fall') 
  
  PD_PD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PD Spring')      
  PD_PD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PD Summer')      
  PD_PD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PD Fall')  
  PD_PD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PD Summer') 
  PD_PD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PD Fall') 
  PD_PD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'PD Fall')
  
  DL_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'DL Spring')      
  DL_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'DL Summer')      
  DL_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'DL Fall')  
  DL_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'DL Summer') 
  DL_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'DL Fall') 
  DL_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'DL Fall')
  
  VL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'VL Spring')      
  VL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'VL Summer')      
  VL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'VL Fall')  
  VL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'VL Summer') 
  VL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'VL Fall') 
  VL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'VL Fall')
  
  ##entre treatments
  ###CD_CL
  CD_CL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Winter')      
  CD_CL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CL Spring')      
  CD_CL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'CL Summer')  
  CD_CL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Fall', sample2 = 'CL Fall') 
  CD_CL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Spring') 
  CD_CL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CD Spring') 
  CD_CL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Summer') 
  CD_CL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CD Summer') 
  CD_CL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'CL Fall') 
  CD_CL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'CD Fall') 
  CD_CL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CL Summer') 
  CD_CL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CD Summer') 
  CD_CL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'CL Fall') 
  CD_CL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'CD Fall') 
  CD_CL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'CL Fall') 
  CD_CL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'CD Fall') 
  
  ##CD_PL
  CD_PD_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Winter')      
  CD_PD_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Spring')      
  CD_PD_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'CL Summer')  
  CD_PD_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'CL Fall') 
  CD_PD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Spring') 
  CD_PD_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Spring') 
  CD_PD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Summer') 
  CD_PD_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Summer') 
  CD_PD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Fall') 
  CD_PD_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Fall') 
  CD_PD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Summer') 
  CD_PD_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Summer') 
  CD_PD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Fall') 
  CD_PD_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Fall') 
  CD_PD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'CL Fall') 
  CD_PD_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PD Fall')
  
  ##CD_PL
  CD_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Winter')      
  CD_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Spring')      
  CD_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Summer')  
  CD_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Fall', sample2 = 'CL Fall') 
  CD_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Spring') 
  CD_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Spring') 
  CD_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Summer') 
  CD_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Summer') 
  CD_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Fall') 
  CD_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Fall') 
  CD_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Summer') 
  CD_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Summer') 
  CD_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Fall') 
  CD_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Fall') 
  CD_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Fall') 
  CD_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Fall') 
  
  ##CD_DL
  CD_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Winter')      
  CD_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'DL Spring')      
  CD_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'DL Summer')  
  CD_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Fall', sample2 = 'DL Fall') 
  CD_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Spring') 
  CD_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CD Spring') 
  CD_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Summer') 
  CD_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CD Summer') 
  CD_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'DL Fall') 
  CD_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CD Fall') 
  CD_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'DL Summer') 
  CD_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CD Summer') 
  CD_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'DL Fall') 
  CD_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CD Fall') 
  CD_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'DL Fall') 
  CD_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'CD Fall') 
  
  ##CD_VL (16)
  CD_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Winter')      
  CD_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'VL Spring')      
  CD_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'VL Summer')  
  CD_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Fall', sample2 = 'VL Fall') 
  CD_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Spring') 
  CD_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CD Spring') 
  CD_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Summer') 
  CD_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CD Summer') 
  CD_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Winter', sample2 = 'VL Fall') 
  CD_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CD Fall') 
  CD_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'VL Summer') 
  CD_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CD Summer') 
  CD_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Spring', sample2 = 'VL Fall') 
  CD_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CD Fall') 
  CD_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CD Summer', sample2 = 'VL Fall') 
  CD_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'CD Fall') 
  
  
  ##CL_PL (16)
  CL_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Winter')      
  CL_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Spring')      
  CL_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Summer')  
  CL_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'PL Fall') 
  CL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Spring') 
  CL_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Spring') 
  CL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Summer') 
  CL_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Summer') 
  CL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Fall') 
  CL_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Fall') 
  CL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Summer') 
  CL_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Summer') 
  CL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Fall') 
  CL_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Fall') 
  CL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Fall') 
  CL_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Fall') 
  
  ##CL_PD (16)
  CL_PD_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Winter')      
  CL_PD_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Spring')      
  CL_PD_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PD Summer')  
  CL_PD_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'PD Fall') 
  CL_PD_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Spring') 
  CL_PD_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Spring') 
  CL_PD_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Summer') 
  CL_PD_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Summer') 
  CL_PD_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PD Fall') 
  CL_PD_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'CL Fall') 
  CL_PD_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Summer') 
  CL_PD_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Summer') 
  CL_PD_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PD Fall') 
  CL_PD_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'CL Fall') 
  CL_PD_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PD Fall') 
  CL_PD_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'CL Fall') 
  
  ##CL_PL (16)
  CL_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Winter')      
  CL_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Spring')      
  CL_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Summer')  
  CL_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'PL Fall') 
  CL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Spring') 
  CL_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Spring') 
  CL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Summer') 
  CL_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Summer') 
  CL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'PL Fall') 
  CL_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'CL Fall') 
  CL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Summer') 
  CL_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Summer') 
  CL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'PL Fall') 
  CL_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'CL Fall') 
  CL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'PL Fall') 
  CL_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'CL Fall') 
  
  ##CL_DL (16)
  CL_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Winter')      
  CL_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Spring')      
  CL_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Summer')  
  CL_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'DL Fall') 
  CL_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Spring') 
  CL_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Spring') 
  CL_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Summer') 
  CL_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Summer') 
  CL_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Fall') 
  CL_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Fall') 
  CL_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Summer') 
  CL_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Summer') 
  CL_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Fall') 
  CL_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Fall') 
  CL_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Fall') 
  CL_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'CL Fall') 
  
  ##CL_DL (16)
  CL_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Winter')      
  CL_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Spring')      
  CL_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Summer')  
  CL_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'DL Fall') 
  CL_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Spring') 
  CL_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Spring') 
  CL_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Summer') 
  CL_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Summer') 
  CL_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'DL Fall') 
  CL_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'CL Fall') 
  CL_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Summer') 
  CL_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Summer') 
  CL_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'DL Fall') 
  CL_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'CL Fall') 
  CL_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'DL Fall') 
  CL_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'CL Fall') 
  
  ##CL_VL (16)
  CL_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Winter')      
  CL_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'VL Spring')      
  CL_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'VL Summer')  
  CL_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Fall', sample2 = 'VL Fall') 
  CL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Spring') 
  CL_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CL Spring') 
  CL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Summer') 
  CL_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CL Summer') 
  CL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Winter', sample2 = 'VL Fall') 
  CL_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'CL Fall') 
  CL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'VL Summer') 
  CL_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CL Summer') 
  CL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Spring', sample2 = 'VL Fall') 
  CL_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'CL Fall') 
  CL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'CL Summer', sample2 = 'VL Fall') 
  CL_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'CL Fall')
  
  ##PD_PL (16)
  PD_PL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Winter')      
  PD_PL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PL Spring')      
  PD_PL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'PL Summer')  
  PD_PL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'PL Fall') 
  PD_PL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Spring') 
  PD_PL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PD Spring') 
  PD_PL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Summer') 
  PD_PL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PD Summer') 
  PD_PL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'PL Fall') 
  PD_PL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'PD Fall') 
  PD_PL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PL Summer') 
  PD_PL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PD Summer') 
  PD_PL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'PL Fall') 
  PD_PL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'PD Fall') 
  PD_PL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'PL Fall') 
  PD_PL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'PD Fall') 
  
  ##PD_DL (16)
  PD_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Winter')      
  PD_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'DL Spring')      
  PD_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'DL Summer')  
  PD_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'DL Fall') 
  PD_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Spring') 
  PD_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PD Spring') 
  PD_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Summer') 
  PD_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PD Summer') 
  PD_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'DL Fall') 
  PD_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PD Fall') 
  PD_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'DL Summer') 
  PD_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PD Summer') 
  PD_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'DL Fall') 
  PD_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PD Fall') 
  PD_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'DL Fall') 
  PD_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'PD Fall') 
  
  ##PD_VL (16)
  PD_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Winter')      
  PD_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'VL Spring')      
  PD_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'VL Summer')  
  PD_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Fall', sample2 = 'VL Fall') 
  PD_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Spring') 
  PD_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PD Spring') 
  PD_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Summer') 
  PD_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PD Summer') 
  PD_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Winter', sample2 = 'VL Fall') 
  PD_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PD Fall') 
  PD_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'VL Summer') 
  PD_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PD Summer') 
  PD_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Spring', sample2 = 'VL Fall') 
  PD_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PD Fall') 
  PD_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PD Summer', sample2 = 'VL Fall') 
  PD_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'PD Fall') 
  
  ##PL_DL (16)
  PL_DL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Winter')      
  PL_DL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'DL Spring')      
  PL_DL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'DL Summer')  
  PL_DL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Fall', sample2 = 'DL Fall') 
  PL_DL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Spring') 
  PL_DL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PL Spring') 
  PL_DL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Summer') 
  PL_DL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PL Summer') 
  PL_DL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'DL Fall') 
  PL_DL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'PL Fall') 
  PL_DL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'DL Summer') 
  PL_DL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PL Summer') 
  PL_DL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'DL Fall') 
  PL_DL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'PL Fall') 
  PL_DL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'DL Fall') 
  PL_DL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'PL Fall') 
  
  ##PL_VL (16)
  PL_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Winter')      
  PL_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'VL Spring')      
  PL_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'VL Summer')  
  PL_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Fall', sample2 = 'VL Fall') 
  PL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Spring') 
  PL_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PL Spring') 
  PL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Summer') 
  PL_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PL Summer') 
  PL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Winter', sample2 = 'VL Fall') 
  PL_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'PL Fall') 
  PL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'VL Summer') 
  PL_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PL Summer') 
  PL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Spring', sample2 = 'VL Fall') 
  PL_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'PL Fall') 
  PL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'PL Summer', sample2 = 'VL Fall') 
  PL_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'PL Fall') 
  
  ##DL_VL
  DL_VL_W <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Winter')      
  DL_VL_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'VL Spring')      
  DL_VL_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'VL Summer')  
  DL_VL_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Fall', sample2 = 'VL Fall') 
  DL_VL_W_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Spring') 
  DL_VL_Sp_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'DL Spring') 
  DL_VL_W_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Summer') 
  DL_VL_Su_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'DL Summer') 
  DL_VL_W_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Winter', sample2 = 'VL Fall') 
  DL_VL_F_W <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Winter', sample2 = 'DL Fall') 
  DL_VL_Sp_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'VL Summer') 
  DL_VL_Su_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'DL Summer') 
  DL_VL_Sp_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Spring', sample2 = 'VL Fall') 
  DL_VL_F_Sp <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Spring', sample2 = 'DL Fall') 
  DL_VL_Su_F <- common.asvs(data = data_for_common_asvs, sample1 = 'DL Summer', sample2 = 'VL Fall') 
  DL_VL_F_Su <- common.asvs(data = data_for_common_asvs, sample1 = 'VL Summer', sample2 = 'DL Fall') 
  
  common_asv_sum <- 
    bind_rows(CD_CD_W_Sp, CD_CD_W_Su, CD_CD_W_F, CD_CD_Sp_Su, CD_CD_Sp_F, CD_CD_Su_F,
              CL_CL_W_Sp, CL_CL_W_Su, CL_CL_W_F, CL_CL_Sp_Su, CL_CL_Sp_F, CL_CL_Su_F,
              PL_PL_W_Sp, PL_PL_W_Su, PL_PL_W_F, PL_PL_Sp_Su, PL_PL_Sp_F, PL_PL_Su_F,
              PD_PD_W_Sp, PD_PD_W_Su, PD_PD_W_F, PD_PD_Sp_Su, PD_PD_Sp_F, PD_PD_Su_F,
              DL_DL_W_Sp, DL_DL_W_Su, DL_DL_W_F, DL_DL_Sp_Su, DL_DL_Sp_F, DL_DL_Su_F,
              VL_VL_W_Sp, VL_VL_W_Su, VL_VL_W_F, VL_VL_Sp_Su, VL_VL_Sp_F, VL_VL_Su_F,
              CD_CL_W, CD_CL_Sp, CD_CL_Su, CD_CL_F, 
              CD_CL_W_Sp, CD_CL_Sp_W, CD_CL_W_Su, CD_CL_Su_W, CD_CL_W_F, CD_CL_F_W, CD_CL_Sp_Su , CD_CL_Su_Sp, CD_CL_Sp_F, CD_CL_F_Sp, CD_CL_Su_F, CD_CL_F_Su, 
              CD_PD_W, CD_PD_Sp, CD_PD_Su, CD_PD_F,
              CD_PD_W_Sp, CD_PD_Sp_W, CD_PD_W_Su, CD_PD_Su_W, CD_PD_W_F, CD_PD_F_W, CD_PD_Sp_Su , CD_PD_Su_Sp, CD_PD_Sp_F, CD_PD_F_Sp, CD_PD_Su_F, CD_PD_F_Su,
              CD_PL_W, CD_PL_Sp, CD_PL_Su, CD_PL_F,
              CD_PL_W_Sp, CD_PL_Sp_W, CD_PL_W_Su, CD_PL_Su_W, CD_PL_W_F, CD_PL_F_W, CD_PL_Sp_Su , CD_PL_Su_Sp, CD_PL_Sp_F, CD_PL_F_Sp, CD_PL_Su_F, CD_PL_F_Su,
              CD_DL_W, CD_DL_Sp, CD_DL_Su, CD_DL_F,
              CD_DL_W_Sp, CD_DL_Sp_W, CD_DL_W_Su, CD_DL_Su_W, CD_DL_W_F, CD_DL_F_W, CD_DL_Sp_Su , CD_DL_Su_Sp, CD_DL_Sp_F, CD_DL_F_Sp, CD_DL_Su_F, CD_DL_F_Su, 
              CD_VL_Sp, CD_VL_Su, CD_VL_F,
              CD_VL_W_Sp, CD_VL_Sp_W, CD_VL_W_Su, CD_VL_Su_W, CD_VL_W_F, CD_VL_F_W, CD_VL_Sp_Su , CD_VL_Su_Sp, CD_VL_Sp_F, CD_VL_F_Sp, CD_VL_Su_F, CD_VL_F_Su, 
              PD_PL_W, PD_PL_Sp, PD_PL_Su, PD_PL_F,
              PD_PL_W_Sp, PD_PL_Sp_W, PD_PL_W_Su, PD_PL_Su_W, PD_PL_W_F, PD_PL_F_W, PD_PL_Sp_Su , PD_PL_Su_Sp, PD_PL_Sp_F, PD_PL_F_Sp, PD_PL_Su_F, PD_PL_F_Su,
              PD_DL_W, PD_DL_Sp, PD_DL_Su, PD_DL_F,
              PD_DL_W_Sp, PD_DL_Sp_W, PD_DL_W_Su, PD_DL_Su_W, PD_DL_W_F, PD_DL_F_W, PD_DL_Sp_Su , PD_DL_Su_Sp, PD_DL_Sp_F, PD_DL_F_Sp, PD_DL_Su_F, PD_DL_F_Su, 
              PD_VL_Sp, PD_VL_Su, PD_VL_F,
              PD_VL_W_Sp, PD_VL_Sp_W, PD_VL_W_Su, PD_VL_Su_W, PD_VL_W_F, PD_VL_F_W, PD_VL_Sp_Su , PD_VL_Su_Sp, PD_VL_Sp_F, PD_VL_F_Sp, PD_VL_Su_F, PD_VL_F_Su,                 
              PL_DL_W, PL_DL_Sp, PL_DL_Su, PL_DL_F,
              PL_DL_W_Sp, PL_DL_Sp_W, PL_DL_W_Su, PL_DL_Su_W, PL_DL_W_F, PL_DL_F_W, PL_DL_Sp_Su , PL_DL_Su_Sp, PL_DL_Sp_F, PL_DL_F_Sp, PL_DL_Su_F, PL_DL_F_Su, 
              PL_VL_Sp, PL_VL_Su, PL_VL_F,
              PL_VL_W_Sp, PL_VL_Sp_W, PL_VL_W_Su, PL_VL_Su_W, PL_VL_W_F, PL_VL_F_W, PL_VL_Sp_Su , PL_VL_Su_Sp, PL_VL_Sp_F, PL_VL_F_Sp, PL_VL_Su_F, PL_VL_F_Su,                 
              DL_VL_Sp, DL_VL_Su, DL_VL_F,
              DL_VL_W_Sp, DL_VL_Sp_W, DL_VL_W_Su, DL_VL_Su_W, DL_VL_W_F, DL_VL_F_W, DL_VL_Sp_Su , DL_VL_Su_Sp, DL_VL_Sp_F, DL_VL_F_Sp, DL_VL_Su_F, DL_VL_F_Su) %>%
    separate(seas_treat.x, '_', into = c('from', 'to')) %>%
    left_join(num_asv_seas_treat, by = c('from' = 'seas_treat')) %>%
    left_join(num_asv_seas_treat, by = c('to' = 'seas_treat')) %>%
    mutate(division = case_when( n.x < n.y ~ n.x,
                                 n.y < n.x ~ n.y,
                                 n.x == n.y ~ n.x)) %>%
    mutate(relative_connectivity = num_common_asvs/division) %>%
    mutate(sample = paste0(from, to))  #%>%
  # select(from, to, relative_connectivity) 
  
  common_asv_sum 
  
  
  ##prepare graph
  edges <- common_asv_sum %>%
    select(from, to) %>%
    filter(from != is.na(from))
  vertices <- common_asv_sum %>%
    select(sample) %>%
    mutate(size = 10) %>%
    as_tibble()
  
  vertices <- data.frame(rbind(as_tibble(common_asv_sum$from), 
                               as_tibble(common_asv_sum$to))) %>%
    group_by(value) %>%
    unique() %>%
    mutate(size = 10) %>%
    separate(value, ' ', into = c('Treatment', 'Season'), remove = F) %>%
    mutate(Treatment = as.factor(Treatment),
           Season = as.factor(Season),
           value = as.factor(value)) %>%
    select(value, Season, Treatment) %>%
    filter(value != is.na(value)) %>%
    as_tibble()
  
  vertices$Season <- vertices$Season %>% 
    factor(levels = c('Winter', 'Spring', 'Summer', 'Fall'))
  
  vertices$Treatment <- vertices$Treatment %>% 
    factor(levels = c('CD', 'CL', 'PD', 'PL', 'DL', 'VL'))
  
  vertices <- vertices %>%
    mutate(value = fct_relevel(value, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                        "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                        "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                        "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
    
    arrange(value) %>%
    as_tibble()
  
  # vertices %$%
  #   value
  connections <- common_asv_sum %>%
    separate(from, ' ', into = c('Treatment', 'Season'), remove = F) %>%
    filter(relative_connectivity != is.na(relative_connectivity)) %>%
    select(from, to, relative_connectivity, Treatment, Season) 
  
  
  num_connections1 <- common_asv_sum %>%
    group_by(from) %>%
    dplyr::summarize(n = n()) %>%
    right_join(vertices, by = c('from' = 'value'))
  
  
  num_connections2 <- common_asv_sum %>%
    group_by(to) %>%
    dplyr::summarize(n = n()) %>%
    right_join(vertices, by = c('to' = 'value'))
  
  num_connections12 <- num_connections1 %>%
    left_join(num_connections2, by=c('Season' = 'Season', 'Treatment' = 'Treatment')) %>%
    mutate(n.x = case_when(#!is.na(n.x) ~ n.x,
      is.na(n.x) ~ '0',
      TRUE ~ as.character(n.x)),
      n.y = case_when(#!is.na(n.x) ~ n.x,
        is.na(n.y) ~ '0',
        TRUE ~ as.character(n.y))
    ) %>%
    mutate(num_connections = (as.numeric(n.x)+as.numeric(n.y)))
  
  common_asv_sum$from
  
  common_asv_sum$to
  
  # # Preparation to draw labels properly:
  # vertices$id=NA
  # myleaves=which(is.na( match(vertices$name, edges$from) ))
  # nleaves=length(myleaves)
  # vertices$id[ myleaves ] = seq(1:nleaves)
  # vertices$angle= 90 - 360 * vertices$id / nleaves
  # vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)
  # vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)
  
  # Build a network object from this dataset:
  mygraph <- graph_from_data_frame(edges, vertices = vertices)
  #graph_test <-  tbl_graph(edges, vertices)
  ##reorder nodes
  # s <- arrange(names(V(mygraph)), vertices$Season)
  # mygraph <- permute(mygraph, match(V(mygraph)$name, s))
  
  # V(mygraph)$node_label <- names(V(mygraph)) 
  
  # The connection object must refer to the ids of the leaves:
  #from = match( connections$Treatment$from, vertices$Treatment$value) (no funciona)
  from = match( connections$from, vertices$value)
  to = match( connections$to, vertices$value)
  
  # Basic dendrogram
  # ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  #   geom_conn_bundle(data = get_con(from = from, to = to, paths = connections$relative_connectivity), alpha = 0.1, colour="#69b3a2") +
  #   #geom_edge_link(size=0.4, alpha=0.1) +
  #   #geom_node_text(aes(x = x*1.01, y=y*1.01, filter = leaf, label=vertices$value, angle = angle, hjust=hjust), size=1.5, alpha=1) +
  #   coord_fixed() +
  #   theme_void() +
  #   theme(
  #     legend.position="none",
  #     plot.margin=unit(c(0,0,0,0),"cm"),
  #   ) +
  #   expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))
  
  ##CIRCULAR CONNECTIONS 
  specific_range <- common_asv_sum %$%
    relative_connectivity %>%
    range()
  
  palete_gradient <- c(#"#005a47",
    #"#678078",
    #"#aaaaaa",
    "#e6e6e6",
             #"#7b7b7b",
             "#003318") 
             
  
  common_asvs_realbloom <- ggraph(mygraph, layout="linear", circular = TRUE) + 
    # geom_conn_bundle(data = get_con(from = from, to = to, alpha = 0.1, colour="#e6e6e6") #alpha=0.2, , tension = 0.9,
    # )+
    geom_edge_arc(aes(width = as.numeric(connections$relative_connectivity), 
                      color = as.numeric(connections$relative_connectivity)
    ),
    fold=FALSE,
    edge_alpha=0.8,
    lineend = 'square',
    linejoin = 'mitre', check_overlap = TRUE, #strength = 0.5
    # title = 'Relative connectivity'
    # edge_width = 1
    ) + #, position = position_jitter(0.5)
    #geom_edge_density(aes(fill = as.numeric(connections$relative_connectivity)))+
    # geom_edge_fan0(aes(#alpha = after_stat(connections$relative_connectivity,
    #                                      color = connections$Treatment))+
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    # scale_color_manual(values = palette_treatments_remei)+
    scale_edge_colour_gradientn(colours = palete_gradient, name = 'Relative connectivity')+
    scale_edge_width_continuous(limits = c(0.3, 0.89), guide = 'none')+
    #geom_edge_diagonal(aes(width = connections$relative_connectivity, color = connections$relative_connectivity))+
    #geom_ar_link(aes(width = connections$relative_connectivity))
    #scale_color_manual(values=mycolor) +
    coord_fixed()+
    geom_node_point(aes(color = Season, size = num_connections12$num_connections), alpha = 0.9)+#, size = 0.75
    scale_size(range=c(8,19), limits = c(8,19)) +
    geom_node_text(aes(label=Treatment), vjust = -0.5, show.legend = TRUE) + #, angle=65, hjust=1, nudge_y = -1.1, size=2.3
    scale_color_manual(values = palette_seasons_4)+
    guides(color = guide_legend(ncol = 1, size = 10,
                                override.aes = aes(label = '')),
           shape = guide_legend(ncol = 1, size = 1),
           # edge = guide_legend(title = 'Relative connectivity',
           #                           override.aes = aes(label = ''),
           #                     name = 'Relative connectivity'),
           size = 'none',
           #width = guide_legend(title = 'Relative connectivity')
    )+
    theme_void() +
    theme(legend.text = element_text(size = 5), 
          legend.title = element_text(size =7))
  
  
  ggsave(filename = 'common_asvs_all_realbloom.pdf', plot = common_asvs_realbloom, 
         path = 'results/figures/corrected_asv_num/',
         width = 188, height = 188, units = 'mm')

  ###Number of EXCLUSIVE ASV growing >1gr/day per sample----------------
  data_for_exclusive_asvs <- reg_all_slopes_chosen_silva_tax %>%
    filter(slope_chosen_days > 1 &
             pvalue_slope_chosen < 0.05) %>%
    dplyr::select(treatment, season, asv_num) %>% #, slope_chosen_days, family_f
    #filter(treatment %in% c('CD', 'CL', 'PD', 'PL')) %>%
    #dplyr::select(-slope_chosen_days, -family_f) %>%
    mutate(seas_treat = paste(treatment, season)) %>%
    group_by(seas_treat, asv_num) %>% ##sembla que PD Winter está repetit al dataset original (trobar a on es duplica, té efecte als gràfics?)
    distinct(seas_treat, asv_num) 
  
  
  total_responding_asvs_sample  <-  data_for_exclusive_asvs %>%
    group_by(seas_treat) %>%
    dplyr::summarize(num_responding_asv = n()) %>%
    mutate(num_responding_asv = as.numeric(num_responding_asv)) %>%
    as_tibble()

# asvs_present_all_other_samples <-  data_for_common_asvs %>%
#     filter(seas_treat != 'CD Winter') %>%
#     group_by(asv_num) %>%
#     distinct()
# 
# CD_winter_exclusive_asvs_gr2 <- data_for_common_asvs %>%
#   filter(seas_treat == 'CD Winter') %>%
#   filter(!asv_num %in% asvs_present_all_other_samples$asv_num)

  ##we filter all ASVs that >1GR/day (other than 0)
  
  exclusive.asvs <- function(data, sample){
    asvs_present_all_other_samples <- data %>%
      dplyr::filter(seas_treat != sample) %>%
      group_by(asv_num) %>%
      distinct()
    exclusive_asv_sample <- data %>%
      dplyr::filter(seas_treat == sample) %>%
      filter(!asv_num %in% asvs_present_all_other_samples$asv_num)
    return(exclusive_asv_sample)
      
  }

  excl_cd_w <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'CD Winter')
  excl_cd_sp <-   exclusive.asvs(data_for_exclusive_asvs, sample = 'CD Spring')
  excl_cd_su <-   exclusive.asvs(data_for_exclusive_asvs, sample = 'CD Summer')
  excl_cd_f <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'CD Fall')
  excl_cl_w <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'CL Winter')
  excl_cl_sp <- exclusive.asvs(data_for_exclusive_asvs, sample = 'CL Spring')
  excl_cl_su <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'CL Summer')
  excl_cl_f <- exclusive.asvs(data_for_exclusive_asvs, sample = 'CL Fall')
  excl_PD_w <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'PD Winter')
  excl_PD_sp <- exclusive.asvs(data_for_exclusive_asvs, sample = 'PD Spring')
  excl_PD_su <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'PD Summer')
  excl_PD_f <- exclusive.asvs(data_for_exclusive_asvs, sample = 'PD Fall')
  excl_PL_w <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'PL Winter')
  excl_PL_sp <- exclusive.asvs(data_for_exclusive_asvs, sample = 'PL Spring')
  excl_PL_su <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'PL Summer')
  excl_PL_f <- exclusive.asvs(data_for_exclusive_asvs, sample = 'PL Fall')
  excl_DL_w <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'DL Winter')
  excl_DL_sp <- exclusive.asvs(data_for_exclusive_asvs, sample = 'DL Spring')
  excl_DL_su <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'DL Summer')
  excl_DL_f <- exclusive.asvs(data_for_exclusive_asvs, sample = 'DL Fall')
  excl_VL_w <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'VL Winter')
  excl_VL_sp <- exclusive.asvs(data_for_exclusive_asvs, sample = 'VL Spring')
  excl_VL_su <-  exclusive.asvs(data_for_exclusive_asvs, sample = 'VL Summer')
  excl_VL_f <- exclusive.asvs(data_for_exclusive_asvs, sample = 'VL Fall')

exclusive_asvs_dataset <-  bind_rows(excl_cd_w ,   excl_cd_sp ,   excl_cd_su ,   excl_cd_f ,   excl_cl_w , 
            excl_cl_sp ,   excl_cl_su ,   excl_cl_f ,   excl_PD_w ,   excl_PD_sp ,
            excl_PD_su ,   excl_PD_f ,   excl_PL_w,   excl_PL_sp ,   excl_PL_su,
            excl_PL_f ,   excl_DL_w,   excl_DL_sp ,   excl_DL_su,   excl_DL_f , 
            excl_VL_w ,   excl_VL_sp ,   excl_VL_su,   excl_VL_f) %>%
  group_by(seas_treat) %>%
  dplyr::summarize(number_exclusive_asvs = n()) %>%
  separate(seas_treat, ' ', into = c('Treatment', 'Season'), remove = F) %>%
  mutate(seas_treat = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                      "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                      "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                      "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall")))
exclusive_asvs_dataset %$%
  number_exclusive_asvs %>%
  range()
exclusive_asvs_absolut <- exclusive_asvs_dataset %>%
  ggplot(aes(seas_treat, number_exclusive_asvs, fill = Season))+
  geom_col()+
  geom_text(aes(label = (Treatment), y = (110)))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_fill_manual(values = palette_seasons_4)+
  annotate('text', x = 0, y = c(20, 40, 60, 80, 100), label = c('20', '40', '60', '80', '100')) +
  scale_y_continuous(limits = c(0, 110))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(filename = 'exclusive_asvs_absolut_1gr.pdf', plot = exclusive_asvs_absolut, 
       path = 'results/figures/corrected_asv_num/',
       width = 88, height = 88, units = 'mm')


##fer el mateix però en relatiu també
exclusive_asvs_relative <-  exclusive_asvs_dataset %>%
  left_join(total_responding_asvs_sample, by = 'seas_treat') %>%
  mutate(number_exclusive_asvs = as.numeric(number_exclusive_asvs)) %>%
  as_tibble() %>%
  mutate(relative_exclusive_asvs = number_exclusive_asvs/num_responding_asv) %>%
  mutate(seas_treat = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall")))

exclusive_asvs_relative_plot <- exclusive_asvs_relative %>%
  ggplot(aes(seas_treat, relative_exclusive_asvs, fill = Season))+
  geom_col()+
  geom_text(aes(label = (Treatment), y = (0.55)))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_fill_manual(values = palette_seasons_4)+
  annotate('text', x = 0, y = c(0.2, 0.4), label = c('0.2', '0.4')) +
  scale_y_continuous(limits = c(0, 0.55))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

ggsave(filename = 'exclusive_asvs_relative_1gr.pdf', plot = exclusive_asvs_relative_plot, 
       path = 'results/figures/corrected_asv_num/',
       width = 88, height = 88, units = 'mm')

####Total growing ASVs (GR > 0)

total_growing_asvs <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0 &
           pvalue_slope_chosen < 0.05) %>%
  dplyr::select(treatment, season, asv_num) %>% #, slope_chosen_days, family_f
  #filter(treatment %in% c('CD', 'CL', 'PD', 'PL')) %>%
  #dplyr::select(-slope_chosen_days, -family_f) %>%
  mutate(seas_treat = paste(treatment, season)) %>%
  group_by(seas_treat, asv_num) %>% ##sembla que PD Winter está repetit al dataset original (trobar a on es duplica, té efecte aals gràfics?)
  distinct(seas_treat, asv_num) %>%
  group_by(seas_treat) %>%
  dplyr::summarize(num_growing_asv = n()) %>%
  #mutate(num_responding_asv = as.numeric(num_responding_asv)) %>%
  as_tibble()

total_growing_asvs <-  total_growing_asvs %>%
 #left_join(, by = 'seas_treat') %>%
  separate('seas_treat', sep = ' ', into = c('Treatment', 'Season'), remove = F) %>%
  #mutate(number_exclusive_asvs = as.numeric(number_exclusive_asvs)) %>%
  as_tibble() %>%
  #mutate(relative_exclusive_asvs = number_exclusive_asvs/num_responding_asv) %>%
  mutate(seas_treat = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall")))
total_growing_asvs %$%
  num_growing_asv %>%
  range()

total_growing_asvs_plot_absolut <- 
  total_growing_asvs %>%
  ggplot(aes(seas_treat, num_growing_asv, fill = Season))+
  geom_col()+
  geom_text(aes(label = (Treatment), y = (300)))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_fill_manual(values = palette_seasons_4)+
  annotate('text', x = 0, y = c(100, 200, 300), label = c('100', '200', '300')) +
  scale_y_continuous(limits = c(0, 300))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

ggsave(filename = 'total_growing_asvs_plot_absolut.pdf', plot = total_growing_asvs_plot_absolut, 
       path = 'results/figures/corrected_asv_num/',
       width = 88, height = 88, units = 'mm')

##  per fer el relatiu comptem nº ASV presents a cada mostra i dividim les growing per aquest valor
rem_relabun_melt <-  read.table("data/rem_relabun_melt.txt", sep="\t", header = TRUE) %>%
filter(season != "Early_fall")

rem_relabun_melt %>%
  head()

total_growing_asvs_relative <- rem_relabun_melt %>%
  dplyr::filter(Abundance != 0 &
                  time == 't0') %>%
  group_by(treatment, season, replicate) %>%
  dplyr::summarize(present_asvs = n()) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_asv_num = mean(present_asvs)) %>%
  left_join(total_growing_asvs, by = c('treatment' = 'Treatment', 'season' = 'Season')) %>%
  mutate(relative_growing_asvs = num_growing_asv/mean_asv_num) #%>%
 # ungroup() %>%
  #mutate(relative_growing_asvs_mean = mean(relative_growing_asvs)) 
total_growing_asvs_relative %$%
  relative_growing_asvs %>%
  range()

total_growing_asvs_plot_relative <- 
  total_growing_asvs_relative %>%
  ggplot(aes(seas_treat, relative_growing_asvs, fill = season))+
  geom_col(alpha = 0.9)+
  geom_text(aes(label = (treatment), y = (0.65)))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_fill_manual(values = palette_seasons_4)+
  annotate('text', x = 0, y = c(0.2, 0.4, 0.6), label = c('0.2', '0.4', '0.6')) +
  scale_y_continuous(limits = c(0, 0.65))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text = element_blank(),
        #axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

ggsave(filename = 'total_growing_asvs_plot_relative.pdf', plot = total_growing_asvs_plot_relative, 
       path = 'results/figures/corrected_asv_num/',
       width = 88, height = 88, units = 'mm')


##constructing relative pannel for those ASVs that are growing >2 GR/day----
# library(gridExtra)
# grid.arrange(exclusive_asvs_relative_plot, common_asvs_2gr)

library(multipanelfigure)
common_exclusive_relative_asvs <- multi_panel_figure(columns = 3, rows = 1, width = 188, height = 120, 
                                                  column_spacing = 0, unit = 'mm',
                                                  panel_label_type = 'upper-alpha') #si no el volem none


common_exclusive_relative_asvs %<>%
  fill_panel(common_asvs_2gr, column = 1:2, row = 1) %<>%
  fill_panel(exclusive_asvs_relative_plot, column = 3, row = 1) 

##CAMBIAR EL NOM ABANS DE GUARDAR!!!!
ggsave('common_exclusive_relative_asvs_2gr.pdf', common_exclusive_relative_asvs,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 188,
       height = 120,
       units = 'mm')

#number of growing ASVs per treatment depending on nº of total ASVs present in that sample-----

##diversity in growth rates (Is the community reacting more )-----
# SHANNON DIVERSITY
# creamos la función de shannon
shannon_div <- function(vector){
  vector <- vector*log(vector)
  # the vector has NA values if the species proportions are 0
  vectorsum <- sum(na.omit(vector))
  return((-1)*vectorsum)
}

##entro la taula d'ASVs------
##This phyloseq object has in sample data a column with FC counts to calculate pseudoabundances
rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_silva_fc.rds")
#rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_gtbd_fc.rds")

##information from my dataset----
rem_fc %>% 
  summarize_phyloseq()

##reorder the environmental data
rem_fc@sam_data$treatment <- rem_fc@sam_data$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
rem_fc@sam_data$season <- rem_fc@sam_data$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Early_fall", "Fall")))

##extract some information from our dataset----
rem_fc %>% 
  ntaxa() ##5670 asv silva // gtbd 4599
rem_fc %>% 
  nsamples() ##312 samples

##subseting by library size----
nsamples(rem_fc) #312
rem_fc_filt <- prune_samples(sample_sums(rem_fc) > 5000, rem_fc)
nsamples(rem_fc_filt) #308 hi ha 4 samples per sota dels 5,000 silva // amb gtbd n'hi ha 7
#write_rds(rem_fc_filt, "data/intermediate_files/rem_fc_filt.rds")



##filter all asv that sum 0-----
# rem_fc_filt@otu_table %>% 
#   colSums()
rem_fc_filt <- prune_taxa(taxa_sums(rem_fc_filt) >0, rem_fc_filt)
rem_fc_filt#4594 asv silva // gtbd 3810

##rel_abund
rem_fc_relabun <- phyloseq::transform_sample_counts(rem_fc_filt, 
                                                    function(x)
                                                    {x / sum(x)})

asv_tab <- rem_fc_relabun@otu_table 
sam_data <- rem_fc_filt@sam_data

#Shannon diversity plot-----
shannon_diversities1 <- apply(asv_tab, 1, shannon_div) #normal
sam_data$shannon <- shannon_diversities1

sam_data %>%
  colnames()

sam_data$season <- sam_data$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))


t0 <- sam_data %>%
  as_tibble() %>%
  subset(time == 't0' &
           season != 'Early_fall') %>%
  group_by(treatment, season, time, hours.x) %>%
  dplyr::summarize(mean_shannon = mean(shannon)) %>%
  mutate(seas_treat = as.factor(paste(treatment,season, sep = ' '))) %>%
  mutate(seas_treat2 = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
  ggplot(aes(seas_treat2, mean_shannon))+
  geom_col(aes(fill = season), alpha = 0.9)+
  scale_fill_manual(values = palette_seasons_4)+
  geom_text(aes(label = (treatment), y = (5.5)))+
   coord_polar()+
  labs(x = '', y = '')+#, title = 'T0'
  annotate('text', x = 0.2, y = c(2, 4), label = c('2', '4')) +
  #annotate('segment', x = 0, xend = 0, y = 0, yend =5.5)+
  scale_y_continuous(limits = c(0, 5.5))+
  #scale_y_continuous(limits = c(0, 0.55))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

t2 <-  sam_data %>%
  as_tibble() %>%
  subset(time == 't2' &
           season != 'Early_fall') %>%
  group_by(treatment, season, time, hours.x) %>%
  dplyr::summarize(mean_shannon = mean(shannon)) %>%
  mutate(seas_treat = as.factor(paste(treatment,season, sep = ' '))) %>%
  mutate(seas_treat2 = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                 "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                 "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                 "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
  ggplot(aes(seas_treat2, mean_shannon))+
  geom_col(aes(fill = season), alpha = 0.9)+
  scale_fill_manual(values = palette_seasons_4)+
  geom_text(aes(label = (treatment), y = (5.5)))+
  coord_polar()+
  labs(x = '', y = '', title = 'T2')+
  annotate('text', x = 0.2, y = c(2, 4), label = c('2', '4')) +
 # annotate('segment', x = 0, xend = 0, y = 0, yend =5)+
  scale_y_continuous(limits = c(0, 5.5))+
  #scale_y_continuous(limits = c(0, 0.55))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

t3 <- sam_data %>%
  as_tibble() %>%
  subset(time == 't3' &
           season != 'Early_fall') %>%
  group_by(treatment, season, time, hours.x) %>%
  dplyr::summarize(mean_shannon = mean(shannon)) %>%
  mutate(seas_treat = as.factor(paste(treatment,season, sep = ' '))) %>%
  mutate(seas_treat2 = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                 "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                 "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                 "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
  ggplot(aes(seas_treat2, mean_shannon))+
  geom_col(aes(fill = season), alpha = 0.9)+
  scale_fill_manual(values = palette_seasons_4)+
  geom_text(aes(label = (treatment), y = (5.5)))+
  coord_polar()+
  labs(x = '', y = '', title = 'T3')+
  annotate('text', x = 0.2, y = c(2, 4), label = c('2', '4')) +
  scale_y_continuous(limits = c(0, 5.5))+
  #annotate('segment', x = 0, xend = 0, y = 0, yend =5)+
  #scale_y_continuous(limits = c(0, 0.55))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

t4 <- sam_data %>%
  as_tibble() %>%
  subset(time == 't4' &
           season != 'Early_fall') %>%
  group_by(treatment, season, time, hours.x) %>%
  dplyr::summarize(mean_shannon = mean(shannon)) %>%
  mutate(seas_treat = as.factor(paste(treatment,season, sep = ' '))) %>%
  mutate(seas_treat2 = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                 "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                 "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                 "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
  ggplot(aes(seas_treat2, mean_shannon))+
  geom_col(aes(fill = season), alpha = 0.9)+
  scale_fill_manual(values = palette_seasons_4)+
  geom_text(aes(label = (treatment), y = (5.5)))+
  coord_polar()+
  labs(x = '', y = '', title = 'T4')+
  annotate('text', x = 0.2, y = c(2, 4), label = c('2', '4')) +
 # annotate('segment', x = 0, xend = 0, y = 0, yend =5)+
  scale_y_continuous(limits = c(0, 5.5))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

shannon_time <- grid.arrange(t0,t2,t3,t4, nrow = 1)

ggsave(filename = 'shannon_time.pdf', plot = shannon_time_cols, 
       path = 'results/figures/corrected_asv_num/',
       width = 288, height = 100, units = 'mm')

##Richness-----
rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_silva_fc.rds")
#rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_gtbd_fc.rds")

##information from my dataset----
rem_fc %>% 
  summarize_phyloseq()

##reorder the environmental data
rem_fc@sam_data$treatment <- rem_fc@sam_data$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
rem_fc@sam_data$season <- rem_fc@sam_data$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Early_fall", "Fall")))

##extract some information from our dataset----
rem_fc %>% 
  ntaxa() ##5670 asv silva // gtbd 4599
rem_fc %>% 
  nsamples() ##312 samples

##subseting by library size----
nsamples(rem_fc) #312
rem_fc_filt <- prune_samples(sample_sums(rem_fc) > 5000, rem_fc)

sam_data <- sample_data(rem_fc_filt) %>%
  as_tibble() %>%
  dplyr::mutate(sample_sum = sample_sums(rem_fc_filt)) %>%
  dplyr::rename('Sample' = '.sample') %>%
  as_tibble() %>%
  dplyr::filter(season != 'Early_fall') %>%
  mutate(Sample = str_replace(Sample, '6642', 'X6642'))

richness <- rem_fc_filt %>%
  estimate_richness() %>%
  rownames_to_column(var = 'Sample') %>%
  #mutate(time_ed = str_replace_all(time, c('t4'= 'tf', 't3' = 'tf')))
  mutate(Sample = str_replace_all(Sample, c('[.]' = '-'))) %>% #'X6642' = '6642',
  #mutate(Sample = str_replace(Sample, '[.]' = '-' )) %>%
 # mutate(Sample = str_replace_all(Sample,'*6642' = '6642'))
  #as_tibble() %>%
  right_join(sam_data, by = 'Sample')

richness$Sample == sam_data$Sample #check just in case

##Richness plot----
richness$season <- richness$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

richness %>%
  group_by(treatment, season, time) %>%
  dplyr::summarize(mean_richness = mean(Observed)) %$%
  mean_richness %>%
  range()

t0 <- richness %>%
  #as_tibble() %>%
  subset(time == 't0' &
           season != 'Early_fall') %>%
  group_by(treatment, season, time) %>%
  dplyr::summarize(mean_richness = mean(Observed)) %>%
  mutate(seas_treat = as.factor(paste(treatment,season, sep = ' '))) %>%
  mutate(seas_treat2 = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                 "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                 "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                 "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
  ggplot(aes(seas_treat2, mean_richness))+
  geom_col(aes(fill = season), alpha = 0.9)+
  scale_fill_manual(values = palette_seasons_4)+
  geom_text(aes(label = (treatment), y = (1000)))+
  coord_polar()+
  labs(x = '', y = '')+#, title = 'T0'
  annotate('text', x = 0, y = c(250, 500, 750, 1000), label = c('250', '500', '750', '1000')) +
  annotate('segment', x = 0, xend = 0, y = 0, yend =5.5)+
  scale_y_continuous(limits = c(0, 1015))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))


t2 <- richness %>%
  #as_tibble() %>%
  subset(time == 't2' &
           season != 'Early_fall') %>%
  group_by(treatment, season, time) %>%
  dplyr::summarize(mean_richness = mean(Observed)) %>%
  mutate(seas_treat = as.factor(paste(treatment,season, sep = ' '))) %>%
  mutate(seas_treat2 = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                 "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                 "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                 "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
  ggplot(aes(seas_treat2, mean_richness))+
  geom_col(aes(fill = season), alpha = 0.9)+
  scale_fill_manual(values = palette_seasons_4)+
  geom_text(aes(label = (treatment), y = (1000)))+
  coord_polar()+
  labs(x = '', y = '')+#, title = 'T0'
  annotate('text', x = 0, y = c(250, 500, 750, 1000), label = c('250', '500', '750', '1000')) +
  #annotate('segment', x = 0, xend = 0, y = 0, yend =5.5)+
  scale_y_continuous(limits = c(0, 1015))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))


t3 <- richness %>%
  #as_tibble() %>%
  subset(time == 't3' &
           season != 'Early_fall') %>%
  group_by(treatment, season, time) %>%
  dplyr::summarize(mean_richness = mean(Observed)) %>%
  mutate(seas_treat = as.factor(paste(treatment,season, sep = ' '))) %>%
  mutate(seas_treat2 = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                 "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                 "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                 "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
  ggplot(aes(seas_treat2, mean_richness))+
  geom_col(aes(fill = season), alpha = 0.9)+
  scale_fill_manual(values = palette_seasons_4)+
  geom_text(aes(label = (treatment), y = (1000)))+
  coord_polar()+
  labs(x = '', y = '')+#, title = 'T0'
  annotate('text', x = 0, y = c(250, 500, 750, 1000), label = c('250', '500', '750', '1000')) +
  #annotate('segment', x = 0, xend = 0, y = 0, yend =5.5)+
  scale_y_continuous(limits = c(0, 1015))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

t4 <- richness %>%
  #as_tibble() %>%
  subset(time == 't4' &
           season != 'Early_fall') %>%
  group_by(treatment, season, time) %>%
  dplyr::summarize(mean_richness = mean(Observed)) %>%
  mutate(seas_treat = as.factor(paste(treatment,season, sep = ' '))) %>%
  mutate(seas_treat2 = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                 "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                 "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                 "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) %>%
  ggplot(aes(seas_treat2, mean_richness))+
  geom_col(aes(fill = season), alpha = 0.9)+
  scale_fill_manual(values = palette_seasons_4)+
  geom_text(aes(label = (treatment), y = (1000)))+
  coord_polar()+
  labs(x = '', y = '')+#, title = 'T0'
  annotate('text', x = 0, y = c(250, 500, 750, 1000), label = c('250', '500', '750', '1000')) +
  #annotate('segment', x = 0, xend = 0, y = 0, yend =1000)+
  scale_y_continuous(limits = c(0, 1015))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

#richness_time <- grid.arrange(t0,t2,t3,t4, nrow = 1)

richness_time <- multi_panel_figure(columns = 2, rows = 2, width = 188, height = 188, 
                                                              column_spacing = 0, unit = 'mm',
                                                              row_spacing = 0,
                                                              panel_label_type = 'upper-alpha')


richness_time %<>%
  fill_panel(t0, column = 1, row = 1) %<>%
  fill_panel(t2, column = 2, row = 1) %<>%
  fill_panel(t3, column = 1, row = 2) %>%
  fill_panel(t4, column = 2, row = 2)

ggsave(filename = 'richnesss_time.pdf', plot = richness_time, 
       path = 'results/figures/corrected_asv_num/',
       width = 200, height = 200, units = 'mm')


###CONSTRUCTING PANNELS FOR THE ANALYSIS COMMON AND EXCLUSIVE ASVS------
####Absolut values
library(multipanelfigure)
common_exclusive_diversity_absolut_asvs <- multi_panel_figure(columns = 2, rows = 2, width = 188, height = 188, 
                                                     column_spacing = 0, unit = 'mm',
                                                     row_spacing = 0,
                                                     panel_label_type = 'upper-alpha')


common_exclusive_diversity_absolut_asvs %<>%
  fill_panel(t0, column = 1, row = 1) %<>%
  fill_panel(total_growing_asvs_plot_absolut, column = 2, row = 1) %<>%
  fill_panel(exclusive_asvs_absolut, column = 1, row = 2) %>%
  fill_panel(common_asvs_1gr_absolut, column = 2, row = 2)


##CAMBIAR EL NOM ABANS DE GUARDAR!!!!
ggsave('common_exclusive_richness_absolut_asv.pdf', common_exclusive_diversity_absolut_asvs,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 210,
       height = 200,
       units = 'mm')


####Relative values----
common_exclusive_diversity_relative_asvs <- multi_panel_figure(columns = 2, rows = 2, width = 188, height = 188, 
                                                              column_spacing = 0, unit = 'mm',
                                                              row_spacing = 0,
                                                              panel_label_type = 'upper-alpha')



common_exclusive_diversity_relative_asvs %<>%
  fill_panel(t0, column = 1, row = 1) %<>%
  fill_panel(total_growing_asvs_plot_relative, column = 2, row = 1) %<>%
  fill_panel(exclusive_asvs_relative_plot, column = 1, row = 2) %>%
  fill_panel(common_asvs_1gr_relative, column = 2, row = 2)


##CAMBIAR EL NOM ABANS DE GUARDAR!!!!
ggsave('common_exclusive_richness_relative_asv.pdf', common_exclusive_diversity_relative_asvs,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 210,
       height = 200,
       units = 'mm')

##SUMMARY PLOT (DISCARDED)
common_asv_sum %>%
  select(from, to, num_common_asvs) %>%
  pivot_longer(cols = c(from, to)) %>%
  filter(name != 'to')

exclusive_asvs_relative %>%
  left_join(total_growing_asvs, by = c('seas_treat', 'Season', 'Treatment')) %>%
  #left_join(common_asv_sum_l, by = c('seas_treat' = 'from')) %>%
  pivot_longer(cols = c('num_growing_asv', 'number_exclusive_asvs', 'num_responding_asv'), 
               names_to = 'category', values_to = 'number') %>%
  ggplot(aes(seas_treat, number))+
  geom_col(aes(fill = category))+
  geom_point(y = 700, aes(color = Season), size = 15, alpha = 0.5)+
  geom_text(aes(label = (Treatment), y = (700)), color = 'black')+
  # geom_col(aes(seas_treat, num_responding_asv))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_color_manual(values = palette_seasons_4)+
  scale_fill_manual(values = c('#08021f', '#E9AF25', '#0F2877'))+
 # annotate('text', x = 0, y = c(100, 200, 300), label = c('100', '200', '300')) +
  #scale_y_continuous(limits = c(0, 300))+
  theme_bw()+
  theme(legend.position = 'right',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))
  

  ##SanketNetwork-----------
# library(networkD3)
#  ##entre seasons i treatments
# reg_all_slopes_chosen_silva_tax_filt <-  
#   reg_all_slopes_chosen_silva_tax %>%
#   filter(slope_chosen_days > 2 &
#            pvalue_slope_chosen < 0.05) %>%
#   #mutate(season_treatment = paste(season,'_',treatment)) %>%
#   mutate(season_treatment = paste(season,'',treatment)) %>%
#   select(season_treatment, asv_num, slope_chosen_days) #phylum_f, class_f, order_f, family_f
# 
# 
# 
# colnames(reg_all_slopes_chosen_silva_tax_filt) <- c("source", "target", "value")
# reg_all_slopes_chosen_silva_tax_filt$target <- paste(reg_all_slopes_chosen_silva_tax_filt$target, " ", sep="")
# # From these flows we need to create a node data frame: it lists every entities involved in the flow
# nodes <- data.frame(name=c(as.character(reg_all_slopes_chosen_silva_tax_filt$source), as.character(reg_all_slopes_chosen_silva_tax_filt$target)) %>% unique())
# 
# # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
# reg_all_slopes_chosen_silva_tax_filt$IDsource=match(reg_all_slopes_chosen_silva_tax_filt$source, nodes$name)-1 
# reg_all_slopes_chosen_silva_tax_filt$IDtarget=match(reg_all_slopes_chosen_silva_tax_filt$target, nodes$name)-1
# 
# # prepare colour scale
# ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'
# 
# # Make the Network
# sankeyNetwork(Links = reg_all_slopes_chosen_silva_tax_filt, Nodes = nodes,
#               Source = "IDsource", Target = "IDtarget",
#               Value = "value", NodeID = "name", 
#               sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20)

x <- c(0.0333, 0.047, 0.139, 0.032) 
xday <-  x*24
mean(xday)
