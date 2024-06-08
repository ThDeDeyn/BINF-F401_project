###
# P1: graph pathways with clusters
###


d.pathway =  c(replicate(8,"Creation of C4 and C2 activators"),
                replicate(4,"Initial triggering of complements"),
                replicate(1,"Complement cascade"),
                replicate(2,"TNFS binds their physiological receptors"),
                replicate(1,"TNFR2 on canonical NF-KB pathway")
)
d.reg_status = "up-regulated"
d.clusters = c( "G4_1", "G4_2", "G4_9", "G4_15", "G4_18", "G4_21", "G4_24", "G4_30",
               "G4_2", "G4_15", "G4_18", "G4_24",
               "G4_2",
               "G4_15",  "G4_30",
               "G4_15")


data.path = data.frame(cluster = d.clusters, reg_status  = d.reg_status,
                        pathway = d.pathway, number = 1)


ggplot(data.path, aes(x= fct_inorder(pathway), y=number, fill = cluster)) +
  geom_bar(position = "fill", stat = "identity",color='black',width=0.9) +
  ggtitle("Up-regulated pathways and related clusters") + 
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(plot.title = element_text(size = 15, face = "bold")) 


###
# P2: most found genes among 320
###

candidates = rownames(top10iN32)
unique(candidates)
counts = table(candidates)
counts = as.data.frame(counts)
counts = counts[order(-counts$Freq),]

write.csv(counts, file = paste0(path_out, "occurence_in_top10.csv"))
