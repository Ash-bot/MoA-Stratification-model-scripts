#Sanity check!!!


setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Identification of a New Chemical Class of Antimalarials ACT213615/Series Matrix File/')

Art1= read.csv(file = 'ACTderiv processed probes.csv',header = T, sep = ';')

Art1[Art1 == "NaN"] = NA
Art1[Art1 == "null"] = NA
sum(is.na(Art1))
# [1] 2784
is.data.frame(Art1)


Art2 = data.frame(Gene.ID = Art1$ID_REF,apply(Art1[,2:11],2,as.numeric))
is.numeric(Art2$ACT.213615_TP1)
# Art2= na.omit(Art2)
l= length(unique(Art2$Gene.ID))
print(l)
# [1] 5462
Art2 = aggregate(.~Gene.ID, data = Art2, FUN = mean, na.rm = T, na.action = na.pass)


# write.table(collapsed_array, "collapsed_Hu_all compounds.txt", sep = "\t", row.names = F)
Art2[Art2 == "NaN"] = NA
sum(is.na(Art2))
# [1] 804
Art3 = data.frame(row.names = Art2$Gene.ID,
                  ArtTP1= Art2$ACT.213615_TP1,
                  ArtTP2= Art2$ACT.213615_TP2,
                  ArtTP3= Art2$ACT.213615_TP3,
                  ArtTP4= Art2$ACT.213615_TP4,
                  ArtTP5= Art2$ACT.213615_TP5,
                  ArtTP1c= Art2$DMSO_TP1,
                  ArtTP2c= Art2$DMSO_TP2,
                  ArtTP3c= Art2$DMSO_TP3,
                  ArtTP4c= Art2$DMSO_TP4,
                  ArtTP5c= Art2$DMSO_TP5)
get = grep('MAL|PF', rownames(Art3))        # | for or                                                                                                                                                                                                                                                                                                               
Art_genes = Art3[get,]
setwd('C:/Users/ashle/Desktop/LAbbook(2019)/')
IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
Art_genes$new_ID = IDs$new_ID[match(rownames(Art_genes), IDs$old_ID)]
new_IDs_filtered = Art_genes[-(which(Art_genes$new_ID %in% NA)),] 
# removes old gene names that no longer exist because there is not a new name for it

new_IDs_filtered=as.data.frame(new_IDs_filtered)
new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
Art=new_IDs_filtered[,-1]                

sum(is.na(Art))
Art = na.omit(Art)


setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')
write.table(Art, 'Act_data.csv', row.names=TRUE, append=FALSE, sep = ';')

###

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Three different histone deacetylase inhibitors/Series Matrix File/')

HDA1= read.delim(file = 'HDAi probes processed.txt',header = T)

HDA1[HDA1 == "NaN"] = NA
HDA1[HDA1 == "null"] = NA
sum(is.na(HDA1))
# [1] 65260
is.data.frame(HDA1)

HDA2 = data.frame(Gene.ID = HDA1$ID_REF,apply(HDA1[,2:11],2,as.numeric))
is.numeric(HDA2$DMSO.control.2h.replicate.1)
l= length(unique(HDA2$Gene.ID))
print(l)
# [1] 5485
HDA2 = aggregate(.~Gene.ID, data = HDA2, FUN = mean, na.rm = T, na.action = na.pass)
HDA2[HDA2 == "NaN"] = NA
sum(is.na(HDA2))

# write.table(collapsed_array, "collapsed_Hu_all compounds.txt", sep = "\t", row.names = F)
HDA3 = data.frame(row.names = HDA2$Gene.ID,
                  HDATP1= HDA2$SAHA.treated.2h.replicate.1,
                  HDATP1.1= HDA2$SAHA.treated.2h..2h..replicate.1,
                  HDATP2= HDA2$TSA.treated.2h..replicate.1,
                  HDATP2.2= HDA2$TSA.treated.2h..2h..replicate.1,
                  HDATP3= HDA2$X2.ASA.9.treated.2h..replicate.1,
                  HDATP3.3= HDA2$X2.ASA.9.treated.2h..2h..replicate.1,
                  HDATP1c= HDA2$DMSO.control.2h.replicate.1,
                  HDATP1c2= HDA2$DMSO.control.2h..2h..replicate.1)
get = grep('MAL|PF', rownames(HDA3))        # | for or                                                                                                                                                                                                                                                                                                               
HDA_genes = HDA3[get,]
setwd('C:/Users/ashle/Desktop/LAbbook(2019)')
IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
HDA_genes$new_ID = IDs$new_ID[match(rownames(HDA_genes), IDs$old_ID)]
new_IDs_filtered = HDA_genes[-(which(HDA_genes$new_ID %in% NA)),] 
# removes old gene names that no longer exist because there is not a new name for it

new_IDs_filtered=as.data.frame(new_IDs_filtered)
new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
HDA=new_IDs_filtered[,-1] 
sum(is.na(HDA))
# [1] 2060
HDAi= na.omit(HDA)

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')
write.table(HDA, 'HDA_data.csv', row.names=TRUE, append=FALSE, sep = ';')
####

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Role of Calcium Signaling in the Transcriptional Regulation of the Apicoplast Genome of Plasmodium falciparum/')

Iono1= read.csv(file = 'Ionomycin probes processed.csv',header = T, sep = ';')

Iono1[Iono1 == "NaN"] = NA
Iono1[Iono1 == "null"] = NA

sum(is.na(Iono1))
# [1] 8231
is.data.frame(Iono1)

Iono2 = data.frame(Gene.ID = Iono1$ID_REF,apply(Iono1[,2:16],2,as.numeric))
is.numeric(Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._Ionomycin.Treatment_TP1_30mins)
l= length(unique(Iono2$Gene.ID))
print(l)
# [1] 5300
Iono2 = aggregate(.~Gene.ID, data = Iono2, FUN = mean, na.rm = T, na.action = na.pass)

# write.table(collapsed_array, "collapsed_Hu_all compounds.txt", sep = "\t", row.names = F)
Iono2[Iono2 == "NaN"] = NA
sum(is.na(Iono2))
# [1] 2866

Iono3 = data.frame(row.names = Iono2$Gene.ID,
                   IonoTP1= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._Ionomycin.Treatment_TP1_30mins,
                   IonoTP2= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._Ionomycin.Treatment_TP2_1hr,
                   IonoTP3= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._Ionomycin.Treatment_TP3_2hr,
                   IonoTP4= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._Ionomycin.Treatment_TP4_4hr,
                   IonoTP5= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._Ionomycin.Treatment_TP5_6hr,
                   IonoTP1c= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._control_TP1_30mins,
                   IonoTP2c= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._control_TP2_1hr,
                   IonoTP3c= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._control_TP3_2hr,
                   IonoTP4c= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._control_TP4_4hr,
                   IonoTP5c= Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._control_TP5_6hr)


get = grep('MAL|PF', rownames(Iono3))        # | for or                                                                                                                                                                                                                                                                                                               
Iono_genes = Iono3[get,]
setwd('C:/Users/ashle/Desktop/LAbbook(2019)/')
IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
Iono_genes$new_ID = IDs$new_ID[match(rownames(Iono_genes), IDs$old_ID)]
new_IDs_filtered = Iono_genes[-(which(Iono_genes$new_ID %in% NA)),] 
# removes old gene names that no longer exist because there is not a new name for it

new_IDs_filtered=as.data.frame(new_IDs_filtered)
new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
Iono=new_IDs_filtered[,-1]                
sum(is.na(Iono))
# [1] 1803

Iono = na.omit(Iono)

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')
write.table(Iono, 'Iono_data.csv', row.names=TRUE, append=FALSE, sep = ';')
####

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Gupta DK,  Gupta AP et al. DNA damage regulation/Series Matrix File/')
dnad1= read.csv(file = 'DNA damaging processed probes.csv',header = T, sep = ';')

dnad1[dnad1 == "NaN"] = NA
dnad1[dnad1 == "null"] = NA
sum(is.na(dnad1))
# [1] 6881
is.data.frame(dnad1)

dnad2 = data.frame(Gene.ID = dnad1$ID_REF,apply(dnad1[,2:16],2,as.numeric))
is.numeric(dnad2$H4K8ac_control_6.hours_1)
l= length(unique(dnad2$Gene.ID))
print(l)
# [1] 5274
dnad2 = aggregate(.~Gene.ID, data = dnad2, FUN = mean, na.rm = T, na.action = na.pass)


# write.table(collapsed_array, "collapsed_Hu_all compounds.txt", sep = "\t", row.names = F)
dnad2[dnad2 == "NaN"] = NA
sum(is.na(dnad2))
# [1] 1700

dnad3 = data.frame(row.names = dnad2$Gene.ID,
                   dnadTP1.1= dnad2$H4K8ac_MMS_6.hours_1,
                   dnadTP1.2= dnad2$H4K8ac_MMS_6.hours_2,
                   dnadTP1.3= dnad2$H4K8ac_MMS_6.hours_3,
                   dnadTP1.1C= dnad2$H4K8ac_control_6.hours_1,
                   dnadTP1.2C= dnad2$H4K8ac_control_6.hours_2,
                   dnadTP1.3C= dnad2$H4K8ac_control_6.hours_3,
                   dnad2TP1.1= dnad2$H3K9ac_MMS_6.hour_1,
                   dnad2TP1.2= dnad2$H3K9ac_MMS_6.hour_2,
                   dnad2TP1.1C= dnad2$H3K9ac_control_6.hour_1,
                   dnad2TP1.2C= dnad2$H3K9ac_control_6.hour_1)


get = grep('MAL|PF', rownames(dnad3))        # | for or                                                                                                                                                                                                                                                                                                               
dnad_genes = dnad3[get,]
setwd('C:/Users/ashle/Desktop/LAbbook(2019)/')
IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
dnad_genes$new_ID = IDs$new_ID[match(rownames(dnad_genes), IDs$old_ID)]
new_IDs_filtered = dnad_genes[-(which(dnad_genes$new_ID %in% NA)),] 
# removes old gene names that no longer exist because there is not a new name for it

new_IDs_filtered=as.data.frame(new_IDs_filtered)
new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
dnad=new_IDs_filtered[,-1]                
sum(is.na(dnad))
# [1] 1803

dnad = na.omit(dnad)
setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')
write.table(dnad, 'MMS_data.csv', row.names=TRUE, append=FALSE, sep = ';')

#Replacing microarray ID numbers with actual gene ID

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Plasmodium falciparum asexual parasites and late gametocytes MMV390048 or MMV642943/Series Matrix File/')

MMV1= read.csv(file = 'MMV compounds need probes.csv',header = T, sep = ';')

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Plasmodium falciparum asexual parasites and late gametocytes MMV390048 or MMV642943/')

Genes= read.csv( file = 'GeneID.csv',header = T, sep = ';')
Genes= Genes[-2]

MMV1$Gene_ID = Genes$ORF[match(MMV1$ID_REF, Genes$ID)]
MMV2= MMV1[,c(12,2:11)]

MMV2[MMV2 == "NaN"] = NA
MMV2[MMV2 == "null"] = NA
sum(is.na(MMV2))
# [1] 34412

is.data.frame(MMV2)

MMV3 = data.frame(Gene.ID = MMV2$Gene_ID,apply(MMV2[,2:11],2,as.numeric))
is.numeric(MMV2$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV642943.and.harvested.24.h.post.expt.initiation..A24.T943.)
l= length(unique(MMV3$Gene.ID))
print(l)
# [1] 5643
MMV3 = aggregate(.~Gene.ID, data = MMV3, FUN = mean, na.rm = T, na.action = na.pass)

MMV3[MMV3 == "NaN"] = NA
MMV3[MMV3 == "null"] = NA
sum(is.na(MMV3))
# [1] 10571

# MMV4 = data.frame(row.names = MMV3$Gene.ID,
#                   MMV_43TP1= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV642943.and.harvested.24.h.post.expt.initiation..A24.T943.,
#                   MMV_43TP2= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV642943.and.harvested.48.h.post.expt.initiation..A48.T943.,
#                   MMV_48TP1= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV390048.and.harvested.48.h.post.expt.initiation..A24.T048.,
#                   MMV_48TP2= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV390048.and.harvested.24.h.post.expt.initiation..A48.T048.,
#                   MMV_TP1c= MMV3$NF54.asexual.ring.stage.parasites.remained.untreated.and.harvested.24.h.post.expt.initiation..A24.UT.,
#                   MMV_TP2c= MMV3$NF54.asexual.ring.stage.parasites.remained.untreated.and.harvested.48.h.post.expt.initiation..A48..UT.,
#                   MMV_43TP1g= MMV3$NF54.stage.IV.V.gametocytes.treated.with.10.x.IC50.of.MMV642943.and.harvested.24.h.post.expt.initiation..LG..48h..T943.,
#                   MMV_48TP1g= MMV3$NF54.stage.IV.V.gametocytes.treated.with.10.x.IC50.of.MMV390048.and.harvested.48.h.post.expt.initiation..LG..48h..T048.,
#                   MMV_48TPgc= MMV3$NF54.stage.IV.V.gametocytes.remained.untreated.and.were.harvested.48.h.post.expt.initiation..LG..48h..UT048.,
#                   MMV_43TPgc= MMV3$NF54.stage.IV.V.gametocytes.remained.untreated.and.were.harvested.48.h.post.expt.initiation..LG..48h..UT943.)

MMV4 = data.frame(row.names = MMV3$Gene.ID,
                  MMV_43TP1= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV642943.and.harvested.24.h.post.expt.initiation..A24.T943.,
                  MMV_43TP2= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV642943.and.harvested.48.h.post.expt.initiation..A48.T943.,
                  MMV_48TP1= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV390048.and.harvested.48.h.post.expt.initiation..A24.T048.,
                  MMV_48TP2= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV390048.and.harvested.24.h.post.expt.initiation..A48.T048.,
                  MMV_TP1c= MMV3$NF54.asexual.ring.stage.parasites.remained.untreated.and.harvested.24.h.post.expt.initiation..A24.UT.,
                  MMV_TP2c= MMV3$NF54.asexual.ring.stage.parasites.remained.untreated.and.harvested.48.h.post.expt.initiation..A48..UT.)

get = grep('MAL|PF', rownames(MMV4))        # | for or                                                                                                                                                                                                                                                                                                               
MMV_genes = MMV4[get,]
setwd('C:/Users/ashle/Desktop/LAbbook(2019)/')
IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
MMV_genes$new_ID = IDs$new_ID[match(rownames(MMV_genes), IDs$old_ID)]
new_IDs_filtered = MMV_genes[-(which(MMV_genes$new_ID %in% NA)),]


new_IDs_filtered=as.data.frame(new_IDs_filtered)
new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
MMV=new_IDs_filtered[,-1]  
sum(is.na(MMV))
# [1] 1811

MMV= na.omit(MMV)
sum(is.na(MMV))
# [1] 0

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')
write.table(MMV, 'MMV_data.csv', row.names=TRUE, append=FALSE, sep = ';')

###Cyclohexamine
# 
# setwd('/Users/ashle/Desktop/Additional datasets/Accepted/Plasmodium falciparum treated with cyclohexylamine/Series Matrix File/')
# 
# Cyc1= read.csv(file = 'Cyclohexamine processed probes.csv',header = T, sep = ';')
# 
# setwd('/Users/ashle/Desktop/Additional datasets/Accepted/Plasmodium falciparum treated with cyclohexylamine/Series Matrix File/')
# 
# Genes= read.delim( file = 'GeneID.txt',header = T)
# Genes= Genes[-3]
# 
# Genes[Genes == "N/A"] = NA
# Genes[Genes == "NULL"] = NA
# 
# 
# Cyc1$Gene_ID = Genes$PlasmoDB.ID[match(Cyc1$ID_REF, Genes$ID)]
# Cyc3= Cyc1[,c(20,2:19)]
# sum(is.na(Cyc3))
# # [1] 3982
# 
# Cyc3[Cyc3 == "NaN"] = NA
# Cyc3[Cyc3 == "Null"] = NA
# sum(is.na(Cyc3))
# # [1] 2394
# 
# is.data.frame(Cyc3)
# Cyc4 = data.frame(Gene.ID = Cyc3$Gene_ID,apply(Cyc3[,2:19],2,as.numeric))
# is.numeric(Cyc4$ctrl.18.1.tech.1)
# l= length(unique(Cyc4$Gene.ID))
# print(l)
# # [1] 2695
# Cyc4 = aggregate(.~Gene.ID, data = Cyc4, FUN = mean, na.rm = T, na.action = na.pass)
# 
# Cyc4[Cyc4 == "NaN"] = NA
# Cyc4[Cyc4 == "null"] = NA
# sum(is.na(Cyc4))
# # [1] 1746
# 
# Cyc5 = data.frame(row.names = Cyc4$Gene.ID,
#                   CycTP1= Cyc4$cyclo.18.1.tech.1,
#                   CycTP1.2=Cyc4$cyclo.18.1.tech.2,
#                   CycTP2= Cyc4$cyclo.18.2,
#                   CycTP3= Cyc4$cyclo.25.1.tech.1,
#                   CycTP3.2= Cyc4$cyclo.25.1.tech.2,
#                   CycTP4= Cyc4$cyclo.25.2,
#                   CycTP5= Cyc4$cyclo.30.2.tech.1,
#                   CycTP5.2= Cyc4$cyclo.30.2.tech.2,
#                   CycTP6= Cyc4$cyclo.30.1,
#                   CycTP1c= Cyc4$ctrl.18.1.tech.1,
#                   CycTP1.2c=Cyc4$ctrl.18.1.tech.2,
#                   CycTP2c= Cyc4$ctrl.18.2,
#                   CycTP3c= Cyc4$ctrl.25.2.tech.1,
#                   CycTP3.2c= Cyc4$ctrl.25.2.tech.2,
#                   CycTP4c= Cyc4$ctrl.25.1,
#                   CycTP5c= Cyc4$ctrl.30.1.tech.1,
#                   CycTP5.2c= Cyc4$ctrl.30.1.tech.2,
#                   CycTP6c= Cyc4$ctrl.30.2)
# 
# get = grep('MAL|PF', rownames(Cyc5))        # | for or                                                                                                                                                                                                                                                                                                               
# Cyc_genes = Cyc5[get,]
# setwd('C:/Users/ashle/Desktop/LAbbook')
# IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
# Cyc_genes$new_ID = IDs$new_ID[match(rownames(Cyc_genes), IDs$old_ID)]
# new_IDs_filtered = Cyc_genes[-(which(Cyc_genes$new_ID %in% NA)),]
# 
# 
# new_IDs_filtered=as.data.frame(new_IDs_filtered)
# new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
# rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
# Cyc=new_IDs_filtered[,-1]  
# sum(is.na(Cyc))
# # [1] 1686
# Cyc= na.omit(Cyc)
###Cyclohexamine

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Plasmodium falciparum treated with cyclohexylamine/Series Matrix File/')

Cyc1= read.csv(file = 'Cyclohexamine processed probes.csv',header = T, sep = ';')

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Plasmodium falciparum treated with cyclohexylamine/Series Matrix File/')

Genes= read.delim( file = 'GeneID.txt',header = T)
Genes= Genes[-3]

Genes[Genes == "N/A"] = NA
Genes[Genes == "NULL"] = NA
Genes= na.omit(Genes)

Cyc1$Gene_ID = Genes$PlasmoDB.ID[match(Cyc1$ID_REF, Genes$ID)]
Cyc4= Cyc1[,c(20,2:19)]
sum(is.na(Cyc4))
# [1] 3982

Cyc4[Cyc4 == "NaN"] = NA
Cyc4[Cyc4 == "Null"] = NA
sum(is.na(Cyc3))
# [1] 13108
is.data.frame(Cyc4)
Cyc4 = data.frame(Gene.ID = Cyc4$Gene_ID,apply(Cyc4[,2:19],2,as.numeric))
is.numeric(Cyc4$ctrl.18.1.tech.1)
Cyc4= na.omit(Cyc4)

l= length(unique(Cyc4$Gene.ID))
print(l)
# [1] 2695
Cyc4 = aggregate(.~Gene.ID, data = Cyc4, FUN = mean, na.rm = T, na.action = na.pass)

Cyc4[Cyc4 == "NaN"] = NA
Cyc4[Cyc4 == "null"] = NA
sum(is.na(Cyc4))
# [1] 1746


Cyc5 = data.frame(row.names = Cyc4$Gene.ID,
                  CycTP1= Cyc4$cyclo.18.1.tech.1,
                  CycTP1.2=Cyc4$cyclo.18.1.tech.2,
                  CycTP2= Cyc4$cyclo.18.2,
                  CycTP3= Cyc4$cyclo.25.1.tech.1,
                  CycTP3.2= Cyc4$cyclo.25.1.tech.2,
                  CycTP4= Cyc4$cyclo.25.2,
                  CycTP5= Cyc4$cyclo.30.2.tech.1,
                  CycTP5.2= Cyc4$cyclo.30.2.tech.2,
                  CycTP6= Cyc4$cyclo.30.1,
                  CycTP1c= Cyc4$ctrl.18.1.tech.1,
                  CycTP1.2c=Cyc4$ctrl.18.1.tech.2,
                  CycTP2c= Cyc4$ctrl.18.2,
                  CycTP3c= Cyc4$ctrl.25.2.tech.1,
                  CycTP3.2c= Cyc4$ctrl.25.2.tech.2,
                  CycTP4c= Cyc4$ctrl.25.1,
                  CycTP5c= Cyc4$ctrl.30.1.tech.1,
                  CycTP5.2c= Cyc4$ctrl.30.1.tech.2,
                  CycTP6c= Cyc4$ctrl.30.2)
get = grep('MAL|PF', rownames(Cyc5))        # | for or
Cyc_genes = Cyc5[get,]
setwd('C:/Users/ashle/Desktop/LAbbook(2019)/')
IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
Cyc_genes$new_ID = IDs$new_ID[match(rownames(Cyc_genes), IDs$old_ID)]
new_IDs_filtered = Cyc_genes[-(which(Cyc_genes$new_ID %in% NA)),]


new_IDs_filtered=as.data.frame(new_IDs_filtered)
new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
Cyc=new_IDs_filtered[,-1]


sum(is.na(Cyc5))
# [1] 1746
Cyc= na.omit(Cyc)

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')
write.table(Cyc, 'Cyc_data.csv',row.names=TRUE, append=FALSE, sep = ';')

###DHB

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/In vitro antiplasmodial activity of Dicoma anomala subsp gerrardii (Asteraceae/Series Matrix File/')

DHB1= read.csv(file = 'GSE29874_series_matrix.csv',header = T, sep = ';')

Genes= read.delim( file = 'GeneID.csv',header = T, sep = ';')
Genes= Genes[-3]

Genes[Genes == "N/A"] = NA
Genes[Genes == "NULL"] = NA

DHB1$Gene_ID = Genes$PlasmoDB.ID[match(DHB1$ID_REF, Genes$ID)]
DHB2= DHB1[,c(21,2:20)]

DHB2[DHB2 == "NaN"] = NA
DHB2[DHB2 == "Null"] = NA
sum(is.na(DHB2))
# [1] 4085

DHB3= DHB2



is.data.frame(DHB3)
DHB4 = data.frame(Gene.ID = DHB3$Gene_ID,apply(DHB3[,2:20],2,as.numeric))
is.numeric(DHB4$DMSO.treated.control.technical.replicate1.after.2.hours.drug.exposure)
l= length(unique(DHB4$Gene.ID))
print(l)
# [1] 2696
DHB4 = aggregate(.~Gene.ID, data = DHB4, FUN = mean, na.rm = T, na.action = na.pass)

DHB4[DHB4 == "NaN"] = NA
DHB4[DHB4 == "null"] = NA
sum(is.na(DHB4))
# [1] 0

DHB5 = data.frame(row.names = DHB4$Gene.ID,
                  DHB4TP1.1= DHB4$Dehydrobrachylaenolide.treated.biological.replicate1.after.2.hours.drug.exposure,
                  DHB4TP1.2= DHB4$Dehydrobrachylaenolide.treated.biological.replicate2.after.2.hours.drug.exposure,
                  DHB4TP1.1t= DHB4$Dehydrobrachylaenolide.treated.technical.replicate1.after.2.hours.drug.exposure,
                  DHB4TP2.1= DHB4$Dehydrobrachylaenolide.treated.biological.replicate1.after.6.hours.drug.exposure,
                  DHB4TP2.2= DHB4$Dehydrobrachylaenolide.treated.biological.replicate2.after.6.hours.drug.exposure,
                  DHB4TP2t= DHB4$Dehydrobrachylaenolide.treated.technical.replicate.after.6.hours.drug.exposure,
                  DHB4TP3.1= DHB4$Dehydrobrachylaenolide.treated.biological.replicate1.after.12.hours.drug.exposure,
                  DHB4TP3.2= DHB4$Dehydrobrachylaenolide.treated.biological.replicate2.after.12.hours.drug.exposure,
                  DHB4TP3t= DHB4$Dehydrobrachylaenolide.treated..technical.replicate.after.12.hours.drug.exposure,
                  DHB4TP1.1c= DHB4$DMSO.treated.control.biological.replicate1.after.2.hours.drug.exposure,
                  DHB4TP1.2c= DHB4$DMSO.treated.control.biological.replicate2.after.2.hours.drug.exposure,
                  DHB4TP1.1tc= DHB4$DMSO.treated.control.technical.replicate1.after.2.hours.drug.exposure,
                  DHB4TP1.2tc= DHB4$DMSO.treated.control.technical.replicate2.after.2.hours.drug.exposure,
                  DHB4TP2.1c= DHB4$DMSO.treated.control.biological.replicate1.after.6.hours.drug.exposure,
                  DHB4TP2.2c= DHB4$DMSO.treated.control.biological.replicate2.after.6.hours.drug.exposure,
                  DHB4TP2tc= DHB4$DMSO.treated.control.technical.replicate.after.6.hours.drug.exposure,
                  DHB4TP3.1c= DHB4$DMSO.treated.control.biological.replicate1.after.12.hours.drug.exposure,
                  DHB4TP3.2c= DHB4$DMSO.treated.control.biological.replicate2.after.12.hours.drug.exposure,
                  DHB4TP3tc= DHB4$DMSO.treated.control.technical.replicate.after.12.hours.drug.exposure)

get = grep('MAL|PF', rownames(DHB5))        # | for or                                                                                                                                                                                                                                                                                                               
DHB_genes = DHB5[get,]
setwd('C:/Users/ashle/Desktop/LAbbook(2019)/')
IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
DHB_genes$new_ID = IDs$new_ID[match(rownames(DHB_genes), IDs$old_ID)]
new_IDs_filtered = DHB_genes[-(which(DHB_genes$new_ID %in% NA)),]


new_IDs_filtered=as.data.frame(new_IDs_filtered)
new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
DHB=new_IDs_filtered[,-1]  
sum(is.na(DHB))
# [1] 0
DHB= na.omit(DHB)

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')
write.table(DHB, 'DHB_data.csv', row.names=TRUE, append=FALSE, sep = ';')
###DFMO

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Functional genomics investigation of PfAdoMetDCODC coinhibition of Plasmodium falciparum/Series Matrix File/')

DFMO1= read.csv(file = 'DMFOseriesmatrix.csv',header = T, sep = ';')

Genes= read.delim( file = 'GeneID.csv',header = T, sep = ';')
Genes= Genes[-3]

Genes[Genes == "N/A"] = NA
Genes[Genes == "NULL"] = NA

DFMO1$Gene_ID = Genes$PLASMODB.ID[match(DFMO1$ID_REF, Genes$ID)]
DFMO2= DFMO1[,c(22,2:21)]

DFMO2[DFMO2 == "NaN"] = NA
DFMO2[DFMO2 == "Null"] = NA
# DFMO3= na.omit(DFMO2)
sum(is.na(DFMO2))
# [1] 37249
DFMO3= DFMO2

is.data.frame(DFMO3)
DFMO4 = data.frame(Gene.ID = DFMO3$Gene_ID,apply(DFMO3[,2:21],2,as.numeric))
is.numeric(DFMO4$X107_Pf_UT_t19_repA)
l= length(unique(DFMO4$Gene.ID))
print(l)
# [1] 5012
DFMO4 = aggregate(.~Gene.ID, data = DFMO4, FUN = mean, na.rm = T, na.action = na.pass)

DFMO4[DFMO4 == "NaN"] = NA
DFMO4[DFMO4 == "null"] = NA
sum(is.na(DFMO4))
# [1] 17237

DFMO5 = data.frame(row.names = DFMO4$Gene.ID,
                   DFMOT1.2= DFMO4$X127_Pf_T_t19_repB,
                   DFMOT1.2t= DFMO4$X116_Pf_T_t19_repB_tech,
                   DFMOT1.2c= DFMO4$X124_Pf_UT_t19_repB,
                   DFMOT1.1c= DFMO4$X107_Pf_UT_t19_repA,
                   DFMOT1.1tc= DFMO4$X095_Pf_UT_t19_repA_tech,
                   DFMOT2.1= DFMO4$X111_Pf_T_t27_repA,
                   DFMOT2.1t= DFMO4$X098_Pf_T_t27_repA_tech,
                   DFMOT2.2= DFMO4$X128_Pf_T_t27_repB,
                   DFMOT2.2t= DFMO4$X102_Pf_T_t27_repB_tech,
                   DFMOT2.1c= DFMO4$X108_Pf_UT_t27_repA,
                   DFMOT2.1tc= DFMO4$X091_Pf_UT_t27_repA_tech,
                   DFMOT2.2c= DFMO4$X125_Pf_UT_t27_repB,
                   DFMOT3.1= DFMO4$X112_Pf_T_t34_repA,
                   DFMOT3.1t= DFMO4$X090_Pf_T_t34_repA_tech,
                   DFMOT3.2= DFMO4$X129_Pf_T_t34_repB,
                   DFMOT3.2t= DFMO4$X094_Pf_T_t34_repB_tech,
                   DFMOT3.1c= DFMO4$X109_Pf_UT_t34_repA,
                   DFMOT3.1tc= DFMO4$X097_Pf_UT_t34_repA_tech,
                   DFMOT3.2cc= DFMO4$X126_Pf_UT_t34_repB,
                   DFMOT3.2c= DFMO4$X129_Pf_T_t34_repB,
                   DFMOT3.2tc= DFMO4$X089_Pf_UT_t34_repB_tech)

get = grep('MAL|PF', rownames(DFMO5))        # | for or                                                                                                                                                                                                                                                                                                               
DFMO_genes = DFMO5[get,]
setwd('C:/Users/ashle/Desktop/LAbbook(2019)/')
IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
DFMO_genes$new_ID = IDs$new_ID[match(rownames(DFMO_genes), IDs$old_ID)]
new_IDs_filtered = DFMO_genes[-(which(DFMO_genes$new_ID %in% NA)),]


new_IDs_filtered=as.data.frame(new_IDs_filtered)
new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
DFMO=new_IDs_filtered[,-1]  
sum(is.na(DFMO))
# [1] 17125
DFMO= na.omit(DFMO)

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')
write.table(DFMO, 'DFMO_data.csv', row.names=TRUE, append=FALSE, sep = ';')

####HU
setwd('C:/Users/ashle/Desktop/LAbbook(2019)/')

hu_Mvalue_genes= read.delim('hu_set.txt')

selected = data.frame(row.names = hu_Mvalue_genes$GeneID,
                      ML_7_H1= hu_Mvalue_genes$ML7_exp7_TP1, 
                      ML_7_H2= hu_Mvalue_genes$ML7_exp7_TP2,
                      ML_7_H4= hu_Mvalue_genes$ML7_exp7_TP3,
                      ML_7_H6= hu_Mvalue_genes$ML7_exp7_TP4,
                      ML_7_H8= hu_Mvalue_genes$ML7_exp7_TP5,
                      ML_7_H10= hu_Mvalue_genes$ML7_exp7_TP6,
                      ML_7_H12= hu_Mvalue_genes$ML7_exp7_TP7,
                      ML_7_H1c=hu_Mvalue_genes$Control_3_TP1,
                      ML_7_H2c=hu_Mvalue_genes$Control_3_TP2,
                      ML_7_H4c=hu_Mvalue_genes$Control_3_TP3,
                      ML_7_H6c=hu_Mvalue_genes$Control_3_TP4,
                      ML_7_H8c=hu_Mvalue_genes$Control_3_TP5,
                      ML_7_H10c=hu_Mvalue_genes$Control_3_TP6,
                      ML_7_H12c=hu_Mvalue_genes$Control_3_TP7,
                      W7_H1= hu_Mvalue_genes$W7_exp8_TP1,
                      W7_H2= hu_Mvalue_genes$W7_exp8_TP2,
                      W7_H4= hu_Mvalue_genes$W7_exp8_TP3,
                      W7_H6= hu_Mvalue_genes$W7_exp8_TP4,
                      W7_H8= hu_Mvalue_genes$W7_exp8_TP5,
                      W7_H10= hu_Mvalue_genes$W7_exp8_TP6,
                      W7_H12= hu_Mvalue_genes$W7_exp8_TP7,
                      W7_H1c=hu_Mvalue_genes$Control_3_TP1,
                      W7_H2c=hu_Mvalue_genes$Control_3_TP2,
                      W7_H4c=hu_Mvalue_genes$Control_3_TP3,
                      W7_H6c=hu_Mvalue_genes$Control_3_TP4,
                      W7_H8c=hu_Mvalue_genes$Control_3_TP5,
                      W7_H10c=hu_Mvalue_genes$Control_3_TP6,
                      W7_H12c=hu_Mvalue_genes$Control_3_TP7,
                      Staurosporine_H1= hu_Mvalue_genes$Staurosporine_exp14_TP1,
                      Staurosporine_H2= hu_Mvalue_genes$Staurosporine_exp14_TP2,
                      Staurosporine_H4= hu_Mvalue_genes$Staurosporine_exp14_TP3,
                      Staurosporine_H6= hu_Mvalue_genes$Staurosporine_exp14_TP4,
                      Staurosporine_H8= hu_Mvalue_genes$Staurosporine_exp14_TP5,
                      Staurosporine_H10= hu_Mvalue_genes$Staurosporine_exp14_TP6,
                      Staurosporine_H12= hu_Mvalue_genes$Staurosporine_exp14_TP7,
                      Staurosporine_H1c=hu_Mvalue_genes$Control_3_TP1,
                      Staurosporine_H2c=hu_Mvalue_genes$Control_3_TP2,
                      Staurosporine_H4c=hu_Mvalue_genes$Control_3_TP3,
                      Staurosporine_H6c=hu_Mvalue_genes$Control_3_TP4,
                      Staurosporine_H8c=hu_Mvalue_genes$Control_3_TP5,
                      Staurosporine_H10c=hu_Mvalue_genes$Control_3_TP6,
                      Staurosporine_H12c=hu_Mvalue_genes$Control_3_TP7,
                      Cyclosporine_H1= hu_Mvalue_genes$CyclosporineA_exp2_TP1,
                      Cyclosporine_H2= hu_Mvalue_genes$CyclosporineA_exp2_TP2,
                      Cyclosporine_H4= hu_Mvalue_genes$CyclosporineA_exp2_TP3,
                      Cyclosporine_H6= hu_Mvalue_genes$CyclosporineA_exp2_TP4,
                      Cyclosporine_H8= hu_Mvalue_genes$CyclosporineA_exp2_TP5,
                      Cyclosporine_H10= hu_Mvalue_genes$CyclosporineA_exp2_TP6,
                      Cyclosporine_H12= hu_Mvalue_genes$CyclosporineA_exp2_TP7,
                      Cyclosporine_H14= hu_Mvalue_genes$CyclosporineA_exp2_TP8,
                      Cyclosporine_H1c= hu_Mvalue_genes$Control_1_TP1,
                      Cyclosporine_H2c= hu_Mvalue_genes$Control_1_TP2,
                      Cyclosporine_H4c= hu_Mvalue_genes$Control_1_TP3,
                      Cyclosporine_H6c= hu_Mvalue_genes$Control_1_TP4,
                      Cyclosporine_H8c= hu_Mvalue_genes$Control_1_TP5,
                      Cyclosporine_H10c= hu_Mvalue_genes$Control_1_TP6,
                      Cyclosporine_H12c= hu_Mvalue_genes$Control_1_TP7,
                      Cyclosporine_H14c= hu_Mvalue_genes$Control_1_TP8,
                      Colchine_H1= hu_Mvalue_genes$Colchicine_exp5_TP1,
                      Colchine_H2= hu_Mvalue_genes$Colchicine_exp5_TP2,
                      Colchine_H4= hu_Mvalue_genes$Colchicine_exp5_TP3,
                      Colchine_H6= hu_Mvalue_genes$Colchicine_exp5_TP4,
                      Colchine_H8= hu_Mvalue_genes$Colchicine_exp5_TP5,
                      Colchine_H1c= hu_Mvalue_genes$Control_2_TP1,
                      Colchine_H2c= hu_Mvalue_genes$Control_2_TP2,
                      Colchine_H4c= hu_Mvalue_genes$Control_2_TP3,
                      Colchine_H6c= hu_Mvalue_genes$Control_2_TP4,
                      Colchine_H8c= hu_Mvalue_genes$Control_2_TP5,
                      PMSF_H1= hu_Mvalue_genes$PMSF_exp16_TP1,
                      PMSF_H2= hu_Mvalue_genes$PMSF_exp16_TP2,
                      PMSF_H4= hu_Mvalue_genes$PMSF_exp16_TP3,
                      PMSF_H6= hu_Mvalue_genes$PMSF_exp16_TP4,
                      PMSF_H8= hu_Mvalue_genes$PMSF_exp16_TP5,
                      PMSF_H1c= hu_Mvalue_genes$Control_5_TP1,
                      PMSF_H2c= hu_Mvalue_genes$Control_5_TP2,
                      PMSF_H4c= hu_Mvalue_genes$Control_5_TP3,
                      PMSF_H6c= hu_Mvalue_genes$Control_5_TP4,
                      PMSF_H8c= hu_Mvalue_genes$Control_5_TP5,
                      Leupeptine_H1= hu_Mvalue_genes$Leupeptine_exp24_TP1,
                      Leupeptine_H2= hu_Mvalue_genes$Leupeptine_exp24_TP2,
                      Leupeptine_H4= hu_Mvalue_genes$Leupeptine_exp24_TP3,
                      Leupeptine_H6= hu_Mvalue_genes$Leupeptine_exp24_TP4,
                      Leupeptine_H8= hu_Mvalue_genes$Leupeptine_exp24_TP5,
                      Leupeptine_H1c= hu_Mvalue_genes$Control_5_TP1,
                      Leupeptine_H2c= hu_Mvalue_genes$Control_5_TP2,
                      Leupeptine_H4c= hu_Mvalue_genes$Control_5_TP3,
                      Leupeptine_H6c= hu_Mvalue_genes$Control_5_TP4,
                      Leupeptine_H8c= hu_Mvalue_genes$Control_5_TP5,
                      Apicidin_H1= hu_Mvalue_genes$Apicidin.troph.5nM._exp17_TP1,
                      Apicidin_H2= hu_Mvalue_genes$Apicidin.troph.5nM._exp17_TP2,
                      Apicidin_H4= hu_Mvalue_genes$Apicidin.troph.5nM._exp17_TP3,
                      Apicidin_H6= hu_Mvalue_genes$Apicidin.troph.5nM._exp17_TP4,
                      Apicidin_H8= hu_Mvalue_genes$Apicidin.troph.5nM._exp17_TP5,
                      Apicidin_H10= hu_Mvalue_genes$Apicidin.troph.5nM._exp17_TP6,
                      Apicidin_H1c= hu_Mvalue_genes$Control_6_TP1,
                      Apicidin_H2c= hu_Mvalue_genes$Control_6_TP2,
                      Apicidin_H4c= hu_Mvalue_genes$Control_6_TP3,
                      Apicidin_H6c= hu_Mvalue_genes$Control_6_TP4,
                      Apicidin_H8c= hu_Mvalue_genes$Control_6_TP5,
                      Apicidin_H10c= hu_Mvalue_genes$Control_6_TP6,
                      Trichostatin_H1= hu_Mvalue_genes$TrichostatinA.IC50._epx21_TP1,
                      Trichostatin_H2= hu_Mvalue_genes$TrichostatinA.IC50._epx21_TP2,
                      Trichostatin_H4= hu_Mvalue_genes$TrichostatinA.IC50._epx21_TP3,
                      Trichostatin_H6= hu_Mvalue_genes$TrichostatinA.IC50._epx21_TP4,
                      Trichostatin_H8= hu_Mvalue_genes$TrichostatinA.IC50._epx21_TP5,
                      Trichostatin_H10= hu_Mvalue_genes$TrichostatinA.IC50._epx21_TP6,
                      Trichostatin_H1c= hu_Mvalue_genes$Control_8_TP1,
                      Trichostatin_H2c= hu_Mvalue_genes$Control_8_TP2,
                      Trichostatin_H4c= hu_Mvalue_genes$Control_8_TP3,
                      Trichostatin_H6c= hu_Mvalue_genes$Control_8_TP4,
                      Trichostatin_H8c= hu_Mvalue_genes$Control_8_TP5,
                      Trichostatin_H10c= hu_Mvalue_genes$Control_8_TP6,
                      Quinine_H1= hu_Mvalue_genes$Quinine_exp9_TP1,
                      Quinine_H2= hu_Mvalue_genes$Quinine_exp9_TP2,
                      Quinine_H4= hu_Mvalue_genes$Quinine_exp9_TP3,
                      Quinine_H6= hu_Mvalue_genes$Quinine_exp9_TP4,
                      Quinine_H8= hu_Mvalue_genes$Quinine_exp9_TP5,
                      Quinine_H10= hu_Mvalue_genes$Quinine_exp9_TP6,
                      Quinine_H1c= hu_Mvalue_genes$Conntrol_4_TP1,
                      Quinine_H2c= hu_Mvalue_genes$Conntrol_4_TP2,
                      Quinine_H4c= hu_Mvalue_genes$Conntrol_4_TP3,
                      Quinine_H6c= hu_Mvalue_genes$Conntrol_4_TP4,
                      Quinine_H8c= hu_Mvalue_genes$Conntrol_4_TP5,
                      Quinine_H10c= hu_Mvalue_genes$Conntrol_4_TP6,
                      Chloroquine_H1= hu_Mvalue_genes$Chloroquine.IC50._exp25_TP1,
                      Chloroquine_H2= hu_Mvalue_genes$Chloroquine.IC50._exp25_TP2,
                      Chloroquine_H4= hu_Mvalue_genes$Chloroquine.IC50._exp25_TP3,
                      Chloroquine_H6= hu_Mvalue_genes$Chloroquine.IC50._exp25_TP4,
                      Chloroquine_H8= hu_Mvalue_genes$Chloroquine.IC50._exp25_TP5,
                      Chloroquine_H10= hu_Mvalue_genes$Chloroquine.IC50._exp25_TP6,
                      Chloroquine_H1c= hu_Mvalue_genes$Control_9_TP1,
                      Chloroquine_H2c= hu_Mvalue_genes$Control_9_TP2,
                      Chloroquine_H4c= hu_Mvalue_genes$Control_9_TP3,
                      Chloroquine_H6c= hu_Mvalue_genes$Control_9_TP4,
                      Chloroquine_H8c= hu_Mvalue_genes$Control_9_TP5,
                      Chloroquine_H10c= hu_Mvalue_genes$Control_9_TP6,
                      Febrifugine_H1= hu_Mvalue_genes$Febrifugine_exp10_TP1,
                      Febrifugine_H2= hu_Mvalue_genes$Febrifugine_exp10_TP2,
                      Febrifugine_H4= hu_Mvalue_genes$Febrifugine_exp10_TP3,
                      Febrifugine_H6= hu_Mvalue_genes$Febrifugine_exp10_TP4,
                      Febrifugine_H8= hu_Mvalue_genes$Febrifugine_exp10_TP5,
                      Febrifugine_H10= hu_Mvalue_genes$Febrifugine_exp10_TP6,
                      Febrifugine_H1c= hu_Mvalue_genes$Conntrol_4_TP1,
                      Febrifugine_H2c= hu_Mvalue_genes$Conntrol_4_TP2,
                      Febrifugine_H4c= hu_Mvalue_genes$Conntrol_4_TP3,
                      Febrifugine_H6c= hu_Mvalue_genes$Conntrol_4_TP4,
                      Febrifugine_H8c= hu_Mvalue_genes$Conntrol_4_TP5,
                      Febrifugine_H10c= hu_Mvalue_genes$Conntrol_4_TP6,
                      Artemisinin_H1= hu_Mvalue_genes$Artemisinin_exp11_TP1,
                      Artemisinin_H2= hu_Mvalue_genes$Artemisinin_exp11_TP2,
                      Artemisinin_H4= hu_Mvalue_genes$Artemisinin_exp11_TP3,
                      Artemisinin_H6= hu_Mvalue_genes$Artemisinin_exp11_TP4,
                      Artemisinin_H8= hu_Mvalue_genes$Artemisinin_exp11_TP5,
                      Artemisinin_H10= hu_Mvalue_genes$Artemisinin_exp11_TP6,
                      Artemisinin_H1c= hu_Mvalue_genes$Conntrol_4_TP1,
                      Artemisinin_H2c= hu_Mvalue_genes$Conntrol_4_TP2,
                      Artemisinin_H4c= hu_Mvalue_genes$Conntrol_4_TP3,
                      Artemisinin_H6c= hu_Mvalue_genes$Conntrol_4_TP4,
                      Artemisinin_H8c= hu_Mvalue_genes$Conntrol_4_TP5,
                      Artemisinin_H10c= hu_Mvalue_genes$Conntrol_4_TP6)
# LogFC

selected[selected == "NaN"] = NA
get = grep('MAL|PF', rownames(selected))        # | for or                                                                                                                                                                                                                                                                                                               
selected_genes = selected[get,]
setwd('C:/Users/ashle/Desktop/LAbbook(2019)/')
IDs=read.delim("new_old_gene_ID.txt")  #2017 updated list
selected_genes$new_ID = IDs$new_ID[match(rownames(selected_genes), IDs$old_ID)]
new_IDs_filtered = selected_genes[-(which(selected_genes$new_ID %in% NA)),] 
# removes old gene names that no longer exist because there is not a new name for it

new_IDs_filtered=as.data.frame(new_IDs_filtered)
new_IDs_filtered = aggregate(.~new_ID, data = new_IDs_filtered, FUN = mean, na.action = na.pass)
rownames(new_IDs_filtered)=new_IDs_filtered$new_ID
new_IDs_filtered=new_IDs_filtered[,-1]  #remove column one -1

library(xlsx)
Hu= na.omit(new_IDs_filtered)
setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')
write.table(Hu, 'Hu_data.csv', row.names=TRUE, append=FALSE, sep = ';')


#####################################MERGE####################################


setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Filtered datasets/Redo/')

Iono= read.csv('Iono_data.csv', header = T, sep = ';')
Iono= na.omit(Iono)

Act= read.csv('Act_data.csv', header = T, sep = ';')
Act= na.omit(Act)

Cyc= read.csv('Cyc_data.csv', header = T, sep = ';')
Cyc= na.omit(Cyc)

DFMO= read.csv('DFMO_data.csv', header = T, sep = ';')
DFMO= na.omit(DFMO)

DHB= read.csv('DHB_data.csv', header = T, sep = ';')
DHB= na.omit(DHB)

HDA= read.csv('HDA_data.csv', header = T, sep = ';')
HDA= na.omit(HDA)

MMS= read.csv('MMS_data.csv', header = T, sep = ';')
MMS= na.omit(MMS)

MMV= read.csv('MMV_data.csv', header = T, sep = ';')
MMV= na.omit(MMV)

Hu= read.csv('Hu_data.csv', sep = ';') #NB go into file and add X by rownames
Hu= na.omit(Hu)

a =merge(Cyc,DFMO,by = 'X')
b= merge(a,DHB, by= 'X')
c= merge(b,MMV, by= 'X')
d= merge(c,Iono, by= 'X')
e= merge(d,DFMO, by= 'X')
f= merge(e,HDA, by= 'X')
g=merge(f,MMS, by= 'X')  
h=merge(g,Hu, by= 'X')
Final= merge(h,Act, by= 'X')

#Merge with datasets that form database


c= merge(DFMO,MMV, by= 'X')
d= merge(c,Iono, by= 'X')
f= merge(d,HDA, by= 'X')
h=merge(f,Hu, by= 'X')
Final= merge(h,Act, by= 'X')

write.table(Final,'Merged_data.csv',  append=FALSE, sep = ';')
df= read.csv('Merged_data.csv', row.names = 1, header = T, sep = ';')
