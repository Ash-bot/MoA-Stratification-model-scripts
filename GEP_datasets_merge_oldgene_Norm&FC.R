
#old ID trail

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Identification of a New Chemical Class of Antimalarials ACT213615/Series Matrix File/')

Art1= read.csv(file = 'ACTderiv processed probes.csv',header = T, sep = ';')

Art1[Art1 == "NaN"] = NA
Art1[Art1 == "null"] = NA
sum(is.na(Art1))
# [1] 2784
is.data.frame(Art1)

Art2 = data.frame(Gene.ID = Art1$ID_REF,apply(Art1[,2:11],2,as.numeric))
is.numeric(Art2$ACT.213615_TP1)
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


sum(is.na(Art3))
Art = na.omit(Art3)
# HDA

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Three different histone deacetylase inhibitors/Series Matrix File/')

HDA1= read.delim(file = 'HDAi probes processed.txt',header = T)

HDA1[HDA1 == "NaN"] = NA
HDA1[HDA1 == "null"] = NA
sum(is.na(HDA1))
# [1] 65260
is.data.frame(HDA1)

HDA2 = data.frame(Gene.ID = HDA1$ID_REF,apply(HDA1[,2:17],2,as.numeric))
is.numeric(HDA2$DMSO.control.2h.replicate.1)
# HDA2= na.omit(HDA2)
# HDA2= HDA2[-which(rowMeans(is.na(HDA2)) > 0.5), ]
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

sum(is.na(HDA3))
# [1] 2060
HDAi= na.omit(HDA3)


setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Role of Calcium Signaling in the Transcriptional Regulation of the Apicoplast Genome of Plasmodium falciparum/')

Iono1= read.csv(file = 'Ionomycin probes processed.csv',header = T, sep = ';')

Iono1[Iono1 == "NaN"] = NA
Iono1[Iono1 == "null"] = NA

sum(is.na(Iono1))
# [1] 8231
is.data.frame(Iono1)

Iono2 = data.frame(Gene.ID = Iono1$ID_REF,apply(Iono1[,2:16],2,as.numeric))
is.numeric(Iono2$P.falciparum.3D7_schizonts_extracellular.calcium.chelation._Ionomycin.Treatment_TP1_30mins)
# Iono2= na.omit(Iono2)
# Iono2= Iono2[-which(rowMeans(is.na(Iono2)) > 0.5), ]
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



sum(is.na(Iono3))
# [1] 1803
# [1] 1955
Iono = na.omit(Iono3)


###############################
setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Plasmodium falciparum asexual parasites and late gametocytes MMV390048 or MMV642943/Series Matrix File/')

MMV1= read.csv(file = 'MMV compounds need probes.csv',header = T, sep = ';')

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Plasmodium falciparum asexual parasites and late gametocytes MMV390048 or MMV642943/')

Genes= read.csv( file = 'GeneID.csv',header = T, sep = ';')
Genes= Genes[-2]
G= na.omit(Genes)
l= length(unique(Genes$ORF))
print(l)

MMV1$Gene_ID = Genes$ORF[match(MMV1$ID_REF, Genes$ID)]
MMV2= MMV1[,c(12,2:11)]

MMV2[MMV2 == "NaN"] = NA
MMV2[MMV2 == "null"] = NA
sum(is.na(MMV2))
# [1] 34412

is.data.frame(MMV2)

MMV3 = data.frame(Gene.ID = MMV2$Gene_ID,apply(MMV2[,2:11],2,as.numeric))
is.numeric(MMV2$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV642943.and.harvested.24.h.post.expt.initiation..A24.T943.)
# MMV3= na.omit(MMV3)
# MMV3= MMV3[-which(rowMeans(is.na(MMV3)) > 0.5), ]
l= length(unique(MMV3$Gene.ID))
print(l)
# [1] 5643
MMV3 = aggregate(.~Gene.ID, data = MMV3, FUN = mean, na.rm = T, na.action = na.pass)

MMV3[MMV3 == "NaN"] = NA
MMV3[MMV3 == "null"] = NA
sum(is.na(MMV3))
# [1] 10571

MMV4 = data.frame(row.names = MMV3$Gene.ID,
                  MMV_43TP1= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV642943.and.harvested.24.h.post.expt.initiation..A24.T943.,
                  MMV_43TP2= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV642943.and.harvested.48.h.post.expt.initiation..A48.T943.,
                  MMV_48TP1= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV390048.and.harvested.48.h.post.expt.initiation..A24.T048.,
                  MMV_48TP2= MMV3$NF54.asexual.ring.stage.parasites.treated.with.10.x.IC50.of.MMV390048.and.harvested.24.h.post.expt.initiation..A48.T048.,
                  MMV_TP1c= MMV3$NF54.asexual.ring.stage.parasites.remained.untreated.and.harvested.24.h.post.expt.initiation..A24.UT.,
                  MMV_TP2c= MMV3$NF54.asexual.ring.stage.parasites.remained.untreated.and.harvested.48.h.post.expt.initiation..A48..UT.)


sum(is.na(MMV4))
# [1] 2515

MMV= na.omit(MMV4)
sum(is.na(MMV))

##########DFMO

setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Accepted/Functional genomics investigation of PfAdoMetDCODC coinhibition of Plasmodium falciparum/Series Matrix File/')

DFMO1= read.csv(file = 'DMFOseriesmatrix.csv',header = T, sep = ';')

Genes= read.delim( file = 'GeneID.csv',header = T, sep = ';')
Genes= Genes[-3]

Genes[Genes == "N/A"] = NA
Genes[Genes == "NULL"] = NA
l= length(unique(Genes$PLASMODB.ID))
print(l)

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
# DFMO4= na.omit(DFMO4)
# DFMO4= DFMO4[-which(rowMeans(is.na(DFMO4)) > 0.5), ]
l= length(unique(DFMO4$Gene.ID))
print(l)
# [1] 5012
DFMO4 = aggregate(.~Gene.ID, data = DFMO4, FUN = mean, na.rm = T, na.action = na.pass)

DFMO4[DFMO4 == "NaN"] = NA
DFMO4[DFMO4 == "null"] = NA
sum(is.na(DFMO4))
# [1] 17237


DFMO5 = data.frame(row.names = DFMO4$Gene.ID,
                   DFMOT19= rowMeans(DFMO4[c('X127_Pf_T_t19_repB', 'X116_Pf_T_t19_repB_tech')], na.rm=TRUE),
                   DFMOT19c= rowMeans(DFMO4[c('X095_Pf_UT_t19_repA_tech', 'X107_Pf_UT_t19_repA','X124_Pf_UT_t19_repB')], na.rm=TRUE),
                   DFMOT27= rowMeans(DFMO4[c('X128_Pf_T_t27_repB', 'X098_Pf_T_t27_repA_tech','X111_Pf_T_t27_repA')], na.rm=TRUE),
                   DFMOT27c= rowMeans(DFMO4[c('X125_Pf_UT_t27_repB', 'X091_Pf_UT_t27_repA_tech','X108_Pf_UT_t27_repA')], na.rm=TRUE),
                   DFMOT34= rowMeans(DFMO4[c('X094_Pf_T_t34_repB_tech','X129_Pf_T_t34_repB', 'X090_Pf_T_t34_repA_tech','X112_Pf_T_t34_repA')], na.rm=TRUE),
                   DFMOT34c= rowMeans(DFMO4[c('X109_Pf_UT_t34_repA','X097_Pf_UT_t34_repA_tech','X126_Pf_UT_t34_repB', 'X129_Pf_T_t34_repB','X089_Pf_UT_t34_repB_tech')], na.rm=TRUE))

DFMO5[DFMO5 == "NaN"] = NA
DFMO5[DFMO5 == "null"] = NA
sum(is.na(DFMO4))

sum(is.na(DFMO5))
# [1] 17793
DFMO= na.omit(DFMO5)
#################HU

setwd('C:/Users/ashle/Desktop/LAbbook(2019)')

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
# 
Hu= na.omit(selected)


setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Old gene ID/1356 genes/')
write.table(DFMO, file= 'DFMO.csv', col.names=TRUE, row.names=TRUE, append=FALSE, sep = ';')
write.table(HDAi, file= 'HDA.csv',  col.names=TRUE, row.names=TRUE, append=FALSE, sep = ';')
write.table(Hu, file= 'Hu.csv', col.names=TRUE, row.names=TRUE, append=FALSE, sep = ';')
write.table(Art, file= 'Art.csv',  col.names=TRUE, row.names=TRUE, append=FALSE, sep = ';')
write.table(Iono, file= 'Iono.csv',col.names=TRUE, row.names=TRUE, append=FALSE, sep = ';')
write.table(MMV, file= 'MMV.csv',col.names=TRUE, row.names=TRUE, append=FALSE, sep = ';')
#Just make sure when recalling these files you go in and make the First row X
setwd('/Users/ashle/Desktop/LAbbook(2019)/Additional datasets/Old gene ID/1356 genes/')
Iono= read.csv('Iono.csv', sep = ',')# , row.names = 1)
Iono1=Iono[!rowSums(Iono[,2:11] >7),]
Iono1= na.omit(Iono1)
# max(Iono) NB!! iono is weird has high values in 88 
#so we just put a limit of no higher than 7 as the other genes look fine, so there may have 
#been some error with these genes during microarray
boxplot(Iono1[2:11], las=2)

Act= read.csv('Art.csv', header = T, sep = ',')
Act= na.omit(Act)

DFMO= read.csv('DFMO.csv', header = T, sep = ',')
DFMO= na.omit(DFMO)

HDA= read.csv('HDA.csv', header = T, sep = ',')
HDA= na.omit(HDA)

MMV= read.csv('MMV.csv', header = T, sep = ',')
MMV= na.omit(MMV)
Hu= read.csv('Hu.csv', header = T, sep = ',')
Hu= na.omit(Hu)

a =merge(Hu,DFMO,by = 'X')
b= merge(a,MMV, by= 'X')
c= merge(b,Iono1, by= 'X')
d= merge(c,HDA, by= 'X')
e= merge(d,Act, by= 'X')
# f= merge(e,Iono1, by= 'X')
# g=merge(f,Cyc, by= 'X')
# h=merge(g,DHB, by= 'X')
# # Final= merge(h,Act, by= 'X')

write.table(e, file= 'Merged data.csv', row.names = FALSE, append=FALSE, sep = ',')
df=read.csv('Merged data.csv', row.names = 1, header = T, sep = ',')
as.numeric(df)


library(limma)
boxplot(df,
        las = 2,
        main= "Unnormalized")
loess = normalizeCyclicLoess(df)
boxplot(loess, las =2,
        main= "Cyclic Loess")

quantile = normalizeQuantiles(df)
boxplot(quantile, las =2,
        main= "Quantile")

scale = normalizeMedianAbsValues(df)
boxplot(scale, las =2,
        main= "Medium scaling (absolute values)")

scale = normalizeMedianValues(df1)
boxplot(scale, las =2,
        main= "Medium scaling (absolute values)")

write.table(loess, 'Norm_MergedData.csv', sep = ',')
################log FC
loess1= read.csv('Norm_MergedData.csv', sep = ',', row.names = 1)
loess= as.data.frame(2^loess1)
FC= data.frame(row.names = rownames(loess),
               HDATP1= ((loess$HDATP1)/(loess$HDATP1c)),
               HDATP1.1= ((loess$HDATP1.1)/(loess$HDATP1c2)),
               HDATP2= ((loess$HDATP2)/(loess$HDATP1c)),
               HDATP2.2= ((loess$HDATP2.2)/(loess$HDATP1c2)),
               HDATP3= ((loess$HDATP3)/(loess$HDATP1c)),
               HDATP3.3= ((loess$HDATP3.3)/(loess$HDATP1c2)),
               ACTTP1= ((loess$ArtTP1)/(loess$ArtTP1c)),
               ACTTP2= ((loess$ArtTP2)/(loess$ArtTP2c)),
               ACTTP3= ((loess$ArtTP3)/(loess$ArtTP3c)),
               ACTTP4= ((loess$ArtTP4)/(loess$ArtTP4c)),
               ACTTP5= ((loess$ArtTP5)/(loess$ArtTP5c)),
               ML_7_H1=((loess$ML_7_H1)/(loess$ML_7_H1c)),
               ML_7_H2=((loess$ML_7_H2)/(loess$ML_7_H2c)),
               ML_7_H4=((loess$ML_7_H4)/(loess$ML_7_H4c)),
               ML_7_H6=((loess$ML_7_H6)/(loess$ML_7_H6c)),
               ML_7_H8=((loess$ML_7_H8)/(loess$ML_7_H8c)),
               ML_7_H10=((loess$ML_7_H10)/(loess$ML_7_H10c)),
               ML_7_H12=((loess$ML_7_H12)/(loess$ML_7_H12c)),
               W7_H1=((loess$W7_H1)/(loess$W7_H1c)),
               W7_H2=((loess$W7_H2)/(loess$W7_H2c)),
               W7_H4=((loess$W7_H4)/(loess$W7_H4c)),
               W7_H6=((loess$W7_H6)/(loess$W7_H6c)),
               W7_H8=((loess$W7_H8)/(loess$W7_H8c)),
               W7_H10=((loess$W7_H10)/(loess$W7_H10c)),
               W7_H12=((loess$W7_H12)/(loess$W7_H12c)),
               Staurosporine_H1=((loess$Staurosporine_H1)/(loess$Staurosporine_H1c)),
               Staurosporine_H2=((loess$Staurosporine_H2)/(loess$Staurosporine_H2c)),
               Staurosporine_H4=((loess$Staurosporine_H4)/(loess$Staurosporine_H4c)),
               Staurosporine_H6=((loess$Staurosporine_H6)/(loess$Staurosporine_H6c)),
               Staurosporine_H8=((loess$Staurosporine_H8)/(loess$Staurosporine_H8c)),
               Staurosporine_H10=((loess$Staurosporine_H10)/(loess$Staurosporine_H10c)),
               Staurosporine_H12=((loess$Staurosporine_H12)/(loess$Staurosporine_H12c)),
               Cyclosporine_H1=((loess$Cyclosporine_H1)/(loess$Cyclosporine_H1c)),
               Cyclosporine_H2=((loess$Cyclosporine_H2)/(loess$Cyclosporine_H2c)),
               Cyclosporine_H4=((loess$Cyclosporine_H4)/(loess$Cyclosporine_H4c)),
               Cyclosporine_H6=((loess$Cyclosporine_H6)/(loess$Cyclosporine_H6c)),
               Cyclosporine_H8=((loess$Cyclosporine_H8)/(loess$Cyclosporine_H8c)),
               Cyclosporine_H10=((loess$Cyclosporine_H10)/(loess$Cyclosporine_H10c)),
               Cyclosporine_H12=((loess$Cyclosporine_H12)/(loess$Cyclosporine_H12c)),
               Cyclosporine_H14=((loess$Cyclosporine_H14)/(loess$Cyclosporine_H14c)),
               Colchine_H1=((loess$Colchine_H1)/(loess$Colchine_H1c)),
               Colchine_H2=((loess$Colchine_H2)/(loess$Colchine_H2c)),
               Colchine_H4=((loess$Colchine_H4)/(loess$Colchine_H4c)),
               Colchine_H6=((loess$Colchine_H6)/(loess$Colchine_H6c)),
               Colchine_H8=((loess$Colchine_H8)/(loess$Colchine_H8c)),
               PMSF_H1=((loess$PMSF_H1)/(loess$PMSF_H1c)),
               PMSF_H2=((loess$PMSF_H2)/(loess$PMSF_H2c)),
               PMSF_H4=((loess$PMSF_H4)/(loess$PMSF_H4c)),
               PMSF_H6=((loess$PMSF_H6)/(loess$PMSF_H6c)),
               PMSF_H8=((loess$PMSF_H8)/(loess$PMSF_H8c)),
               Leupeptine_H1=((loess$Leupeptine_H1)/(loess$Leupeptine_H1c)),
               Leupeptine_H2=((loess$Leupeptine_H2)/(loess$Leupeptine_H2c)),
               Leupeptine_H4=((loess$Leupeptine_H4)/(loess$Leupeptine_H4c)),
               Leupeptine_H6=((loess$Leupeptine_H6)/(loess$Leupeptine_H6c)),
               Leupeptine_H8=((loess$Leupeptine_H8)/(loess$Leupeptine_H8c)),
               Apicidin_H1=((loess$Apicidin_H1)/(loess$Apicidin_H1c)),
               Apicidin_H2=((loess$Apicidin_H2)/(loess$Apicidin_H2c)),
               Apicidin_H4=((loess$Apicidin_H4)/(loess$Apicidin_H4c)),
               Apicidin_H6=((loess$Apicidin_H6)/(loess$Apicidin_H6c)),
               Apicidin_H8=((loess$Apicidin_H8)/(loess$Apicidin_H8c)),
               Apicidin_H10=((loess$Apicidin_H10)/(loess$Apicidin_H10c)),
               Trichostatin_H1=((loess$Trichostatin_H1)/(loess$Trichostatin_H1c)),
               Trichostatin_H2=((loess$Trichostatin_H2)/(loess$Trichostatin_H2c)),
               Trichostatin_H4=((loess$Trichostatin_H4)/(loess$Trichostatin_H4c)),
               Trichostatin_H6=((loess$Trichostatin_H6)/(loess$Trichostatin_H6c)),
               Trichostatin_H8=((loess$Trichostatin_H8)/(loess$Trichostatin_H8c)),
               Trichostatin_H10=((loess$Trichostatin_H10)/(loess$Trichostatin_H10c)),
               Quinine_H1=((loess$Quinine_H1)/(loess$Quinine_H1c)),
               Quinine_H2=((loess$Quinine_H2)/(loess$Quinine_H2c)),
               Quinine_H4=((loess$Quinine_H4)/(loess$Quinine_H4c)),
               Quinine_H6=((loess$Quinine_H6)/(loess$Quinine_H6c)),
               Quinine_H8=((loess$Quinine_H8)/(loess$Quinine_H8c)),
               Quinine_H10=((loess$Quinine_H10)/(loess$Quinine_H10c)),
               Chloroquine_H1=((loess$Chloroquine_H1)/(loess$Chloroquine_H1c)),
               Chloroquine_H2=((loess$Chloroquine_H2)/(loess$Chloroquine_H2c)),
               Chloroquine_H4=((loess$Chloroquine_H4)/(loess$Chloroquine_H4c)),
               Chloroquine_H6=((loess$Chloroquine_H6)/(loess$Chloroquine_H6c)),
               Chloroquine_H8=((loess$Chloroquine_H8)/(loess$Chloroquine_H8c)),
               Chloroquine_H10=((loess$Chloroquine_H10)/(loess$Chloroquine_H10c)),
               Febrifugine_H1=((loess$Febrifugine_H1)/(loess$Febrifugine_H1c)),
               Febrifugine_H2=((loess$Febrifugine_H2)/(loess$Febrifugine_H2c)),
               Febrifugine_H4=((loess$Febrifugine_H4)/(loess$Febrifugine_H4c)),
               Febrifugine_H6=((loess$Febrifugine_H6)/(loess$Febrifugine_H6c)),
               Febrifugine_H8=((loess$Febrifugine_H8)/(loess$Febrifugine_H8c)),
               Febrifugine_H10=((loess$Febrifugine_H10)/(loess$Febrifugine_H10c)),
               Artemisinin_H1=((loess$Artemisinin_H1)/(loess$Artemisinin_H1c)),
               Artemisinin_H2=((loess$Artemisinin_H2)/(loess$Artemisinin_H2c)),
               Artemisinin_H4=((loess$Artemisinin_H4)/(loess$Artemisinin_H4c)),
               Artemisinin_H6=((loess$Artemisinin_H6)/(loess$Artemisinin_H6c)),
               Artemisinin_H8=((loess$Artemisinin_H8)/(loess$Artemisinin_H8c)),
               Artemisinin_H10=((loess$Artemisinin_H10)/(loess$Artemisinin_H10c)),
               MMV43TP1= ((loess$MMV_43TP1)/(loess$MMV_TP1c)),
               MMV43TP2= ((loess$MMV_43TP2)/(loess$MMV_TP2c)),
               MMV48TP1= ((loess$MMV_48TP1)/(loess$MMV_TP1c)),
               MMV48TP2= ((loess$MMV_48TP2)/(loess$MMV_TP2c)),
               DFMO_TP19= ((loess$DFMOT19)/(loess$DFMOT19c)),
               DFMO_TP27= ((loess$DFMOT27)/(loess$DFMOT27c)),
               DFMO_TP34= ((loess$DFMOT34)/(loess$DFMOT34c)),
               # CycTP25= ((loess$CycT25)/(loess$CycT25c)),
               # CycTP30= ((loess$CycT30)/(loess$CycT30c)),
               # DHBTP2= ((loess$DHBT2)/(loess$DHBT2c)),
               # DHBTP6= ((loess$DHBT6)/(loess$DHBT30)),
               # DHBTP12= ((loess$DHBT12)/(loess$DHBT30c)),
               IonoTP1= ((loess$IonoTP1)/(loess$IonoTP1c)),
               IonoTP2= ((loess$IonoTP2)/(loess$IonoTP2c)),
               IonoTP3= ((loess$IonoTP3)/(loess$IonoTP3c)),
               IonoTP4= ((loess$IonoTP4)/(loess$IonoTP4c)),
               IonoTP5= ((loess$IonoTP5)/(loess$IonoTP5c)))
# CisTP1= ((loess$)/(loess$)),
# CisTP2= ((loess$)/(loess$)),
# CisTP3= ((loess$)/(loess$)),
# MMSTP4= ((loess$)/(loess$)),
# MMSTP5= ((loess$)/(loess$)))

FC1= log2(FC)

write.table(FC1, "FC_LoessNorm_mergdata.csv", sep=';')

