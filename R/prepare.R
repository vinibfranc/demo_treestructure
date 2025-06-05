library(devtools)
#install_github(repo="emvolz-phylodynamics/treestructure")
libs_load <- c("data.table", "ape", "dplyr", "glue", "lubridate", "treedater", "treestructure", "ggtree")
invisible( lapply(libs_load, library, character.only=TRUE) )

# 1) Download msa, nwk and metadata files here; https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2//2022/04/02/
# 2) Move files to data/; unzip nwk and metadata files manually; use xz --decompress public-2022-04-02.all.msa.fa.xz for fasta (xz-utils library)

#columns are: strain	genbank_accession	date	country	host	completeness	length	Nextstrain_clade	pangolin_lineage
md <- fread(file="data/metadata.tsv", sep="\t") #>4.637M seqs
md <- md %>% dplyr::select(strain, date, country, host, completeness, pangolin_lineage)

#tr <- read.tree(file="data/public-2021-05-31.all.nwk")

md_reduced <- md[md$date != "?" & (md$host  %in% c("","Homo sapiens")) & md$completeness != "Partial",]
md_reduced$date <- as.Date(md_reduced$date, format='%Y-%m-%d', na.rm = TRUE)
md_reduced <- md_reduced[complete.cases(md_reduced$date), ] #keep only complete dates: >4.623M seqs

md1 <- md_reduced[md_reduced$date <= "2020-02-29",] #929
md2 <- md_reduced[md_reduced$date <= "2020-03-15",] #mid-march: 6781; end of march: >26k seqs

write.table(md1$strain, file="data/md1.txt", col.names = F, row.names = F, quote = F)
write.table(md2$strain, file="data/md2.txt", col.names = F, row.names = F, quote = F)

# fasta file is huge to load in R, so writing txt with strains and extracting those with seqkit
system("seqkit grep -n -f data/md1.txt data/public-2022-04-02.all.msa.fa -o data/seqs1_up_2020_02_29.fa")
system("seqkit grep -n -f data/md2.txt data/public-2022-04-02.all.msa.fa -o data/seqs2_up_2020_03_15.fa")
# removed unzipped fasta after as >100 GB after these filters

# seqs1: removed manually from fasta the seqs below (>100 gaps in beginning or end); keeping USA/WA1/2020|MZ433205.1|2020-01-25 as first in the US,
# leaving the ones with Ns in the middle, in the end mask all cols from 29695 (209 final cols): 917 remaining
# ALSO removed this one (USA/MA-MASPHL-08626/2019|OM322554.1|2019-10-07 listed below) because from 2019-10-07?
seqs1_to_rm <- c("DP0058|LC570965.1|2020-02-15", "Iceland/13/2020|OU161584.1|2020-02-27", "Netherlands/Berlicum_1363564/2020|LR878228.1|2020-02-24",
																	"CHN/WHUHnCoV020/2020|MT079853.1|2020-01-22", "TWN/CGMH-CGU-04/2020|MT370517.1|2020-02-27", "USA/WA-S8/2020|MT598638.1|2020-02-24",
																	"USA/WA-S113/2020|MT627239.1|2020-02-29", "CHN/Meizhou/2020|MT856370.1|2020-01-26",
																	"DP0462|LC571000.1|2020-02-16", "DP0690|LC571014.1|2020-02-17", "DP0700|LC571020.1|2020-02-17", "DP0779|LC571032.1|2020-02-17", 
																	"USA/MA-MASPHL-08626/2019|OM322554.1|2019-10-07")
md1 <- md1[!(md1$strain %in% seqs1_to_rm),]
md1$strain[md1$strain == "Wuhan/WH19008/2019?|NMDC60013002-06|2019-12-30"] <- "Wuhan/WH19008/2019_|NMDC60013002-06|2019-12-30"
# Wuhan/WH19008/2019?|NMDC60013002-06|2019-12-30 changed ? for _ (aln and md1 file after estimating ML tree)

# iqtree ML tree run for first dataset (up to end onf february 2020)
HOME <- Sys.getenv("HOME")
IQTREE2_BIN <- glue("{HOME}/tools/iqtree-2.2.2.6-Linux/bin/iqtree2")
# Using -keep-ident to make sure identical seqs have bootstrap support assigned
system(glue("{IQTREE2_BIN} -s data/seqs1_up_2020_02_29.fa -m GTR+F+R3 -nt AUTO -ntmax 3 -B 1000 -nm 500 -keep-ident > /dev/null 2>&1 &"))

# look at ML tree in figtree: Japan sequences "rooting" tree and actual root as internal node?
# inspect in tempest
md1_tempest <- md1 %>% dplyr::select(strain, date)
colnames(md1_tempest) <- c("sequence_name","sample_date")
write.table(md1_tempest, file=glue("data/seqs1_up_2020_02_29_dates.tsv"), quote = F, sep="\t", row.names = F)
# some overlapping with outliers from treedater (below), but tmrca dates go to early 2019 or 2018 here

mltr <- ape::read.tree(glue("data/seqs1_up_2020_02_29.fa.treefile"))
# will use Wuhan/WH01/2019|LR757998.1|2019-12-26 as root
mltr_root <- ape::root(mltr, outgroup="Wuhan/WH01/2019|LR757998.1|2019-12-26", resolve.root=TRUE)
plot(mltr_root, show.tip.label = F)
sts <- decimal_date(as.Date(md1$date))
names(sts) <- md1$strain

NCPU <- 3
timetr <- dater(mltr_root, sts, s=29903, minblen=1/365, quiet=FALSE, clock="strict", omega0=0.0008, meanRateLimits=c(0.0002, 0.0010), parallel_foreach = TRUE, ncpu=NCPU) #clock="additive", maxit=10
plot(timetr, show.tip.label = T)

RES_RTT <- "results/rtt"
system(glue("mkdir -p {RES_RTT}"))

pdf(glue("{RES_RTT}/timetr_rtt.pdf"))
rootToTipRegressionPlot(timetr, show.tip.labels = F)
dev.off()

outliers <- outlierTips(timetr, alpha=0.05)
pdf(glue("{RES_RTT}/q_outliers.pdf"))
hist(outliers$q, breaks=200)
dev.off()
# Remove all tips that don't have a high q-value
to_rm <- outliers[(outliers$q < 0.05) ,]
print("TO RM:")
print(nrow(to_rm)) #25

mltr2_outl_rm <- drop.tip(mltr_root, rownames(to_rm))
timetr2 <- dater(mltr2_outl_rm, sts, s=29903, minblen=1/365, quiet=FALSE, clock="strict", omega0=0.0008, meanRateLimits=c(0.0002, 0.0010), parallel_foreach = TRUE, ncpu=NCPU)
plot(timetr2, show.tip.label = F)

pdf(glue("{RES_RTT}/timetr2_rtt.pdf"))
rootToTipRegressionPlot(timetr2, show.tip.labels = F)
dev.off()

timetr2_boot <- as.integer(mltr2_outl_rm$node.label)
timetr2_boot[is.na(timetr2_boot)] <- 95
print(length(timetr2_boot))

timetr2_phylo <- timetr2
class(timetr2_phylo) <- "phylo"

RDS_VIG <- "results/export_vignette"
system(glue("mkdir -p {RDS_VIG}"))
saveRDS(timetr2_phylo, glue("{RDS_VIG}/timetr2_phylo_sc2_feb2020.rds"))
saveRDS(mltr2_outl_rm, glue("{RDS_VIG}/mltr2_outl_rm_sc2_feb2020.rds"))

#mcs50 and BT90 (2 phylotypes)
#mcs25 and BT80 (4 phylotypes)
trestruct_res <- trestruct(timetr2_phylo, minCladeSize=30, nodeSupportValues=timetr2_boot, nodeSupportThreshold=80, level=0.01, ncpu=NCPU) #4 pts
saveRDS(trestruct_res, glue("{RDS_VIG}/trestruct_res.rds"))
trestruct_res_df <- as.data.frame(trestruct_res)

trestruct_res_nobt <- trestruct(timetr2_phylo, minCladeSize=30, nodeSupportValues=FALSE, level=0.01, ncpu=NCPU) #13 pts
saveRDS(trestruct_res_nobt, glue("{RDS_VIG}/trestruct_res_nobt.rds"))

RES_TRESTR <- "results/trestruct"
system(glue("mkdir -p {RES_TRESTR}"))


pdf(glue("{RES_TRESTR}/trestruct_sup.pdf"))
plot(trestruct_res, use_ggtree = T) + ggtree::geom_tippoint()
dev.off()

pdf(glue("{RES_TRESTR}/trestruct_nosup.pdf"))
plot(trestruct_res_nobt, use_ggtree = T) + ggtree::geom_tippoint()
dev.off()

# for adding tips to existing tree will work with data up to mid-March 2020
seqs2_fa <- read.dna("data/seqs2_up_2020_03_15.fa", format="fasta") #6780
seqs2_fa_filt1 <- seqs2_fa[!(rownames(seqs2_fa) %in% c(seqs1_to_rm,to_rm$taxon)),] #6742, remove previous seqs with too many gaps and outliers from treedater
write.dna(seqs2_fa_filt1, file="data/seqs2_adj1_up_2020_03_15.fa", format="fasta") #quite slow

# manually adjust ends of seqs2_adj1_up_2020_03_15.fa (masking now the first 300 nucleotides instead of 100 and from 29694 to end as before) without removing manually potentially problematic as would take a long time

md2$strain[md2$strain == "Wuhan/WH19008/2019?|NMDC60013002-06|2019-12-30"] <- "Wuhan/WH19008/2019_|NMDC60013002-06|2019-12-30"
md2 <- md2[md2$strain %in% rownames(seqs2_fa_filt1),] #6741

RDS_F <- "results/rds"
system(glue("mkdir -p {RDS_F}"))
save.image(file=glue("{RDS_F}/all_before_mltr_add.RData"))

# running
system(glue("{IQTREE2_BIN} -s data/seqs2_adj1_up_2020_03_15.fa -m GTR+F+R3 -nt AUTO -ntmax 3 -B 1000 -nm 500 -keep-ident > /dev/null 2>&1 &"))

load(file=glue("{RDS_F}/all_before_mltr_add.RData"))

# look at ML tree in figtree
mltr_add <- ape::read.tree(glue("data/seqs2_adj1_up_2020_03_15.fa.treefile"))
mltr_add_root <- ape::root(mltr_add, outgroup="Wuhan/WH01/2019|LR757998.1|2019-12-26", resolve.root=TRUE)
# remove a few seqs with long branches
rm_long_br <- c("KOR/CNUHV03/2020|MT678839.1|2020-03-11", "POL/IHG_PAS_25_P7726/2020|MZ047087.1|2020-03-15", "JPN/JP_HiroFH173c/2020|OM816720.1|2020-03-02",
																"USA/MA-UMASSMED-P007D03/2020|OM485489.1|2020-03-15", "USA/MA-UMASSMED-P007D04/2020|OM485490.1|2020-03-15", "USA/WV091040/2020|MZ500330.1|2020-01-21",
																"USA/MA-UMASSMED-P007C06/2020|OM485484.1|2020-03-15", "USA/MA-UMASSMED-P007D08/2020|OM485491.1|2020-03-15", "USA/MA-UMASSMED-P007E10/2020|OM485497.1|2020-03-15",
																"POL/IHG_PAS_25_P7908/2020|MZ047082.1|2020-03-15", "USA/MA-UMASSMED-P007D11/2020|OM485492.1|2020-03-15", "USA/MA-UMASSMED-P007C03/2020|OM485483.1|2020-03-15",
																"Switzerland/VD-CHUV-GEN2474/2020|OU575951.1|2020-03-09", "USA/MA-UMASSMED-P007C10/2020|OM485486.1|2020-03-15", "USA/MA-UMASSMED-P007C11/2020|OM485487.1|2020-03-15",
																"USA/MA-UMASSMED-P007B12/2020|OM485481.1|2020-03-15", "USA/MA-UMASSMED-P007E06/2020|OM485495.1|2020-03-15", "USA/MA-UMASSMED-P007C01/2020|OM485482.1|2020-03-15",
																"POL/IHG_PAS_25_7726/2020|MZ047085.1|2020-03-15","Japan/Donner381/2020|BS001050.1|2020-02-15","Japan/Donner380/2020|BS001049.1|2020-02-01","USA/MA-UMASSMED-P007C07/2020|OM485485.1|2020-03-15",
																"USA/MA-UMASSMED-P007D01/2020|OM485488.1|2020-03-15","USA/WV064569/2020|MZ500497.1|2020-01-25","USA/WV064576/2020|MZ500751.1|2020-01-25","Japan/Donner379/2020|BS001048.1|2020-02-15",
																"Japan/Donner382/2020|BS001051.1|2020-02-24", "JPN/JP_Hiro84244c/2020|OL636004.1|2020-03-02", "USA/MA-UMASSMED-P007E08/2020|OM485496.1|2020-03-15",
																"USA/FL-BPHL-2024/2020|MW286530.1|2020-03-13")
# rm_long_br <- c("Switzerland/ZH-UZH-1000477806/2020|FR990353.1|2020-02-29")

mltr_add_root <- drop.tip(phy=mltr_add_root, tip=rm_long_br) # 6712 tips
write.tree(mltr_add_root, file="data/seqs2_adj1_up_2020_03_15_rooted.fa.treefile")

saveRDS(mltr_add_root, glue("{RDS_VIG}/mltr_addtips_mar2020.rds"))

trestruct_add_tips <- addtips(trst=trestruct_res, tre=mltr_add_root)
pdf(glue("{RES_TRESTR}/trestruct_added_tips.pdf"))
plot(trestruct_add_tips, use_ggtree = T) + ggtree::geom_tippoint()
dev.off()
