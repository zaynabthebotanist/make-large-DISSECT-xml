
library(seqinr)
library(ShortRead)
library(stringr)
library(DescTools)
library(XML)
library(schoolmath)

## import individual alignments

files <- list.files(pattern = "*.fasta")
n <- c(1:length(files))
file <- list()
for (i in n) {
  file[[i]] <- readDNAStringSet(files[[i]], "fasta")
}

# 1.) sequence data 

data_id <- strsplit(files, ".fasta")

seq_id <- vector("list", length = 524)
n <- c(1:length(seq_id))
for (i in n) {
  seq_id[[i]] <- paste0("seq_", names(file[[i]]))
}

append <- as.list(c("a", c(1:523)))
seq_id1 <- list()
for(i in n) {
  seq_id1[[i]] <- paste0(seq_id[[i]], append[[i]])
}
seq_id1[[1]] <- str_replace(seq_id1[[1]], "_hap1a", "_hap1")
seq_id1[[1]] <- str_replace(seq_id1[[1]], "_hap2a", "_hap2")

n <- c(1:length(files))
taxon <- list()
for (i in n) {
  taxon[[i]] <- names(file[[i]])
}  
  
value <- list() 
for (i in n) {
  value[[i]] <- as.character(file[[i]])
}  

dataline_beg <- list()
for (i in n) {
  dataline_beg[[i]] <- paste0('<data id=\"', data_id[[i]], '\" name="alignment">')  ##remove the \ later on
}

seqs <- list()
for (i in n) {
  seqs[[i]] <- paste0('<sequence id=\"', seq_id1[[i]], '\" taxon=\"', taxon[[i]], '\" totalcount=\"4\" value=\"', value[[i]], '\"/>')
}

dataline_end <- as.list(rep("</data>", 524))

combined <- list()
for (i in n) {
  combined[[i]] <- c(dataline_beg[[i]], seqs[[i]], dataline_end[[i]])
}

fileConn <- file("seqdata.txt")
writeLines(unlist(combined), fileConn)
close(fileConn)


# 2.) clock rate models

#<prior id="bdcGrowthRatePrior.t:Species" name="distribution" x="@bdcGrowthRate.t:Species">
#  <LogNormal id="LogNormalDistributionModel.0" name="distr">
#  <parameter id="RealParameter.10" estimate="false" name="M">1.0</parameter>
#  <parameter id="RealParameter.11" estimate="false" lower="0.0" name="S" upper="5.0">1.25</parameter>
#  </LogNormal>
#  </prior>
#  <prior id="ClockPrior.c:C6422483" name="distribution" x="@clockRate.c:C6422483">
#  <LogNormal id="LogNormalDistributionModel.2" name="distr">
#  <parameter id="RealParameter.14" estimate="false" name="M">1.0</parameter>
#  <parameter id="RealParameter.15" estimate="false" lower="0.0" name="S" upper="5.0">1.25</parameter>
#  </LogNormal>
#  </prior>
#  <prior id="ClockPrior.c:C6423145" name="distribution" x="@clockRate.c:C6423145">
#  <LogNormal id="LogNormalDistributionModel.3" name="distr">
#  <parameter id="RealParameter.16" estimate="false" name="M">1.0</parameter>
#  <parameter id="RealParameter.17" estimate="false" lower="0.0" name="S" upper="5.0">1.25</parameter>
#  </LogNormal>
#  </prior>

#A few things need to be set up:
#(1) insert the chromosome ID
#(2) change the LogNormaDistributionModel number, start at 14
#(3) change the RealParameter number

data_id

a <- as.list(as.character(1:524))  # (for 2)

14 + (523*2)    ##how many numbers are needed for (3)
even <- seq(12, 1060, 2)
odd <- seq(13, 1061, 2)

##construct line-by-line

line1 <- list()
n <- c(1:length(data_id))
for (i in n) {
  line1[[i]] <- paste0('<prior id="ClockPrior.c:', data_id[[i]], '" name="distribution" x="@clockRate.c:', data_id[[i]], '">')
}

line2 <- list()
for (i in n) {
  line2[[i]] <- paste0('<LogNormal id="LogNormalDistributionModel.', a[[i]],'" name="distr">')
}

line3 <- list()
for (i in n) {
  line3[[i]] <- paste0('<parameter id="RealParameter.', even[[i]], '" estimate="false" name="M">1.0</parameter>')
}

line4 <- list()
for (i in n) {
  line4[[i]] <- paste0('<parameter id="RealParameter.', odd[[i]], '" estimate="false" lower="0.0" name="S" upper="5.0">1.25</parameter>')
}

line5 <- as.list(rep("</LogNormal>", 524))
line6 <- as.list(rep("</prior>", 524))

clockrates <- list()
for (i in n) {
  clockrates[[i]] <- c(line1[[i]], line2[[i]], line3[[i]], line4[[i]], line5[[i]], line6[[i]])
}

fileConn <- file("clockratedata.txt")
writeLines(unlist(clockrates), fileConn)
close(fileConn)


# 3.) tree rates need to be set relative to the first alignment

#<tree id="Tree.t:C6432679" name="stateNode">
#  <taxonset id="TaxonSet.C6432679" spec="TaxonSet">
#  <alignment idref="C6432679"/>
#  </taxonset>
#  </tree>
#  <parameter id="clockRate.c:C6432679" name="stateNode">1.0</parameter>

data_id

n <- c(1:length(data_id))
line1 <- list()
for (i in n) {
  line1[[i]] <- paste0('<tree id="Tree.t:', data_id[[i]], '" name="stateNode">')
}

line2 <- list()
for (i in n) {
  line2[[i]] <- paste0('<taxonset id="TaxonSet.', data_id[[i]], '" spec="TaxonSet">')
}

line3 <- list()
for (i in n) {
  line3[[i]] <- paste0('<alignment idref="', data_id[[i]], '"/>')
}

line4 <- as.list(rep("</taxonset>", 524))
line5 <- as.list(rep("</tree>", 524))

line6 <- list()
for (i in n) {
  line6[[i]] <- paste0('<parameter id="clockRate.c:', data_id[[i]], '" name="stateNode">1.0</parameter>')
}
  
  
n <- c(1:length(data_id))
clockratesparams <- list()
for (i in n) {
  clockratesparams[[i]] <- c(line1[[i]], line2[[i]], line3[[i]], line4[[i]],
                             line5[[i]], line6[[i]])
}

fileConn <- file("clockrateparamdata.txt")
writeLines(unlist(clockratesparams), fileConn)
close(fileConn)


# 4.) RandomGeneTree additions

# <init id="RandomGeneTree.t:C6432679" spec="beast.evolution.speciation.RandomGeneTree" initial="@Tree.t:C6432679" speciesTree="@Tree.t:Species" taxa="@C6432679">
# <populationModel id="RGTPopulationModel.t:C6432679" spec="ConstantPopulation">
# <parameter id="RGTPopSize.t:C6432679" name="popSize">0.05</parameter>
# </populationModel>
# </init>

n <- c(1:length(data_id))

line1 <- list()
for (i in n) {
  line1[[i]] <- paste0('<init id="RandomGeneTree.t:', data_id[[i]], 
                       '" spec="beast.evolution.speciation.RandomGeneTree" ', 'initial="@Tree.t:',
                       data_id[[i]], '" speciesTree="@Tree.t:Species" taxa="@', data_id[[i]], '">')
}

line2 <- list()
for (i in n) {
  line2[[i]] <- paste0('<populationModel id="RGTPopulationModel.t:', data_id[[i]], '" spec="ConstantPopulation">')
}

line3 <- list()
for (i in n) {
  line3[[i]] <- paste0('<parameter id="RGTPopSize.t:', data_id[[i]], '" name="popSize">0.05</parameter>')
}

line4 <- as.list(rep("</populationModel>", 524))
line5 <- as.list(rep("</init>", 524))

randomgenetree <- list()
for (i in n) {
  randomgenetree[[i]] <- c(line1[[i]], line2[[i]], line3[[i]], line4[[i]], line5[[i]])
}

fileConn <- file("randomgenetree.txt")
writeLines(unlist(randomgenetree), fileConn)
close(fileConn)


# 5.) gtreeandcoalfactor

#<geneTree id="gTreeCF.t:C6432679" spec="stacey.GtreeAndCoalFactor" tree="@Tree.t:C6432679"/>

gtreeandcoalfactor <- list()
for (i in n) {
  gtreeandcoalfactor[[i]] <- paste0('<geneTree id="gTreeCF.t:', data_id[[i]], '" spec="stacey.GtreeAndCoalFactor" tree="@Tree.t:', data_id[[i]], '"/>')
}

fileConn <- file("gtreeandcoalfactor.txt")
writeLines(unlist(gtreeandcoalfactor), fileConn)
close(fileConn)


# 6.) treelikelihood:

#  <distribution id="treeLikelihood.C6422483" spec="TreeLikelihood" data="@C6422483" tree="@Tree.t:C6422483">
#  <siteModel id="SiteModel.s:C6422483" spec="SiteModel">
#  <parameter id="mutationRate.s:C6422483" estimate="false" name="mutationRate">1.0</parameter>
#  <parameter id="gammaShape.s:C6422483" estimate="false" name="shape">1.0</parameter>
#  <parameter id="proportionInvariant.s:C6422483" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>
#  <substModel id="JC69.s:C6422483" spec="JukesCantor"/>
#  </siteModel>
#  <branchRateModel id="StrictClock.c:C6422483" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:C6422483"/>
#  </distribution>

line1 <- list()
for (i in n) {
  line1[[i]] <- paste0('<distribution id="treeLikelihood.', data_id[[i]], '" spec="TreeLikelihood" data="@', data_id[[i]], '" tree="@Tree.t:', data_id[[i]], '">')
}

line2 <- list()
for (i in n) {
  line2[[i]] <- paste0('<siteModel id="SiteModel.s:', data_id[[i]], '" spec="SiteModel">')
}

line3 <- list()
for (i in n) {
  line3[[i]] <- paste0('<parameter id="mutationRate.s:', data_id[[i]], '" estimate="false" name="mutationRate">1.0</parameter>')
}

line4 <- list()
for (i in n) {
  line4[[i]] <- paste0('<parameter id="gammaShape.s:', data_id[[i]], '" estimate="false" name="shape">1.0</parameter>')
}

line5 <- list()
for (i in n) {
  line5[[i]] <- paste0('<parameter id="proportionInvariant.s:', data_id[[i]], '" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>')
}

line6 <- list()
for (i in n) {
  line6[[i]] <- paste0('<substModel id="JC69.s:', data_id[[i]], '" spec="JukesCantor"/>')
}

line7 <- as.list(rep("</siteModel>", 524))

line8 <- list()
for (i in n) {
  line8[[i]] <- paste0('<branchRateModel id="StrictClock.c:', data_id[[i]], '" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:', data_id[[i]], '"/>')
}

line9 <- as.list(rep("</distribution>", 524))


treelikelihood <- list()
for (i in n) {
  treelikelihood[[i]] <- c(line1[[i]], line2[[i]], line3[[i]], line4[[i]], line5[[i]],
                           line6[[i]], line7[[i]], line8[[i]], line9[[i]])
}

fileConn <- file("treelikelihood.txt")
writeLines(unlist(treelikelihood), fileConn)
close(fileConn)

# 7.) genetreeshorts

#<geneTree idref="Tree.t:C6432679"/>

genetreeshorts <- list()
for (i in n) {
  genetreeshorts[[i]] <- paste0('<geneTree idref="Tree.t:', data_id[[i]], '"/>')
}

fileConn <- file("genetreeshorts.txt")
writeLines(unlist(genetreeshorts), fileConn)
close(fileConn)

# 8.) up

#<up idref="clockRate.c:C6432679"/>

up <- list()
for (i in n) {
  up[[i]] <- paste0('<up idref="clockRate.c:', data_id[[i]], '"/>')
}

fileConn <- file("up.txt")
writeLines(unlist(up), fileConn)
close(fileConn)


# 9.) down

#<down idref="Tree.t:C6432679"/>

down <- list()
for (i in n) {
  down[[i]] <- paste0('<down idref="Tree.t:', data_id[[i]], '"/>')
}

fileConn <- file("down.txt")
writeLines(unlist(down), fileConn)
close(fileConn)



# 10.) clockscaler

#<operator id="StrictClockRateScaler.c:C6432679" spec="ScaleOperator" parameter="@clockRate.c:C6432679" scaleFactor="0.75" weight="3.0"/>
#<operator id="treeScaler.t:C6432679" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:C6432679" weight="3.0"/>
#<operator id="treeRootScaler.t:C6432679" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:C6432679" weight="3.0"/>
#<operator id="UniformOperator.t:C6432679" spec="Uniform" tree="@Tree.t:C6432679" weight="30.0"/>
#<operator id="SubtreeSlide.t:C6432679" spec="SubtreeSlide" tree="@Tree.t:C6432679" weight="15.0"/>
#<operator id="narrow.t:C6432679" spec="Exchange" tree="@Tree.t:C6432679" weight="15.0"/>
#<operator id="wide.t:C6432679" spec="Exchange" isNarrow="false" tree="@Tree.t:C6432679" weight="3.0"/>
#<operator id="WilsonBalding.t:C6432679" spec="WilsonBalding" tree="@Tree.t:C6432679" weight="3.0"/>
#<operator id="updown.C6432679" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">
#<up idref="clockRate.c:C6432679"/>
#<down idref="Tree.t:C6432679"/>
#</operator>

line1 <- list()
for (i in n) {
  line1[[i]] <- paste0('<operator id="StrictClockRateScaler.c:', data_id[[i]], '" spec="ScaleOperator" parameter="@clockRate.c:', data_id[[i]], '" scaleFactor="0.75" weight="3.0"/>')
}

line2 <- list()
for (i in n) {
  line2[[i]] <- paste0('<operator id="treeScaler.t:', data_id[[i]], '" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:', data_id[[i]], '" weight="3.0"/>')
}

line3 <- list()
for (i in n) {
  line3[[i]] <- paste0('<operator id="treeRootScaler.t:', data_id[[i]], '" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:', data_id[[i]], '" weight="3.0"/>')
}

line4 <- list()
for (i in n) {
  line4[[i]] <- paste0('<operator id="UniformOperator.t:', data_id[[i]], '" spec="Uniform" tree="@Tree.t:', data_id[[i]], '" weight="30.0"/>')
}

line5 <- list()
for (i in n) {
  line5[[i]] <- paste0('<operator id="SubtreeSlide.t:', data_id[[i]], '" spec="SubtreeSlide" tree="@Tree.t:', data_id[[i]], '" weight="15.0"/>')
}

line6 <- list()
for (i in n) {
  line6[[i]] <- paste0('<operator id="narrow.t:', data_id[[i]], '" spec="Exchange" tree="@Tree.t:', data_id[[i]], '" weight="15.0"/>')
}

line7 <- list()
for (i in n) {
  line7[[i]] <- paste0('<operator id="wide.t:', data_id[[i]], '" spec="Exchange" isNarrow="false" tree="@Tree.t:', data_id[[i]], '" weight="3.0"/>')
}

line8 <- list()
for (i in n) {
  line8[[i]] <- paste0('<operator id="WilsonBalding.t:', data_id[[i]], '" spec="WilsonBalding" tree="@Tree.t:', data_id[[i]], '" weight="3.0"/>')
}

line9 <- list()
for (i in n) {
  line9[[i]] <- paste0('<operator id="updown.', data_id[[i]], '" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">')
}

line10 <- list()
for (i in n) {
  line10[[i]] <- paste0('<up idref="clockRate.c:', data_id[[i]], '"/>')
}

line11 <- list()
for (i in n) {
  line11[[i]] <- paste0('<down idref="Tree.t:', data_id[[i]], '"/>')
}

line12 <- as.list(rep("</operator>", 524))

clockscaler <- list()
for (i in n) {
  clockscaler[[i]] <- c(line1[[i]], line2[[i]], line3[[i]], line4[[i]],
                        line5[[i]], line6[[i]], line7[[i]], line8[[i]],
                        line9[[i]], line10[[i]], line11[[i]], line12[[i]])
}

fileConn <- file("clockscaler.txt")
writeLines(unlist(clockscaler), fileConn)
close(fileConn)


# 11.) updown

#  <operator id="strictClockUpDownOperator.c:C6432679" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">
#  <up idref="clockRate.c:C6432679"/>
#  <down idref="Tree.t:C6432679"/>
#  </operator>

line1 <- list()
for (i in n) {
  line1[[i]] <- paste0('<operator id="strictClockUpDownOperator.c:', data_id[[i]], '" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">')
}

line2 <- list()
for (i in n) {
  line2[[i]] <- paste0('<up idref="clockRate.c:', data_id[[i]], '"/>')
}

line3 <- list()
for (i in n) {
  line3[[i]] <- paste0('<down idref="Tree.t:', data_id[[i]], '"/>')
}

line4 <- as.list(rep("</operator>", 524))

operator <- list()
for (i in n) {
  operator[[i]] <- c(line1[[i]], line2[[i]], line3[[i]], line4[[i]])
}

fileConn <- file("operator.txt")
writeLines(unlist(operator), fileConn)
close(fileConn)


# 12.) likelogger

# <log idref="treeLikelihood.C6432679"/>
# <log id="TreeHeight.t:C6432679" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.t:C6432679"/>  
# <log idref="clockRate.c:C6432679"/>

line1 <- list()
for (i in n) {
  line1[[i]] <- paste0('<log idref="treeLikelihood.', data_id[[i]], '"/>')
}

line2 <- list()
for (i in n) {
  line2[[i]] <- paste0('<log id="TreeHeight.t:', data_id[[i]], 
                       '" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.t:',
                       data_id[[i]], '"/>')
}

line3 <- list()
for (i in n) {
  line3[[i]] <- paste0('<log idref="clockRate.c:', data_id[[i]], '"/>')
}


likelogger <- list()
for (i in n) {
  likelogger[[i]] <- c(line1[[i]], line2[[i]], line3[[i]])
}

fileConn <- file("likelogger.txt")
writeLines(unlist(likelogger), fileConn)
close(fileConn)

# 13.) logger

#  <logger id="treelog.t:C6432679" fileName="$(tree).trees" logEvery="5000" mode="tree">
#  <log id="TreeWithMetaDataLogger.t:C6432679" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:C6432679"/>
#  </logger>

line1 <- list()
for (i in n) {
  line1[[i]] <- paste0('<logger id="treelog.t:', data_id[[i]], '" fileName="$(tree).trees" logEvery="5000" mode="tree">')
}

line2 <- list()
for (i in n) {
  line2[[i]] <- paste0('<log id="TreeWithMetaDataLogger.t:', data_id[[i]], '" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:', data_id[[i]], '"/>')
}

line3 <- as.list(rep("</logger>", 524))

logger <- list()
for (i in n) {
  logger[[i]] <- c(line1[[i]], line2[[i]], line3[[i]])
}

fileConn <- file("logger.txt")
writeLines(unlist(logger), fileConn)
close(fileConn)
