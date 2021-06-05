#Step 62: Install bioconductor and other auxiliary packages
install.packages("BiocManager")
install.packages("tidyverse")
library("tidyverse")

#Step 63: Install and load the package TPP
BiocManager::install("TPP")
library("TPP")


#Step 64: Load the "proteinGroups.txt" file and extract the relevant columns.
proteinGroups <- read.delim("proteinGroups.txt")
proteinGroups <-  proteinGroups[!proteinGroups$Reverse=="+",]
proteinGroups <- proteinGroups[!proteinGroups$Potential.contaminant=="+",]
proteinGroups <-  proteinGroups[!proteinGroups$Only.identified.by.site=="+",]
proteinGroups <- select(proteinGroups, Protein.IDs, contains("Reporter.intensity.corrected"))
proteinGroups$qssm <- 6
proteinGroups$qupm <- 6
control_a <- select(proteinGroups, Protein.IDs, contains("control_a"), qupm, qssm)
control_b <- select(proteinGroups, Protein.IDs, contains("control_b"), qupm, qssm)
treated_a <- select(proteinGroups, Protein.IDs, contains("treated_a"), qupm, qssm)
treated_b <- select(proteinGroups, Protein.IDs, contains("treated_b"), qupm, qssm)
colnames(control_a) <- c("gene_name",	"rel_fc_131L", "rel_fc_130H", "rel_fc_130L",	"rel_fc_129H",	"rel_fc_129L",	"rel_fc_128H",	"rel_fc_128L",	"rel_fc_127H",	"rel_fc_127L",	"rel_fc_126", "qssm", "qupm")
colnames(control_b) <- c("gene_name",	"rel_fc_131L", "rel_fc_130H", "rel_fc_130L",	"rel_fc_129H",	"rel_fc_129L",	"rel_fc_128H",	"rel_fc_128L",	"rel_fc_127H",	"rel_fc_127L",	"rel_fc_126", "qssm", "qupm")
colnames(treated_a) <- c("gene_name",	"rel_fc_131L", "rel_fc_130H", "rel_fc_130L",	"rel_fc_129H",	"rel_fc_129L",	"rel_fc_128H",	"rel_fc_128L",	"rel_fc_127H",	"rel_fc_127L",	"rel_fc_126", "qssm", "qupm")
colnames(treated_b) <- c("gene_name",	"rel_fc_131L", "rel_fc_130H", "rel_fc_130L",	"rel_fc_129H",	"rel_fc_129L",	"rel_fc_128H",	"rel_fc_128L",	"rel_fc_127H",	"rel_fc_127L",	"rel_fc_126", "qssm", "qupm")


#Step 65: Plot the protein abundance per temperature in a boxplot.
windowsFonts(A = windowsFont("Arial"))
par(family = "A")
#Control A
boxplot(control_a[2:11], ylim=c(0, 1.2), ylab="Protein relative abundance", xlab="Temperature (°C)", names=c(33,37,41,45,49,53,57,61,65,69))
#Control B
boxplot(control_b[2:11], ylim=c(0, 1.2), ylab="Protein relative abundance", xlab="Temperature (°C)", names=c(33,37,41,45,49,53,57,61,65,69))
#Treated A
boxplot(treated_a[2:11], ylim=c(0, 1.2), ylab="Protein relative abundance", xlab="Temperature (°C)", names=c(33,37,41,45,49,53,57,61,65,69))
#Treated B
boxplot(treated_b[2:11], ylim=c(0, 1.2), ylab="Protein relative abundance", xlab="Temperature (°C)", names=c(33,37,41,45,49,53,57,61,65,69))

#Step 66
n <- length(proteinGroups$Protein.IDs)
p_coverage <- n/8023*100
print(paste("The protein coverage is",p_coverage,"%"))

#Step 67: Load the experiment files in the environment.
table_Control_1  <- control_a
Control_1 = as.data.frame.matrix(table_Control_1) 
table_Control_2  <- control_b
Control_2 = as.data.frame.matrix(table_Control_2) 
table_drug_1  <- treated_a
Drug_1 = as.data.frame.matrix(table_drug_1) 
table_drug_2  <- treated_b
Drug_2 = as.data.frame.matrix(table_drug_2) 
hdacTR_data <- list(Control_1 = Control_1, Control_2 = Control_2, Drug_1=Drug_1, Drug_2=Drug_2)
resultPath = file.path(getwd(), 'analysis')

#Step 68: Load the temperature setup
nnhdacTR_config = read.csv("TPP_config.txt", sep=' ', header=TRUE, check.names=FALSE)

#Step 69: Start the workflow 
TPPresults <- analyzeTPPTR(configTable=nnhdacTR_config, 
                          data=hdacTR_data, 
                          normalize = FALSE, 
                          nCores=2, 
                          resultPath=resultPath,
                          plotCurves=TRUE,
                          pValFilter = list(minR2 = 0.65, maxPlateau = 0.4))

#Step 70: Generate a list of Tm targets
Tm_targets <- filter(TPPresults, fulfills_all_4_requirements==TRUE)
write.csv(Tm_targets, "Tm_targets.csv")

#Step 71: Generate a list of top NPARC targets
NPARC_targets <- filter(TPPresults, p_adj_NPARC<0.01)
write.csv(NPARC_targets, "NPARC_targets.csv")


