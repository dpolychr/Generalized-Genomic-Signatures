library(plyr);

# Comment: Read the "All.txt" file, which contains the accession numbers of the bacteria that are included in the final collection. Store the data in a variable, i.e. "bacteria".
bacteria <- read.table( file="bacteria_last.txt", header=T, sep="\t" );

# Comment: Set the variables where the results of the classification will be stored.
eval_j48 <- eval_LMT <- eval_SMO <- eval_Log <- phylum.class <- phyla <- list();

# Comment: "n" defines the odds (differences or ratios) of the oligonucleotides that are going to be used.
for (n in c( "DR_GGS_w2", "OR_GGS_w2", "DR_GGS_w3", "OR_GGS_w3", "DR_GGS_w4", "OR_GGS_w4", "DR_GGS_w5", "OR_GGS_w5" )) {

# n = "DR_GGS_w2"

# DATA MANIPULATION
# Comment: "..." is the path of the corresponding file, i.e. "DR_GGS_w2.txt".
odds <- read.table( file=paste("./", n, ".txt", sep=""), header=T, sep="\t" );
# Comment: Keep only the entries of the file (i.e. "DR_GGS_w2.txt") that correspond to the bacteria of the final collection (stored in "All.txt").
odds <- odds[ odds[,"accession"] %in% rownames(bacteria), ];
# Comment: Remove all duplicate rows that are present in the file (i.e. "DR_GGS_w2.txt").
odds <- odds[ !duplicated(odds[,"accession"]), ];
# Comment: Set the accession numbers of the remaining bacteria as names of the rows of the object "odds".
rownames(odds) <- odds[ ,"accession" ];

# Comment: Store the phylum/class of each entry represented in "odds" in a separate variable.
phylum.class[[n]] <- bacteria[ rownames(odds), "phylum.class" ];

# Convert to factor
phylum.class[[n]] <- as.factor(phylum.class[[n]])

# Comment: "mapvalues" is a function of the package "plyr": Do not use separate classes in the case of Proteobacteria. Instead, classify the bacteria according to their phylum.
# phyla[[n]] <- plyr::mapvalues(phylum.class[[n]], from=c("Alphaproteobacteria", "Betaproteobacteria", "Deltaproteobacteria", #"Epsilonproteobacteria", "Gammaproteobacteria"), to=rep("Proteobacteria", 5) );
# To keep all the classes
phyla[[n]] = phylum.class[[n]]

# Comment: Remove the "accession" column from the variable "odds". Accession numbers are now soted as "rownames" of "odds".
odds <- odds[, setdiff( colnames(odds), "accession") ];
# Comment: Replace "NAs", if any, with zeros.
odds[ is.na(odds) ] <- 0;


# RWeka CLASSIFICATION

library(RWeka);

# Comment: Create Weka classifier. "J48" generates unpruned or pruned C4.5 decision trees (Quinlan, 1993).
# odds_j48 <- J48(phyla[[n]]~ ., data=odds);
# Evaluate Weka classifier (10-fold validation).
# eval_j48[[n]] <- evaluate_Weka_classifier(odds_j48, numFolds=10, complexity=FALSE, seed=1, class=TRUE);

# Comment: Create Weka classifier. "LMT" implements ?Logistic Model Trees? (Landwehr, 2003; Landwehr et al., 2005).
 odds_LMT <- LMT(phyla[[n]]~ ., data=odds);
# Evaluate Weka classifier (10-fold validation).
 eval_LMT[[n]] <- evaluate_Weka_classifier(odds_LMT, numFolds=10, complexity=FALSE, seed=1, class=TRUE);

# Comment: Create Weka classifier. "SMO" implements John C. Platt?s sequential minimal optimization algorithm for training a support vector classifier using polynomial or RBF kernels.
# Evaluate Weka classifier (10-fold validation).
# odds_SMO <- SMO(phyla[[n]]~ ., data=odds);
# eval_SMO[[n]] <- evaluate_Weka_classifier(odds_SMO, numFolds=10, complexity=FALSE, seed=1, class=TRUE);

# Comment: Create Weka classifier. "Log" builds multinomial logistic regression models based on ridge estimation (le Cessie and van Houwelingen, 1992).
 odds_Log <- Logistic(phyla[[n]]~ ., data=odds);
# Evaluate Weka classifier (10-fold validation).
 eval_Log[[n]] <- evaluate_Weka_classifier(odds_Log, numFolds=10, complexity=FALSE, seed=1, class=TRUE);
}

#Save objects and then comment what you have run and run SMO

#  COMMENT: Each object with the name "eval_..." is a list (its names correspond to the odds, differences of ratios, of the oligonucleotides) with the results of the 10 Fold Cross Validation of the Weka classifiers. To retrieve the percentage of the correctly or incorrectly classified bacteria, among with other statistics, see "details" (i.e. eval_j48[[ "DR_GGS_w2" ]][[ "details" ]]). For the Confusion Matrices, see "confusionMatrix" (i.e. eval_j48[[ "DR_GGS_w2" ]][[ "confusionMatrix" ]])