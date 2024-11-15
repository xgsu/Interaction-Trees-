
####################################################
####################################################
# PART I: ILLUSTRATION: ONE SINGLE IT TREE ANALYSIS
####################################################
####################################################

# rm(list=ls(all=TRUE))
source("Functions-IT.R")

# ====================
# GENERATE DATA
# ====================

set.seed(123)
beta0 <- rep(2, 6)
dat <- rdat(n=600, beta=beta0, p= 3, digits=1) 
dim(dat); head(dat)

# =============================================================
# CONSTRUCT ONE SINGLE IT VIA TEST SAMPLE METHOD (METHOD I) 
# =============================================================

# GENERATE TEST SAMPLE DATA
test.dat <- rdat(n=400, beta=beta0, p= 3, digits=1)  

# OBTAIN A LARGE INITIAL TREE
split.var <- 3:6  # COLUMNS OF COVARIATES
cols.ctg <- c(6)  # COLUMNS OF CATEGORICAL COVARIATES
tree0 <- grow.INT(data=dat, test=test.dat, min.ndsz=20, n0=5, 
	split.var=split.var, ctg=cols.ctg, max.depth=5)
# PLOT THE INITIAL TREE
plot.tree(tree=tree0, textDepth=10, leaf.label=FALSE, lines="rectangle")

# PRUNE AND SELECT BEST TREE SIZE
prn <- prune.size.testsample(tree0)
# PRUNING INFORMATION, LOOK FOR THE MAXIMUM Ga.2, Ga.3, Ga.4, OR Ga.BIC
prn$result 

# PLOT OF G.a VS. TREE SIZE (# TERMINALS)
tmp <- prn$result 
col.Ga <- 7:10
for (j in c(4, col.Ga)) tmp[,j] <- as.numeric(as.character(tmp[,j]))
y.min <- min(tmp[, col.Ga]); y.max <- max(tmp[, 7:10])
x.min <- min(tmp[, 4]); x.max <- max(tmp[, 4])
plot(c(x.min, x.max), c(y.min, y.max), type="n", xlab="tree size", ylab="G.a Values")
for (j in col.Ga) lines(tmp$size.tmnl, tmp[,j], col=j, lty=1)
legend(x.max-5, y.max, col=col.Ga, lty=1, legend=colnames(tmp)[col.Ga])    


# THE BEST TREE SIZES (USING BIC)
prn$size
bsize <- prn$size[4]; bsize 
# OBTAIN THE FINAL TREE STRUCTURE
btree <- obtain.btree(tree0, bsize=bsize); # btree

# PLOT THE TREE STRUCTURE
source("Functions-IT.R")
plot.tree(tree=btree, textDepth=10, 
          leaf.label=TRUE, lines="rectangle")
title(main="The Final IT Structure")


# GENERATE LaTeX CODES FOR PACKAGE {PCTricks}
plot.tree.latex(btree, file="tree-code.tex", digits=5)



# =============================================================
# CONSTRUCT ONE SINGLE IT VIA BOOTSTRAP METHOD (METHOD II) 
# =============================================================

# GROW AND PRUNE
boot.result <- boottrap.grow.prune(B=25, data=dat, N0=20, n0=5, 
	split.var=split.var, ctg=cols.ctg, 
	max.depth=8, LeBlanc=T, min.boot.tree.size=1)  
# THE INITIAL LARGE TREE CAN BE FOUND FROM TEH RESULTS
tree0 <- boot.result$initial.tree; tree0; 
# DETERMINE THE BEST TREE SIZE AND OBTAIN THE BEST TREE MODELS
boot.prune <- boot.result$boot.prune
OUT <- bootstrap.size(boot.prune, n=nrow(dat), 
	plot.it=TRUE, filename=NULL, horizontal=T);
names(OUT)
OUT$G.a 	# BEST G.a 
OUT$bsize	# BEST TREE SIZES
OUT$btree  	# BEST TREE MODELS

# BEST TREE USING BIC
btree <- OUT$btree$`G.log(n)`;  btree
plot.tree(tree=btree, textDepth=6, lines="rectangle")

# IF WE WANT A TREE MODEL WITH ARBITRARY NUMBER OF TERMINAL NODES 
btree <- obtain.btree(tree0, bsize=4); btree


# =======================
# AMALGAMATION
# =======================

n.groups <- 3
result.amalgamation <- amalagamation(dat, btree, char.var=cols.ctg, n.groups=n.groups, details=F)
result.amalgamation$merge.info
dat0 <- result.amalgamation$data; dat0$node 

# =============================
# SUMMARIZE THE FINAL SUBGROUPS 
# =============================

dat$subgroup <- dat0$node
col.grp <- which(colnames(dat)=="subgroup")
summary.grp <- summary.subgroups(dat, var.equal=T, digits=4, col.grp=col.grp)
summary.grp







#################################
#################################
# PART II: RANDOM FORESTS OF IT 
#################################
#################################

source("Functions-randomForests-IT.R")

# =====================================
# BUILT A NUMBER OF RANDOM FORESTS
# =====================================

names(dat)
# GROW AND PRUNE
# ctg <- NA  #### DONOT USE "NULL" IF NO NOMINAL COVARIATES
col.y <- 1; col.trt <- 2
ntree <- 200
RF.fit1 <- Build.RF.IT(dat=dat, col.y=col.y, col.trt=col.trt, split.var=split.var, ctg=cols.ctg, 
	N0=20, n0=5,  max.depth=10,
	ntree=ntree, mtry = max(floor(length(split.var)/3), 1),
	avoid.nul.tree=T)
ntree <- 200
RF.fit2 <- Build.RF.IT(dat=dat, col.y=col.y, col.trt=col.trt, split.var=split.var, ctg=cols.ctg, 
	N0=20, n0=5,  max.depth=10,
	ntree=ntree, mtry = max(floor(length(split.var)/3), 1),
	avoid.nul.tree=T)

# COMBINE TWO RANDOM FORESTS
RF.fit <- combine.RFIT(RF.fit1, RF.fit2)
names(RF.fit)

# =========================
# VARIABLE IMPORTANCE
# =========================

VI <- Variable.Importance(RF.fit, n0=5, sort=T, details=F, truncate.zeros=T); VI

# BAR PLOT OF VARIABLE IMPORTANCE 
plot.VI(VI, filename="fig-VI.eps", horizontal=T, rain.bow=F)
plot.VI(VI, filename=NULL, horizontal=T)



################## THE END ####################### 
