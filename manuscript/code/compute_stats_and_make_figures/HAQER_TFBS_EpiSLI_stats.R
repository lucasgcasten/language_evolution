### load packages
library(tidyverse)
library(Biostrings)
library(TFBSTools)
library(doMC)
library(JASPAR2020)

### prepare the PWMs for the TFs we will investigate
opts <- list()
opts[["species"]] <- 9606
opts[["all_versions"]] <- FALSE
pwm_set = getMatrixSet(JASPAR2020, opts)
pwm = lapply(pwm_set,function(x) toPWM(TFBSTools::Matrix(x),pseudocounts=colSums(TFBSTools::Matrix(x))))
names(pwm) = sapply(pwm_set,name)

### function for calculating delta TF motif scores
tfdiff = function(x1,x2,m){
	l = ncol(m)
	w = nchar(x1)[1]
	midpoint = ((w-1)/2)+1
	st = midpoint-l
	ed = midpoint+l
	library(genomation)
	s1 = patternMatrix(m,subseq(DNAStringSet(x1),start=st,end=ed),asPercentage=T,min.score=NULL)
	s1rc = patternMatrix(m,subseq(reverseComplement(DNAStringSet(x1)),start=st,end=ed),asPercentage=T,min.score=NULL)
	s2 = patternMatrix(m,subseq(DNAStringSet(x2),start=st,end=ed),asPercentage=T,min.score=NULL)
	s2rc = patternMatrix(m,subseq(reverseComplement(DNAStringSet(x2)),start=st,end=ed),asPercentage=T,min.score=NULL)

	s1max = pmax(apply(s1,1,max),apply(s1rc,1,max))
	s2max = pmax(apply(s2,1,max),apply(s2rc,1,max))

	return(list(refmax=s1max/100,altmax=s2max/100))
}


### read in the variants
har = readRDS("/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR_all_variants.hg19_normed.vcf.rds")
rand = readRDS("/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/RAND_all_variants.hg19_normed.vcf.rds")
haq = readRDS("/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER_v2_all_variants.hg19_normed.vcf.rds")

### fill in HCA alleles when it's obvious (human, neanderthal, chimp, and gorrilla all have same allele)
har$HCA_allele = case_when(is.na(har$HCA_allele) & har$neanderthal_allele %in% c(har$ref, har$alt) & har$neanderthal_allele == har$chimp_allele & har$neanderthal_allele == har$bonobo_allele & har$neanderthal_allele == har$gorilla_allele ~ har$neanderthal_allele,
                           TRUE ~ har$HCA_allele)
rand$HCA_allele = case_when(is.na(rand$HCA_allele) & rand$neanderthal_allele %in% c(rand$ref, rand$alt) & rand$neanderthal_allele == rand$chimp_allele & rand$neanderthal_allele == rand$bonobo_allele & rand$neanderthal_allele == rand$gorilla_allele ~ rand$neanderthal_allele,
                           TRUE ~ rand$HCA_allele)
haq$HCA_allele = case_when(is.na(haq$HCA_allele) & haq$neanderthal_allele %in% c(haq$ref, haq$alt) & haq$neanderthal_allele == haq$chimp_allele & haq$neanderthal_allele == haq$bonobo_allele & haq$neanderthal_allele == haq$gorilla_allele ~ haq$neanderthal_allele,
                           TRUE ~ haq$HCA_allele)

cols = list(rand=colnames(rand),har=colnames(har),haq=colnames(haq))
ctb = unclass(table(unlist(cols),rep(names(cols),sapply(cols,length))))

### a function to filter sites
fsites = function(x){
	x = as.data.frame(x)
	allele_cols = c("ref","alt","neanderthal_allele","chimp_allele",
	"bonobo_allele","gorilla_allele","orangutan_allele")
	alleles = as.data.frame(x[,allele_cols])
	alt_is = as.matrix(alleles[,-c(1:2)] == alleles[,2])
	ref_is = as.matrix(alleles[,-c(1:2)] == alleles[,1])
	mode(ref_is)="integer"
	mode(alt_is)="integer"

	### drop variants that are not within feature boundaries
	keep = x$distance_to_nearest==0
	### only rare variants
	keep = keep & x$max_AF < 0.01
	### drop variants for which we don't have an indicator of great ape alleles
	keep = keep & rowSums(is.na(ref_is)) ==0
	### drop variants where either ref or alt are not 1 in length (i.e., SNVs only)
	keep = keep & nchar(x$ref) ==1 & nchar(x$alt)==1

	### sample names 
	sample_names = c("sample1","sample10","sample100","sample101","sample102","sample103","sample104","sample105",
	"sample106","sample107","sample109","sample11","sample110","sample111","sample113","sample114",
	"sample115","sample116","sample117","sample118","sample119","sample12","sample120","sample121",
	"sample122","sample123","sample124","sample125","sample126","sample127","sample128","sample129",
	"sample13","sample130","sample131","sample132","sample133","sample134","sample135","sample136",
	"sample138","sample14","sample140","sample141","sample142","sample143","sample144","sample145",
	"sample146","sample147","sample148","sample149","sample15","sample150","sample151","sample152",
	"sample153","sample154","sample155","sample156","sample157","sample158","sample159","sample16",
	"sample160","sample161","sample163","sample164","sample165","sample166","sample168","sample169",
	"sample17","sample170","sample171","sample172","sample173","sample174","sample175","sample176",
	"sample177","sample178","sample179","sample18","sample180","sample181","sample182","sample183",
	"sample184","sample185","sample186","sample187","sample188","sample189","sample19","sample190",
	"sample2","sample20","sample201","sample202","sample203","sample204","sample205","sample206",
	"sample207","sample208","sample209","sample210","sample211","sample212","sample213","sample214",
	"sample215","sample217","sample218","sample219","sample220","sample221","sample222","sample224",
	"sample225","sample226","sample227","sample229","sample23","sample230","sample231","sample232",
	"sample233","sample234","sample235","sample236","sample237","sample238","sample239","sample240",
	"sample241","sample242","sample243","sample244","sample245","sample246","sample247","sample25",
	"sample259","sample26","sample260","sample261","sample262","sample263","sample264","sample265",
	"sample266","sample267","sample268","sample269","sample27","sample270","sample271","sample272",
	"sample273","sample274","sample275","sample276","sample277","sample279","sample28","sample280",
	"sample281","sample282","sample284","sample285","sample286","sample287","sample288","sample289",
	"sample29","sample291","sample292","sample293","sample294","sample295","sample298","sample299",
	"sample3","sample30","sample300","sample301","sample302","sample304","sample305","sample306",
	"sample307","sample309","sample31","sample310","sample311","sample312","sample313","sample315",
	"sample316","sample317","sample318","sample319","sample32","sample321","sample322","sample323",
	"sample324","sample325","sample326","sample327","sample328","sample329","sample33","sample330",
	"sample331","sample332","sample333","sample334","sample335","sample336","sample337","sample338",
	"sample339","sample34","sample340","sample341","sample342","sample343","sample344","sample347",
	"sample348","sample349","sample35","sample350","sample351","sample354","sample355","sample357",
	"sample358","sample359","sample36","sample360","sample361","sample362","sample363","sample364",
	"sample365","sample366","sample367","sample368","sample369","sample37","sample370","sample371",
	"sample374","sample375","sample376","sample377","sample378","sample379","sample38","sample380",
	"sample381","sample382","sample383","sample384","sample386","sample387","sample388","sample389",
	"sample39","sample390","sample391","sample392","sample393","sample397","sample4","sample40",
	"sample401","sample403","sample404","sample405","sample407","sample41","sample411","sample413",
	"sample414","sample415","sample416","sample417","sample42","sample43","sample44","sample45",
	"sample46","sample47","sample48","sample5","sample50","sample51","sample52","sample53","sample54",
	"sample55","sample56","sample57","sample58","sample59","sample6","sample60","sample61","sample62",
	"sample63","sample64","sample65","sample66","sample67","sample68","sample69","sample7","sample70",
	"sample71","sample72","sample73","sample74","sample75","sample76","sample77","sample79","sample8",
	"sample80","sample81","sample82","sample84","sample85","sample86","sample88","sample89","sample9",
	"sample90","sample91","sample92","sample93","sample94","sample95","sample96","sample97","sample99")

	### annotation column names
	ann_names = c("ID","chromosome","pos","ref","alt","rsid","gene","AF_EpiSLI","max_AF",
	"max_AF_pop","cadd","fathmm","impact","motif_name","motif_score_change",
	"HCA_allele","neanderthal_allele","chimp_allele",
	"bonobo_allele","gorilla_allele","orangutan_allele")

	### get genotype matrix
	gg = as.matrix(x[keep,sample_names])
	mode(gg) = "integer"
	gg[is.na(gg)] = 0
	rownames(gg) = x[keep,]$ID

	### organize annotation matrix
	ann = x[keep,ann_names]
	return(list(ann,gg,ref_is[keep,],alt_is[keep,]))
}

ll = lapply(list(rand,har,haq),fsites)
sapply(ll,function(x) dim(x[[1]]))

sources=c("random","HAR","HAQER")
for(i in 1:3){
	ll[[i]][[1]]$source = sources[i]
}

v = do.call('rbind',lapply(ll,function(x) x[[1]]))

##############################################################
### impute the reversion status where missing
##############################################################
library(glmnet)
rev = as.integer(v$alt == v$HCA_allele)
sum(rev, na.rm = TRUE)
sum(as.integer(v$ref == v$HCA_allele), na.rm = TRUE)
xrev = do.call('rbind',lapply(ll,function(x) cbind(x[[3]],x[[4]])))
xrev = cbind(xrev,as.matrix(v[,17:21]=="-"))

### weight each group equally
table(rev,v$source)
sum(table(rev,v$source))/6 ## ~ 391

wt = 391/table(as.factor(rev):as.factor(v$source))
wt = wt[as.character(as.factor(rev):as.factor(v$source))]

tr = !is.na(rev)
fit = cv.glmnet(y=rev[tr],x=xrev[tr,],alpha=0.9,family="binomial",weights=wt[tr]) # ,weights=wt[tr]

prd = predict(fit,xrev,s=fit$lambda.min,type="response")[,1]
rvr = as.integer(prd >= 0.95) ## classify reversions as all unlabeled variants with a predicted reversion probability >= 95%
rvr[tr] = rev[tr] ## relabel the training data according to their original labels

v$reversion = rvr
# v$reversion = prd
table(rvr, v$source)

##############################################################
### fetch the genomic sequence around the variant sites,
### insert the alternate base onto the ref background
##############################################################
library(GenomicRanges)
gr = v[,c(2,3,3,1)]
colnames(gr)[3] = "end"
colnames(gr)[2] = "start"
gr[[1]] = paste("chr",gr[[1]],sep='')
gr = as(gr,"GRanges")
max(sapply(pwm,ncol))
gr = resize(gr,width=51,fix="center")

library(BSgenome.Hsapiens.UCSC.hg19)
sref = getSeq(BSgenome.Hsapiens.UCSC.hg19,gr)
salt = as.character(sref)
substr(salt,26,26) = v$alt


### get the maximal motif scores for ref and alt sequences
### (only need to run once, this takes a while to run)
library(doMC)
registerDoMC(cores=10)
n = length(pwm)
print(n)
## takes a while to run - so commenting out and just reading in the results
s = foreach(i=1:n) %dopar% tfdiff(sref,salt,pwm[[i]])
names(s) = names(pwm)[1:n]
write_rds(s, '/wdata/lcasten/sli_wgs/HCA_reversion/data/TFBS_rare_variant_PWM_deltas_v2.rds')
# s <- read_rds('/wdata/lcasten/sli_wgs/HCA_reversion/data/TFBS_rare_variant_PWM_deltas.rds')
s <- read_rds('/wdata/lcasten/sli_wgs/HCA_reversion/data/TFBS_rare_variant_PWM_deltas_v2.rds')

score_ref = do.call('rbind',lapply(s,function(x) x[[1]]))
score_alt = do.call('rbind',lapply(s,function(x) x[[2]]))

### we see that increased motif match scores for certain TF families are
### favored along the arc of human evolution. Is individual variation
### in TFBS motif match scores related to core language ability (F1)?
load("/wdata/jmichaelson/SLI-seq/2024/Nature/phenotypes.Rdata.pxz", verbose = TRUE)
ph = p0[,2] ## 2nd col is F1 (core language)

### -log MAF (not used in final analysis)
# lmaf = -log10(pmax(v$max_AF,1e-5)) 
# lmaf = lmaf-min(lmaf)

### this order makes sense from a perspective, 
### but we flip the sign of the beta in a later analysis such that higher scores indicate increased motif integrity in the human ref allele
si = score_alt-score_ref

### this is just a variable that tracks 
### the sequence context type (HAQER, HAR, random)
seq_context = as.factor(v$source)
seq_context=relevel(seq_context,"random")


###################################################
### calculate per-motif burden scores and fit
### a linear model that associates with F1 core
### language scores
###################################################
gg = do.call('rbind',lapply(ll,function(x) x[[2]]))

## RAND HCA reversion burden per sample
si2 = scale(si)
brd_random = scale((t(gg * rvr * as.integer(v$source=="random")* ifelse(v$AF_EpiSLI < .05 & v$ref == v$neanderthal_allele, 1, 0)) %*% t(si2))) #  * ifelse(v$max_AF < .05, 1, 0))
cfl_random = t(sapply(1:ncol(brd_random), function(i) summary(lm(ph~brd_random[,i]))$coef[2,]))
rownames(cfl_random) = colnames(brd_random)
head(cfl_random[order(cfl_random[,3],decreasing=T),],20)

## HARs HCA reversion burden per sample 
# v$ref == v$neanderthal_allele & 
brd_har = scale(t(gg * rvr * as.integer(v$source=="HAR") * ifelse(v$AF_EpiSLI < .05 & v$ref == v$neanderthal_allele, 1, 0)) %*% t(si2)) #  * ifelse(v$AF_EpiSLI < .05, 1, 0))
cfl_har = t(sapply(1:ncol(brd_har), function(i) summary(lm(ph~brd_har[,i]))$coef[2,]))
rownames(cfl_har) = colnames(brd_har)
head(cfl_har[order(cfl_har[,3],decreasing=T),],20)

## HAQERs HCA reversion burden per sample
brd_haq = scale(t(gg * rvr * as.integer(v$source=="HAQER") * ifelse(v$AF_EpiSLI < .05 & v$ref == v$neanderthal_allele, 1, 0)) %*% t(si2)) #  * ifelse(v$AF_EpiSLI < .05, 1, 0))
cfl_haq = t(sapply(1:ncol(brd_haq), function(i) summary(lm(ph~brd_haq[,i]))$coef[2,]))
rownames(cfl_haq) = colnames(brd_haq)
head(cfl_haq[order(cfl_haq[,3],decreasing=T),],20)

### get the TF family annotations and arrange colors
### for plotting (this is gory code but works)
pr = read.table("manuscript/supplemental_materials/enrichment.InterPro.tsv",comment.char="",stringsAsFactors=F,header=T,sep="\t")
pl = strsplit(pr[[10]],split=",")
names(pl) = pr[[2]]
ptb = unclass(table(rep(names(pl),sapply(pl,length)),unlist(pl)))
ptb = ptb[c(1,3,8,15,19),]
lang = rownames(cfl_haq)
names(lang) = sapply(strsplit(lang,split="\\:|\\("),'[',1)
library(RColorBrewer)
set.seed(65666)
pal = sample(brewer.pal(name="Dark2",n=nrow(ptb)))
tfcols = structure(rep("grey",nrow(cfl_haq)),names=lang)
for(i in 1:5){
	tfcols[lang[colnames(ptb)[ptb[i,]==1]]] = pal[i]
	print(i)
}


##########################################################################
### here we fit models that account for avg. motif score per sequence
### source (e.g., HAQER, HAR, random); the coefficient of interest
### is the source-specific effect of reversions (the seq_context:rev term)
##########################################################################
betas = t(apply(si2,1,function(y) summary(lm(scale(y)~0+seq_context+seq_context:ifelse(rev == 1 & v$neanderthal_allele == v$ref, 1, 0)))$coef[,1]))
beta_se = t(apply(si2,1,function(y) summary(lm(scale(y)~0+seq_context+seq_context:ifelse(rev == 1 & v$neanderthal_allele == v$ref, 1, 0)))$coef[,2]))
pvals = t(apply(si2,1,function(y) summary(lm(scale(y)~0+seq_context+seq_context:ifelse(rev == 1 & v$neanderthal_allele == v$ref, 1, 0)))))

## look at human-gained TF binding for forkhead genes of interest
pvals[[which(rownames(si2) == 'FOXP2')]]

###########################################################################
### run the York regression; note that we flip the sign of the beta so that
### positive regression coefficients indicate better motif score with the 
### REF allele (which is the most common allele)
###########################################################################
###---------------------------------------------------------------
### York regression functions (uncertainty in X and Y directions)
###---------------------------------------------------------------
york_regression <- function(x, y, sx, sy, r = 0, tol = 1e-8, max_iter = 50) {
  # Input validation
  if (length(x) != length(y) || length(x) != length(sx) || length(x) != length(sy)) {
    stop("All input vectors must have the same length")
  }
  if (any(sx <= 0) || any(sy <= 0)) {
    stop("Standard errors must be positive")
  }
  if (abs(r) > 1) {
    stop("Correlation coefficient must be between -1 and 1")
  }
  
  n <- length(x)
  
  # Initial guess for slope using weighted mean of y/x and x/y regressions
  wx <- 1/sx^2
  wy <- 1/sy^2
  xbar <- sum(wx * x) / sum(wx)
  ybar <- sum(wy * y) / sum(wy)
  u <- x - xbar
  v <- y - ybar
  b <- sum(wy * u * v) / sum(wy * u^2)
  
  # Iterative solution
  b_old <- b + 2 * tol  # Ensure first iteration
  iter <- 0
  
  while (abs(b - b_old) > tol && iter < max_iter) {
    b_old <- b
    
    # Calculate weights
    w <- (sx^2 + b^2 * sy^2 - 2 * b * r * sx * sy)^(-1)
    
    # Calculate weighted means
    xbar <- sum(w * x) / sum(w)
    ybar <- sum(w * y) / sum(w)
    
    # Adjusted variables
    u <- x - xbar
    v <- y - ybar
    
    # Update slope
    b <- sum(w * u * v) / sum(w * u^2)
    
    iter <- iter + 1
  }
  
  if (iter >= max_iter) {
    warning("Maximum iterations reached without convergence")
  }
  
  # Calculate intercept
  a <- ybar - b * xbar
  
  # Calculate uncertainties
  x_adj <- x - xbar
  y_adj <- y - ybar
  
  # Standard errors of slope and intercept
  var_b <- 1 / sum(w * x_adj^2)
  se_b <- sqrt(var_b)
  se_a <- sqrt((1/sum(w)) + (xbar^2 * var_b))
  
  # Goodness of fit
  chi_sq <- sum(w * (y - (a + b * x))^2)
  df <- n - 2  # degrees of freedom
  p_value <- 1 - pchisq(chi_sq, df)
  
  # Return results
  list(
    intercept = a,
    slope = b,
    se_intercept = se_a,
    se_slope = se_b,
    iterations = iter,
    chi_square = chi_sq,
    p_value = p_value,
    xbar = xbar,
    ybar = ybar,
    n = n,
    weights = w
  )
}

add_york_line <- function(fit, xlim = NULL, n_points = 100, conf_level = 0.95, 
                           line_col = "blue", ci_col = rgb(0, 0, 1, 0.2)) {
  # If xlim not provided, use the range of observed x values +/- 10%
  if (is.null(xlim)) {
    x_range <- diff(range(x))
    xlim <- range(x) + c(-0.1, 0.1) * x_range
  }
  
  # Create sequence of x values for plotting
  x_seq <- seq(xlim[1], xlim[2], length.out = n_points)
  
  # Calculate predicted y values
  y_seq <- fit$intercept + fit$slope * x_seq
  
  # Calculate standard error of prediction
  # For each x point: var(a + bx) = var(a) + x²var(b) + 2x*cov(a,b)
  # Note: cov(a,b) = -xbar * var(b) for York regression
  var_pred <- fit$se_intercept^2 + 
              x_seq^2 * fit$se_slope^2 + 
              2 * x_seq * (-fit$xbar * fit$se_slope^2)
  
  # Calculate confidence intervals
  t_crit <- qt((1 + conf_level)/2, fit$n - 2)
  ci_width <- t_crit * sqrt(var_pred)
  
  # Add regression line
  lines(x_seq, y_seq, col = line_col, lwd = 2)
  
  # Add confidence intervals
  polygon(c(x_seq, rev(x_seq)), 
         c(y_seq + ci_width, rev(y_seq - ci_width)),
         col = ci_col, border = NA)
  
  # Return invisibly the computed values
  invisible(list(
    x = x_seq,
    y = y_seq,
    ci_lower = y_seq - ci_width,
    ci_upper = y_seq + ci_width
  ))
}

pt_error = function(x,y,sx,sy){
 arrows(x - sx, y, x + sx, y, code = 3, angle = 90, length = 0.05,col=rgb(0,0,0,0.2))
 arrows(x, y - sy, x, y + sy, code = 3, angle = 90, length = 0.05,col=rgb(0,0,0,0.2))
}

add_york_stats <- function(fit, x_pos = "topleft", y_pos = NULL) {
  # Create formatted strings using expression
  beta_text <- bquote(beta == .(format(round(fit$slope, 2), nsmall = 2)) %+-% 
                      .(format(round(fit$se_slope, 2), nsmall = 2)))
  
  chi_text <- bquote(chi^2 == .(format(round(fit$chi_square, 2), nsmall = 2)))
  
  # Format p-value in scientific notation if very small
  if (fit$p_value < 0.001) {
    # Extract base and exponent from scientific notation
    p_parts <- strsplit(format(fit$p_value, scientific = TRUE, digits = 2), "e")[[1]]
    base <- as.numeric(p_parts[1])
    exponent <- as.numeric(p_parts[2])
    p_text <- bquote(italic(p) == .(round(base, 2)) %*% 10^.(exponent))
  } else {
    p_text <- bquote(italic(p) == .(format(round(fit$p_value, 3), nsmall = 3)))
  }
  
  # Combine all expressions
  stats_expr <- as.expression(bquote(
    atop(
      .(beta_text),
      atop(.(chi_text), .(p_text))
    )
  ))
  
  # Add text to plot
  if (is.character(x_pos)) {
    legend(x_pos, legend = stats_expr, bty = "n")
  } else {
    legend(x_pos, y_pos, legend = stats_expr, bty = "n")
  }
}

##
cor.test(-1*betas[,4],cfl_random[,1])
cor.test(-1*betas[,5],cfl_haq[,1])
cor.test(-1*betas[,6],cfl_har[,1])

##
res_r = york_regression(-1*betas[,4],cfl_random[,1],beta_se[,4],cfl_random[,2])
res_haq = york_regression(-1*betas[,5],cfl_haq[,1],beta_se[,5],cfl_haq[,2])
res_har = york_regression(-1*betas[,6],cfl_har[,1],beta_se[,6],cfl_har[,2])

###########################################################################
### Make TFBS figures
###########################################################################
## function to add stats to left hand figures
add_york_stats_size <- function(fit, x_pos = "topleft", y_pos = NULL, size = 1) {
  # Create formatted strings using expression
  beta_text <- bquote(beta == .(format(round(fit$slope, 2), nsmall = 2)) %+-% 
                      .(format(round(fit$se_slope, 2), nsmall = 2)))
  
  chi_text <- bquote(chi^2 == .(format(round(fit$chi_square, 2), nsmall = 2)))
  
  # Format p-value in scientific notation if very small
  if (fit$p_value < 0.001) {
    # Extract base and exponent from scientific notation
    p_parts <- strsplit(format(fit$p_value, scientific = TRUE, digits = 2), "e")[[1]]
    base <- as.numeric(p_parts[1])
    exponent <- as.numeric(p_parts[2])
    p_text <- bquote(italic(p) == .(round(base, 2)) %*% 10^.(exponent))
  } else {
    p_text <- bquote(italic(p) == .(format(round(fit$p_value, 3), nsmall = 3)))
  }
  
  # Combine all expressions
  stats_expr <- as.expression(bquote(
    atop(
      .(beta_text),
      atop(.(chi_text), .(p_text))
    )
  ))
  
  # Add text to plot
  if (is.character(x_pos)) {
    legend(x_pos, legend = stats_expr, bty = "n", cex = size)
  } else {
    legend(x_pos, y_pos, legend = stats_expr, bty = "n", cex = size)
  }
}

## sequence class colors for consistency
haq_col = "#762776"
har_col = "#e04468"
rand_col = "#dcc699"
### the left figure panels
## plot set up
png("manuscript/figures/TFBS_left_panels.png",4,12,res=300,units="in")
par(mfrow=c(3,1))
xrg = c(-1.5,1.5)
yrg = c(-0.2,0.2)
pt_col = rgb(0,0,0,0.3)
cxl = 1.5

## add HAQER panel
plot(-1*betas[,5],cfl_haq[,1],
	xlim=xrg,ylim=yrg,cex=res_haq$weights/median(res_haq$weights),col=pt_col,cex.lab=cxl,
	xlab="Human/Neanderthal gained motif integrity",ylab="Language-motif integrity association",main="HAQERs", cex.main = 1.5, font.main = 2)
pt_error(-1*betas[,5],cfl_haq[,1],beta_se[,5],cfl_haq[,2])
add_york_line(res_haq,xlim=xrg,line_col=haq_col,ci_col=rgb(0.27,0.51,0.71,0.2))
add_york_stats_size(res_haq, size = 1.1)
## add HAR panel
plot(-1*betas[,6],cfl_har[,1],
	xlim=xrg,ylim=yrg,cex=res_har$weights/median(res_har$weights),col=pt_col,cex.lab=cxl,
	xlab="Human/Neanderthal gained motif integrity",ylab="Language-motif integrity association",main="HARs", cex.main = 1.5, font.main = 2)
pt_error(-1*betas[,6],cfl_har[,1],beta_se[,6],cfl_har[,2])
add_york_line(res_har,xlim=xrg,line_col='grey',ci_col=rgb(0.2,0.2,0.2,0.2))
add_york_stats_size(res_har, size = 1.1)
## add RAND panel
plot(-1*betas[,4],cfl_random[,1],
	xlim=xrg,ylim=yrg,cex=res_r$weights/median(res_r$weights),col=pt_col,cex.lab=cxl,
	xlab="Human/Neanderthal gained motif integrity",ylab="Language-motif integrity association",main="RAND", cex.main = 1.5, font.main = 2)
pt_error(-1*betas[,4],cfl_random[,1],beta_se[,4],cfl_random[,2])
add_york_line(res_r,xlim=xrg,line_col='grey',ci_col=rgb(0.2,0.2,0.2,0.2))
add_york_stats_size(res_r, size = 1.1)
dev.off()

### the center figure panel zooming in on the TF effects seen in HAQER sequences
## plot set up
png("manuscript/figures/TFBS_center_panel.png",7,7,res=300,units="in")
grp = unique(tfcols)
xrg = c(-0.55,0.4)
yrg = c(-0.2,0.17)
x1 = -1*betas[,5]
y1 = cfl_haq[,1]
sig = (-1*betas[,5]/beta_se[,5]) > 1.96 & cfl_haq[,3]>1.96
tf_map <- read.csv("manuscript/supplemental_materials/TFBS_data/TF_family_color_mapping.csv")
top_per_type <- data.frame(gene = names(sig), stat = -1*betas[,5]/beta_se[,5], sig = sig) %>% 
  inner_join(tf_map) %>% 
  arrange(desc(stat)) %>%
  filter(sig == TRUE) %>%
  group_by(tf) %>% 
  slice_head(n = 1)
sig = sig == TRUE & names(sig) %in% top_per_type$gene  
## make figure
plot(x1,y1,
	xlab="Human/Neanderthal gained motif integrity",ylab="Language-motif integrity association", main="HAQERs", cex.main = 1.5, font.main = 2,
	xlim=xrg,ylim=yrg,cex = 0,col=tfcols, cex.lab = 1.25, type = 'n',
	pch=ifelse(sig,16,1))
rect(xrg[1] - .0375, yrg[1] - .025, 0, yrg[2] + .025, col = "grey95", border = NA)  # Left half
rect(xrg[1]- .025, yrg[1] - .025, xrg[2] + .035, 0, col = "grey95", border = NA)  # Bottom-right
## make figure cex = (res_haq$weights/median(res_haq$weights)) / 2
points(x1,y1,
	cex = (res_haq$weights/median(res_haq$weights)) / 3,col=tfcols,
	pch=ifelse(sig,16,1))
xlab("Human gained motif integrity", cex.lab = 1.25)    
sig = sig | rownames(betas)=="FOXP2"
set.seed(8)
# abline(v=seq(-0.4,0.4,0.1),h=seq(-0.2,0.2,0.1),col=grey(0.8))
segments(c(0,0),c(0,0),c(1,0),c(0,1),lty=2,col=grey(0.15),lwd=2)
for(i in 1:6){
	hch = chull(x1[tfcols==grp[i]],y1[tfcols==grp[i]])
	idx = c(hch,hch[1])
	lines(x1[tfcols==grp[i]][idx],y1[tfcols==grp[i]][idx],col=grp[i],lwd=2)
}
txt = rownames(cfl_haq)[sig]
# txt[1] = "FOXP2"
text(x1[sig] + .01,y1[sig],txt,
	col=tfcols[sig],pos=sample(1:4,sum(sig),replace=T))
text(0.225,0.165,"Convergence of selection\nand language effects")
legend("topleft",legend=c("ETS domain","Forkhead","Homeobox","Basic Helix-Loop-Helix","Zinc Finger C2H2","other"),lty=1,col=c(pal,"grey"),bty='n',border=NA,lwd=3)
legend("bottomleft",legend=c("p < 0.05 for both positive language and positive selection effects","all others"),pch=c(16,1),col='grey',border=NA,bty='n',cex=0.8,pt.cex=1)
dev.off()


### the analysis of enrichment in quadrant I
fam = structure(c(rownames(ptb),"Other"),names=c(pal,"grey"))
fam[grp]

ft = lapply(grp,function(z) fisher.test(table(x1 > 0 & y1 > 0,tfcols%in%z) + 1))

names(ft) = fam[grp]
enr = data.frame(OR=log2(sapply(ft,function(x) x$est)),
	log2(t(sapply(ft,function(x) x$conf[1:2]))),
	p=sapply(ft,function(x) x$p.value))
rownames(enr) = names(ft)
enr = enr[order(enr[,1],decreasing=T),]

### bottom center panel - a barplot showing the enrichment results
## plot set up
enr_conv <- enr %>% 
	as.data.frame() %>% 
	rownames_to_column('TF_family') %>% 
	as_tibble() %>% 
	dplyr::rename(log2_OR = OR, ci_lower = X1, ci_upper = X2) %>% 
  mutate(TF_family_clean = case_when(TF_family == 'Homeobox domain' ~ 'Homeobox',
                               TF_family == 'Fork head domain' ~ 'Forkhead',
                               TF_family == 'Other' ~ 'Other',
                               TF_family == 'Myc-type, basic helix-loop-helix (bHLH) domain' ~ 'Basic\nHelix-Loop-Helix',
                               TF_family == 'Ets domain' ~ 'ETS domain',
                               TF_family == 'Zinc finger C2H2 superfamily' ~ 'Zinc Finger\nC2H2'))
png("manuscript/figures/TFBS_barplot.png",5,5,res=300,units="in")
par(mar=c(8,6,1,2))
col2 = apply(col2rgb(names(fam)[match(rownames(enr),fam)])/255,2,function(x) rgb(x[1],x[2],x[3],0.5))
col1 = names(fam)[match(rownames(enr),fam)]
bp = barplot(enr[,1],las=2,col=col2,
	ylab="Convergence of selection and\nlanguage effects (log2 OR)",
	border=NA,ylim=c(min(enr[,2]),max(enr[,3])),
	names.arg=enr_conv$TF_family_clean)
abline(h=0,col='grey',lty=2,lwd=2)
points(bp[,1],enr[,1],pch=ifelse(enr[,4]<0.05,16,1),cex=1.5,col=col1)
arrows(bp[,1], enr[,2], bp[,1], enr[,3], code = 3, angle = 90, length = 0.05,col=col1, lwd = 2.5)
ofst = c(0.5,0.5,-0.5,-0.7,-0.9,-1)
# legend("topright",legend=c("p > 0.05","p < 0.05"),pch=c(1,16),col='grey',border=NA,bty='n')
dev.off()


#################################################
## gather and save statistics
#################################################
## clean up TFBS reversion x F1 results
res_haq_clean <- res_haq %>% 
	as.data.frame() %>% 
	as_tibble() %>% 
	select(beta = slope, std.error = se_slope, iter = iterations, chi_square, p.value = p_value)  %>% 
	distinct() %>% 
	mutate(y = 'core_language_F1', x = 'HAQER_TFBS_selection') %>% 
	relocate(y, x)
res_har_clean <- res_har %>% 
	as.data.frame() %>% 
	as_tibble() %>% 
	select(beta = slope, std.error = se_slope, iter = iterations, chi_square, p.value = p_value)  %>% 
	distinct() %>% 
	mutate(y = 'core_language_F1', x = 'HAR_TFBS_selection') %>% 
	relocate(y, x)
res_rand_clean <- res_r %>% 
	as.data.frame() %>% 
	as_tibble() %>% 
	select(beta = slope, std.error = se_slope, iter = iterations, chi_square, p.value = p_value)  %>% 
	distinct() %>% 
	mutate(y = 'core_language_F1', x = 'RAND_TFBS_selection') %>% 
	relocate(y, x)
## merge TFBS selection data
bind_rows(res_haq_clean, res_har_clean, res_rand_clean) %>% 
	write_csv('manuscript/supplemental_materials/stats/TFBS_reversion_core_language_selection_results.csv')

## clean up enrichment of human-divergent selection and language results
enr_conv <- enr %>% 
	as.data.frame() %>% 
	rownames_to_column('TF_family') %>% 
	as_tibble() %>% 
	dplyr::rename(log2_OR = OR, ci_lower = X1, ci_upper = X2)
enr_conv %>% 
	write_csv('manuscript/supplemental_materials/stats/TFBS_TF_family_binding_convergence_results.csv')
