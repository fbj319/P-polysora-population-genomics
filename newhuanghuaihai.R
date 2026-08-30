#####################################################################start analyse########
rm(list =ls())
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(png)
library(reshape2)
library(ggplot2) 
library(cowplot)

getwd()
setwd("D:/fbj/huanghuaihai")
vcfc <- read.vcfR("SCR_Huanghuaihai_qc.vcf.gz")

pop.data <- read.table("pop_group_huanghuaihai1.txt", sep = "\t", header = TRUE)
#colnames(vcf96@gt)[-1] <-pop.data$AccessID
colnames(vcfc@gt)[-1] == pop.data$AccessID


#remove ".fa" from vcf96_2@gt
#unlist(lapply(strsplit(colnames(vcf96_2@gt), split="\\."), function(x){x[1]}))
colnames(vcfc@gt) <- sub("\\.fa", "", colnames(vcfc@gt))

#rearrangement of txt
match(colnames(vcfc@gt)[-1], pop.data$AccessID)
pop.data <- pop.data[match(colnames(vcfc@gt)[-1], pop.data$AccessID),]
colnames(vcfc@gt)[-1] == pop.data$AccessID

#Converting the dataset to a genlight object. convert the data set into an object that is usable by poppr, adegenet, or any of the other population genetics packages in R. 
#The vcfR2genlight function subsets the data to filter loci that are not bi-allelic, returning an object that contains only loci with two alleles. The warning is to make sure we are aware that this action has taken place.
#vcfcgl <- vcfR2genlight(vcfc)
vcfcgl  <- vcfR2genind(vcfc)
#specify a ploidy of 2 for the entire population
ploidy(vcfcgl) <- 2



# 方法1：用read.delim（默认制表符分隔，更稳定）
pop.data <- read.delim("pop_group_Huanghuaihai1.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# 方法2：明确指定分隔符为制表符，并跳过可能的空行
pop.data <- read.table("pop_group_Huanghuaihai1.txt", sep = "\t", header = TRUE, 
                           +                        stringsAsFactors = FALSE, skip = 0,  # 若有1行空行，可设skip=1
                           +                        na.strings = c("", "NA"))  # 只将空字符串和"NA"视为缺
 # 检取结果
head(pop.data, 5)  # 查看前5行是否正确包含AH211-1、AH211-4等
#############################richness diversity evenness
#added the State column from our pop.data data frame to the pop slot of our genlight object
#pop(gl.rubi) <- pop.data$Country
pop(vcfcgl) <- pop.data$Province
pop(vcfcgl) <- pop.data$city
  

#######克隆矫正
##### # 嵌套分层：location（上层）包含city（下层）
###strata(vcfcgl) <- pop.data[, c("Province", "city")]  # 提取这两列作为分层
###rust <- clonecorrect(vcfcgl, strata = ~Province/city, keep = 1:2)
  
  #genotyptic richness
rustpop_diversity <- poppr(vcfcgl, total = T, sample = 999)
rustpop_diversity



#############################tree
library("ggplot2")
library("phangorn")
geo_tree <- aboot(x = vcfcgl, sample = 999, tree = "nj",  missing = "ignore", quiet = TRUE, showtree = FALSE, distance = bitwise.dist, cutoff = 75)
## branches according to Kuhner and Felsenstein (1994)
#spectral_pal <- spectral(length(levels(pop(vcfcgl))))
#cols <- spectral_pal[pop(rust)]
#svg("figures/nj_phylogeny.svg", width = 10, height = 10)

#draw the tree 
pop(vcfcgl) <- pop.data$Province
pop(vcfcgl) <- pop.data$city
#nursery_cols <- spectral_pal[pop(vcfcgl)]

tip_pal <-c("deeppink","seagreen2","gold")
nursery_cols <- tip_pal[pop(vcfcgl)]
#par(mar = c(0.1, 0.1, 0.1, 0.1))
#axis(side = 1, line = -2)
#axis(side = 2, line = -2)
#par(pin = c(8,8))
#par(mar = c(0.01, 0.01, 0.01, 0.01))

#q <- plot.phylo(ladderize(midpoint(geo_tree)), cex = 0.5, font = 2, adj = 0,
#                label.offset = 0.0125, underscore = TRUE, tip.color = nursery_cols,
#                type = "fan", align.tip.label = TRUE, open.angle = 10, x.lim = c(-0.6, 0.6), y.lim = c(-0.6, 0.6)) 

q <- plot.phylo(ladderize(midpoint(geo_tree)), cex = 0.5, font = 2, adj = 0,
                label.offset = 0.001, underscore = TRUE, tip.color = nursery_cols,
                type = "fan", align.tip.label = TRUE, open.angle = 10, x.lim = c(-0.04, 0.04), y.lim = c(-0.048, 0.048)) 

# draw boostraps
nodelabels(round(geo_tree$node.label, 1),  adj = c(1.3, -0.5),
           frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#axis(1)
#axis(2)
#axis(side = 1)
# add legend
#legend(-0.2,0.2, legend = levels(pop(rust)), pch = 22, pt.bg = spectral_pal, pt.cex = 2.5, title = "area")
#legend(-0.6,0.6, legend = levels(pop(vcfcgl)), pch = 22, pt.bg = spectral_pal, pt.cex = 2.5, title = "area", text.width = 0.01)# add scale
#legend(-0.056,0.05, legend = levels(pop(vcfcgl)),cex = 0.8, pch = 22, pt.bg = spectral_pal, pt.cex = 1.5, title = "Province", text.width = 0.005)# add scale
legend(-0.074,0.05, legend = levels(pop(vcfcgl)),cex = 0.8, pch = 22, pt.bg = tip_pal, pt.cex = 1.5, title = "Province", text.width = 0.005)# add scale

axisPhylo(xaxp= c(2, 9, 7))


pop(vcfcgl) <- pop.data$city
#tip_pal <- virid(length(levels(pop(vcfcgl))))

#tip_pal <-c("blue","cyan2","green","yellow","darkblue","red","purple")
#tip_pal <-c("red3","royalblue3","green4","purple4","tan2","yellow","tan4")
#tip_pal <-c("firebrick1","royalblue3","green4","purple3","darkorange","yellow","tan4")
tip_pal <-c("purple3","darkorange","yellow","tan4","red3","tan2","green","blue","cyan2")
tip_cols <- tip_pal[pop(vcfcgl)]
tiplabels(frame = "n", cex = 1, font = 3, xpd = NA, offset = 0.008, pch = 19, col = tip_cols) #0.07
#legend(-0.05,0.2, legend = levels(pop(rust)), pch = 22, pt.bg = tip_pal, pt.cex = 2.5, title = "month")
#legend(-0.72,0.6, legend = levels(pop(vcfcgl)), pch = 22, pt.bg = tip_pal, pt.cex = 2.5, title = "City", text.width = 0.02)
legend(-0.090,0.05, legend = levels(pop(vcfcgl)), cex = 0.8, pch = 22, pt.bg = tip_pal, pt.cex = 1.5, title = "City", text.width = 0.005)






###########################Minimum spanning networks#

library(igraph)
pop(vcfcgl) <- pop.data$Province
pop(vcfcgl) <- pop.data$city
vcfc.dist <- bitwise.dist(vcfcgl )
vcfc.msn <- poppr.msn(vcfcgl, vcfc.dist, showplot = FALSE, include.ties = T)
node.size <- rep(1, times = nInd(vcfcgl))
names(node.size) <- indNames(vcfcgl)
#vertex.attributes(rubi.msn$graph)$size <- node.size
myCol <- c("deeppink","seagreen2","gold")
set.seed(9)
plot_poppr_msn(vcfcgl, vcfc.msn , palette = brewer.pal(n = nPop(vcfcgl),
                                                       name = "Dark2"), gadj = 9, nodebase = 1.15, nodelab = 6, mlg = F, wscale  = F, inds = "xxx")

plot_poppr_msn(vcfcgl, vcfc.msn, inds = "none", palette = myCol, nodebase = 2.25, nodescale = 0.1, wscale = 0.1)




###############################################################DAPC_supply###################################################
#####BIC
###GBS data
library(vcfR)
getwd()
setwd("D:/fbj/huanghuaihai")
vcfc <- read.vcfR("SCR_Huanghuaihai_qc.vcf.gz")

pop.data <- read.table("pop_group_huanghuaihai1.txt", sep = "\t", header = TRUE)
#colnames(vcf96@gt)[-1] <-pop.data$AccessID
colnames(vcfc@gt)[-1] == pop.data$AccessID

#remove ".fa" from vcf96_2@gt
#unlist(lapply(strsplit(colnames(vcf96_2@gt), split="\\."), function(x){x[1]}))
colnames(vcfc@gt) <- sub("\\.fa", "", colnames(vcfc@gt))

#rearrangement of txt
match(colnames(vcfc@gt)[-1], pop.data$AccessID)
pop.data <- pop.data[match(colnames(vcfc@gt)[-1], pop.data$AccessID),]
colnames(vcfc@gt)[-1] == pop.data$AccessID

vcfcgl <- vcfR2genlight(vcfc)

library(adegenet)
maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(vcfcgl, n.pca = 40, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)
##   Group K      BIC
## 1     1 1 205.8632
## 2     2 1 205.8632
## 3     3 1 205.8632
## 4     4 1 205.8632
## 5     5 1 205.8632
## 6     6 1 205.8632
p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1











######PCA分析

#Principal components analysis
vcfc.pca <- glPca(vcfcgl, nf = 3)
barplot(vcfc.pca$eig, col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Proportion of variance explained")
title(xlab="Eigenvalue")


vcfc.pca.scores <- as.data.frame(vcfc.pca$scores)
vcfc.pca.scores$pop <- pop(vcfcgl)

#theme_bw, The classic dark-on-light ggplot2 theme. May work better for presentations displayed with a projector.
library(ggplot2)
set.seed(9)
p <- ggplot(vcfc.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()






#############交叉验证确定PCA值
#Discriminant analysis of principal components (DAPC)
#Cross validation: DAPC analysis
#The argument NA.method allows to replace missing data (NAs)
set.seed(999)
pop(vcfcgl) <- pop.data$city
pramx <- xvalDapc(tab(vcfcgl, NA.method = "mean"), pop(vcfcgl))
pramx[2:6]

######绘制DAPC
cols <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
          "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")  # 9种颜色
scatter(pramx$DAPC, col = cols, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

#note that scatter can also represent a single discriminant function, which is especially useful when only one of these has been 
#retained (e.g. in the case k = 2). This is achieved by plotting the densities of individuals on a given discriminant function with
#diferent colors for diferent groups:
scatter(pramx$DAPC,1,1, col=cols, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

####估算最佳遗传聚类数量 (K)
library("adegenet")
#table.value(table(pred.sup$assign, pop(x.sup)), col.lab=levels(pop(x.sup)))
grp <- find.clusters(vcfcgl, max.n.clust=40)
#####群体划分结果验证
table(pop(vcfcgl), grp$grp)
table.value(table(pop(vcfcgl), grp$grp), col.lab=paste("inf", 1:3),
            row.lab=paste("ori", 1:9))

####DAPC核心分析
dapc1 <- dapc(vcfcgl, grp$grp)
#######LD结果可视化
myCol <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")

# 绘制线性判别(LD)散点图
scatter(
  dapc1,
  col = myCol,        # 使用自定义颜色
  cex = 2,            # 点的大小
  pch = 17:22,        # 点形状(三角形到正方形)
  bg = "white",       # 背景色
  cstar = 0,          # 禁用群体星形连线
  scree.da = FALSE,   # 隐藏DA碎石图
  scree.pca = FALSE,  # 隐藏PCA碎石图
  legend = TRUE,      # 显示图例
  posi.leg = "bottomleft",  # 图例位置
  clabel = FALSE,     # 禁用数据点标签
  solid = 0.8,        # 填充不透明度
  cell = 0,           # 禁用椭圆
  posi.da = "bottomright"  # DA轴标签位置
)

# 添加坐标轴标题
title(xlab = "Linear Discriminant 1 (LD1)", 
      ylab = "Linear Discriminant 2 (LD2)")



#####DAPC结果可视化
scatter(dapc1)
myCol <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")
scatter(dapc1, col = myCol, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "bottomleft", cleg = 0.75)
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)
myCol <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")
scatter(dapc1, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft")

########DAPC核心分析（设置最佳PCA=30，K=3对应的判别函数数量）
mydapc <- dapc(vcfcgl, n.pca = 30, n.da = 2)
scatter(mydapc, col = myCol, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)
#mydapc <- dapc(vcfcgl, n.pca = 3, n.da = 2)
#mydapc <- dapc(vcfcgl, var.contrib = TRUE, scale = FALSE, n.pca = 3, n.da = 2)
mydapc <- dapc(vcfcgl, pop = pop(vcfcgl), parallel = T,  n.pca = 3, n.da = 2)
scatter(mydapc, col = myCol, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)




#######群体成员概率可视化
#####堆叠条形图 (ggplot2)
#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(mydapc,col = cols, posi = 'top')
# separate the samples by population
dapc.p <- as.data.frame(mydapc$posterior)
dapc.p$pop <- pop(vcfcgl)
dapc.p$indNames <- rownames(dapc.p)
library(reshape2)
dapc.p <- melt(dapc.p)
colnames(dapc.p) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
p <- ggplot(dapc.p, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

######密度图 (adegenet)
#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(dapc1,col = myCol, posi = 'top')
# separate the samples by population
dapc.p <- as.data.frame(dapc1$posterior)
dapc.p$pop <- pop(vcfcgl)
dapc.p$indNames <- rownames(dapc.p)
library(reshape2)
dapc.p <- melt(dapc.p)
colnames(dapc.p) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
p <- ggplot(dapc.p, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = myCol) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

#####群体遗传分化评估
#Population structure
library("mmod")
#Hendrick's standardized GSTGST to assess population structure

Gst_Hedrick(vcfc)










### Fst
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
#### Figure 3
library("hierfstat")
##Write the function
setwd("D:/fbj/huanghuaihai")
vcfc <- read.vcfR("SCR_Huanghuaihai_qc.vcf.gz")
pop.data <- read.table("pop_group_huanghuaihai1.txt", sep = "\t", header = TRUE)
colnames(vcfc@gt) <- sub("\\.fa", "", colnames(vcfc@gt))
match(colnames(vcfc@gt)[-1], pop.data$AccessID)
pop.data <- pop.data[match(colnames(vcfc@gt)[-1], pop.data$AccessID),]
colnames(vcfc@gt)[-1] == pop.data$AccessID
pop.data[1,]
vcfcgl  <- vcfR2genind(vcfc)
ploidy(vcfcgl) <- 2




# And change the names so we know what they are
strata(vcfcgl) <- data.frame(pop.data)[,-c(1,5)]
nameStrata(vcfcgl) <- ~Province/city
setPop(vcfcgl) <- ~Province/city
setPop(vcfcgl) <- ~Province
setPop(vcfcgl) <- ~city

pop(vcfcgl)

##############AMOVA_pbulish


vcfcgl
monpop <- as.genclone(vcfcgl)




# And change the names so we know what they are
nameStrata(monpop) <- ~Province/city

setPop(monpop) <- ~Province/city
setPop(monpop) <- ~Province
setPop(monpop) <- ~city



pop(monpop)
splitStrata(monpop) <- ~Province/city
setPop(monpop) <- ~city
amova.result <- poppr.amova(monpop, ~city) # gets a warning of zero distances
amova.result 
poppr.amova(monpop, ~Province/city)
poppr.amova(monpop, ~Province/city, filter = TRUE, threshold = 0.1) # no warning
amova.test <- randtest(amova.result, nrepet = 999) # Test for significance
plot(amova.test)
amova.test

setPop(monpop) <- ~Province
amova.result <- poppr.amova(monpop, ~Province) # gets a warning of zero distances
amova.result 
amova.test <- randtest(amova.result, nrepet = 999) # Test for significance
plot(amova.test)
amova.test


######################################Fst
pairwise.neifst(monpop,diploid=TRUE)

data(gtrunchier)
pairwise.neifst(gtrunchier[,-2],diploid=TRUE)





#############################################Fst_pbulish
library("hierfstat")
library("dplyr")
library("tidyr")
library("rlang")
library("ggplot2")
library("poppr")
library("adegenet")
vcfcgl
monpop <- as.genclone(vcfcgl)

# And change the names so we know what they are
nameStrata(monpop) <- ~Province/City/Month/City_month
setPop(monpop) <- ~Province/City/Month
setPop(monpop) <- ~Province
setPop(monpop) <- ~City
setPop(monpop) <- ~Month
setPop(monpop) <- ~City_month
pop(monpop)
splitStrata(monpop) <- ~Province/City/Month
setPop(monpop) <- ~City

rust <- monpop
setPop(rust) <- ~City_month
setPop(rust) <- ~Month
setPop(rust) <- ~Province
setPop(rust) <- ~City
setPop(rust) <- ~Group
pop(rust)
rust <- genind2hierfstat(rust )
x.fst <- pairwise.WCfst(rust)
x.fst 
pairwise.neifst(rust)

heatmap (x.fst, Rowv = NA, Colv = NA, symm = TRUE)

heatmap (x.fst, Rowv = NA, Colv = NA, symm = TRUE, name = "mat")
#heatmap (x.fst, Rowv = NA, Colv = NA)


library("ComplexHeatmap")
library("corrplot")

col3 <- colorRampPalette(c("red", "white", "blue"))
corrplot (x.fst, type = 'lower', order = 'hclust', col.lim=c(-0.005, 0.02), tl.col = 'black',
          cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))

corrplot (x.fst, type = 'lower', order = 'hclust', col.lim=c(-0.003, 0.0150),  col = col3(2000))

corrplot (x.fst, type = 'lower', order = 'hclust', method = 'color',col.lim=c(-0.003, 0.015), tl.col = 'black',
          cl.ratio = 0.2, tl.srt = 45)
corrplot(x.fst, method = 'color',type = 'lower', order = 'hclust',col.lim=c(-0.005, 0.02))

corrplot(x.fst,col.lim=c(-0.005, 0.02),  col = col3(2000))




corrplot.mixed(x.fst)



























setwd("D:/rust/publish_papers/4. population genetics in the pathogen's winter-reproduction regions")
vcfc <- read.vcfR("addreplicates_vcf_qualitycontrol_min2allelecount.vcf.gz")
pop.data <- read.table("winter_reproduction_locations-POP.txt", sep = "\t", header = TRUE)
colnames(vcfc@gt) <- sub("\\.fa", "", colnames(vcfc@gt))
match(colnames(vcfc@gt)[-1], pop.data$AccessID)
pop.data <- pop.data[match(colnames(vcfc@gt)[-1], pop.data$AccessID),]
colnames(vcfc@gt)[-1] == pop.data$AccessID
pop.data[1,]
vcfcgl  <- vcfR2genind(vcfc)
ploidy(vcfcgl) <- 2

# And change the names so we know what they are
strata(vcfcgl) <- data.frame(pop.data)[,-c(1,5)]
nameStrata(vcfcgl) <- ~Province/City/Month/Gro
nameStrata(vcfcgl) <- ~Province/City/Month/City_month
setPop(vcfcgl) <- ~Gro
setPop(vcfcgl) <- ~Province/City/Month/City_month
setPop(vcfcgl) <- ~Province/City/Month/Gro
setPop(vcfcgl) <- ~Province
setPop(vcfcgl) <- ~City
setPop(vcfcgl) <- ~Month
setPop(vcfcgl) <- ~City_month
pop(vcfcgl)

##############AMOVA_pbulish


vcfcgl
monpop <- as.genclone(vcfcgl)

# And change the names so we know what they are
nameStrata(monpop) <- ~Province/City/Month/City_month
nameStrata(monpop) <- ~Province/City/Month/Gro
setPop(monpop) <- ~Gro
setPop(monpop) <- ~Province/City/Month
setPop(monpop) <- ~Province
setPop(monpop) <- ~City
setPop(monpop) <- ~Month
setPop(monpop) <- ~City_month
pop(monpop)
splitStrata(monpop) <- ~Province/City/Month
setPop(monpop) <- ~City

vcfcgl
monpop <- as.genclone(vcfcgl)

# And change the names so we know what they are
nameStrata(monpop) <- ~Province/City/Month/City_month
setPop(monpop) <- ~Province/City/Month
setPop(monpop) <- ~Province
setPop(monpop) <- ~City
setPop(monpop) <- ~Month
setPop(monpop) <- ~City_month
pop(monpop)
splitStrata(monpop) <- ~Province/City/Month
setPop(monpop) <- ~City

rust <- monpop
setPop(rust) <- ~City_month
setPop(rust) <- ~Month
setPop(rust) <- ~Province
setPop(rust) <- ~City
setPop(rust) <- ~Group
pop(rust)
rust <- genind2hierfstat(rust )
x.fst <- pairwise.WCfst(rust)
x.fst 
pairwise.neifst(rust)

heatmap (x.fst, Rowv = NA, Colv = NA, symm = TRUE)

heatmap (x.fst, Rowv = NA, Colv = NA, symm = TRUE, name = "mat")
#heatmap (x.fst, Rowv = NA, Colv = NA)








































library("ggcorrplot")

ggcorrplot(x.fst, hc.order=TRUE,outline.color="white",
           type="lower",colors = c("#6D9EC1", "white", "#E46726"),
           ggtheme = ggplot2::theme_void())


ggcorrplot(x.fst)

ggcorrplot(x.fst, hc.order=TRUE,outline.color="white",
           type="lower")


library(gplots)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)#?????ÿ?????ɫ
hM <- format(round(x.fst, 3))#?????ݱ???2λС??

heatmap.2(x.fst,
          trace="none",#????ʾtrace
          #col=coul,#?޸???ͼ??ɫ
          density.info = "none",#ͼ??ȡ??density
          key.xlab ='Fst',
          key.title = "",
          cexRow = 1,cexCol = 1,#?޸ĺ???????????
          Rowv = F,Colv = F, #ȥ??????
          margins = c(6, 6),
          cellnote = hM,notecol='black'#????????ϵ????ֵ???޸???????ɫ
)



























attach(rust)
loci <- rust[,c(2:9)]
levels <- rust[,1]
varcomp.glob(levels, loci)
test.within(rust, test.lev = pop,  )











require(graphics); require(grDevices)
x  <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
hv <- heatmap(x, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              xlab = "specification variables", ylab =  "Car Models",
              main = "heatmap(<Mtcars data>, ..., scale = \"column\")")

utils::str(hv) # the two re-ordering index vectors

## no column dendrogram (nor reordering) at all:
heatmap(x, Colv = NA, col = cm.colors(256), scale = "column",
        RowSideColors = rc, margins = c(5,10),
        xlab = "specification variables", ylab =  "Car Models",
        main = "heatmap(<Mtcars data>, ..., scale = \"column\")")

## "no nothing"
heatmap(x, Rowv = NA, Colv = NA, scale = "column",
        main = "heatmap(*, NA, NA) ~= image(t(x))")


round(Ca <- cor(attitude), 2)
symnum(Ca) # simple graphic
heatmap(Ca,               symm = TRUE, margins = c(6,6)) # with reorder()
heatmap(Ca, Rowv = FALSE, symm = TRUE, margins = c(6,6)) # _NO_ reorder()
## slightly artificial with color bar, without and with ordering:
cc <- rainbow(nrow(Ca))
heatmap(Ca, Rowv = FALSE, symm = TRUE, RowSideColors = cc, ColSideColors = cc,
        margins = c(6,6))
heatmap(Ca,		symm = TRUE, RowSideColors = cc, ColSideColors = cc,
        margins = c(6,6))

## For variable clustering, rather use distance based on cor():
symnum( cU <- cor(USJudgeRatings) )

hU <- heatmap(cU, Rowv = FALSE, symm = TRUE, col = topo.colors(16),
              distfun = function(c) as.dist(1 - c), keep.dendro = TRUE)
## The Correlation matrix with same reordering:
round(100 * cU[hU[[1]], hU[[2]]])













data(humDNAm)
amovahum <- amova(humDNAm$samples, sqrt(humDNAm$distances), humDNAm$structures)
amovahum

data(Aeut)
strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
agc <- as.genclone(Aeut)
agc
amova.result <- poppr.amova(agc, ~Pop/Subpop)
amova.result
amova.test <- randtest(amova.result) # Test for significance
plot(amova.test)
amova.test

## Not run: 

# You can get the same results with the pegas implementation
amova.pegas <- poppr.amova(agc, ~Pop/Subpop, method = "pegas")
amova.pegas
amova.pegas$varcomp/sum(amova.pegas$varcomp)

# Clone correction is possible
amova.cc.result <- poppr.amova(agc, ~Pop/Subpop, clonecorrect = TRUE)
amova.cc.result
amova.cc.test <- randtest(amova.cc.result)
plot(amova.cc.test)
amova.cc.test


# Example with filtering
data(monpop)
splitStrata(monpop) <- ~Tree/Year/Symptom
poppr.amova(monpop, ~Symptom/Year) # gets a warning of zero distances
poppr.amova(monpop, ~Symptom/Year, filter = TRUE, threshold = 0.1) # no warning



## End(Not run)
