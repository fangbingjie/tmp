# 生信多元统计分析（R语言）

[TOC]

## 一、线性代数基础-矩阵

```r
rm(list=ls())

#安装R包
library("BiocManager")
BiocManager::install("pheatmap")
BiocManager::install("ggplot2")
BiocManager::install("ggtree")

setwd("C:/Users/rexki/.r_data")

a <- c(1,2,3,4)
b <- c(4,3,2,1)

#向量点乘
crossprod(a,b)

#创建一个测试矩阵test_mat
test_mat <- matrix(c(1,1,3,1,2,1,3,1,1),nrow = 3,ncol = 3)

#求矩阵的逆
test_mat_inverse <- solve(test_mat)

#单位矩阵=矩阵*矩阵逆
test_mat_identity <- round(test_mat %*% test_mat_inverse)

#矩阵特征值和特征向量
test_mat_eigen <- eigen(test_mat)
λ <- test_mat_eigen$values
t_mat <- test_mat_eigen$vectors

#矩阵谱分解 t*diag(λ)*t转置或逆
t_mat %*% diag(λ) %*% solve(t_mat)
t_mat %*% diag(λ) %*% t(t_mat)

#GTEx
load("GTEx.RData")

GTEx_mat <- GTEx.TPM.gene.mat[,c(-1,-2)]
dim(GTEx_mat)

#协方差矩阵
GTEx_mat_cov <- cov(GTEx_mat)

#相关系数矩阵
GTEx_mat_cor <- cor(GTEx_mat)

#热图查看相关系数
library(pheatmap)
pheatmap(GTEx_mat_cor)
```

## 二、主成分分析-PCA

```r
rm(list=ls())

#安装R包
library("BiocManager")
BiocManager::install("FactoMineR")
BiocManager::install("corrplot")
BiocManager::install("factoextra")
BiocManager::install("backports")
BiocManager::install("data.table")
BiocManager::install("RColorBrewer")

#加载iris数据
head(iris)
iris_mat <- iris[,-5]

#使用PCA包，得到特征值
library("FactoMineR")
iris_pca <- FactoMineR::PCA(iris_mat,graph = F)
iris_pca$eig
iris_pca$ind$coord[,1]

########################  手动运算PCA  #########################
#scale归一化
iris_mat_scale <- scale(iris_mat,center = T)
#得到协方差矩阵
iris_mat_cov <- cov(iris_mat_scale)
#得到特征值和特向量征向量
iris_mat_cov_eigen <- eigen(iris_mat_cov)
iris_mat_cov_eigen$values
iris_mat_cov_eigen$vectors

pc1 <- iris_mat_scale %*% iris_mat_cov_eigen$vectors[,1]
pc2 <- iris_mat_scale %*% iris_mat_cov_eigen$vectors[,2]

####对比包算的和手动算的,在一根直线上即相等
pc1
iris_pca$ind$coord[,1]

plot(x =pc1, y=iris_pca$ind$coord[,1] )

##########################   PCA画图   #########################
rm(list=ls())

# --------------------------------------------------------->>>>>>>>>>>>>
# 1.PCA part
# --------------------------------------------------------->>>>>>>>>>>>>
library(corrplot)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)

#去掉第5列
iris_mat <- iris[,-5]

#plot correlation  查看变量间是否有相关关系
corrplot(cor(iris_mat),diag = F,type = "lower")

#run pca
iris_pca <- FactoMineR::PCA(iris_mat,graph = F)

########################   PCA 贡献率   #######################
fviz_eig(iris_pca,addlabels = T)

########################    PCA 点图    #######################
# version 1
fviz_pca_ind(iris_pca,pointsize="coord",pointshape=21,fill=iris$Species)

# version 2  
#选3个颜色
RColorBrewer::display.brewer.all()
color_list <- RColorBrewer::brewer.pal(3,"Set1")

#加阴影，加椭圆，加调色板
fviz_pca_ind(iris_pca,
             label="none",
             habillage = iris$Species,
             palette = color_list,
             addEllipses = TRUE)

#########################   PCA 载荷图    #####################
fviz_pca_var(iris_pca,col.var = "coord")

# plot PC1
fviz_contrib(iris_pca, choice = "var", axes = 1)

# plot PC2
fviz_contrib(iris_pca, choice = "var", axes = 2)

#########################  使用GTEx数据 PCA  ######################
#GTEx
rm(list = ls())

library(corrplot)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)

#加载数据
load("GTEx.Rdata")

accession_df <- read.csv(file = "GTEx_accession_table.csv")
annotation_df <- read.csv(file = "GTEx_v7_Annotations_SampleAttributesDS.txt",header = T,sep = "\t")

#得到GTEx 矩阵，进行log2 +1,行加基因名
GTEx_mat <- GTEx.TPM.gene.mat[,c(-1,-2)]
GTEx_mat_log2 <- log2(GTEx.TPM.gene.mat[,c(-1,-2)] + 1)
rownames(GTEx_mat_log2) <- as.character(GTEx.TPM.gene.mat[,1])

#查看整体分布情况
hist(as.matrix(GTEx_mat_log2))

#计算 gene 标准差 
SD_value <- apply(GTEx_mat,1,FUN = function(x){sd(x)})

#筛选变化最大的基因top 500或1000,2000,3000
GTEx_mat_log2_sort <- GTEx_mat_log2[order(SD_value,decreasing = T),]
GTEx_mat_log2_sort_top <- GTEx_mat_log2_sort[1:1000,]

# get index
col_index = limma::strsplit2(colnames(GTEx_mat_log2_sort_top),split = "\\.")[,5]

# 增加组织信息,转换基因名.为-
accession_vec <- gsub(x = colnames(GTEx_mat_log2_sort_top),pattern = "\\.",replacement = "-") 
tissue_info <- as.character(annotation_df$SMTS[match(accession_vec,annotation_df$SAMPID)])

#增加tissue一列
new_colname <- sprintf("%s_%s",tissue_info, col_index)
GTEx_mat_log2_sort_top_fix <- GTEx_mat_log2_sort_top
colnames(GTEx_mat_log2_sort_top_fix) = new_colname

# show cor 展示相关关系
# uniqe tissue
GTE_mat_log2_uniqe_tissue <- GTEx_mat_log2_sort_top_fix[,!duplicated(tissue_info)]
colnames(GTE_mat_log2_uniqe_tissue) <- tissue_info[!duplicated(tissue_info)]
corrplot(cor(GTE_mat_log2_uniqe_tissue),diag = F,type = "lower")

# all case
GTE_mat_log2_all_case <- GTEx_mat_log2_sort_top_fix
colnames(GTE_mat_log2_all_case) <- sprintf("%s.%d",tissue_info,1:length(tissue_info))
corrplot(cor(GTE_mat_log2_all_case),diag = F,type = "lower")

#run PCA
GTEx_pca <- FactoMineR::PCA(t(GTE_mat_log2_all_case),graph = F,scale.unit = T)
GTEx_pca_uniqe <- FactoMineR::PCA(t(GTE_mat_log2_uniqe_tissue),graph = F,scale.unit = T)

# 特征值大小
GTEx_pca$eig

# var info
summary(GTEx_pca$var)

##################### PCA 贡献率 #####################
fviz_eig(GTEx_pca,addlabels = T)

##################### PCA 载荷图 #####################
fviz_pca_var(GTEx_pca, col.var = "cos2")

# plot PC1
fviz_contrib(GTEx_pca, choice = "var", axes = 1)

# plot PC2
fviz_contrib(GTEx_pca, choice = "var", axes = 2)
##################### PCA 点图 #######################
#### plot directly
fviz_pca_ind(GTEx_pca,pointsize="cos2",pointshape=21, fill=tissue_info)

fviz_pca_ind(GTEx_pca,pointshape=21, fill=tissue_info)

#### plot with shadow
#选8个颜色
RColorBrewer::display.brewer.all()
color_list = colorRampPalette(brewer.pal(8,"Set2"))(length(unique(tissue_info)))

#加阴影，加椭圆，加调色板
fviz_pca_ind(GTEx_pca,
             label = "none", # hide individual labels
             habillage = factor(tissue_info), # color by groups
             palette = color_list,
             addEllipses = TRUE # Concentration ellipses
)
```

![](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107102358112.png)

![image-20210107102425688](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107102425688.png)

![image-20210107102450592](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107102450592.png)

![image-20210107102522137](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107102522137.png)

![image-20210107102607559](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107102607559.png)

## 三、**层次聚类与k-means聚类**

```r
rm(list=ls())


# ---------------------------------------------------------------------------->>>>>>>>
# 0.make data
# ---------------------------------------------------------------------------->>>>>>>>

setwd("C:/Users/rexki/.r_data")

#安装R包
library("BiocManager")
BiocManager::install("fpc")
BiocManager::install("limma")
BiocManager::install("modeltools")

# ---------------------------------------------------------------------------->>>>>>>>
# 1.hcluster 层次聚类
# ---------------------------------------------------------------------------->>>>>>>>

rm(list=ls())

load("GTEx.RData")

# tissue的注释文件
accession_df <- read.csv(file="GTEx_accession_table.csv")
annotation_df <- read.csv(file="GTEx_v7_Annotations_SampleAttributesDS.txt",header = T,sep = "\t")

# make matrix  得到GTEx 矩阵，进行log2 +1,行加基因名
GTEx_mat <- GTEx.TPM.gene.mat[,c(-1,-2)]
GTEx_mat_log2 <- log2(GTEx_mat + 1)
rownames(GTEx_mat_log2) <- as.character(GTEx.TPM.gene.mat$Name)

#计算 gene var 
SD_value <- apply(GTEx_mat,1,FUN = function(x){sd(x)})

#筛选变化最大的基因top 500或1000,2000,3000
GTEx_mat_log2_sort <- GTEx_mat_log2[order(SD_value,decreasing = T),]
GTEx_mat_log2_sort_top <- GTEx_mat_log2_sort[1:1000,]

# get index
col_index = limma::strsplit2(colnames(GTEx_mat_log2_sort_top),split = "\\.")[,5]

# 增加tissue信息,转换基因名.为-
accession_vec <- gsub(x = colnames(GTEx_mat_log2_sort_top),pattern = "\\.",replacement = "-") 
tissue_info <- as.character(annotation_df$SMTS[match(accession_vec,annotation_df$SAMPID)])

#增加tissue一列
new_colname <- sprintf("%s_%s",tissue_info, col_index)
GTEx_mat_log2_sort_top_fix <- GTEx_mat_log2_sort_top
colnames(GTEx_mat_log2_sort_top_fix) <- new_colname

# run hclust
dist_mat <- dist(t(GTEx_mat_log2_sort_top_fix))
plot(hclust(dist_mat))

# plot with cor 
library(pheatmap)
pheatmap(cor(GTEx_mat_log2_sort_top_fix))


# ---------------------------------------------------------------------------->>>>>>>>
# 2.kernel clustring
# ---------------------------------------------------------------------------->>>>>>>>
library(fpc)

#加载数据
iris_mat <- iris[,-5]

# set seed
set.seed(2021)

# k-means
par(mfrow=c(2,3))
for(k_num in c(2,3,4)){
  kmeans_res = kmeans(iris_mat, k_num, iter.max = 10000)
  plot(iris_mat, pch=kmeans_res$cluter, col=kmeans_res$cluster)
  title(paste("k-means: k=", k_num))
}
```

![image-20210107102153641](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107102153641.png)

![image-20210107102232505](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107102232505.png)

## 四、多维标度法（MDS）

```r
rm(list = ls())

# ---------------------------------------------------------------------------->>>>>>>>
# 0.make data
# ---------------------------------------------------------------------------->>>>>>>>

setwd("C:/Users/rexki/.r_data")

#安装R包
library("BiocManager")
BiocManager::install("Rtsne")
BiocManager::install("plotly")


# ---------------------------------------------------------------------------->>>>>>>>
# MDS
# ---------------------------------------------------------------------------->>>>>>>>

rm(list=ls())

load("GTEx.RData")

# tissue的注释文件
accession_df <- read.csv(file="GTEx_accession_table.csv")
annotation_df <- read.csv(file="GTEx_v7_Annotations_SampleAttributesDS.txt",header = T,sep = "\t")

# make matrix  得到GTEx 矩阵，进行log2 +1,行加基因名
GTEx_mat <- GTEx.TPM.gene.mat[,c(-1,-2)]
GTEx_mat_log2 <- log2(GTEx_mat + 1)
rownames(GTEx_mat_log2) <- as.character(GTEx.TPM.gene.mat$Name)

#计算 gene var 
SD_value <- apply(GTEx_mat,1,FUN = function(x){sd(x)})

#筛选变化最大的基因top 500或1000,2000,3000
GTEx_mat_log2_sort <- GTEx_mat_log2[order(SD_value,decreasing = T),]
GTEx_mat_log2_sort_top <- GTEx_mat_log2_sort[1:1000,]

# get index
col_index = limma::strsplit2(colnames(GTEx_mat_log2_sort_top),split = "\\.")[,5]

# 增加tissue信息,转换基因名.为-
accession_vec <- gsub(x = colnames(GTEx_mat_log2_sort_top),pattern = "\\.",replacement = "-") 
tissue_info <- as.character(annotation_df$SMTS[match(accession_vec,annotation_df$SAMPID)])

#增加tissue一列
new_colname <- sprintf("%s_%s",tissue_info, col_index)
GTEx_mat_log2_sort_top_fix <- GTEx_mat_log2_sort_top
colnames(GTEx_mat_log2_sort_top_fix) <- new_colname

# run MDS

dist_mat = dist(t(GTEx_mat_log2_sort_top), method = "minkowski", p=5)
dist_mat = sqrt(2*(1 - abs(cor(GTEx_mat_log2_sort_top,method = "kendall"))))

MDS_res = cmdscale(dist_mat,eig=TRUE, k=2)

MDS_xy_df = data.frame(MDS_1 = MDS_res$points[,1],
                       MDS_2 = MDS_res$points[,2],
                       tissue = factor(tissue_info))

# plot with ggplot2
library(ggplot2)

ggplot(data = MDS_xy_df, aes(MDS_1,MDS_2,color=tissue)) + 
  geom_point()  + 
  geom_text(aes(label=tissue)) + 
  theme_bw()
```

![image-20210107102043226](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107102043226.png)

## 五、**tSNE**

```r
rm(list = ls())

# ---------------------------------------------------------------------------->>>>>>>>
# 0.make data
# ---------------------------------------------------------------------------->>>>>>>>

setwd("C:/Users/rexki/.r_data")

#安装R包
library("BiocManager")
BiocManager::install("Rtsne")
BiocManager::install("plotly")

# ---------------------------------------------------------------------------->>>>>>>>
# t-SNE
# ---------------------------------------------------------------------------->>>>>>>>
# load package 
rm(list=ls())

library(Rtsne)

load("GTEx.RData")

# tissue的注释文件
accession_df <- read.csv(file="GTEx_accession_table.csv")
annotation_df <- read.csv(file="GTEx_v7_Annotations_SampleAttributesDS.txt",header = T,sep = "\t")

# make matrix  得到GTEx 矩阵，进行log2 +1,行加基因名
GTEx_mat <- GTEx.TPM.gene.mat[,c(-1,-2)]
GTEx_mat_log2 <- log2(GTEx_mat + 1)
rownames(GTEx_mat_log2) <- as.character(GTEx.TPM.gene.mat$Name)

#计算 gene var 
SD_value <- apply(GTEx_mat,1,FUN = function(x){sd(x)})

#筛选变化最大的基因top 500或1000,2000,3000
GTEx_mat_log2_sort <- GTEx_mat_log2[order(SD_value,decreasing = T),]
GTEx_mat_log2_sort_top <- GTEx_mat_log2_sort[1:1000,]

# 增加tissue信息,转换基因名.为-
accession_vec <- gsub(x = colnames(GTEx_mat_log2_sort_top),pattern = "\\.",replacement = "-") 
tissue_info <- as.character(annotation_df$SMTS[match(accession_vec,annotation_df$SAMPID)])

#run t-SNE
set.seed(2021)

tSNE_res = Rtsne(t(GTEx_mat_log2_sort_top),dims = 3,perplexity = 10,pca = T)


#################### plot region #######################
#三维坐标
tSNE_res_df = as.data.frame(tSNE_res$Y)

colnames(tSNE_res_df) = c("tSNE1","tSNE2","tSNE3")

tSNE_res_df$tissue = tissue_info

ggplot(data = tSNE_res_df, aes(tSNE1,tSNE2,color=tissue)) + 
  geom_point()  + 
  geom_text(aes(label=tissue)) + 
  theme_bw()

ggplot(data = tSNE_res_df, aes(tSNE1,tSNE3,color=tissue)) + 
  geom_point()  + 
  geom_text(aes(label=tissue)) + 
  theme_bw()

ggplot(data = tSNE_res_df, aes(tSNE2,tSNE3,color=tissue)) + 
  geom_point()  + 
  geom_text(aes(label=tissue)) + 
  theme_bw()

# use plotly
library(plotly)

p <- plot_ly(tSNE_res_df, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, color = ~tissue) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'tSNE1'),
                      yaxis = list(title = 'tSNE2'),
                      zaxis = list(title = 'tSNE3')))

```

![](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107105357879.png)

## 六、矩阵的奇异值分解 (SVD)

```r
rm(list = ls())

# ---------------------------------------------------------------------------->>>>>>>>
# 0.make data
# ---------------------------------------------------------------------------->>>>>>>>
setwd("C:/Users/rexki/.r_data")

#安装R包
library("BiocManager")
BiocManager::install("")

# ---------------------------------------------------------------------------->>>>>>>>
# SVD
# ---------------------------------------------------------------------------->>>>>>>>
rm(list=ls())

load("GTEx.RData")

# tissue的注释文件
accession_df <- read.csv(file="GTEx_accession_table.csv")
annotation_df <- read.csv(file="GTEx_v7_Annotations_SampleAttributesDS.txt",header = T,sep = "\t")

# make matrix  得到GTEx 矩阵，进行log2 +1,行加基因名
GTEx_mat <- GTEx.TPM.gene.mat[,c(-1,-2)]
GTEx_mat_log2 <- log2(GTEx_mat + 1)
rownames(GTEx_mat_log2) <- as.character(GTEx.TPM.gene.mat$Name)

#计算 gene var 
SD_value <- apply(GTEx_mat,1,FUN = function(x){sd(x)})

#筛选变化最大的基因top 500或1000,2000,3000
GTEx_mat_log2_sort <- GTEx_mat_log2[order(SD_value,decreasing = T),]
GTEx_mat_log2_sort_top <- GTEx_mat_log2_sort[1:1000,]

# get index
col_index = limma::strsplit2(colnames(GTEx_mat_log2_sort_top),split = "\\.")[,5]

# 增加tissue信息,转换基因名.为-
accession_vec <- gsub(x = colnames(GTEx_mat_log2_sort_top),pattern = "\\.",replacement = "-") 
tissue_info <- as.character(annotation_df$SMTS[match(accession_vec,annotation_df$SAMPID)])

#增加tissue一列
new_colname <- sprintf("%s_%s",tissue_info, col_index)
GTEx_mat_log2_sort_top_fix <- GTEx_mat_log2_sort_top
colnames(GTEx_mat_log2_sort_top_fix) <- new_colname


# run SVD

GTEx_mat_log2_fit <- as.matrix(GTEx_mat_log2)

svd_res = svd(GTEx_mat_log2)

U = as.matrix(svd_res$u)
D = diag(svd_res$d)
V = as.matrix(svd_res$v)

back_mat <- U %*% D %*% t(V)

plot(x=GTEx_mat_log2_fit[1:3000], y=back_mat[1:3000])

#压缩

i <- 20
back_mat_part <- svd_res$u[,1:i] %*% diag(svd_res$d[1:i]) %*% t(svd_res$v[,1:i])

plot(x=GTEx_mat_log2_fit[1:3000], y=back_mat_part[1:3000])

```

![](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210107113706129.png)

## 七、**多元回归**

```r
# ---------------------------------------------------------------------------->>>>>>>>
# linear regression
# ---------------------------------------------------------------------------->>>>>>>>
rm(list = ls())

#加载数据
diabetes_df <- read.csv(file = "diabetes.csv")

# ---------------------------------------------------------------------------->>>>>>>>
# 1. run lm 基础回归
# ---------------------------------------------------------------------------->>>>>>>>
# lm 
lm.res = lm(y ~ ., diabetes_df)
summary(lm.res)

# kappa
kappa(diabetes_df[,-1])

# vif  一般vif>5 ,kappa > 15,代表存在共线性
library(car)
sort(car::vif(lm.res),decreasing = T)[1:10]

# ---------------------------------------------------------------------------->>>>>>>>
# 2. run step 逐步回归
# ---------------------------------------------------------------------------->>>>>>>>
lm.res.step = step(lm(y ~ ., diabetes_df))
summary(lm.res.step)

# vif 
sort(vif(lm.res.step),decreasing = T)

# ---------------------------------------------------------------------------->>>>>>>>
# 3. ridge 岭回归
# ---------------------------------------------------------------------------->>>>>>>>
library(ridge)

lm.res.ridge = ridge::linearRidge(y ~ . , diabetes_df)
summary(lm.res.ridge)

# ---------------------------------------------------------------------------->>>>>>>>
# 4. lasso 回归
# ---------------------------------------------------------------------------->>>>>>>>
library(lars)

lasso.x = as.matrix(diabetes_df[,-1])
lasso.y = as.matrix(diabetes_df[,1])

lasso.res = lars::lars(x=lasso.x, y=lasso.y)
summary(lasso.res)

# get coef
min(lasso.res$Cp)
lasso.coef = lars::coef.lars(lasso.res, mode="step",s=15)

# ---------------------------------------------------------------------------->>>>>>>>
# 5. partial least squares 偏最小二乘回归
# ---------------------------------------------------------------------------->>>>>>>>
library(pls)

pls.x = as.matrix(diabetes_df[,-1])
pls.y = as.matrix(diabetes_df[,1])

pls.res = pls::plsr(pls.y ~ pls.x, 64, validation="CV")

pls.res$loadings
pls.res$coefficients[,,6]

plot(R2(pls.res))
```

