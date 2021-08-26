# 生信R语言实战

[TOC]

## 1.根据R包`org.Hs.eg.db`找到下面ensembl 基因ID 对应的基因名(symbol)

```r
ENSG00000000003.13
ENSG00000000005.5
ENSG00000000419.11
ENSG00000000457.12
ENSG00000000460.15
ENSG00000000938.11
```

```r
rm(list = ls())
options(stringsAsFactors = F)
library(org.Hs.eg.db)
library(limma)

#加载数据
a1 = read.table('e1.txt')

g2s = toTable(org.Hs.egSYMBOL)
g2e = toTable(org.Hs.egENSEMBL)

#新增一列ensembl_id，为去掉点后的对应基因ID
a1$ensembl_id = unlist(lapply(a1$V1,function(x){
  limma::strsplit2(x,split = '\\.')[,1]
}))

#merge
tmp = merge(a1,g2e,by='ensembl_id')
tmp = merge(tmp,g2s,by='gene_id')
```

## 2.根据R包`hgu133a.db`找到下面探针对应的基因名(symbol)

```r
1053_at
117_at
121_at
1255_g_at
1316_at
1320_at
1405_i_at
1431_at
1438_at
1487_at
1494_f_at
1598_g_at
160020_at
1729_at
177_at
```

```r
rm(list = ls())
options(stringsAsFactors = F)
library(hgu133a.db)

#加载数据
a2 = read.table('e2.txt')
colnames(a2)='probe_id'
ids=toTable(hgu133aSYMBOL)

#merge
tmp1=merge(ids,a2,by='probe_id')
tmp2=ids[match(a2$probe_id,ids$probe_id),]
```

## 3.找到R包`CLL`内置的数据集的表达矩阵里面的TP53基因的表达量，并且绘制在 `progres.-stable`分组的boxplot图

```r
rm(list = ls())
options(stringsAsFactors = F)

suppressPackageStartupMessages(library(CLL))
data(sCLLex)
sCLLex
exprSet=exprs(sCLLex) 
pd=pData(sCLLex)
library(hgu95av2.db)

ids=toTable(hgu95av2SYMBOL)

boxplot(exprSet['1939_at',] ~ pd$Disease) ## sig
boxplot(exprSet['1974_s_at',] ~ pd$Disease)
boxplot(exprSet['31618_at',] ~ pd$Disease)
```

## 4.找到BRCA1基因在TCGA数据库的乳腺癌数据集([Breast Invasive Carcinoma (TCGA, PanCancer Atlas)](http://www.cbioportal.org/study?id=brca_tcga_pan_can_atlas_2018))的表达情况

```
http://www.cbioportal.org/index.do
```

```r
rm(list = ls())
options(stringsAsFactors = F)
a4=read.table('e4-plot.txt',sep = '\t',fill = T,header = T)

colnames(a4)=c('id','subtype','expression','mut')
dat=a4
library(ggstatsplot)
ggbetweenstats(data =dat, x = subtype,  y = expression)
library(ggplot2)
ggsave('plot-again-BRCA1-TCGA-BRCA-cbioportal.png')
```

## 5.找到TP53基因在TCGA数据库的乳腺癌数据集的表达量分组看其是否影响生存

```r
http://www.oncolnc.org/
```

```r
rm(list = ls())
options(stringsAsFactors = F)
a5=read.table('BRCA_7157_50_50.csv',sep = ',',fill = T,header = T)

dat=a5
library(ggplot2)
library(survival)
library(survminer) 
table(dat$Status)
dat$Status=ifelse(dat$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=dat)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsave('survival_TP53_in_BRCA_TCGA.png')


###################################
head(a5)
b=read.table('e4-plot.txt',sep = '\t',fill = T,header = T)
colnames(b)=c('Patient','subtype','expression','mut')
head(b)
b$Patient=substring(b$Patient,1,12)
tmp=merge(a5,b,by='Patient')

table(tmp$subtype)

type='BRCA_LumB'
x=tmp[tmp$subtype==type,] 
library(ggplot2)
library(survival)
library(survminer) 
#table(x$Status)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)  

table(tmp$subtype)

type='BRCA_Normal'
x=tmp[tmp$subtype==type,] 
library(ggplot2)
library(survival)
library(survminer) 
#table(x$Status)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE) 

table(tmp$subtype)
type='BRCA_Basal'

x=tmp[tmp$subtype==type,] 
library(ggplot2)
library(survival)
library(survminer) 
#table(x$Status)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE) 

table(tmp$subtype)

type='BRCA_Her2'
x=tmp[tmp$subtype==type,] 
library(ggplot2)
library(survival)
library(survminer) 
#table(x$Status)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE) 

table(tmp$subtype)

type='BRCA_LumA'
x=tmp[tmp$subtype==type,] 
library(ggplot2)
library(survival)
library(survminer) 
#table(x$Status)
x$Status=ifelse(x$Status=='Dead',1,0)
sfit <- survfit(Surv(Days, Status)~Group, data=x)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE) 
```

## 6.下载数据集GSE17215的表达矩阵并且提取下面的基因画热图

```r
ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDCA1 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KNTC2 KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 ORC6L PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T
```

根据基因名拿到探针ID，缩小表达矩阵绘制热图，没有检查到的基因直接忽略即可

```r
rm(list = ls()) 
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据 
f='GSE17215_eSet.Rdata'

library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE17215', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE17215_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)


library(hgu133a.db)
ids=toTable(hgu133aSYMBOL)
head(ids)
dat=dat[ids$probe_id,]
dat[1:4,1:4] 
ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol
dat[1:4,1:4]  
dim(dat)

ng='ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDCA1 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KNTC2 KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 ORC6L PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T'
ng=strsplit(ng,' ')[[1]]
table(ng %in%  rownames(dat))
ng=ng[ng %in%  rownames(dat)]
dat=dat[ng,]
dat=log2(dat)
pheatmap::pheatmap(dat,scale = 'row')
```

## 7.下载数据集GSE24673的表达矩阵计算样本的相关性并且绘制热图，需要标记上样本分组信息

```r
rm(list = ls())  
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据 
f='GSE24673_eSet.Rdata'

library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE24673', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE24673_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)
pd=pData(a)
group_list=c('rbc','rbc','rbc',
             'rbn','rbn','rbn',
             'rbc','rbc','rbc',
             'normal','normal')
dat[1:4,1:4]
M=cor(dat)
pheatmap::pheatmap(M)
tmp=data.frame(g=group_list)
rownames(tmp)=colnames(M)
pheatmap::pheatmap(M,annotation_col = tmp)
```

![image-20210108171444736](C:%5CUsers%5Crexki%5CAppData%5CRoaming%5CTypora%5Ctypora-user-images%5Cimage-20210108171444736.png)

## 8.找到 GPL6244 platform of Affymetrix Human Gene 1.0 ST Array 对应的R的bioconductor注释包，并且安装它

```r
options()$repos
options()$BioC_mirror 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
BiocManager::install("hugene10sttranscriptcluster.db",ask = F,update = F)
options()$repos
options()$BioC_mirror
```

## 9.下载数据集GSE42872的表达矩阵，并且分别挑选出 所有样本的(平均表达量/sd/mad/)最大的探针，并且找到它们对应的基因。

```r
rm(list = ls())  
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据 
f='GSE42872_eSet.Rdata'

library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE42872_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)
pd=pData(a)
# (平均表达量/sd/mad/)最大的探针
boxplot(dat)
sort(apply(dat,1,mean),decreasing = T)[1]
sort(apply(dat,1,sd),decreasing = T)[1]
sort(apply(dat,1,mad),decreasing = T)[1]
```

## 10.下载数据集GSE42872的表达矩阵，并且根据分组使用limma做差异分析，得到差异结果矩阵

```r
rm(list = ls())  
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据 
f='GSE42872_eSet.Rdata'

library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE42872_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)
pd=pData(a)
# (平均表达量/sd/mad/)最大的探针
boxplot(dat)
group_list=unlist(lapply(pd$title,function(x){
  strsplit(x,' ')[[1]][4]
}))


exprSet=dat
exprSet[1:4,1:4]
# DEG by limma 
suppressMessages(library(limma)) 
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix<-makeContrasts("progres.-stable",levels = design)
contrast.matrix 
##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
```

# 统计学-表达矩阵(R语言）

```r
rm(list = ls())
options(stringsAsFactors = F)

#加载数据
library(airway)
data(airway)
RNAseq_expr=assay(airway)
RNAseq_gl=colData(airway)[,3] #分组矩阵

```

## 1.把RNAseq_expr第一列全部加1后取log2后计算平均值和标准差

```r
####################################################################
# 1.把RNAseq_expr第一列全部加1后取log2后计算平均值和标准差
####################################################################
tmp=log2(RNAseq_expr[,1]+1)
mean(tmp)
sd(tmp)
```

2.根据上一步得到平均值和标准差生成同样个数的随机的正态分布数值

```r
####################################################################
# 2.根据上一步得到平均值和标准差生成同样个数的随机的正态分布数值
####################################################################
a=rnorm(length(tmp),mean = mean(tmp),sd = sd(tmp))
a=sort(a)
plot(a)
points(sort(tmp))
```

## 3.删除RNAseq_expr第一列低于5的数据后，重复Q1和Q2

```r
####################################################################
# 3.删除RNAseq_expr第一列低于5的数据后，重复Q1和Q2
####################################################################
tmp=RNAseq_expr[,1]
tmp=tmp[tmp>5]
tmp=log2(tmp)
a=rnorm(length(tmp),mean = mean(tmp),sd = sd(tmp))
a=sort(a)
plot(a)
points(sort(tmp))
```

## 4.基于Q3对RNAseq_expr的第一列和第二列进行T检验

```r
####################################################################
# 4.基于Q3对RNAseq_expr的第一列和第二列进行T检验
####################################################################
x=RNAseq_expr[,1]
x=x[x>5]
x=log2(x)

y=RNAseq_expr[,2]
y=y[y>5]
y=log2(y)

t.test(x,y)
library(ggpubr)
df=data.frame(value=c(x,y),
              group=c(rep('x',length(x)),rep('y',length(y))))
ggboxplot(df, y = "value", x = "group")
```

## 5.取RNAseq_expr行之和最大的那一行根据分组矩阵进行T检验

```r
####################################################################
# 5.取RNAseq_expr行之和最大的那一行根据分组矩阵进行T检验
####################################################################
pos=which.max(rowSums(RNAseq_expr))
t.test(RNAseq_expr[pos,]~RNAseq_gl)
pos
```

## 6.取RNAseq_expr的MAD最大的那一行根据分组矩阵进行T检验

```r
####################################################################
# 6.取RNAseq_expr的MAD最大的那一行根据分组矩阵进行T检验
####################################################################
pos=which.max(apply(RNAseq_expr,1,mad))
t.test(RNAseq_expr[pos,]~RNAseq_gl)
pos
```

## 7.对RNAseq_expr全部加1后取log2后重复Q5和Q6

```r
####################################################################
# 7.对RNAseq_expr全部加1后取log2后重复Q5和Q6
####################################################################
RNAseq_expr=log2(RNAseq_expr+1)
pos=which.max(rowSums(RNAseq_expr))
pos
t.test(RNAseq_expr[pos,]~RNAseq_gl)
pos=which.max(apply(RNAseq_expr,1,mad))
pos
t.test(RNAseq_expr[pos,]~RNAseq_gl)
```

## 8.取RNAseq_expr矩阵的MAD最高的100行，对列和行分别进行层次聚类

```r
####################################################################
# 8.取RNAseq_expr矩阵的MAD最高的100行，对列和行分别进行层次聚类
####################################################################
cg=names(tail(sort(apply(RNAseq_expr,1,mad)),100))
dat=RNAseq_expr[cg,]
plot(hclust(dist(t(dat))))
colnames(dat)
RNAseq_gl
plot(hclust(dist( dat )))
```

## 9.对Q8矩阵按照行和列分别归一化并且热图可视化

```r
####################################################################
#  9.对Q8矩阵按照行和列分别归一化并且热图可视
####################################################################
cg=names(tail(sort(apply(RNAseq_expr,1,mad)),100))
dat=RNAseq_expr[cg,]
pheatmap::pheatmap(scale(dat))
pheatmap::pheatmap(t(scale(t(dat))))
```

![image-20210109151601140](D:%5CDesktop%5C%E7%94%9F%E4%BF%A1%E6%8A%80%E8%83%BD%E6%A0%91%5C4%E3%80%81%E7%94%9F%E4%BF%A1%E4%BA%BA%E5%BA%94%E8%AF%A5%E8%BF%99%E6%A0%B7%E5%AD%A6R%E8%AF%AD%E8%A8%80%5C%E7%94%9F%E4%BF%A1R%E8%AF%AD%E8%A8%80%E5%AE%9E%E6%88%98.assets%5Cimage-20210109151601140.png)

![image-20210109151617562](D:%5CDesktop%5C%E7%94%9F%E4%BF%A1%E6%8A%80%E8%83%BD%E6%A0%91%5C4%E3%80%81%E7%94%9F%E4%BF%A1%E4%BA%BA%E5%BA%94%E8%AF%A5%E8%BF%99%E6%A0%B7%E5%AD%A6R%E8%AF%AD%E8%A8%80%5C%E7%94%9F%E4%BF%A1R%E8%AF%AD%E8%A8%80%E5%AE%9E%E6%88%98.assets%5Cimage-20210109151617562.png)

# 统计学-统计检验（R语言）

```r
options(stringsAsFactors = F)
rm(list=ls())
library(airway) 
RNAseq_expr=assay(airway)
e1=RNAseq_expr[apply(RNAseq_expr,1,function(x) sum(x>0)>1),] #过滤每一列都为0的行
RNAseq_gl=colData(airway)[,3] #分组矩阵
```

## 1.对e1每一行独立根据分组矩阵进行T检验，检查为什么有些行T检验失败

```r
####################################################################
#  1.对e1每一行独立根据分组矩阵进行T检验，检查为什么有些行T检验失败
####################################################################
apply(e1, 1, function(x){
  t.test(x~RNAseq_gl)$p.value
})
```

## 2.找出T检验失败的行并且从e1矩阵剔除掉

```r
####################################################################
# 2.找出T检验失败的行并且从e1矩阵剔除掉
####################################################################
e1_a=e1[,RNAseq_gl=='trt']
e1_b=e1[,RNAseq_gl=='untrt']
a_filter=apply(e1_a, 1,function(x) sd(x)>0)
b_filter=apply(e1_b, 1,function(x) sd(x)>0)
table(a_filter,b_filter)
e1=e1[a_filter | b_filter,]
```

## 3.对e1矩阵进行加1后log2的归一化命名为e2,对e1,e2的T检验P值做相关性分析

```r
####################################################################
# 3.对e1矩阵进行加1后log2的归一化命名为e2,对e1,e2的T检验P值做相关性分析
####################################################################
p1=apply(e1, 1, function(x){
  t.test(x~RNAseq_gl)$p.value
}) 
e2=log(e1+1)
p2=apply(e2, 1, function(x){
  t.test(x~RNAseq_gl)$p.value
}) 
plot(p1,p2)
cor(p1,p2)

```

![image-20210109153237414](D:%5CDesktop%5C%E7%94%9F%E4%BF%A1%E6%8A%80%E8%83%BD%E6%A0%91%5C4%E3%80%81%E7%94%9F%E4%BF%A1%E4%BA%BA%E5%BA%94%E8%AF%A5%E8%BF%99%E6%A0%B7%E5%AD%A6R%E8%AF%AD%E8%A8%80%5C%E7%94%9F%E4%BF%A1R%E8%AF%AD%E8%A8%80%E5%AE%9E%E6%88%98.assets%5Cimage-20210109153237414.png)