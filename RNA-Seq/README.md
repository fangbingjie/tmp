# GSE56010-m6A RNA-Seq转录分析

> N6-methyladenosine-dependent RNA structural switches regulate RNA–protein interactions	
>
> www.nature.com/doifinder/10.1038/nature14234


## 一、环境准备

```shell
#安装miniconda
SRC_PATH = ~/src  # 软件下载目录
cd $SRC_PATH  # 进入目录
# 安装aspera，之后下载sra或fastq文件使用
wget -c http://download.asperasoft.com/download/sw/connect/3.7.4/aspera-connect-3.7.4.147727-linux-64.tar.gz

#安装miniconda，之后使用各种软件可以集成
wget -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-py37_4.9.2-Linux-x86_64.sh
bash Miniconda3-py37_4.9.2-Linux-x86_64.sh
source ~/.bashrc  # 加载环境变量
conda --help  # 有提示证明安装成功

#anaconda添加源
conda config --add channels http://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels http://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels http://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels http://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/  
conda info  # 查看配置信息

#创建一个conda虚拟环境并激活
conda create -n RNA-Seq python=2
conda activate -n RNA-Seq

#安装对应软件，多安装了几个
conda install -y mamba   # 安装mamba，之后用mamba下载更快
mamba install -y sra-tools  # sra下载（prefetch），sra转换fastq（fastq-dump）
mamba install -y trimmomatic cutadapt trim-galore # 去接头工具
mamba install -y fastqc multiqc # qc质控、多qc查看
mamba install -y star hisat2 bowtie2  # 比对工具
mamba install -y salmon bwa # 比对工具（转录比对、基因组比对）
mamba install -y samtools # 压缩排序
mamba install -y subread htseq  # 计数RNA-seq（featureCounts）
mamba install -y bedtools deeptools  # 功能注释
```



## 二、下载sra文件或fastq文件

```shell
#创建工作目录,设为WORK_DIR变量
mkdir -p ~/data/RNA-Seq/GSE56010-m6A/1.raw_fq 
WORK_DIR=~/data/RNA-Seq/GSE56010-m6A 
cd $WORK_DIR/1.raw_fq

# ENA数据库获取fq.txt,里面文件如图所示，aspera批量下载
# 以下代码保存为step1-aspera-fq.sh
cat fq.txt |while read id 
do
	~/.aspera/connect/bin/ascp -QT -l 300m -P33001  \
	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh   \
	--host=fasp.sra.ebi.ac.uk --user=era-fasp --mode=recv \
	$id  ./
done
nohup bash step1-aspera-fq.sh 1> step1-aspera-fq.log 2>&1 &
```

​	**下载的若是sra文件，需要转换成fastq文件**

```shell
#对文件夹下的所有sra文件批量处理,以下代码保存为sra2fq.sh
ls *sra | xargs -n1 -i fastq-dump --split-3 {} --gzip

nohup bash sra2fq.sh 1> sra2fq.log 2>&1 &
```

![image-20210113234527271](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210113234527271.png)

![image-20210113234414700](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210113234414700.png)

```shell
    SRR			对应样本处理
SRR1573494	RNA-seq-control-rep1
SRR1573495	RNA-seq-control-rep2
SRR1573498	RNA-seq-METTL3KD-rep1
SRR1573499	RNA-seq-METTL3KD-rep2
SRR1573500	RNA-seq-METTL14KD-rep1
SRR1573501	RNA-seq-METTL14KD-rep2
```



## 三、对fastq数据进行质控

```shell
# 创建fastqc输出目录
mkdir -p $WORK_DIR/2.fastqc_output

# 使用fastqc进行质控,此处使用6线程,以下代码保存为step2-fastqc.sh
ls *gz | xargs fastqc -t 6 -o $WORK_DIR/2.fastqc_output

nohup  bash step2-fastqc.sh 1> step2-fastqc.log 2>&1 &
```

```shell
# 通过mutliqc 查看
cd $WORK_DIR/2.fastqc_output
multiqc ./   # 需要python3环境

#multiqc_report.html是所有样本的整合报告
```

![image-20210112115119348](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210112115119348.png)

![image-20210112115139752](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210112115139752.png)

## 四、接头处理后再次质控

```shell
# 创建一个去接头后的数据存放目录
mkdir -p $WORK_DIR/3.clean_fq/unpaired
cd $WORK_DIR/3.clean_fq

ls $WORK_DIR/1.raw_fq/*_1.fastq.gz  >1
ls $WORK_DIR/1.raw_fq/*_2.fastq.gz  >2
paste 1 2  > config

# 保存以下代码为step3-cutadapt.sh
cat $1 |while read id
do
        arr=(${id})
        fq1=${arr[0]}
        fq2=${arr[1]} 
  trimmomatic PE \
  -threads 4 -phred33  $fq1 $fq2 \
  clean.$(basename $fq1) unpaired.$(basename $fq1) \
  clean.$(basename $fq2) unpaired.$(basename $fq2) \
  ILLUMINACLIP:/home/fbj/miniconda3/envs/RNA-Seq/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 \
  SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50
done

nohup  bash step3-cutadapt.sh config 1> step3-cutadapt.log 2>&1 &  # config是传递进去的参数
ls unpaired* | xargs -I {} mv {} $WORK_DIR/3.clean_fq/unpaired  # 将unpaired文件转移至相应目录
```

![image-20210113235144391](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210113235144391.png)

```shell
# 创建clean_fastqc输出目录
mkdir -p $WORK_DIR/4.clean_fastqc_output
cd $WORK_DIR/4.clean_fastqc_output

# 使用fastqc进行质控,此处使用6线程,以下代码保存为step4-fastqc.sh
ls $WORK_DIR/3.clean_fq/*gz | xargs fastqc -t 6 -o $WORK_DIR/4.clean_fastqc_output

nohup  bash step4-fastqc.sh 1> step4-fastqc.log 2>&1 &
```

```shell
# 通过mutliqc 查看
cd $WORK_DIR/4.clean_fastqc_output
multiqc ./   # 需要python3环境

#multiqc_report.html是所有样本的整合报告
```

![image-20210112203107429](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210112203107429.png)

![image-20210112203253753](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210112203253753.png)

## 五、下载比对参考基因组文件，为HISAT2配置index

HISAT2是TopHat2/Bowti2的继任者，使用改进的BWT算法，实现了更快的速度和更少的资源占用

HISAT2官网有已经配置好的人类和小鼠基因组索引，可以直接下载，下载最好用迅雷，4个多G

```shell
# 创建用于放置hiseq2索引的目录,设为REF_DIR变量
mkdir -p ~/data/reference/GRCh38 
REF_DIR=~/data/reference/GRCh38
cd $REF_DIR

# 从hiseq2官网下载人类参考基因组索引
wget -c https://genome-idx.s3.amazonaws.com/hisat/grch38_snptran.tar.gz  
```

```shell
# 从ensemble中下载最新版本的人类基因组注释文件（gtf格式）
wget -c ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.chr.gtf.gz 

# 从ensemble中下载人类基因组序列
wget -c ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

```

#### 需要自己配置hiseq2索引时使用以下方式：

```shell
# 配置HISAT2的index-只有基因组（可以下载，不建议建索引）
hisat2-build –p 6 Homo_sapiens.GRCh38.dna.toplevel.fa GRCh38_ensembl_dna > hisat2_build.log 2>&1 &  

# 建立基因组+转录组+SNP索引
# 需要下载snp信息文件
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz

# HISAT2提供两个Python脚本将GTF文件转换成hisat2-build能使用的文件
extract_exons.py Homo_sapiens.GRCh38.102.chr.gtf > genome.exon &  # 外显子信息文件
extract_splice_sites.py Homo_sapiens.GRCh38.102.chr.gtf > genome.ss &  # 可变剪切位点信息文件
extract_snps.py snp142Common.txt > genome.snp &  # snp信息文件

# 开始构建索引
hisat2-build -p 6 Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--snp genome.snp \
	--ss genome.ss \
	--exon genome.exon \
	GRCh38_ensembl_dna_tran > hisat2_build_tran.log 2>&1 & 
```



## 六、序列比对

```shell
# 创建用于存放比对结果的目录
mkdir -p ~/data/RNA-Seq/GSE56010-m6A/5.align_output/bam
ALIGN_DIR=~/data/RNA-Seq/GSE56010-m6A/5.align_output
cd $ALIGN_DIR/bam


# 批量生成Sam文件再压缩为Bam文件,以下代码保存为step5-fq2sam.sh (python3)
ls $WORK_DIR/3.clean_fq/*gz | xargs -I {} basename {} |cut -d "_" -f 1 | sort -u | while read id
do
	hisat2 -p 6 -x $REF_DIR/grch38_snp_tran/genome_snp_tran \
	-1 $WORK_DIR/3.clean_fq/${id}_1.fastq.gz \
	-2 $WORK_DIR/3.clean_fq/${id}_2.fastq.gz \
	-S $ALIGN_DIR/bam/${id}_hisat2.sam 
done

nohup bash step5-fq2sam.sh > step5-fq2sam.log 2>&1 &
# -x 指定基因组索引前缀
# -1 指定第一个fastq文件
# -2 指定第二个fastq文件
# -S 指定输出的SAM文件

```



## 七、Sam转Bam文件并排序

```shell
#以下代码保存为step6-sam2bam.sh
ls *.sam | cut -d "_" -f 1 |while read id
do 
	samtools sort -O bam -@ 6 -o ${id}_sorted.bam ${id}_hisat2.sam
done

nohup bash step6-sam2bam.sh  > step6-sam2bam.log 2>&1 &

rm -rf *.sam
```

![image-20210114153336671](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210114153336671.png)

## 八、Bam构建索引，multiqc检查flagstat

```shell
# 批量创建索引
ls *.bam |xargs -i samtools index -@ 6 {} &

# 批量创建falgstat
ls *.bam | cut -d "_" -f 1 |while read id ;do samtools flagstat -@ 6 ${id}_sorted.bam > ${id}.flagstat ;done &

multiqc ./*.flagstat  #  python3
```

![image-20210114172434754](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210114172434754.png)

![image-20210114172841563](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210114172841563.png)

## 八、featureCounts得到count值

```shell
# 创建存放count矩阵的目录
mkdir -p $ALIGN_DIR/count
cd $ALIGN_DIR/count

# 以下代码保存为step7-countdata.sh
gtf="$REF_DIR/Homo_sapiens.GRCh38.102.chr.gtf"
featureCounts  -t exon  -g gene_id  -T 6  -p  -Q 20  -a $gtf  -o all.id.txt $ALIGN_DIR/bam/*.bam  1>counts.id.log 2>&1 

nohup bash step7-countdata.sh > step7-countdata.log 2>&1 &
```

```shell
# featureCounts参数

-a 输入GTF/GFF基因组注释文件
-p 这个参数是针对paired-end数据
-F 指定-a注释文件的格式，默认是GTF
-g 从注释文件中提取Meta-features信息用于read count，默认是gene_id
-t 跟-g一样的意思，其是默认将exon作为一个feature
-o 输出文件
-T 多线程数
-Q 最小mapping质量值
```



## 九、使用R及count矩阵进行差异分析（DEG）

### 9-1 构造一个run_DEG_RNAseq函数

```R
#################################################################################################
# Author Bingjie Fang
# E-mail fbj1001@126.com
# run_DEG_RNA-seq
#################################################################################################
rm(list = ls())
options(stringsAsFactors = F)
library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(pheatmap)

load(file ='exprSet.Rdata')

draw_h_v <- function(exprSet,need_DEG,n='DEseq2'){
  ## 差异表达基因分析需要两列数据 一个是log2FoldChange，一个是padj
    
  ## heatmap 热图
  choose_gene=head(rownames(need_DEG),50) ## 50 maybe better
  choose_matrix=exprSet[choose_gene,]  ## 提取前50的表达矩阵
  choose_matrix=t(scale(t(choose_matrix)))  ##归一化
  pheatmap(choose_matrix,filename = paste0(n,'_need_DEG_top50_heatmap.png'))
  
  
  ## volcano火山图
  logFC_cutoff <- with(need_DEG,mean(abs( log2FoldChange)) + 2*sd(abs( log2FoldChange)) )
  # logFC_cutoff=1
  
  need_DEG$change = as.factor(ifelse(need_DEG$padj < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,ifelse(need_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(need_DEG[need_DEG$change =='DOWN',])
  )
  library(ggplot2)
  g = ggplot(data=need_DEG, 
             aes(x=log2FoldChange, y=-log10(padj), 
                 color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
  print(g)
  ggsave(g,filename = paste0(n,'_volcano.png'))
}


run_DEG_RNAseq <- function(exprSet,group_list,
                           g1="ctrl",g2="KD",
                           pro='test'){
  print(table(group_list))
  colnames(exprSet)
  cat(paste0('Now process the project : ',pro))
  
### ------------------------------------------------------------------>>>>>>>>>>
###
### Firstly run DEseq2 
###
### ------------------------------------------------------------------>>>>>>>>>>

# 过滤,每行count数总和大于20且都不为0
exprSet <- exprSet[rowSums(exprSet) > 20 & apply(exprSet,1,function(x){ all(x > 0) }),]

# 创建DESeq2对象dds
group_list$condition <- as.factor(group_list$condition)  # condition转换为因子
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = group_list,
                              design = ~condition)

dds <- DESeq(dds)  # 直接用DESeq函数处理
dds.res <- results(dds)   # 提取想要的差异分析结果，g2，g1相比，g1是对照组
dds.res.sorted <- dds.res[order(dds.res$padj),]  # 按照padj排序
DEG =as.data.frame(dds.res.sorted)


## MAplot画图
png(filename = 'DESeq2_MAplot.png')
DESeq2::plotMA(dds.res,alpha=0.001)
dev.off()

nrDEG=DEG
DEG_DEseq2=nrDEG
nrDEG=DEG_DEseq2[,c(2,6)]  #提取lodFC，padj两列
colnames(nrDEG)=c('log2FoldChange','padj') 
draw_h_v(exprSet,nrDEG,paste0(pro,'_DEseq2'))  #根据函数画图

### ------------------------------------------------------------------>>>>>>>>>>
###
### Then run edgeR 
###
### ------------------------------------------------------------------>>>>>>>>>>

dge.edgeR.obj <- DGEList(counts = exprSet, group = group_list$condition)

dge.edgeR.obj <- calcNormFactors(dge.edgeR.obj,method = "TMM")  #选择Normalization方式

dge.edgeR.obj$samples

# design matrix
design.mat <- model.matrix(~group_list$condition)
rownames(design.mat) <- colnames(dge.edgeR.obj)
colnames(design.mat) <- levels(group_list$condition)

dge.edgeR.obj <- estimateDisp(dge.edgeR.obj,design.mat)  # 函数直接运行

# plot dispersion
plotBCV(dge.edgeR.obj, cex = 0.8)

# plot var and mean
plotMeanVar(dge.edgeR.obj, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)

# -------------------------------------------------------->>>>>>>>>>
# test with likelihood ratio test
# -------------------------------------------------------->>>>>>>>>>
dge.edgeR.obj.res <- exactTest(dge.edgeR.obj)
DEGs.edgeR.res <- as.data.frame(topTags(dge.edgeR.obj.res,n=nrow(exprSet),sort.by = "logFC"))

# MA plot
select.sign.gene = decideTestsDGE(dge.edgeR.obj.res, p.value = 0.001) 
select.sign.gene_id = rownames(dge.edgeR.obj.res)[as.logical(select.sign.gene)]
plotSmear(dge.edgeR.obj.res, de.tags = select.sign.gene_id, cex = 0.5,ylim=c(-4,4)) 
abline(h = c(-2, 2), col = "blue")


# -------------------------------------------------------->>>>>>>>>>
# test with likelihood ratio test   LR越接近1越不显著
# -------------------------------------------------------->>>>>>>>>>
fit <- glmFit(dge.edgeR.obj, design.mat)
lrt <- glmLRT(fit, coef=2)
nr_DEG <- as.data.frame(topTags(lrt,n=nrow(dge.edgeR.obj),sort.by = "logFC"))

DEG_edgeR=nr_DEG
nrDEG=DEG_edgeR[,c(1,5)]
colnames(nrDEG)=c('log2FoldChange','padj')
draw_h_v(exprSet,nrDEG,paste0(pro,'_edgeR'))


### ------------------------------------------------------------------>>>>>>>>>>
###
### Then run limma 
###
### ------------------------------------------------------------------>>>>>>>>>> 

voom.res <- limma::voom(exprSet, design = design.mat, normalize.method = "quantile")
voom.fit <- lmFit(voom.res, design.mat)
voom.eBayes <- eBayes(voom.fit)
voom.DEG <- topTable(voom.eBayes,coef = 2, number = Inf)


##################################################################################
save(DEG_DEseq2,DEG_edgeR,
     dds,exprSet,group_list,
     file = paste0(pro,'_DEG_results.Rdata')) 

allg=intersect(rownames(DEG_DEseq2),rownames(DEG_edgeR))

nrDEG=cbind(DEG_edgeR[allg,c(1,5)],
            DEG_DEseq2[allg,c(2,6)]) 
colnames(nrDEG)
print(cor(nrDEG[,c(1,3)]))
write.csv(nrDEG,file =  paste0(pro,'_DEG_results.csv'))
return(nrDEG)
}

getDEGs <- function(DEG_DEseq2,DEG_edgeR,thre_logFC=1,thre_p=0.05){
  head(DEG_DEseq2)
  head(DEG_edgeR)
  # head(DEG_limma_voom)
  # thre_logFC=1;thre_p=0.05
  u1=rownames(DEG_DEseq2[with(DEG_DEseq2,log2FoldChange>thre_logFC & padj<thre_p),])
  u2=rownames(DEG_edgeR[with(DEG_edgeR,logFC>thre_logFC & FDR<thre_p),])
  # u3=rownames(DEG_limma_voom[with(DEG_limma_voom,logFC>thre_logFC & adj.P.Val<thre_p),])
  
  d1=rownames(DEG_DEseq2[with(DEG_DEseq2,log2FoldChange < -thre_logFC & padj<thre_p),])
  d2=rownames(DEG_edgeR[with(DEG_edgeR,logFC < -thre_logFC & FDR<thre_p),])
  # d3=rownames(DEG_limma_voom[with(DEG_limma_voom,logFC < -thre_logFC & adj.P.Val<thre_p),])
  
  u=intersect(u1,u2)
  d=intersect(d1,d2)
  
  return(list(up=u,down=d))
}
```

### 9-2 提取表达矩阵，创建分组信息

```R
#################################################################################################
# Author Bingjie Fang
# E-mail fbj1001@126.com
# exprSet
#################################################################################################
rm(list = ls())
options(stringsAsFactors = F)

# 安装相应的包
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
}

# 加载R包
library(org.Hs.eg.db) 
library(limma)

# -------------------------------------------------------------------->>>>>>>>>>
# 1.转换ensemblID和geneID
# -------------------------------------------------------------------->>>>>>>>>>
raw_df <- read.table(file = "all.id.txt",header = T)

# 若ensembl有点号后面的，运行{}里代码去掉版本号
{
  raw_df$ensembl_id = unlist(lapply(raw_df$V1,function(x){
    limma::strsplit2(x,split = '\\.')[,1]
  }))
}

# 提取包中的基因注释信息
g2s <- toTable(org.Hs.egSYMBOL)
g2e <- toTable(org.Hs.egENSEMBL)

# 关联ensembl和symbol，得到过滤后的矩阵raw_fix_df
colnames(raw_df)[1] <- 'ensembl_id'  # 将第一列的列名geneid改为ensembl_id
raw_median_df <- merge(raw_df,g2e,by='ensembl_id',all.x = T)  # 通过ensembl_id关联raw_df与g2e
raw_fix_df <- merge(raw_median_df,g2s,by='gene_id',all.x = T)  # 通过gene_id关联raw_median_df与g2s

# 过滤矩阵raw_fix_df中重复行，NA缺省值的行，转换行名
raw_fix_df <- raw_fix_df[!duplicated(raw_fix_df$ensembl_id),]  # 删除ensembl_id重复的条目
raw_fix_df <- raw_fix_df[match(raw_df$ensembl_id,raw_fix_df$ensembl_id),]  # 通过match()函数将raw_df中ensembl_id的顺序放到raw_fix_df的ensembl_id中去
raw_fix_df <- na.omit(raw_fix_df)  # 删除NA所在的行
raw_fix_df <- raw_fix_df[!duplicated(raw_fix_df$symbol),]  # 删除symbol重复的条目
rownames(raw_fix_df) <- as.character(raw_fix_df[,ncol(raw_fix_df)])  # 转换行名为symbol


# -------------------------------------------------------------------->>>>>>>>>>
# 2.提取表达矩阵,创建分组信息
# -------------------------------------------------------------------->>>>>>>>>>
raw_fix_df <- as.data.frame(raw_fix_df) 
exprSet <- raw_fix_df[,c(8:11)]  # 提取表达矩阵
colnames(exprSet) <- c("Ctrl_rep1","Ctrl_rep2","METTL3_KD_rep1","METTL3_KD_rep2")  # 转换列名
# write.table(exprSet,file = 'all.id.fix.tsv',col.names = T,row.names = T,quote = F,sep = "\t")
group_list <- data.frame(condition = c(rep("ctrl",2), rep("KD",2)))

rownames(group_list) <- colnames(exprSet)
save(exprSet,group_list,file = 'exprSet.Rdata')
```

### 9-3 DEG run

```R
#################################################################################################
# Author Bingjie Fang
# E-mail fbj1001@126.com
# DEG
#################################################################################################

rm(list=ls())
options(stringsAsFactors = F)
load(file = 'exprSet.Rdata')
group_list
#group_list=relevel(group_list,ref = 'ctrl')
source('run_DEG_RNA-seq.R')
run_DEG_RNAseq(exprSet,group_list,
               g1="ctrl",g2="KD",
               pro='m6A')

```

#### DESeq2-MAplot

![image-20210116150433648](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210116150433648.png)

#### DESeq2-heatmap

![image-20210116171828115](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210116171828115.png)

#### DEseq2_volcano

![image-20210116171901964](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210116171901964.png)

#### edgeR-BCV

![image-20210116162919370](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210116162919370.png)

#### edgeR-MeanVar

![image-20210116163449435](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210116163449435.png)

#### edgeR-MAplot

![image-20210116164216454](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210116164216454.png)

#### edgeR-heatmap

![image-20210116172005924](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210116172005924.png)

#### edgeR-volcano

![image-20210116172040220](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210116172040220.png)

## 十、GO和KEGG分析

### 10-1 构造一个run_go，run_kegg函数

```r
################################################################################################
# Author Bingjie Fang
# E-mail fbj1001@126.com
# run-GO  run-KEGG run-GSEA
################################################################################################
rm(list=ls())
options(stringsAsFactors = F)

### GO database analysis 
### 做GO数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
run_go <- function(gene_up,gene_down,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)
  # 因为go数据库非常多基因集，所以运行速度会很慢。
  if(T){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',ont ))
        ego <- enrichGO(gene          = gene,
                        #universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = "SYMBOL" ,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file =paste0(pro, '_go_enrich_results.Rdata'))
    
  }
  # 只有第一次需要运行，就保存结果，下次需要探索结果，就载入即可。
  load(file=paste0(pro, '_go_enrich_results.Rdata'))
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0(pro, '_dotplot_',n1[i],'_',n2[j],'.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }  
}

## KEGG pathway analysis
### 做KEGG数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
run_kegg <- function(gene_up,gene_down,geneList=F,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  ###   over-representation test
  # 下面把3个基因集分开做超几何分布检验
  # 首先是上调基因集。
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk=kk.up
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keytype='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.up.csv'))
  
  # 然后是下调基因集。
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  kk=kk.down
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keytype='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.down.csv'))
  
  # 最后是上下调合并后的基因集。
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  kk=kk.diff
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keytype='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.diff.csv'))
  
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.01,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.01,];up_kegg$group=1
  #画图设置
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg,filename = paste0(pro,'_kegg_up_down.png') )
  
  if(geneList){
    ###  GSEA 
    ## GSEA算法跟上面的使用差异基因集做超几何分布检验不一样。
    kk_gse <- gseKEGG(geneList     = geneList,
                      organism     = 'hsa',
                      nPerm        = 1000,
                      minGSSize    = 20,
                      pvalueCutoff = 0.9,
                      verbose      = FALSE)
    head(kk_gse)[,1:6]
    gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
    gseaplot(kk_gse, 'hsa04110',title = 'Cell cycle') 
    kk=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keytype='ENTREZID')
    tmp=kk@result
    write.csv(kk@result,paste0(pro,'_kegg.gsea.csv'))
    
    
    
    down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
    up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
    
    g_kegg=kegg_plot(up_kegg,down_kegg)
    print(g_kegg)
    ggsave(g_kegg,filename = paste0(pro,'_kegg_gsea.png'))
    # 
  }
  
}

kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}
```

### 10-2 run GO/KEGG

```r
#################################################################################################
# Author Bingjie Fang
# E-mail fbj1001@126.com
# GO/KEGG/GSEA
#################################################################################################

rm(list=ls())
options(stringsAsFactors = F)

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

source('run_DEG_RNA-seq.R')
load(file = 'm6A_DEG_results.Rdata')

deg=getDEGs(DEG_DEseq2,DEG_edgeR,thre_logFC=1,thre_p=0.05)

gene_up= bitr(unique(deg$up), fromType = "SYMBOL",
              toType = c( "ENTREZID"),
              OrgDb = org.Hs.eg.db)[,2] 

gene_down= bitr(unique(deg$down), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)[,2] 
gene_diff=c(deg$up,deg$down) 

save(dds,deg,DEG_DEseq2,DEG_edgeR,exprSet,group_list,gene_diff,gene_down,gene_up,
     draw_h_v,getDEGs,run_DEG_RNAseq,file = 'kegg_go_pre.Rdata')

source('kegg_and_go_up_and_down.R', spaced = F)
load('kegg_go_pre.Rdata')
#######################################################################################
#######################################################################################

run_go(gene_up,gene_down,pro='m6A')


run_kegg(gene_up,gene_down,pro='m6A')

go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
ggsave('gene_up_GO_all_barplot.png')
go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_down_GO_all_barplot.png')
```

#### m6A_dotplot_gene_up_BP

![image-20210117153210395](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117153210395.png)

#### m6A_dotplot_gene_up_CC

![image-20210117153327148](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117153327148.png)

#### m6A_dotplot_gene_up_MF

![image-20210117153351048](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117153351048.png)

#### m6A_dotplot_gene_down_BP

![image-20210117153414162](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117153414162.png)

#### m6A_dotplot_gene_down_CC

![image-20210117153444226](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117153444226.png)

#### m6A_dotplot_gene_down_MF

![image-20210117153505751](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117153505751.png)

#### m6A_dotplot_gene_diff_BP

![image-20210117153527664](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117153527664.png)

#### m6A_dotplot_gene_diff_CC

![image-20210117153552719](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117153552719.png)

#### m6A_dotplot_gene_diff_MF

![image-20210117153615542](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117153615542.png)

#### up_barplot

![image-20210117155528722](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117155528722.png)

#### down_barplot

![image-20210117155632096](https://github.com/fangbingjie/gitbio/blob/master/figure/image-20210117155632096.png)
