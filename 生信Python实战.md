# 生信Python实战

[TOC]

## 1、对FASTQ的操作

- 5,3段截掉几个碱基
- 序列长度分布统计
- FASTQ 转换成 FASTA
- 统计碱基个数及GC%

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File  : Fastq_analysis.py
@Author: BJ_Fang
@Date  : 2021/1/10 
@Desc  : fastq文件数据处理工具
#   1.cut off 5',3' bases(input the num)
#   2.the length distribution of reads
#   3.FASTQ to FASTA
#   4.the number of bases and count the ratio of GC
"""


from optparse import OptionParser
import sys
import filetype
import gzip

args = sys.argv


def cut_fastq(filename, cut_5=0, cut_3=0):
    """
    从5'端,3'端截掉几个碱基
    :param:filename  -- the file which will be dealed with
    :param:cut_5  -- 5' bases number which will be cutted off
    :param:cut_3  -- 3' bases number which will be cutted off
    """

    if cut_3 < 0:
        raise ValueError('cut_3: %s < 0' % cut_3)
    if cut_5 < 0:
        raise ValueError('cut_5: %s < 0' % cut_5)

    kind = filetype.guess(filename)
    if kind is None:
        with open(filename) as f:
            while f:
                line_1 = f.readline().strip("\n")
                if not (line_1 and f):
                    break
                line_2 = f.readline().strip("\n")
                line_3 = f.readline().strip("\n")
                line_4 = f.readline().strip("\n")
                line2 = line_2[cut_5:-cut_3]
                line4 = line_4[cut_5:-cut_3]
                print(line_1+"\n"+line2+"\n"+line_3+"\n"+line4)
    elif kind.extension == "gz":
        with gzip.open(filename, 'rb') as gf:
            while gf:
                line_1 = gf.readline().decode("utf8").strip("\n")
                if not (line_1 and gf):
                    break
                line_2 = gf.readline().decode("utf8").strip("\n")
                line_3 = gf.readline().decode("utf8").strip("\n")
                line_4 = gf.readline().decode("utf8").strip("\n")
                line2 = line_2[cut_5:-cut_3]
                line4 = line_4[cut_5:-cut_3]
                print(line_1+"\n"+line2+"\n"+line_3+"\n"+line4)


def reads_dis(filename):
    """
    统计fastq文件reads长度分布
    parma:filename -- a fastq format file which will be statistics
    """

    read_dis = {}
    kind = filetype.guess(filename)
    if kind is None:
        with open(filename) as f:
            while f:
                line_1 = f.readline().strip("\n")
                if not (line_1 and f):
                    break
                line_2 = f.readline().strip("\n")
                line_3 = f.readline().strip("\n")
                line_4 = f.readline().strip("\n")
                length = len(line_2)
                if not length in read_dis:
                    read_dis[length] = 0
                read_dis[length] += 1
    elif kind.extension == "gz":
        with gzip.open(filename, 'rb') as gf:
            while gf:
                line_1 = gf.readline().decode("utf8").strip("\n")
                if not (line_1 and gf):
                    break
                line_2 = gf.readline().decode("utf8").strip("\n")
                line_3 = gf.readline().decode("utf8").strip("\n")
                line_4 = gf.readline().decode("utf8").strip("\n")
                length = len(line_2)
                if not length in read_dis:
                    read_dis[length] = 0
                read_dis[length] += 1
    for length, num in read_dis.items():
        print("length : %s \t Number : %s" % (length, num))
    return reads_dis


def fastq2fasta(filename):
    """
    将fastq格式文件转换为fasta格式
    parma:filename  -- input the FASTQ format file
    the line_1 is ID and the line_2 is bases
    """

    kind = filetype.guess(filename)
    if kind is None:
        with open(filename) as f:
            while f:
                line_1 = f.readline().strip("\n")
                if not (line_1 and f):
                    break
                line_2 = f.readline().strip("\n")
                line_3 = f.readline().strip("\n")
                line_4 = f.readline().strip("\n")
                print(">", line_1)
                print(line_2)
    elif kind.extension == "gz":
        with gzip.open(filename,'rb') as gf:
            while gf:
                line_1 = gf.readline().decode("utf8").strip("\n")
                if not (line_1 and gf):
                    break
                line_2 = gf.readline().decode("utf8").strip("\n")
                line_3 = gf.readline().decode("utf8").strip("\n")
                line_4 = gf.readline().decode("utf8").strip("\n")
                print(">", line_1)
                print(line_2)


def fq_length(filename):
    """
    统计fastq文件中碱基数量(reads总长度)
    parma:the fastq filename
    """

    length = 0
    kind = filetype.guess(filename)
    if kind is None:
        with open(filename) as f:
            while f:
                line_1 = f.readline().strip("\n")
                if not (line_1 and f):
                    break
                line_2 = f.readline().strip("\n")
                line_3 = f.readline().strip("\n")
                line_4 = f.readline().strip("\n")
                length += len(line_2)
    elif kind.extension == "gz":
        with gzip.open(filename, 'rb') as gf:
            while gf:
                line_1 = gf.readline().decode('utf8').strip("\n")
                # line_1 = [x.decode('utf8').strip() for x in gf.readlines()]
                if not (line_1 and gf):
                    break
                line_2 = gf.readline().decode("utf8").strip("\n")
                line_3 = gf.readline().decode("utf8").strip("\n")
                line_4 = gf.readline().decode("utf8").strip("\n")
                length += len(line_2)
    print(length)
    return length


def gc_count(filename):
    """
    统计fastq文件中reads的GC含量
    param:filename  -- the input FASTQ format file

    """
    count_gc = 0
    kind = filetype.guess(filename)
    if kind is None:
        with open(filename) as f:
            while f:
                line_1 = f.readline().strip("\n")
                if not (line_1 and f):
                    break
                line_2 = f.readline().strip("\n")
                line_3 = f.readline().strip("\n")
                line_4 = f.readline().strip("\n")
                count_gc += line_2.count("G") + line_2.count("C")
    elif kind.extension == "gz":
        with gzip.open(filename, 'rb') as gf:
            while gf:
                line_1 = gf.readline().decode('utf8').strip("\n")
                if not (line_1 and gf):
                    break
                line_2 = gf.readline().decode("utf8").strip("\n")
                line_3 = gf.readline().decode("utf8").strip("\n")
                line_4 = gf.readline().decode("utf8").strip("\n")
                count_gc += line_2.count("G") + line_2.count("C")
    ratio = 0
    total = fq_length(filename)
    ratio = (count_gc * 1.0) / total
    print("GC : %s" % (ratio))
    return ratio


def main(args):
    parser = OptionParser()
    parser.add_option("-f", "--fastq", dest="filename", help=".fastq format file or .fastq.gz format file",
                      metavar="FILE")
    parser.add_option("--cut", "--cut-fastq", dest="cut", help="cut fastq file", action="store_true", default=False)
    parser.add_option("-5", "--cut-5", dest="cut_5", type=int, help="cut fastq N(N>0) bases from 5'", metavar="INT",
                      default=0)
    parser.add_option("-3", "--cut-3", dest="cut_3", type=int, help="cut fastq N(N>0) bases from 3'", metavar="INT",
                      default=0)
    parser.add_option("--len-dis", dest="len_dis", help="sequence length distribution", action="store_true",
                      default=False)
    parser.add_option("--fasta", "--fastq2fasta", dest="fasta", help="convert fastq to fasta", action="store_true",
                      default=False)
    parser.add_option("--gc", "--count-gc", dest="gc", help="count the GC%", action="store_true", default=False)
    parser.add_option("--len", "--fastq-length", dest="fq_len", help="total reads number", action="store_true",
                      default=False)
    (options, args) = parser.parse_args()
    if not options.filename:
        parser.print_help()
    filename = options.filename
    if options.cut:
        cut_5 = options.cut_5
        cut_3 = options.cut_3
        cut_fastq(filename, cut_5, cut_3)
    if options.len_dis:
        reads_dis(filename)
    if options.fq_len:
        fq_length(filename)
    if options.fasta:
        fastq2fasta(filename)
    if options.gc:
        gc_count(filename)


if __name__ == '__main__':
    main(args)

```

## 2、对FASTA的操作

- 取互补序列
- 取反向序列
- 取反向互补序列
- DNA to RNA
- 大小写字母形式输出
- 按照序列长度/名字排序

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File  : Fasta_analysis.py
@Author: BJ_Fang
@Date  : 2021/1/10 
@Desc  : fasta文件数据处理工具
"""


import sys
import re
from collections import OrderedDict
from optparse import OptionParser

args = sys.argv


def complement_seq(filename):
    """
        取互补序列
        :param filename: fastq filename
        :return:

    """

    tmp_dic = {'A': 'T','G': 'C','C': 'G','T': 'A','N': 'N'}
    seq = OrderedDict()

    with open(filename) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                seq_id = line
                seq[seq_id] = ''
            else:
                line = line.upper()
                for i in line:
                    seq[seq_id] += tmp_dic[i]

        for ID,sequence in seq.items():
            print('%s\n%s' % (ID,sequence))


def reverse_seq(filename):
    """
        取反向序列
        :param filename: fastq filename
        :return:

    """

    with open(filename) as f:
        while True:
            line_1 = f.readline().rstrip()
            if not line_1:
                break
            line_2 = f.readline().rstrip()
            line_2 = line_2[::-1]
            print(line_1)
            print(line_2)


def dna_2_rna(filename):
    """
        DNA to RNA
        :param filename: fastq filename
        :return:

    """

    with open(filename) as f:
        while True:
            line_1 = f.readline().rstrip()
            if not line_1:
                break
            line_2 = f.readline().rstrip().upper()
            line_2 = line_2.replace('T','U')
            print(line_1)
            print(line_2)


def format_output(filename,**kwargs):
    """
        大小写字母形式输出
        :param filename: fastq filename
        :param out_filename: output
        :param kwargs: variable args
        :return:

    """
    print("format_output")
    capital = kwargs["capital"] if "capital" in kwargs else True
    # 将if...else语句缩减为单一的条件表达式，语法为expression1 if A else expression2，如果A为True，条件表达式的结果为expression1，否则为expression2
    with open(filename) as f:
        while True:
            line_1 = f.readline().rstrip()
            if not line_1:
                break
            line_2 = f.readline().rstrip()
            if capital:
                line_2 = line_2.upper()
            else:
                line_2 = line_2.lower()
            print(line_2)


def sort_seq(filename,sort_type):
    with open(filename) as f:
        fasta = {}

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                ID = line
                fasta[ID] = ''
            else:
                fasta[ID] += line

        if sort_type == 'id':
            fasta = sorted(fasta.items(),key=lambda i: i[0])
        elif sort_type == 'len':
            fasta = sorted(fasta.items(),key=lambda i: len(i[1]))
        else:
            fasta = fasta.items()

        for k,v in fasta:
            print('%s\n%s' % (k,v))


def main(args):
    parser = OptionParser()

    parser.add_option("-f","--fasta",dest="filename",help="fasta filename",metavar="FILE")


    parser.add_option("--com","--complent_seq",dest="complement",help="complent fasta sequence",action="store_true",
                      default=False)

    parser.add_option("--rev","--reverse_seq",dest="reverse",help="reverse fasta sequence",action="store_true",
                      default=False)

    parser.add_option("--rna","--dna2rna",dest="rna",help="convert DNA to RNA",action="store_true",default=False)

    parser.add_option("--upper",dest="format_upper",help="output upper fasta sequence",action="store_true",
                      default=False)

    parser.add_option("--lower",dest="format_lower",help="output lower fasta sequence",action="store_true",
                      default=False)

    parser.add_option("--sort_type",dest="sort_type",help="output fasta by sort =id or =len",default=False)

    (options,args) = parser.parse_args()

    if not options.filename:
        parser.print_help()

    filename = options.filename

    if options.complement:
        complement_seq(filename)

    if options.reverse:
        reverse_seq(filename)

    if options.rna:
        dna_2_rna(filename)

    if options.format_upper:
        format_output(filename,capital=True)

    if options.format_lower:
        format_output(filename,capital=False)

    if options.sort_type:
        sort_seq(filename,options.sort_type)


if __name__ == '__main__':
    main(args)
    
```

### 2-1 fasta操作补充

#### 2-1-1 fasta文件读取

```python
#FASTA文件读取
def read_fasta(file_path=r""):
    """
    Loading FASTA file and return a iterative object
    """

    line = ""

    try:
        fasta_handle = open(file_path,"r")
    except:
        raise IOError("Your input FASTA file is not right!")

    # make sure the file is not empty
    while True:
        line = fasta_handle.readline()
        if line == "":
            return
        if line[0] == ">":
            break

    # when the file is not empty, we try to load FASTA file
    while True:
        if line[0] != ">":
            raise ValueError("Records in Fasta files should start with '>' character")
        title = line[1:].rstrip()
        lines = []
        line = fasta_handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = fasta_handle.readline()

        yield title,"".join(lines).replace(" ","").replace("\r","")

        if not line:
            return

    fasta_handle.close()
    assert False, "Your input FASTA file have format problem."

for header,seq in read_fasta(file_path=r"…………"):
    print (header+"\n"+seq)
    
```

#### 2-1-2 统计碱基数量

```python
# 以chrM.fa为例，定义变量用来保存chrM的所有序列
chrM_genome = ""

with open("./chrM.fa","r") as input_file:	 
	for line in input_file: 				  # 使用for循环进行读取文件，每次读取1行
		if line[0] != ">": 					  # 当读取到的内容不是标题的时候，就是序列内容
			line = line.strip().upper() 	  # strip去掉行末的换行符，upper是把所有的序列都转换成大写 
			chrM_genome = chrM_genome + line  # 将新读取的行添加到之前的序列

# 统计A，T，G，C的数量
print("A count is {0}".format(chrM_genome.count("A")))
print("T count is {0}".format(chrM_genome.count("T")))
print("C count is {0}".format(chrM_genome.count("C")))
print("G count is {0}".format(chrM_genome.count("G")))

# 统计线粒体参考基因组的长度
print("chrM total length is: {0} bp".format(len(chrM_genome)))

```

#### 2-1-3 人类hg19基因组统计

```python
#生成hg19染色体的名字和序列对应的字典，调用上述的read_fasta函数，UCSC下载hg19_only_chromosome.fa 
hg19_genome = {}
for chr_name,chr_seq in read_fasta(file_path=r"D:/data_file/hg19_only_chromosome.fa"):
    hg19_genome[chr_name] = chr_seq

#求各染色体长度,获取全基因组长度
for index in range(1,26):
    if index <= 22:
        chr_name = "chr{0}".format(index)
    elif index == 23:
        chr_name = "chr X"
    elif index == 24:
        chr_name = "chr Y"
    elif index == 25:
        chr_name = "chr M"
    print (chr_name，len(hg19_genome.get(chr_name )))  # hg19 各染色体长度

#统计hg19基因组全长
hg19_total_len = 0 
for index in range(1,26):
    if index <= 22:
        chr_name = "chr{0}".format(index)
    elif index == 23:
        chr_name = "chr X"
    elif index == 24:
        chr_name = "chr Y"
    elif index == 25:
        chr_name = "chr M"
    hg19_total_len += len(hg19_genome.get(chr_name ))  #基因组全长
print(hg19_total_len)

#统计染色体碱基数量
hg19_genome["chr22"].count("A")  # 第22条染色体碱基A的含量

#提取基因组特定区域序列，输出，支持正负链
chr_name = "chr1"
chr_start = 10000000
chr_end = 10000100
hg19_genome[chr_name][chr_start-1:chr_end-1].upper() #取特定区域序列（正链）

#取反向互补序列（string库）
import string
translate_table = string.maketrans("ATGC","TACG")  # 制作转义表
a = hg19_genome[chr_name][chr_start-1:chr_end-1].upper()
a.translate(translate_table)  #互补序列
a.translate(translate_table)[::-1]  # 反向互补序列(负链)

#反向互补函数
def com_rev(input_seq):
    translate_table = string.maketrans("ATGC","TACG")
    return input_seq.upper().translate(translate_table)[::-1]
```



#### 2-1-4 fastq转换为fasta，按指定长度输出

```python
output_file = open(r"D:\data_file\test.fa","w")

with open(r"D:\data_file\test.fastq","r") as input_fastq:
    for index,line in enumerate(input_fastq):
        if index % 4 == 0:
            output_file.write(">" + line.strip()[1:] + "\n")
        elif index % 4 == 1:
            for i in range(0,len(line.strip()),40):  # 40个字符串为一行输出
                output_file.write(line.strip()[i:i+40] + "\n")
        elif index % 4 == 2:
            continue
        elif index % 4 == 3:
            continue
output_file.close()
```



#### 2-1-5 统计人类基因组多少个基因及exon长度（UCSC注释文件）

```python
# UCSC下载人类参考转录组注释 hg19_refseq_gene_table.txt
# 一个基因多个转录本，选取最长的转录本

gene_exon_len_dict= {}
exon_total_len = 0

with open(r"D:\data_file\hg19_refseq_gene_table.txt","r") as input_file
	header =  input_file.readline()  # 第一行读进,循环在第二行开始
    for line in input_file:
        line = line.strip().split("\t")
        exon_count = int(line[8])  # exon数量在第9列
        exon_start = map(int,line[9].strip(",").split(","))  # exon_start在第10列，逗号分隔
        exon_end = map(int,line[10].strip(",").split(","))  # exon_end在第11列，逗号分隔
        for index in range(exon_count):
            exon_total_len = exon_total_len + exon_end[index] - exon_start[index]
        gene_id = line[12]
        if gene_exon_len_dict.get(gene_id) is None:
            gene_exon_len_dict[gene_id] = exon_total_len
        else:
            if gene_exon_len_dict.get(gene_id) < exon_total_len
            	gene_exon_len_dict.get[gene_id] = exon_total_len
        
print len(gene_exon_len_dict)  # 总共多少基因  2万7
#############################
exon_total = 0
for gene_id in gene_exon_len_dict.keys():
    exon_total = exon_total + gene_exon_len_dict[gene_id]
print (exon_total) # exon全长
print (exon_total * 1.0 / hg19_total_len * 100)  # exon占基因组比例（2.54%）
```

## 3、实战

#### 3-1 统计人类外显子长度

```python
#坐标文件（NCBI注释文件）
wget ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
    
import sys
import re

args=sys.argv


filename=args[1]
exon_length=0
aDict={}
with open(filename) as f:
      for line in f:
          if line.startswith("#"):
              continue
          lineL=line.strip().split("\t")
          exon_position=lineL[-2] #取出倒数第二列，坐标列
          if exon_position=="-":  #有的基因没有外显子的坐标，用-代替的，所以这行就要除掉，不然会报错
              continue
          exon_position=re.sub("\[|\]","",exon_position) #把坐标列的[]去除，注意正则表达式的用法
          exonL=exon_position.split(",")  
          for exon in exonL:
              exonS=lineL[0]+":"+exon    #有点基因会有相同坐标的外显子，所以要去除这一部分，注意要加上染色体的编号，染色体不同而坐标一样就没事
              if exonS not in aDict:    #如果坐标没有在字典，即第一次出现，就将其放入字典，并继续操作。
                   aDict[exonS]=1
                   exon_pL=exon.split("-")
                   exon_start=int(exon_pL[0].strip())
                   exon_end=int(exon_pL[1].strip())
                   exon_length+=exon_end-exon_start
print(exon_length)
```

#### 3-2 统计hg38 每条染色体长度，N的含量，GC含量

```python
import sys
import collections

args=sys.argv


filename=args[1]
count_ATCG=collections.OrderedDict() #构建有顺序的字典
Bases=["a","t","g","c","n"]
for line in open(filename):
   if line.startswith(">"):
      chr_id = line[1:]
      count_ATCG[chr_id] = {}
      for base in Bases:
        count_ATCG[chr_id][base] = 0 ##字典中，chr_id属于key，base和0共同构成#value，而value中，base又成了key，0成了value
   else:
      line = line.lower()
      for base in Bases:
        count_ATCG[chr_id][base] += line.count(base)

for chr_id, ATCG_count in count_ATCG.items():
    count_sum = sum(ATCG_count.values())
    count_GC = ATCG_count['g'] + ATCG_count['c']
    print(chr_id)
    for base in Bases:
        print("%s: %s" % (base, ATCG_count[base])) #统计碱基数量

    print("Sum: %s" % count_sum)  # 统计染色体总长度
    print("N %%: %.2f%%" % (ATCG_count['n']*100.00/count_sum))  # 统计染色体N的含量
    print("GC %%: %.2f%%" % (count_GC*100.00/count_sum))  # 统计GC含量

```

#### 3-3 统计hg38 每条染色体基因、转录本的分布

```python
#注释文件下载（ensembl）
wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.87.chr.gtf.gz   
grep -v "#" Homo_sapiens.GRCh38.87.chr.gtf > hg38.gtf


#统计各个染色体基因数量的分布，包括总的基因和各种类型基因数量的分布


import sys
import collections

args=sys.argv


filename=args[1]
count_gene=collections.OrderedDict()
with open(filename) as f:
   for line in f:
      lineL=line.strip().split("\t")
      chr_num="chr"+lineL[0]
      type_info=lineL[2]
      if type_info=="gene":
         All_gene="all"+"_"+type_info
         if chr_num not in count_gene:
            count_gene[chr_num]={}
            count_gene[chr_num][All_gene]=0
         count_gene[chr_num][All_gene]+=1    #计算all gene个数
         descrip_info=lineL[8]
         descrip_info=descrip_info[:-1]      #主要是有的biotype就是最后一部分了，这时候后面得到的b_type包含“;”,所以一开始就把最后的";"去掉
         descrip_infoL=descrip_info.split("; ")
         biotype=descrip_infoL[4]
         biotypeL=biotype.split(" ")
         b_type=biotypeL[1]
         b_type=eval(b_type)
         if b_type not in count_gene[chr_num]:
              count_gene[chr_num][b_type]=0
         count_gene[chr_num][b_type]+=1
for chr_num,countAll in count_gene.items():
   for k,v in countAll.items():
      print(chr_num,k,v)
```

```r
#r作图统计 all_gene 分布
rm(list = ls())
library(ggplot2)
data <- read.table("count_gene.txt")
chr_name <- data$V1
type <- data$V2 
number <- data$V3
p <- ggplot(data,aes(x=chr_name,y=number))+
            theme_bw()+
            theme(axis.text.x=element_text(angle=90,size =8,vjust = 0.5,hjust = 0.5),
            panel.grid = element_blank())+
            geom_bar(aes(fill=chr_name),stat = "identity",position="dodge",width=0.8)
p
```

![image-20210110193810058](D:%5CDesktop%5C%E7%94%9F%E4%BF%A1%E6%8A%80%E8%83%BD%E6%A0%91%5C6.%E7%94%9F%E4%BF%A1%E5%BA%94%E8%AF%A5%E8%BF%99%E6%A0%B7%E5%AD%A6python%5C%E7%94%9F%E4%BF%A1Python%E5%AE%9E%E6%88%98.assets%5Cimage-20210110193810058.png)

```python
#统计一下基因转录本数量的分布
import sys
import collections

args=sys.argv


filename=args[1]
count_transcript=collections.OrderedDict()
with open(filename) as f:
   for line in f:
      lineL=line.strip().split("\t")
      chr_num="chr"+lineL[0]
      type_info=lineL[2]
      if type_info=="transcript":
         if chr_num not in count_transcript:
            count_transcript[chr_num]={}
         descrip_info=lineL[8]
         descrip_infoL=descrip_info.split("; ")
         gene_ID=descrip_infoL[0]
         gene_IDL=gene_ID.split(" ")
         gene_ID=gene_IDL[1]
         gene_ID=eval(gene_ID)
         if gene_ID not in count_transcript[chr_num]:
              count_transcript[chr_num][gene_ID]=0
         count_transcript[chr_num][gene_ID]+=1
for chr_num,countAll in count_transcript.items():
   for k,v in countAll.items():
      print(chr_num,k,v)
```

#### 3-4 多个同样的行列式文件合并

```python
import glob
import collections

mydict=collections.OrderedDict()
list_dirs=glob.glob("./*.txt")  # 当前文件下的文件路径读取到列表中
for i in list_dirs:
   for line in open (i):
      array=line.strip().split("\t")
      if array[0] not in mydict:
           mydict[array[0]]=[array[1]]    #注意array[1]外面的中括号
      else:
           mydict[array[0]].append(array[1])    #字典中，对一个key增加多个value的方法

for gene_name in mydict:
    print("%s\t%s"%(gene_name,"\t".join(mydict[gene_name])))
```

#### 3-5 根据GTF画基因的多个转录本结构

```python
#先从hg38的gtf中提取"ANXA1"基因
grep '"ANXA1"' hg38.gtf >ANXA1.gtf


import collections
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches

matplotlib.style.use("ggplot")
gene_cord=collections.OrderedDict()
gene_exon=collections.OrderedDict()
gene_utr=collections.OrderedDict()
gene_start_codon=collections.OrderedDict()
gene_stop_codon=collections.OrderedDict()

with open("ANXA1.gtf") as f:
   for line in f:
      lineL=line.strip().split("\t")
      if lineL[2]=="gene":
          gene_info=lineL[-1]
          gene_infoL=gene_info.split("; ")
          gene_name=gene_infoL[2]
          gene_nameL=gene_name.split(" ")
          gene_name=eval(gene_nameL[1])
          gene_cord[gene_name]=[int(lineL[3]),int(lineL[4])]  #对基因构建字典,key是基因名，value是基因坐标
          gene_exon[gene_name]={}
          gene_utr[gene_name]={}
      elif lineL[2]=="transcript":
          trans_info=lineL[-1]
          trans_infoL=trans_info.split("; ")
          trans_name=trans_infoL[9]
          trans_nameL=trans_name.split(" ")
          trans_name=eval(trans_nameL[1])
          gene_exon[gene_name][trans_name]=[]                #字典中套字典,最终的value先建立为空list
          gene_exon[gene_name][trans_name].append(lineL[6])  #往空list中添加正负链信息的元素
          gene_utr[gene_name][trans_name]=[]
      elif lineL[2]=="exon":
          gene_exon[gene_name][trans_name].extend([int(lineL[3]),int(lineL[4])]) #往list中加exon坐标
      elif lineL[2] in ["five_prime_utr","three_prime_utr"]:
          gene_utr[gene_name][trans_name].extend([int(lineL[3]),int(lineL[4])])  #往list中加UTR坐标

fig= plt.figure()  #建一个画板
ax = fig.add_subplot(111)  #将画板分为1行1列，图像画在从左到右从上到下第一块

trans_num=0
for trans_name,exon in gene_exon[gene_name].items():
    trans_num+=1
    if exon[0]=="+":
      col="limegreen"
    if exon[0]=="-":
      col="magenta"
    for i in range(1,len(exon),2):   #注意！！！python中for i in range(n,m)  范围包括n,不包括m
        #对每个exon画矩形图
        #ax.add_patch(patches.Rectangle((x,y),width,height))
        rect=patches.Rectangle((exon[i],trans_num-0.1),exon[i+1]-exon[i],0.2,color=col,fill=True)
        ax.add_patch(rect)
        #对每个intron画线图
        if i < len(exon)-2:
            #ax.plot([x1,x2],[y1,y2])
            ax.plot([exon[i+1],exon[i+2]],[trans_num,trans_num],color="red",linewidth="0.5")

#对UTR画图
trans_num=0
for trans_name,utr in gene_utr[gene_name].items():
    trans_num+=1
    if utr==[]:
        continue
    else:
        for i in range(0,len(utr),2):
            rect=patches.Rectangle((utr[i],trans_num-0.1),utr[i+1]-utr[i],0.2,color="black")
            ax.add_patch(rect)



ax.set_xlabel(gene_name,color="blue",fontsize=14,fontweight="bold")
ax.yaxis.set_major_locator(ticker.FixedLocator(range(1,trans_num+1)))
ax.set_yticklabels(gene_exon[gene_name].keys())
plt.xlim(gene_cord[gene_name])
plt.ylim([0.5,trans_num+0.5])
plt.show()

fig.savefig('rect.png', dpi=90, bbox_inches='tight')
```

![image-20210110200830100](D:%5CDesktop%5C%E7%94%9F%E4%BF%A1%E6%8A%80%E8%83%BD%E6%A0%91%5C6.%E7%94%9F%E4%BF%A1%E5%BA%94%E8%AF%A5%E8%BF%99%E6%A0%B7%E5%AD%A6python%5C%E7%94%9F%E4%BF%A1Python%E5%AE%9E%E6%88%98.assets%5Cimage-20210110200830100.png)

#### 3-6 下载最新版KEGG信息，格式化好，统计pathway信息

```python
import sys
import collections
import re

args=sys.argv


filename=args[1]
keggD=collections.OrderedDict()
with open(filename) as f:
   for line in f:
     if line.startswith("A"):
         mth=re.search('A(\S+)(\s+)(.+)',line)  #"A(\S+)(\s+)(.+)" 指的是匹配以A开头，后面跟的是一个或多个非空白字符，
         cate1=mth.group(3)                     #再跟着的是一个或多个空白字符，再后面跟一个或多个除换行符外的字符。/表示转义，A没必要转义
         keggD[cate1]=collections.OrderedDict()

     if line.startswith("B"):
         if len(line)== 2:
            continue
         mth=re.search('(\d+)(\s+)(.+)',line)  #'(\d+)(\s+)(.+)'指的是匹配一个或多个数字，后面跟着的是一个或多个空白字符，再后面跟着的是一个或多个
         cate2=mth.group(3)                    # 除换行符外的任意字符
         keggD[cate1][cate2]=collections.OrderedDict()

     if line.startswith("C"):
         mth=re.search('(\d+)(\s+)(.+)',line)  # '(\d+)(\s+)(.+)'同上
         pathID="hsa"+mth.group(1)
         pathname=mth.group(3)
         pathname=re.sub("(\s+)\[(\S+)\]","",pathname)  # (\s+)\[(\S+)\] 匹配一个或多个空白字符，后面跟着的是[,再后面跟着的是一个或多个非空白字符，再跟]
         pathway=pathID+":"+pathname
         keggD[cate1][cate2][pathway]=[[],[]]
     if line.startswith("D"):
         lst=line.split(";")
         geneinfo=lst[0].split("\t")[0]
         mth=re.search("(\d+)(\s+)(\S+)",geneinfo)  #"(\d+)(\s+)(\S+)" 匹配的是一个或多个数字字符，后面跟一个或多个空白字符，再后面跟的是一个或多个非空白字符  gene_ID=mth.group(1)
         gene_symbol=mth.group(3)
         keggD[cate1][cate2][pathway][0].append(gene_ID)
         keggD[cate1][cate2][pathway][1].append(gene_symbol)

for cate1,value1 in keggD.items():
    for cate2,value2 in value1.items():
        for pathway,value3 in value2.items():
            geneIDs=";".join(value3[0])
            genes=";".join(value3[1])
            print(cate1,cate2,pathway,geneIDs,genes,sep="\t")
```

```python
# 统计pathway

import sys

args=sys.argv
filename=args[1]
cnt_pathway=0
all_gene=[]
for line in open(filename):
    if len(line)>3 :
       cnt_pathway += 1

    lineL=line.split("\t")
    gene_symL=lineL[-1].split(";")
    for gene in gene_symL:
        if gene not in all_gene:
           all_gene.append(gene)

print(cnt_pathway,len(all_gene)) 
```

#### 3-7 超几何分布检验

```python
import os
import re
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from collections import OrderedDict

kegg = OrderedDict()

#计算m和k

pop=[]
for line in open("hsa.kegg.txt"):
    lineL=line.strip().split("\t")
    if len(line)>3:     #只有len（line）>3才表示这条是有基因的pathway
        gene_id=lineL[-2]
        pathway=lineL[2]
        gene_idL=gene_id.strip().split(";")
        kegg[pathway]=gene_idL
        for i in gene_idL:
            if i not in pop:
                pop.append(i)
DEG=[]
for line in open("ehbio.DESeq2.all.DE.entrez.txt") :
    line=line.strip()
    if line in pop:
        DEG.append(line)
pop_number=len(pop)
DEG_number=len(DEG)

print("the number of all gene in pathway is :%d" %pop_number )#注意%d而不是d%
print("the number of diff gene in pathway is :%d" %DEG_number)


#超几何分布检验

keggEnrichment=OrderedDict()
for k,v in kegg.items():
    cnt=[]
    pathway_gene_cnt=len(v)
    for i in DEG:
        if i in v:
            cnt.append(i)
    dif_in_pathway=len(cnt)
    if dif_in_pathway == 0:
        print(k)
    else:
        pValue=hypergeom.sf(dif_in_pathway-1,pop_number,pathway_gene_cnt,DEG_number)
        keggEnrichment[k]=[dif_in_pathway-1,pop_number,pathway_gene_cnt,DEG_number,pValue,";".join(cnt)]

#用pandas的dataframe便于操作

keggOUT=pd.DataFrame.from_dict(keggEnrichment,orient="columns",dtype=None)
keggOUT=pd.DataFrame.transpose(keggOUT) #行列的转置
keggOUT.columns=["a","m","n","k","pValue","gene"]
keggOUT=keggOUT.sort_values(by="pValue",axis=0)


#FDR
def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

pValue=keggOUT["pValue"]
fdr=p_adjust_bh(pValue)
fdr=pd.DataFrame(fdr,index=keggOUT.index)
keggOUT.insert(5,"FDR",fdr)

keggOUT.loc[keggOUT["FDR"]<0.05] #筛选FDR<0.05
```

#### 3-8 ID转换

```python
import collections

geneDict=collections.OrderedDict()
with open("geneid2symbol.txt") as fh:
  for line in fh:
    lineL=line.strip().split("\t")
    gene_id=lineL[0]
    symbol=lineL[1]
    geneDict[gene_id]=symbol

for line in open("my_geneID.txt"):
   line=line.strip()
   print(geneDict[line])
```

#### 3-9 根据指定染色体及坐标得到参考碱基

```python
import sys
args=sys.argv
filename=args[1]
import collections
chr_num=args[2]
base_num=args[3]
up_stream=int(base_num)-1
down_stream=int(base_num)+1
num=0
stop=0
with open(filename) as fh:
  chrome=">"+chr_num
  for line in fh:
     line=line.strip()
     if line == chrome :
          num=1
          continue
     if num:
        if line.startswith(">"):
           break
        for i in line:
           if num == up_stream:
              print(i)
              stop=stop+1
           if num == int(base_num):
              print(i)
              stop=stop+1
           if num == down_stream:
              print(i)
              stop=stop+1
           num+=1
     if stop == 3:
         break
```

#### 3-10 根据指定染色体及坐标得到参考碱基

```python
# 下载tss文件
wget http://www.biotrainee.com/jmzeng/tmp/hg38.tss
    
import sys
args=sys.argv
filename1=args[1]
filename2=args[2]
aDict={}

for line in open (filename1):
   lineL=line.strip().split("\t")
   gene_ID=lineL[0]
   chr_name=lineL[1]
   start=lineL[2]
   end=lineL[3]
   if chr_name not in aDict:
     aDict[chr_name]={}
   aDict[chr_name][start,end]=gene_ID


for line in open (filename2):
   lineList=line.strip().split("\t")
   Chr_name=lineList[0]
   Start=int(lineList[1])
   End=int(lineList[2])
   for k1,v1 in aDict[Chr_name].items():
     if int(k1[0])<=Start <= int(k1[1]) or int(k1[0])<=End <= int(k1[1]) or Start <=int(k1[0])<=End or Start <=int(k1[1])<=End :
        print(Chr_name,Start,End,v1)    
```

#### 3-11 文件内容按染色体分开写出

```python
import sys
args=sys.argv
filename=args[1]
aList=[]
for i in range(1,23):
  i=str(i)     #这个转换特别重要，如果不换，存到列表中的将是数字，不是字符串
  aList.append(i)

with open(filename) as fh:
  for line in fh:
     if line.startswith("#"):
         continue
     lineL=line.strip().split("\t")
     chr_num=lineL[0]            #当然这里也是可以考虑转成int，上面的不用转str，但是！！，并不是所有的chr_num都是数字形式的str，所以不能转int，会报错
     if chr_num in aList:
        name=chr_num+".txt"
        f=open(name,"a")     #a很重要，表示写入文件，如果文件存在，则在结尾追加
        f.write(line)    #写入文件
        f.close()   #文件关闭
     else :
        name="else_chr.txt"
        f=open(name,"a")
        f.write(line)
        f.close()
```

#### 3-12 json格式数据的格式化

```python
import sys
import re
args=sys.argv
filename=args[1]
aDict={}
with open (filename) as fh :
   pattern = re.compile("\s+\"(.+)\"\s:\s\"(.+)\",")
   """
      一个或多个空格+"+一个或多个除换行符之外的任意字符+"+一个空格+:+一个空格+"+一个或多个除换行符之外的任意字符+"+,
   这样我只需要将每个line右边的空格换行符去掉，然后对每行去匹配这个pattern就能得到需要整理的行
   """
      for line in fh:
     line=line.rstrip()
     mth=pattern.search(line)
     if mth:           #如果匹配上，mth就不是空
        a = mth.group(1)  #取出第一个括号里的内容
        b = mth.group(2)  #取出第二个括号里的内容
        if a not in aDict:
           aDict[a]=[b]
        else:
           aDict[a].append(b)
for k,v in aDict.items():
    z='\t'.join(v)            #因为v是列表，这一步就是把列表转换成字符串，方便查看
    print(k,z,sep="\t")
```

