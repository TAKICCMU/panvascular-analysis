####更新后的步骤10
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#加载包
library(Fast2TWAS)
library(data.table)

#输入参数-----------------------------------------------------------------------
#读入index文件
setwd("/path/to/panvas_analysis/Genomic_SEM/")
index <- fread("/path/to/panvas_analysis/Genomic_SEM/gwasfile_info.txt")
################################################################################
# RUN LD-SCORE REGRESSION-----------------------------------------------
# 运行多变量 LD-Score 回归，以获得遗传协方差矩阵 (S) 和相应的抽样协方差矩阵 (V)
################################################################################
gwas_files<-paste0("/path/to/panvas_analysis/Genomic_SEM/",index$phenotype,".sumstats.gz")
sample.prev <- c()
for (file in index$filename) {
  # 判断GWAS数据是否是多个GWAS数据META后得到的数据，以及GWAS数据类型是病例对照
  casecontrol <- ifelse(!is.na(index$case[which(file==index$filename)]) & !is.na(index$controls[which(file==index$filename)]),TRUE,FALSE)
  ismeta <- index$ismeta[which(file==index$filename)]
  
  #sample.prev
  if (casecontrol) {
    if (ismeta) {
      single_sample.prev <- 0.5
    }else{
      # 疾病流行率的计算方法是病例数除以参与者总数（即病例 + 对照组）。可能的范围 = 0-1，保留四位有效数
      # single_sample.prev <- index$case[which(file==index$filename)]/(index$case[which(file==index$filename)]+index$controls[which(file==index$filename)])
      single_sample.prev <- round(index$case[which(file==index$filename)]/(index$case[which(file==index$filename)]+index$controls[which(file==index$filename)]), 4)
    }
    if (single_sample.prev < 0 | single_sample.prev > 1) {
      stop("sample.prev计算错误，请确认数据后重新分析")
    }
  }else{
    single_sample.prev <- NA
  }
  sample.prev <- c(sample.prev,single_sample.prev)
}
print(sample.prev)
population.prev<-index$population.prev
print(population.prev)
ld<-"/path/to/panvas_analysis/Genomic_SEM/eur_w_ld_chr/"
wld<-"/path/to/panvas_analysis/Genomic_SEM/eur_w_ld_chr/"
phen<-index$phenotype
print(phen)
compiler::loadcmp("/path/to/panvas_analysis/Genomic_SEM/10.ld-score_regression_pro.Rc", env = environment())

####更新后的步骤11
library(Fast2TWAS)
library(data.table)
################################################################################
# ESTIMATE THE COMMON FACTOR MODEL -------------------------------------
################################################################################
#加载LDSCoutput.RData数据
load("/path/to/panvas_analysis/Genomic_SEM/LDSCoutput.RData")
covstruc<-LDSCoutput
estimation<-"DWLS"
compiler::loadcmp("./11.common_factor_pro.Rc", env = environment())

####更新后的步骤12
library(Fast2TWAS)
compiler::loadcmp("./12.EFA_pro_lib.Rc", env = environment())

#加载LDSCoutput.RData数据
#导入12步的数据
load("/path/to/panvas_analysis/Genomic_SEM/LDSCoutput.RData")

# Exploratory Factor Analysis --------------------------------------------------
#smooth the S matrix for EFA using the nearPD function in the Matrix package. 
Ssmooth<-as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)
# 计算最大可行因子数
p <- ncol(Ssmooth)  
max_factors <- floor((2*p + 1 - sqrt(8*p + 1))/2)  # 最大因子数公式
factors_num <- 1:max_factors
#lapply编写函数factors_num设置为最大因子+1
all <-   lapply(factors_num, function(factors_num) {
  EFA_factor <- try(factanal(covmat = Ssmooth, factors = factors_num, rotation = "promax"),silent = TRUE)
  print(EFA_factor)
  return(EFA_factor)
})

print(all)

#保存数据
save(all,file = "EFA.RData")

# 开始KMO分析-------------------------------------------------------------------
compiler::loadcmp("/path/to/panvas_analysis/Genomic_SEM/EFA_pro.Rc", env = environment())
# KMO 值的解释标准（Kaiser, 1974） 
# KMO ≥ 0.90：非常适合（marvelous） 
# 0.80 ≤ KMO < 0.90：好（meritorious） 
# 0.70 ≤ KMO < 0.80：一般（middling） 
# 0.60 ≤ KMO < 0.70：勉强可以（mediocre） 
# 0.50 ≤ KMO < 0.60：较差（miserable） 
# KMO < 0.50：不适合因子分析

# 开始并行分析---###结果强烈支持二因子模型
paLDSC(
  S = LDSCoutput$S_Stand,
  V = LDSCoutput$V_Stand,
  r = 500,
  p = .95,
  diag = FALSE,
  fa = TRUE,
  fm = "minres",
  nfactors = 1,
  save.pdf = TRUE
)

#Specify the Genomic confirmatory factor model----------------------------------
CFAofEFA <- '
F1 =~ NA*AA + AAA + TAAD
F2 =~ NA*CAD + Stroke + PAD + CA + CAVS
G =~ NA*F1 + l2*F2
G ~~ 1*G
F1 ~~ u*F1
F2 ~~ u*F2
F1 ~~ 0*F2

AAA ~~ 1e-6*AAA

'

#run the model
CFA <- usermodel(
  LDSCoutput,
  estimation = "DWLS",
  model = CFAofEFA,
  CFIcalc = TRUE,
  std.lv = TRUE,
  imp_cov = TRUE,
  fix_resid = FALSE,
  toler = 1e-20
)

# 保存两个表格
write.csv(CFA$modelfit, "EFA_pro_STable4a.txt", row.names = FALSE)
write.csv(CFA$results, "EFA_pro_STable4b.txt", row.names = FALSE)
# 绘图 -------------------------------------------------------------------------
compiler::loadcmp("./usermodel_path_diagram_plot.Rc", env = environment())

####更新后的步骤13
library(GenomicSEM)

#输入参数-----------------------------------------------------------------------
#读入index文件
index <- data.table::fread("/path/to/panvas_analysis/Genomic_SEM/gwasfile_info.txt")

################################################################################
# Prepare the summary statistics for GWAS---------------------------------------
# 为了运行多变量GWAS, SNP和相应的SNP效应需要跨表型矫正，以便在所有情况下相同的等位基因是参考等位基因。
################################################################################

#GWAS数据路径
gwas_files<-paste0("/path/to/panvas_analysis/Genomic_SEM/",index$phenotype,".txt")

#用于计算跨性状 SNP 方差的参考文件
#下载地址：https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v
ref <- "/path/to/panvas_analysis/Genomic_SEM/reference.1000G.maf.0.005.txt.gz"

#性状名称 
phen<-index$phenotype

#SE 是否采用对数标度：简单说就是判断SE是否是Beta值的SE
# 1. 请注意，对于许多性状，即使回归效应是以几率比（ORs）的形式报告的，SE 也将以对数标度的形式报告。
# 2. 这一信息通常可以在汇总统计的 README 文件中找到。
# 3. 如果该参数被列为 FALSE，且所有其他附加参数均为 FALSE，则函数期望得到 OR 的标准误差。
# 4. 函数将根据效应大小中位数自动判断效应列是 OR 还是 logistic 标度。
se.logit=rep(TRUE, length(phen))

#表型是否是使用观察到的最小二乘 (OLS；即线性) 估计法进行分析的连续性状：
# 1. 如果没有提供任何信息，则该函数默认为 FALSE；
# 2. 我们为这个参数和所有其他默认为 NULL 的参数写上 “OLS=NULL”，以便在需要时清楚参数的名称；
# 3. 如果该参数实际上对所有性状都为 FALSE，则无需将其写出。
OLS=NULL

#逻辑向量
# 1. 表型是否为二分结果，且在GWAS汇总统计文件中只有Z统计量；
# 2. 或表型是否为使用OLS分析的二分结果；
# 3. 如果没有提供任何信息，函数默认为 NULL。
linprob=NULL

# 用户提供的 N，按文件参数的性状顺序排列：
# 1. 如果没有提供任何信息，函数默认为 NULL；
# 2. 在使用 linprob 参数反推出逻辑贝塔值时，需要提供参与 GWAS 的队列的有效样本量总和；
# 3. 对于 OLS 变换，这应该是总样本量；
# 4. 如果提供的是某些性状的样本量，而不是其他性状的样本量，那么对于用户不希望提供样本量的性状，可以使用 NA；
# 5. 如果汇总统计文件中包含样本大小列，那么如果用户希望使用特定于 SNP 的样本大小来执行重定标，也可以列出 NA；
# 6. 不过，我们要再次指出，这应该是二分性状的有效样本量之和；
# 7. 在本例中，所有性状都是病例/对照性状，并且已经报告了几率比或 logistic betas，或者在酒精使用障碍的情况下，只提供了 Z 统计量，但在 GWAS 摘要数据中存在有效 N 的总和；
# 8. 因此，对于这些性状，该参数被列为 NULL。
N=NULL

#如果用户知道 GWAS 是在一个已经标准化的表型上运行的，那么这个参数可以用来直接提供 GWAS 统计摘要中连续性性状的 betas 的列名:
# 1. 在这种情况下，软件包将直接使用这些 betas，除了将符号与参考等位基因对齐之外，不会执行任何进一步的转换；
# 2. 用户使用该参数的主要原因是连续性状的汇总统计文件中没有 MAF；
# 3. 更具体地说，当样本内 MAF 不可用时，函数会使用参考 MAF 来反推出部分标准化的贝塔值，而样本内 MAF 和参考 MAF 之间的不匹配可能会导致汇总统计中贝塔值的相对缩放出现极小量的偏差。
betas=NULL

#用作输入质量下限的数值（INFO），即imputation后SNP推算质量的纳入标准，默认为0.6
info.filter=0.6

#作为次要等位基因频率下限的数值，默认为0.01
maf.filter=0.01

#输出中是否包含插入缺失（indels），默认 = FALSE
keep.indel=FALSE

#表示 sumstats 是以并行还是串行方式处理汇总统计文件。默认为 “FALSE”，表示串行处理。 
parallel=TRUE

#表示并行运行时要使用多少个内核：
# 1. 默认值为 NULL，在这种情况下，sumstats 将使用比本地环境可用内核总数少 1 的内核；
# 2. 或者当可用内核数超过性状数时，使用与性状数相同的内核cores。
cores = NULL

#表示是否应从输出中删除strand不明确的 SNP，默认 = FALSE。
ambig=FALSE

# direct.filter	指示是否应删除半数以上贡献队列中信息缺失的 SNP（如方向列中的缺失信息所示），默认 = FALSE。
direct.filter=FALSE

#运行Prepare the summary statistics for GWAS
# 1.sumstats 函数要求变量的排列顺序与上述 ldsc 函数相同。
# 2.sumstats 函数需要 4 个必要参数，以及 6 个默认为 NULL 的附加参数：
compiler::loadcmp("/path/to/panvas_analysis/CKM/sumstats.Rc", env = environment())

####更新后的步骤14
library(Fast2TWAS)
################################################################################
# Combine the summary statistics and LDSC output and run the common factor GWAS#
################################################################################
#导入LDSCoutput和PSYCH_sumstats数据
load("/path/to/panvas_analysis/Genomic_SEM/LDSCoutput.RData")
load("/path/to/panvas_analysis/Genomic_SEM/PSYCH_sumstats.RData")
print(LDSCoutput)
print(PSYCH_sumstats)

#ldsc()函数的输出结果
covstruc=LDSCoutput
estimation="DWLS"
SNPs=PSYCH_sumstats
cores= ifelse( Sys.info()["sysname"]!="Linux",parallel::detectCores(),parallel::detectCores()/2)
toler=1e-20
SNPSE=FALSE
parallel=TRUE
GC = "standard"
MPI=FALSE
TWAS=FALSE
smooth_check=FALSE
printwarn=T
sub = c("G~SNP")
std.lv=TRUE
fix_measurement=TRUE
Q_SNP = FALSE
trait = "userGWAS"
overwrite = TRUE
split_SNP = 10000
split_CHR = FALSE
model <- '

# 一阶因子
F1 =~ NA*AA + AAA + TAAD
F2 =~ NA*CAD + Stroke + PAD + CA + CAVS

# 二阶因子
G =~ NA*F1 + l2*F2
G ~~ 1*G

# 等值约束（保证识别）
F1 ~~ u*F1
F2 ~~ u*F2

# 因子独立
F1 ~~ 0*F2

# 防止 Heywood（关键）
AAA ~~ 1e-6*AAA

# GWAS（正式分析时打开）
G ~ SNP

'
#########注意！注意！注意！注意！注意！慎重！慎重！慎重！慎重！慎重！###########
# 先测试能否运行，避免出错 -----------------------------------------------------
#########注意！注意！注意！注意！注意！慎重！慎重！慎重！慎重！慎重！###########
# 这里写你的代码，先从PSYCH_sumstats中随机挑选10000行测试一下能否成功运行并生成结果
random_numbers <- sample(1:nrow(SNPs), 1000)
compiler::loadcmp("./try_userGWAS.Rc",env = environment())

#########注意！注意！注意！注意！注意！慎重！慎重！慎重！慎重！慎重！###########
#########只有上一步成功运行了之后才能运行这一步，否则会出错，且耗费大量时间#####
####这一步非常消耗时间(根据电脑性能，GWAS数量和cores等不同，时长可达到数天)#####
compiler::loadcmp("./userGWAS.Rc",env = environment())

# 保存结果到文件
save(PSYCH_factor, file = "userGWA.RData")
allresult <- do.call(rbind,PSYCH_factor)
#循环保存结果
for (factor in sub) {
  fa <- trimws(strsplit(factor, "~")[[1]])
  #筛选
  allresult_sub <- subset(allresult,lhs==fa[1] & rhs==fa[2])
  table(allresult_sub$CHR)
  dim(allresult_sub)
  
  ################################################################################
  # Calculating sample size for factors. -----------------------------------------
  ################################################################################
  #restrict to MAF >= 10%
  allresult_sub_limit_MAF<-subset(allresult_sub, allresult_sub$MAF >= 0.1)
  #calculate expected sample size (N_hat)
  N_hat<-mean(1/((2*allresult_sub_limit_MAF$MAF*(1-allresult_sub_limit_MAF$MAF))*allresult_sub_limit_MAF$SE^2),na.rm = TRUE)
  
  #显示N_hat
  print(N_hat)
  
  #添加expected sample size (N_hat)列
  allresult_sub$N_hat <- N_hat
  
  print(head(allresult_sub))
  
  #保存数据
  save(allresult_sub,N_hat,
       file=paste0(trait,"_",factor,"_userGWAS_N_hat.RData"))
  data.table::fwrite(allresult_sub,
                     file = paste0(trait,"_",factor,"_userGWAS_N_hat.txt"),
                     sep = "\t")
}
