####Genomic SEM
library(devtools)
library(GenomicSEM)
library(dplyr)
library(ggplot2)
library(data.table)
rm(list = ls())
setwd("/path/to/panvas_analysis/Genomic_SEM")

####CAVS清洗
trait="CAVS"
tes <-fread("calcific aortic valve stenosis_38.gz", select = c("rs_id","chromosome", "base_pair_location","effect_allele", "other_allele", "beta", "standard_error", "p_value","effect_allele_frequency","n")) 
colnames(tes)= c("SNP","CHR","BP","effect_allele", "other_allele", "beta", "se", "pval","eaf","N")
options(timeout = 600)
GWAS_hg19=gwas_liftover(gwasfile=tes,
                        build_from = "hg38",
                        build_to = "hg19",
                        savefile = paste0(trait,"hg19.txt"))
head(GWAS_hg19)

data=gwas2format(gwas_dat  = GWAS_hg19,
                 chr_col = "CHR",                      #GWAS数据中SNPs的染色体列名，默认CHR
                 pos_col = "BP",                      #GWAS数据中SNPs的基因位置列名，默认BP
                 snp_col = "SNP",                      #GWAS数据中SNPs的rsid列名，默认SNP
                 effect_allele_col = "effect_allele",  #GWAS数据中SNPs的效应等位基因列名，默认effect_allele
                 other_allele_col = "other_allele",    #GWAS数据中SNPs的非效应等位基因列名，默认other_allele 
                 eaf_col = "eaf",                      #GWAS数据中SNPs的效应等位基因频率列名，默认eaf
                 beta_col = "beta",                    #GWAS数据中SNPs的效应值列名，默认beta
                 se_col = "se",                        #GWAS数据中SNPs的标准误差列名，默认se
                 pval_col = "pval",                    #GWAS数据中SNPs的pval值列名，默认pval
                 log_pval = FALSE,                     #GWAS数据中SNPs的pval值是否为-log10(pval)，默认FALSE
                 ncase_col = "ncase",                  #GWAS数据中SNPs的病例数列名，默认ncase
                 ncontrol_col = "ncontrol",            #GWAS数据中SNPs的对照数列名，默认ncontrol
                 samplesize_col = "N",        #GWAS数据中SNPs的总样本数列名，默认N
                 min_pval = 1e-200,                    #强制设置pval的最小值，默认为1e-200
                 samplesize = NULL)                    #强制设置GWAS数据中SNPs的总样本数，默认NULL
#查看数据
head(data)
#保存数据
fwrite(data, file = paste0(trait,".txt"),sep = "\t")

####TAA清洗
trait="TAA"
tes <- fread("2021_Roychowdhury_et_al_TAA.txt.gz")
tes$EAF <- ifelse(tes$MAF <= 0.5,tes$MAF,1 - tes$MAF)
head(tes)
data=gwas2format(gwas_dat  = tes,
                 chr_col = "CHROM",                      #GWAS数据中SNPs的染色体列名，默认CHR
                 pos_col = "POS",                      #GWAS数据中SNPs的基因位置列名，默认BP
                 snp_col = "SNPID",                      #GWAS数据中SNPs的rsid列名，默认SNP
                 effect_allele_col = "ALT",  #GWAS数据中SNPs的效应等位基因列名，默认effect_allele
                 other_allele_col = "REF",    #GWAS数据中SNPs的非效应等位基因列名，默认other_allele 
                 eaf_col = "EAF",                      #GWAS数据中SNPs的效应等位基因频率列名，默认eaf
                 beta_col = "BETA",                    #GWAS数据中SNPs的效应值列名，默认beta
                 se_col = "SE",                        #GWAS数据中SNPs的标准误差列名，默认se
                 pval_col = "P",                    #GWAS数据中SNPs的pval值列名，默认pval
                 log_pval = FALSE,                     #GWAS数据中SNPs的pval值是否为-log10(pval)，默认FALSE
                 ncase_col = "1351",                  #GWAS数据中SNPs的病例数列名，默认ncase
                 ncontrol_col = "18295 ",            #GWAS数据中SNPs的对照数列名，默认ncontrol
                 samplesize_col = "N",        #GWAS数据中SNPs的总样本数列名，默认N
                 min_pval = 1e-200,                    #强制设置pval的最小值，默认为1e-200
                 samplesize = NULL)                    #强制设置GWAS数据中SNPs的总样本数，默认NULL
#查看数据
head(data)
#保存数据
fwrite(data, file = paste0(trait,".txt"),sep = "\t")


####PAD清洗
trait="PAD"
tes <-fread("pad_primary.txt.gz", select = c("chr", "bp","or_allele","nonor_allele","or", "or_se", "p","eaf","ncases","ncontrols")) #注意这里我设置的是制表符分隔符，如果你的文件不是制表符的话，需要修改成对应的分隔符
colnames(tes)= c("chromosome","pos","Allele1", "Allele2", "OR", "StdErr", "P-value","eaf","Ncases","Ncontrols")
tes$N=tes$Ncases+tes$Ncontrols
tes$Beta=log(tes$OR)
head(tes)
colnames(tes)= c("chromosome","pos","Allele1", "Allele2", "OR", "StdErr", "P-value","eaf","Ncases","Ncontrols","N","Effect")
tes$chromosome.start <- paste(tes$chromosome, tes$pos, sep = ":")
head(tes)
match = read.table("/path/to/panvas_analysis/CVD/snp150_hg19.txt.gz",header=T,check.names=F,sep="\t") 
colnames(match) <- c("chromosome.start","rsid")
need=dplyr::left_join(tes,match,by="chromosome.start") #如果snp150_hg19.txt文件中有对应的RS号，则比对到test.txt文件中，如果没有的话，就变为NA
need[need == "NA"] <- NA
need$rsid[need$rsid == ""] <- NA
test <- na.omit(need)
write.table(test, file = "PAD_anno.txt",sep = "\t", row.names = FALSE, quote = FALSE)
test$N_total=test$Ncases+test$Ncontrols

##对PAD尝试新的注释方法
# ========= 1) 安装与加载 =========
#GitHub安装
#先设置一下令牌
gitcreds::gitcreds_set() #输入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
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("nicksunderland/genepi.utils")
library(genepi.utils)
# ========= 2) 清洗字段 =========
library(data.table)
library(genepi.utils)
trait="PAD"
tes <-fread("pad_primary.txt.gz", select = c("chr", "bp","or_allele","nonor_allele","or", "or_se", "p","eaf","ncases","ncontrols")) #注意这里我设置的是制表符分隔符，如果你的文件不是制表符的话，需要修改成对应的分隔符
tes$N=tes$ncases+tes$ncontrols
tes$Beta=log(tes$OR)
head(tes)
setDT(res)
# -------------------------
# 0) 统一字段 & hg19(=b37)
# -------------------------
# 你的列：chr, bp, or_allele, nonor_allele, snp
res[, chr := sub("^chr", "", as.character(chr), ignore.case = TRUE)]
res[, bp  := as.integer(bp)]
res[, `:=`(
  or_allele    = toupper(as.character(or_allele)),
  nonor_allele = toupper(as.character(nonor_allele))
)]

# -------------------------
# 1) 从 snp 列提取 I/D，并做一致性检查（可选但强烈推荐）
#    snp 形如 chr10:100024611:I
# -------------------------
res[, snp_type := fifelse(grepl(":(I|D)$", snp), sub(".*:(I|D)$", "\\1", snp), NA_character_)]

# 用等位基因长度差推断 I/D（对 indel 很有用）
# 插入：alt 比 ref 长；缺失：alt 比 ref 短
res[, inferred_type := {
  lo <- nchar(or_allele)
  ln <- nchar(nonor_allele)
  fifelse(lo > ln, "I",
          fifelse(lo < ln, "D", "S"))  # S=SNV 或长度相同的替换
}]

# 标记“snp列类型”和“长度推断类型”不一致的（这类最容易对不上dbSNP）
res[, type_mismatch := !is.na(snp_type) & snp_type %chin% c("I","D") & inferred_type %chin% c("I","D") & snp_type != inferred_type]

# -------------------------
# 2) 基础清洗：去掉明显非法等位基因（否则映射会失败/误匹配）
# -------------------------
valid_allele <- function(x) grepl("^[ACGT]+$", x)

res[, bad_allele := !(valid_allele(or_allele) & valid_allele(nonor_allele))]
# 你可以选择直接过滤掉 bad_allele 的行（推荐）
res_clean <- res[!is.na(chr) & !is.na(bp) & !bad_allele]

# -------------------------
# 3) 准备 genepi.utils 输入：把列名对齐为 chr/pos/EA/NEA
#    注意：genepi.utils 的 ea_col/nea_col 只是“用于匹配”的两列，不等于必须是效应/非效应，
#    它会在 flip="report" 时尝试对调匹配并报告 flip。
# -------------------------
dt_map <- res_clean[, .(
  chromosome = chr,
  pos= bp,
  EA= or_allele,
  NEA= nonor_allele,
  snp,
  snp_type,
  inferred_type,
  type_mismatch
)]

# -------------------------
# 4) dbSNP 目录：必须存在，且里面要有 vcf.gz + .tbi
# -------------------------
dbsnp_dir <- "/path/to/panvas_analysis/Genomic_SEM/dbsnp/hg19/"
dbsnp_dir <- path.expand(dbsnp_dir)
stopifnot(dir.exists(dbsnp_dir))

# （可选）简单检查目录下是否有 vcf.gz/tbi
vcf_gz <- list.files(dbsnp_dir, pattern = "\\.vcf\\.gz$", full.names = TRUE)
tbi    <- list.files(dbsnp_dir, pattern = "\\.tbi$", full.names = TRUE)
if (length(vcf_gz) == 0) warning("dbsnp_dir 里没发现 .vcf.gz 文件：", dbsnp_dir)
if (length(tbi) == 0)    warning("dbsnp_dir 里没发现 .tbi 文件：", dbsnp_dir)

# -------------------------
# 5) 分批映射（按染色体/或按块），更稳更省内存
# -------------------------
setkey(dt_map, chromosome, pos)

by_chr <- split(dt_map, dt_map$chromosome)

anno_list <- lapply(by_chr, function(x) {
  # 避免小批次中出现奇怪类型
  x[, chromosome := as.character(chromosome)]
  x[, pos := as.integer(pos)]
  
  genepi.utils::chrpos_to_rsid(
    dt        = x,
    chr_col   = "chromosome",
    pos_col   = "pos",
    ea_col    = "EA",
    nea_col   = "NEA",
    build     = "b37",          # hg19/GRCh37
    dbsnp_dir = dbsnp_dir,      # 关键：你原代码漏了这个
    flip      = "report",
    verbose   = FALSE
  )
})

res_anno <- rbindlist(anno_list, use.names = TRUE, fill = TRUE)
setDT(res_anno)

# -------------------------
# 6) 合并回原表（保持原字段不丢）
# -------------------------
# res_anno 里应该新增 rsid 和 rsid_flip_match 等列（取决于你安装的 genepi.utils 版本）
# 为了稳妥，只取需要的列再 merge
keep_cols <- intersect(names(res_anno), c("chromosome","pos","EA","NEA","rsid","rsid_flip_match"))
map_out <- unique(res_anno[, ..keep_cols])

# 回填到原数据：按 chr/bp + alleles 合并（更严格）；如果你想宽松可只按 chr/bp 合并
setnames(map_out, c("chromosome","pos","EA","NEA"), c("chr","bp","or_allele","nonor_allele"))

res_final <- merge(
  res, map_out,
  by = c("chr","bp","or_allele","nonor_allele"),
  all.x = TRUE
)

# -------------------------
# 7) 结果汇总与质量检查
# -------------------------
setDT(res_final)

res_final[, .(
  n_total = .N,
  n_rsid  = sum(!is.na(rsid)),
  n_flip  = sum(isTRUE(rsid_flip_match), na.rm = TRUE),
  n_type_mismatch = sum(type_mismatch, na.rm = TRUE),
  n_bad_allele = sum(bad_allele, na.rm = TRUE)
)]

# 看看没注释上的，优先关注 indel + type_mismatch
head(res_final[is.na(rsid)][order(type_mismatch, decreasing = TRUE)], 20)
##调整数据格式
data=gwas2format(gwas_dat  = test,
                 chr_col = "chromosome",                      #GWAS数据中SNPs的染色体列名，默认CHR
                 pos_col = "pos",                      #GWAS数据中SNPs的基因位置列名，默认BP
                 snp_col = "rsid",                      #GWAS数据中SNPs的rsid列名，默认SNP
                 effect_allele_col = "Allele1",  #GWAS数据中SNPs的效应等位基因列名，默认effect_allele
                 other_allele_col = "Allele2",    #GWAS数据中SNPs的非效应等位基因列名，默认other_allele 
                 eaf_col = "eaf",                      #GWAS数据中SNPs的效应等位基因频率列名，默认eaf
                 beta_col = "Effect",                    #GWAS数据中SNPs的效应值列名，默认beta
                 se_col = "StdErr",                        #GWAS数据中SNPs的标准误差列名，默认se
                 pval_col = "P-value",                    #GWAS数据中SNPs的pval值列名，默认pval
                 log_pval = FALSE,                     #GWAS数据中SNPs的pval值是否为-log10(pval)，默认FALSE
                 ncase_col = "Ncases",                  #GWAS数据中SNPs的病例数列名，默认ncase
                 ncontrol_col = "Ncontrols",            #GWAS数据中SNPs的对照数列名，默认ncontrol
                 samplesize_col = "N",        #GWAS数据中SNPs的总样本数列名，默认N
                 min_pval = 1e-200,                    #强制设置pval的最小值，默认为1e-200
                 samplesize = NULL)                    #强制设置GWAS数据中SNPs的总样本数，默认NULL
#查看数据
head(data)
#保存数据
fwrite(data,
       file = paste0(trait,".txt"),
       sep = "\t")


####颈动脉硬化CA清洗
trait="CA"
tes <- fread("GCST007435.tsv.gz")
head(tes)
data=gwas2format(gwas_dat  = tes,
                 chr_col = "chromosome",                      #GWAS数据中SNPs的染色体列名，默认CHR
                 pos_col = "base_pair_location",                      #GWAS数据中SNPs的基因位置列名，默认BP
                 snp_col = "rs_id",                      #GWAS数据中SNPs的rsid列名，默认SNP
                 effect_allele_col = "effect_allele",  #GWAS数据中SNPs的效应等位基因列名，默认effect_allele
                 other_allele_col = "other_allele",    #GWAS数据中SNPs的非效应等位基因列名，默认other_allele 
                 eaf_col = "effect_allele_frequency",                      #GWAS数据中SNPs的效应等位基因频率列名，默认eaf
                 beta_col = "beta",                    #GWAS数据中SNPs的效应值列名，默认beta
                 se_col = "standard_error",                        #GWAS数据中SNPs的标准误差列名，默认se
                 pval_col = "p_value",                    #GWAS数据中SNPs的pval值列名，默认pval
                 log_pval = FALSE,                     #GWAS数据中SNPs的pval值是否为-log10(pval)，默认FALSE
                 ncase_col = "21540",                  #GWAS数据中SNPs的病例数列名，默认ncase
                 ncontrol_col = "26894 ",            #GWAS数据中SNPs的对照数列名，默认ncontrol
                 samplesize_col = "n",        #GWAS数据中SNPs的总样本数列名，默认N
                 min_pval = 1e-200,                    #强制设置pval的最小值，默认为1e-200
                 samplesize = NULL)                    #强制设置GWAS数据中SNPs的总样本数，默认NULL
#查看数据
head(data)
#保存数据
fwrite(data, file = paste0(trait,".txt"),sep = "\t")

####AAA清洗
#head(tes)
trait="AAA"
tes <- fread("AAA_MVP.gz")
tes$se=abs(log(tes$odds_ratio)/qnorm(tes$p_value/2))
tes$beta = log(tes$odds_ratio)
tes$N=tes$num_cases+tes$num_controls
head(tes)
library(dplyr)
tes <- as.data.frame(tes)
tes <- tes %>%
  select(rsid, chromosome, base_pair_location, effect_allele, 
         other_allele, beta, se, p_value, effect_allele_frequency, 
         num_cases, num_controls, N)
colnames(tes)= c("SNP","CHR","BP","effect_allele", "other_allele", "beta", "se", "pval","eaf","ncase","ncontrol","N")
tes <- na.omit(tes)
options(timeout = 600)
GWAS_hg19=gwas_liftover(gwasfile=tes,
                        build_from = "hg38",
                        build_to = "hg19",
                        savefile = paste0(trait,"hg19.txt"))
head(GWAS_hg19)

data=gwas2format(gwas_dat  = GWAS_hg19,
                 chr_col = "CHR",                      #GWAS数据中SNPs的染色体列名，默认CHR
                 pos_col = "BP",                      #GWAS数据中SNPs的基因位置列名，默认BP
                 snp_col = "SNP",                      #GWAS数据中SNPs的rsid列名，默认SNP
                 effect_allele_col = "effect_allele",  #GWAS数据中SNPs的效应等位基因列名，默认effect_allele
                 other_allele_col = "other_allele",    #GWAS数据中SNPs的非效应等位基因列名，默认other_allele 
                 eaf_col = "eaf",                      #GWAS数据中SNPs的效应等位基因频率列名，默认eaf
                 beta_col = "beta",                    #GWAS数据中SNPs的效应值列名，默认beta
                 se_col = "se",                        #GWAS数据中SNPs的标准误差列名，默认se
                 pval_col = "pval",                    #GWAS数据中SNPs的pval值列名，默认pval
                 log_pval = FALSE,                     #GWAS数据中SNPs的pval值是否为-log10(pval)，默认FALSE
                 ncase_col = "ncase",                  #GWAS数据中SNPs的病例数列名，默认ncase
                 ncontrol_col = "ncontrol",            #GWAS数据中SNPs的对照数列名，默认ncontrol
                 samplesize_col = "N",        #GWAS数据中SNPs的总样本数列名，默认N
                 min_pval = 1e-200,                    #强制设置pval的最小值，默认为1e-200
                 samplesize = NULL)                    #强制设置GWAS数据中SNPs的总样本数，默认NULL
head(data)
fwrite(data, file = paste0(trait,".txt"),sep = "\t")

####AA清洗
trait="AA"
#finngen数据预处理生成建议的GWAS数据格式
N=464256
dat_fin_hg38=finngenClean(finngenfile = "finngen_R12_I9_aortic dissection.gz",
                          samplesize = N,
                          savefile = paste0(trait,"_QC_hg38.txt") )
head(dat_fin_hg38)

options(timeout = 600)
GWAS_hg19=gwas_liftover(gwasfile=paste0(trait,"_QC_hg38.txt") ,
                        build_from = "hg38",
                        build_to = "hg19",
                        savefile = paste0(trait,"_QC_hg19.txt"))
head(GWAS_hg19)

data=gwas2format(gwas_dat  = GWAS_hg19,
                 chr_col = "CHR",                      #GWAS数据中SNPs的染色体列名，默认CHR
                 pos_col = "BP",                      #GWAS数据中SNPs的基因位置列名，默认BP
                 snp_col = "SNP",                      #GWAS数据中SNPs的rsid列名，默认SNP
                 effect_allele_col = "effect_allele",  #GWAS数据中SNPs的效应等位基因列名，默认effect_allele
                 other_allele_col = "other_allele",    #GWAS数据中SNPs的非效应等位基因列名，默认other_allele 
                 eaf_col = "eaf",                      #GWAS数据中SNPs的效应等位基因频率列名，默认eaf
                 beta_col = "beta",                    #GWAS数据中SNPs的效应值列名，默认beta
                 se_col = "se",                        #GWAS数据中SNPs的标准误差列名，默认se
                 pval_col = "pval",                    #GWAS数据中SNPs的pval值列名，默认pval
                 log_pval = FALSE,                     #GWAS数据中SNPs的pval值是否为-log10(pval)，默认FALSE
                 ncase_col = "ncase",                  #GWAS数据中SNPs的病例数列名，默认ncase
                 ncontrol_col = "ncontrol",            #GWAS数据中SNPs的对照数列名，默认ncontrol
                 samplesize_col = N,        #GWAS数据中SNPs的总样本数列名，默认N
                 min_pval = 1e-200,                    #强制设置pval的最小值，默认为1e-200
                 samplesize = NULL)              
head(data)
data$N=N
fwrite(data, file = paste0(trait,".txt"),sep = "\t")


####index文件制作
#性状名
filename=c("AA.txt",
           "stroke.txt",
           "AAA.txt",
           "CAVS.txt",
           "CAD.txt",
           "TAA.txt",
           "PAD.txt",
           "CA.txt")

#phenotype是GWAS数据的性状名称
phenotype=c("AA",
            "Stroke",
            "AAA",
            "CAVS",
            "CAD",
            "TAAD",
            "PAD",
            "CA")
#cases和controls分别是GWAS数据的病例数和对照数
#若GWAS数据的性状是连续性，则cases和controls设置为NA
cases = c(1150,
          73822,
          16289,
          14819,
          122733,
          1351,
          12086,
          21540
          )
#cases和controls分别是GWAS数据的病例数和对照数
#若GWAS数据的性状是连续性，则cases和controls设置为NA
controls = c(463106,
             1223086,
             427111,
             941863,
             424528,
             18295,
             437462,
             26894)

#samplesize指的是GWAS数据的总样本量
samplesize=c(464256,
             1296908,
             443400,
             956682,
             547261,
             19646,
             449548,
             48434)
####Ntotal=4226135
#ismeta判断GWAS数据是否是多个GWAS数据META后得到的数据
ismeta = c(FALSE,
           TRUE, 
           FALSE, 
           TRUE, 
           TRUE,
           FALSE,
           TRUE,
           TRUE
           )

#population.prev：群体流行率向量：
# 1.这些估计值可从多种来源（如大规模流行病学研究）获得，可能范围 = 0-1。
# 2.同样，如果性状是连续的，其值应等于NA。
population.prev <- c(0.0025, 
                     0.012, 
                     0.0367, 
                     0.01549, 
                     0.02,
                     0.012,
                     0.1,
                     0.2)

input_info_data <- data.frame(filename=filename,
                              phenotype=phenotype,
                              cases=cases,
                              controls=controls,
                              samplesize=samplesize,
                              ismeta=ismeta,
                              population.prev=population.prev)

#保存文件
data.table::fwrite(input_info_data,
                   file = "gwasfile_info.txt",
                   na = "NA",
                   quote = FALSE,
                   sep = "\t")
#读入index文件
setwd("/path/to/panvas_analysis/Genomic_SEM/")
index <- data.table::fread("/path/to/panvas_analysis/Genomic_SEM/gwasfile_info.txt")
for (file in index$filename) {
  #读入数据
  dat <- data.table::fread(paste0("/path/to/panvas_analysis/Genomic_SEM/",file),showProgress = F)
  dat_2 <- dat
  
  # 判断GWAS数据是否是多个GWAS数据META后得到的数据，且GWAS数据类型是病例对照--------
  casecontrol <- ifelse(!is.na(index$case[which(file==index$filename)]) & !is.na(index$controls[which(file==index$filename)]),TRUE,FALSE)
  ismeta <- index$ismeta[which(file==index$filename)]
  if (casecontrol & ismeta ) {
    if (!all(is.na(dat$eaf))) {
      # 若GWAS数据中有eaf则将计算MAF，并修改列名
      dat <- dat %>% na.omit()
      dat$MAF <- ifelse(1-dat$eaf>0.5,dat$eaf,1-dat$eaf)
      dat$eaf <- NULL
      dat_ref <- dat
    }else{
      # 计算有效样本量
      dat$eaf <- NULL
      ref<-data.table::fread("/path/to/panvas_analysis/Genomic_SEM/reference.1000G.maf.0.005.txt.gz",data.table=FALSE)
      dat_ref <-merge(dat,ref,by="SNP",suffixes = c("", ".y"))
      dat_ref <- subset(dat_ref,select = c(colnames(dat_ref)[which(colnames(dat_ref) %in% colnames(Fast2TWAS::show_demo_gwas()))],"MAF"))
    }
    #calculate effective sample size implied by GWAS summary statistics 
    dat_ref$N<-4/((2*dat_ref$MAF*(1-dat_ref$MAF))*dat_ref$se^2)
    
    #calculate total effective N to cap backed out Neff
    Ncases<-index$case[which(file==index$filename)] 
    Ncontrols<-index$controls[which(file==index$filename)] 
    v<-Ncases/(Ncases+Ncontrols)
    
    TotalNeff<-4*v*(1-v)*(Ncases+Ncontrols)
    
    #cap at 1.1 of total effective N
    dat_ref$N<-ifelse(dat_ref$N > 1.1*TotalNeff, 1.1*TotalNeff, dat_ref$N)
    
    #lower limit of 0.5 of total effective N
    dat_ref$N<-ifelse(dat_ref$N < 0.5*TotalNeff, 0.5*TotalNeff, dat_ref$N)
    
    #若GWAS数据中原本的eaf列为空则去掉MAF列
    if (all(is.na(dat_2$eaf))) {
      dat_ref$MAF<-NULL
    }
    
  }else{
    #若GWAS数据不是多个GWAS数据META后得到的数据，和/或GWAS数据类型不是病例对照
    if (!all(is.na(dat$eaf))) {
      # 若GWAS数据中有eaf则将计算MAF，并修改列名
      dat <- dat %>% na.omit()
      dat$MAF <- ifelse(1-dat$eaf>0.5,dat$eaf,1-dat$eaf)
      dat$eaf <- NULL
      dat_ref <- dat
    }else{
      dat_ref <- dat
    }
  }
  
  #保存数据
  data.table::fwrite(dat_ref,
                     file = paste0(index$phenotype[which(file==index$filename)],".txt"),
                     sep = "\t")
  
}

library(Fast2TWAS)
index <- data.table::fread("/path/to/panvas_analysis/Genomic_SEM/gwasfile_info.txt")
setwd("/path/to/panvas_analysis/Genomic_SEM/")
# 开始调用munge函数转化数据
#GWAS数据路径
gwas_files<-paste0("/path/to/panvas_analysis/Genomic_SEM/",index$phenotype,".txt")
hm3<-"/path/to/panvas_analysis/Genomic_SEM/w_hm3.noMHC.snplist"
phen<-index$phenotype
N=rep(NA, length(phen))
info.filter=0.9
maf.filter=0.01
parallel = TRUE
cores = ifelse(length(gwas_files)<parallel::detectCores()/2,length(gwas_files),parallel::detectCores()/2)
compiler::loadcmp("/path/to/panvas_analysis/CKM/munge.Rc", env = environment())
gwas_files<-paste0("/path/to/panvas_analysis/Genomic_SEM/",index$phenotype,".sumstats.gz")

