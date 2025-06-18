library(data.table)
library(ggplot2)
library(ggsci)
library(dplyr)
library(multcompView)
calpercent = function(x,y){
  x = round(x,2)
  y = round(y,2)
  return(round(((x - y)/y)*100,2))
}

classModel = function(x){
x[mname%flike%"gblup",MODE:="GBLUP"]
if(any(x$mname == "pblup")) x[mname%flike%"pblup",MODE:="PBLUP"]
x[mname%flike%"cnn",MODE:="CNN"]
x[mname%flike%"fc",MODE:="FCN"]
x[mname%flike%"mul",MODE:="MNNDR"]

x[mname%flike%"linear",type:="MT_continuous"]
x[mname%flike%"mt",type:="MT_binary"]
x[mname%flike%"blup_mt",type:="MT_continuous"]
x[is.na(x$type),type:="ST"]

x[mname%flike%"_a",combine:="DL_gblup"]
x[!mname%flike%"_a",combine:="original"]
return(x)
}

checking = function(wac_mean_sd,pac_mean_sd,cal){

wdt = wac_mean_sd[combine == "original",]
pdt = pac_mean_sd[combine == "original",]

# if(cal == "ACC") model = c("PBLUP","GBLUP","FCN","CNN","MNN") else model = c("GBLUP","FCN","CNN","MNN")
model = c("PBLUP","GBLUP","FCN","CNN","MNNDR")
type = c("ST","MT_continuous","MT_binary")
if(cal == "COV") value_type = c("accov","fcov","wcov","fwcov","wfcov")
if(cal == "VAR") value_type = c("va","fva","wva","fwcova")
if(cal == "ACC") value_type = c("ac","famacc","wfac")
output = data.table()
for(v in value_type){
  for(m in model){
    for(t in type){

      if(v %in% c("wcov","wva","fwcova","wfac") & m == "PBLUP") next

      if(v %in% c("accov","ac","va")) temp_dt = pdt else temp_dt = wdt

      #向量化实现，可以节省关于model和type的两个for循环
      # if(v %in% c("wcov","wva","fwcova","wfac")) temp_dt = temp_dt[MODE!="PBLUP",]

      # vsd = paste(v,"_sd",sep = "")

      # temp_value = temp_dt[,..v]

      # temp_value_sd = temp_dt[,..vsd]

      # temp_output = data.table(Vtype = v,MODE = temp_dt$MODE,type = temp_dt$type,value = temp_value, value_sd = temp_value_sd)

      #colnames(temp_output) = c("Vtype","MODE","type","value","value_sd")

      # output = rbind(output, temp_output)

      vsd = paste(v,"_sd",sep = "")

      temp_value = temp_dt[MODE ==m & type==t,..v]

      temp_value_sd = temp_dt[MODE ==m & type==t,..vsd]

      temp_output = data.table(Vtype = v,MODE = m,type = t,value = temp_value[[1]], value_sd = temp_value_sd[[1]])

      output = rbind(output,temp_output)

    }
  }
}
return(output)
}


checking_rep = function(wac_mean_sd,pac_mean_sd,cal){

wdt = wac_mean_sd[combine == "original",]
pdt = pac_mean_sd[combine == "original",]

# if(cal == "ACC") model = c("PBLUP","GBLUP","FCN","CNN","MNN") else model = c("GBLUP","FCN","CNN","MNN")
model = c("PBLUP","GBLUP","FCN","CNN","MNNDR")
type = c("ST","MT_continuous","MT_binary")
if(cal == "COV") value_type = c("accov","fcov","wcov","fwcov","wfcov")
if(cal == "VAR") value_type = c("va","fva","wva","fwcova")
if(cal == "ACC") value_type = c("ac","famacc","wfac")
output = data.table()
for(v in value_type){
  for(m in model){
    for(t in type){

      if(v %in% c("wcov","wva","fwcova","wfac") & m == "PBLUP") next

      if(v %in% c("accov","ac","va")) temp_dt = pdt else temp_dt = wdt

      temp_value = temp_dt[MODE ==m & type==t,..v]

      temp_output = data.table(Vtype = v,MODE = m,type = t,value = temp_value[[1]])

      output = rbind(output,temp_output)

    }
  }
}
return(output)
}

theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 12),#坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          #legend.title=element_blank(),
          #legend.position=c(0.5, 0.95),#图例在绘图区域的位置
          #legend.position="none",
          legend.position="right",
          #legend.direction = "horizontal",
          legend.direction = "vertical",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}


pac = fread("pac.csv",sep = ",")
pac_mean_sd = fread("pac_mean_sd.csv",sep = ",")
setorder(pac_mean_sd,ac)
pac_mean_sd
pac_mean_sd[!mname%flike%"_a",]
t.test(x = pac[mname =="mul_linear_a",ac],y = pac[mname=="gblup_mt",ac])

#检测模型之间的差异
anova_model<- aov(ac~MODE,data = pac[combine == "original"&type == "ST",] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]


anova_model<- aov(ac~MODE,data = pac[combine == "original"&type == "MT_continuous",] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]
a = 0.74
b=0.7
calpercent(x=a,y=b)


ckck = rbind(pac[combine == "original"&type == "MT_binary",],pac[combine == "original"&MODE=="GBLUP"& type== "MT_continuous",])
anova_model<- aov(ac~MODE,data = ckck )

anova_model<- aov(ac~MODE,data = pac[combine == "original"&type == "MT_binary",] )

summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]


#检测单双性状之间的差异
Me = "CNN"
anova_model<- aov(ac~type,data = pac[combine == "original"&MODE==Me,] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"type",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$type, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]
a = 0.73
b=0.74
calpercent(x=a,y=b)


Me = "GBLUP"
t.test(x = pac[combine == "original"&MODE ==Me&type == "ST",ac],
y = pac[combine == "original"&MODE ==Me&type == "MT_continuous",ac])

Me = "PBLUP"
t.test(x = pac[combine == "original"&MODE ==Me&type == "ST",ac],
y = pac[combine == "original"&MODE ==Me&type == "MT_continuous",ac])
# pac_mean_sd[mname%flike%"_a"&MODE =="CNN",MODE:="CNN_a"]
# pac_mean_sd[mname%flike%"_a"&MODE =="FCN",MODE:="FCN_a"]
# pac_mean_sd[mname%flike%"_a"&MODE =="MNN",MODE:="MNN_a"]
# zt = pac_mean_sd
# zt$value = zt[,..pv]
# zt$valuese = zt[,..pvse]
# zt$type = factor(zt$type,level = c("ST","MT_continuous","MT_binary"))
# zt$MODE = factor(zt$MODE,level = c("pblup","gblup","fc","cnn","mul","fc_a","cnn_a","mul_a"))
# zt$MODE = factor(zt$MODE,level = c("pblup","gblup","fc","fc_a","cnn","cnn_a","mul","mul_a"))
# P <- ggplot(data = zt, aes(x = MODE, y = value, group = type, fill = type)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
#   geom_text(aes(label = round(value,2)), position = position_dodge(width = 0.9),vjust = -1, size = 4) +  # 添加数值标签
#   xlab("Model") +
#   ylab("Within-family accuracy") +
#   theme_zg() +
#   theme(legend.position ="right", legend.title = element_blank()) +
#   geom_errorbar(aes(ymin = value - valuese, ymax = value + valuese),
#                 position = position_dodge(width = 0.9),  # 与条形图对齐
#                 width = 0.2, alpha = 0.5) +
#   labs(fill = "Mat") +
#   scale_fill_aaas()

# ggsave("Figrue_pac_medium_complete.pdf", P , width = 15, height = 5, dpi = 300)

wac = fread("wac.csv",sep = ",")
wac_mean_sd = fread("wac_mean_sd.csv",sep = ",")
setorder(wac_mean_sd,wac)
wac_mean_sd
t.test(x = wac[mname =="gblup_mt",wac],y = wac[mname=="mul_mt",wac])

wac_mean_sd = classModel(wac_mean_sd)

pv = "famacc"
pvse = "famacc_sd"
zt = wac_mean_sd[combine == "original",]
zt$value = zt[,..pv]
zt$valuese = zt[,..pvse]
zt$type = factor(zt$type,level = c("ST","MT_continuous","MT_binary"))
zt$MODE = factor(zt$MODE,level = c("PBLUP","GBLUP","FCN","CNN","MNNDR"))

pv = "wfac"
pvse = "wfac_sd"
zt = wac_mean_sd[combine == "original"&MODE!="PBLUP",]
zt$value = zt[,..pv]
zt$valuese = zt[,..pvse]
zt$type = factor(zt$type,level = c("ST","MT_continuous","MT_binary"))
zt$MODE = factor(zt$MODE,level = c("GBLUP","FCN","CNN","MNNDR"))

P <- ggplot(data = zt, aes(x = MODE, y = value, group = type, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
  geom_text(aes(label = round(value,2)), position = position_dodge(width = 0.9),vjust = -2, size = 4) +  # 添加数值标签
  xlab("Model") +
  ylab("Within-family accuracy") +
  theme_zg() +
  theme(legend.position ="right", legend.title = element_blank()) +
  geom_errorbar(aes(ymin = value - valuese, ymax = value + valuese),
                position = position_dodge(width = 0.9),  # 与条形图对齐
                width = 0.2, alpha = 0.5) +
  labs(fill = "Mat") +
  scale_fill_aaas()

# ggsave("Figrue_wac_pure.pdf", P , width = 15, height = 6, dpi = 300)
# ggsave("Figrue_wac_low.pdf", P , width = 15, height = 6, dpi = 300)
# ggsave("Figrue_wac_medium.pdf", P , width = 15, height = 6, dpi = 300)
# ggsave("Figrue_wac_high.pdf", P , width = 15, height = 6, dpi = 300)


ggsave("Figrue_wfac_pure.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_wfac_low.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_wfac_medium.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_wfac_high.pdf", P , width = 15, height = 6, dpi = 300)

ggsave("Figrue_fac_pure.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_fac_low.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_fac_medium.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_fac_high.pdf", P , width = 15, height = 8, dpi = 300)



wac = classModel(wac)
#检测模型之间的差异
anova_model<- aov(wfac~MODE,data = wac[combine == "original"&type == "ST"&MODE!="PBLUP",] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]

anova_model<- aov(wfac~MODE,data = wac[combine == "original"&type == "MT_continuous"&MODE!="PBLUP",] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]
a = 0.33
b=0.31
calpercent(x=a,y=b)


ckck = rbind(wac[combine == "original"&type == "MT_binary",],wac[combine == "original"&MODE=="GBLUP"& type== "MT_continuous",])
anova_model<- aov(wfac~MODE,data = ckck )

anova_model<- aov(wfac~MODE,data = wac[combine == "original"&type == "MT_binary",] )


summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]


#检测单双性状之间的差异
Me = "CNN"
anova_model<- aov(wfac~type,data = wac[combine == "original"&MODE==Me,] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"type",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$type, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]
a = 0.2
b=0.31
calpercent(x=a,y=b)


Me = "GBLUP"
t.test(x = wac[combine == "original"&MODE ==Me&type == "ST",wfac],
y = wac[combine == "original"&MODE ==Me&type == "MT_continuous",wfac])

Me = "PBLUP"
t.test(x = wac[combine == "original"&MODE ==Me&type == "ST",wfac],
y = wac[combine == "original"&MODE ==Me&type == "MT_continuous",wfac])
# wac_mean_sd[mname%flike%"_a"&MODE =="cnn",MODE:="cnn_a"]
# wac_mean_sd[mname%flike%"_a"&MODE =="fc",MODE:="fc_a"]
# wac_mean_sd[mname%flike%"_a"&MODE =="mul",MODE:="mul_a"]
# zt = wac_mean_sd
# zt$value = zt[,..pv]
# zt$valuese = zt[,..pvse]
# zt$type = factor(zt$type,level = c("ST","MT_linear","MT_binary"))
# zt$MODE = factor(zt$MODE,level = c("gblup","fc","cnn","mul","fc_a","cnn_a","mul_a"))
# zt$MODE = factor(zt$MODE,level = c("gblup","fc","fc_a","cnn","cnn_a","mul","mul_a"))
# P <- ggplot(data = zt, aes(x = MODE, y = value, group = type, fill = type)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
#   geom_text(aes(label = round(value,2)), position = position_dodge(width = 0.9),vjust = -2, size = 4) +  # 添加数值标签
#   xlab("Model") +
#   ylab("Within-family accuracy") +
#   theme_zg() +
#   theme(legend.position ="right", legend.title = element_blank()) +
#   geom_errorbar(aes(ymin = value - valuese, ymax = value + valuese),
#                 position = position_dodge(width = 0.9),  # 与条形图对齐
#                 width = 0.2, alpha = 0.5) +
#   labs(fill = "Mat") +
#   scale_fill_aaas()
  
# ggsave("Figrue_wac_medium_complete.pdf", P , width = 15, height = 5, dpi = 300)




#检测模型之间的差异
anova_model<- aov(famacc~MODE,data = wac[combine == "original"&type == "ST",] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]

anova_model<- aov(famacc~MODE,data = wac[combine == "original"&type == "MT_continuous",] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]
a = 0.89
b=0.86
calpercent(x=a,y=b)


ckck = rbind(wac[combine == "original"&type == "MT_binary",],wac[combine == "original"&MODE%in%c("GBLUP","PBLUP")& type== "MT_continuous",])
anova_model<- aov(famacc~MODE,data = ckck )

anova_model<- aov(famacc~MODE,data = wac[combine == "original"&type == "MT_binary",] )

summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]


#检测单双性状之间的差异
Me = "MNN"
anova_model<- aov(famacc~type,data = wac[combine == "original"&MODE==Me,] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"type",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$type, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]
a = 0.87
b=0.89
calpercent(x=a,y=b)


Me = "GBLUP"
t.test(x = wac[combine == "original"&MODE ==Me&type == "ST",famacc],
y = wac[combine == "original"&MODE ==Me&type == "MT_continuous",famacc])

Me = "PBLUP"
t.test(x = wac[combine == "original"&MODE ==Me&type == "ST",famacc],
y = wac[combine == "original"&MODE ==Me&type == "MT_continuous",famacc])

#=============================================================================
pac = fread("pac.csv",sep = ",")
wac = fread("wac.csv",sep = ",")
pac = classModel(pac)
wac = classModel(wac)

pac_mean_sd = fread("pac_mean_sd.csv",sep = ",")
wac_mean_sd = fread("wac_mean_sd.csv",sep = ",")
pac_mean_sd <- classModel(pac_mean_sd)
wac_mean_sd = classModel(wac_mean_sd)

#checking
zt = checking(wac_mean_sd,pac_mean_sd,cal = "ACC")
zt1 = zt
zt1$MODE = factor(zt1$MODE,level = c("PBLUP","GBLUP","FCN","CNN","MNNDR"))
zt1$type = factor(zt1$type,level = c("ST","MT_continuous","MT_binary"))
zt1[Vtype == "ac",Vtype:="Individual"]
zt1[Vtype == "famacc",Vtype:="Family"]
zt1[Vtype == "wfac",Vtype:="Within-family"]
zt1$Vtype = factor(zt1$Vtype,level = c("Individual","Family","Within-family"))
P <- ggplot(data = zt1, aes(x = MODE, y = value, group = Vtype, fill = Vtype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
  geom_text(aes(label = round(value,2)), position = position_dodge(width = 0.9),vjust = -2, size = 4) +  # 添加数值标签
  xlab("Model") +
  ylab("Accuracy") +
  theme_zg() +
  theme(legend.position ="right", legend.title = element_blank()) +
  geom_errorbar(aes(ymin = value - value_sd, ymax = value + value_sd),
                position = position_dodge(width = 0.9),  # 与条形图对齐
                width = 0.2, alpha = 0.5) +
  labs(fill = "Mat") +
  scale_fill_aaas()+facet_wrap(~type, scales = "free_x")
ggsave("Figrue_acc_pure.pdf", P , width = 17, height = 8, dpi = 300)
ggsave("Figrue_acc_low.pdf", P , width = 17, height = 8, dpi = 300)
ggsave("Figrue_acc_medium.pdf", P , width = 17, height = 8, dpi = 300)
ggsave("Figrue_acc_high.pdf", P , width = 17, height = 8, dpi = 300)





zt = checking(wac_mean_sd,pac_mean_sd,cal = "VAR")
zt1 = zt[Vtype!="fwcova",]
zt1$MODE = factor(zt1$MODE,level = c("PBLUP","GBLUP","FCN","CNN","MNNDR"))
zt1$type = factor(zt1$type,level = c("ST","MT_continuous","MT_binary"))
zt1[Vtype == "va",Vtype:="Individual"]
zt1[Vtype == "fva",Vtype:="Family"]
zt1[Vtype == "wva",Vtype:="Within-family"]
zt1$Vtype = factor(zt1$Vtype,level = c("Individual","Family","Within-family"))
P <- ggplot(data = zt1, aes(x = MODE, y = value, group = Vtype, fill = Vtype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
  geom_text(aes(label = round(value,2)), position = position_dodge(width = 0.9),vjust = -2, size = 4) +  # 添加数值标签
  xlab("Model") +
  ylab("Variance") +
  theme_zg() +
  theme(legend.position ="right", legend.title = element_blank()) +
  geom_errorbar(aes(ymin = value - value_sd, ymax = value + value_sd),
                position = position_dodge(width = 0.9),  # 与条形图对齐
                width = 0.2, alpha = 0.5) +
  labs(fill = "Mat") +
  scale_fill_aaas()+facet_wrap(~type, scales = "free_x")
ggsave("Figrue_var_pure.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_var_low.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_var_medium.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_var_high.pdf", P , width = 15, height = 6, dpi = 300)




zt = checking(wac_mean_sd,pac_mean_sd,cal = "COV")
zt1 = zt[!Vtype%in%c("fwcov","wfcov"),]
zt1$MODE = factor(zt1$MODE,level = c("PBLUP","GBLUP","FCN","CNN","MNNDR"))
zt1$type = factor(zt1$type,level = c("ST","MT_continuous","MT_binary"))
zt1[Vtype == "accov",Vtype:="Individual"]
zt1[Vtype == "fcov",Vtype:="Family"]
zt1[Vtype == "wcov",Vtype:="Within-family"]
zt1$Vtype = factor(zt1$Vtype,level = c("Individual","Family","Within-family"))
P <- ggplot(data = zt1, aes(x = MODE, y = value, group = Vtype, fill = Vtype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
  geom_text(aes(label = round(value,2)), position = position_dodge(width = 0.9),vjust = -2, size = 4) +  # 添加数值标签
  xlab("Model") +
  ylab("Variance") +
  theme_zg() +
  theme(legend.position ="right", legend.title = element_blank()) +
  geom_errorbar(aes(ymin = value - value_sd, ymax = value + value_sd),
                position = position_dodge(width = 0.9),  # 与条形图对齐
                width = 0.2, alpha = 0.5) +
  labs(fill = "Mat") +
  scale_fill_aaas()+facet_wrap(~type, scales = "free_x")

ggsave("Figrue_cov_pure.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_cov_low.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_cov_medium.pdf", P , width = 15, height = 6, dpi = 300)
ggsave("Figrue_cov_high.pdf", P , width = 15, height = 8, dpi = 300)




pac_mean_sd = fread("pac_mean_sd.csv",sep = ",")
wac_mean_sd = fread("wac_mean_sd.csv",sep = ",")
pac_mean_sd <- classModel(pac_mean_sd)
wac_mean_sd = classModel(wac_mean_sd)

pac = fread("pac.csv",sep = ",")
wac = fread("wac.csv",sep = ",")
pac = classModel(pac)
wac = classModel(wac)

value_type = c("ac","famacc","wfac")
type = c("ST","MT_continuous","MT_binary")
output_test = data.table()

for(t in type){

for(v in value_type){

  if(v == "ac") temp_dt = pac[type == t,] else temp_dt = wac[type == t,]

  if(v == "wfac") temp_dt = temp_dt[MODE!="PBLUP",]
  
  aov_res <- aov(as.formula(paste(v, "~ MODE")), data = temp_dt)
  
  tukey_res <- TukeyHSD(aov_res,conf.level = 0.95)
  
  letters <- multcompLetters4(aov_res, tukey_res)
  
  output_test = rbind(output_test,data.table(MODE = names(letters$MODE$Letters),Letter = letters$MODE$Letters,type = t,Vtype = v))
   
}

}

# df = checking_rep(wac,pac,cal="ACC")

# output_test <- data.table()
# for(t in unique(df$type)){
# for (g in unique(df$Vtype)) {

#   sub <- df %>% filter(Vtype == g, type == t)
  
#   # Run ANOVA + Tukey
#   aov_res <- aov(value ~ MODE, data = sub)
#   tukey_res <- TukeyHSD(aov_res)
#   letters <- multcompLetters4(aov_res, tukey_res)
  
#   letter_df <- data.frame(
#     MODE = names(letters$MODE$Letters),
#     Letter = letters$MODE$Letters,
#     Vtype = g,
#     type = t
#   )
  
#   output_test <- bind_rows(output_test, letter_df)
# }
# }
zt = checking(wac_mean_sd,pac_mean_sd,cal = "ACC")
zt1 <- merge(output_test, zt, by = c("MODE", "type","Vtype"), all = FALSE)
zt1$MODE = factor(zt1$MODE,level = c("PBLUP","GBLUP","FCN","CNN","MNNDR"))
zt1$type = factor(zt1$type,level = c("ST","MT_continuous","MT_binary"))
zt1[Vtype == "ac",Vtype:="Individual"]
zt1[Vtype == "famacc",Vtype:="Family"]
zt1[Vtype == "wfac",Vtype:="Within-family"]
zt1$Vtype = factor(zt1$Vtype,level = c("Individual","Family","Within-family"))
zt1[, letter_y_value := value + 1.1 * value_sd, by = Vtype]
zt1[, letter_y_test := value + 3.1 * value_sd, by = Vtype]
P <- ggplot(data = zt1[type == "ST",], aes(x = MODE, y = value, group = MODE, fill = MODE)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
  geom_text(aes(label = round(value,2),y = letter_y_value), position = position_dodge(width = 0.9),vjust = 0, size = 4) +  # 添加数值标签
  xlab("Model") +
  ylab("Accuracy") +
  theme_zg() +
  theme(legend.position ="right", legend.title = element_blank()) +
  geom_errorbar(aes(ymin = value - value_sd, ymax = value + value_sd),
                position = position_dodge(width = 0.9),  # 与条形图对齐
                width = 0.2, alpha = 0.5) +
  labs(fill = "Mat") +
  scale_fill_aaas()+geom_text(aes(label = Letter, y = letter_y_test), position = position_dodge(width = 0.9),vjust = 0) + facet_wrap(~type + Vtype, scales = "free")+
  ggtitle(label = "ST and MT prediction across models",subtitle = "Simulated dataset (high non-additive)") +
  theme(
    plot.title = element_text(
      size = 20,          # 大字号
      face = "bold",      # 加粗
      hjust = 0.5,        # 居中
      color = "black"     # 颜色（可选）
    ),
    # 调整x轴标签字体
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    # 调整y轴标签字体
    axis.title.y = element_text(size = 16, face = "bold", color = "black"),
    #调整副标题
    plot.subtitle = element_text(size = 12, hjust = 0.5, face = "italic")
  )
ggsave("Figrue_ac_high.pdf", P , width = 15, height = 6, dpi = 300)
