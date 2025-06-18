library(data.table)
library(ggplot2)
library(ggsci)
library(dplyr)
library(multcompView)
calpercent = function(x,y){
  x = round(x,2)
  y = round(y,2)
  return(round(((x - y)/y)*100,1))
}
a = 0.72
b=0.7
calpercent(x=a,y=b)

classModel = function(x){
x[mname%flike%"gblup",MODE:="GBLUP"]
if(any(x$mname == "pblup")) x[mname%flike%"pblup",MODE:="PBLUP"]
x[mname%flike%"cnn",MODE:="CNN"]
x[mname%flike%"fc",MODE:="FCN"]
x[mname%flike%"mul",MODE:="MNNDR"]

x[mname%flike%"linear",type:="MT"]
# x[mname%flike%"mt",type:="MT_binary"]
x[mname%flike%"blup_mt",type:="MT"]
x[is.na(x$type),type:="ST"]

# x[mname%flike%"_a",combine:="DL_gblup"]
x[!mname%flike%"_a",combine:="original"]
return(x)
}

checking = function(pac_mean_sd,cal){

pdt = pac_mean_sd[combine == "original",]

# if(cal == "ACC") model = c("PBLUP","GBLUP","FCN","CNN","MNN") else model = c("GBLUP","FCN","CNN","MNN")
model = c("GBLUP","FCN","CNN","MNNDR")
type = c("ST","MT")
if(cal == "ACC") value_type = c("ac","auc")
output = data.table()
for(v in value_type){
  for(m in model){
    for(t in type){

      temp_dt = pdt
      
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
# pac_mean_sd[!mname%flike%"_a",]
# t.test(x = pac[mname =="mul_linear_a",ac],y = pac[mname=="gblup_mt",ac])


pac_mean_sd <- classModel(pac_mean_sd)


pac = classModel(pac)

#检测模型之间的差异
anova_model<- aov(ac~MODE,data = pac[combine == "original"&type == "ST",] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]


anova_model<- aov(ac~MODE,data = pac[combine == "original"&type == "MT",] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"MODE",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$MODE, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]
a = 0.35
b=0.27
calpercent(x=a,y=b)


#检测单双性状之间的差异
Me = "MNNDR"
anova_model<- aov(ac~type,data = pac[combine == "original"&MODE==Me,] )
summary(anova_model)
tukey_results <- TukeyHSD(anova_model,"type",conf.level = 0.95)
tukey_results_dt <- as.data.table(tukey_results$type, keep.rownames = "comparison")
colnames(tukey_results_dt)[5] = "pvalue"
tukey_results_dt[pvalue<=0.05,]
a = 0.4
b=0.38
calpercent(x=a,y=b)


Me = "MNNDR"
t.test(x = pac[combine == "original"&MODE ==Me&type == "ST",ac],
y = pac[combine == "original"&MODE ==Me&type == "MT",ac])



#checking
# zt = checking(wac_mean_sd,pac_mean_sd,cal = "ACC")
# zt1 = zt
zt1 = pac_mean_sd
zt1$MODE = factor(zt1$MODE,level = c("PBLUP","GBLUP","FCN","CNN","MNNDR"))
zt1$type = factor(zt1$type,level = c("ST","MT"))
# zt1[Vtype == "ac",Vtype:="Population-wide"]
# zt1[Vtype == "famacc",Vtype:="Between-family"]
# zt1[Vtype == "wfac",Vtype:="Within-family"]
# zt1$Vtype = factor(zt1$Vtype,level = c("Population-wide","Between-family","Within-family"))
P <- ggplot(data = zt1, aes(x = MODE, y = ac, group = type, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
  geom_text(aes(label = round(ac,2)), position = position_dodge(width = 0.9),vjust = -2, size = 4) +  # 添加数值标签
  xlab("Model") +
  ylab("Accuracy") +
  theme_zg() +
  theme(legend.position ="right", legend.title = element_blank()) +
  geom_errorbar(aes(ymin = ac - ac_sd, ymax = ac + ac_sd),
                position = position_dodge(width = 0.9),  # 与条形图对齐
                width = 0.2, alpha = 0.5) +
  labs(fill = "Mat") +
  scale_fill_aaas()#+facet_wrap(~type, scales = "free_x")
ggsave("Figrue_acc_crap.pdf", P , width = 17, height = 8, dpi = 300)
ggsave("Figrue_acc_trout.pdf", P , width = 17, height = 8, dpi = 300)


crap_geno = fread("geno.txt",header = TRUE)
nrow(crap_geno)
ncol(crap_geno)-1

trout_geno = fread("myGD.txt",header = TRUE)
nrow(trout_geno)
ncol(trout_geno)-1


zt = checking(pac_mean_sd,cal = "ACC")
zt1 = zt
zt1$MODE = factor(zt1$MODE,level = c("GBLUP","FCN","CNN","MNNDR"))
zt1$type = factor(zt1$type,level = c("ST","MT"))
zt1[Vtype == "ac",Vtype:="Individual accuracy"]
zt1[Vtype == "auc",Vtype:="AUC"]
zt1$Vtype = factor(zt1$Vtype,level = c("Individual accuracy","AUC"))
P <- ggplot(data = zt1, aes(x = MODE, y = value, group = Vtype, fill = Vtype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
  geom_text(aes(label = round(value,2)), position = position_dodge(width = 0.9),vjust = -2, size = 4) +  # 添加数值标签
  xlab("Model") +
  ylab("Value") +
  theme_zg() +
  theme(legend.position ="right", legend.title = element_blank()) +
  geom_errorbar(aes(ymin = value - value_sd, ymax = value + value_sd),
                position = position_dodge(width = 0.9),  # 与条形图对齐
                width = 0.2, alpha = 0.5) +
  labs(fill = "Mat") +
  scale_fill_aaas()+facet_wrap(~type, scales = "free_x")
ggsave("Figrue_acc_crap_all.pdf", P , width = 17, height = 8, dpi = 300)
ggsave("Figrue_acc_trout_all.pdf", P , width = 17, height = 8, dpi = 300)


value_type = c("ac","auc")
type = c("ST","MT")
output_test = data.table()

for(t in type){

for(v in value_type){

  temp_dt = pac[type == t,]
  
  aov_res <- aov(as.formula(paste(v, "~ MODE")), data = temp_dt)
  
  tukey_res <- TukeyHSD(aov_res,conf.level = 0.95)
  
  letters <- multcompLetters4(aov_res, tukey_res)
  
  output_test = rbind(output_test,data.table(MODE = names(letters$MODE$Letters),Letter = letters$MODE$Letters,type = t,Vtype = v))
   
}
}
zt = checking(pac_mean_sd,cal = "ACC")
zt1 <- merge(output_test, zt, by = c("MODE", "type","Vtype"), all = FALSE)
zt1$MODE = factor(zt1$MODE,level = c("GBLUP","FCN","CNN","MNNDR"))
zt1$type = factor(zt1$type,level = c("ST","MT"))
zt1[Vtype == "ac",Vtype:="Individual accuracy"]
zt1[Vtype == "auc",Vtype:="AUC"]
zt1$Vtype = factor(zt1$Vtype,level = c("Individual accuracy","AUC"))
zt1[, letter_y_value := value + 1.1 * value_sd, by = Vtype]
zt1[, letter_y_test := value + 3.5 * value_sd, by = Vtype]
zt1[Vtype == "AUC", letter_y_test := value + 5.1 * max(value_sd)]
P <- ggplot(data = zt1, aes(x = MODE, y = value, group = MODE, fill = MODE)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
  geom_text(aes(label = round(value,2),y = letter_y_value), position = position_dodge(width = 0.9),vjust = 0, size = 4) +  # 添加数值标签
  xlab("Model") +
  ylab("Value") +
  theme_zg() +
  theme(legend.position ="right", legend.title = element_blank()) +
  geom_errorbar(aes(ymin = value - value_sd, ymax = value + value_sd),
                position = position_dodge(width = 0.9),  # 与条形图对齐
                width = 0.2, alpha = 0.5) +
  labs(fill = "Mat") +
  scale_fill_aaas()+geom_text(aes(label = Letter, y = letter_y_test), position = position_dodge(width = 0.9),vjust = 0) + facet_wrap(~type+Vtype, scales = "free")+
  ggtitle(label = "ST prediction",subtitle = "Simulated dataset (high non-additive)") +
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
ggsave("Figrue_acc_crap_all.pdf", P , width = 17, height = 8, dpi = 300)


# zt = checking(wac_mean_sd,pac_mean_sd,cal = "VAR")
# zt1 = zt[Vtype!="fwcova",]
# zt1$MODE = factor(zt1$MODE,level = c("PBLUP","GBLUP","FCN","CNN","MNN"))
# zt1$type = factor(zt1$type,level = c("ST","MT_continuous","MT_binary"))
# zt1[Vtype == "va",Vtype:="Population-wide"]
# zt1[Vtype == "fva",Vtype:="Between-family"]
# zt1[Vtype == "wva",Vtype:="Within-family"]
# zt1$Vtype = factor(zt1$Vtype,level = c("Population-wide","Between-family","Within-family"))
# P <- ggplot(data = zt1, aes(x = MODE, y = value, group = Vtype, fill = Vtype)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
#   geom_text(aes(label = round(value,2)), position = position_dodge(width = 0.9),vjust = -2, size = 4) +  # 添加数值标签
#   xlab("Model") +
#   ylab("Variance") +
#   theme_zg() +
#   theme(legend.position ="right", legend.title = element_blank()) +
#   geom_errorbar(aes(ymin = value - value_sd, ymax = value + value_sd),
#                 position = position_dodge(width = 0.9),  # 与条形图对齐
#                 width = 0.2, alpha = 0.5) +
#   labs(fill = "Mat") +
#   scale_fill_aaas()+facet_wrap(~type, scales = "free_x")
# ggsave("Figrue_var_pure.pdf", P , width = 15, height = 6, dpi = 300)
# ggsave("Figrue_var_low.pdf", P , width = 15, height = 6, dpi = 300)
# ggsave("Figrue_var_medium.pdf", P , width = 15, height = 6, dpi = 300)
# ggsave("Figrue_var_high.pdf", P , width = 15, height = 6, dpi = 300)




# zt = checking(wac_mean_sd,pac_mean_sd,cal = "COV")
# zt1 = zt[!Vtype%in%c("fwcov","wfcov"),]
# zt1$MODE = factor(zt1$MODE,level = c("PBLUP","GBLUP","FCN","CNN","MNN"))
# zt1$type = factor(zt1$type,level = c("ST","MT_continuous","MT_binary"))
# zt1[Vtype == "accov",Vtype:="Population-wide"]
# zt1[Vtype == "fcov",Vtype:="Between-family"]
# zt1[Vtype == "wcov",Vtype:="Within-family"]
# zt1$Vtype = factor(zt1$Vtype,level = c("Population-wide","Between-family","Within-family"))
# P <- ggplot(data = zt1, aes(x = MODE, y = value, group = Vtype, fill = Vtype)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.9) +  # 确保条形图分组不重叠
#   geom_text(aes(label = round(value,2)), position = position_dodge(width = 0.9),vjust = -2, size = 4) +  # 添加数值标签
#   xlab("Model") +
#   ylab("Variance") +
#   theme_zg() +
#   theme(legend.position ="right", legend.title = element_blank()) +
#   geom_errorbar(aes(ymin = value - value_sd, ymax = value + value_sd),
#                 position = position_dodge(width = 0.9),  # 与条形图对齐
#                 width = 0.2, alpha = 0.5) +
#   labs(fill = "Mat") +
#   scale_fill_aaas()+facet_wrap(~type, scales = "free_x")

# ggsave("Figrue_cov_pure.pdf", P , width = 15, height = 6, dpi = 300)
# ggsave("Figrue_cov_low.pdf", P , width = 15, height = 6, dpi = 300)
# ggsave("Figrue_cov_medium.pdf", P , width = 15, height = 6, dpi = 300)
# ggsave("Figrue_cov_high.pdf", P , width = 15, height = 8, dpi = 300)
