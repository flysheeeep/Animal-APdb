### get trait-related ATSS

# tissue_dir = as.character(args[1])
# loadpath <- as.character(args[2])
# trait_file <- as.character(args[3])

# tissue_dir: promoter usage dir
# tissue_info: 'sample tissue information from SRA'
# trait_file: 'sample matadata from SRA'


species = gsub('_trait1','', tools::file_path_sans_ext(basename(trait_file)))
trait_file = read.table(trait_file,sep = '\t',header = T)

tissue_info$sample <- tools::file_path_sans_ext(basename(tissue_info$run)) 



##'@trait_analysis


# creat list tissue2sex to save tissues samlpes having sex trait
sex_trait = dplyr::filter(trait_file, sex != 'Unclear')
n=1
sex_trait$number = NA
for (i in unique(sex_trait$tissue)) {
  sex_trait[which(sex_trait$tissue == i),]$number = length(unique(sex_trait[which(sex_trait$tissue == i),]$sex))
}

sex_trait = sex_trait[which(sex_trait$number == 2),]

get_sex = function(sex_trait){
  for (i in unique(sex_trait$tissue)) {
    
    #filter tissue
    trait = sex_trait[which(sex_trait$tissue == i),]
    tissue_sex = trait$sex
    if (length(tissue_sex[which(tissue_sex == 'M')]) < 3 | length(tissue_sex[which(tissue_sex == 'F')]) < 3) {
      print(paste0('no sex trait in tissue: ',i))
      next
    }else {
      n = n+1
      samples = trait$run
      sex_M = trait[trait$sex == 'M',]$run
      sex_F = trait[trait$sex == 'F',]$run
      tissue_rel = paste0(tissue_dir,'/',i,'_modi.txt')
      tissue_rel = read.table(file = tissue_rel,header = T)
      sheet1 = tissue_rel[,1:7]
      sheet1$tissue = i
      sheet1$trait = 'sex'
      sheet1$FC = NA
      sheet1$FDR = NA
      for (j in 1:nrow(tissue_rel)) {
        male = as.numeric(tissue_rel[j,][sex_M])
        female = as.numeric(tissue_rel[j,][sex_F])
        if (length(which(male != 0 )) < 3 | length(which(female != 0 )) < 3 ) {
          print(paste0('no sex trait in tissue 1: ',i))
          next
        }else if (setequal(sort(male), sort(female))) {
          print(paste0('no sex trait in tissue 2: ',i))
          next
        }else {
          test = wilcox.test(male,female)
          FDR = test$p.value
          FC = mean(male)/mean(female)
          sheet1[j,]$FC = FC
          sheet1[j,]$FDR =FDR
          
        }
      }
      
      sheet1 = na.omit(sheet1)
      sheet1$FDR<-p.adjust(sheet1$FDR,method="fdr",length(sheet1$FDR))
      #filter by FDR < 0.05, FC >1.5
      sheet1[abs(sheet1$FC)<=1.5,'FC']<-NA
      sheet1[sheet1$FDR>=0.05,'FDR']<-NA
      sheet1 = na.omit(sheet1)
      if (nrow(sheet1) > 0) {
        write.table(sheet1,paste0('/path/to/save',species,'/',i,'_Gender.txt'),sep = '\t',quote = F,row.names = F)
      }
    }
    
    
  }
}
get_sex(sex_trait = sex_trait)


# Trait2: Embryo vs Postnatal
tissue2age = as.character()
m=1

age_trait = dplyr::filter(trait_file, ageFlag != 'Unclear')

get_age = function(age_trait){
  for (i in unique(age_trait$tissue)) {
    
    #filter tissue
    trait2 = age_trait[which(age_trait$tissue == i),]
    tissue_age = trait2$ageFlag
    if (length(tissue_age[which(tissue_age == 'Postnatal')]) < 3 | length(tissue_age[which(tissue_age == 'Embryo')]) < 3) {
      print(paste0('no age trait in tissue: ',i))
      next
    }else {
      m = m+1
      samples = trait2$run
      age_Em = trait2[trait2$ageFlag == 'Embryo',]$run
      age_Po = trait2[trait2$ageFlag == 'Postnatal',]$run
      tissue_rel = paste0(tissue_dir,'/',i,'_modi.txt')
      tissue_rel = read.table(file = tissue_rel,header = T)
      sheet2 = tissue_rel[,1:7]
      sheet2$tissue = i
      sheet2$trait = 'age'
      sheet2$Age_FC = NA
      sheet2$Age_FDR = NA
      for (j in 1:nrow(tissue_rel)) {
        Embryo = as.numeric(tissue_rel[j,][age_Em])
        Postnatal = as.numeric(tissue_rel[j,][age_Po])
        if (length(which(Embryo != 0 )) < 3 | length(which(Postnatal != 0 )) < 3 ) {
          print(paste0('no age trait in tissue: ',i))
          next
        }else if (setequal(sort(Embryo), sort(Postnatal))) {
          print(paste0('no age trait in tissue: ',i))
          next
        }else {
          test = wilcox.test(Embryo,Postnatal)
          FDR = test$p.value
          FC = mean(Embryo)/mean(Postnatal)
          sheet2[j,]$Age_FC = FC
          sheet2[j,]$Age_FDR =FDR
        }
      }
      sheet2 = na.omit(sheet2)
      sheet2$Age_FDR<-p.adjust(sheet2$Age_FDR,method="fdr",length(sheet2$Age_FDR))
      #filter by FDR < 0.05, FC >1.5
      sheet2[abs(sheet2$Age_FC)<=1.5,'Age_FC']<-NA
      sheet2[sheet2$Age_FDR>=0.05,'Age_FDR']<-NA
      sheet2 = na.omit(sheet2)
      if (nrow(sheet2)>0) {
        write.table(sheet2,paste0('/path/to/save',species,'/',i,'_Age.txt'),sep = '\t',quote = F,row.names = F)
      }
    }
  }
}
get_age(age_trait = age_trait)

#stage=2:wilcox.test, stage=3|4: ANNOVA, stage >= 5: Spearman correlation
#Trait3: Embryo

tissue2Embryo = as.character()

get_ageE = function(age_trait){
  for (i in unique(age_trait$tissue)) {
    #filter tissue
    trait3 = age_trait[which(age_trait$tissue == i),]
    trait3 = trait3[which(trait3$ageFlag == 'Embryo'),]
    tissue_Emb = unique(trait3$ageFlagE)
    tissue_Emb = tissue_Emb[tissue_Emb != -1]
    
    tissue_rel = paste0(tissue_dir,'/',i,'_modi.txt')
    tissue_rel = read.table(file = tissue_rel,header = T)
    if (length(tissue_Emb) == 1) {
      print(paste0('no Embryo trait in tissue 1: ',i))
      next
    }else if (length(tissue_Emb) ==2 ) {
      sheet3 = tissue_rel[,1:7]
      sheet3$tissue = i
      sheet3$trait = 'embryo'
      sheet3$FC = NA
      sheet3$FDR = NA
      sheet3$Flag = 1
      
      #wilcox.test
      Emb_1 = trait3[which(trait3$ageFlagE == sort(tissue_Emb)[1]),]$run
      Emb_2 = trait3[which(trait3$ageFlagE == sort(tissue_Emb)[2]),]$run
      
      for (j in 1:nrow(tissue_rel)) {
        Emb_1.t = as.numeric(tissue_rel[j,][Emb_1])
        Emb_2.t = as.numeric(tissue_rel[j,][Emb_2])
        if (length(which(Emb_1.t != 0 )) < 3 | length(which(Emb_2.t != 0 )) < 3 ) {
          print(paste0('no Embryo trait in tissue 2: ',i))
          next
        }else if (setequal(sort(Emb_1.t), sort(Emb_2.t))) {
          print(paste0('no Embryo trait in tissue 2: ',i))
          next
        }else {
          test = wilcox.test(Emb_1.t,Emb_2.t)
          FDR = test$p.value
          FC = mean(Emb_1.t)/mean(Emb_2.t)
          sheet3[j,]$FC = FC
          sheet3[j,]$FDR =FDR
        }
      }
      
      
      sheet3 = na.omit(sheet3)
      sheet3$FDR<-p.adjust(sheet3$FDR,method="fdr",length(sheet3$FDR))
      #filter by FDR < 0.05, FC >1.5
      sheet3[abs(sheet3$FC)<=1.5,'FC']<-NA
      sheet3[sheet3$FDR>=0.05,'FDR']<-NA
      sheet3 = na.omit(sheet3)
      if (nrow(sheet3)>0) {
        print(paste0(paste0('Embryo trait 2 in tissue: ',i)))
        write.table(sheet3,paste0('/path/to/save',species,'/',i,'_Embryo.txt'),sep = '\t',quote = F,row.names = F)
      }
    }else if (length(tissue_Emb) >= 3 ) {
      
      #Spearman
      sheet3 = tissue_rel[,1:7]
      sheet3$tissue = i
      sheet3$trait = 'embryo'
      sheet3$Rho = NA
      sheet3$FDR = NA
      sheet3$Flag = 3
      
      
      tmp = data.frame(matrix(nrow = nrow(trait3[which(trait3$ageFlagE %in% tissue_Emb),]), ncol = 3))
      colnames(tmp) = c('run','ageFlagE','usage')
      tmp$run = trait3[which(trait3$ageFlagE %in% tissue_Emb),]$run
      for (h in 1:nrow(tmp)) {
        tmp[h,]$ageFlagE = trait3[which(trait3$run == tmp[h,]$run),]$ageFlagE
      }
      
      for (j in 1:nrow(tissue_rel)) {
        tmp$usage = as.numeric(tissue_rel[j,tmp$run])
        if (length(unique(tmp[which(tmp$usage != 0),]$ageFlagE)) < 3) {
          sheet3$Rho[j] = NA
          sheet3$FDR[j] = NA
          next
        }else if (length(tmp[which(tmp$usage != 0),]$run) < 5) {
          sheet3$Rho[j] = NA
          sheet3$FDR[j] = NA
          next
        }else{
          correlation = cor.test(x = as.numeric(tmp$usage), y = as.numeric(tmp$ageFlagE), alternative = "two.sided",method = "spearman")
          sheet3$Rho[j] = correlation[["estimate"]][["rho"]]
          sheet3$FDR[j] = correlation[["p.value"]]
        }
      }
      sheet3 = na.omit(sheet3)
      sheet3$FDR<-p.adjust(sheet3$FDR,method="fdr",length(sheet3$FDR))
      #filter by FDR < 0.05, FC >1.5
      sheet3[sheet3$FDR>=0.05,'FDR']<-NA
      sheet3[sheet3$Rho<0.3,'Rho']<-NA
      sheet3 = na.omit(sheet3)
      if (nrow(sheet3)>0) {
        print(paste0(paste0('Embryo trait 3 in tissue: ',i)))
        write.table(sheet3,paste0('/path/to/save',species,'/',i,'_Embryo.txt'),sep = '\t',quote = F,row.names = F)
      }
    }
    
  }
}

get_ageE(age_trait = age_trait)


#Trait4: Postnatal
get_ageP = function(age_trait){
  for (i in unique(age_trait$tissue)) {
    
    #filter tissue
    trait3 = age_trait[which(age_trait$tissue == i),]
    trait3 = trait3[which(trait3$ageFlag == 'Postnatal'),]
    tissue_Post = unique(trait3$ageFlagP)
    tissue_Post = tissue_Post[tissue_Post != -1]
    
    tissue_rel = paste0(tissue_dir,'/',i,'_modi.txt')
    tissue_rel = read.table(file = tissue_rel,header = T)
    
    if (length(tissue_Post) <= 1) {
      print(paste0('no postnatal trait in tissue: ',i))
      next
    }else if (length(tissue_Post) ==2 ) {
      sheet3 = tissue_rel[,1:7]
      sheet3$tissue = i
      sheet3$trait = 'postnatal'
      sheet3$FC = NA
      sheet3$FDR = NA
      sheet3$Flag = 1
      
      #wilcox.test
      Post_1 = trait3[which(trait3$ageFlagP == sort(tissue_Post)[1]),]$run
      Post_2 = trait3[which(trait3$ageFlagP == sort(tissue_Post)[2]),]$run
      if (length(Post_1) < 3 | length(Post_2) <3 ) {
        print(paste0('no postnatal trait in tissue: ',i))
      }else {
        for (j in 1:nrow(tissue_rel)) {
          Post_1.t = as.numeric(tissue_rel[j,][Post_1])
          Post_2.t = as.numeric(tissue_rel[j,][Post_2])
          if (length(which(Post_1.t != 0 )) < 3 | length(which(Post_2.t != 0 )) < 3 ) {
            next
          }else if (setequal(sort(Post_1.t), sort(Post_2.t))) {
            next
          }else {
            test = wilcox.test(Post_1.t,Post_2.t)
            FDR = test$p.value
            FC = mean(Post_1.t)/mean(Post_2.t)
            sheet3[j,]$FC = FC
            sheet3[j,]$FDR =FDR
          }
        }
        
        
        sheet3 = na.omit(sheet3)
        sheet3$FDR<-p.adjust(sheet3$FDR,method="fdr",length(sheet3$FDR))
        #filter by FDR < 0.05, FC >1.5
        sheet3[abs(sheet3$FC)<=1.5,'FC']<-NA
        sheet3[sheet3$FDR>=0.05,'FDR']<-NA
        sheet3 = na.omit(sheet3)
        if (nrow(sheet3)>0) {
          print(paste0(paste0('Post trait 1 in tissue: ',i)))
          write.table(sheet3,paste0('/path/to/save',species,'/',i,'_Postnatal.txt'),sep = '\t',quote = F,row.names = F)
        }
      }
      
    }else if (length(tissue_Post) > 2) {
      
      #Spearman
      sheet3 = tissue_rel[,1:7]
      sheet3$tissue = i
      sheet3$trait = 'postnatal'
      sheet3$Rho = NA
      sheet3$FDR = NA
      sheet3$Flag = 3
      
      
      tmp = data.frame(matrix(nrow = nrow(trait3[which(trait3$ageFlagP %in% tissue_Post),]), ncol = 3))
      colnames(tmp) = c('run','ageFlagP','usage')
      tmp$run = trait3[which(trait3$ageFlagP %in% tissue_Post),]$run
      for (h in 1:nrow(tmp)) {
        tmp[h,]$ageFlagP = trait3[which(trait3$run == tmp[h,]$run),]$ageFlagP
      }
      
      for (j in 1:nrow(tissue_rel)) {
        tmp$usage = as.numeric(tissue_rel[j,tmp$run])
        if (length(unique(tmp[which(tmp$usage != 0),]$ageFlagP)) < 3) {
          sheet3$Rho[j] = NA
          sheet3$FDR[j] = NA
          next
        }else if (length(tmp[which(tmp$usage != 0),]$run) < 5) {
          sheet3$Rho[j] = NA
          sheet3$FDR[j] = NA
          next
        }else {
          correlation = cor.test(x = as.numeric(tmp$usage), y = as.numeric(tmp$ageFlagP), alternative = "two.sided",method = "spearman")
          sheet3$Rho[j] = correlation[["estimate"]][["rho"]]
          sheet3$FDR[j] = correlation[["p.value"]]
        }
        
      }
      sheet3 = na.omit(sheet3)
      sheet3$FDR<-p.adjust(sheet3$FDR,method="fdr",length(sheet3$FDR))
      #filter by FDR < 0.05, FC >1.5
      sheet3[sheet3$FDR>=0.05,'FDR']<-NA
      sheet3[sheet3$Rho<0.3,'Rho']<-NA
      sheet3 = na.omit(sheet3)
      if (nrow(sheet3)>0) {
        print(paste0(paste0('Post trait 3 in tissue: ',i)))
        write.table(sheet3,paste0('/path/to/save',species,'/',i,'_Postnatal.txt'),sep = '\t',quote = F,row.names = F)
      }
    }
    
  }
}

get_ageP(age_trait = age_trait)
