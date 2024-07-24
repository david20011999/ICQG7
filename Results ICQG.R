rm(list = ls())
devtools::load_all()


#load necessary packages. They will be
pacman::p_load(pedigree, pedigreemm,ggplot2)

nIndfounder=1000 #individual founder population

nChrfounder=2 #chromosomes in founder population

segSites=1000 

nQTL=300

nSNP=segSites-nQTL

mean_tot=c(0,0,0)

var_tot= c(10,10,10)

cor_tot = matrix(data = c(1.00, 0.8, 0.40,
                          0.80, 1.00, 0.40,
                          0.40, 0.40, 1.00),
                 byrow = T, nrow=3,ncol=3)

mean_ID=c(0,0,0)

var_ID= c(9,9,9)

cor_ID = matrix(data = c(1.00, 0.8, 0.90,
                         0.80, 1.00, 0.90,
                         0.90, 0.90, 1.00),
                byrow = T, nrow=3,ncol=3)




#obtaining funder genomes
founderGenomes = quickHaplo(nInd =nIndfounder,
                         nChr = nChrfounder,
                         segSites= segSites)

#defining simulation parameters
SP=SimParam$new(founderGenomes)
SP$setSexes("yes_sys")
SP$addTraitAI(nQtlPerChr= nQTL ,mean=mean_tot,var=var_tot,meanID=mean_ID,
              varID=var_ID,corA = cor_tot,corID = cor_ID)
SP$addSnpChip(nSnpPerChr = nQTL)

a <- SP$traits[[1]]@addEff
d <- rep(0,600)
i <- SP$traits[[1]]@impEff

#defining base population
basePop <- newPop(founderGenomes)

parameters <- genParam(basePop)

#creating phenotypes
basePop <- setPheno(basePop,varE = c(10, 10, 10))

#first cross for separate types

#creating paternal type
PatMales = selectInd(basePop, 100, trait=1, use = "pheno", sex = "M")
PatDams = selectInd(basePop, 100, trait=1, use = "pheno", sex = "F")

MatMales = selectInd(basePop, 100, trait=2, use = "pheno", sex = "M")
MatDams = selectInd(basePop, 100, trait=2, use = "pheno", sex = "F")

for (generation in (1:50)){
  
  Patpop = randCross2(males = PatMales, females = PatDams, nCrosses = 1000, nProgeny = 1)
  Matpop = randCross2(males = MatMales, females = MatDams, nCrosses = 1000, nProgeny = 1)
  
  Patpop <- setPheno(Patpop, varE = c(10, 10, 10)) #obtain their phenotypes
  Matpop <- setPheno(Matpop, varE = c(10, 10, 10))
  
  # paternal type
  PatMales = selectInd(Patpop, 100, trait=1, use = "pheno", sex = "M")
  PatDams = selectInd(Patpop, 100, trait=1, use = "pheno", sex = "F")
  
  
  # maternal type
  MatMales = selectInd(Matpop, 100, trait=2, use = "pheno", sex = "M")
  MatDams = selectInd(Matpop, 100, trait=2, use = "pheno", sex = "F")
  
}
df <- NULL

for (reply in 1:10) {
  
  print(reply)
  
BreedA = randCross2(males = PatMales, females = PatDams, nCrosses = 1000, nProgeny = 1)
BreedB = randCross2(males = MatMales, females = MatDams, nCrosses = 1000, nProgeny = 1)


#BV with imprinting 

for (generation in (1:10)){

  #Purebreed selection
  
  BreedA_males_pure <- selectInd(BreedA, 50, use = "bv",trait = 1, sex = "M")
  BreedA_females_pure <- selectInd(BreedA, 50, use = "bv",trait = 1, sex = "F")
  
  BreedB_males_pure <- selectInd(BreedA, 50, use = "bv",trait = 1, sex = "M")
  BreedB_females_pure <- selectInd(BreedA, 50, use = "bv",trait = 1, sex = "F")
  
  #Crossbreeding selection
  
  BreedA_males_cross <- BreedA[!(BreedA@iid  %in% BreedA_males_pure@iid) & BreedA@sex=="M"]

  BreedB_females_cross <- BreedB[!(BreedB@iid  %in% BreedB_females_pure@iid) & BreedB@sex=="F"]
  
  BreedA_males_cross@ebv <- bvP(BreedA_males_cross)
  BreedB_females_cross@ebv <- bvM(BreedB_females_cross)
  
  BreedA_males_cross <- selectInd(BreedA_males_cross, 50, use = "ebv",trait = 1, sex = "M")
  
  BreedB_females_cross <- selectInd(BreedB_females_cross, 50, use = "ebv",trait = 1, sex = "F")

  #Purebreed reproduction
  BreedA = randCross2(males = BreedA_males_pure, females = BreedA_females_pure, nCrosses = 1000, nProgeny = 1)
  BreedB = randCross2(males = BreedB_males_pure, females = BreedB_females_pure, nCrosses = 1000, nProgeny = 1)

  #Crossbreed production
  
  Terminal = randCross2(males = BreedA_males_cross, females = BreedB_females_cross, nCrosses = 1000, nProgeny = 1)
  Geno <- pullQtlGeno(Terminal)
  Dom <- Geno
  Dom[Dom==2] <- 0
  Imp <- Dom + -2*pullQtlHaplo(Terminal,haplo = 1)*Dom
  
  performance <- Geno %*% matrix(a,ncol = 1)  + Dom %*% matrix(d,ncol = 1) + Imp %*% matrix(i,ncol = 1) 
  df <- rbind(df,data.frame(Generation = generation, GV = performance, Estimation = "AI"))
  
}


#BV Without imprinting

BreedA = randCross2(males = PatMales, females = PatDams, nCrosses = 1000, nProgeny = 1)
BreedB = randCross2(males = MatMales, females = MatDams, nCrosses = 1000, nProgeny = 1)

for (generation in (1:10)){
  
  #Purebreed selection
  
  BreedA_males_pure <- selectInd(BreedA, 50, use = "bv",trait = 1, sex = "M")
  BreedA_females_pure <- selectInd(BreedA, 50, use = "bv",trait = 1, sex = "F")
  
  BreedB_males_pure <- selectInd(BreedA, 50, use = "bv",trait = 1, sex = "M")
  BreedB_females_pure <- selectInd(BreedA, 50, use = "bv",trait = 1, sex = "F")
  
  #Crossbreeding selection
  
  BreedA_males_cross <- BreedA[!(BreedA@iid  %in% BreedA_males_pure@iid) & BreedA@sex=="M"]
  
  BreedB_females_cross <- BreedB[!(BreedB@iid  %in% BreedB_females_pure@iid) & BreedB@sex=="F"]
  
  BreedA_males_cross <- selectInd(BreedA_males_cross, 50, use = "bv",trait = 1, sex = "M")
  
  BreedB_females_cross <- selectInd(BreedB_females_cross, 50, use = "bv",trait = 1, sex = "F")
  
  #Purebreed reproduction
  BreedA = randCross2(males = BreedA_males_pure, females = BreedA_females_pure, nCrosses = 1000, nProgeny = 1)
  BreedB = randCross2(males = BreedB_males_pure, females = BreedB_females_pure, nCrosses = 1000, nProgeny = 1)
  
  #Crossbreed production
  
  Terminal = randCross2(males = BreedA_males_cross, females = BreedB_females_cross, nCrosses = 1000, nProgeny = 1)
  
  Geno <- pullQtlGeno(Terminal)
  Dom <- Geno
  Dom[Dom==2] <- 0
  Imp <- Dom + -2*pullQtlHaplo(Terminal,haplo = 1)*Dom
  
  performance <- Geno %*% matrix(a,ncol = 1)  + Dom %*% matrix(d,ncol = 1) + Imp %*% matrix(i,ncol = 1) 
  df <- rbind(df,data.frame(Generation = generation, GV = performance, Estimation = "A"))
  
}

}

aux_mean <- as.data.frame(t(tapply(df$GV, list(df$Estimation, df$Generation), mean)))
aux_sd <- as.data.frame(t(tapply(df$GV, list(df$Estimation, df$Generation), sd)))

df <- rbind(data.frame(Generation=rownames(aux_mean), mean=aux_mean$A, sd = aux_sd$A, Strategy="Additive"),
            data.frame(Generation=rownames(aux_mean), mean=aux_mean$AI, sd = aux_sd$AI, Strategy="Additive + \nImprinting"))


ggplot(data=df, aes(x=as.numeric(Generation), y=mean, group = Strategy, col=Strategy, fill=Strategy)) + 
  geom_line(linewidth = 1.5) + 
  geom_ribbon(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd),alpha = 0.4,color = "transparent")+
  scale_x_continuous(breaks=seq(1, 10, 1))+
  labs(y ="GV", x = "Generation")

  

