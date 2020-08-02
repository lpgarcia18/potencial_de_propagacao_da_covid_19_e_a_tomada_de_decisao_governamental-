# Ambiente ----------------------------------------------------------------
options(scipen=999)
gc()
set.seed(1)


# Pacotes -----------------------------------------------------------------
library(readr)
library(tidyverse)
library(EpiEstim)
#devtools::install_github("tidyverse/googlesheets4")
library(googlesheets4)
library(caTools)
library(reshape2)

# Key functions - Funções extraídas de  https://www.biorxiv.org/content/10.1101/835181v2.full.pdf  
# Código: https://github.com/kpzoo/model-selection-for-epidemic-renewal-models
source('apeEstim.R')
source('apePredPost.R')


# Total de infectados por data de início de sintomas. Dados com nowcasting, extraídos de https://github.com/geinfosms/covid_geinfo
covid <- read_csv("dados/covid.csv")


res_base_diarios <- list()

for(i in seq_along(covid$INICIO_SINTOMAS)){
  
	# Cortando a base em janelas de 31 dias ----------------------------------------------------------
	incidencia <- covid[i:(31+i),] #Utilizando dados de um mês
	
	# Truncagem à esquerda
	trunc <- 2
	
	Icovid <- incidencia$MEDIANA_CASOS
	gencovid <- EpiEstim::discr_si(c(0:max(incidencia$MEDIANA_CASOS)), mu = 4.8, sigma = 2.3)
	Lcovid = overall_infectivity(Icovid, gencovid)
	
	
	# Priors
	Rprior = c(1, 5); a = 0.025 #Confidence interval level
	# Completando valores NA
	Lcovid[is.na(Lcovid)] = 0# <------ importante
	
	
	# Melhor janela e predições do Rt 
	Rmodcovid = apeEstim(Icovid, gencovid, Lcovid, Rprior, a, trunc, 'covid')
	Rcovid <- Rmodcovid[[2]][[4]]
	RcovidCI_025 <- Rmodcovid[[2]][[5]][1,]
	RcovidCI_975 <- Rmodcovid[[2]][[5]][2,]
	DATA <- tail(incidencia$INICIO_SINTOMAS,-2)
	res_base <- data.frame(DATA = DATA, Rt = Rcovid, Rt_025 = RcovidCI_025, Rt_975 = RcovidCI_975)

	res_base_diarios[[i]] <- res_base
        names(res_base_diarios)[[i]] <-  as.character(tail(res_base$DATA,1))
 }


#Unindo com base do decreto
decretos <- read_csv("dados/classificacao_decretos.csv")
anterior_publicacao <- (as.Date(decretos$PUBLICACAO, format = "%d/%m/%Y") - 1) %>% as.character() #Selecionando Rt anterior a data da publicação do decreto
res_base_diarios_sel <- subset(res_base_diarios, as.character(names(res_base_diarios)) %in% as.character(anterior_publicacao)) 
res_base_diarios_sel <- res_base_diarios_sel %>% 
	lapply(function(x)tail(x,14)) %>%
	lapply(function(x)mutate(x,DEFASAGEM = c("D_14","D_13","D_12","D_11",
						 "D_10","D_9","D_8","D_7","D_6",
						 "D_5","D_4","D_3","D_2","D_1"))) %>%
	lapply(function(x)mutate(x,DATA = NULL)) %>%
	lapply(function(x)melt(x, id.var = "DEFASAGEM"))  %>%
	lapply(function(x)dcast(x,formula = . ~ DEFASAGEM + variable , fun.aggregate = NULL))

	
names(res_base_diarios_sel) <- (as.Date(names(res_base_diarios_sel))+1) #Fazendo com que o nome das bases seja igual à data de publicação dos decretos para juntar as bases
res_base_diarios_sel <- do.call(rbind,res_base_diarios_sel)
res_base_diarios_sel$DATA <- rownames(res_base_diarios_sel)
decretos$PUBLICACAO <- as.character(as.Date(decretos$PUBLICACAO, format = "%d/%m/%Y"))
decretos_analies <- merge(decretos, res_base_diarios_sel, by.x = "PUBLICACAO", by.y = "DATA", all= T)
decretos_analies$CLASSIFICACAO <- as.character(decretos_analies$CLASSIFICACAO) 
for(i in seq_along(decretos_analies$PUBLICACAO)){
	  if(
	     ((decretos_analies$D_1_Rt_975[i] >= 1 | #R>=1 ao menos 1 dia e Aumenta medidas de distanciamento = Consonante
	       decretos_analies$D_2_Rt_975[i] >= 1 |
	       decretos_analies$D_3_Rt_975[i] >= 1 |
	       decretos_analies$D_4_Rt_975[i] >= 1 |
	       decretos_analies$D_5_Rt_975[i] >= 1 |
	       decretos_analies$D_6_Rt_975[i] >= 1 |
	       decretos_analies$D_7_Rt_975[i] >= 1 |
	       decretos_analies$D_8_Rt_975[i] >= 1 |
	       decretos_analies$D_9_Rt_975[i] >= 1 |
	       decretos_analies$D_10_Rt_975[i] >= 1 |
	       decretos_analies$D_11_Rt_975[i] >= 1 |
	       decretos_analies$D_12_Rt_975[i] >= 1 |
	       decretos_analies$D_13_Rt_975[i] >= 1 |
	       decretos_analies$D_14_Rt_975[i] >= 1 ) & (
	       decretos_analies$CLASSIFICACAO[i] == "Aumenta medidas de distanciamento")) |
	     ((decretos_analies$D_1_Rt_975[i] < 1 & #R<1 todos os dias e Flexibiliza medidas de distanciamento = Consonante
	       decretos_analies$D_2_Rt_975[i] < 1 &
	       decretos_analies$D_3_Rt_975[i] < 1 &
	       decretos_analies$D_4_Rt_975[i] < 1 &
	       decretos_analies$D_5_Rt_975[i] < 1 &
	       decretos_analies$D_6_Rt_975[i] < 1 &
	       decretos_analies$D_7_Rt_975[i] < 1 &
	       decretos_analies$D_8_Rt_975[i] < 1 &
	       decretos_analies$D_9_Rt_975[i] < 1 &
	       decretos_analies$D_10_Rt_975[i] < 1 &
	       decretos_analies$D_11_Rt_975[i] < 1 &
	       decretos_analies$D_12_Rt_975[i] < 1 &
	       decretos_analies$D_13_Rt_975[i] < 1 &
	       decretos_analies$D_14_Rt_975[i] < 1) & (
	       decretos_analies$CLASSIFICACAO[i] == "Flexibiliza medidas de distanciamento"))
	     ){decretos_analies$ANALISE[i] <- "Consonante"
	      }else{
	       decretos_analies$ANALISE[i] <- "Dissonante"
	}	
}
	
decr <- decretos_analies
decr$D_1_Rt <- round(decr$D_1_Rt,3)
decr$D_1_Rt_025 <- round(decr$D_1_Rt_025,3)
decr$D_1_Rt_975 <- round(decr$D_1_Rt_975,3)

decr$D_2_Rt <- round(decr$D_2_Rt,3)
decr$D_2_Rt_025 <- round(decr$D_2_Rt_025,3)
decr$D_2_Rt_975 <- round(decr$D_2_Rt_975,3)

decr$D_3_Rt <- round(decr$D_3_Rt,3)
decr$D_3_Rt_025 <- round(decr$D_3_Rt_025,3)
decr$D_3_Rt_975 <- round(decr$D_3_Rt_975,3)

decr$D_4_Rt <- round(decr$D_4_Rt,3)
decr$D_4_Rt_025 <- round(decr$D_4_Rt_025,3)
decr$D_4_Rt_975 <- round(decr$D_4_Rt_975,3)

decr$D_5_Rt <- round(decr$D_5_Rt,3)
decr$D_5_Rt_025 <- round(decr$D_5_Rt_025,3)
decr$D_5_Rt_975 <- round(decr$D_5_Rt_975,3)

decr$D_6_Rt <- round(decr$D_6_Rt,3)
decr$D_6_Rt_025 <- round(decr$D_6_Rt_025,3)
decr$D_6_Rt_975 <- round(decr$D_6_Rt_975,3)

decr$D_7_Rt <- round(decr$D_7_Rt,3)
decr$D_7_Rt_025 <- round(decr$D_7_Rt_025,3)
decr$D_7_Rt_975 <- round(decr$D_7_Rt_975,3)

decr$D_8_Rt <- round(decr$D_8_Rt,3)
decr$D_8_Rt_025 <- round(decr$D_8_Rt_025,3)
decr$D_8_Rt_975 <- round(decr$D_8_Rt_975,3)

decr$D_9_Rt <- round(decr$D_9_Rt,3)
decr$D_9_Rt_025 <- round(decr$D_9_Rt_025,3)
decr$D_9_Rt_975 <- round(decr$D_9_Rt_975,3)

decr$D_10_Rt <- round(decr$D_10_Rt,3)
decr$D_10_Rt_025 <- round(decr$D_10_Rt_025,3)
decr$D_10_Rt_975 <- round(decr$D_10_Rt_975,3)

decr$D_11_Rt <- round(decr$D_11_Rt,3)
decr$D_11_Rt_025 <- round(decr$D_11_Rt_025,3)
decr$D_11_Rt_975 <- round(decr$D_11_Rt_975,3)

decr$D_12_Rt <- round(decr$D_12_Rt,3)
decr$D_12_Rt_025 <- round(decr$D_12_Rt_025,3)
decr$D_12_Rt_975 <- round(decr$D_12_Rt_975,3)

decr$D_13_Rt <- round(decr$D_13_Rt,3)
decr$D_13_Rt_025 <- round(decr$D_13_Rt_025,3)
decr$D_13_Rt_975 <- round(decr$D_13_Rt_975,3)

decr$D_14_Rt <- round(decr$D_14_Rt,3)
decr$D_14_Rt_025 <- round(decr$D_14_Rt_025,3)
decr$D_14_Rt_975 <- round(decr$D_14_Rt_975,3)


decr$D1 <- paste0(decr$D_1_Rt," (", decr$D_1_Rt_025, "-", decr$D_1_Rt_975, ")")
decr$D2 <- paste0(decr$D_2_Rt," (", decr$D_2_Rt_025, "-", decr$D_2_Rt_975, ")")
decr$D3 <- paste0(decr$D_3_Rt," (", decr$D_3_Rt_025, "-", decr$D_3_Rt_975, ")")
decr$D4 <- paste0(decr$D_4_Rt," (", decr$D_4_Rt_025, "-", decr$D_4_Rt_975, ")")
decr$D5 <- paste0(decr$D_5_Rt," (", decr$D_5_Rt_025, "-", decr$D_5_Rt_975, ")")
decr$D6 <- paste0(decr$D_6_Rt," (", decr$D_6_Rt_025, "-", decr$D_6_Rt_975, ")")
decr$D7 <- paste0(decr$D_7_Rt," (", decr$D_7_Rt_025, "-", decr$D_7_Rt_975, ")")
decr$D8 <- paste0(decr$D_8_Rt," (", decr$D_8_Rt_025, "-", decr$D_8_Rt_975, ")")
decr$D9 <- paste0(decr$D_9_Rt," (", decr$D_9_Rt_025, "-", decr$D_9_Rt_975, ")")
decr$D10 <- paste0(decr$D_10_Rt," (", decr$D_10_Rt_025, "-", decr$D_10_Rt_975, ")")
decr$D11 <- paste0(decr$D_11_Rt," (", decr$D_11_Rt_025, "-", decr$D_11_Rt_975, ")")
decr$D12 <- paste0(decr$D_12_Rt," (", decr$D_12_Rt_025, "-", decr$D_12_Rt_975, ")")
decr$D13 <- paste0(decr$D_13_Rt," (", decr$D_13_Rt_025, "-", decr$D_13_Rt_975, ")")
decr$D14 <- paste0(decr$D_14_Rt," (", decr$D_14_Rt_025, "-", decr$D_14_Rt_975, ")")


#Escrevendo o resultado
write.csv(decr, "resultados/resultado.csv", row.names = F)


 

