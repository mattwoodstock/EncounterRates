library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(viridis)
library(ggridges)

# Create a color-blind friendly palette with 9 colors
color_blind_friendly_palette_8 <- viridis_pal()(8)
color_blind_friendly_palette_12 <- viridis_pal()(12)

## Add in visual ranges
setwd("D:/pz_encounter/outputs") #Set your working directory

dat <- read.csv("Prey visual range results.csv")

dat <- dat %>% mutate(Vis_rat = Pred_vis/Prey_vis,Abun_rat = PredN/PreyN) %>% rename(Velo_rat = Speed) %>% mutate(Vis_dist = round((3*(100000*Prey_vis)/4*3.14)^(1/3),0),Velo_rat = round(Velo_rat,2))

encounter_density <- ggplot(dat,aes(x=PreyN,y=Encounter/300,color=as.factor(Velo_rat)))+
  geom_errorbar(aes(ymin=dat$Encounter/300-dat$Enc_sd/300,ymax=dat$Encounter/300+dat$Enc_sd/300))+
  geom_point(size=1.5)+
  geom_smooth(fill=NA)+
  facet_wrap(.~as.factor(Vis_dist))+
  scale_color_manual(values = color_blind_friendly_palette_8, name = "Prey Abundance (n)")+
  labs(x="Prey abundance (n)", y = "Encounter Rate (N per timestep)")+
  theme_classic()+
  xlim(0,900)+
  theme(legend.title = element_blank())


ggsave("Encounter_rate_prey_abundance.png",dpi=300,encounter_density,height = unit(7,"in"),width=unit(12,"in"))

dens_velo <- ggplot(dat,aes(x=Velo_rat,y=Encounter/300,color=as.factor(PreyN)))+
  geom_errorbar(aes(ymin=dat$Encounter/300-dat$Enc_sd/300,ymax=dat$Encounter/300+dat$Enc_sd/300))+
  geom_vline(xintercept = 1)+
  geom_point(size=1.5)+
  geom_smooth(fill=NA)+
  facet_wrap(.~as.factor(Vis_dist))+
  labs(x="Swimming Speed Ratio (Predator/Prey)", y = "Encounter Rate (N per timestep)")+
  theme_classic()+
  scale_x_continuous(
    trans = "log10",
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  scale_color_manual(values = color_blind_friendly_palette_12, name = "Prey Abundance (n)")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),legend.title = element_blank())+
  guides(fill=guide_legend(title="Preys"))

ggsave("Encounter_rate_prey_velo.png",dpi=300,dens_velo,height = unit(7,"in"),width=unit(12,"in"))


ggplot(dat,aes(x=Velo_rat,y=Encounter/300,color=as.factor(PreyN)))+
  geom_vline(xintercept = 1)+
  geom_point(size=1.5)+
  #geom_smooth(fill=NA)+
  facet_wrap(.~as.factor(Vis_dist))+
  labs(x="Swimming Speed Ratio (Predator/Prey)", y = "Encounter Rate (N per timestep)")+
  theme_classic()+
  #scale_x_continuous(
  #  trans = "log10",
  #  labels = scales::trans_format("log10", scales::math_format(10^.x))
  #)+
  #scale_color_manual(values = color_blind_friendly_palette_12, name = "Prey Abundance (n)")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),legend.title = element_blank())+
  guides(fill=guide_legend(title="Preys"))

ggsave("Encounter rate_prey velo.png",dpi=300,dens_velo,height = unit(7,"in"),width=unit(12,"in"))


## Histogram of prey densities
setwd("D:/pz_encounter/outputs")
dat1 <- read.csv("Prey visual range results.csv")
dat1 <- dat1 %>% mutate(Simulation = 1:nrow(dat1))

#Set condition for histogram plot
all_dat <- c()

for (a in 1:nrow(dat1)){
  dat2 <- read.csv("Prey Densities.csv",header=F)[1:max(dat1$Simulation)][,a] %>% as.data.frame()
  colnames(dat2) <- a
  dat2_pivot <- dat2 %>% pivot_longer(cols = c(colnames(dat2)),values_to = "Density",names_to = "Simulation") %>% filter(Density != 0)
  new_dat2_pivot <- merge(dat2_pivot,dat1,by=c("Simulation"))
  all_dat <- rbind(all_dat,new_dat2_pivot)
}
all_dat2 <- all_dat %>% mutate(Prop_dens = Density/PreyN)

unique(all_dat2$Prey_vis)
ggplot(all_dat,aes(x=Density,fill=as.factor(Prey_vis)))+
  geom_density(alpha = 0.5)+
  #geom_histogram(binwidth = 5,color="black")+
  theme_classic()+
  facet_wrap(.~PreyN)+
  xlim(0,max(all_dat$PreyN))+
  labs(x="Percievable Preys (n)",y="Density")

all_dat2[(all_dat2$Prey_vis == 0.4) & (all_dat2$Prey_velo == 1e-4),]

count <- 1
setwd("D:/pz_encounter/outputs/Densities")

library(cowplot)
for (vis in 1:n_distinct(all_dat2$Prey_vis)){
  all_dat4 <- all_dat2 %>% filter(Prey_vis == unique(all_dat2$Prey_vis)[vis])
  p <- list()
  for (velo in 1:n_distinct(all_dat2$Prey_velo)){
    dens <- ggplot(all_dat2[all_dat2$Prey_velo == unique(all_dat2$Prey_velo)[velo],],aes(x=Density,y=as.factor(PreyN),fill=as.factor(PreyN)))+
      geom_density_ridges(scale = 3, rel_min_height = 0.01)+
      #geom_histogram(bins=30)+
      theme_classic()+
      #facet_wrap(Prey_velo~.,scales = "free")+
      #xlim(0,max(all_dat2$PreyN))+
      scale_fill_manual(values = color_blind_friendly_palette_17, name = "Prey Abundance (n)")+
      theme(legend.title =  element_blank(),legend.position = "none")+
      labs(x="Percievable Preys Per Predator (n)",y="Prey Abundance (n)",title=paste("Prey_velo: ",unique(all_dat2$Prey_velo)[velo],sep=""))+
      coord_cartesian(xlim=c(0,NA))
    
    p[velo] <- list(dens = dens)
  }
  
  combined_plot <- cowplot::plot_grid(plotlist = p)
  
  ggsave(paste("Prey Vis_",unique(all_dat2$Prey_vis)[vis],".png",sep=""),dpi=300,combined_plot,height = unit(6,"in"),width=unit(10,"in"))
}

setwd("D:/pz_encounter/outputs")

prey_density <- ggplot(all_dat2[(all_dat2$Prey_vis == 0.4) & (all_dat2$Prey_velo == 1e-4),],aes(x=Density,y=as.factor(PreyN),fill=as.factor(PreyN)))+
  geom_density_ridges(scale = 3, rel_min_height = 0.01)+
  #geom_histogram(bins=30)+
  theme_classic()+
  #facet_wrap(Prey_velo~.,scales = "free")+
  #xlim(0,max(all_dat2$PreyN))+
  scale_fill_manual(values = color_blind_friendly_palette_17, name = "Prey Abundance (n)")+
  theme(legend.title =  element_blank(),legend.position = "none")+
  labs(x="Percievable Preys Per Predator (n)",y="Prey Abundance (n)")+
  coord_cartesian(xlim=c(0,NA))

prop_density <- ggplot(all_dat2[(all_dat2$Prey_vis == 0.4) & (all_dat2$Prey_velo == 1e-4),],aes(x=Prop_dens,y=as.factor(PreyN),fill=as.factor(PreyN)))+
  geom_vline(xintercept = 0.5)+
  geom_density_ridges()+
  #geom_histogram(bins=30)+
  theme_classic()+
  #facet_wrap(Prey_velo~.,scales = "free")+
  #xlim(0,max(all_dat2$PreyN))+
  labs(x="Proportion of Preys Percievable",y="Prey Abundance (n)")+
  scale_fill_manual(values = color_blind_friendly_palette_17, name = "Prey Abundance (n)")+
  theme(legend.title =  element_blank(),legend.position = c(0.9,0.5))+
  guides(color = guide_legend(title.position = "bottom", title.hjust = 0.5))+  # Flip the legend
  xlim(0,1)

densities <-ggarrange(prey_density,prop_density,labels=c("A","B"))
ggsave("Combined_Prey_densities.png",dpi=300,densities,height = unit(6,"in"),width=unit(8,"in"))
  
ggsave("Prop_densities.png",dpi=300,prop_density,height = unit(7,"in"),width=unit(12,"in"))

all_dat3 <- all_dat2 %>% group_by(Prey_velo,PreyN,Prey_vis) %>% summarize(sd(Density)) %>% rename(SD = `sd(Density)`) %>% mutate(Speed = round(1e-4/Prey_velo,2),Vis_dist = round((3*(100000*Prey_vis)/4*3.14)^(1/3),0)) %>% filter(Vis_dist %in% c(0,3,46,62))

density_var <- ggplot(all_dat3,aes(x=PreyN,y=SD,color=as.factor(Speed)))+
  geom_point(size=1.5)+
  geom_smooth(method = "lm",se=F,fill=NA)+
  facet_wrap(.~as.factor(Vis_dist),nrow=1)+
  theme_classic()+
  theme(legend.title = element_blank())+
  labs(x="Prey abundance (n)", y = "SD in Percievable Preys Per Predator Distribution")+
  scale_color_manual(values = color_blind_friendly_palette_8, name = "Prey Abundance (n)")

full_density <- ggarrange(densities,density_var,ncol=1,labels=c(NA,"C"))

ggsave("Full_Prey_Density_Variance.png",dpi=300,full_density,height = unit(10,"in"),width=unit(8,"in"))

## Estimate attack rates from data ----

# Define your functional response curve equation
functional_response <- function(a, P, h) {
  enc <- a * P / (1 + a * h * P)
  return(enc)
}

# Define a function to minimize (sum of squared differences between observed and predicted values)
objective_function <- function(a, observed_enc, P, h) {
  predicted_enc <- functional_response(a, P, h)
  sum_squared_diff <- sum((observed_enc - predicted_enc)^2)
  return(sum_squared_diff)
}

df <- data.frame(Speed = numeric(),Vis = numeric(),a = numeric(),Prey_vis = numeric(), Prey_velo = numeric())

for (i in 1:n_distinct(dat$Velo_rat)){
  for (j in 1:n_distinct(dat$Vis_dist)){
    sub_dat <- dat %>% filter(Velo_rat == unique(dat$Velo_rat)[i], Vis_dist == unique(dat$Vis_dist)[j])
    observed_enc <- sub_dat$Encounter/300  # Your actual encounter rates
    P <- sub_dat$PreyN  # Your prey density vector
    h <- 8/60
    
    initial_guess_of_a <- 0
    
    result <- optim(par = initial_guess_of_a, fn = objective_function, 
                    observed_enc = observed_enc, P = P, h = h,lower=0,upper=5000,method = "L-BFGS-B")
    
    # Optimal value for 'a'
    optimal_a <- result$par

    # Generate data for plotting
    plot_data <- data.frame(P = seq(0, max(P), length.out = 100))
    plot_data$EncounterRates <- functional_response(optimal_a, plot_data$P, h)
    
    df <- df %>% add_row(Speed = unique(dat$Velo_rat)[i],Vis = unique(dat$Vis_dist)[j],a = optimal_a,Prey_vis = unique(dat$Prey_vis)[j],Prey_velo = unique(dat$Prey_velo)[i])
}
}
df$a[df$a == 0] <- NA

ggplot(df)+
  geom_tile(aes(x=as.factor(Speed),y=as.factor(Vis),fill=a))+
  scale_fill_gradient2(low="white",high="blue")+
  labs(x= "Relative Swim Speed (Predator/Prey)",y = "Prey Visual Range (m)")+
  annotate("text", x = 2.5, y = 8, label = "Random Encounters\nLikely Unrealistic Interaction")+
  annotate("text",x=as.factor(df$Speed),y=as.factor(df$Vis),label=round(df$a,3))

## Consumption plot ----

first <- all_dat2 %>% group_by(PreyN,Prey_vis,Prey_velo) %>% summarize(quantile(Density,0.025)) %>% rename(first = `quantile(Density, 0.025)`)

mid <- all_dat2 %>% group_by(PreyN,Prey_vis,Prey_velo) %>% summarize(quantile(Density,0.5)) %>% rename(mid = `quantile(Density, 0.5)`)

third <- all_dat2 %>% group_by(PreyN,Prey_vis,Prey_velo) %>% summarize(quantile(Density,0.975)) %>% rename(third = `quantile(Density, 0.975)`)

new_dat <- merge(first,third,by=c("PreyN","Prey_vis","Prey_velo"))
new_dat2 <- merge(new_dat,mid,by=c("PreyN","Prey_vis","Prey_velo"))
head(df)
head(new_dat2)
new_dat3 <- merge(new_dat2,df,by=c("Prey_vis","Prey_velo"))

consumption = new_dat3 %>% mutate(HighQ = 0.7*5*5000*(a*(third/500))/(1+(8*60)*a*(third/500)),LowQ = 0.7*5*5000*(a*(first/500))/(1+(8*60)*a*(first/500))) %>% mutate(Grid = 0.7*5*5000*(a*((PreyN/1000)))/(1+(8*60)*a*((PreyN/1000))), MidQ = 0.7*5*5000*(a*(mid/500))/(1+(8*60)*a*(mid/500))) %>% filter(Vis != 0,Speed >= 2) 

head(onsumption)
consumption_plot <- ggplot(consumption,aes(x=PreyN))+
  geom_line(aes(x=PreyN,y=consumption$LowQ))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_line(aes(x=PreyN,y=consumption$Grid),color="red",linewidth=1.5)+
  geom_line(aes(x=PreyN,y=consumption$MidQ),color="black",linewidth=1.5)+
  geom_ribbon(aes(ymin = consumption$LowQ,ymax=consumption$HighQ),alpha = 0.2)+
  facet_grid(Vis~Speed)+
  labs(x="Prey Abundance",y="Consumption (J per timestep)")+
  theme_classic()

ggsave("Consumption_Rates.png",dpi=300,consumption_plot,height = unit(7,"in"),width=unit(12,"in"))
### Working area ----

