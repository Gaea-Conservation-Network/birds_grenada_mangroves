### Grenada Mangrove Isotopes

data<-read.csv("C:/Users/crrem/Dropbox/documents/grenada moangroves/2019 Grenada isotopes.csv")
head(data)
library(ggplot2)

data2<-(data[-c(24:28),])
plants<-subset(data2, Type=="plant")
data2_<-c()

ggplot(data2, aes(x=X18O, y=X2H, col=species))+
  geom_point(size=4)+
  scale_color_brewer(palette="Pastel1")+
  geom_abline(aes(intercept = 10, slope = 8), linetype = "dashed")+
  facet_wrap(~Site)+
  theme_bw()



