# champion grants

library(tidyverse)

champs = read.csv("~/Desktop/qcbs champions 2026/names.csv", col.names = "name")

df = table(champs) |> as.data.frame()
df = df[order(df$Freq, decreasing = T),]

choice = df[which(df$Freq>1),]
write.csv(choice, "~/Desktop/qcbs2026-champions.csv")
