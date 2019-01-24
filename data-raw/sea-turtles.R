### sea-turtles
### Example datasets on sea turtles

# Loggerhead sea turtle state-based life tables from Crouse et al. (1987) and Crowder et
# al. (1994)
loggerhead <- list(data.frame(stage = c("eggs", "small juveniles", "large juveniles",
                                        "subadults", "novice breeders",
                                        "1st-yr remigrants", "mature breeders"),
                              survival = c(0.6747, 0.7857, 0.6758, 0.7425,
                                           0.8091, 0.8091, 0.8091),
                              maternity = c(0, 0, 0, 0, 127, 4, 80),
                              duration = c(1, 7, 8, 6, 1, 1, Inf)),
                   data.frame(stage = c("eggs", "small juveniles", "large juveniles",
                                        "subadults", "adults"),
                              survival = c(0.6747, 0.75, 0.6758, 0.7425,
                                           0.8091),
                              maternity = c(0, 0, 0, 0, 76.5),
                              duration = c(1, 7, 8, 6, Inf)))
save(loggerhead, file = "data/loggerhead.rda")