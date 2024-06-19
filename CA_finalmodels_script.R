library(ggplot2)
library(deSolve)
library(reshape2)


plot(raster(nec_abund_mat))
nec.grid <- nec_abund_mat
nec.grid.bin <- nec_abund_mat
nec.grid.bin[nec.grid.bin<2] <- 0
nec.grid.bin[nec.grid.bin>=2] <- 1
plot(raster(nec.grid))
plot(raster(nec.grid.bin))

plot(raster(oligo_abund_mat))
oligo.grid <- oligo_abund_mat
oligo.grid.bin <- oligo_abund_mat
oligo.grid.bin[oligo.grid.bin<2] <- 0
oligo.grid.bin[oligo.grid.bin>=2] <- 1
plot(raster(oligo.grid.bin))


setwd('D:/pesquisas/Zoonotic modelling/rasters/hanta')
hanta_cases <- raster('casos_res10km.tif')
plot(hanta_cases)
casos.matrix <- as.matrix(hanta_cases)
casos.matrix[is.na(casos.matrix)] <- 0

plot(raster(casos.matrix))
max(casos.matrix)

hist(casos.matrix)
casos.mat.bin <- casos.matrix
casos.mat.bin[casos.mat.bin<4.5] <- 0
casos.mat.bin[casos.mat.bin>=4.5] <- 1
plot(raster(casos.mat.bin))
casos.bin.raster <- raster(casos.mat.bin,
                           xmn=-1523073.7859432739205658,
                           xmx=3070100.7613228280097246,
                           ymn=-214981.2677185528154951,
                           ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
setwd('D:/pesquisas/Zoonotic modelling/resultados')
writeRaster(casos.bin.raster,
            'casos.bin.raster.tif',
            format='GTiff',
            overwrite = T)

library(sf)
setwd('D:/pesquisas/Zoonotic modelling/shapefiles/hanta')
casos.teste <- raster(casos.mat.bin,
                      xmn=-1523073.7859432739205658,
                      xmx=3070100.7613228280097246,
                      ymn=-214981.2677185528154951,
                      ymx=4209259.7714078342542052,
                      crs="ESRI:102033")
casos.points <- st_read('centroids.shp')
setwd('D:/pesquisas/Zoonotic modelling/shapefiles/brasil_adm')
brasil.adm <- st_read('brasil_adm_albers.shp')
brasil.mun <- st_read('mun_hanta.shp')
new.bb <- c(xmin=-1523073.7859432739205658,
            xmax=3070100.7613228280097246,
            ymin=-214981.2677185528154951,
            ymax=4209259.7714078342542052)
attr(new.bb,"class") <- "bbox"
attr(st_geometry(brasil.mun),"bbox") <- new.bb

setwd('D:/pesquisas/Zoonotic modelling/shapefiles/brasil_biomas')
brasil.biomas <- st_read('brasil_biomas_albers.shp')
plot(casos.teste)
plot(brasil.adm[3],  col=rgb(1, 0, 0,0), add=T)
plot(brasil.mun[1], col=rgb(1, 0, 0,0), add=T)
plot(casos.points,add=T)


nec_casos <- nec.grid.bin+casos.mat.bin #matriz de onde houve óbitos e necromys está presente
nec_casos[nec_casos==1] <- 0
nec_casos[is.na(nec_casos)] <- 0
nec_casos[nec_casos==2] <- 1
plot(raster(nec_casos)) 

oligo_casos <- oligo.grid.bin+casos.mat.bin #matriz de onde houve óbitos e necromys está presente
oligo_casos[oligo_casos==1] <- 0
oligo_casos[is.na(oligo_casos)] <- 0
oligo_casos[oligo_casos==2] <- 1
plot(raster(oligo_casos)) 

#primeiro calcular a contaminação em cada célula no momento t
sir_equations <- function(time, variables, parameters) 
{
  with(as.list(c(variables, parameters)), 
       {
         dS <- -beta * I * S
         dI <-  beta * I * S - gamma * I
         dR <-  gamma * I
         return(list(c(dS, dI, dR)))
       })
}



#Valores de beta estimados a partir de Figueiredo et al 2010 para Necromys
##---
##Considerando os resultados desses autores, três taxas de infecção
##anual foram calculadas: 
##2005-2006-> 0,124938737
##2006-2007-> 0,176091259
##2007-2008-> 0,56427143
##Para cada taxa foi calculada uma curva de crescimento logístico, e uma 
##taxa de infecção média anual foi então calculada em 0,2884
##Isso equivale a 0,02404 ind/mês ou 0,000801 ind/dia


parameters_values <- c(beta  = 0.02404, # infectious contact rate (/ind/month)
                       gamma = 0.0001      # recovery rate (/ind/month) #não temos dados suficientes
                       )


#Agora deixando prontas as funções que vão rodar os modelos
#de automatos celulares
num_neighbors <- function(grid, i, j)
{
  sum(grid[max(1, i-1):min(nrow(grid), i+1),
           max(1, j-1):min(ncol(grid), j+1)]) - grid[i, j]
}


update_grid <- function(grid)
{
  new_grid <- grid
  for (i in 1:nrow(grid))
  {
    for (j in 1:ncol(grid))
    {
      neighbors <- num_neighbors(grid, i, j)
      new_grid[i, j] <- ifelse(grid[i, j] == 1, #verifique se a célula está contaminada
                               ifelse(neighbors <= 3, 0, 1), #se estiver, cheque se tem mais de 2 células contaminadas ao redor
                               ifelse(neighbors > 3, 1, 0)) #se não, cheque se tem 3 ou mais células contaminadas ao redor
    }
  }
  new_grid
}


#MODELOS CONJUGADOS
##Necromys lasiurus
time_values <- seq(0,120) #tempo (em meses)

nec.teste <- nec_casos
nec.teste.list <- list(nec.teste)
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 121,    # Maximum value of the progress bar
                     style = 3,    
                     width = 50,   
                     char = "=")   
nec.sir_values<-list()
for(i in 1:nrow(nec.grid))
{
  nec.sir_values[[i]] <- list()
  for(j in 1:ncol(nec.grid))
  {
    for(t in 1)
    {
      if(nec.teste[i,j]==1)
      {
        initial_values <- c(S = nec.grid[i,j]/2-0.1,   # number of susceptibles at time = 0
                            I = nec.grid[i,j]/2+0.1,  # number of infectious at time = 0
                            R = 0)   # number of recovered (and immune) at time = 0
        nec.sir_values[[i]][[j]] <- ode(y = initial_values,
                                    times = time_values,
                                    func = sir_equations,
                                    parms = parameters_values)
      }
      else{nec.sir_values[[i]][[j]]<-0}
    }
  }
}
plot(raster(nec.teste.list[[1]]))

nec.grid_zeros <- nec.grid
nec.grid_zeros[is.na(nec.grid_zeros)] <- 0

for(t in 1)
{
  for(i in 1:nrow(nec.grid))
  {
    for(j in 1:ncol(nec.grid))
    {
      if(nec.teste[i,j]!=1){next}
      else
      {
        if(nec.grid[i,j]>max(nec.grid_zeros)/2) #se a abundância local for maior que metade da abundância máxima esperada
          {
          ifelse((nec.sir_values[[i]][[j]][t,2]+nec.sir_values[[i]][[j]][t,4])<nec.sir_values[[i]][[j]][t,3],
               nec.teste[i,j] <- 1,
               nec.teste[i,j] <- 0)
          }
      }
    }
  }
  teste.temp <- update_grid(nec.teste)
  nec.teste.list[[t+1]] <- teste.temp
  plot(raster(nec.teste.list[[t+1]]))
}

for(t in 2:(max(time_values)+1))
{
  for(i in 1:length(nec.sir_values))
  {
    for(j in 1:length(nec.sir_values[[i]]))
    {
      temp <- nec.teste.list[[t]][i,j]-nec.teste.list[[t-1]][i,j]
      if(temp==1)
      {
        initial_values <- c(S = nec.grid[i,j]-nec.grid[i,j]*0.003,   # number of susceptibles at time = 0
                            I = nec.grid[i,j]*0.003,  # number of infectious at time = 0
                            R = 0)   # number of recovered (and immune) at time = 0
        nec.sir_values[[i]][[j]] <- ode(y = initial_values,
                                    times = time_values,
                                    func = sir_equations,
                                    parms = parameters_values)
      }
      if(is.list(nec.sir_values[[i]][[j]]))
      {
        if(nec.grid[i,j]>max(nec.grid_zeros)/2)
          {
        ifelse(nec.sir_values[[i]][[j]][t,2]+nec.sir_values[[i]][[j]][t,4]<nec.sir_values[[i]][[j]][t,3],
               nec.teste[i,j] <- 1,
               nec.teste[i,j] <- 0)
          }
      }
    }
  }
  nec.teste <- update_grid(nec.teste)
  nec.teste.list[[t+1]] <- nec.teste
  plot(raster(nec.teste.list[[t+1]]))
  setTxtProgressBar(pb,t)
}
close(pb)

nec.diff_120_0 <- nec.teste.list[[121]]-nec.teste.list[[1]]
nec.diff_48_0 <- nec.teste.list[[49]]-nec.teste.list[[1]]
nec.dif_teste <- nec.teste.list[[121]]+nec_casos
nec_diff <- nec.dif_teste + nec.grid.bin
plot(raster(nec.diff_120_0))
plot(raster(nec.diff_48_0))
plot(raster(nec_diff))
plot(raster(nec.dif_teste))

sum(nec.teste.list[[1]]==1)
sum(nec.teste.list[[121]]==1)
nrow(nec.teste.list[[1]])*ncol(nec.teste.list[[1]])

nec.teste.sum.list <- data.frame()
for (i in 1:length(nec.teste.list))
{
  nec.teste.sum.list[i,1] <- i
  nec.teste.sum.list[i,2] <- sum(nec.teste.list[[i]]==1)
}
colnames(nec.teste.sum.list) <- c('timesteps','pixels')
ggplot(data=nec.teste.sum.list[-3,],aes(x=timesteps, y=pixels))+
  geom_line(col = 'blue')+
  geom_point()+
  labs(title = 'Necromys lasiurus')+
  ylab("Number of infected pixels")+
  xlab("Timesteps (months)")+
  scale_x_continuous(breaks = seq(0,120,20))+
  geom_vline(xintercept = 23, col = 'red')+
  theme_classic()


#xmn=-1522695.0816999999806285,
#xmx=3069574.6892999997362494,
#ymn=-214549.7362999999895692,
#ymx=4209011.1821999996900558

##EXTENSÃO DOS MAPAS DO BRASIL NO QGIS
#xmn=-1523073.7859432739205658,
#xmx=3070100.7613228280097246,
#ymn=-214981.2677185528154951,
#ymx=4209259.7714078342542052
nec.start.raster <- raster(nec_casos,
                           xmn=-1523073.7859432739205658,
                           xmx=3070100.7613228280097246,
                           ymn=-214981.2677185528154951,
                           ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
nec.final.raster <- raster(nec.teste.list[[121]],
                           xmn=-1523073.7859432739205658,
                           xmx=3070100.7613228280097246,
                           ymn=-214981.2677185528154951,
                           ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
nec.48.raster <- raster(nec.teste.list[[49]],
                        xmn=-1523073.7859432739205658,
                        xmx=3070100.7613228280097246,
                        ymn=-214981.2677185528154951,
                        ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
nec.diff.raster1 <- raster(nec.diff_120_0,
                           xmn=-1523073.7859432739205658,
                           xmx=3070100.7613228280097246,
                           ymn=-214981.2677185528154951,
                           ymx=4209259.7714078342542052,
                          crs="ESRI:102033")
nec.diff.raster2 <- raster(nec.diff_48_0,
                           xmn=-1523073.7859432739205658,
                           xmx=3070100.7613228280097246,
                           ymn=-214981.2677185528154951,
                           ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
nec.diff.raster3 <- raster(nec_diff,
                           xmn=-1523073.7859432739205658,
                           xmx=3070100.7613228280097246,
                           ymn=-214981.2677185528154951,
                           ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
nec.diff.raster4 <- raster(nec.dif_teste,
                           xmn=-1523073.7859432739205658,
                           xmx=3070100.7613228280097246,
                           ymn=-214981.2677185528154951,
                           ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
plot(nec.diff.raster1)
plot(nec.diff.raster2)
plot(nec.diff.raster3)
plot(nec.diff.raster4)
plot(brasil.adm[3],  col=rgb(1, 0, 0,0), add=T)
plot(brasil.mun[1], col=rgb(1, 0, 0,0), add=T)
#plot(brasil.biomas,add=T,col=rgb(1,0,0,0))
#plot(nec.sir_values[[76]][[233]][,1],nec.sir_values[[76]][[233]][,2])
#plot(nec.sir_values[[76]][[233]][,1],nec.sir_values[[76]][[233]][,3],add=T)
setwd('D:/pesquisas/Zoonotic modelling/resultados/necromys')

writeRaster(nec.start.raster,
            't0.tif',
            format='GTiff',
            overwrite = T)
writeRaster(nec.48.raster,
            't48.tif',
            format='GTiff',
            overwrite = T)
writeRaster(nec.final.raster,
            't120.tif',
            format='GTiff',
            overwrite = T)
writeRaster(nec.diff.raster1,
            'expansao.tif',
            format='GTiff',
            overwrite = T)
writeRaster(nec.diff.raster2,
            'dif_inicial.tif',
            format='GTiff',
            overwrite = T)
writeRaster(nec.diff.raster3,
            'dif_inicial.tif',
            format='GTiff',
            overwrite = T)
writeRaster(nec.diff.raster4,
            'nec_t0_expansao.tif',
            format='GTiff',
            overwrite = T)

##Oligoryzomys nigripes
time_values <- seq(0,120) #tempo (em meses)

oligo.teste <- oligo_casos
plot(raster(oligo.teste))
oligo.teste.list <- list(oligo.teste)
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 121, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
oligo.sir_values<-list()
for(i in 1:nrow(oligo.grid))
{
  oligo.sir_values[[i]] <- list()
  for(j in 1:ncol(oligo.grid))
  {
    for(t in 1)
    {
      if(oligo.teste[i,j]==1)
      {
        initial_values <- c(S = oligo.grid[i,j]/2-0.1,   # number of susceptibles at time = 0
                            I = oligo.grid[i,j]/2+0.1,  # number of infectious at time = 0
                            R = 0)   # number of recovered (and immune) at time = 0
        oligo.sir_values[[i]][[j]] <- ode(y = initial_values,
                                    times = time_values,
                                    func = sir_equations,
                                    parms = parameters_values)
      }
      else{oligo.sir_values[[i]][[j]]<-0}
    }
  }
}
plot(raster(oligo.teste.list[[1]]))

oligo.grid_zeros <- oligo.grid
oligo.grid_zeros[is.na(oligo.grid_zeros)] <- 0

for(t in 1)
{
  for(i in 1:nrow(oligo.grid))
  {
    for(j in 1:ncol(oligo.grid))
    {
      if(oligo.teste[i,j]!=1){next}
      else
      {
        if(oligo.grid[i,j]>max(oligo.grid_zeros)/2)
        {
        ifelse((oligo.sir_values[[i]][[j]][t,2]+oligo.sir_values[[i]][[j]][t,4])<oligo.sir_values[[i]][[j]][t,3],
               oligo.teste[i,j] <- 1,
               oligo.teste[i,j] <- 0)
        }
      }
    }
  }
  teste.temp <- update_grid(oligo.teste)
  oligo.teste.list[[t+1]] <- teste.temp
  plot(raster(oligo.teste.list[[t+1]]))
}

for(t in 2:(max(time_values)+1))
{
  for(i in 1:length(oligo.sir_values))
  {
    for(j in 1:length(oligo.sir_values[[i]]))
    {
      temp <- oligo.teste.list[[t]][i,j]-oligo.teste.list[[t-1]][i,j]
      if(temp==1)
      {
        initial_values <- c(S = oligo.grid[i,j]-oligo.grid[i,j]*0.02404,   # number of susceptibles at time = 0
                            I = oligo.grid[i,j]*0.02404,  # number of infectious at time = 0
                            R = 0)   # number of recovered (and immune) at time = 0
        oligo.sir_values[[i]][[j]] <- ode(y = initial_values,
                                    times = time_values,
                                    func = sir_equations,
                                    parms = parameters_values)
      }
      if(is.list(oligo.sir_values[[i]][[j]]))
      {
        if(oligo.grid[i,j]>max(oligo.grid_zeros)/2)
        {
        ifelse(oligo.sir_values[[i]][[j]][t,2]+oligo.sir_values[[i]][[j]][t,4]<oligo.sir_values[[i]][[j]][t,3],
               oligo.teste[i,j] <- 1,
               oligo.teste[i,j] <- 0)
        }
      }
    }
  }
  oligo.teste <- update_grid(oligo.teste)
  oligo.teste.list[[t+1]] <- oligo.teste
  plot(raster(oligo.teste.list[[t+1]]))
  setTxtProgressBar(pb,t)
}
close(pb)

oligo.diff_120_0 <- oligo.teste.list[[121]]-oligo.teste.list[[1]]
oligo.diff_48_0 <- oligo.teste.list[[49]]-oligo.teste.list[[1]]
oligo.dif_teste <- oligo.teste.list[[121]]+oligo_casos
oligo_diff <- oligo.dif_teste + oligo.grid.bin
plot(raster(oligo.diff_120_0))
plot(raster(oligo.diff_48_0))
plot(raster(oligo_diff))

sum(oligo.teste.list[[1]]==1)
sum(oligo.teste.list[[121]]==1)
nrow(oligo.teste.list[[1]])*ncol(oligo.teste.list[[1]])

oligo.teste.sum.list <- data.frame()
for (i in 1:length(oligo.teste.list))
{
  oligo.teste.sum.list[i,1] <- i
  oligo.teste.sum.list[i,2] <- sum(oligo.teste.list[[i]]==1)
}
colnames(oligo.teste.sum.list) <- c('timesteps','pixels')
ggplot(data=oligo.teste.sum.list[-3,],aes(x=timesteps, y=pixels))+
  geom_line(col='blue')+
  geom_point()+
  labs(title = "Oligoryzomys nigripes")+
  ylab("Number of infected pixels")+
  xlab("Timesteps (months)")+
  scale_x_continuous(breaks = seq(0,120,20))+
  geom_vline(xintercept = 12, col = 'red')+
  theme_classic()

oligo.start.raster <- raster(oligo_casos,
                             xmn=-1523073.7859432739205658,
                             xmx=3070100.7613228280097246,
                             ymn=-214981.2677185528154951,
                             ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
oligo.final.raster <- raster(oligo.teste.list[[121]],
                             xmn=-1523073.7859432739205658,
                             xmx=3070100.7613228280097246,
                             ymn=-214981.2677185528154951,
                             ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
oligo.diff.raster1 <- raster(oligo.diff_120_0,
                             xmn=-1523073.7859432739205658,
                             xmx=3070100.7613228280097246,
                             ymn=-214981.2677185528154951,
                             ymx=4209259.7714078342542052,
                             crs="ESRI:102033")
oligo.diff.raster2 <- raster(oligo.diff_48_0,
                             xmn=-1523073.7859432739205658,
                             xmx=3070100.7613228280097246,
                             ymn=-214981.2677185528154951,
                             ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
oligo.diff.raster3 <- raster(oligo_diff,
                             xmn=-1523073.7859432739205658,
                             xmx=3070100.7613228280097246,
                             ymn=-214981.2677185528154951,
                             ymx=4209259.7714078342542052,
                           crs="ESRI:102033")
oligo.diff.raster4 <- raster(oligo.dif_teste,
                           xmn=-1523073.7859432739205658,
                           xmx=3070100.7613228280097246,
                           ymn=-214981.2677185528154951,
                           ymx=4209259.7714078342542052,
                           crs="ESRI:102033")

setwd('D:/pesquisas/Zoonotic modelling/resultados/oligoryzomys')
writeRaster(oligo.start.raster,
            't0.tif',
            format='GTiff',
            overwrite = T)
writeRaster(oligo.final.raster,
            't120.tif',
            format='GTiff',
            overwrite = T)
writeRaster(oligo.diff.raster1,
            'expansao_0_120.tif',
            format='GTiff',
            overwrite = T)
writeRaster(oligo.diff.raster2,
            'expansao_4anos_0_48.tif',
            format='GTiff',
            overwrite = T)
writeRaster(oligo.diff.raster3,
            'dif_inicial.tif',
            format='GTiff',
            overwrite = T)
writeRaster(oligo.diff.raster4,
            'dif_inicial.tif',
            format='GTiff',
            overwrite = T)

plot(oligo.diff.raster1)
plot(oligo.diff.raster2)
plot(oligo.diff.raster3)


#--------
## Ambas as espécies agregadas

nec_oligo_bin <- nec.grid.bin+oligo.grid.bin
plot(raster(nec_oligo_bin))
nec_oligo_bin[nec_oligo_bin>0] <- 1
nec_oligo_casos <- nec_oligo_bin + casos.mat.bin
nec_oligo_casos[nec_oligo_casos==1] <- 0
nec_oligo_casos[nec_oligo_casos==2] <- 1
nec_oligo_casos[is.na(nec_oligo_casos)] <- 0
plot(raster(nec_oligo_casos))

nec_oligo_grid <- nec.grid + oligo.grid
plot(raster(nec_oligo_grid))
#nrow(nec_oligo_casos)==nrow(nec_oligo_grid)
#ncol(nec_oligo_casos)==ncol(nec_oligo_grid)
time_values <- seq(0,120) #tempo (em meses)

nec.oligo.teste <- nec_oligo_casos
nec.oligo.teste.list <- list(nec.oligo.teste)
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 121, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
nec.oligo.sir_values<-list()
for(i in 1:nrow(nec_oligo_grid))
{
  nec.oligo.sir_values[[i]] <- list()
  for(j in 1:ncol(nec_oligo_grid))
  {
    for(t in 1)
    {
      if(nec.oligo.teste[i,j]==1)
      {
        initial_values <- c(S = nec_oligo_grid[i,j]/2-0.001,   # number of susceptibles at time = 0
                            I = nec_oligo_grid[i,j]/2+0.001,  # number of infectious at time = 0
                            R = 0)   # number of recovered (and immune) at time = 0
        nec.oligo.sir_values[[i]][[j]] <- ode(y = initial_values,
                                          times = time_values,
                                          func = sir_equations,
                                          parms = parameters_values)
      }
      else{nec.oligo.sir_values[[i]][[j]]<-0}
    }
  }
}
plot(raster(nec.oligo.teste.list[[1]]))

for(t in 1)
{
  for(i in 1:nrow(nec_oligo_grid))
  {
    for(j in 1:ncol(nec_oligo_grid))
    {
      if(nec.oligo.teste[i,j]!=1){next}
      else
      {
        ifelse((nec.oligo.sir_values[[i]][[j]][t,2]+nec.oligo.sir_values[[i]][[j]][t,4])<nec.oligo.sir_values[[i]][[j]][t,3],
               nec.oligo.teste[i,j] <- 1,
               nec.oligo.teste[i,j] <- 0)
      }
    }
  }
  teste.temp <- update_grid(nec.oligo.teste)
  nec.oligo.teste.list[[t+1]] <- teste.temp
  plot(raster(nec.oligo.teste.list[[t+1]]))
}

for(t in 2:(max(time_values)+1))
{
  for(i in 1:length(sir_values))
  {
    for(j in 1:length(nec.oligo.sir_values[[i]]))
    {
      temp <- nec.oligo.teste.list[[t]][i,j]-nec.oligo.teste.list[[t-1]][i,j]
      if(temp==1)
      {
        initial_values <- c(S = nec_oligo_grid[i,j]-nec_oligo_grid[i,j]*0.0205,   # number of susceptibles at time = 0
                            I = nec_oligo_grid[i,j]*0.0205,  # number of infectious at time = 0
                            R = 0)   # number of recovered (and immune) at time = 0
        nec.oligo.sir_values[[i]][[j]] <- ode(y = initial_values,
                                          times = time_values,
                                          func = sir_equations,
                                          parms = parameters_values)
      }
      if(is.list(nec.oligo.sir_values[[i]][[j]]))
      {
        ifelse(nec.oligo.sir_values[[i]][[j]][t,2]+nec.oligo.sir_values[[i]][[j]][t,4]<nec.oligo.sir_values[[i]][[j]][t,3],
               nec.oligo.teste[i,j] <- 1,
               nec.oligo.teste[i,j] <- 0)
      }
    }
  }
  nec.oligo.teste <- update_grid(nec.oligo.teste)
  nec.oligo.teste.list[[t+1]] <- nec.oligo.teste
  plot(raster(nec.oligo.teste.list[[t+1]]))
  setTxtProgressBar(pb,t)
}
close(pb)

plot(raster(nec.oligo.teste.list[[121]]+nec.oligo.teste.list[[1]]))

plot(raster(nec.oligo.teste.list[[121]]+nec.grid.bin+oligo.grid.bin))
