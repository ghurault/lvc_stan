# Data --------------------------------------------------------------------

process_data <- function(dfA, dfB){
  # Process data
  #
  # Args:
  # dfA: dataframe corresponding to experiment A
  # dfB: dataframe corresponding to experiment B
  #
  # Returns:
  # Dataframe of the data
  
  df1 <- na.omit(rbind(data.frame(Age = dfA[, "Age"], Volume = dfA[, "Saccharomyces"], Experiment = "A"),
                       data.frame(Age = dfB[, "Age"], Volume = dfB[, "Saccharomyces"], Experiment = "B")))
  
  df2 <- na.omit(rbind(data.frame(Age = dfA[, "Age"], Volume = dfA[, "Schixosachararomyces"], Experiment = "A"),
                       data.frame(Age = dfB[, "Age"], Volume = dfB[, "Schixosachararomyces"], Experiment = "B")))
  
  df1$Species <- "Saccharomyces"
  df2$Species <- "Schixosachararomyces"
  
  df <- rbind(df1, df2)
  df$Condition <- "Single"
  
  df12 <- na.omit(rbind(data.frame(Age = dfA[, "Age"], Sa = dfA[, "Saccharomyces.mixed"], Sc = dfA[, "Schixosacharomyces.mixed"], Experiment = "A"),
                        data.frame(Age = dfB[, "Age"], Sa = dfB[, "Saccharomyces.mixed"], Sc = dfB[, "Schixosacharomyces.mixed"], Experiment = "B")))
  
  df12 <- reshape2::melt(df12, id.vars = c("Age", "Experiment"), variable.name = "Species", value.name = "Volume")
  df12$Species <- factor(df12$Species, levels = c("Sa", "Sc"), labels = c("Saccharomyces", "Schixosachararomyces"))
  df12$Condition <- "Mixed"
  
  df <- rbind(df, df12)
  df$Condition <- factor(df$Condition, levels = c("Single", "Mixed"))
  return(df)
}

plot_data <- function(df){
  # Plot data
  #
  # Args:
  # df: dataframe of data
  #
  # Returns:
  # Ggplot
  
  library(ggplot2)
  
  ggplot(data = df,
         aes(x = Age, y = Volume, colour = Experiment)) +
    geom_point(size = 2) + geom_line() +
    facet_grid(rows = vars(Species), cols = vars(Condition)) +
    scale_colour_manual(values = c("#000000", "#E69F00")) +
    theme_bw(base_size = 20) + theme(legend.position = "top")
}

# Posterior Predictive Checks ---------------------------------------------

process_replications_density <- function(data_stan, fit, maxVolume = 20){
  # Summarise replications as densities
  #
  # Args:
  # fit: stanfit object
  # data_stan: data input to stan
  # maxVolume: maximum volume value to consider
  #
  # Returns:
  # Dataframe of replications
  
  t <- c(0, data_stan$t_rep)
  s_name <- c("Saccharomyces", "Schixosachararomyces")
  
  extract_density <- function(tmp, s, cond){
    do.call("rbind",
            lapply(1:ncol(tmp), function(x){
              dA <- density(tmp[, x, 1], kernel = "gaussian", from = 0, to = maxVolume, n = 128)
              dB <- density(tmp[, x, 2], kernel = "gaussian", from = 0, to = maxVolume, n = 128)
              
              rbind(data.frame(Age = t[x], Volume = dA$x, Density = dA$y, Species = s, Condition = cond, Experiment = "A"),
                    data.frame(Age = t[x], Volume = dB$x, Density = dB$y, Species = s, Condition = cond, Experiment = "B"))
            })
    )
  }
  
  yrep <- extract(fit, pars = c("y1_rep", "y2_rep"))
  rep_single <- lapply(1:2,
                       function(i){extract_density(yrep[[i]], s_name[i], "Single")})
  
  yrep <- extract(fit, pars = c("y12_rep"))[[1]]
  rep_mixed <- lapply(1:2, function(i){extract_density(yrep[, , i, ], s_name[i], "Mixed")})
  
  do.call(rbind, c(rep_single, rep_mixed))
}

process_replications_spaghetti <- function(data_stan, fit, draws = 100){
  # Summarise replications as densities
  #
  # Args:
  # fit: stanfit object
  # data_stan: data input to stan
  # draws: number of replications to show
  #
  # Returns:
  # Dataframe of replications
  
  t <- c(0, data_stan$t_rep)
  s_name <- c("Saccharomyces", "Schixosachararomyces")
  
  extract_draws <- function(tmp, s , cond){
    do.call(rbind,
            lapply(1:2,
                   function(ex){
                     out <- melt(tmp[sample(1:nrow(tmp), draws), , ex], varnames = c("Draw", "Age"), value.name = "Volume")
                     out$Age <- t[out$Age]
                     out$Species <- s
                     out$Condition <- cond
                     out$Experiment <- c("A", "B")[ex]
                     return(out)
                   }))
  }
  
  yrep <- extract(fit, pars = c("y1_rep", "y2_rep"))
  rep_single <- lapply(1:2,
                       function(i){extract_draws(yrep[[i]], s_name[i], "Single")})
  
  yrep <- extract(fit, pars = c("y12_rep"))[[1]]
  rep_mixed <- lapply(1:2, function(i){extract_draws(yrep[, , i, ], s_name[i], "Mixed")})
  
  do.call(rbind, c(rep_single, rep_mixed))
}

plot_PPC <- function(rep, df){
  # Plot PPC
  #
  # Args:
  # rep: replications dataframe (either as density or draws)
  # df1: dataframe of species 1
  # df2: dataframe of species 2
  # df12: dataframe of mixed experiment
  #
  # Returns:
  # List of ggplots
  
  library(ggplot2)
  library(cowplot)
  
  # If else statement for density or spaghetti
  if ("Density" %in% colnames(rep)){
    palette <- c("#FFFFFF", "#DEEBF7", "#9ECAE1", "#3182BD", "#000000") # blue
    
    plot_fun <- function(rep, df, s, ex){
      df <- subset(df, Species == s)
      rep <- subset(rep, Species == s & Experiment == ex)
      rep <- subset(rep, Age < 10 * ceiling(max(df$Age) / 10))
      df <- subset(df, Experiment == ex)
      
      ggplot() +
        geom_tile(data = rep, aes(x = Age, y = Volume, fill = Density)) +
        scale_fill_gradientn(colours = palette) + #, trans = "sqrt") + # problem multiplicate error lead to big density and log is horrible
        geom_point(data = df, aes(x = Age, y = Volume), colour = "#D55E00", size = 2) +
        geom_line(data = df, aes(x = Age, y = Volume), colour = "#D55E00", size = 1) +
        scale_x_continuous(limits = c(0, NA), expand = c(0.01, 0)) +
        scale_y_continuous(limits = c(0, NA), expand = c(0.01, 0)) +
        labs(subtitle = paste(s, " (", ex, ")", sep = "")) +
        theme_classic(base_size = 20) + theme(legend.position = "bottom")
    }
  }else{
    plot_fun <- function(rep, df, s, ex){
      df <- subset(df, Species == s)
      rep <- subset(rep, Species == s & Experiment == ex)
      rep <- subset(rep, Age < 10 * ceiling(max(df$Age) / 10))
      df <- subset(df, Experiment == ex)
      
      ggplot() +
        geom_point(data = rep, aes(x = Age, y = Volume, group = Draw), size = 1 , colour = "grey", alpha = 0.1) + # point or line
        geom_smooth(data =  aggregate(Volume ~ Age, rep, mean), aes(x = Age, Volume), size = 2, colour = "blue", method = "loess", span = 0.2) +
        geom_point(data = df, aes(x = Age, y = Volume), size = 2) +
        geom_line(data = df, aes(x = Age, y = Volume), size = 1) +
        scale_x_continuous(limits = c(0, NA), expand = c(0.01, 0)) +
        scale_y_continuous(limits = c(0, NA), expand = c(0.01, 0)) +
        labs(subtitle = paste(s, " (", ex, ")", sep = "")) +
        theme_classic(base_size = 20) + theme(legend.position = "bottom")
    }
  }
  
  tmp <- expand.grid(Species = c("Saccharomyces", "Schixosachararomyces"), Experiment = c("A", "B"))
  l <- lapply(c("Single", "Mixed"),
              function(cond){
                pl <- lapply(1:nrow(tmp),
                             function(x){
                               plot_fun(subset(rep, Condition == cond),
                                        subset(df, Condition == cond),
                                        s = tmp$Species[x], ex = tmp$Experiment[x])
                             })
                plot_grid(plotlist = pl, nrow = 2, ncol = 2)
              })
  names(l) <- c("Single", "Mixed")
  
  return(l)
}
