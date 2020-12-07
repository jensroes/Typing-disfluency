server <- function(input, output, session) {
  
  source("disfluency-app/environment.R")
  
  data <- reactive(d %>% filter(subj %in% as.numeric(input$N_example)))
 
  output$example <- renderDataTable({
    datatable(data(), 
              options = list(searching = F, pageLength = 5,
                             autoWidth = F), rownames = F) })
  
  
  # by participant density plot with normal distribution (can be removed)
  # allow to select ppts
  # which between log and normal scale
  
  d_example <- reactive({
    
    transform <- input$scaled
    grouping <- input$group
    
    if(!grouping){
      d_example <- data() %>% mutate(logIKI = log(IKI)) %>% group_by(subj) 
    }
    if(grouping){
      d_example <- data() %>% mutate(logIKI = log(IKI)) %>% group_by(subj, component) 
    }
    
    d_example <- d_example %>% 
      mutate(component = recode(component, LF = 'LF bigrams'),
             mean = round(mean(IKI),0),
             logmean = mean(logIKI),
             median = median(IKI),
             mode = dmode(IKI),
             sd = round(sd(IKI),0),
             logsd = sd(logIKI),
             logsd = sd(log(IKI)),
             n = n(),
             lognorm = dnorm(logIKI, logmean, logsd),
             norm = dnorm(IKI, mean, sd)) %>% ungroup() %>%
      mutate(maxnorm = max(norm),
             maxlognorm = max(lognorm),
             maxIKI = max(IKI),
             maxlogIKI = max(logIKI),
             subj = paste0("Participant id: ", subj))
    
  })
  
  d_summary <- reactive({
    transform <- input$scaled
    grouping <- input$group
    
    d_example <- d_example()
    
    if(!grouping){
      d_summary <- d_example %>% group_by(subj) %>% select(mean:n, maxnorm:maxlogIKI) %>% unique()
    }
    if(grouping){
      d_summary <- d_example %>% group_by(subj, component) %>% select(mean:n, maxnorm:maxlogIKI) %>% unique()
    }
    
    
  })
  
  output$density <- renderPlot({
    
    transform <- input$scaled
    grouping <- input$group
    
    d_example <- d_example()
    d_summary <- d_summary()
    
    if(!grouping){comp_group <- c("subj")}
    if(grouping){ comp_group <- c("subj","component")}
    
    
    grid <- with(d_example, seq(min(IKI), max(IKI), length = 100))
    normaldens <- plyr::ddply(d_example, comp_group, function(df) {
      data.frame( 
        IKI = grid,
        logIKI = log(grid),
        density = dnorm(grid, mean(df$IKI), sd(df$IKI)),
        logdensity = dnorm(log(grid), mean(log(df$IKI)), sd(log(df$IKI)))
      )
    })
    
    size = 1
    if(transform == "Normal"){
      x <- "IKIs [in msecs]"
      if(!grouping){
        plot <- ggplot(d_example, aes(x = IKI))
      }
      else if(grouping){
        plot <- ggplot(d_example, aes(x = IKI, colour = component))
      }
      plot <- plot + 
        geom_vline(data = filter(d_summary, subj == unique(subj)), aes(xintercept = mean), 
                  linetype = "dotted", colour = "grey30") +
        geom_line(aes(y = density, linetype = 'Normal'), size = size, data = normaldens) 
        
    }
    if(transform == "Log-scaled"){
      x <- "log-scaled IKIs [in msecs]"
      if(!grouping){
        plot <- ggplot(d_example, aes(x = logIKI)) 
      }
      else if(grouping){
        plot <- ggplot(d_example, aes(x = logIKI, colour = component, fill = component))
      }
      plot <- plot + 
        geom_vline(data = filter(d_summary, subj == unique(subj)), aes(xintercept = logmean), 
                  linetype = "dotted", colour = "grey30") +
        geom_line(aes(y = logdensity, linetype = 'Normal'), size = size, data = normaldens) 
    }
      
    plot <- plot + 
      geom_line(aes(y = ..density.., linetype = 'Empirical'), size = size, stat = 'density') +
      scale_linetype_manual("Density", values = c("dashed", "solid")) +
      labs(y = "Density", x = x,caption = "Vertical line shows empirical mean") +
      facet_wrap(~subj, ncol = 4, scales = "free_y") +
      theme(strip.text = element_text(hjust = 0)) +
      scale_colour_viridis_d("")

    plot

  }) #, height = 600, width = 800s
  
  
  mixedmodel <- reactive({
    d <- data() #%>% filter(component == "LF")
    d <- d %>% filter(component == "LF")
    
    #d %>% summarise(max = max(IKI))
    
    subs <- unique(d$subj)
    s_diff <- .01
    sigma <- sd(log(d$IKI))
    mu <- mean(log(d$IKI))
    mu_diff <- .1
    mixs <- map(subs, ~tryCatch({c(mixtools::normalmixEM(log(d[d$subj == .x,]$IKI), k = 2,
                                             lambda = c(.5,.5),
                                             mu = c(mu-mu_diff,mu+mu_diff), 
                                             sigma = c(sigma-s_diff,sigma+s_diff),
                                             fast=FALSE,
                                             maxit=10000,
                                             epsilon = 1e-16), subj = .x)}, error=function(e){} ))

#    mixs %>% enframe() %>% unlist()
    
    subjs_converged <- which(lengths(mixs)!=0)
    mus <- map(subjs_converged, ~mixs[[.x]]$mu)
    lambdas <- map(subjs_converged, ~mixs[[.x]]$lambda)
#    sigmas <- map(subjs_converged, ~mixs[[.x]]$sigma)

    mus <- do.call(rbind,lapply(mus,function(x) x[1:2])) %>%
      as.tibble() %>%
      rename(fluent_typing = V1,
             disfluent_typing = V2) %>%
      mutate_all(exp) 

    lambdas <- do.call(rbind,lapply(lambdas,function(x) x[2])) %>%
      as.tibble() %>%
      rename(proportion_of_disfluencies = V1) 

#    sigmas <- do.call(rbind,lapply(sigmas,function(x) x[1:2])) %>%
#      as.tibble() %>%
#      rename(typing = V1,
#             disfluency = V2) %>%
#      mutate(parameter = "sigma",
#             subj = subjs_converged)
    
    
    mixs <- bind_cols(mus, lambdas) %>%
      mutate(subj = subjs_converged)
  
    mixs  
    #table(mixs$parameter)
    
    #ggplot(mixs, aes(x = values, colour = writing, fill = writing)) +
    #  geom_histogram() +
    #  facet_wrap(~parameter, scales = "free")
    
        
    #  stat_function(aes(colour = 'Log-Normal'), fun = dlnorm, args = list(meanlog = d_summary$logmean, sdlog = d_summary$logsd)) 
    #   stat_function(fun = dnorm, 
    #                args = list(mean = my_mix[["mu"]][[1]], 
    #                            sd = my_mix[["sigma"]][[1]])) +
    #  stat_function(fun = dnorm, 
    #                args = list(mean = my_mix[["mu"]][[2]], 
    #                            sd = my_mix[["sigma"]][[2]]))
    
  })

  output$mog <- renderDataTable({
    datatable(mixedmodel(), # %>% group_by(subj, writing, ), 
              options = list(searching = F, pageLength = 5,
                             autoWidth = F), rownames = F) })
   
}



