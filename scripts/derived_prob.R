
# source("~/00-mount.R")

source("scripts/setup.R")
source("scripts/write_nc.R")
# source("scripts/functions_derived.R")

"~/bucket_mine/results/global_spei_ww/new_derived" -> dir_derived
# dir.create(dir_derived)

read_delim("~/bucket_mine/misc_data/CMIP5_model_temp_thresholds.csv") -> thresholds

thresholds %>% 
  select(1:6) %>% 
  pivot_longer(-Model, names_to = "Warm") %>% 
  mutate(Warm = str_sub(Warm, 3)) -> thresholds

thresholds %>% 
  mutate(Warm = ifelse(str_length(Warm) == 1, str_glue("{Warm}.0"), Warm)) -> thresholds



# ********************

for(dom in (c("AFR", "AUS", "CAM", "CAS", "EAS", "EUR", "NAM", "SAM", "SEA", "WAS"))){      # ***************
  
  print(str_glue("********** PROCESSING DOM {dom} **********"))
  
  str_glue("~/bucket_mine/results/global_spei_ww/new_raw/{dom}") %>%
    list.files(full.names = T) %>%
    str_split("_", simplify = T) %>%
    {str_glue("{.[,7]}_{.[,8]}")} %>%
    unique() -> mods
  
  for(mod in mods){                                                                    # **********************
    
    plan(multicore, gc = T)
    
    print(str_glue("PROCESSING MODEL {mod}"))
    
    for(acc in c("03", "06", "12")){   
      
      print(str_glue("   PROCESSING ACC {acc}"))
      
      # plan(sequential)
      # plan(multicore, workers = 10)
      
      print(str_glue("      Importing into stars"))
      tic("      --Done")
      
      f <- character()
      while(length(f) == 0){
        
        str_glue("~/bucket_mine/results/global_spei_ww/new_raw/{dom}") %>% 
          list.files(full.names = T) %>%
          .[str_detect(., mod)] %>% 
          .[str_detect(., str_glue("spei-{acc}"))] -> f
        
      }
      
      f %>% 
        future_map(read_ncdf, .options = furrr_options(seed = NULL)) %>% 
        suppressMessages() %>% 
        do.call(c, .) -> s
      
      toc()
      
      
      mod %>%
        str_split("_", simplify = T) %>% 
        .[,ncol(.)] %>% 
        str_split("-", simplify = T) -> mod_sh
      
      if(str_detect(mod, "MPI")){
        mod_sh %>% 
          {str_glue("{.[,3]}-{.[,4]}-{.[,5]}")} -> mod_sh
      } else if(str_detect(mod, "_MIROC-")){
        mod_sh %>% 
          .[,2] -> mod_sh
      } else {
        mod_sh %>% 
          {str_glue("{.[,(ncol(.)-1)]}-{.[,ncol(.)]}")} -> mod_sh
      }
      
      
      
      print(str_glue("      Processing 6 warming levels"))
      tic(str_glue("      --Done"))
      
      future_walk(c("0.5", "1.0", "1.5", "2.0", "2.5", "3.0"), function(warm){
        
        # print(str_glue("      Processing warming level {warm}"))
        
        if(warm == 0.5){
          c(1971, 2000) -> start_end
          
        } else {
          
          thresholds %>% 
            filter(Model == mod_sh) %>% 
            filter(Warm == warm) %>% 
            pull(value) %>% 
            {c(.-10, .+10)} -> start_end
          
        }
        
        s %>%
          filter(year(time) >= start_end[1],
                 year(time) <= start_end[2]) -> ss
        
        # PROBABILITY D3 OVER
        {
          fname <- str_glue("~/bucket_mine/results/global_spei_ww/new_derived/spei-{acc}_probabilityD3over_{dom}_{str_split(mod, '_', simplify = T)[,1]}_{mod_sh}_{warm}C.nc")

          ss %>%
            st_apply(c(1,2), function(x){

              if(all(is.na(x))){
                NA
              } else {
                mean(x <= -1.6, na.rm = T)
              }

            },
            FUTURE = F,
            .fname = "prob") -> s_f

          func_write_nc_notime(s_f,
                               fname)
        }
        
        
        
        # PROBABILITY D0 UNDER
        {
          fname <- str_glue("~/bucket_mine/results/global_spei_ww/new_derived/spei-{acc}_probabilityD0under_{dom}_{str_split(mod, '_', simplify = T)[,1]}_{mod_sh}_{warm}C.nc")
          
          ss %>%
            st_apply(c(1,2), function(x){
              
              if(all(is.na(x))){
                NA
              } else {
                mean(x > -0.5, na.rm = T)
              }
              
            },
            FUTURE = F,
            .fname = "prob") -> s_f
          
          func_write_nc_notime(s_f,
                               fname)
        }
        
        
        
        # MEAN
        {
          fname <- str_glue("~/bucket_mine/results/global_spei_ww/new_derived/spei-{acc}_mean_{dom}_{str_split(mod, '_', simplify = T)[,1]}_{mod_sh}_{warm}C.nc")

          ss %>%
            st_apply(c(1,2), function(x){

              if(all(is.na(x))){
                NA
              } else {
                mean(x, na.rm = T)
              }

            },
            FUTURE = F,
            .fname = "prob") -> s_f

          func_write_nc_notime(s_f,
                               fname)
        }
        
        
        
        # PERCENTILES
        {
          probs = c("0.05", "0.10", "0.50", "0.90", "0.95")

          ss %>%
            st_apply(c(1,2), function(x){

              if(all(is.na(x))){
                rep(NA, length(probs))
              } else {
                quantile(x, prob = as.numeric(probs), na.rm = T)
              }

            },
            FUTURE = F,
            .fname = "perc") %>%
            aperm(c(2,3,1)) -> s_f

          imap(probs, function(pp, i){

            s_f %>%
              slice(perc, i) -> ss_f

            fname <- str_glue("~/bucket_mine/results/global_spei_ww/new_derived/spei-{acc}_{str_sub(pp,3)}perc_{dom}_{str_split(mod, '_', simplify = T)[,1]}_{mod_sh}_{warm}C.nc")

            func_write_nc_notime(ss_f,
                                 fname)

          })
        }
        
        
      }) # end of WL loop
      toc()
      
      
    } # end of acc loop
    
    plan(sequential)
    gc()
    
  } # end of mod loop
  
} # end of dom loop




