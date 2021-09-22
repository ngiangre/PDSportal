library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)
library(RSQLite)
library(data.table)
library(pool)
library(cowplot)
library(shinybusy)
library(shinyjqui)
library(shinyjs)

con <- dbPool(
    drv = RSQLite::SQLite(),
    dbname = "data/effect_peds_19q2_v0.1_20210831.sqlite"
)
onStop(function() poolClose(con))

tmp <- 
  tbl(con,"ade_raw") %>% 
  select(safetyreportid,nichd,sex) %>% 
  distinct()

total_stage_reports <- 
  tmp %>% 
  group_by(nichd) %>% 
  count(name = "total") %>% 
  collect() %>% 
  data.table()
total_stage_sex_reports <- 
  tmp %>% 
  group_by(nichd,sex) %>% 
  count(name = "total") %>% 
  collect() %>% 
  data.table()

min_date <- 
  tbl(con,"ade_raw") %>% 
  select(receive_date) %>% 
  arrange(receive_date) %>% 
  head(1) %>% 
  collect() %>% unlist %>% unname
max_date <- 
  tbl(con,"ade_raw") %>% 
  select(receive_date) %>% 
  arrange(desc(receive_date)) %>% 
  head(1) %>% 
  collect() %>% unlist %>% unname

drug_table <-
  tbl(con,"drug") %>%
  collect() %>%
  data.table() %>%
  na.omit() %>%
  .[,.(atc_concept_id,N = ndrugreports,
       code = atc_concept_code,
       ATC5=atc_concept_name,
       ATC4=atc4_concept_name,
       ATC3=atc3_concept_name,
       ATC2=atc2_concept_name,
       ATC1=atc1_concept_name)] %>% 
  .[order(N,decreasing = T)]

drugNames <- 
  drug_table[order(N,decreasing = T),.(atc_concept_id,code,N,ATC5)][,paste0(ATC5," [",code,"] (N=",scales::comma(N,accuracy = 1),")")]
drugIDs <- 
  drug_table[order(N,decreasing = T),.(atc_concept_id,ATC5)][,atc_concept_id]
names(drugIDs) <- drugNames

event_table <-
  tbl(con,"event") %>%
  collect() %>%
  data.table() %>%
  na.omit() %>%
  .[,.(meddra_concept_id,N = neventreports,
       code = meddra_concept_code_1,
       PT=meddra_concept_name_1,
       HLT=meddra_concept_name_2,
       HLGT=meddra_concept_name_3,
       SOC=meddra_concept_name_4)]

null_dist <- 
  tbl(con,"ade_null_distribution") %>% 
  collect() %>% 
  data.table()

null_dist_summary <- 
  null_dist %>% 
  .[,
    .(
      null_lwr = quantile(gam_score,c(0.01)),
      null_mean = mean(gam_score),
      null_upr = quantile(gam_score,c(0.99))
    ),
    nichd
  ]

theme_big <- 
  theme_classic(base_size=16) + 
  theme(
    strip.text = element_text(color="black",face="bold"),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(color="black",face="bold",
                               angle=45,vjust=1,hjust=1),
    axis.text.y = element_blank(),
    axis.title.y = element_text(color="black",face="bold"),
    legend.position = "top",
    legend.text = element_text(color="black",face="bold"),
    legend.key.size = unit(0.4,"cm"),
    legend.box.margin = margin(-0.4,-0.4,0,-0.4,unit = "cm")
  )

theme_small <- 
  theme_classic(base_size=12) + 
  theme(
    strip.text = element_text(color="black",face="bold",size=8),
    axis.title.x = element_text(color="black",face="bold",size=8),
    axis.title.y = element_text(color="black",face="bold",size=8),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1,color="black",face="bold",size=6),
    axis.text.y = element_text(color="black",face="bold",size=8),
    legend.text = element_text(color="black",face="bold",size=8),
    legend.key.size = unit(0.4,"cm"),
    legend.box.margin = margin(-0.4,-0.4,-0.4,-0.4,unit = "cm")
  )

stages <-
  c("term_neonatal","infancy",
    "toddler","early_childhood",
    "middle_childhood","early_adolescence",
    "late_adolescence")
stages_split <- str_replace(stages,"_","\n")

ade_cohort <- function(drugs=c(),events=c(),database=con){
  
  if(length(events)==0 | length(drugs)==0){errorCondition("no drugs or events given")}
  
  tmp <-
    expand.grid(
      atc_concept_id = drugs,
      meddra_concept_id = events
    ) %>%
    merge(
      tbl(con,"ade_nichd") %>%
        filter(
          atc_concept_id %in% drugs &
            meddra_concept_id %in% events) %>%
        collect() %>% data.table(),
      by=c("atc_concept_id","meddra_concept_id")
    ) %>%
    data.table() %>%
    merge(
      tbl(con,"ade") %>%
        filter(
          atc_concept_id %in% drugs &
            meddra_concept_id %in% events) %>%
        collect() %>% data.table(),
      by=c("ade","atc_concept_id","meddra_concept_id")
    ) %>%
    merge(
      tbl(con,"ade_null") %>% 
        collect() %>% 
        data.table(),
      by="nichd"
    ) %>% 
    merge(
      tbl(con,"event") %>%
        filter(
          meddra_concept_id %in% events
        ) %>%
        select(meddra_concept_id,meddra_concept_name_1) %>%
        rename(meddra_concept_name = meddra_concept_name_1) %>%
        collect() %>% data.table() %>% na.omit() %>% unique(),
      by=c("meddra_concept_id")
    ) %>%
    merge(
      tbl(con,"drug") %>%
        filter(
          atc_concept_id %in% drugs
        ) %>%
        select(atc_concept_id,atc_concept_name) %>%
        collect() %>% data.table() %>% na.omit() %>% unique(),
      by=c("atc_concept_id")
    )
  
  tmp$significance <-
    ifelse(tmp$gt_null_statistic==1,"Nominal","NA")
  tmp[gt_null_99==1,"significance"]="Null model"
  
  tmp
  
}

plot_ade_risks_null_shade <- function(x,color="red",theme=theme_big,ts=5){
  x$nichd_split <- str_replace(x$nichd,"_","\n")
  x$NICHD = factor(x$nichd_split,levels=str_replace(stages,"_","\n"))
  x %>% 
    merge(
      null_dist_summary,
      by="nichd"
    ) %>%
    ggplot(aes(factor(nichd_split,levels=stages_split),gam_score,group=ade)) +
    geom_ribbon(aes(ymin=null_lwr,ymax=null_upr),
                fill="lightgray",alpha=0.3) +
    geom_bar(stat="identity",aes(y=log10(as.numeric(DE)+1)),fill="lightgray",color="black") +
    ggrepel::geom_label_repel(
      aes(
        label=DE,
        y=log10(as.numeric(DE)+1)
      ),
      vjust=1.3,fontface="bold",size=ts) +
    geom_point(color=color,size=2) +
    geom_errorbar(aes(ymin=gam_score_90mse,ymax=gam_score_90pse),
                  width=0.1,color=color,size=1) +
    geom_hline(yintercept = 0,color="black",size=0.5) +
    facet_wrap(~ade_name,labeller = label_wrap_gen(width = 20)) +
    xlab("") +
    ylab("Risk of ADE (GAM log odds)") + 
    theme +
    theme(
      axis.line = element_line(),
      axis.ticks = element_line(),
      axis.text.y = element_text(color="black",face="bold",size=12)
    )
  
}

population_summary_stage <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
  pop <- 
    tbl(con,"ade_raw") %>% 
    filter(
      atc_concept_id %in% drugs &
        meddra_concept_id %in% events
    ) %>% 
    collect() %>% 
    data.table()
  pop$nichd_split <- str_replace(pop$nichd,"_","\n")
  
  tmp <- 
    merge(
      total_stage_reports,
      pop[,
          .(safetyreportid,nichd)
      ] %>% 
        unique() %>% 
        .[,.N,.(nichd)],
      all.x=T
    ) %>%
    .[,.(nichd,nichd_split = str_replace(nichd,"_","\n"),N = ifelse(is.na(N),0,N))]
  
  merge(
    tmp,
    total_stage_reports,
    by="nichd",
    all.y = T
  ) %>% 
    .[,.(nichd_split,N,total,
         prop = (N/total),
         prop_norm = (((N/total) - min((N/total)))/(max((N/total)) - min((N/total)))),
         prop_norm_mm = (((N/total) - min((N/total)))/(max((N/total)) - min((N/total))))*max(N)
    )
    ] %>% 
    ggplot(aes(factor(nichd_split,levels=stages_split),N)) +
    geom_bar(stat="identity",color="black",fill='lightgray') +
    geom_point(aes(y=prop_norm_mm),color="yellowgreen",pch=21) +
    geom_line(aes(y=prop_norm_mm,group=1),color="yellowgreen",linetype="dotted") +
    ggrepel::geom_label_repel(aes(label=paste0(round(prop*100,2),"% of ",scales::comma(total)),y=prop_norm_mm),
                              color="yellowgreen",size=ts,fontface="bold",
                              fill = "white") +
    xlab("") +
    coord_cartesian(clip="off") +
    scale_y_continuous(labels=scales::comma) +
    geom_text(aes(label=N,y=N),vjust=-1,fontface="bold",size=ts) +
    ylab("Number of reports") +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    theme
}

population_summary_sex <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
  pop <- 
    tbl(con,"ade_raw") %>% 
    filter(
      atc_concept_id %in% drugs &
        meddra_concept_id %in% events
    ) %>% 
    collect() %>% 
    data.table()
  pop$nichd_split <- str_replace(pop$nichd,"_","\n")
  
  merge(
    total_stage_sex_reports,
    pop[,
        .(safetyreportid,nichd,sex)
    ] %>% 
      unique() %>% 
      .[,.N,.(sex,nichd)],
    all.x=T
  ) %>%
    .[,.(sex,nichd,nichd_split = str_replace(nichd,"_","\n"),N = ifelse(is.na(N),0,N))] %>%
    ggplot(aes(factor(nichd_split,levels=stages_split),N,fill=sex)) +
    geom_bar(stat="identity",position = position_dodge(width=0.9),color="black") +
    xlab("") +
    coord_cartesian(clip="off") +
    scale_y_continuous(labels=scales::percent) +
    geom_text(aes(label=N,y=N),vjust=-1,position = position_dodge(width=0.9),fontface="bold",size=ts) +
    colorspace::scale_fill_discrete_qualitative(breaks=c("Female","Male"),labels=c("Female","Male")) +
    guides(fill=guide_legend(title = NULL)) +
    ylab("Number of reports") +
    theme
}

population_summary_poly <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
  pop <- 
    tbl(con,"ade_raw") %>% 
    filter(
      atc_concept_id %in% drugs &
        meddra_concept_id %in% events
    ) %>% 
    collect() %>% 
    data.table()
  pop$nichd_split <- str_replace(pop$nichd,"_","\n")
  
  pop[,
      .(safetyreportid,polypharmacy,nichd_split)
  ] %>% 
    merge(
      data.table(nichd_split = stages_split),
      all.y=T
    ) %>%
    ggplot(aes(factor(nichd_split,levels=stages_split),polypharmacy)) +
    ggbeeswarm::geom_quasirandom(groupOnX = T) +
    xlab("") +
    ylab("Number of drugs") +
    theme +
    theme(
      axis.text.y = element_text(color="black",face="bold"),
      axis.line = element_line()
    )
}

population_summary_reporter <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
  pop <- 
    tbl(con,"ade_raw") %>% 
    filter(
      atc_concept_id %in% drugs &
        meddra_concept_id %in% events
    ) %>% 
    collect() %>% 
    data.table()
  pop$nichd_split <- str_replace(pop$nichd,"_","\n")
  pop$reporter_qualification <- str_replace(pop$reporter_qualification," ","\n")
  
  merge(
    expand.grid(
      nichd_split = stages_split,
      reporter_qualification = pop[,unique(reporter_qualification)]
    ) %>% data.table(),
    pop[
      ,.(safetyreportid,reporter_qualification,nichd_split)
    ] %>% 
      unique() %>% 
      .[,
        .(N = .N),
        .(nichd_split,reporter_qualification)
      ],
    all.x=T
  ) %>%
    .[,.(nichd_split,reporter_qualification,N = ifelse(is.na(N),0,N))] %>% 
    ggplot(aes(factor(nichd_split,levels=stages_split),N,fill=reporter_qualification)) +
    geom_bar(stat="identity",position = position_dodge(width=0.9),color="black") +
    xlab("") +
    scale_y_continuous(labels=scales::percent) +
    coord_cartesian(clip="off") +
    geom_text(aes(label=N,y=N),vjust=-1,position = position_dodge(width=0.9),fontface="bold",size=ts) +
    colorspace::scale_fill_discrete_qualitative() +
    guides(fill=guide_legend(title=NULL,ncol=2)) +
    ylab("Number of reports") +
    theme
}

population_summary_atc <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
  pop <- 
    tbl(con,"ade_raw") %>% 
    filter(
      atc_concept_id %in% drugs &
        meddra_concept_id %in% events
    ) %>% 
    collect() %>% 
    data.table()
  pop$nichd_split <- str_replace(pop$nichd,"_","\n")
  reports=pop[,unique(safetyreportid)]
  
  atc_names <- 
    tbl(con,"atc_raw_map") %>% 
    collect() %>% 
    data.table() %>% 
    .[,unique(atc1_concept_name[order(atc1_concept_name)])] %>% 
    str_replace_all(" ","\n")
  tmp <- 
    merge(
      expand.grid(
        nichd_split = stages_split,
        atc1_concept_name = atc_names
      ) %>% data.table() %>% 
        merge(
          tbl(con,"drug") %>% 
            select(atc1_concept_name,atc1_concept_code) %>% 
            distinct() %>% 
            collect() %>% 
            data.table() %>% 
            .[,.(code = atc1_concept_code,
                 atc1_concept_name = str_replace_all(atc1_concept_name," ","\n"))],
          by="atc1_concept_name"
        ),
      merge(
        tbl(con,"ade_raw") %>% 
          filter(safetyreportid %in% reports &
                   !(atc_concept_id %in% drugs)) %>% 
          select(safetyreportid,atc_concept_id,nichd) %>% 
          collect() %>% 
          data.table() %>% 
          .[,.(safetyreportid,atc_concept_id,nichd_split = str_replace_all(nichd,"_","\n"))] %>% 
          unique(),
        tbl(con,"drug") %>% 
          collect() %>% 
          data.table() %>% 
          .[,.(atc_concept_id,code = atc1_concept_code,
               atc1_concept_name = str_replace_all(atc1_concept_name," ","\n"))],
        by="atc_concept_id"
      ) %>% 
        .[,.(safetyreportid,code,atc1_concept_name,code,nichd_split)] %>% 
        unique() %>% 
        .[,.N,
          .(atc1_concept_name,
            nichd_split,code)
        ] %>% na.omit(),
      all.x=T,
      by=c("atc1_concept_name","nichd_split","code")
    )
  tmp %>% 
    .[atc1_concept_name!="VARIOUS",
      .(nichd_split,
        atc1_concept_name = paste0(atc1_concept_name,"\n(ATC 1st code: ",code,")"),
        N = ifelse(is.na(N),0,N))] %>% 
    ggplot(aes(factor(nichd_split,levels=stages_split),N)) +
    geom_bar(stat="identity",color="black",fill='lightgray') +
    xlab("") +
    coord_cartesian(clip="off") +
    scale_y_continuous(labels=scales::comma) +
    geom_text(aes(label=N,y=N),vjust=-1,fontface="bold",size=ts) +
    facet_wrap(~atc1_concept_name~.,labeller = label_wrap_gen(width=20)) +
    ylab("Number of reports") +
    theme +
    theme(
      strip.text = element_text(face="bold",size=12),
      axis.text.x = element_text(face="bold",size=10)
    )
}

population_summary_date <- function(drugs=c(),events=c(),database=con,ts=5,theme=theme_big){
  pop <- 
    tbl(con,"ade_raw") %>% 
    filter(
      atc_concept_id %in% drugs &
        meddra_concept_id %in% events
    ) %>% 
    collect() %>% 
    data.table()
  pop$nichd_split <- str_replace(pop$nichd,"_","\n")
  pop$receive_date <- as.Date(pop$receive_date)
  
  range_ <- seq.Date(as.Date(min_date),as.Date(max_date),by = 1)
  date_range <- data.table(receive_date = range_)
  date_range$interval <- cut(range_,breaks = "quarters")
  
  pop[,
      .(receive_date,nichd_split,safetyreportid)
  ] %>% 
    unique() %>% 
    merge(
      expand.grid(
        nichd_split = stages_split,
        receive_date = range_
      ) %>% data.table(),
      all.y=T
    ) %>%
    merge(
      date_range,
      by="receive_date",
      all.x=T
    ) %>% 
    .[,.(N = length(unique(na.omit(safetyreportid)))),.(interval,nichd_split)] %>% 
    .[,.(interval = as.Date(interval),N,nichd_split)] %>% 
    ggplot(aes(interval,N,group=1)) +
    geom_point() +
    geom_line() +
    scale_x_date(breaks="3 years",date_labels = "%Y") +
    scale_y_continuous(labels=scales::number_format(accuracy=1),lim=c(0,NA)) +
    facet_wrap(~factor(nichd_split,stages_split)) +
    xlab("Time") +
    ylab("Number\nof reports") +
    theme +
    theme(
      axis.text.y = element_text(color="black",face="bold"),
      axis.line = element_line()
    )
}

population_summary <- function(drugs=c(),events=c(),database=con,theme=theme_small){
  
  pop <- 
    tbl(con,"ade_raw") %>% 
    filter(
      atc_concept_id %in% drugs &
        meddra_concept_id %in% events
    ) %>% 
    collect() %>% 
    data.table()
  pop$receive_date <- as.Date(pop$receive_date)
  pop$nichd_split <- str_replace(pop$nichd,"_","\n")
  pop$reporter_qualification <- str_replace(pop$reporter_qualification," ","\n")
  reports=pop[,unique(safetyreportid)]
  bm <- 0.2
  ts <- 2.5
  ##stage
  tmp <- 
    merge(
      total_stage_reports,
      pop[,
          .(safetyreportid,nichd)
      ] %>% 
        unique() %>% 
        .[,.N,.(nichd)],
      all.x=T
    ) %>%
    .[,.(nichd,nichd_split = str_replace(nichd,"_","\n"),N = ifelse(is.na(N),0,N))]
  
  gstage <- 
    merge(
      tmp,
      total_stage_reports,
      by="nichd",
      all.y = T
    ) %>% 
    .[,.(nichd_split,N,total,
         prop = (N/total),
         prop_norm = (((N/total) - min((N/total)))/(max((N/total)) - min((N/total)))),
         prop_norm_mm = (((N/total) - min((N/total)))/(max((N/total)) - min((N/total))))*max(N)
    )
    ] %>% 
    ggplot(aes(factor(nichd_split,levels=stages_split),N)) +
    geom_bar(stat="identity",color="black",fill='lightgray') +
    geom_point(aes(y=prop_norm_mm),color="yellowgreen",pch=21) +
    geom_line(aes(y=prop_norm_mm,group=1),color="yellowgreen",linetype="dashed") +
    geom_line(aes(y=prop_norm_mm,group=1),color="yellowgreen",linetype="dotted") +
    ggrepel::geom_label_repel(aes(label=paste0(round(prop*100,2),"% of ",scales::comma(total)),y=prop_norm_mm),
                              color="yellowgreen",size=ts,fontface="bold",
                              fill = "white") +
    xlab("") +
    coord_cartesian(clip="off") +
    scale_y_continuous(labels=scales::comma) +
    geom_text(aes(label=N,y=N),vjust=-1,fontface="bold",size=ts) +
    ylab("Number\nof reports") +
    theme +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      plot.margin = unit(c(0.4,0,0,0), "cm")
    )
  ##Sex
  gsex <- 
    merge(
      total_stage_sex_reports,
      pop[,
          .(safetyreportid,nichd,sex)
      ] %>% 
        unique() %>% 
        .[,.N,.(sex,nichd)],
      all.x=T
    ) %>%
    .[,.(sex,nichd,nichd_split = str_replace(nichd,"_","\n"),N = ifelse(is.na(N),0,N))] %>%
    ggplot(aes(factor(nichd_split,levels=stages_split),N,fill=sex)) +
    geom_bar(stat="identity",position = position_dodge(width=0.9),color="black") +
    xlab("") +
    coord_cartesian(clip="off") +
    scale_y_continuous(labels=scales::percent) +
    geom_text(aes(label=N,y=N),vjust=-1,position = position_dodge(width=0.9),fontface="bold",size=ts) +
    colorspace::scale_fill_discrete_qualitative(breaks=c("Female","Male"),labels=c("Female","Male")) +
    guides(fill=guide_legend(title = NULL)) +
    ylab("Number\nof reports") +
    theme +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.box.margin = margin(-0.4,0.6,0,-0.4,unit = "cm")
    ) 
  ##Polypharmacy
  gpoly <- 
    pop[,
        .(safetyreportid,polypharmacy,nichd_split)
    ] %>% 
    merge(
      data.table(nichd_split = stages_split),
      all.y=T
    ) %>%
    ggplot(aes(factor(nichd_split,levels=stages_split),polypharmacy)) +
    ggbeeswarm::geom_quasirandom(groupOnX = T,size=0.5) +
    xlab("") +
    ylab("Polypharmacy") +
    theme +
    theme(
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      axis.line = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) 
  
  ##Reporter
  grep <- merge(
    expand.grid(
      nichd_split = stages_split,
      reporter_qualification = pop[,unique(reporter_qualification)]
    ) %>% data.table(),
    pop[
      ,.(safetyreportid,reporter_qualification,nichd_split)
    ] %>% 
      unique() %>% 
      .[,
        .(N = .N),
        .(nichd_split,reporter_qualification)
      ],
    all.x=T
  ) %>%
    .[,.(nichd_split,reporter_qualification,N = ifelse(is.na(N),0,N))] %>% 
    ggplot(aes(factor(nichd_split,levels=stages_split),N,fill=reporter_qualification)) +
    geom_bar(stat="identity",position = position_dodge(width=0.9),color="black") +
    xlab("") +
    scale_y_continuous(labels=scales::percent) +
    coord_cartesian(clip="off") +
    geom_text(aes(label=N,y=N),vjust=-1,position = position_dodge(width=0.9),fontface="bold",size=ts) +
    colorspace::scale_fill_discrete_qualitative() +
    guides(fill=guide_legend(title=NULL,ncol=2)) +
    ylab("Number\nof reports") +
    theme +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.box.margin = margin(-0.4,0.6,0,-0.4,unit = "cm")
    )
  
  ##ATC class proportion
  atc_names <- 
    tbl(con,"atc_raw_map") %>% 
    collect() %>% 
    data.table() %>% 
    .[,unique(atc1_concept_name[order(atc1_concept_name)])] %>% 
    str_replace_all(" ","\n")
  tmp <- 
    merge(
      expand.grid(
        nichd_split = stages_split,
        atc1_concept_name = atc_names
      ) %>% data.table() %>% 
        merge(
          tbl(con,"drug") %>% 
            select(atc1_concept_name,atc1_concept_code) %>% 
            distinct() %>% 
            collect() %>% 
            data.table() %>% 
            .[,.(code = atc1_concept_code,
                 atc1_concept_name = str_replace_all(atc1_concept_name," ","\n"))],
          by="atc1_concept_name"
        ),
      merge(
        tbl(con,"ade_raw") %>% 
          filter(safetyreportid %in% reports &
                   !(atc_concept_id %in% drugs)) %>% 
          select(safetyreportid,atc_concept_id,nichd) %>% 
          collect() %>% 
          data.table() %>% 
          .[,.(safetyreportid,atc_concept_id,nichd_split = str_replace_all(nichd,"_","\n"))] %>% 
          unique(),
        tbl(con,"drug") %>% 
          collect() %>% 
          data.table() %>% 
          .[,.(atc_concept_id,code = atc1_concept_code,
               atc1_concept_name = str_replace_all(atc1_concept_name," ","\n"))],
        by="atc_concept_id"
      ) %>% 
        .[,.(safetyreportid,code,atc1_concept_name,code,nichd_split)] %>% 
        unique() %>% 
        .[,.N,
          .(atc1_concept_name,
            nichd_split,code)
        ] %>% na.omit(),
      all.x=T,
      by=c("atc1_concept_name","nichd_split","code")
    )
  gatc <- tmp %>% 
    .[atc1_concept_name!="VARIOUS",
      .(nichd_split,
        atc1_concept_name,
        N = ifelse(is.na(N),0,N))] %>% 
    ggplot(aes(factor(nichd_split,levels=stages_split),N)) +
    geom_bar(stat="identity",color="black",fill='lightgray') +
    xlab("") +
    coord_cartesian(clip="off") +
    scale_y_continuous(labels=scales::comma) +
    geom_text(aes(label=N,y=N),vjust=-1,fontface="bold",size=ts) +
    facet_wrap(~atc1_concept_name~.,labeller = label_wrap_gen(width=20)) +
    ylab("Number of reports") +
    theme +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face="bold",size=12),
      axis.text.x = element_text(face="bold",size=10)
    )
  
  #Date
  range_ <- seq.Date(as.Date(min_date),as.Date(max_date),by = 1)
  date_range <- data.table(receive_date = range_)
  date_range$interval <- cut(range_,breaks = "quarters")
  
  gdate <- 
    pop[,
        .(receive_date,nichd_split,safetyreportid)
    ] %>% 
    unique() %>% 
    merge(
      expand.grid(
        nichd_split = stages_split,
        receive_date = range_
      ) %>% data.table(),
      all.y=T
    ) %>%
    merge(
      date_range,
      by="receive_date",
      all.x=T
    ) %>% 
    .[,.(N = length(unique(na.omit(safetyreportid)))),.(interval,nichd_split)] %>% 
    .[,.(interval = as.Date(interval),N,nichd_split)] %>% 
    ggplot(aes(interval,N,group=1)) +
    geom_point() +
    geom_line() +
    scale_x_date(breaks="3 years",date_labels = "%Y") +
    scale_y_continuous(labels=scales::number_format(accuracy=1),lim=c(0,NA)) +
    facet_wrap(~factor(nichd_split,stages_split)) +
    xlab("Time") +
    ylab("Number\nof reports") +
    theme +
    theme(
      strip.background = element_blank(),
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      axis.line = element_blank(),
      legend.position = "none"
    )
  
  #https://wilkelab.org/cowplot/articles/plot_grid.html
  plot_grid(
    plot_grid(
      plot_grid(gstage,gsex,gpoly,ncol=1,nrow=3),
      plot_grid(gatc),
      nrow=1,ncol=2,
      rel_widths = c(2,4.5)
    ),
    plot_grid(grep,gdate,nrow=1,ncol=2,rel_widths = c(2.5,2)),
    nrow=2,ncol=1,
    rel_heights = c(4.5,2)
  )
  
}

get_drug_name <- function(x){
  drug_table %>%
    .[atc_concept_id %in% x,
      unique(ATC5)]
}

get_event_name <- function(x){
  event_table %>%
    .[meddra_concept_id %in% x,
      unique(PT)]
}


body <-
  dashboardBody(
    useShinyjs(),
    fluidRow(
      column(#start left side
        width=3,
        add_busy_spinner(spin="self-building-square",position="top-left",color="red"),
        box(
          id="select-container",
          height="40%",
          width=NULL,
          title="Quantify population-level evidence for an adverse drug event (ADE) in the pediatric population",
          solidHeader = T,
          status = "primary",
          shinyjs::hidden(
            div(
              id="welcome",
              h3(strong("Start selecting ADEs"))
            )
          ),
          h5(strong("Select a drug(s):")),
          h5(em("Drug name [ATC 5th code] Number of total drug reports")),
          selectInput("drugInput", 
                      label=NULL,
                      choices=drugIDs,
                      multiple=T
          ),
          shinyjs::hidden(
            div(
              id="eventInput_div",
              h5(strong("Select an adverse event(s):")),
              h5(em("Event name [MedDRA SOC name] Number of drug reports of event")),
              selectizeInput("eventInput",
                             label=NULL,
                             choices=c(),
                             multiple=T
              )
            )
          ),
          shinyjs::hidden( 
            div(id="calculate-button",
                style="display: flex; justify-content: center;",
                actionButton("goButton", 
                             HTML("Get population risk<br>of all ADE pairs"),
                             style="font-weight: bold; word-wrap: break-word;"
                )
            )
          )
        ),
        box(
          id="nichd-to-age-container",
          title="Child development stages",
          solidHeader = T,
          height="30%",
          width=NULL,
          collapsible = T,
          collapsed = T,
          status="primary",
          div(
            h4(em("Defined by NICHD "),tags$a(href="https://doi.org/10.1542/peds.2012-0055I",target="_blank","here")),
            h5(em("Term neonatal : Birth to 1 month old")),
            h5(em("Infancy : 1 month to 1 year old")),
            h5(em("Toddler: 1 year to 2 years old")),
            h5(em("Early childhood: 2 years to 6 years old")),
            h5(em("Middle childhood: 6 years old to 12 years old")),
            h5(em("Early adolescence: 12 years old to 18 years old")),
            h5(em("Late adolescence: 18 years old to 21 years old"))
          )
        ),
        box(
          id="tips-container",
          title = "Tips when using the app",
          height="30%",
          width = NULL,
          solidHeader = T,
          collapsible = T,
          collapsed = T,
          status = "primary",
          div(
            tags$ul(
              tags$li(
                "The drug name provided should be the main ingredient. You can search for the drug name from the Brand name using",tags$a(href="https://athena.ohdsi.org/search-terms/terms?query",target="_blank","OHDSI's Athena")
              ),
              tags$li(
                "Each time you add or remove drugs and events, press the calculate button to update the plots."
              ),
              tags$li(
                "You can save the plots by right clicking on the images."
              ),
              tags$li(
                "The estimated ADE risks should be interpreted in light of the ADE report demographics."
              ),
              tags$li(
                "The ADE reporting information are all drug reports contained in the FDA's Adverse Event Reporting System up through 2019Q2. See the download tab for the pediatric-specific reporting we call Pediatric FAERS."
              )
            )
          )
        ),
        box(
          id="cite-container",
          height="30%",
          width = NULL,
          solidHeader=T,
          collapsible=T,
          collapsed = T,
          status="primary",
          title="Please cite us!",
          p(
            "Giangreco, Nicholas P. and Tatonetti, Nicholas P., A Database of Pediatric Drug Effects to Evaluate Ontogenic Mechanisms From Child Growth and Development. Available at SSRN: https://ssrn.com/abstract=3898786 or http://dx.doi.org/10.2139/ssrn.3898786"
          ),
          p(
            "Giangreco, N.P., Tatonetti, N.P. Evaluating risk detection methods to uncover ontogenic-mediated adverse drug effect mechanisms in children. BioData Mining 14, 34 (2021). https://doi.org/10.1186/s13040-021-00264-9."
          )
        ),
        box(
          id="disclaimer-container",
          height="30%",
          width = NULL,
          solidHeader=T,
          collapsible=T,
          collapsed = T,
          status="primary",
          title="Disclaimer",
          p(
            "This application is for research purposes only - please consult a healthcare professional for advice on appropriate pediatric drug treatments and potential adverse events."
          )
        )
      ),#end left side
      column(#start right side
        width=9,
        tabBox(
          id="container",
          title="",
          height="100%",
          width=NULL,
          tabPanel(
            "ADE report demographics",
            width=NULL,
            shinyjs::hidden(
              div(id="population_summary_text_div",
                  h4(strong("For children taking")),
                  em(textOutput("drug_text")),
                  h4(strong("we summarize the demographics for the reports of these adverse event(s):")),
                  em(textOutput("event_text")),
                  h4(strong("These data are used to estimate ADE risk at all child development stages, shown in the next tab."))
              )
            ),
            hr(),
            tabBox(
              id="population_summary_tabs",
              height=NULL,
              width=NULL,
              tabPanel(
                "Reporting across childhood",
                height=NULL,
                width=NULL,
                div(id="population_summary_stage_text_div",
                    h4(em("The percent of reports at each stage is shown in yellow-green. For example, 1 report of the adverse drug event(s) out of 100 total reports at that stage would give a reporting rate of 1%."))
                ),
                jqui_resizable(plotOutput("population_summary_stage",height="400px"))
              ),
              tabPanel(
                "By sex",
                height=NULL,
                width=NULL,
                jqui_resizable(plotOutput("population_summary_sex",height="300px"))
              ),
              tabPanel(
                "By reporter type",
                height=NULL,
                width=NULL,
                jqui_resizable(plotOutput("population_summary_reporter",height="300px"))
              ),
              tabPanel(
                "By polypharmacy",
                height=NULL,
                width=NULL,
                jqui_resizable(plotOutput("population_summary_poly",height="300px"))
              ),
              tabPanel(
                "By co-medications",
                height=NULL,
                width=NULL,
                h4(em("All other co-medication types")),
                jqui_resizable(plotOutput("population_summary_atc",height="600px"))
              ),
              tabPanel(
                "By reporting date",
                height=NULL,
                width=NULL,
                jqui_resizable(plotOutput("population_summary_date",height="400px"))
              ),
              tabPanel(
                "All graphs",
                height=NULL,
                width=NULL,
                jqui_resizable(plotOutput("population_summary",height="800px"))
              )
            )
          ),
          tabPanel(
            "ADE risk through childhood",
            width=NULL,
            shinyjs::hidden(
              div(id="ade_risk_text_div",
                  h4(em("We estimate ADE risk, and the 90% interval of risk variation, by sharing reporting information across all the child development stages.")),
                  h4(em("The height of the bars represent the number of adverse drug event reports at that stage reported to the FDA.")),
                  h4(em("The lightgray region indicates risk estimated from random reporting - risk estimates above this region indicate a non-random risk.")),
                  jqui_resizable(plotOutput("ade_risk_null_shade"))
              )
            ),
            DT::dataTableOutput("ade_risk_table")
          ),
          tabPanel(
            "About this app",
            width=NULL,
            div(
              h2("Overview",align = "center"),
              h4(
                "Drug side effects in children were mined from openFDA, the publicly available repository of the FDA. This is the first resource of 460,837 adverse drug events (ADEs) with estimated risk across child development stages. This R shinydashboard application is the publically accessible web application of our resource. Please see our publications on the left under 'Please cite us!'"
              ),
              h4(
                "The goal of the Pediatric Drug Safety portal (PDSportal) is to access our generated database to easily identify and evaluate evidence for adverse drug events in the pediatric population. See our publications in the sidebar for details."
              ),
              h4(
                "The intended usage is for summarizing the risk demographics and the estimated risk of an adverse drug event(s) occuring across child development stages. Note that ADE risk was estimated for drug-event pairs and not for a drug and multiple adverse events. In contrast, the risk demographics is shown in aggregate across drug-event pairs.
                "
              ),
              h4("Typical questions this app can evaluate are 'What adverse events have been reported with this drug?', 'What adverse events are not likely risks from a drug treatment?','What is the estimated risk of an ADE at and relative to other child development stages?', 'What is the significance of the estimated ADE risk at a particular stage?', and 'What are the population characteristics for an ADE risk?'"
              )
            ) 
          ), 
          tabPanel(
            "Download data",
            width=NULL,
            div(id="dump-div",
                tags$a(href='https://pds-database.s3.amazonaws.com/effect_peds_19q2_v0.1_20210831.sql.gz', 
                       target='_blank', 
                       h3('MySQL dump v0.1')
                ),
                tags$a(href='https://pds-database.s3.amazonaws.com/effect_peds_19q2_v0.2_20210920.sql.gz', 
                       target='_blank', 
                       h3('MySQL dump v0.2')
                ),
                tags$a(href='https://pds-database.s3.amazonaws.com/effect_peds_19q2_v0.1_20210831.sqlite.gz', 
                       target='_blank', 
                       h3('SQLite database v0.1')
                ),
                tags$a(href='https://pds-database.s3.amazonaws.com/effect_peds_19q2_v0.2_20210920.sqlite.gz', 
                       target='_blank', 
                       h3('SQLite database v0.2')
                ),
                h3("Select zipped CSV files:"),
                tags$a(href='https://pds-database.s3.amazonaws.com/database_generation_er_tables/dictionary.csv.gz', 
                       target='_blank', 
                       h5('Dictionary table (Descriptions of all fields from the 19 tables in the database)')
                ),
                tags$a(href='https://pds-database.s3.amazonaws.com/database_generation_er_tables/ade_raw.csv.gz', 
                       target='_blank', 
                       h5('Pediatric FAERS')
                ),
                tags$a(href='https://pds-database.s3.amazonaws.com/database_generation_er_tables/ade_nichd.csv.gz', 
                       target='_blank', 
                       h5('All 460,837 ADE risk estimates')
                ),
                tags$a(href='https://pds-database.s3.amazonaws.com/database_generation_er_tables/drug.csv.gz', 
                       target='_blank', 
                       h5('Drug vocabulary')
                ),
                tags$a(href='https://pds-database.s3.amazonaws.com/database_generation_er_tables/event.csv.gz', 
                       target='_blank', 
                       h5('Adverse event vocabulary')
                )
            )
          )
        )#end tabBox
      )#end right side
    )#end fluidRow
  )#end dashboardBody

sidebar <- dashboardSidebar(disable = T)

ui <-
  dashboardPage(
    header = dashboardHeader(
      title='PDSportal'
    ),
    sidebar = sidebar,
    body = body,
    skin="blue"
  )

server <- function(input,output,session) {
  
  #update event selection after drug selection
  observeEvent(input$drugInput,{
    
    isolate({
      #input <- list("drugInput" = 21603356)
      drugs_ <- input$drugInput
    })
    
    event_count <-
      tbl(con,"ade_raw") %>%
      filter(atc_concept_id %in% drugs_) %>%
      distinct(safetyreportid,meddra_concept_id) %>%
      group_by(meddra_concept_id) %>%
      count() %>%
      collect() %>%
      data.table()
    
    event_table_update <-
      merge(
        event_table[,.(meddra_concept_id,code,PT,SOC)],
        event_count[,.(meddra_concept_id,N=n)],
        by="meddra_concept_id"
      )
    
    eventNames <-
      event_table_update[
        order(N,decreasing = T),
        .(meddra_concept_id,code,N,PT,SOC)
      ][,
        paste0(stringr::str_replace_all(PT,' ','_'),
               " [",
               stringr::str_replace_all(SOC,' ','_'),
               "] (N =",
               scales::comma(N,accuracy = 1),
               ")")
      ]
    eventIDs <-
      event_table_update[order(N,decreasing = T),.(meddra_concept_id,code,N,PT)][,meddra_concept_id]
    names(eventIDs) <- eventNames
    
    
    updateSelectizeInput(session, "eventInput",
                         label = NULL,
                         choices = eventIDs
    )
    
  })
  
  #show event selection if drug(s) is selected
  observe({
    shinyjs::toggleElement(id="eventInput_div",
                           condition={!is.null(input$drugInput)})
  })
  
  #show calculate button after event selection
  observe({
    shinyjs::toggleElement(id="calculate-button",
                           condition={!is.null(input$eventInput)})
  })
  
  plots <- reactiveValues()
  text <- reactiveValues()
  
  #start plotting after pressing calculate button
  observeEvent(input$goButton,{
    
    show_modal_progress_line() # show the modal window
    
    drug_vec <- input$drugInput %>% as.integer()
    event_vec <- input$eventInput %>% as.integer()
    
    text$drugs <- 
      get_drug_name(drug_vec)
    
    text$events <- 
      get_event_name(event_vec)
    
    plots$ade_risk_null_shade <- 
      plot_ade_risks_null_shade(
        ade_cohort(
          drug_vec,
          event_vec
        )
      )
    
    update_modal_progress(0.2) # update progress bar value
    
    plots$population_summary_stage <- 
      population_summary_stage(
        drug_vec,
        event_vec
      )
    
    update_modal_progress(0.4) # update progress bar value
    
    plots$population_summary_sex <- 
      population_summary_sex(
        drug_vec,
        event_vec
      )
    
    plots$population_summary_poly <- 
      population_summary_poly(
        drug_vec,
        event_vec
      )
    
    plots$population_summary_reporter <- 
      population_summary_reporter(
        drug_vec,
        event_vec
      )
    
    update_modal_progress(0.6) # update progress bar value
    
    plots$population_summary_atc <- 
      population_summary_atc(
        drug_vec,
        event_vec
      )
    
    update_modal_progress(0.8) # update progress bar value
    
    plots$population_summary_date <- 
      population_summary_date(
        drug_vec,
        event_vec
      )
    
    update_modal_progress(1) # update progress bar value
    
    plots$population_summary <- 
      population_summary(
        drug_vec,
        event_vec
      )
    
    update_modal_progress(0.2) # update progress bar value
    remove_modal_progress()
  })
  
  #show population summary text after pressing calculate button
  observe({
    shinyjs::toggleElement(id="welcome",
                           condition=is.null(input$drugInput))
    shinyjs::toggleElement(id="population_summary_text_div",
                           condition=input$goButton>0)
    shinyjs::toggleElement(id="population_summary_tabs",
                           condition=input$goButton>0)
    shinyjs::toggleElement(id="population_summary_atc_text_div",
                           condition=input$goButton>0)
    shinyjs::toggleElement(id="population_summary_stage_text_div",
                           condition=input$goButton>0)
  })
  
  output$drug_text <- 
    renderText({
      ifelse(
        length(text$drugs)<2,
        paste0(text$drugs),
        ifelse(
          length(text$drugs)==2,
          paste0(text$drugs[1]," or ",text$drugs[2]),
          paste0(
            paste0(text$drugs[1:(length(text$drugs)-1)],collapse=", "),
            ", or ",text$drugs[length(text$drugs)]
          )
        )
      )
    })
  output$event_text <- 
    renderText({
      ifelse(
        length(text$events)<2,
        paste0(text$events),
        ifelse(
          length(text$events)==2,
          paste0(text$events[1]," or ",text$events[2]),
          paste0(
            paste0(text$events[1:(length(text$events)-1)],collapse=", "),
            ", or ",text$events[length(text$events)]
          )
        )
      )
    })
  
  #show ade risk text after pressing calculate button
  observe({
    shinyjs::toggleElement(id="ade_risk_text_div",condition=input$goButton>0)
  })
  
  output$ade_risk_null_shade <- renderPlot({plots$ade_risk_null_shade})
  output$population_summary_stage <- renderPlot({plots$population_summary_stage})
  output$population_summary_sex <- renderPlot({plots$population_summary_sex})
  output$population_summary_poly <- renderPlot({plots$population_summary_poly})
  output$population_summary_reporter <- renderPlot({plots$population_summary_reporter})
  output$population_summary_atc <- renderPlot({plots$population_summary_atc})
  output$population_summary_date <- renderPlot({plots$population_summary_date})
  output$population_summary <- renderPlot({plots$population_summary})
  
  output$ade_risk_table <- DT::renderDataTable({
    
    # Take a dependency on input$goButton
    if (input$goButton == 0)
      return()
    
    isolate({
      drug_vec <- input$drugInput %>% as.integer()
      event_vec <- input$eventInput %>% as.integer()
      ade_cohort(
        drug_vec,
        event_vec
      ) %>% 
        .[,.(Drug = atc_concept_name,
             Event = meddra_concept_name,
             Stage = str_replace(nichd,"_"," "),
             N = DE,
             `Lower` = gam_score_90mse %>% signif(digits = 3),
             `Log Odds` = gam_score %>% signif(digits = 3),
             `Upper` = gam_score_90pse %>% signif(digits = 3),
             `Random threshold` = null_99 %>% signif(digits = 3)
        )]
    })
  },
  server=F,
  extensions = c('Buttons', 'Scroller'),
  options = list(
    paging = TRUE,
    searching = TRUE,
    scrollY=T,
    fixedColumns = TRUE,
    autoWidth = TRUE,
    ordering = TRUE,
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel')
  ),
  class = "display"
  )
}

shinyApp(ui=ui,server=server)
