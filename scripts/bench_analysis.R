library(ggplot2)

criterion_dir <- commandArgs(TRUE)[1]

d <- data.frame()
for (name in list.dirs(criterion_dir, F, F)) {
  d <- rbind(d, read.csv(file.path(criterion_dir, name, 'new', 'raw.csv')))
}

s <- strsplit(as.character(d$group), ' ')
cfg <- as.data.frame(t(as.data.frame(lapply(s, '[', 1:5))))
names(cfg) <- c('format', 'reader', 'seqlen', 'config', 'data_len')
cfg$data_len <- as.numeric(as.character(cfg$data_len))
cfg$seqlen <- as.numeric(as.character(cfg$seqlen))
cfg$cap_test <- grepl('_cap', as.character(cfg$reader))

d <- cbind(d, cfg)


d$gb_per_s <- (as.numeric(d$data_len) / 1e9) / (d$sample_time_nanos / 1e9) * d$iteration_count


reader_plot <- function(data, facets) {
  ggplot(data, aes(reader, gb_per_s, fill=reader)) +
    stat_summary(fun.y=mean, geom='bar', width=1, colour='#222222', size=0.2) + 
    stat_summary(fun.data=mean_se, geom = 'errorbar', width=0.2, alpha=0.5) +
    facet_grid(paste(paste(facets, collapse= "+"), " ~ ."),
               space='free_y', scale='free_y', switch='y') +
    expand_limits(x=1) +
    scale_fill_grey(start=0.05, end=0.95) +
    labs(x='reader', y='GB/s') +
    theme_bw() + theme(
      strip.text.y = element_text(angle=180),
      strip.background = element_rect(fill='white', colour='gray50'),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      #panel.border = element_rect(color=NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing = unit(-.5, 'pt')
    ) +
    coord_flip(expand=F, ylim=c(0, max(data$gb_per_s, na.rm=T)))
}

reader_cmp <- subset(d, !cap_test)
cap_cmp <- subset(d, cap_test)

outdir <- file.path(dirname(dirname(criterion_dir)), 'bench_results')
dir.create(outdir, F)

# comparison of readers
for (fmt in levels(reader_cmp$format)) {
  data <- subset(reader_cmp, format == fmt & !is.na(data_len))
  png(file.path(outdir, sprintf('bench_%s.png', fmt)), width=1400, height=30*length(unique(data$group))+150, res=200)
  print(reader_plot(data, c('seqlen', 'config')))
  dev.off()
  
  sub <- subset(data, seqlen == 500 & !grepl('(records|seq|iter)', config, perl=T))
  png(file.path(outdir, sprintf('bench_%s_simple.png', fmt)), width=1400, height=30*length(unique(sub$group))+150, res=200)
  print(reader_plot(sub, c('config')))
  dev.off()
}

# test different buffer capacities
lapply(split(cap_cmp, cap_cmp$format), function(data) {
  fmt <- data$format[1]
  png(file.path(outdir, sprintf('bench_cap_%s.png', fmt)), width=1400, height=1000, res=200)
  data$bufsize <- as.numeric(gsub('([0-9]+)ki', '\\1', data$config))
  print(ggplot(data, aes(bufsize, gb_per_s, color=as.factor(seqlen), group=seqlen)) +
          stat_summary(fun.y=mean, geom='point') +
          stat_summary(fun.y=mean, geom='line') +
          stat_summary(fun.data=mean_se, geom = 'errorbar', width=0.1, alpha=0.5) +
          scale_x_continuous(trans='log1p', breaks=2^(0:25)) +
          labs(x='Buffer size (Kib)', y='GB/s', color='Sequence length') +
          theme_bw())
  dev.off()
})
