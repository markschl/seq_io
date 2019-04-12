library(ggplot2)

criterion_dir <- commandArgs(T)[1]

#### read data ####

d <- data.frame()
for (subdir_name in list.dirs(criterion_dir, F, F)) {
  subdir = file.path(criterion_dir, subdir_name)
  for (bench_name in list.dirs(subdir, F, F)) {
    f <- file.path(subdir, bench_name, 'new', 'raw.csv')
    if (!file.exists(f))
      f = file.path(subdir, bench_name, '0', 'new', 'raw.csv')
    if (file.exists(f))
      d <- rbind(d, read.csv(f, stringsAsFactors=F))
  }
}

s <- strsplit(as.character(d$function.), ' ')
d$format = sapply(s, '[', 1)
d$name = sapply(s, '[', 2)
d$kind = sapply(s, '[', 3)
d$throughput_gb_s = with(d, throughput_num * 1024^3 / (sample_measured_value / iteration_count * 1e9))

d$kind[is.na(d$kind)] = 'borrow' # TODO: remove
d = subset(d, kind != 'discard')
d = split(d, ifelse(grepl('_cap_', d$group), 'cap', 'readers'))
d$readers$this_crate = grepl('seqio', d$readers$name)

#### setup output directory ####

outdir <- file.path(dirname(dirname(criterion_dir)), 'bench_results')
dir.create(outdir, F)

#### better labels ####

group_labels = c(
  fasta = 'FASTA',
  fasta_multiline = 'FASTA\n(lines wrapped)',
  fastq = 'FASTQ'
)

kind_labels = c(
  borrow = 'Iterate over records\nborrowing data',
  seq_iter = '… + iterate over\nsequence lines',
  seq = '… + access contigouous\nsequence',
  id = '… + access record ID', 
  owned = 'Iterate over / copy to\n owned (allocated) records',
  recordset = 'read_record_set() + iterate',
  parallel = 'Parallel processing\n(just iteration)'
)

name_labels = c(
  seqio = '(seq_io standard)',
  seqio_bytes = 'byte slice',
  seqio_str = 'str slice',
  seqio_multi = 'fastq::multiline::Reader',
  seqio_single = 'fasta::single_line::Reader',
  seqio_single_linestore = 'fasta::single_line::Reader (LineStore)',
  seqio_fastx = 'seq_io::fastx',
  seqio_fastx_dynamic = 'seq_io::fastx::dynamic',
  seqio_rset = 'per record set',
  seqio_record = 'per record',
  seqio_records = 'seq_io Reader::records()',
  seqio_clone_into = 'RefRecord::clone_into_owned()',
  seqio_seq_given = 'RefRecord::full_seq_given()',
  fastq_rs = 'fastq-rs',
  bio = 'Rust-Bio'
)

# for simple overview
kind_chosen = c('borrow', 'owned')
name_chosen = c('seqio', 'seqio_fastx', 'fastq_rs', 'bio')
d$readers$chosen = (d$readers$kind %in% kind_chosen) & (d$readers$name %in% name_chosen)

stopifnot(all(d$readers$group %in% names(group_labels)))
d$readers$group = factor(d$readers$group, names(group_labels))
levels(d$readers$group) = group_labels

stopifnot(all(d$readers$kind %in% names(kind_labels)))
d$readers$kind = factor(d$readers$kind, names(kind_labels))
levels(d$readers$kind) = kind_labels

stopifnot(all(d$readers$name %in% names(name_labels)))
d$readers$name = factor(d$readers$name, names(name_labels))
levels(d$readers$name) = name_labels

#### comparison of readers ####

reader_plot <- function(data, facets) {
  ggplot(data, aes(name, throughput_gb_s, fill=this_crate)) +
    stat_summary(fun=mean, geom='bar', width=0.8, colour='#222222', size=0.2) +
    stat_summary(fun.data=mean_se, geom = 'errorbar', width=0.2, alpha=0.5) +
    facet_grid(kind ~ group, space='free_y', scales='free_y') +
    coord_flip() +
    scale_x_discrete(limits=rev) +
    scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
    scale_fill_grey(start=0.5, end=0.9) +
    labs(y='Throughput [GiB/s]') +
    theme_bw() + theme(
      strip.text.y = element_text(angle=0, hjust=0),
      strip.background = element_rect(fill='white', colour='gray50'),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing.x = unit(2, 'pt'),
      panel.spacing.y = unit(1, 'pt'),
      legend.position = 'none'
    )
}

png(file.path(outdir, 'reader_comparison.png'), width=1400, height=700, res=150)
reader_plot(d$readers)
dev.off()

png(file.path(outdir, 'reader_comparison_simple.png'), width=1400, height=250, res=150)
reader_plot(subset(d$readers, chosen))
dev.off()


#### test different buffer capacities and sequence lengths ####

d$cap$capacity = as.integer(gsub('k', '', d$cap$kind)) / 1024
d$cap$seq_length = as.integer(d$cap$name)

png(file.path(outdir, 'bench_cap.png'), width=1400, height=1000, res=200)
ggplot(d$cap, aes(capacity, throughput_gb_s, color=as.factor(seq_length), linetype=toupper(format))) +
    stat_summary(fun=mean, geom='point') +
    stat_summary(fun=mean, geom='line') +
    stat_summary(fun.data=mean_se, geom = 'errorbar', width=0.1, alpha=0.5) +
    expand_limits(y=0) +
    scale_x_continuous(trans='log1p', breaks=2^(0:25)) +
    labs(x='Buffer size [KiB]', y='Throughput [GiB/s]', color='Sequence length', linetype='Format') +
    scale_color_brewer(palette='Set1') +
    theme_bw()
dev.off()
