library (dplyr)
library (tidyr)
library (ggplot2)


Coverage <- read.table("Results/Coverage.all", stringsAsFactors = FALSE, header=T)
Coverage <- rename ( Coverage, cov=Coverage, pos=Pos, ID=sample)
###### Functions #####
#coverage

slide <- function(cov.df, setup.df)
{
  cov = rep(NA, nrow(setup.df))
  for(i in 1:nrow(setup.df))
  {
    s = setup.df$starts[i]
    e = setup.df$ends[i]
    subset(cov.df, pos >= s & pos < e, select = c(cov)) -> position
    mean(position$cov) -> cov[i]
  }
  out <- data.frame(mean =cov, pos = setup.df$pos)
  out$ID = unique(cov.df$ID)
  return(out)
}

cov_plot <- function(cov.df, title)
{
  cov.df  %>% summarize(first = min(pos), last = max(pos)) %>% plyr::adply(1,function(x) data.frame(starts = seq(x$first,x$last,by=400))) %>% mutate(ends = ifelse(starts + 400 < last, starts + 400, last)) -> setup
  setup %>% select(starts, ends) -> setup_means
  setup$pos <- apply(setup_means, 1, function(x) mean(x))
  plyr::ddply(cov.df, ~ID, slide, setup) -> cov.slid.df
  
  
  cov.plot <- ggplot(cov.slid.df, mapping = aes(x = as.factor(pos), y = mean)) + geom_boxplot(fill="white")
  cov.plot <- cov.plot + ggtitle(title) + ylab("Read depth")  + xlab(" Genome Position")
  cov.plot <- cov.plot + theme(axis.title.y = element_text(vjust=1.2))
  cov.plot <- cov.plot + theme(legend.position = "none") + theme_classic()
  cov.plot <- cov.plot + theme(text = element_text(size = 15), axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 10))
  return(cov.plot)
}

# ======================== Plot ===================


cov_plot(Coverage, title = "") -> coverage.plot.all


ggsave(filename = "../Results/Plots/coverage_plot.pdf", plot = coverage.plot.all, device = "pdf", width = 6, height = 4)

