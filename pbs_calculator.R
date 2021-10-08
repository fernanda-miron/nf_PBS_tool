# This script can take Fst Data from VCFTools
# and calculate PBS values

# Charging library
library("dplyr")
library("ggplot2")
library("ggrepel")
library("tidyr")
library("cowplot")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "./results"
#args[2] <- "./test/data/pbs.png"

## Place args into named object
file_dir <- args[1]
png_file <- args[2]

## Get all the fst.csv files in path
csv.names <- list.files(path = file_dir, pattern = "*.csv", full.names = T)

# Define a function to import and process fst output from vcftools
fst_reader<-function(filename,pops){
  fst<-read.csv(file = filename, header = T, stringsAsFactors = F) #read in fst output from vcftools
  colnames(fst)[1]<-paste(pops,".Fst",sep = "") #change fst column name to identify pops
  fst[,1][which(is.na(fst[,1]))]<-0 #change NA's to 0, the NA's are produced by vcftools when there is no variation at a site
  assign(pops,fst,envir = .GlobalEnv) #export dataframe with pop names as variable names
}

## Getting dataframe names
pop1 <- as.list(strsplit(csv.names[1], "/") [[1]])
pop1_saved <- as.character(pop1[[3]])

pop2 <- as.list(strsplit(csv.names[2], "/") [[1]])
pop2_saved <- as.character(pop2[[3]])

pop3 <- as.list(strsplit(csv.names[3], "/") [[1]])
pop3_saved <- as.character(pop3[[3]])

## Uploading fst files
df1 <- fst_reader(filename = csv.names[1], pops = pop1_saved )
df2 <- fst_reader(filename = csv.names[2], pops = pop2_saved )
df3 <- fst_reader(filename = csv.names[3], pops = pop3_saved )


# Make a list of your dataframes and join them all together
predata <- left_join(x = df1,
                            y = df2,
                            by = c("Initial_position", "Final_Position", "SNP_Number",
                                   "gene_lenght", "SNP_by_1KB"))

all_fst <- left_join(x = predata,
                     y = df3,
                     by = c("Initial_position", "Final_Position","SNP_Number",
                            "gene_lenght", "SNP_by_1KB"))

# Making a function for PBS calculation
dopbs<-function(pop1_out,pop1_pop2,pop2_out) {
  Tpop1_out= -log(1-pop1_out)
  Tpop1_pop2= -log(1-pop1_pop2)
  Tpop2_out= -log(1-pop2_out)
  pbs= (Tpop1_out + Tpop1_pop2 - Tpop2_out)/2
  pbs
}

# Running for my data
pbsresults<-dopbs(pop1_out = all_fst[1],
              pop1_pop2 = all_fst[7],
              pop2_out = all_fst[8])

# Turn PBS into data frame and merge
pbsresults_2 <-all_fst[,2:6]
pbsresults_2 <-cbind(pbsresults_2 ,pbsresults)
pbsresults_2 <- pbsresults_2 %>% 
  rename(PBS_value = paste0(pop1_saved,".Fst"))

#set negative pbs values to 0 for convenience
pbsresults_2$PBS_value[which(pbsresults_2$PBS_value<0)]<-0

# Arranging data to see values
arreglado <- pbsresults_2[order(-pbsresults_2$PBS_value),]

# Uploading diccionary
# Getting dictionary name
dictionary <- list.files(path = file_dir, pattern = "*.txt", full.names = T)

diccionary <- read.table(file = dictionary, header = T,
                         sep = "\t", stringsAsFactors = F)
colnames(diccionary) <- c("ID", "Initial_position", "Final_Position", "Gene_name")

# Merging to know the gen name
final_dataset <- left_join(x = arreglado,
                     y = diccionary,
                     by = c("Initial_position", "Final_Position"))
# Plotting
p1 <- ggplot(data = final_dataset, aes(x = SNP_by_1KB , y = PBS_value )) +
  geom_point(color = "#F95738", size = 2)
p1

# Obteining better plot
p2 <- p1 + labs(title = "PBS by Genetic Region") +
  xlab("Number of SNPs by 10KB")+
  ylab("PBS") +
  scale_x_continuous(breaks = seq(0, 
                                  ceiling(max(final_dataset[5])), 
                                  by = 100)) +
  geom_hline( yintercept = 0.060, lty = "dashed" ) +
  theme_light()
p2

# Labbeling
p3 <- p2 + geom_label_repel(data= final_dataset %>% filter(PBS_value > 0.060),
                      aes(label=Gene_name),
                          color = 'black',
                          size = 4)

p3

## For the gen with major PBS, we will make two graphs
## Getting initial position of gen with higher PBS
major_PBS <- final_dataset[1,1]
gene_name <- final_dataset[1,8]

## Using function for reading
fst_reader_snp<-function(filename,pops){
  fst<-read.table(file = filename, header = T, stringsAsFactors = F, sep = "\t") #read in fst output from vcftools
  colnames(fst)[3]<-paste(pops,".Fst",sep = "") #change fst column name to identify pops
  fst[,3][which(is.na(fst[,3]))]<-0 #change NA's to 0, the NA's are produced by vcftools when there is no variation at a site
  assign(pops,fst,envir = .GlobalEnv) #export dataframe with pop names as variable names
}

## Reading data
snp.name <- list.files(path = file_dir, pattern = paste0("*",major_PBS,"*"), full.names = T)

## Preparing data
pop1_snp <- as.list(strsplit(snp.name[1], "/") [[1]])
pop1s_saved <- as.character(pop1[[3]])

pop2_snp <- as.list(strsplit(snp.name[2], "/") [[1]])
pop2s_saved <- as.character(pop2[[3]])

pop3_snp <- as.list(strsplit(snp.name[3], "/") [[1]])
pop3s_saved <- as.character(pop3[[3]])

## Reading data
psnp1 <- fst_reader_snp(filename = snp.name[1], pops = pop1s_saved)
psnp2 <- fst_reader_snp(filename = snp.name[2], pops = pop2s_saved)
psnp3 <- fst_reader_snp(filename = snp.name[3], pops = pop3s_saved)

# Make a list of your dataframes and join them all together
predata <- left_join(x = psnp1,
                     y = psnp2,
                     by = c("CHROM", "POS"))

all_fst <- left_join(x = predata,
                     y = psnp3,
                     by = c("CHROM", "POS"))

# Running for my data
my_pbs <-dopbs(pop1_out = all_fst[3],
              pop1_pop2 = all_fst[5],
              pop2_out = all_fst[4])

# Turn PBS into data frame and merge
pbsresults<-all_fst[,1:2]
pbsresults<-cbind(pbsresults,my_pbs)
pbsresults <- pbsresults %>% 
  rename(PBS_value = paste0(pop1_saved,".Fst"))

#set negative PBS values to 0 for convenience
pbsresults$PBS_value[which(pbsresults$PBS_value<0)]<-0

# Arranging data to see values
arreglado <- pbsresults[order(-pbsresults$PBS_value),]

# Plotting
p1 <- ggplot(data = arreglado, mapping = aes(x = PBS_value)) +
  geom_histogram( fill="#EE964B") +
  coord_cartesian(ylim=c(0,nrow(arreglado))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = paste("PBS", gene_name),
       y = "Number of SNPs",
       x = "PBS value") +
  theme_bw()
p1

## Making spiderplot graph
## Import and process AF output from vcftools
## Reading data
frq.names <- list.files(path = file_dir, pattern = "*.frq", full.names = T)

## Reading first data frame
AF_1 <- read.table(file = frq.names[1], sep = "\t",
                        header = T, stringsAsFactors = F, 
                        fill = T, row.names = NULL)

## Changing df format
colnames(AF_1) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "AF1_1", "AF2_1")
AF_1 <- transform(AF_1, CHROM = as.numeric(CHROM))
AF_1 <- transform(AF_1, POS = as.numeric(POS))

## Reading second dataframe
AF_2 <- read.table(file = frq.names[2], sep = "\t",
                        header = T, stringsAsFactors = F, 
                        fill = T, row.names = NULL)

## Changing df format
colnames(AF_2) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "AF1_2", "AF2_2")
AF_2 <- transform(AF_2, CHROM = as.numeric(CHROM))
AF_2 <- transform(AF_2, POS = as.numeric(POS))

## Reading third dataframe
AF_3 <- read.table(file = frq.names[3], sep = "\t",
                        header = T, stringsAsFactors = F, 
                   fill = T, row.names = NULL)

## Changing df format
colnames(AF_3) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "AF1_3", "AF2_3")
AF_3 <- transform(AF_3, CHROM = as.numeric(CHROM))
AF_3 <- transform(AF_3, POS = as.numeric(POS))

## Merge data
merged_data.df <- arreglado %>% left_join(AF_1,
                                            by = c("CHROM"="CHROM", "POS"="POS")) %>%
  left_join(AF_2, by = c("CHROM"="CHROM","POS"="POS","N_ALLELES" = "N_ALLELES")) %>%
  left_join(AF_3, by = c("CHROM"="CHROM","POS"="POS","N_ALLELES" = "N_ALLELES")) %>%
  select(2:3,6:7,9:10,12:13)

## Changing formats
fixed_data.df <- merged_data.df %>%
  mutate(AF1_1 = gsub(x = AF1_1,
                        pattern = ".*:",
                        replacement = "")) %>%
  mutate(AF2_1 = gsub(x = AF2_1,
                        pattern = ".*:",
                        replacement = "")) %>%
  mutate(AF1_2 = gsub(x = AF1_2,
                        pattern = ".*:",
                        replacement = "")) %>%
  mutate(AF2_2 = gsub(x = AF2_2,
                        pattern = ".*:",
                        replacement = "")) %>%
  mutate(AF1_3 = gsub(x = AF1_3,
                        pattern = ".*:",
                        replacement = "")) %>%
  mutate(AF2_3 = gsub(x = AF2_3,
                        pattern = ".*:",
                        replacement = ""))

## Add SNP´s name
nombres <- sprintf("SNP%s",seq(1:nrow(fixed_data.df)))
fixed_data.df$SNP <- nombres

## Pivot dataframe
fixed_data.df <- pivot_longer(fixed_data.df, cols = c("AF1_1", "AF1_2", "AF1_3",
                                                      "AF2_1", "AF2_2", "AF2_3" ),
                              names_to = "AF", values_to = "valor")

## FIltering
row_name <- "SNP1"
fst_values <- c("AF2_1", "AF2_2", "AF2_3")
fixed_data.df <- transform(fixed_data.df, valor = as.numeric(valor))

## Ploting first spider plot
spider_uno.p <- fixed_data.df %>%
  filter(SNP == row_name) %>%
  filter(AF == fst_values) %>%
  ggplot( mapping = aes(x = AF, y = valor) ) +
  geom_point( size = 3, color = "#49697F" )
spider_uno.p

## Using geom segment geometry
spider_dos.p <- spider_uno.p +
  geom_segment(
    aes( x = AF, xend = AF,
         y = 0.0, yend = valor), color = "#49697F", size = 1
  )
spider_dos.p

## Coord polar
spider_tres.p <- spider_dos.p + 
  coord_polar()
spider_tres.p

## Adding value
spider_cuatro.p <- spider_tres.p +
  geom_text( aes(label = valor))
spider_cuatro.p

## Improving geometry
spider_cuatro.p <- spider_tres.p +
  geom_text( aes(label = valor),position = position_nudge(y = 0.3))
spider_cuatro.p

# Cleaning plot
spider_cinco.p <- spider_cuatro.p +
  theme_light() +                         
  theme(panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),       
        axis.text.y = element_blank(),       
        axis.title = element_blank()         
  )
spider_cinco.p

## Adding dinamic title
spider_seis.p <- spider_cinco.p +
  labs(title = paste("AF of gene",
                       gene_name, "SNP with higher PBS value")) +
  xlab("Number of SNPs by 10KB")+
  ylab("PBS")
spider_seis.p

## Merging 3 plots
grid1 <- plot_grid(p1, spider_seis.p, align = "h", labels = c('B', 'C'))
grid2 <- plot_grid(p3, grid1, nrow = 2, labels = 'A' )
grid2

## Saving
ggsave(filename = png_file, 
       plot = grid2, 
       device = "png", 
       width = 15, height = 10, units = "in", 
       dpi = 300)
