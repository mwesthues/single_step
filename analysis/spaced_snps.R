if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "rvest", "magrittr")

## -- DATA SELECTION -----------------------------------------------------
genos <- readRDS("./data/processed/common_genotypes.RDS")

# snp
snp <- readRDS("./data/processed/snp_mat.RDS")
snp_nms <- genos %>%
  filter(Data_Type == "snp") %>% 
  .$G
snp <- snp %>%
  .[match(snp_nms, rownames(.)), ]

# Load data frame with meta information on marker names and their positions in 
# the genome.
snp_map <- fread("./data/raw/marker_map.txt", sep = "\t")
snp_map <- snp_map[complete.cases(snp_map), ]
snp_map <- snp_map[snp_map$chromosome != 0, ]
map_marker_names <- as.character(snp_map[["markername"]])
map_marker_names <- gsub(pattern = "[.]", replacement = "_",
                         x = map_marker_names)
map_marker_names <- gsub(pattern = "[-]", replacement = "_", 
                         x = map_marker_names)
snp_map[, markername := map_marker_names]

# Get the length of each of the ten maize chromosomes. We need this information
# so that we can align the chromosomes sequentially and treat them like an 
# entire maize genome with a global length. Thereby, we can select equally 
# spaced SNPs with respect to the genome instead of individual chromosomes.
ensembl <- read_html("http://ensembl.gramene.org/Zea_mays/Info/Annotation/")
chr_lengths <- ensembl %>%
  html_nodes(., "#chromosome_table td:nth-child(2)") %>%
  html_text(.) %>%
  .[seq_len(10)] %>%
  c(0, .) %>%
  .[seq_len(10)] %>%
  as.list() %>%
  map(as.data.frame) %>%
  bind_rows(.id = "chromosome") %>%
  rename(Chr_Length = `.x[[i]]`) %>%
  mutate(Chr_Length = gsub(",", replacement = "", x = Chr_Length),
         Chr_Length = as.numeric(Chr_Length),
         Chr_Length = cumsum(Chr_Length),
         chromosome = as.integer(chromosome))

# Create global SNP positions
global_map <- snp_map %>%
  left_join(y = chr_lengths, by = "chromosome") %>%
  as_tibble() %>%
  rename(Chr_position = position) %>%
  mutate(Genome_position = Chr_position + Chr_Length) %>%
  arrange(Genome_position)

# Function to choose m evenly spaced elements from a sequence of length n
# http://stackoverflow.com/a/9873804/2323832
space_m_equally <- function(m, n) {
  i <- seq(from = 0, to = m - 1)
  i * n %/% m + n %/% (2 * m)
}
# Return the SNP loci that are physically closest to the previously selected,
# artificially equi-distant loci.
select_closest_matches <- function(x, m) {
  y <- space_m_equally(m, n = max(x))
  vapply(y, FUN = function(i) {
    x[which.min(abs(i - x))]
  }, FUN.VALUE = double(1))
}
# Compute the distance between two adjacent loci.
dist_neighbor <- function(x) {
  y <- x - c(0, x[seq(from = 1, to = length(x) - 1)])
  y[1] <- x[2] - x[1]
  y
}

# Select different sets of equi-distant loci:
# 125, 250, 500, 1000, 2500, 5000, 10000 SNPs
# The actual numbers in 'n_snps' differ from the specified set sizes because
# it is possible that the same locus is the closest neighbor to more than one
# breakpoint.
n_snps <- c(125, 250, 500, 1007, 2548, 5144, 10682)
even_snp_lst <- lapply(n_snps, FUN = function(snp_number) {
  global_map %>%
    distinct(Genome_position, .keep_all = TRUE) %>%
    filter(Genome_position %in% select_closest_matches(Genome_position,
                                                       m = snp_number))
})
names(even_snp_lst) <- n_snps

# For each number of equally spaced loci, summarize the distance between each
# locus and its next neighbor.
even_snp_lst %>%
  map("Genome_position") %>%
  map(dist_neighbor) %>%
  map(summary)
# Plot distance-distribution of neighbors.
even_snp_lst %>%
  bind_rows(.id = "SNP_Number") %>%
  mutate(SNP_Number = as.factor(SNP_Number)) %>%
  group_by(SNP_Number) %>%
  mutate(Neighbor_dist = dist_neighbor(Genome_position)) %>%
  ungroup() %>%
  select(SNP_Number, Neighbor_dist) %>%
  ggplot(., aes(x = Neighbor_dist)) +
  geom_histogram(fill = "white", color = "black") +
  facet_wrap(~ SNP_Number, scales = "free")
