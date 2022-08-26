# This script will download all 3' UTRs for a given 
# species using Ensembl's BioMart function 
#================================================

# Load required packages into library 
library(biomaRt)
library(dplyr)
library(glue)

# Create an Ensembl Mart object 
ensmart = useMart('ENSEMBL_MART_ENSEMBL')

# See all available datasets within Ensembl 
datasets <- listDatasets(ensmart)

#Function to download all 3' UTRs from Ensembl
utr3_ensembl <- function(species_ID, species_dataset){
  # Create species specific Mart object
  species_mart <- useMart('ENSEMBL_MART_ENSEMBL', dataset = species_dataset)
  # Download all genes from species specific Mart object 
  genes = getBM('ensembl_gene_id', mart = species_mart)
  # Retrieve all 3' UTRs from species specific Mart object
  utr3s <- getSequence(id=genes[, 1], type='ensembl_gene_id', seqType='3utr', mart=species_mart)
  # Filter the dataset to remove all missing sequences, insufficiently long, and duplicate entries
  final_seqs <- subset(utr3s, `3utr` != 'Sequence unavailable')
  final_seqs <- subset(final_seqs, nchar(`3utr`)>=100)
  final_seqs <- unique(final_seqs)
  # Calculate the percentage of genes that have sufficiently long 3' UTR sequences 
  seq_count <- utr3s %>%
    count(`3utr`, sort = TRUE) 
  NA_seqs <- seq_count[1,2] # Identify all missing sequences 
  num_seqs <- length(utr3s$`3utr`) # Calculate total number of sequences 
  missing_utr_per <- (NA_seqs/num_seqs)*100 # Calculate missing UTR %
  over_100_utr <- length(final_seqs$`3utr`) # Count number of UTRs over 100nt
  over_100_per <- (over_100_utr/num_seqs*100) # Calculate 100nt+ UTR %
  output_str_1 <- glue('% of genes that are missing UTR sequences: {missing_utr_per}')
  output_str_2 <- glue('% of UTRs over 100nt: {over_100_per}')
  print(species_ID)
  print(output_str_1)
  print(output_str_2)
  # Export 3' UTRs to FASTA file 
  file_name <- glue('{species_ID}_3utrs.fasta')
  exportFASTA(final_seqs, file_name)
}

#Download and assess the quality of UTR data in each mammalian species available in Ensembl
utr3_ensembl('Erinaceus_europaeus','eeuropaeus_gene_ensembl')
utr3_ensembl('Sorex_araneus','saraneus_gene_ensembl')
utr3_ensembl('Pteropus_vampyrus','pvampyrus_gene_ensembl')
utr3_ensembl('Myotis_lucifugus','mlucifugus_gene_ensembl')
utr3_ensembl('Vicugna_pacos','vpacos_gene_ensembl')
utr3_ensembl('Camelus_dromedarius','cdromedarius_gene_ensembl')
utr3_ensembl('Delphinapterus_leucas','dleucas_gene_ensembl')
utr3_ensembl('Monodon_monoceros','mmonoceros_gene_ensembl')
utr3_ensembl('Capra_hircus','chircus_gene_ensembl')
utr3_ensembl('Ovis_aries','oaries_gene_ensembl')
utr3_ensembl('Moschus_moschiferus','mmoschiferus_gene_ensembl')
utr3_ensembl('Bos_taurus','btaurus_gene_ensembl')
utr3_ensembl('Bos_indicus','bihybrid_gene_ensembl')
utr3_ensembl('Bos_mutus','bmutus_gene_ensembl')
utr3_ensembl('Bison_bison','bbbison_gene_ensembl')
utr3_ensembl('Catagonus_wagneri','cwagneri_gene_ensembl')
utr3_ensembl('Sus_scrofa','sscrofa_gene_ensembl')
utr3_ensembl('Panthera_pardus','ppardus_gene_ensembl')
utr3_ensembl('Felis_catus','fcatus_gene_ensembl')
utr3_ensembl('Canis_lupus_familiaris','clfamiliaris_gene_ensembl')
utr3_ensembl('Ailuropoda_melanoleuca','amelanoleuca_gene_ensembl')
utr3_ensembl('Ursus_maritimus','umaritimus_gene_ensembl')
utr3_ensembl('Equus_caballus','ecaballus_gene_ensembl')
utr3_ensembl('Equus_asinus','eaasinus_gene_ensembl')
utr3_ensembl('Otolemur_garnettii','ogarnettii_gene_ensembl')
utr3_ensembl('Propithecus_coquereli','pcoquereli_gene_ensembl')
utr3_ensembl('Microcebus_murinus','mmurinus_gene_ensembl')
utr3_ensembl('Callithrix_jacchus','cjacchus_gene_ensembl')
utr3_ensembl('Aotus_nancymaae','anancymaae_gene_ensembl')
utr3_ensembl('Saimiri_boliviensis','sbboliviensis_gene_ensembl')
utr3_ensembl('Cebus_capucinus','ccapucinus_gene_ensembl')
utr3_ensembl('Mandrillus_leucophaeus','mleucophaeus_gene_ensembl')
utr3_ensembl('Cercocebus_atys','catys_gene_ensembl')
utr3_ensembl('Papio_anubis','panubis_gene_ensembl')
utr3_ensembl('Macaca_mulatta','mmulatta_gene_ensembl')
utr3_ensembl('Macaca_fascicularis','mfascicularis_gene_ensembl')
utr3_ensembl('Macaca_nemestrina','mnemestrina_gene_ensembl')
utr3_ensembl('Chlorocebus_sabaeus','csabaeus_gene_ensembl')
utr3_ensembl('Rhinopithecus_roxellana','rroxellana_gene_ensembl')
utr3_ensembl('Rhinopithecus_bieti','rbieti_gene_ensembl')
utr3_ensembl('Nomascus_leucogenys','nleucogenys_gene_ensembl')
utr3_ensembl('Gorilla_gorilla','ggorilla_gene_ensembl')
utr3_ensembl('Pan_troglodytes','ptroglodytes_gene_ensembl')
utr3_ensembl('Pan_paniscus','ppaniscus_gene_ensembl')
utr3_ensembl('Homo_sapiens','hsapiens_gene_ensembl')
utr3_ensembl('Pongo_abelii','pabelii_gene_ensembl')
utr3_ensembl('Oryctolagus_cuniculus','ocuniculus_gene_ensembl')
utr3_ensembl('Ochotona_princeps','oprinceps_gene_ensembl')
utr3_ensembl('Mesocricetus_auratus','mauratus_gene_ensembl')
utr3_ensembl('Microtus_ochrogaster','mochrogaster_gene_ensembl')
utr3_ensembl('Nannospalax_galili','ngalili_gene_ensembl')
utr3_ensembl('Jaculus_jaculus','jjaculus_gene_ensembl')
utr3_ensembl('Dipodomys_ordii','dordii_gene_ensembl')
utr3_ensembl('Ictidomys_tridecemlineatus','itridecemlineatus_gene_ensembl')
utr3_ensembl('Marmota_marmota','mmmarmota_gene_ensembl')
utr3_ensembl('Cavia_porcellus','cporcellus_gene_ensembl')
utr3_ensembl('Octodon_degus','odegus_gene_ensembl')
utr3_ensembl('Chinchilla_lanigera','clanigera_gene_ensembl')
utr3_ensembl('Rattus_norvegicus','rnorvegicus_gene_ensembl')
utr3_ensembl('Mus_musculus','mmusculus_gene_ensembl')
utr3_ensembl('Mus_spretus','mspretus_gene_ensembl')
utr3_ensembl('Mus_caroli','mcaroli_gene_ensembl')
utr3_ensembl('Dasypus_novemcinctus','dnovemcinctus_gene_ensembl')
utr3_ensembl('Choloepus_hoffmanni','choffmanni_gene_ensembl')
utr3_ensembl('Procavia_capensis','pcapensis_gene_ensembl')
utr3_ensembl('Loxodonta_africana','lafricana_gene_ensembl')
utr3_ensembl('Echinops_telfairi','etelfairi_gene_ensembl')