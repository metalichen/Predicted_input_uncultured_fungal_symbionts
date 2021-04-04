# Predicted_input_uncultured_fungal_symbionts

Scripts used in the analysis for the "Predicted input of uncultured fungal symbionts to a lichen symbiosis from metagenome-assembled genomes" paper, published in [Genome Biology and Evolution](https://academic.oup.com/gbe/advance-article/doi/10.1093/gbe/evab047/6163286) in 2021.

## Abstract
Basidiomycete yeasts have recently been reported as stably associated secondary fungal symbionts (SFSs) of many lichens, but their role in the symbiosis remains unknown. Attempts to sequence their genomes have been hampered both by the inability to culture them and their low abundance in the lichen thallus alongside two dominant eukaryotes (an ascomycete fungus and chlorophyte alga). Using the lichen Alectoria sarmentosa, we selectively dissolved the cortex layer in which SFSs are embedded to enrich yeast cell abundance and sequenced DNA from the resulting slurries as well as bulk lichen thallus. In addition to yielding a near-complete genome of the filamentous ascomycete using both methods, metagenomes from cortex slurries yielded a 36- to 84-fold increase in coverage and near-complete genomes for two basidiomycete species, members of the classes Cystobasidiomycetes and Tremellomycetes. The ascomycete possesses the largest gene repertoire of the three. It is enriched in proteases often associated with pathogenicity and harbours the majority of predicted secondary metabolite clusters. The basidiomycete genomes possess ∼35% fewer predicted genes than the ascomycete and have reduced secretomes even compared to close relatives, while exhibiting signs of nutrient limitation and scavenging. Furthermore, both basidiomycetes are enriched in genes coding for enzymes producing secreted acidic polysaccharides, representing a potential contribution to the shared extracellular matrix. All three fungi retain genes involved in dimorphic switching, despite the ascomycete not being known to possess a yeast stage. The basidiomycete genomes are an important new resource for exploration of lifestyle and function in fungal-fungal interactions in lichen symbioses.


## List of scripts:
### coverage.sh  
Script producing coverage files used to make GC-coverage plots  

### gc-coverage-plot.R  
R script used to make the GC-coverage plots  

### funannotate_annotation.sh  
Wrapper script used to produce genome annotations using the funannotate pipeline  

### get_secretome.R
R scrip used to combine annotations produced by SignalP, TMHMM, and WolfPSORT to produce secretome
