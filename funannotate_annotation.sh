#Genome annotation using funannotate pipeline

ASSEMBLYPATH=./rintermedia_PEKF01.1.fsa_nt
SPECIES="Ramalina intermedia"
CPUS=12
BASENAME="ram_int"
SP=${SPECIES// /_}
echo $SP

ASSEMBLY=$(basename "$ASSEMBLYPATH" .fna)

funannotate clean -i $ASSEMBLYPATH -o "$ASSEMBLY"_clean.fa

funannotate sort -i "$ASSEMBLY"_clean.fa -o "$ASSEMBLY"_clean_sort.fa -b $BASENAME

funannotate mask -i "$ASSEMBLY"_clean_sort.fa -o "$ASSEMBLY"_clean_sort_masked.fa --cpus $CPUS

#the next line needs to be changes dependening if masking was performed or not!!!
funannotate predict -i "$ASSEMBLY"_clean_sort_masked.fa -o "$BASENAME"_preds -s "Ramalina intermedia" --optimize_augustus --cpus $CPUS --name "$BASENAME"_predicted_gene

funannotate iprscan -i "$BASENAME"_preds -m local --iprscan_path /scratch/1/gulnara/interproscan/interproscan-5.34-73.0/interproscan.sh -c $CPUS

funannotate remote -i "$BASENAME"_preds -m phobius -e tagirdzh@ualberta.ca

#the next command is run in the predictions directory:
cd "$BASENAME"_preds
SP=${SPECIES// /_}
echo $SP
emapper.py  -i predict_results/"$SP".proteins.fa --output eggnog_results -d euk --data_dir /scratch/1/gulnara/annotations/own_libraries/bryoria_tortuosa/bry_tor_preds/data --cpu $CPUS --override -m diamond
cd ..

funannotate annotate -i "$BASENAME"_preds --sbt ../bryoria_tortuosa/template.sbt.txt --eggnog "$BASENAME"_preds/eggnog_results.emapper.annotations --cpus $CPUS
