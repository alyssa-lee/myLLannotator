#!/usr/bin/env bash

echo "filtering AMR proteins in background..."
(head -n1 ../data/Maddamsetti2024/all-proteins.csv && cat ../data/Maddamsetti2024/all-proteins.csv | grep -E "chloramphenicol|Chloramphenicol|tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating|macrolide|lincosamide|streptogramin|Multidrug resistance|multidrug resistance|antibiotic resistance|lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\S*|glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase|bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase|trimethoprim|dihydrofolate reductase|dihydropteroate synthase|sulfonamide|Sul1|sul1|sulphonamide|quinolone|Quinolone|oxacin|qnr|Qnr|Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt|macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA|QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\S*") \
> ../results/filtered-AMR-proteins.csv &

echo "filtering MGE proteins in foreground..."
(head -n1 ../data/Maddamsetti2024/all-proteins.csv && cat ../data/Maddamsetti2024/all-proteins.csv | grep -E "IS|transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|tra[A-Z]|conjugate transposon|Transpos\S*|Tn[0-9]|tranposase|Tnp|Ins|ins|relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation|Mob\S*|Plasmid|Rep|Conjug\S*|capsid|phage|Tail|tail|head|tape measure|antiterminatio|Phage|virus|Baseplate|baseplate|coat|entry exclusion|Integrase|integrase|excision\S*|exonuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip|intron|antitoxin|toxin|Toxin|Reverse transcriptase|hok|Hok|competence|addiction") \
> ../results/filtered-MGE-proteins.csv




