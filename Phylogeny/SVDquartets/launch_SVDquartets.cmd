#nexus

begin paup;
 
 log file=svdq_paup_36.log;

 execute /mnt/DATA/B_JUN_reseq/I_phylogeny/A_SVDq_PAUP/JUN_36_best_FILT_edit_POP.min4.nexus;

 outgroup 1;

 svdq nquartets=500000 taxpartition=juncoforms nthreads=20 seed=0 bootstrap=standard nrep=500 evalq=all showScores=yes;

 SaveTrees file=svd.tre;

end;
