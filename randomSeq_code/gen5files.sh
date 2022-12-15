for nm in `seq 1 5`
do ./rand.out -t /home/mleone2/ml_class_project/FASTA_CV/Mo2015_EXCpos_Ctx_fold${nm}_testPos.fa -n 40000 -o /home/mleone2/ml_class_project/results/random_seqs_test${nm}.fa
done