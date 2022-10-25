for fnm in Neg10x Pos
do python seq_functions.py --fasta_name /home/mleone2/ml_class_project/FASTA_CV/Mo2015_EXCpos_Ctx_fold1_test${fnm}.fa --prediction_file /home/mleone2/ml_class_project/results/${fnm}.txt
done
python seq_functions.py --fasta_name /home/mleone2/ml_class_project/results/random_seqs_test.fa --prediction_file /home/mleone2/ml_class_project/results/random_seqs_test_predictions.txt
