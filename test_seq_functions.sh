for idx in `seq 1 5`
    # do for fnm in Neg10x Pos
    #     do python seq_functions.py --fasta_name /home/mleone2/ml_class_project/FASTA_CV/Mo2015_EXCpos_Ctx_fold${idx}_train${fnm}.fa --prediction_file /home/mleone2/ml_class_project/results/fold${idx}_train${fnm}_10000result.txt
    # done
    do python seq_functions.py --fasta_name /home/mleone2/ml_class_project/results/random_seqs_test${idx}.fa --prediction_file /home/mleone2/ml_class_project/results/random_seqs_test${idx}_predictions_10000result.txt
done