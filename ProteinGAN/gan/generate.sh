#!/bin/bash
#SBATCH -p moh1
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH -t 04-00:00:00
#SBATCH -o /home/mleone2/ml_class_project/ProteinGAN/src/%j.out
#SBATCH -e /home/mleone2/ml_class_project/ProteinGAN/src/%j.err

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "start at $tic\n\n";


python -u -m generate --batch_size 16 --name test_run --steps 3000 -shuffle_buffer_size 100000 --loss_type non_saturating --discriminator_learning_rate 0.0001 --generator_learning_rate 0.0001 --dilation_rate 2 --n_seqs 1000 --gf_dim 44 --df_dim 30 --dataset protein/our_data --nouse_cpu --architecture gumbel --pooling conv


toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "finish at $toc\n\n";
