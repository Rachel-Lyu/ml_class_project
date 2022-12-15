#!/bin/bash
#SBATCH -p moh1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH -t 04-00:00:00
#SBATCH -o /home/mleone2/ml_class_project/ProteinGAN/src/%j.out
#SBATCH -e /home/mleone2/ml_class_project/ProteinGAN/src/%j.err

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "start at $tic\n\n";

python -u -m train_gan --batch_size 16 --name moh1_run_16_10000_step_2 --steps 10000 -shuffle_buffer_size 100000 --loss_type non_saturating --discriminator_learning_rate 0.0001 --generator_learning_rate 0.0001 --dilation_rate 1 --gf_dim 44 --df_dim 30 --dataset protein/our_data --architecture gumbel --pooling conv 

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "finish at $toc\n\n";
