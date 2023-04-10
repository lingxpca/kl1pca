# Kernelized ℓ1-norm PCA for Denoising

This repository contains the implementation of the Kernel ℓ1-norm Principal Component Analysis (KPCA) algorithm for denoising, as described in the paper:

Kernel ℓ1-norm Principal Component Analysis for Denoising

# Introduction
The proposed method focuses on denoising data using KPCA, which combines a projection-free preimage estimation algorithm with an ℓ1-norm KPCA. This approach is insensitive to outliers and computationally efficient, providing better performance in terms of mean squared error compared to the ℓ2-norm KPCA. The algorithm can be applied to a range of unsupervised learning tasks, such as denoising and clustering.

Getting Started
Prerequisites
* gcc
* Intel-MKL

# Instructions:
Prerequeist: You need install intel one-API mkl library
* make clean && make
* go to exec folder
* ./krpca filename #ofcomponents rbf l1 variance
