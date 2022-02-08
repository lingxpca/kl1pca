# Kernelize L1 norm PCA for Denoising

In this paper we describe a method for denoising data using kernelized principal component analysis (PCA) that is able to recover preimages of the intrinsic variables in the feature space using a single line search along the gradient descent direction of its squared projection error. This method combines a projection-free preimage estimation algorithm with an L1-norm kernel PCA (KPCA). Those two stages provide distinct dvantages to other KPCA preimage methods in the sense that it is insensitive to outliers and computationally efficient. The method can enhance the results of a range of unsupervised learning tasks such as denoising,
clustering, and dimensionality reduction. Numerical experiments on Amsterdam Library of Object Images nd the MNIST handwritten digits demonstrate the proposed method performs better in terms of mean quared error than L2-norm analogue as well as on toy synthetic data. The proposed method is applied to different data sets and the performances are reported.

# Instructions:
* make clean && make
* go to exec folder
* ./krpca filename #ofcomponents rbf l1 variance
