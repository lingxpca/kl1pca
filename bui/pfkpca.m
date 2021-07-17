load t2;
variance = 0; seed = 1;
rng(seed);
X = true + normrnd(0,sqrt(variance),size(true));
data2.X=X;
rho_m=1;
L = 4;
myKPCA(data2, rho_m, L, 1:size(X,1))