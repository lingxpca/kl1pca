function [Xq] = myKPCA( data, rho_m, nL_vector, whidx )
%myKPCA KPCA with the Gaussian RBF kernel
%   
% 2 methods; (1) Scholkopf, and (3) Contour gradient
format shortG 
[N, n]=size(data.X); % N: sample size; n: dimension
oneOverN = 1/N;
ones1n = ones(1,n);
Ntest=length(whidx);

% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_exh=500; % value of I in the paper (algorithm part)
h=1.2;

% standardize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=data.X; %


clear data;
tic;
% inter-distances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm2=zeros(N,N);
for i=1:N,
    for j=(i+1):N,
        norm2(i,j)=(X(i,:)'-X(j,:)')'*(X(i,:)'-X(j,:)');
        norm2(j,i) =norm2(i,j);
    end
end 

% rho equals mean component variance times rho_m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_comp_var=mean(var(X));
rho=(mean_comp_var)*size(X,2)*rho_m;% Bui: the same with rho = sum(Var(X))*rho_m;

% compute the kernel matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=zeros(N,N);
for i=1:N,
    for j=i:N,
        K(i,j) = 1+(X(i,:)')'*(X(j,:)');
        K(j,i) = K(i,j);
    end
end 

%for i=1:N,
%    for j=i:N,
%        K(i,j) = exp(-norm2(i,j)/rho);
%         K(j,i) = K(i,j);
%     end
% end 

% centering the kernel matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unit = ones(N, 1);
unit2=ones(N,N); % buffer
Kone = K*unit;
oneKone = unit'*Kone; % buffer
oneKoneUnit2 = oneOverN^2*oneKone*unit;


K_tilde = K - unit2*K/N - K*unit2/N + unit*oneKone*unit'/N^2;
K_tilde=min(K_tilde, K_tilde')

% find eigenvalues and eigenvectors of the centered kernel matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[evecs,evals] = eig(K_tilde)
evals = real(diag(evals));


% scree plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(1,1,1)
% plot(1:N, evals(N:(-1):1), 'g')
% xlim([1,length(X)])
% ylim([0, evals(N)]);
% pause;
scree_val=evals(N:(-1):1);
% denoising by 2 different methods
% loop over different choices of the number of retained eigenvectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for vv=1:length(nL_vector),
    
    vv;
    nL=nL_vector(vv);
    
    large_idx=N-nL+(1:nL);
    alpha=evecs(:,large_idx)
    levals=evals(large_idx);
    diag(levals.^(-0.5));
    alpha=alpha*diag(levals.^(-0.5));%normaiization of v(=e-vector in feature space)
    alpha*alpha'
    Xp1=zeros(Ntest,n);
%     Xp2=Xp1;
    Xq=Xp1;
   
    for ii=1:Ntest,
        i=whidx(ii);
        if i==1,
            idx_lst=2:N;
        else
            idx_lst=[1:(i-1), (i+1):N];
        end
        
        x=X(i,:)';%2 to i
        
      
        % Step S1:
        [g0, v0]=obj_rbf_dirS1(x, X, Kone, alpha, rho, oneKone, unit, oneOverN, ones1n, unit2, oneKoneUnit2); % compute gradient g0 and steepest direction v0 
      
        v=-v0/norm(v0)% normalize direction v to unit vector
       
        
        if g0>10^(-9),
            % Step S2:
            pjxn_v=X*v;
            min_pjxn_v=min(pjxn_v(idx_lst))- pjxn_v(i);
            max_pjxn_v=max(pjxn_v(idx_lst))- pjxn_v(i);

   
            clear y_exh; % b/c we use length(y_exh) below
            delta = max_pjxn_v/(n_exh-1);
            y_exh(1)=delta;

            k=2;
            while y_exh(length(y_exh)) < max_pjxn_v && delta > 0
                delta =  delta*h;
                y_exh(k)=y_exh(k-1) + delta;
                k=k+1;
            end
            k=k-1;
            if k>1
                y_exh(k)= min(max_pjxn_v, y_exh(k));   
            end

            % Step S3:
            buff = g0;
            for i_exh=1:k, 
                f_exh=obj_rbfS3(x+y_exh(i_exh)*v, X, Kone, alpha, rho, unit, oneOverN, unit2, oneKoneUnit2,oneKone);
                if f_exh > buff,
                    if i_exh>1
                        gp_exh=rbf_dirS3(x+y_exh(i_exh-1)*v, X, Kone, alpha, rho, v, unit, oneOverN, ones1n, unit2, oneKoneUnit2);
                        if gp_exh>0
                            i_exh=i_exh-1; 
                        else
                            i_exh =i_exh; 
                        end
                    end
%                     i_exh
                    break
                end
                buff = f_exh;
            end
            
           
            % Step S4:           
            if i_exh>1
                Xq(ii,:)=(x + (y_exh(i_exh-1)+y_exh(i_exh))/2*v)';
            else       
                Xq(ii,:)= (x + y_exh(i_exh)*v)';
            end
                 
        else
            Xq(ii,:)=x';
        end
    end
%         Xq=Xq.*(ones(Ntest,1)*stdevX)+ones(Ntest,1)*meanX; %Unstandardize   
end
end


function [f, gradf]=stage2_rbf(x, X, K, c, alpha, rho)

N=size(X,1);
n=size(X,2);
Kx=zeros(N,1);
for i=1:N,
    Kx(i,1)=exp(-(x-X(i,:)')'*(x-X(i,:)')/rho);
end
unit=ones(N,1);
JKx=-2/rho*(Kx*ones(1,n)).*(unit*x'-X);


Kx_tilde=Kx-1/N*K*unit-1/N*unit*unit'*Kx+1/N^2*(unit'*K*unit)*unit;
JKx_tilde=JKx-1/N*unit*unit'*JKx;


f=1-2/N*unit'*Kx+1/N^2*unit'*K*unit-2*(Kx_tilde'*alpha)*c;

if nargout>1,
    gradf=-2/N*(unit'*JKx)'-2*(JKx_tilde'*alpha)*c;
end

end


function [f, gradf]=obj_rbf_dirS1(x, X, Kone, alpha, rho, oneKone, unit, oneOverN, ones1n, unit2, oneKoneUnit2)
tg = unit*x'-X; %(Nxn)

    Kx=exp(-sum(tg.^2,2)/rho) ;% (Nx1)
    JKx=-2/rho*(Kx*ones1n).*tg; % derivative of Kx with respect to x (N*n)
   
    Kx_tilde=Kx-oneOverN*(Kone-unit2*Kx)+ oneKoneUnit2;
    JKx_tilde=JKx-oneOverN*unit2*JKx;% centralized JKx
    f=1-2*oneOverN*unit'*Kx+oneOverN^2*oneKone-(Kx_tilde'*alpha)*(alpha'*Kx_tilde);
    gradf=-2*oneOverN*(JKx'*unit)-2*(JKx_tilde'*alpha)*(alpha'*Kx_tilde); 
end

function f=obj_rbfS3(x, X, Kone, alpha, rho, unit, oneOverN, unit2, oneKoneUnit2, oneKone)
tg = unit*x'-X;

    Kx=exp(-sum(tg.^2,2)/rho);

    Kx_tilde = Kx - oneOverN*(Kone - unit2*Kx) + oneKoneUnit2;
    
    f=1-2*oneOverN*unit'*Kx+oneOverN^2*oneKone-(Kx_tilde'*alpha)*(alpha'*Kx_tilde);
    
end

function gradf=rbf_dirS3(x, X, Kone, alpha, rho, v, unit, oneOverN, ones1n, unit2, oneKoneUnit2)
tg = unit*x'-X;
    
    Kx=exp(-sum(tg.^2,2)/rho);
    JKx=-2/rho*(Kx*ones1n).*tg;

    Kx_tilde = Kx - oneOverN*(Kone - unit2*Kx) + oneKoneUnit2;
    JKx_tilde=JKx - oneOverN*unit2*JKx;
    2*oneOverN*(JKx'*unit)+(JKx_tilde'*alpha)*(alpha'*Kx_tilde);
    gradf=-2*(oneOverN*(JKx'*unit)+(JKx_tilde'*alpha)*(alpha'*Kx_tilde))'*v;
    
end