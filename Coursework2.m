%% Leticia Viada Campos
% Mathematical Methods Computational homework 2
%
clear all, close all, clc, format long, format compact
%%
%
% Exercise 1
%
% (a)
% We need A to be strictly symmetric diagonal dominant matrix. Then
% $a_{ii}>\Sigma(a_{ij})$ for $i!=j$. As $\epsilon \in [0,1]$, we need $1>
% 2\epsilon + 2\epsilon^2$. Thus, $0<\epsilon<\frac{\sqrt(3)-1}{2}$. The
% upperbound is approximately 0.37
%%
% (b)
tol=1e-10;
nmax=1e3;
x0=zeros(5,1);
[A,b]=matrix(5,0.3);
[xJac, niterJac, relresiterJac2, xiterJac2] = itermeth(A,b,x0,nmax,tol,'J');
[xGS, niterGS2, relresiterGS, xiterGS] = itermeth(A,b,x0,nmax,tol,'G');
niterJac % Number of iterations for the Jacobi method
niterGS2 % Number of iterations for the Gauss-Seidel method
%%
%
% (c)
spec_BJ=[]; % Array to store all the spectral values for the different BJ
spec_BG=[]; % Array to store all the spectral values for the different BGS
for epsi=0:0.01:1 % for loop to find the values
    [A,b]=matrix(5,epsi);
    D=diag(diag(A)); L=tril(A)-D;U=triu(A)-D;
    BJ=-inv(D)*(L+U); eigBJ=eig(BJ); spectralBJ=max(abs(eigBJ)); % find the spectral value
    spec_BJ(end+1)=spectralBJ; % add the value at the end of the array spec_BJ
    BG=-inv(L+D)*U; eigBG=eig(BG); spectralBG=max(abs(eigBG)); % find the spectral value
    spec_BG(end+1)=spectralBG; % add the value at the end of the array spec_BG
end
eps=[0:0.01:1];
figure
plot(eps,spec_BJ)
hold on
plot(eps,spec_BG)
yline(1,'LineWidth',1.5)
grid on
xlabel('\epsilon')
legend('specBJ','specBG')
%%
%
% Theorem 4.6.1 states that if the method is consistent then it converges for
% all initial guesses iff the spectral radius is less than 1. For the Jacobi
% method, we need $\epsilon<0.45$ and for the Gauss-Seidel method we need
% $\epsilon<0.71$. Putting both of them together we get that for both methods 
% to converge $\epsilon<0.434$ needs to be satisfied. 
%
% In part (a) we have found a more restricted interval, but we can see
% that for some other values for which A is not strictly diagonal dominant,
% the Jacobi method and the Gauss-Seidel method converge.
% Moreover, the smaller the spectral value, the faster a method converges. 
% For the stated values 
% of $\epsilon$ we can expect the Gaus-Seidel Method to converge faster.
% would recommend using the Gauss-Seidel Method if $\epsilon=0.5$ as the
% Jacobi method will not converge for that value.
[A,b]=matrix(5,0.5);
[xJac2, niterJac2, relresiterJac2, xiterJac2] = itermeth(A,b,x0,nmax,tol,'J');
[xGS2, niterGS2, relresiterGS2, xiterGS2] = itermeth(A,b,x0,nmax,tol,'G');
niterJac2 % Number of iterations for the Jacobi method when \epsilon=0.5
niterGS2 % Number of iterations for the Gauss-Seidel method \epsilon=0.5