%% Leticia Viada Campos
% Mathematical Methods Computational homework 2
%
clear all, close all, clc, format long, format compact
%%
%
% Exercise 2
%
% (a)
rhoBJ=[];
rhoBG=[];
for N=[5,10,20,40,80]
    h=1\N;
    A=(2/h^2)*diag(ones(N-1,1))-(1/h^2)*diag(ones(N-2,1),1)...
        -(1/h^2)*diag(ones(N-2,1),-1);
    D=diag(diag(A)); L=tril(A)-D;U=triu(A)-D;
    BJ=-inv(D)*(L+U); eigBJ=eig(BJ); spectralBJ=max(abs(eigBJ)); % find the spectral value
    rhoBJ(end+1)=spectralBJ; % add the value at the end of the array rhoBJ
    BG=-inv(L+D)*U; eigBG=eig(BG); spectralBG=max(abs(eigBG)); % find the spectral value
    rhoBG(end+1)=spectralBG; % add the value at the end of the array rhoBG
end
%%
%
% We have seen in theorem 4.6.1 that a consistent method is convergent iff
% its spectral value is less than 1. In the arrays above we have that
% condition satisfied. Thus, for these values of N I expect both methods to
% be convergent. Moreover, the smaller the spectral radius the faster it
% converges, so we can see that the Gauss-Seidel method converges faster
% than the Jacobi method.
Nvec=[5,10,20,40,80];
loglog(Nvec,1-rhoBJ)
hold
grid on
loglog(Nvec,1-rhoBG)
legend('BJ','BS')
%%
%
% We can expect the relation to be $\rho(B)=1-CN^{m}$. We have
% sketched the graphs $1-\rho(B)$ and obtained straight lines for both methods.
% The slopes are $m$ and offset $D=\log(C)$, thus $C=e^{D}$. For the Jacobi
% method we get
x=linspace(5,80);
mbj=(log(1-rhoBJ(end))-log(1-rhoBJ(1)))/(log(Nvec(end))-log(Nvec(1)))
dbj=log(1-rhoBJ(1))-mbj*log(Nvec(1));
cbj=exp(dbj)
f=@(x)cbj*x.^mbj;
%%
%
% Similarly for the Gauss-Seidel method we obtain
mbg=(log(1-rhoBG(end))-log(1-rhoBG(1)))/(log(Nvec(end))-log(Nvec(1)))
dbg=log(1-rhoBG(1))-mbg*log(Nvec(1));
cbg=exp(dbg)
g=@(x)cbg*x.^mbg;
%%
% 
% We therefore get that for the Jacobi method an approximation is
% $\rho(B)=1-7.99N^{(-1.95)}$ and for the Gauss-Seidel method an approximation is
% $\rho(B)=1-4.68N^{(-1.98)}$
% We have that the error is proportional to $\rho(B)^k$. To attain a
% tolerance, tol, $\rho(B)^k\leq tol$. Applying logarithms at boths sides
% we get that k increases as $\frac{\log(tol)}{\log(1-CN^{m})}$. Thus it has
% order $O(\frac{1}{\log(1-CN^{m})})$, making it converge faster than direct
% methods for large N. These direct methods have $O(N^3)$
%%
% 
% (b)
tol=1e-10;
nmax=1e5;
for N=[5,10,20,40,80]
    h=1/N;
    u0=zeros(N-1,1);
    A=(2/h^2)*diag(ones(N-1,1))-(1/h^2)*diag(ones(N-2,1),1)...
        -(1/h^2)*diag(ones(N-2,1),-1);
    b=transpose(sin(pi*h*(1:N-1)));
    if N==5
        u5 = itermeth(A,b,u0,nmax,tol,'G'); % find the solution for N=5
        u5=[0;u5;0]; % add u_0 and u_N
        x5=linspace(0,1,N+1); % find the range to sketch the x_N
    elseif N==10
        u10 = itermeth(A,b,u0,nmax,tol,'G');% find the solution for N=10
        u10=[0;u10;0]; % add u_0 and u_N
        x10=linspace(0,1,N+1); % find the range to sketch the x_N
    elseif N==20
        u20=itermeth(A,b,u0,nmax,tol,'G'); % find the solution for N=20
        u20=[0;u20;0]; % add u_0 and u_N
        x20=linspace(0,1,N+1); % find the range to sketch the x_N
    elseif N==40
        u40=itermeth(A,b,u0,nmax,tol,'G'); % find the solution for N=40
        u40=[0;u40;0]; % add u_0 and u_N
        x40=linspace(0,1,N+1); % find the range to sketch the x_N
    else
        u80=itermeth(A,b,u0,nmax,tol,'G'); % find the solution for N=80
        u80=[0;u80;0]; % add u_0 and u_N
        x80=linspace(0,1,N+1); % find the range to sketch the x_N
    end
end
x1=linspace(0,1);
y=@(x1)pi.^(-2)*sin(pi*x1); % solution to the differential equation
figure
plot(x1,y(x1))
hold on
plot(x5,u5); plot(x10,u10);plot(x20,u20);plot(x40,u40);plot(x80,u80)
legend('y(x)=\pi^{-2}sin(\pi x)','N=5','N=10','N=20','N=40','N=80')
%%
% 
% We can see that the solutions converge to the solution of the
% differential equation as N grows
error_vect=[]; hvect=1./Nvec;
e5=max(abs(transpose(u5)-y(x5))); error_vect(end+1)=e5; % calculate the maximum error for N=5
e10=max(abs(transpose(u10)-y(x10)));error_vect(end+1)=e10; % same but for N=10
e20=max(abs(transpose(u20)-y(x20)));error_vect(end+1)=e20; % N=20
e40=max(abs(transpose(u40)-y(x40)));error_vect(end+1)=e40; % N=40
e80=max(abs(transpose(u80)-y(x80)));error_vect(end+1)=e80; % N=80
figure
loglog(hvect,error_vect) % plot of the error
% find an approximation for the error_vect plot
p=(log(error_vect(end))-log(error_vect(1)))/(log(hvect(end))-log(hvect(1))) % gradient gives 1.989 which is approx. 2
d=log(error_vect(1))-p*log(hvect(1)); % intercept with the axis gives log(c)
C=exp(d);
figure
loglog(hvect,error_vect) 
hold
loglog(hvect,C*hvect.^p) % plot of the error and the approximation of the error
legend('error_vect', 'C*hvect^p')
%%
%
% We can see that the numerical results support the theoretical result.
%%
% 
% (c)
%
% Now consider the RHS of the differential equation to be 1. We have to
% adjust b to get an 1xN-1 vector of ones. As before, we solve the system
% for the different values of N and plot it in the same graph as the
% solution.

for N=[5,10,20,40,80]
    h=1/N;
    u0=zeros(N-1,1);
    A=(2/h^2)*diag(ones(N-1,1))-(1/h^2)*diag(ones(N-2,1),1)...
        -(1/h^2)*diag(ones(N-2,1),-1);
    b2=transpose(ones(1,N-1)); % adjust b
    if N==5
        u52 = itermeth(A,b2,u0,nmax,tol,'G'); % find the solution for N=5
        u52=[0;u52;0]; % add u_0 and u_N
        x52=linspace(0,1,N+1); % find the range to sketch the x_N
    elseif N==10
        u102 = itermeth(A,b2,u0,nmax,tol,'G'); % find the solution for N=10
        u102=[0;u102;0];% add u_0 and u_N
        x102=linspace(0,1,N+1);% find the range to sketch the x_N
    elseif N==20
        u202=itermeth(A,b2,u0,nmax,tol,'G');  % find the solution for N=20
        u202=[0;u202;0];% add u_0 and u_N
        x202=linspace(0,1,N+1); % find the range to sketch the x_N
    elseif N==40
        u402=itermeth(A,b2,u0,nmax,tol,'G'); % find the solution for N=40
        u402=[0;u402;0];% add u_0 and u_N
        x402=linspace(0,1,N+1); % find the range to sketch the x_N
    else
        u802=itermeth(A,b2,u0,nmax,tol,'G'); % find the solution for N=80
        u802=[0;u802;0];% add u_0 and u_N
        x802=linspace(0,1,N+1); % find the range to sketch the x_N

    end
end
x1=linspace(0,1);
y2=@(x)(1./2).*x.*(1-x); % solution to the differential equation
figure
plot(x1,y2(x1))
hold on
plot(x52,u52);plot(x102,u102);plot(x202,u202);plot(x402,u402);plot(x802,u802)
legend('y(x)=x(1-x)/2','N=5','N=10','N=20','N=40','N=80') 
%%
%
% Now, we calculate the maximum error that we have for each one of the
% values of N and plot them. 
error_vector=[]; hvect=1./Nvec;
e52=max(abs(transpose(u52)-y2(x52))); error_vector(end+1)=e52; % calculate the maximum error for N=5
e102=max(abs(transpose(u102)-y2(x102)));error_vector(end+1)=e102; % same but for N=10
e202=max(abs(transpose(u202)-y2(x202)));error_vector(end+1)=e202; % N=20
e402=max(abs(transpose(u402)-y2(x402)));error_vector(end+1)=e402; % N=40
e802=max(abs(transpose(u802)-y2(x802)));error_vector(end+1)=e802; % N=80
figure
loglog(hvect,error_vector) % plot of the error
%%
% We have seen that $e(N)\leq CN^{-2}$ as N gets largerr, where C
% is proportional to the fourth derivative of $y(x)$. 
% We have that $y(x)=\frac{x(x-1)}{2}$, so the fourth derivative will be zero. Giving that
% the error,e(N), as N gets larger will go to zero. Because
% $e(N)\leq 0$ and $e(N)\geq 0$, since we have absolute values, $e(N)=0$ for
% N sufficiently large.