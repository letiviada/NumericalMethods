%% Leticia Viada Campos
% Mathematical Methods Computational homework 3
%
clear all, close all, clc, format long, format compact
%%
%
% Exercise 1
%(a) 
f=@(t,y)1-y.^2;
y0=0;
y_ex=@(t)(exp(2*t)-1)./(exp(2*t)+1);
tmax=20;t=tmax*(0:0.01:1); ytmax=y_ex(tmax); y=y_ex(t);
%%
% We use the forward Euler method to compute the numerical solution of the
% plot
%%

[t_fe19, u_fe19] = feuler(f,[0,tmax],y0,19);
[t_fe21, u_fe21] = feuler(f,[0,tmax],y0,21);
[t_fe40, u_fe40] = feuler(f,[0,tmax],y0,40);
[t_fe80, u_fe80] = feuler(f,[0,tmax],y0,80);
[t_fe160, u_fe160] = feuler(f,[0,tmax],y0,160); % t values for plotting exact solution
figure
plot(t,y,t_fe19,u_fe19,t_fe21,u_fe21,t_fe40,u_fe40,t_fe80,u_fe80,t_fe160,u_fe160);
legend('Exact solution','Forward Euler N=19','Forward Euler N=21','Forward Euler N=40','Forward Euler N=80','Forward Euler N=160');
title('y^\prime(t) =1-y^2');
%%
% We compute the maximal error of each method and store it in err_fe
err_fe=[];
err_fe(end+1) = max(abs(y_ex(t_fe19) - u_fe19)); err_fe(end+1) = max(abs(y_ex(t_fe21) - u_fe21));
err_fe(end+1) = max(abs(y_ex(t_fe40) - u_fe40)); err_fe(end+1) = max(abs(y_ex(t_fe80) - u_fe80));
err_fe(end+1) = max(abs(y_ex(t_fe160) - u_fe160));
err_fe
%%
% Similarly we compute the error at the end of the interval
err_feend=[];
err_feend(end+1)=abs(ytmax - u_fe19(end)); err_feend(end+1)=abs(ytmax - u_fe21(end));
err_feend(end+1)=abs(ytmax - u_fe40(end)); err_feend(end+1)=abs(ytmax - u_fe80(end));
err_feend(end+1)=abs(ytmax - u_fe160(end));
err_feend
%%
%
% (b)
[t_he19, u_he19] = heun(f,[0,tmax],y0,19);
[t_he21, u_he21] = heun(f,[0,tmax],y0,21);
[t_he40, u_he40] = heun(f,[0,tmax],y0,40);
[t_he80, u_he80] = heun(f,[0,tmax],y0,80);
[t_he160, u_he160] = heun(f,[0,tmax],y0,160);
figure
plot(t,y,t_he19,u_he19,t_he21,u_he21,t_he40,u_he40,t_he80,u_he80,t_he160,u_he160);
legend('Exact solution','Heun N=19','Heun N=21','Heun N=40','Heun N=80','Heun N=160');
title('y^\prime(t) =1-y^2');

err_he=[];
err_he(end+1) = max(abs(y_ex(t_he19) - u_he19)); err_he(end+1) = max(abs(y_ex(t_he21) - u_he21));
err_he(end+1) = max(abs(y_ex(t_he40) - u_he40)); err_he(end+1) = max(abs(y_ex(t_he80) - u_he80));
err_he(end+1) = max(abs(y_ex(t_he160) - u_he160));
err_he

err_heend=[];
err_heend(end+1)=abs(ytmax - u_he19(end)); err_heend(end+1)=abs(ytmax - u_he21(end));
err_heend(end+1)=abs(ytmax - u_he40(end)); err_heend(end+1)=abs(ytmax - u_he80(end));
err_heend(end+1)=abs(ytmax - u_he160(end));
err_heend
%%
%
% (c)
Nvec=[19,21,40,80,160];
loglog(20./Nvec,err_fe,20./Nvec, err_he)
legend('err_{fe}','err_{he}')
%% 
% Using rules of logs, we get that the error is approximately E(N)=CN^p.
% And so we find an approximation for the err_fe and err_he plots.
p_fe=(log(err_fe(end))-log(err_fe(1)))/(log(Nvec(end))-log(Nvec(1))) % gradient gives -1.13625 which is approx. -1
d_fe=log(err_fe(1))-p_fe*log(Nvec(1)); % intercept with the axis gives log(c)
C_fe=exp(d_fe);
% find an approximation for the err_he plot
p_he=(log(err_he(end))-log(err_he(1)))/(log(Nvec(end))-log(Nvec(1))) % gradient gives -2.50204 which is approx. 2
d_he=log(err_he(1))-p_he*log(Nvec(1)); % intercept with the axis gives log(c)
C_he=exp(d_he);
%%
%
% Therefore, we can see that the numerical results coincide with the
% theoretical results as forward Euler's method is of first order of convergence and Heun's
% method is of second order of convergence. 

%%
% 
% (d)
% For forward Euler's method we have that for N=19, the approximation does
% not tend to 1 as $t\rightarrow \infty$. Looking at err_fe,end we can see
% that as N gets larger, the error tends to zero. In fact, for N=19
% err_fe,end gives approximately 0.211 and for N=21,40,80,160 err_fe,end=0.
% Similarly, for Heun's method we have that for N=19 the approximation
% does not tend to 1 either with err_he,end=0.36336. For N=21 it tends to a value between the N=19 case
% and N=40,80,160; with err_he,end approximately 0.03596. For N=40,80,160
% we get that err_he,end=0.
%
%%
% We define h=20/N. The method is stable for h<1, in our numerical results
% this is satisfied for N=21,40,80,160 but not for N=19.
% For forward Euler's method the lack of stability is shown through the numerical
% approximation increasingly oscillating, as seen in the first graph.
% For Heun's method we see instability for N=19 as the final point from our calculation is
% approximately 0.6366 and the numerical approximation will tend to a value
% close to this. 
