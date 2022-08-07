clear;
n = 3;
m = [1,1,.5].*1e2; c = [5,2,1].*1e3; k = [1,1,2].*1e7; 
J = diag(m);
C = diag(c+[c(2),c(3),0])-diag([c(2),c(3)],1)-diag([c(2),c(3)],-1);
K = diag(k+[k(2),k(3),0])-diag([k(2),k(3)],1)-diag([k(2),k(3)],-1);
x0 = [.1;0;0]; v0 = [0;0;0]; 
f = @(t)[.3;1;-1].*1e5.*sin(160*pi*t);
ti = linspace(0,.5,501);
%%
if length(J)~=n||length(C)~=n||length(K)~=n
    disp('Warning: check matrix dimensions!');
elseif ~issymmetric(J)||~issymmetric(C)||~issymmetric(K)
    disp('Warning: assymmetric matrix!');
else
    [l1,u1,xt1,vt1] = draw(n,J,zeros(n),K,x0,v0,@(t)zeros(n,1),ti);
    [l2,u2,xt2,vt2] = draw(n,J,C,K,x0,v0,@(t)zeros(n,1),ti);
    [l3,u3,xt3,vt3] = draw(n,J,zeros(n),K,x0,v0,f,ti);
    [l4,u4,xt4,vt4] = draw(n,J,C,K,x0,v0,f,ti);
end
%%
clc;
disp('undamped eigenvalues:'); disp(l1); 
disp('undamped eigenvectors:'); disp(u1);
disp('damped eigenvalues:'); disp(l2); 
disp('damped eigenvectors:'); disp(u2);
for i = 1:n
    figure('windowstate','maximize');
    subplot(2,1,1); hold on; grid on;
    title('Free vibrations');
    plot(ti,xt1(i,:),'linewidth',1.5);
    plot(ti,xt2(i,:),'linewidth',1.5);
    legend('undamped','damped','location','ne');
    subplot(2,1,2); hold on; grid on;
    title('Forced vibrations');
    plot(ti,xt3(i,:),'linewidth',1.5);
    plot(ti,xt4(i,:),'linewidth',1.5);
    legend('undamped','damped','location','ne');
end
%%
function [l,u,xt,vt] = draw(n,J,C,K,x0,v0,f,ti)
A = [C,J;J,zeros(n)]; B = [K,zeros(n);zeros(n),-J]; 
y0 = [x0;v0]; q = @(t)[f(t);zeros(n,1)]; 
[U,l] = eig(-B,A,'v');  u = U(1:n,:)./U(1,:);
z0 = U.'*A*y0; a = diag(U.'*A*U);
g = @(t) integral(@(tj)U.'*q(tj).*exp(-l*tj),0,t,'arrayvalued',true);
y = @(t) U*((z0+g(t)).*exp(l*t)./a); yt = zeros(2*n,length(ti));
for i = 1:length(ti)
    yt(:,i) = y(ti(i));
end
xt = yt(1:n,:); vt = yt(n+1:2*n,:);
end