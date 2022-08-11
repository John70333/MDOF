clear;
n = 3;
m = [1,1,.5].*1e2; c = [5,2,1].*1e3; k = [1,1,2].*1e7; 
J = diag(m);
C = diag(c+[c(2),c(3),0])-diag([c(2),c(3)],1)-diag([c(2),c(3)],-1);
K = diag(k+[k(2),k(3),0])-diag([k(2),k(3)],1)-diag([k(2),k(3)],-1);
x0 = [.1;0;0]; v0 = [0;0;0]; 
f = @(t)[.3;1;-1].*1e5.*sin(160*pi*t);  p = [1,1/80,10];
ti = linspace(0,.5,501);
%%
if length(J)~=n||length(C)~=n||length(K)~=n
    disp('Warning: check matrix dimensions!');
elseif ~issymmetric(J)||~issymmetric(C)||~issymmetric(K)
    disp('Warning: assymmetric matrix!');
else
    [l1,u1,xt1,vt1] = draw(n,J,zeros(n),K,x0,v0,@(t)zeros(n,1),p,ti);
    [l2,u2,xt2,vt2] = draw(n,J,C,K,x0,v0,@(t)zeros(n,1),p,ti);
    [l3,u3,xt3,vt3] = draw(n,J,zeros(n),K,x0,v0,f,p,ti);
    [l4,u4,xt4,vt4] = draw(n,J,C,K,x0,v0,f,p,ti);
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
    xlabel('time');
    ylabel('displacement');
    plot(ti,xt1(i,:),'linewidth',1.5);
    plot(ti,xt2(i,:),'linewidth',1.5);
    legend('undamped','damped','location','ne');
    subplot(2,1,2); hold on; grid on;
    title('Forced vibrations');
    xlabel('time');
    ylabel('displacement');
    plot(ti,xt3(i,:),'linewidth',1.5);
    plot(ti,xt4(i,:),'linewidth',1.5);
    legend('undamped','damped','location','ne');
    saveas(gcf,sprintf('%ddof%d.png',n,i));
end
%%
function [l,u,xt,vt] = draw(n,J,C,K,x0,v0,f,p,ti)
A = [C,J;J,zeros(n)]; B = [K,zeros(n);zeros(n),-J];
y0 = [x0;v0]; q = @(t)[f(t);zeros(n,1)];
[U,l] = eig(-B,A,'v'); u = U(1:n,:)./U(1,:);
z0 = U.'*A*y0; a = diag(U.'*A*U);
if p(1) == 1
    [a0,an,bn,wn] = fss(@(t)U.'*q(t),0,p(2),2*n,p(3));
    g = @(x)-a0./l;
    ai = (-an.*l-bn.*wn)./(l.^2+wn.^2);
    bi = (-bn.*l+an.*wn)./(l.^2+wn.^2);
    for i = 1:p(3)
        g = @(x)g(x)+ai(:,i).*cos(wn(i).*x)+bi(:,i).*sin(wn(i).*x);
    end
    g = @(x)g(x).*exp(-l.*x)-g(0);
elseif p(1) == 0
    g = @(t)integral(@(tj)U.'*q(tj).*exp(-l*tj),0,t,'arrayvalued',true);
end
y = @(t) U*((z0+g(t)).*exp(l*t)./a); yt = zeros(2*n,length(ti));
for i = 1:length(ti)
    yt(:,i) = y(ti(i));
end
xt = yt(1:n,:); vt = yt(n+1:2*n,:);
end
function [a0,a,b,w] = fss(f,l,u,nf,nl)
a0 = 1/(u-l)*integral(f,l,u,'arrayvalued',true);
a = zeros(nf,nl); b = zeros(nf,nl); w = zeros(1,nl);
for n = 1:nl
    w(n) = 2*n*pi/(u-l);
    a(:,n) = 2/(u-l)*integral(@(x)f(x).*cos(w(n)*x),l,u,'arrayvalued',true);
    b(:,n) = 2/(u-l)*integral(@(x)f(x).*sin(w(n)*x),l,u,'arrayvalued',true);
end
end
