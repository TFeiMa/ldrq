function [x0,lambda,g] = RQ()
% clear;clc;
% A1,b1,c1�ֱ�Ϊ�����ϵľ��������
% A2,c2,ΪԼ������ͱ���
% m:������ɾ����ά��
global A1 A2 b1 c1 c2 m epislon

% ������ɵ�һЩ����
% A1 = [1,0;0,0];A2= [0,1];A2 = A2'*A2;b1=[0;0];c1=2;c2=1;
% m = 22;
% rng(1);A1 = rand(m,m)*100;A1 = A1'+A1;
% rng(2);A2 = rand(m,m)*10;A2 = A2'*A2;
% rng(3);b1 = rand(m,1)*100;c1 = 50; c2 = 100;
% ȷ��lanmbda_max
A = [A1,b1;
        b1',c1];
mu = eigs(A,1,'SA');
lambda_min = 0;lambda_max = (c1-mu)/c2;
% epislon = 1e-4;

f = []; g = []; x0 = [];
k = 0; u0 = 10;d_lambda = 2;
% dx = 10; min_gev = rand(n,1);
% tic
while abs(d_lambda) > epislon && abs(u0) > epislon
    k = k + 1;
    lambda(k) = (lambda_min+lambda_max)/2;
    A = [A1+lambda(k)*A2,b1
        b1',c1-lambda(k)*c2];
%     min_ev ������min_e ����ֵ
    [min_ev,min_e] = eigs(A,1,'SA');
%     min_ev = -ev(:,1);
%     [V,D] = eig(A);
%     min_e = min(diag(V));
%     min_ev = V(:,1);
    
    y = min_ev./min_ev(end);
    
    x0 = y(1:end-1);
    f1 = y'*[A1,b1;b1',c1]*y ;
    f2 = 1+norm(x0)^2 ;
    f(k) = (f1)/(f2);
%     ��΢���ж�
%     u(k) = (x0'*A2*x0 - c2)/(f2);
    u0 = (x0'*A2*x0 - c2)/(f2);;
    if u0 > 0
        lambda_min = lambda(k);
    elseif u0 < 0
        lambda_max = lambda(k);
    else
        break
    end
    g(k) = min_e;
    gap = f(k) - g(k);
    d_lambda = lambda_max-lambda_min;
end
% toc
% figure(1)
% subplot(2,2,1)
% plot(1:k,f,'-*',1:k,g,'-^')
% legend('f','g')
% subplot(2,2,2)
% plot(1:k,gap,'-*',1:k,lambda,'-^')
% legend('gap','lambda')
% subplot(2,2,3)
% plot(x0(1,:),x0(2,:),'-*')
% legend('x')
% subplot(2,2,4)
% plot(1:k,u,'-+')
% legend('���ݶ�u')
% if abs(gap) < 1e-2
if x0'*A2*x0 - c2 < 1e-2
    fprintf('��ý�:x0 \n')
    x0
    fprintf('ԭ��������ֵf��%f \n',f(end))
    fprintf('��ż��������ֵg��%f \n',g(end))
    fprintf('x0*A2*x0 - c2 = %f ]\n',x0'*A2*x0 - c2)
    fprintf('lambda��%f \n',lambda(end))
    
else
    fprintf('���ַ�x0*A2*x0 - c2 = %f ',x0'*A2*x0 - c2)
    fprintf('���ַ���õĽ�Ϊ��%f \n')
    x0
    f1 = [x0;1]'*[A1,b1;b1',c1]*[x0;1] ;  f2 = 1+norm(x0)^2 ;
    f = (f1)/(f2); 
    fprintf('���ַ����ԭ��������ֵf��%f \n',f)
    fprintf('ʹ��nesterov�����ݶȷ����: \n')
    lambda_star=lambda(end); value=g(end);
    x0 = nesterov_grad_constraint(lambda_star,value)
    f1 = [x0;1]'*[A1,b1;b1',c1]*[x0;1] ;  f2 = 1+norm(x0)^2 ;
    f = (f1)/(f2); 
    fprintf('ԭ��������ֵf��%f \n',f)
    fprintf('��ż��������ֵg��%f \n',g(end))
    fprintf('nesterov x0*A2*x0 - c2 = %f ',x0'*A2*x0 - c2)
end
lambda = lambda(end);
g = g(end);
end