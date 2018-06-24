function [x1] = nesterov_grad_constraint(lambda1,value1)
global A1 A2 b1 c1 c2 epislon rho m lambda value dh y
% 最优乘子和对偶最优值
lambda = lambda1; value = value1;

% epsilon = 1e-6; m=2;

% A1 = [1,0;0,0];A2= [0,1];A2 = A2'*A2;b1=[0;0];c1=2;c2=1;
% rng(1);A1 = rand(m,m)*100;A1 = A1'+A1;
% rng(2);A2 = rand(m,m)*10;A2 = A2'*A2;
% rng(3);b1 = rand(m,1)*100;c1 = 50; c2 = 100;

% 计算目标函数的二次项和0次项
H = A1 + lambda*A2 -value*eye(m); c = c1 - lambda*c2 - value;
mu = eigs(H,1,'SA'); L = eigs(H,1,'LA');

% eig_v = eig(H);
% mu = min(eig_v); L = max(eig_v);
q = mu/L; rho = L;
% 初始化
x1 = rand(m,1)*10; y = x1;
alpha0 = 0.9; alpha1 = 0.9;
dx = 5;

% x = phr('fun', 'cfun','dfun', 'dcfun',x1);
% y = (x'*A1*x + 2*b1'*x - c1)/(1+norm(x)^2);

while dx > epislon
    dh = 2*H*y + 2*b1;
    x0 = x1;
%    求解子问题
    x1 = nesterov_subproblem(dh,y);    
    alpha0 = alpha1;
    alpha1 = 1/2*(sqrt((alpha0^2 - q)^2 + 4*alpha0^2) - (alpha0^2 - q));
    beta = (alpha0*(1-alpha0))/(alpha0^2 + alpha1);
    dx = norm(x1-x0);
    y = x1 + beta*(x1 - x0);
end

end