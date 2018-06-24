function [x] = nesterov_subproblem(dh,y)
global A1 A2 b1 c1 c2 epislon rho

% 初始化
t0 = dh'*y - 1/(4*rho)*norm(dh)^2;
eta_l = 0; eta_r = (rho*norm(y)^2 - t0)/c2;
d_eta = eta_r - eta_l; m = size(A2,1); cx = 6;k=0;
while d_eta > epislon && abs(cx) > epislon
   eta = (eta_l + eta_r)/2;
%        子问题是一个无约束的强凸二次规划，计算二次项和一次项
   H = rho*eye(m) + eta*A2; b = dh - 2*rho*y;
%    使用加速梯度法计算无约束子问题
    k=k+1;
   x = nesterov_grad_unconstraint(H,b);
%        约束函数的值，也就是次微分
   cx = x'*A2*x - c2;
   if cx >0
       eta_l = eta;
   elseif cx<0
       eta_r = eta;
   end
   d_eta = eta_r - eta_l;
end
% x = 2*x;
end