function [x2] = nesterov_grad_unconstraint(H,b)
global A1 A2 b1 c1 c2 epislon
% 初始化
rho1 = eigs(H,1); m = size(H,1);
y1 = zeros(m,1);x1 = zeros(m,1);theta1 = 1;
df = 6;
grad_f = 2*H*y1+b;
% 写得不对，得改
while abs(df) > epislon
    grad_f = 2*H*y1+b;
    x2 = y1 - 1/rho1 * grad_f;
    theta2 = 1/2*(1+sqrt(1+4*theta1^2));
    gama = (1-theta1)/theta2;
    y1 = x2 + (1-theta1)/theta2*(x2-x1); 

    f = x2'*H*x2 + b'*x2;
    df = x2'*H*x2 + b'*x2 - x1'*H*x1 - b'*x1;
    x1 = x2;theta1= theta2;
end

end