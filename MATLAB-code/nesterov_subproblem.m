function [x] = nesterov_subproblem(dh,y)
global A1 A2 b1 c1 c2 epislon rho

% ��ʼ��
t0 = dh'*y - 1/(4*rho)*norm(dh)^2;
eta_l = 0; eta_r = (rho*norm(y)^2 - t0)/c2;
d_eta = eta_r - eta_l; m = size(A2,1); cx = 6;k=0;
while d_eta > epislon && abs(cx) > epislon
   eta = (eta_l + eta_r)/2;
%        ��������һ����Լ����ǿ͹���ι滮������������һ����
   H = rho*eye(m) + eta*A2; b = dh - 2*rho*y;
%    ʹ�ü����ݶȷ�������Լ��������
    k=k+1;
   x = nesterov_grad_unconstraint(H,b);
%        Լ��������ֵ��Ҳ���Ǵ�΢��
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