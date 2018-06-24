function [] = rq_main()
% 对比二分搜索和sdp问题的主程序
global A1 A2 b1 c1 c2 m epislon
% A1 = [1,3;3,0];A2= [2,1];A2 = A2'*A2;b1=[0;0];c1=2;c2=1;
% A2=[2,1;1,3]
% A1 = [83.4044   72.0439
%    72.0439   60.46650];A2= [10   25
%    25  62.5 ];c1=50;c2=100;

epislon = 1e-6;
time_binary = []; time_sdpt3 = []; tiem_sedumi = [];
solve_binary = {}; solve_sdpt3 = {}; solve_sedumi = {};
lambda_binary = {}; lambda_sdpt3 = {}; lambda_sedumi = {};
t_binary = {}; t_sdpt3 = {}; t_sedumi = {};
% m=1000;
dimision = [10,100,500,1000,1500,2000];
% dimision = [10];
for k = 1:size(dimision,2)
    m = dimision(k);
    fprintf('m=%d\n',m)
    A1 = sprandn(m,m,0.1)*10;A1 = A1'+A1;
    A2 = sprandn(m,m,0.1)*10;A2 = A2'+A2; 
    A2 = A2 - eigs(A2,1,'SA')*eye(m);
    b1 = rand(m,1)*10;c1 = rand()*10; c2 = rand()*10;
    
%     二分法-
    fprintf('now,binary method is runing \n')
    t1 = clock();
    [solve_binary{k},lambda_binary{k},t_binary{k}] = RQ();
    time_binary(k) = etime(clock,t1);
    
%     cvx调用sdpt3
    fprintf('now,sdpt3 is runing \n')
    t1 = clock();
    [solve_sdpt3{k},lambda_sdpt3{k},t_sdpt3{k}] = cvx_sdp(1);
    time_sdpt3(k) = etime(clock,t1);
%     cvx调用sedumi
    fprintf('now,sedumi is runing \n')
    t1 = clock();
    [solve_sedumi{k},lambda_sedumi{k},t_sedumi{k}] = cvx_sdp(0);
    time_sedumi(k) = etime(clock,t1);
save('RQ_test.mat','time_binary','time_sdpt3','time_sedumi','solve_binary',...
   'solve_sdpt3', 'solve_sedumi','lambda_binary','t_binary','lambda_sdpt3','t_sdpt3',...
   'lambda_sedumi','t_sedumi')

end
save('RQ_test.mat','time_binary','time_sdpt3','time_sedumi','solve_binary',...
   'solve_sdpt3', 'solve_sedumi','lambda_binary','t_binary','lambda_sdpt3','t_sdpt3',...
   'lambda_sedumi','t_sedumi')

end