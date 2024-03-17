
%%% LJ potential generate data %%%%%
A=1;A_star=1.1;B=2;B_star=2;
P=[A,A_star,B,B_star];
xg=[1,2];
fun= @(x)LJ(x,P);
x= fsolve(fun,xg);
x1_exp = x(1); x2_exp = x(2); 
x0=0; x3=3;

%%% estimate the parameters %%%%
% initial guess
A=0.9;A_star=1.0;B=2.0;B_star=2.0;
P_old=[A,A_star,B,B_star];
xg = [1.1,2.1];

%initialize
loss_arr=[];
norm_ld = 1.0;
alpha=-1e-0;
% alpha = -0.01;
options = optimset('Display','off');

%% loop starts %%
for ii=1:20000
fun=@(x)LJ(x,[A,A_star,B,B_star]);
x= fsolve(fun,xg);
x1 = x(1); x2 = x(2);

% loss function
ls1 = x1-x1_exp; ls2 = x2 - x2_exp;
loss = 1/2*((ls1)^2 + (ls2)^2);

%% adjoint equation %%
Kadj = zeros(2,2); 
% right hand side
rhs = [ls1;ls2];

Kadj(1,1) = 156*A./(x1-x0)^14 - 42*B./(x1 - x0)^8 + 156*A_star./(x2-x1)^14 - 42*B_star./(x2 - x1)^8;
Kadj(1,2) = -156*A_star./(x2-x1)^14 + 42*B_star./(x2-x1)^8;
Kadj(2,1) = -156*A_star./(x2-x1)^14 + 42*B_star./(x2-x1)^8;
Kadj(2,2) = 156*A_star./(x2-x1)^14 - 42*B_star./(x2 - x1)^8 + 156*A./(x3-x2)^14 - 42*B./(x3 - x2)^8;

lamda = Kadj'\rhs;

%% derivative of loss function %%

ld_A = lamda(1)*12./(x1-x0)^13 - lamda(2)*12./(x3-x2)^13;
ld_Astar = -lamda(1)*12./(x2-x1)^13 + lamda(2)*12./(x2-x1)^13;
ld_B = -lamda(1)*6./(x1-x0)^7 + lamda(2)*6./(x3-x2)^7;
ld_Bstar = lamda(1)*6./(x2-x1)^7 - lamda(2)*6./(x2-x1)^7;

ld=[ld_A;ld_Astar;ld_B;ld_Bstar];

%% initializing things for MMA update
% if ii == 1
%     norm_ld = norm(ld);
% end
% ld = ld/norm_ld;


%%%%% update the parameters %%%%
%     if ii >= 2
%         if loss/old_loss < 0.95
%             alpha = alpha*10.0;
%         elseif loss/old_loss > 1.05
%             alpha = alpha/10.0;
%         end
%     end

P_old=[A;A_star;B;B_star];
P_new = P_old + alpha*ld;
A=P_new(1);
A_star=P_new(2);
B=P_new(3);
B_star = P_new(4);

%% store the loss funtion at every iteration %%%
loss_arr=[loss_arr;loss];
old_loss = loss;

end

%% plot loss function %%
figure(2)
hold on
% plot(loss_arr,'ok')
semilogy(loss_arr(1:20000)/loss_arr(1))


%%% computing the equilibrium equations by solving the non-linear equations for LJ potential

function F = LJ(x,P)
x0=0;x3=3;
A=P(1);A_star =P(2);B=P(3);B_star =P(4);
F(1) = -12*A/(x(1)-x0)^13 + 6*B/(x(1)-x0)^7 + 12*A_star/(x(2)-x(1))^13 - 6*B_star/(x(2)-x(1))^7;
F(2) = -12*A/(x(2)-x(1))^13 + 6*B/(x(2)-x(1))^7 + 12*A_star/(x3-x(2))^13 - 6*B_star/(x3-x(2))^7;

end