
%%% generate data %%%%
k=1;k_star=1;a=1;a_star=1.2;
x1_exp = -(-a*k + a_star*k_star - k_star*x3)./(k + 2*k_star);
x2_exp = -(a*k - a_star*k_star - k*x3 - k_star*x3)./(k + 2*k_star);
x0=0;x3=3;
% x1_exp = 28/30; x2_exp = 62/30; x0=0; x3=3;


%%%% calculate the positions of atoms given the parameters %%%%
% intial guess
k=0.9; k_star = 1.1; a=1.1; a_star = 1.1;
P_old=[a;a_star;k;k_star];

%initialize
loss_arr=[];

for ii=1:10000
% solve for positions
x1 = -(-a*k + a_star*k_star - k_star*x3)./(k + 2*k_star); 
x2 = -(a*k - a_star*k_star - k*x3 - k_star*x3)./(k + 2*k_star);

% loss function
ls1 = x1-x1_exp; ls2 = x2 - x2_exp;
loss = 1/2*((ls1)^2 + (ls2)^2);

%% adjoint equation %%
Kadj = zeros(2,2); 
% right hand side
rhs = [ls1;ls2];

Kadj(1,1) = k+k_star;
Kadj(1,2) = -k_star;
Kadj(2,1) = -k_star;
Kadj(2,2) = k_star + k;

lamda = Kadj'\rhs;

%% derivative of loss function %%

ld_a = (lamda(1)-lamda(2))*k;
ld_astar =-(lamda(1)-lamda(2))*k_star;
ld_k = -(lamda(1)*(x1 - x0 - a) - lamda(2)*(x3- x2 - a));
ld_kstar = -(lamda(1)*(x2 - x1 - a_star) + lamda(2)*(x2 - x1 - a_star));

ld=[ld_a;ld_astar;ld_k;ld_kstar];

%%%%% update the parameters %%%%
alpha=-1e-2;
P_old = [a;a_star;k;k_star];
P_new = P_old + alpha*ld;
a=P_new(1);
a_star=P_new(2);
k=P_new(3);
k_star = P_new(4);

%% store the loss funtion at every iteration %%%
loss_arr=[loss_arr;loss];

end

%% plot loss function %%
figure(1)
hold on
plot(loss_arr,'ok')
