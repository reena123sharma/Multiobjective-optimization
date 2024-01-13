function y = ga_multiobjective1(x,P_l,i,itr)
% Initialize for three objectives 
load('PP_L.mat');   load('PP_solar.mat');   load('PP_W.mat');
PP_l=P_l;   PP_solar=P_solar;   PP_W=P_W;
clear P_l P_solar P_W
P_l=PP_l(itr:itr+24,1); P_solar=PP_solar(itr:itr+24,1); P_W=PP_W(itr:itr+24,1);
% load('P_L.mat');
% pll(:,1)=P_l;   plt(:,1)=Tim;       %****when Tim vector is same for all loads****%
% clear Tim
% load('P_s.mat');
% psl(:,1)=P_solar;   pst(:,1)=Tim;   %****
% clear Tim
% load('P_w2.mat');
% pwl(:,1)=P_W;   pwt(:,1)=Tim;       %****
% clear Tim P_l P_solar P_W
% 
% P_l=zeros(25,1);    P_solar=zeros(25,1);    P_W=zeros(25,1);
% for i=1:25
%     [temp1,temp2]=find(plt(:,1)<=i);
%     P_l(i,1)=pll(temp1(end,1),1);
%     clear temp1 temp2
%     [temp1,temp2]=find(pst(:,1)<=i);
%     P_solar(i,1)=psl(temp1(end,1),1);
%     clear temp1 temp2
%     [temp1,temp2]=find(pwt(:,1)<=i);
%     P_W(i,1)=pwl(temp1(end,1),1);
%     clear temp1 temp2    
% end
% P_solar= 5*P_solar*0.95*0.95;
% P_W= 0.025*P_W;
% P_l = P_l;
% % P_solar= 5*P_solar;
% % P_W= 0.110*P_W;
% % P_l =  P_l;
y = zeros(3,1);
Ts = 1;
dt = Ts*24;
Ef = 0.5;
e1 = 0.9;
A_dg = 0.246 ;
P_dgr = 40;
C_fuel = 65;
e2 = 0.95;
%  global P_load
%P_max = 50;
Cost_grid = 2;
y(1) = (sum(Cost_grid*(x(4) + x(5))) +C_fuel*A_dg*sum(x(1))+C_fuel*A_dg*sum(x(2)));
y(2) = 1 - sum((P_solar(i) + x(3)+ P_W(i))/sum(P_l(i)));
y(3) = sum(e1*P_dgr + e2*(x(1)+x(2)))*dt*Ef;
end

