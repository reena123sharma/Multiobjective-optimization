clc;  clear all;
load('PP_L.mat');   load('PP_solar.mat');   load('PP_W.mat');
PP_l=P_l;   PP_solar=P_solar;   PP_W=P_W;
clear P_l P_solar P_W
%initiate loop************************************************
ic=zeros(8,25);
zmin_avg=zeros(24,12);
for itr=1:24
    P_l=PP_l(itr:itr+24,1); P_solar=PP_solar(itr:itr+24,1); P_W=PP_W(itr:itr+24,1);
    %till here P_l, P_solar, P_W [48x1]
    Ts =1;          %sampling time
    m =24;          %iteration
    dt =Ts*24;
    Ef = 0.5;
    e1 = 0.9;
    P_dgr = 40;
    e2 = 0.95;
    A_dg = 0.246 ;
    Cost_grid = 2;
    C_fuel = 65;
    %global P_load;
    global i;
    % P_load = 60.*ones(m,1);
    P_buy = zeros(m,1);
    P_sell = zeros(m,1);
    nvars = 8;
    x = zeros(3,8);

    %iteration-1 start---------------------------------------------
    i=1;
    fitnessfcn = @(x) ga_multiobjective1(x,P_l,i,itr);
    rng default % For reproducibility
    % lb = [0 0 0 0 0 0 0 0]; % Lower bound
    % ub = [inf inf 0 inf inf inf inf inf];
    % lb = []; % Lower bound
    % ub = []; % Upper bound
    %options = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto);
    options = optimoptions(@gamultiobj,'PopulationSize',100);
    options = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto,'PopulationSize',100);
    options = optimoptions(options,'Tolfun',1e-3,'MaxstallGenerations',50); 
        objectivesToPlot = [1,2,3]; % or whatever two objectives you want
         plotfn = @(options,state,flag)gaplotpareto(options,state,flag,objectivesToPlot);
         options = gaoptimset('PlotFcns',plotfn);

    %change for each loop**********************************************
    for j =1:3
        P_dg1(i,j) = ic(1,itr);
        P_dg2(i,j) = ic(2,itr);
        P_fc(i,j) = ic(3,itr);
        P_grid1(i,j) = ic(4,itr);
        P_grid2(i,j) = ic(5,itr);
        P_s1(i,j) = ic(6,itr);
        P_s2(i,j)=ic(7,itr);
        P_fcmax(i,j)= ic(8,itr);
        index(i,j) =1;
        P_dg(i,j) = P_dg1(i,j)+P_dg2(i,j);
        P_grid(i,j) = P_grid1(i,j)+P_grid2(i,j);
    end
    fprintf("1st section has been successfully executed \n")
    %*******************************************************************

     Y1 = [];
     i=2;
     h=1;
     A = [1,1,0,0,0,0,0,0];  %;0,0,-1,0,0,0,0,0];
     b = P_dgr;   %;0];
     lb = [0 0 0 0 0 0 0 0]; % Lower bound
     ub = [inf inf 0 inf inf inf inf inf];
     Aeq =[0,-1,1,0,-1,0,-1,1; 0,1,1,1,0,1,0,0; 0,0,0,0,0,1,1,0];
     beq =[P_fcmax(i-1) P_l(i) P_solar(i)+ P_W(i)] ;  
     [x,Fval,exitFlag,Output] = gamultiobj(fitnessfcn,nvars,A, b,Aeq,beq,lb,ub,options);
        am=1;
     if exitFlag==1
        y = [];
        ind=h.*ones(length(x),1);

    for j=1:length(x)
       y(j,1) = (sum(Cost_grid*(x(j,4) + x(j,5))) +C_fuel*A_dg*sum(x(j,1))+C_fuel*A_dg*sum(x(j,2)));
       y(j,2) = 1 -sum((P_solar(i) + x(j,3)+ P_W(i))/sum(P_l(i)));
       y(j,3) = sum(e1*P_dgr + e2*(x(j,1)+x(j,2)))*dt*Ef;
    end
     Y = [x,y,ind];


        if am==1
            Y1=Y;
            am=am+1;
        else
            Y1=cat(1,Y1,Y);
            am=am+1;
        end
     end
    clear am

    if  isempty (Y1)
        disp('empty')
       else
           disp('non-epmty')
    end
    sol =sortNB (Y1);
    for j =1:3
        P_dg1(i,j) = sol(j,1);
        P_dg2(i,j) = sol(j,2);
        P_fc(i,j) = sol(j,3);
        P_grid1(i,j) = sol(j,4);
        P_grid2(i,j) = sol(j,5);
        P_s1(i,j) =sol(j,6);
        P_s2(i,j) = sol(j,7);
        P_fcmax(i,j) =sol(j,8);
        index(i,j) =j;
        P_dg(i,j) = (P_dg1(i,j) + P_dg2(i,j));
        P_grid(i,j) = P_grid1(i,j) + P_grid2(i,j);
    end

     fprintf("2nd section has been successfully executed \n")
     zmin =zeros(25 ,12);
     zmin(2,:)=(Y1(1,:));
     for i=3:25
    % lb = [0 0 0 0 0 0 0 0]; % Lower bound
    % ub = [1300 1300 1300 1300 1300 1300 1300 1300];
     lb = [0 0 0 0 0 0 0 0]; % Lower bound
     ub = [inf inf P_fcmax(i-1) inf inf inf inf inf];
     h = 1;
     fprintf(num2str(i))
     fprintf('\n')
     Y1 = [];
     for j=1:3
     A = [1,1,0,0,0,0,0,0];  %;0,0,-1,0,0,0,0,0];
     b = P_dgr ; %;0]; 
     Aeq =[0,-1,1,0,-1,0,-1,1; 0,1,1,1,0,1,0,0; 0,0,0,0,0,1,1,0];
     beq =[P_fcmax(i-1) P_l(i) P_solar(i)+ P_W(i)] ;  
      x = [];
      [x,Fval,exitFlag,Output] = gamultiobj(fitnessfcn,nvars,A, b,Aeq,beq,lb,ub,options);
     if exitFlag==1
            y = [];
            Y = [];
            ind=h.*ones(length(x),1);
            h = h+1;
     for p=1:length(x)
       y(p,1) =  (sum(Cost_grid*(x(p,4) + x(p,5))) +C_fuel*A_dg*sum(x(p,1))+C_fuel*A_dg*sum(x(p,2)));
       y(p,2) = 1 - sum((P_solar(i) + x(p,3)+ P_W(i))/sum(P_l(i)));
       y(p,3) =sum(e1*P_dgr + e2*(x(p,1)+x(p,2)))*dt*Ef;
     end
         Y = [x,y,ind];
    if i==25
        if j==1
            par=Y;
        end
    end
        %Y1 = zeros(105,9);

        if j==1
            Y1=Y;
        else
            Y1=cat(1,Y1,Y);
        end
        %   Y1(:,5) =[];
     end
     end

        if  isempty (Y1)
            disp('empty')
           else
               disp('non-epmty')
        end

     sol = sortNB(Y1);
     for j =1:3
        P_dg1(i,j) = sol(j,1);
        P_dg2(i,j) = sol(j,2);
        P_fc(i,j) = sol(j,3);
        P_grid1(i,j) = sol(j,4);
        P_grid2(i,j) = sol(j,5);
        P_s1(i,j) =sol(j,6);
        P_s2(i,j) = sol(j,7);
        P_fcmax(i,j) =sol (j,8);
        index(i,j) =sol(j,12);
        P_dg(i,j) = (P_dg1(i,j) + P_dg2(i,j));
        P_grid(i,j) = P_grid1(i,j) + P_grid2(i,j);
     end 
     Y1(:,12) =[];
     [N, m] = size(Y1);
        clear m
        % Initialize the Rank number to 1.
        Rank = 1;

        % There is nothing to this assignment, used only to manipulate easily in
        % MATLAB.

         F(Rank).f = [];
        indv = [];


        for l = 1 : N 

            % Number of individuals that dominate this individual
            indv(l).n = 0;
            % Individuals which this individual dominate
            indv(l).p = [];
            for r = 1 : N
                dom_less = 0;
                dom_equal = 0;
                dom_more = 0;
                for k = 1 : 3
                    if (Y1(l,8+ k) < Y1(r,8+ k))
                        dom_less = dom_less + 1;
                    elseif (Y1(l,8 + k) == Y1(r,8 + k))
                        dom_equal = dom_equal + 1;
                    else
                        dom_more = dom_more + 1;
                    end
                end
                if dom_less == 0 && dom_equal ~= 3
                    indv(l).n = indv(l).n + 1;
                elseif dom_more == 0 && dom_equal ~= 3
                    indv(l).p = [indv(l).p r];
                end
            end 

            if indv(l).n == 0
                Y1(l,3+ 8 + 1) = 1;
                 F(Rank).f = [F(Rank).f l];
            end

        end

        % calculate the subsequent fronts%%%%%%%%%%%%%%%%%%%%%%%%%%
        while ~isempty(F(Rank).f)
           Q = [];
           for l = 1 : length(F(Rank).f)
               if ~isempty(indv(F(Rank).f(l)).p)
                    for j = 1 : length(indv(F(Rank).f(l)).p)
                        indv(indv(F(Rank).f(l)).p(j)).n = ...
                            indv(indv(F(Rank).f(l)).p(j)).n - 1;
                        if indv(indv(F(Rank).f(l)).p(j)).n == 0
                            Y1(indv(F(Rank).f(l)).p(j),3 + 8 +1) = ...
                                Rank + 1;
                            Q = [Q indv(F(Rank).f(l)).p(j)];
                        end
                    end
               end
           end
           Rank =  Rank + 1;
           F(Rank).f = Q;
        end
    %     a=F(Rank).f;
    %     if length(a)<=2
    %         a(1,length(a)+1)=a(1,length(a));
    %     end
    %     F(Rank).f=a;
    %     clear a

        [~,IOR] = sort(Y1(:,3 +8+1));
        for l = 1 : length(IOR)
            sorted_based_on_Rank(l,:) = Y1(IOR(l),:);
        end
        current_index = 0;

        %crowding distance for each individual in each Rank%%%%%%%%%%%%%%%%%%%
    for Rank = 1 : (length(F) - 1)
        %     objective = [];
             distance = 0;
              w=[];
    %        if i==3
    %            w=zmin(i-1,:);
    %            w(:,12)=[];
    %        else
    %            clear aab1 aab2
    %            [aab1,aab2] = min(zmin(3:i,11));
    %            w=zmin(aab2,:);
    %           w(:,12)=[];
    %        end
               w=zmin(i-1,:);
    %            w(:,12)=[];
            previous_index = current_index + 1;
            a=length(F(Rank).f);    %********correction****************
            if a==0
                a=3;
            end
            for l = 1 : a
              w(l,:) = sorted_based_on_Rank(current_index+ l ,:);
            end
             current_index = current_index + l ; %%%%%%%%%%%%% k %%%%%%%%%%%%%%
    %%%%%%%%   Sort each individual based on the objective%%%%%%%%%%
               SBO= [];
            for l = 1:3
              [SBO,IO] = sort(w(:,8+l));
                   SBO = [];
                for j = 1 : length(IO)
                    SBO(j,:) = w(IO(j),:);
                end
                f_max = SBO(length(IO), 8 + l);
                f_min = SBO(1, 8 + l);
                 w(IO(length(IO)),3 + 8 + 1 + l) = inf;
                 w(IO(1),3 + 8+ 1  + l) = inf;
                 for j = 2 : length(IO)-1
                    next_obj  = SBO(j+1,8 + l);
                    previous_obj  = SBO(j-1 ,8 +l);
                    if (f_max - f_min == 0)
                        w(IO(j),3 + 8 + 1 +l) = inf;
                    else
                        w(IO(j),3 + 8 + 1 +l) = ...
                             (next_obj - previous_obj)/(f_max - f_min);
                    end
                 end

            end

            distance = [];
            distance(:,1) = zeros(a,1); %*********correction**************
            for o = 1 : 3
                distance(:,1) = distance(:,1) + w(:,3 + 8 +1+ o);
            end
            w(:,3 + 8 + 2) = distance;
            w = w(:,1 : 3 + 8 + 2);
             z(previous_index:current_index,:) = w;


                %********************\\\///\\\///*****************%
             clear previous_index
             current_index=0;


    end
    clear a
    %[ab1,ab2] = min(z(:,11));
    % zmin(i,:)=z(ab2,:);
    g1=Y1(:,9:11);
    
    tempp(1,:)=min(g1);
    tempp(2,:)=max(g1);
    for wq=1:size(g1,1)
        g1(wq,1)=(g1(wq,1)-tempp(1,1))/(tempp(2,1)-tempp(1,1));
        g1(wq,2)=(g1(wq,2)-tempp(1,2))/(tempp(2,2)-tempp(1,2));
        g1(wq,3)=(g1(wq,3)-tempp(1,3))/(tempp(2,3)-tempp(1,3));
    end
    
    for ii=1:size(g1,1)
    gg(ii,1)= (sum(g1(ii,:))/3);
    end
    [am1 am2]=find(gg(:,1)==min(gg));
    zmin(i,:)=Y1(am1(1,1),:);
    clear ii g1 gg am2
     end
     
     if itr==1
         zmin_avg=zmin(2:25,:);
     else         
         ac1=0; ac2=1;
         zzmin=zeros(24,12);
         for bb=1:24
             if (itr+ac1)<=24
                 zzmin(itr+ac1,:)=zmin(bb+1,:);
                 ac1=ac1+1;
             else
                 zzmin(ac2,:)=zmin(bb+1,:);
                 ac2=ac2+1;
             end                 
         end
         clear ac1 ac2
         for aa=1:24
             zmin_avg(aa,:)=(zmin_avg(aa,:)+zzmin(aa,:))/2;
         end
         clear zzmin
     end     
     ic(:,itr+1)=zmin_avg(itr+1,1:8)';
     save(['results_initial_hr' num2str(itr)])
     clearvars -except PP_l PP_solar PP_W itr ic zmin_avg
end
save('ic.mat','ic');
save('zmin_avg.mat','zmin_avg');

 
 