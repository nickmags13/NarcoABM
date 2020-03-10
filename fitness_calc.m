function [t_firstmovcorr,t_firstmovpval,deptmovcorr,deptmovpval]=...
    fitness_calc(tfirstmov,deptflowsts,cntryflowsts)


%calculate correlation of first movement
load('tmovref_156.mat')

% tmov=120-linspace(1,120,length(tfirstmov));
[tmov_rho,tmov_pval]=corr(tfirstmov,tmov','rows','complete');
tcheck=(max(tfirstmov) >= 100);
t_firstmovcorr=tmov_rho;
% t_firstmovcorr=tmov_rho*tcheck;
t_firstmovpval=tmov_pval+(1-tcheck);

deptrefvecs{1}=[5500 6300 2250 3450 1050 1375 625 625 1600];    %Peten, 2005-2013
deptrefvecs{2}=[29750 26812 37034 49043 41399 41535];   %Darien, 2009-2014
% deptrefvecs{3}=[4000 886 32782 26225.6 11400 25765 28684 34926];    %Puntarenas, 2007-2014
deptrefvecs{3}=[770 1275 3450 5000 23450 23779 92661 159368 114905 92472 ...
    73836 191770 209586 269160];    %Costa Rica, 2001-2014
deptrefvecs{4}=[3930 12800 4463 6906 63395 109875 219354 214087 107219 55268]; %gracias a dios, 2005-2014
deptrefvecs{5}=[560 6766 7045 17690 9850 3715 6100 25258];  %colon, 2006-2013
deptrefvecs{6}=[5575 775 32750 37650 48038 9400 9117 3970]; %atlantico norte, 2007-2014

tspan=[5 13; 9 14; 1 14; 5 14; 6 13; 7 14];
fl_order=[1 2 4 5 6 7];
deptmovcorr=zeros(1,length(fl_order));
deptmovpval=zeros(1,length(fl_order));
for gg=1:6
    if gg == 2
        [dept_rho,dept_pval]=corr(deptflowsts(fl_order(gg),tspan(gg,1):...
            tspan(gg,2))'+deptflowsts(3,tspan(gg,1):tspan(gg,2))',deptrefvecs{gg}');
    elseif gg == 3
        [dept_rho,dept_pval]=corr(cntryflowsts(fl_order(gg),tspan(gg,1):...
            tspan(gg,2))',deptrefvecs{gg}');
    else
        [dept_rho,dept_pval]=corr(deptflowsts(fl_order(gg),tspan(gg,1):...
            tspan(gg,2))',deptrefvecs{gg}');
    end
    if isnan(dept_rho) == 1
        dept_rho=0;
        dept_pval=1;
    end
    deptmovcorr(gg)=dept_rho;
    deptmovpval(gg)=dept_pval;
end
end