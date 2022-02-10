%%%%%%%%%%% Batch Run Results %%%%%%%%%%%%%%%%%%%%%%%

% cd D:\Narcoscapes\Batch_ABM_Interdiction_Results\batch_pint
% load batch_pint_results

TSTART=1;
TMAX=180;
MRUNS=30;
ERUNS=11;
ndto=2;

BRUNS=24;   %batch runs of force package levels [batch, mrun, erun]
subbatchind=[reshape(repmat(1:MRUNS,1,ERUNS),MRUNS*ERUNS,1) ...
    reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1)];
batchind=[reshape(repmat(1:BRUNS,ERUNS*MRUNS,1),BRUNS*MRUNS*ERUNS,1) ...
    repmat(subbatchind,BRUNS,1)];
brunlist=1:BRUNS;
frcpkgs=[18 19 20 21 22 23 19 20 21 22 23 24 20 21 22 23 24 25 21 22 23 ...
    24 25 26];
forcepkgs=[1 1; 2 1; 3 1; 4 1; 5 1; 6 1; 1 2; 2 2; 3 2; 4 2; 5 2; 6 2; 1 3; ...
    2 3; 3 3; 4 3; 5 3; 6 3; 1 4; 2 4; 3 4; 4 4; 5 4; 6 4]; %[Pac Carib]
fplevels=sortrows([brunlist' frcpkgs'],2);

MDN_INTRT=zeros(BRUNS,3);   % lowest, baseline, highest probability of interception
QNT_INTRT=zeros(2,3,BRUNS);
mdn_intrt=zeros(BRUNS,ERUNS);
qnt_intrt=zeros(5,ERUNS,BRUNS);
AVG_SL=zeros(BRUNS,3);
STD_SL=zeros(BRUNS,3);
avg_sl=zeros(BRUNS,ERUNS);
std_sl=zeros(BRUNS,ERUNS);
MDN_DLVR=zeros(BRUNS,3);
QNT_DLVR=zeros(2,3,BRUNS);
mdn_dlvr=zeros(BRUNS,ERUNS);
qnt_dlvr=zeros(5,ERUNS,BRUNS);
AVG_SLVAL=zeros(BRUNS,3);
STD_SLVAL=zeros(BRUNS,3);
avg_slval=zeros(BRUNS,ERUNS);
std_slval=zeros(BRUNS,ERUNS);
MDN_RETURN=zeros(BRUNS,3);
QNT_RETURN=zeros(2,3,BRUNS);
mdn_return=zeros(BRUNS,ERUNS);
qnt_return=zeros(5,ERUNS,BRUNS);
MDN_COST=zeros(BRUNS,3);
QNT_COST=zeros(2,3,BRUNS);
mdn_cost=zeros(BRUNS,ERUNS);
qnt_cost=zeros(5,ERUNS,BRUNS);
MDN_DIST=zeros(BRUNS,3);
QNT_DIST=zeros(2,3,BRUNS);
mdn_dist=zeros(BRUNS,ERUNS);
qnt_dist=zeros(5,ERUNS,BRUNS);
MDN_PRIMDIST=zeros(BRUNS,3);
QNT_PRIMDIST=zeros(2,3,BRUNS);
mdn_primdist=zeros(BRUNS,ERUNS);
qnt_primdist=zeros(5,ERUNS,BRUNS);
MDN_KSVAL=zeros(BRUNS,3);
QNT_KSVAL=zeros(2,3,BRUNS);
mdn_ksval=zeros(BRUNS,ERUNS);
qnt_ksval=zeros(5,ERUNS,BRUNS);
MDN_KSSIG=zeros(BRUNS,3);
QNT_KSSIG=zeros(2,3,BRUNS);
mdn_kssig=zeros(BRUNS,ERUNS);
qnt_kssig=zeros(5,ERUNS,BRUNS);

%%% Binning for intra-run variability analysis
hist_intrt=zeros(10,ERUNS,BRUNS);
intedges=linspace(min(min(min(totintrt))),max(max(max(totintrt))),11);


sub_slval=zeros(MRUNS,1);
sub_rttime=zeros(MRUNS,1);
for br=1:BRUNS
    for er=1:ERUNS
        qnt_intrt(:,er,br)=quantile(totintrt(:,er,br),[0.025,0.25,0.5,0.75,0.975])';
        mdn_intrt(br,er)=qnt_intrt(3,er,br);
        avg_sl(br,er)=mean(totsl(:,er,br));
        std_sl(:,er,br)=std(totsl(:,er,br));
        qnt_dlvr(:,er,br)=quantile(totdlvr(:,er,br),[0.025,0.25,0.5,0.75,0.975])';
        mdn_dlvr(br,er)=qnt_dlvr(3,er,br);
        
        %%% Intra-run variability
        hist_intrt(:,er,br)=histcounts(totintrt(:,er,br),intedges);
        
        ind=find(batchind(:,1) == br & batchind(:,3) == er);
        for mr=1:MRUNS
            sub_slval(mr)=sum(cell2mat(slmetrics(ind(mr)).slval));
            qnt_return(:,er,br)=quantile(slmetrics(ind(mr)).returntime,[0.025,0.25,0.5,0.75,0.975])';
            qnt_cost(:,er,br)=quantile(slmetrics(ind(mr)).costroutes,[0.025,0.25,0.5,0.75,0.975])';
            qnt_dist(:,er,br)=quantile(slmetrics(ind(mr)).distroutes,[0.025,0.25,0.5,0.75,0.975])';
            qnt_primdist(:,er,br)=quantile(slmetrics(ind(mr)).primdisplace,[0.025,0.25,0.5,0.75,0.975])';
            qnt_ksval(:,er,br)=quantile(slmetrics(ind(mr)).ksnetwork{:},[0.025,0.25,0.5,0.75,0.975])';
            qnt_kssig(:,er,br)=quantile(slmetrics(ind(mr)).kspval{:},[0.025,0.25,0.5,0.75,0.975])';
        end
        mdn_return(br,er)=qnt_return(3,er,br);
        mdn_cost(br,er)=qnt_cost(3,er,br);
        mdn_dist(br,er)=qnt_dist(3,er,br);
        mdn_primdist(br,er)=qnt_primdist(5,er,br);
        mdn_ksval(br,er)=qnt_ksval(3,er,br);
        mdn_kssig(br,er)=qnt_kssig(3,er,br);
        avg_slval(br,er)=mean(sub_slval);
        
    end
    MDN_INTRT(br,:)=mdn_intrt(br,[1 6 ERUNS]);
    QNT_INTRT(:,:,br)=qnt_intrt([2 4],[1 6 ERUNS],br);
    AVG_SL(br,:)=avg_sl(br,[1 6 ERUNS]);
    MDN_DLVR(br)=mean(mdn_dlvr(br,[1 6 ERUNS]));
    QNT_DLVR(:,:,br)=qnt_dlvr([2 4],[1 6 ERUNS],br);
    AVG_SLVAL(br,:)=avg_slval(br,[1 6 ERUNS]);
    MDN_RETURN(br,:)=mdn_return(br,[1 6 ERUNS]);
    QNT_RETURN(:,:,br)=qnt_return([2 4],[1 6 ERUNS],br);
    MDN_COST(br,:)=mdn_cost(br,[1 6 ERUNS]);
    QNT_COST(:,:,br)=qnt_cost([2 4],[1 6 ERUNS],br);
    MDN_DIST(br,:)=mdn_dist(br,[1 6 ERUNS]);
    QNT_DIST(:,:,br)=qnt_dist([2 4],[1 6 ERUNS],br);
    MDN_PRIMDIST(br,:)=mdn_primdist(br,[1 6 ERUNS]);
    QNT_PRIMDIST(:,:,br)=qnt_primdist([2 4],[1 6 ERUNS],br);
    MDN_KSVAL(br,:)=mdn_ksval(br,[1 6 ERUNS]);
    QNT_KSVAL(:,:,br)=qnt_ksval([2 4],[1 6 ERUNS],br);
    MDN_KSSIG(br,:)=mdn_kssig(br,[1 6 ERUNS]);
    QNT_KSSIG(:,:,br)=qnt_kssig([2 4],[1 6 ERUNS],br);
end


cd D:\Narcoscapes\Batch_ABM_Interdiction_Results\batch_pint\Results
%%%%%% Surface Pareto plots %%%%%%
rows=repmat(fplevels(:,2),1,ERUNS);
low_pint=linspace(0.1,10,6);
high_pint=linspace(10,20,6);
p_int=10.*[low_pint high_pint(2:6)];
cols=repmat(p_int,BRUNS,1);

% Interdiction rate
hh1=figure;
set(hh1,'Color','white');
surface(cols,rows,mdn_intrt(fplevels(:,1),:))
zlabel('Total Interdiction Rate (% of total flow)')
ylabel('Total Force Packages')
xlabel('Relative Prob. of Interception (%)')
title('Median Interdiction Rate')
view(3)
saveas(hh1,'IntRt','png')

% Seizure volume
hh2=figure;
set(hh2,'Color','white');
surface(cols,rows,avg_sl(fplevels(:,1),:)./1000)
zlabel('Avg. Total Volume Seized (MT)')
ylabel('Total Force Packages')
xlabel('Relative Prob. of Interception (%)')
title('Total Volume Seized')
view(3)
saveas(hh2,'SLvol','png')

% Seizure value
hh3=figure;
set(hh3,'Color','white');
surface(cols,rows,avg_slval(fplevels(:,1),:)./1000000)
zlabel('Avg. Total Value Seized ($M)')
ylabel('Total Force Packages')
xlabel('Relative Prob. of Interception (%)')
title('Total Value Seized')
view(3)
saveas(hh3,'SLval','png')

% Return time
hh4=figure;
set(hh4,'Color','white');
surface(cols,rows,mdn_return(fplevels(:,1),:))
zlabel('Median Trafficking Return Time')
ylabel('Total Force Packages')
xlabel('Relative Prob. of Interception (%)')
title('Trafficking Node Return Time')
view(3)
saveas(hh4,'ReturnTime','png')

% Displacement costs
hh5=figure;
set(hh5,'Color','white');
surface(cols,rows,mdn_cost(fplevels(:,1),:)./1000000)
zlabel('Median Network Displacement Costs ($M)')
ylabel('Total Force Packages')
xlabel('Relative Prob. of Interception (%)')
title('Displacement Costs')
view(3)
saveas(hh5,'RouteCost','png')

% Displacement distances
hh6=figure;
set(hh6,'Color','white');
surface(cols,rows,mdn_dist(fplevels(:,1),:))
zlabel('Median Network Displacement Dist (km)')
ylabel('Total Force Packages')
xlabel('Relative Prob. of Interception (%)')
title('Displacement Distance')
view(3)
saveas(hh6,'RouteDist','png')

% Displacement primary movements
hh7=figure;
set(hh7,'Color','white');
surface(cols,rows,mdn_primdist(fplevels(:,1),:))
zlabel('Median Primary Movement Displacement Dist (km)')
ylabel('Total Force Packages')
xlabel('Relative Prob. of Interception (%)')
title('Primary Movement Displacement Distance (km)')
view(3)
saveas(hh7,'PrimeDist','png')

% Portrait Divergence test stat
hh8=figure;
set(hh8,'Color','white');
surface(cols,rows,mdn_ksval(fplevels(:,1),:))
zlabel('Portrait Divergence Test Stat')
ylabel('Total Force Packages')
xlabel('Relative Prob. of Interception (%)')
title('Portrait Divergence')
view(3)
saveas(hh8,'NetDiv','png')

% Portrait Divergence p-value
hh9=figure;
set(hh9,'Color','white');
surface(cols,rows,mdn_kssig(fplevels(:,1),:))
zlabel('Portrait Divergence p-value')
ylabel('Total Force Packages')
xlabel('Relative Prob. of Interception (%)')
title('Portrait Divergence Significance')
view(3)
saveas(hh8,'NetDivPval','png')

%%%%%% Versus force packages plots %%%%%%%%
% Interdiction rate
h1=figure;
set(h1,'Color','white')
plot(fplevels(:,2),MDN_INTRT(fplevels(:,1)),'.-','MarkerSize',10)
ylabel('Total Interdiction Rate (% of total flow)')
xlabel('Total Force Packages')

% S&L volume
h2=figure;
set(h2,'Color','white')
plot(fplevels(:,2),AVG_SL(fplevels(:,1))./1000,'.-','MarkerSize',10)
ylabel('Total S&L Volume (MT)')
xlabel('Total Force Packages')

% S&L value
h3=figure;
set(h3,'Color','white')
plot(fplevels(:,2),AVG_SLVAL(fplevels(:,1))./1000000,'.-','MarkerSize',10)
ylabel('Total S&L Value (Million $)')
xlabel('Total Force Packages')

% Primary Displacement
h4=figure;
set(h4,'Color','white')
plot(fplevels(:,2),MDN_PRIMDIST(fplevels(:,1)),'.','MarkerSize',10)
ylabel('Primary Movement Displacement (km)')
xlabel('Total Force Packages')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Intra-Run Variability Analysis   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intcenters=intedges(1:10)+diff(intedges)./2;
histcols=repmat(p_int,10,1);
h_1=figure;
set(h_1,'Color','white')
Sintrt=surf(histcols,intcenters,hist_intrt(:,:,1),'FaceColor','interp',...
    'FaceAlpha','0.75','EdgeAlpha',0);
% Sintrt.EdgeColor='none';
xlabel('Relative Interception Probability (%)')
ylabel('Interdiction Rate')
zlabel('Model Replicant Count')
title('Interdiction Rate: Batch Run 1')
view(3)

h_2=figure;
set(h_2,'Color','white')
Sintrt=surf(histcols,intcenters,hist_intrt(:,:,9),'FaceColor','interp',...
    'FaceAlpha','0.75','EdgeAlpha',0);
% Sintrt.EdgeColor='none';
xlabel('Relative Interception Probability (%)')
ylabel('Interdiction Rate')
zlabel('Model Replicant Count')
title('Interdiction Rate: Batch Run 9 (baseline)')

h_3=figure;
set(h_3,'Color','white')
Sintrt=surf(histcols,intcenters,hist_intrt(:,:,24),'FaceColor','interp',...
    'FaceAlpha','0.75','EdgeAlpha',0);
% Sintrt.EdgeColor='none';
xlabel('Relative Interception Probability (%)')
ylabel('Interdiction Rate')
zlabel('Model Replicant Count')
title('Interdiction Rate: Batch Run 24')
