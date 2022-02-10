%%%%%%%% Node and State Space Analysis %%%%%%%%%%%%%%

cd D:\Narcoscapes\Batch_ABM_Interdiction_Results\batch_pint
load \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM\GLOBALDIST
% load \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM\dept_nei_binning.mat
% load D:\Narcoscapes\Batch_ABM_Interdiction_Results\batch_pint\batch_pint_results.mat

[CAfull,CAfullattr]=shaperead('D:\CentralAmerica\GADM\CA_full_theater_0_simple_01.shp','UseGeoCoords',...
    true);

[CAadm1,CAattr1]=shaperead('D:\CentralAmerica\GADM\g2015_2014_1\CAadm1.shp',...
    'UseGeoCoords',true);
nnodes=length(GLOBDIST);
TSTART=1;
TMAX=180;
MRUNS=30;
ERUNS=11;
ndto=2;
endnodeset=[156 161 162 163];
nodelist=1:163;
% nodenames=strcat(num2str(nodelist'));
actnodeset=nodelist(~ismember(nodelist,endnodeset));
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

fnames=dir;
fnamescell=struct2cell(fnames);

%%%% Flow Metrics
prim_num=zeros(MRUNS,TMAX);
prim_cost=zeros(MRUNS,TMAX);
prim_dist=zeros(MRUNS,TMAX);
slsuc_edges=zeros(MRUNS,TMAX,BRUNS);
slsuc_vol=zeros(MRUNS,TMAX,BRUNS);
slsuc_val=zeros(MRUNS,TMAX,BRUNS);
slsuc_nodeind=cell(MRUNS,TMAX);
sl_nodeind=zeros(length(GLOBDIST),TMAX,MRUNS);
sl_nodevol=zeros(length(GLOBDIST),TMAX,MRUNS);
sl_nodeval=zeros(length(GLOBDIST),TMAX,MRUNS);
ntwk_nodes=zeros(MRUNS,TMAX,BRUNS);
ntwk_edges=zeros(MRUNS,TMAX,BRUNS);
ntwk_degree=zeros(MRUNS,TMAX,BRUNS);
ntwk_clstr=zeros(MRUNS,TMAX,BRUNS);
ntwk_spath=zeros(MRUNS,TMAX,BRUNS);
ntwk_lpath=zeros(MRUNS,TMAX,BRUNS);
ntwk_ddist=cell(MRUNS,TMAX,BRUNS);

%%% Node metrics
flowprob=zeros(length(GLOBDIST),length(GLOBDIST),MRUNS);
tflowprob=zeros(length(GLOBDIST),length(GLOBDIST),TMAX);
tnodestate=zeros(length(GLOBDIST),TMAX);
patheven=zeros(length(GLOBDIST),TMAX,MRUNS);
flowdiv=zeros(MRUNS,TMAX,BRUNS);
floweven=zeros(MRUNS,TMAX,BRUNS);

tsind=ceil((1:TMAX)./12);
for br=1:BRUNS
    for er=6
        for mr=MRUNS*ERUNS*(br-1)+MRUNS*(er-1)+1:MRUNS*ERUNS*(br-1)+MRUNS*(er-1)+MRUNS
            imr=mr-(MRUNS*ERUNS*(br-1)+MRUNS*(er-1));
            
            h=strcmp(sprintf('supplychain_results_optint_batch_%d_%d_%d.mat',...
                batchind(mr,1),batchind(mr,2),batchind(mr,3)),fnamescell(1,:));
            load(fnamescell{1,h==1})
            
            G=digraph(EdgeTable.EndNodes(:,1),EdgeTable.EndNodes(:,2),...
                GLOBDIST(sub2ind(size(GLOBDIST),EdgeTable.EndNodes(:,1),...
                EdgeTable.EndNodes(:,2))));
            
            %%% Flow metrics
            slsuc_edges(imr,:,br)=slcount_edges;
            slsuc_vol(imr,:,br)=slcount_vol;
            slsuc_val(imr,:,br)=sum(slval,1);
            
%             ntwk_nodes(imr,:,br)=sum(nactnodes,1);
            
            for t=2:TMAX
%                 prim_num(imr,t)=length(activeroute{1,t});
%                 prim_cost(imr,t)=sum(CTRANS(1,activeroute{1,t},t));
%                 prim_dist(imr,t)=sum(GLOBDIST(1,activeroute{1,t}));
%                 islnodes=slnodes{t};
%                 sl_nodeind(slnodes{t},t,imr)=1;
%                 sucind=unique(islnodes(sum(slsuccess(:,slnodes{t},t),1)>0));
%                 slsuc_nodeind(imr,t)=mat2cell(sucind,length(sucind),1);
%                 sl_nodevol(slnodes{t},t,imr)=sum(slsuccess(:,slnodes{t},t),1);
%                 sl_nodeval(slnodes{t},t,imr)=sum(repmat(PRICE(slnodes{t},t)',...
%                     length(GLOBDIST),1).*slsuccess(:,slnodes{t},t),1);
                if isempty(find(FLOW(:,:,t)>0,1)) == 0
                    [snd,rcv]=ind2sub(size(FLOW(:,:,t)),find(FLOW(:,:,t)>0));
                    ntwk_nodes(imr,t,br)=length(find(OUTFLOW(2:nnodes,t)>0));
                    ntwk_edges(imr,t,br)=length(snd);
                    actG=subgraph(G,unique([snd; rcv]));
                    actnodelist=1:length(unique(rcv))+1;
                    iendnodes=1+find(ismember(unique(rcv),endnodeset)==1);
                    [TR,D]=shortestpathtree(actG,1,'OutputForm','cell');
                    ntwk_spath(imr,t,br)=min(D(iendnodes));
                    ntwk_lpath(imr,t,br)=max(D(iendnodes));
                    inD=indegree(actG);
                    outD=outdegree(actG);
                    degrees=inD+outD;
                    %                 ntwk_degree(imr,t,br)=mean(degrees(~ismember(actnodelist,[1; iendnodes])));
                    ntwk_degree(imr,t,br)=mean(degrees);
                    ntwk_ddist(imr,t,br)=mat2cell(degrees,length(degrees),1);
                end
                
            
                tflowprob(:,:,t)=FLOW(:,:,t)./OUTFLOW(1,t);
                sub_tflowprob=tflowprob(:,:,t);
                flowdiv(imr,t,br)=-1*sum(sub_tflowprob(sub_tflowprob > 0).*...
                    log(sub_tflowprob(sub_tflowprob > 0)));
            end
            floweven(imr,:,br)=flowdiv(imr,:,br)./max(flowdiv(imr,:,br));
%             %%% Node metrics
%             flowprob(:,:,imr)=mean(tflowprob,3);
%             for t=2:TMAX
%                 for i=1:length(actnodeset)
%                     inode=actnodeset(i);
%                     isend=EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==inode,2);
%                     p_flow=FLOW(inode,isend,t)./OUTFLOW(inode,t);
%                     ikeep=p_flow > 0;
%                     tnodestate(inode,t)=-1*sum(p_flow(ikeep).*log(p_flow(ikeep)));
%                 end
%             end
%             t_hrzn=1;
%             for tt=2:TMAX-t_hrzn
%                 for i=1:length(actnodeset)
%                     inode=actnodeset(i);
%                     patheven(inode,tt,imr)=-1*sum((tnodestate(inode,tt)+...
%                         t_hrzn/tnodestate(inode,tt))*log(tnodestate(inode,tt)+...
%                         t_hrzn/tnodestate(inode,tt)));
%                 end
%             end
        end
    end
end

%%
cd D:\Narcoscapes\Batch_ABM_Interdiction_Results\batch_pint\Results
%%%%%%% Single scenario, single experimental setting %%%%%%%%%%%%%%%%%%%
%%% State-space plots %%%
% volspace=linspace(floor(min(min(slsuc_vol(:,2:TMAX)))./1000),ceil(max(max(slsuc_vol(:,2:TMAX)))./1000),100);
% valspace=linspace(floor(min(min(slsuc_val(:,2:TMAX)))./1000000),ceil(max(max(slsuc_val(:,2:TMAX)))./1000000),100);
% divspace=linspace(min(min(flowdiv(:,2:TMAX))),ceil(max(max(flowdiv(:,2:TMAX)))),100);
% intspace=linspace(floor(min(min(slsuc_edges(:,2:TMAX)))),ceil(max(max(slsuc_edges(:,2:TMAX)))),100);
% div_vol=zeros(length(divspace),length(volspace));
% div_val=zeros(length(divspace),length(valspace));
% div_int=zeros(length(divspace),length(intspace));
% vol_val=zeros(length(volspace),length(valspace));
mrun=18;
% brun=12;
Bruns=[9 1 6 19 24];
volmax=max(max(max(slsuc_vol)))./1000;
valmax=max(max(max(slsuc_val)))./1000000;
intmax=max(max(max(slsuc_edges)));
divmax=ceil(max(max(max(flowdiv))));
divchangemax=ceil(max(max(max(diff(flowdiv,1,2)))));
divchangemin=floor(min(min(min(diff(flowdiv,1,2)))));
dgrchangemax=ceil(max(max(max(diff(ntwk_degree,1,2)))));
dgrchangemin=floor(min(min(min(diff(ntwk_degree,1,2)))));
degreemax=ceil(max(max(max(ntwk_degree))));
nodemax=max(max(max(ntwk_nodes)));
nndchangemax=ceil(max(max(max(diff(ntwk_nodes,1,2)))));
nndchangemin=floor(min(min(min(diff(ntwk_nodes,1,2)))));
spathmax=max(max(max(ntwk_spath)));
lpathmax=max(max(max(ntwk_lpath)));
% Flow volume
cmap=[zeros(length(4:TMAX),1) 0.6*ones(length(4:TMAX),1) linspace(0,1,length(4:TMAX))'];

for bb=1:length(Bruns)
    brun=Bruns(bb);
    
    h1=figure;
    h1.Visible='off';
    set(h1,'Color','white')
    colormap(cmap)
    plot(slsuc_vol(mrun,3:TMAX-1,brun)./1000,flowdiv(mrun,4:TMAX,brun),'-k')
    hold on
    plot(slsuc_vol(mrun,3,brun)./1000,flowdiv(mrun,4,brun),'-ok',...
        'MarkerSize',8,'MarkerEdgeColor','k')
    for tt=4:TMAX-2
        rgb=[0 0.6 tt/TMAX];
        plot(slsuc_vol(mrun,tt,brun)./1000,flowdiv(mrun,tt+1,brun),'-ok',...
            'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    end
    rgb=[0 0.6 (TMAX-1)/TMAX];
    plot(slsuc_vol(mrun,TMAX-1,brun)./1000,flowdiv(mrun,TMAX,brun),'-xk',...
        'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    xlim([0 volmax])
    ylim([0 divmax])
    ylabel('Flow Diversity, t+1')
    xlabel('Interdiciton Volume (MT), t')
    c=colorbar('Ticks',[0,0.5,1],'TickLabels',{'3','90','179'});
    c.Label.String='Time Step';
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h1,sprintf('div_vol_%d_pac%d_carib%d.png',mrun,forcepkgs(brun,1),forcepkgs(brun,2)))
    
    hh1=figure;
    hh1.Visible='off';
    set(hh1,'Color','white')
    colormap(cmap)
    plot(slsuc_vol(mrun,3:TMAX-1,brun)./1000,ntwk_nodes(mrun,4:TMAX,brun),'-k')
    hold on
    plot(slsuc_vol(mrun,3,brun)./1000,ntwk_nodes(mrun,4,brun),'-ok',...
        'MarkerSize',8,'MarkerEdgeColor','k')
    for tt=4:TMAX-2
        rgb=[0 0.6 tt/TMAX];
        plot(slsuc_vol(mrun,tt,brun)./1000,ntwk_nodes(mrun,tt+1,brun),'-ok',...
            'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    end
    rgb=[0 0.6 (TMAX-1)/TMAX];
    plot(slsuc_vol(mrun,TMAX-1,brun)./1000,ntwk_nodes(mrun,TMAX,brun),'-xk',...
        'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    xlim([0 volmax])
    ylim([0 nodemax])
    ylabel('Number of Active Nodes, t+1')
    xlabel('Interdiciton Volume (MT), t')
    c=colorbar('Ticks',[0,0.5,1],'TickLabels',{'3','90','179'});
    c.Label.String='Time Step';
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(hh1,sprintf('nodes_vol_%d_pac%d_carib%d.png',mrun,forcepkgs(brun,1),forcepkgs(brun,2)))
    
    % Number of successful interdictions
    cmap=[zeros(length(4:TMAX),1) 0.6*ones(length(4:TMAX),1) linspace(0,1,length(4:TMAX))'];
    h2=figure;
    h2.Visible='off';
    set(h2,'Color','white')
    colormap(cmap)
    plot(slsuc_edges(mrun,3:TMAX-1,brun),flowdiv(mrun,4:TMAX,brun),'-k')
    hold on
    plot(slsuc_edges(mrun,3,brun),flowdiv(mrun,4,brun),'-ok',...
        'MarkerSize',8,'MarkerEdgeColor','k')
    for tt=4:TMAX-2
        rgb=[0 0.6 tt/TMAX];
        plot(slsuc_edges(mrun,tt,brun),flowdiv(mrun,tt+1,brun),'-ok',...
            'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    end
    rgb=[0 0.6 (TMAX-1)/TMAX];
    plot(slsuc_edges(mrun,TMAX-1,brun),flowdiv(mrun,TMAX,brun),'-xk',...
        'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    xlim([0 intmax])
    ylim([0 divmax])
    ylabel('Flow Diversity, t+1')
    xlabel('Number of Successful Interdictions, t')
    c=colorbar('Ticks',[0,0.5,1],'TickLabels',{'3','90','179'});
    c.Label.String='Time Step';
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h2,sprintf('div_int_%d_pac%d_carib%d.png',mrun,forcepkgs(brun,1),forcepkgs(brun,2)))
    
    %flow value
    cmap=[zeros(length(4:TMAX),1) 0.6*ones(length(4:TMAX),1) linspace(0,1,length(4:TMAX))'];
    h3=figure;
    h3.Visible='off';
    set(h3,'Color','white')
    colormap(cmap)
    plot(slsuc_val(mrun,3:TMAX-1,brun)./1000000,flowdiv(mrun,4:TMAX,brun),'-k')
    hold on
    plot(slsuc_val(mrun,3,brun)./1000000,flowdiv(mrun,4,brun),'-ok',...
        'MarkerSize',8,'MarkerEdgeColor','k')
    for tt=4:TMAX-2
        rgb=[0 0.6 tt/TMAX];
        plot(slsuc_val(mrun,tt,brun)./1000000,flowdiv(mrun,tt+1,brun),'-ok',...
            'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    end
    rgb=[0 0.6 (TMAX-1)/TMAX];
    plot(slsuc_val(mrun,TMAX-1,brun)./1000000,flowdiv(mrun,TMAX,brun),'-xk',...
        'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    xlim([0 valmax])
    ylim([0 divmax])
    ylabel('Flow Diversity, t+1')
    xlabel('Interdiciton Value ($Million), t')
    c=colorbar('Ticks',[0,0.5,1],'TickLabels',{'3','90','179'});
    c.Label.String='Time Step';
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h3,sprintf('div_val_%d_pac%d_carib%d.png',mrun,forcepkgs(brun,1),forcepkgs(brun,2)))
    
    %%% Network metrics
    h10_1=figure;
    h10_1.Visible='off';
    set(h10_1,'Color','white')
    colormap(cmap)
    plot(ntwk_degree(mrun,3:TMAX-1,brun),flowdiv(mrun,3:TMAX-1,brun),'-k')
    hold on
    plot(ntwk_degree(mrun,3,brun),flowdiv(mrun,3,brun),'-ok',...
        'MarkerSize',8,'MarkerEdgeColor','k')
    for tt=4:TMAX-2
        rgb=[0 0.6 tt/TMAX];
        plot(ntwk_degree(mrun,tt,brun),flowdiv(mrun,tt,brun),'-ok',...
            'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    end
    rgb=[0 0.6 (TMAX-1)/TMAX];
    plot(ntwk_degree(mrun,TMAX,brun),flowdiv(mrun,TMAX,brun),'-xk',...
        'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    xlim([0 degreemax])
    ylim([0 divmax])
    xlabel('Avg. Degree, t')
    ylabel('Flow Diversity, t')
    c=colorbar('Ticks',[0,0.5,1],'TickLabels',{'3','90','179'});
    c.Label.String='Time Step';
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h10_1,sprintf('degree_div_%d_pac%d_carib%d.png',mrun,forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h10_2=figure;
    h10_2.Visible='off';
    set(h10_2,'Color','white')
    colormap(cmap)
    plot(slsuc_edges(mrun,3:TMAX-1,brun),ntwk_degree(mrun,4:TMAX,brun),'-k')
    hold on
    plot(slsuc_edges(mrun,3,brun),ntwk_degree(mrun,4,brun),'-ok',...
        'MarkerSize',8,'MarkerEdgeColor','k')
    for tt=4:TMAX-2
        rgb=[0 0.6 tt/TMAX];
        plot(slsuc_edges(mrun,tt,brun),ntwk_degree(mrun,tt+1,brun),'-ok',...
            'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    end
    rgb=[0 0.6 (TMAX-1)/TMAX];
    plot(slsuc_edges(mrun,TMAX-1,brun),ntwk_degree(mrun,TMAX,brun),'-xk',...
        'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    xlim([0 intmax])
    ylim([0 degreemax])
    xlabel('Number of Successful Interdictions, t')
    ylabel('Average Degree, t+1')
    c=colorbar('Ticks',[0,0.5,1],'TickLabels',{'3','90','179'});
    c.Label.String='Time Step';
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h10_2,sprintf('int_degree_%d_pac%d_carib%d.png',mrun,forcepkgs(brun,1),forcepkgs(brun,2)))
    
    
    % Volume vs. value
    % cmap=[zeros(length(4:TMAX),1) 0.6*ones(length(4:TMAX),1) linspace(0,1,length(4:TMAX))'];
    h4=figure;
    h4.Visible='off';
    set(h4,'Color','white')
    colormap(cmap)
    plot(slsuc_vol(mrun,3:TMAX-1,brun)./1000,slsuc_val(mrun,3:TMAX-1,brun)./1000000,'-k')
    hold on
    plot(slsuc_vol(mrun,3,brun)./1000,slsuc_val(mrun,3,brun)./1000000,'-ok',...
        'MarkerSize',8,'MarkerEdgeColor','k')
    for tt=4:TMAX-1
        rgb=[0 0.6 tt/TMAX];
        plot(slsuc_vol(mrun,tt,brun)./1000,slsuc_val(mrun,tt,brun)./1000000,'-ok',...
            'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    end
    rgb=[0 0.6 (TMAX-1)/TMAX];
    plot(slsuc_vol(mrun,TMAX,brun)./1000,slsuc_val(mrun,TMAX,brun)./1000000,'-xk',...
        'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    xlim([0 volmax])
    ylim([0 valmax])
    xlabel('Interdiction Volume (MT), t')
    ylabel('Interdiciton Value ($Million), t')
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h1,sprintf('vol_val_%d_pac%d_carib%d.png',mrun,forcepkgs(brun,1),forcepkgs(brun,2)))
    
    hh4=figure;
    hh4.Visible='off';
    set(hh4,'Color','white')
    colormap(cmap)
        plot(flowdiv(mrun,3:TMAX-1,brun),ntwk_nodes(mrun,3:TMAX-1,brun),'-k')
    hold on
    plot(flowdiv(mrun,3,brun),ntwk_nodes(mrun,3,brun),'-ok',...
        'MarkerSize',8,'MarkerEdgeColor','k')
    for tt=4:TMAX-1
        rgb=[0 0.6 tt/TMAX];
        plot(flowdiv(mrun,tt,brun),ntwk_nodes(mrun,tt,brun),'-ok',...
            'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    end
    rgb=[0 0.6 (TMAX-1)/TMAX];
    plot(flowdiv(mrun,TMAX,brun),ntwk_nodes(mrun,TMAX,brun),'-xk',...
        'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',rgb)
    xlim([0 divmax])
    ylim([0 nodemax])
    xlabel('Flow Diversity, t')
    ylabel('Number of Active Nodes, t')
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(hh4,sprintf('div_nodes_%d_pac%d_carib%d.png',mrun,forcepkgs(brun,1),forcepkgs(brun,2)))

    %%%%%%%%%%   Density surfaces   %%%%%%%%%%%%%%%%%%%
    
%     for r=1:MRUNS
%         for tt=3:TMAX-1
%             ix=find(divspace <= flowdiv(r,tt+1,brun),1,'last');
%             iy=find(volspace <= slsuc_vol(r,tt,brun)./1000,1,'last');
%             div_vol(ix,iy)=div_vol(ix,iy)+1;
%     
%             iy=find(valspace <= slsuc_val(r,tt)./1000000,1,'last');
%             div_val(ix,iy)=div_val(ix,iy)+1;
%     
%             iy=find(intspace <= slsuc_edges(r,tt),1,'last');
%             div_int(ix,iy)=div_int(ix,iy)+1;
%     
%             ix=find(volspace <= slsuc_vol(r,tt)./1000,1,'last');
%             iy=find(valspace <= slsuc_val(r,tt)./1000000,1,'last');
%             vol_val(ix,iy)=vol_val(ix,iy)+1;
%         end
%     end
    
    h1_1=figure;
    h1_1.Visible='off';
    set(h1_1,'Color','white')
    scatter(reshape(slsuc_vol(:,2:TMAX-1,brun)./1000,1,MRUNS*(TMAX-2)),...
        reshape(flowdiv(:,3:TMAX,brun),1,MRUNS*(TMAX-2)))
%     surf(divspace,volspace,flipud(div_vol./(MRUNS*length(3:TMAX-1))),'FaceColor','interp',...
%         'FaceAlpha','0.85','EdgeAlpha',0);
    xlim([0 volmax])
    ylim([0 divmax])
    ylabel('Flow Diversity, t+1','FontSize',12)
    xlabel('Interdiciton Volume (MT), t','FontSize',12)
%     zlabel('Probability','FontSize',12)
%     xlim([min(divspace) max(divspace)])
%     zlim([0.00001 0.01])
    set(gca,'FontSize',12)
%     cb1=colorbar;
%     cb1.Label.String='Probability';
%     cb1.Label.FontSize=12;
%     view(2)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h1_1,sprintf('div_vol_t1_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h1_1_1=figure;
    h1_1_1.Visible='off';
    set(h1_1_1,'Color','white')
    scatter(reshape(slsuc_vol(:,2:TMAX,brun)./1000,1,MRUNS*(TMAX-1)),...
        reshape(flowdiv(:,2:TMAX,brun),1,MRUNS*(TMAX-1)))
    xlim([0 volmax])
    ylim([0 divmax])
    ylabel('Flow Diversity, t','FontSize',12)
    xlabel('Interdiciton Volume (MT), t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h1_1_1,sprintf('div_vol_t_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h1_1_2=figure;
    h1_1_2.Visible='off';
    set(h1_1_2,'Color','white')
    scatter(reshape(slsuc_vol(:,2:TMAX-1,brun)./1000,1,MRUNS*(TMAX-2)),...
        reshape(flowdiv(:,3:TMAX,brun)-flowdiv(:,2:TMAX-1,brun),1,MRUNS*(TMAX-2)))
    xlim([0 volmax])
    ylim([divchangemin divchangemax])
    ylabel('Change in Flow Diversity, t+1','FontSize',12)
    xlabel('Interdiciton Volume (MT), t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h1_1_2,sprintf('div_vol_change_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h2_1=figure;
    h2_1.Visible='off';
    set(h2_1,'Color','white')
    scatter(reshape(slsuc_val(:,2:TMAX,brun)./1000000,1,MRUNS*(TMAX-1)),...
        reshape(flowdiv(:,2:TMAX,brun),1,MRUNS*(TMAX-1)))
    xlim([0 valmax])
    ylim([0 divmax])
    ylabel('Flow Diversity, t','FontSize',12)
    xlabel('Interdiciton Value ($Million), t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h2_1,sprintf('div_val_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))

    h2_1_1=figure;
    h2_1_1.Visible='off';
    set(h2_1_1,'Color','white')
    scatter(reshape(slsuc_val(:,2:TMAX-1,brun)./1000000,1,MRUNS*(TMAX-2)),...
        reshape(flowdiv(:,3:TMAX,brun),1,MRUNS*(TMAX-2)))
    xlim([0 valmax])
    ylim([0 divmax])
    ylabel('Flow Diversity, t+1','FontSize',12)
    xlabel('Interdiciton Value ($Million), t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h2_1_1,sprintf('div_val_t_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h2_1_2=figure;
    h2_1_2.Visible='off';
    set(h2_1_2,'Color','white')
    scatter(reshape(slsuc_val(:,2:TMAX-1,brun)./1000000,1,MRUNS*(TMAX-2)),...
        reshape(flowdiv(:,3:TMAX,brun)-flowdiv(:,2:TMAX-1,brun),1,MRUNS*(TMAX-2)))
    xlim([0 valmax])
    ylim([divchangemin divchangemax])
    ylabel('Change in Flow Diversity, t+1','FontSize',12)
    xlabel('Interdiciton Value ($M), t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h2_1_2,sprintf('div_val_change_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h3_1=figure;
    h3_1.Visible='off';
    set(h3_1,'Color','white')
    scatter(reshape(slsuc_edges(:,2:TMAX,brun),1,MRUNS*(TMAX-1)),...
        reshape(flowdiv(:,2:TMAX,brun),1,MRUNS*(TMAX-1)))
    xlim([0 intmax])
    ylim([0 divmax])
    ylabel('Flow Diversity, t','FontSize',12)
    xlabel('Number of Successful Interdictions, t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h3_1,sprintf('div_int_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))

    h3_1_1=figure;
    h3_1_1.Visible='off';
    set(h3_1_1,'Color','white')
    scatter(reshape(slsuc_edges(:,2:TMAX-1,brun),1,MRUNS*(TMAX-2)),...
        reshape(flowdiv(:,3:TMAX,brun),1,MRUNS*(TMAX-2)))
    xlim([0 intmax])
    ylim([0 divmax])
    ylabel('Flow Diversity, t+1','FontSize',12)
    xlabel('Number of Successful Interdictions, t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h3_1_1,sprintf('div_int_t1_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h3_1_2=figure;
    h3_1_2.Visible='off';
    set(h3_1_2,'Color','white')
    scatter(reshape(slsuc_edges(:,2:TMAX-1,brun),1,MRUNS*(TMAX-2)),...
        reshape(flowdiv(:,3:TMAX,brun)-flowdiv(:,2:TMAX-1,brun),1,MRUNS*(TMAX-2)))
    xlim([0 intmax])
    ylim([divchangemin divchangemax])
    ylabel('Change in Flow Diversity, t+1','FontSize',12)
    xlabel('Number of Successful Interdiciton Events, t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h3_1_2,sprintf('div_int_change_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h20_1=figure;
    h20_1.Visible='off';
    set(h20_1,'Color','white')
    scatter(reshape(ntwk_degree(:,3:TMAX,brun),1,MRUNS*(TMAX-2)),...
        reshape(flowdiv(:,3:TMAX,brun),1,MRUNS*(TMAX-2)))
    xlim([0 degreemax])
    ylim([0 divmax])
    xlabel('Avg. Degree, t','FontSize',12)
    ylabel('Flow Diversity, t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h20_1,sprintf('degree_div_scat_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h20_2=figure;
    h20_2.Visible='off';
    set(h20_2,'Color','white')
    scatter(reshape(slsuc_vol(:,2:TMAX-1,brun)./1000,1,MRUNS*(TMAX-2)),...
        reshape(ntwk_degree(:,3:TMAX,brun),1,MRUNS*(TMAX-2)))
    xlim([0 volmax])
    ylim([0 degreemax])
    ylabel('Avg. Degree, t+1','FontSize',12)
    xlabel('Interdiction Volume (MT), t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h20_2,sprintf('degree_vol_scat_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h20_3=figure;
    h20_3.Visible='off';
    set(h20_3,'Color','white')
    scatter(reshape(slsuc_edges(:,2:TMAX-1,brun),1,MRUNS*(TMAX-2)),...
        reshape(ntwk_degree(:,3:TMAX,brun),1,MRUNS*(TMAX-2)))
    xlim([0 intmax])
    ylim([0 degreemax])
    ylabel('Avg. Degree, t+1','FontSize',12)
    xlabel('Number of Successful Interdictions, t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h20_3,sprintf('degree_int_scat_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))

    h4_1=figure;
    h4_1.Visible='off';
    set(h4_1,'Color','white')
    scatter(reshape(slsuc_vol(:,2:TMAX,brun)./1000,1,MRUNS*(TMAX-1)),...
        reshape(slsuc_val(:,2:TMAX,brun)./1000000,1,MRUNS*(TMAX-1)))
    xlim([0 volmax])
    ylim([0 valmax])
    xlabel('Interdiciton Volume (MT), t')
    ylabel('Interdiciton Value ($Million), t')
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h4_1,sprintf('vol_val_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
    
    h5_1=figure;
    h5_1.Visible='off';
    set(h5_1,'Color','white')
    scatter(reshape(flowdiv(:,2:TMAX,brun),1,MRUNS*(TMAX-1)),...
        reshape(ntwk_nodes(:,2:TMAX,brun),1,MRUNS*(TMAX-1)))
    xlabel('Flow Diversity, t','FontSize',12)
    ylabel('Number of Active Nodes, t','FontSize',12)
    % zlabel('Probability','FontSize',12)
    % xlim([min(divspace) max(divspace)])
    % zlim([0.00001 0.01])
    set(gca,'FontSize',12)
    % cb1=colorbar;
    % cb1.Label.String='Probability';
    % cb1.Label.FontSize=12;
    % view(2)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(brun,1),forcepkgs(brun,2)))
    saveas(h5_1,sprintf('div_nodes_pac%d_carib%d.png',forcepkgs(brun,1),forcepkgs(brun,2)))
end

%%%%% Pareto plots %%%%%%
% total volume
h5=figure;
set(h5,'Color','white')
vol_runsum=reshape(sum(slsuc_vol(:,2:TMAX,:),2)./1000,MRUNS,BRUNS);
vol_pareto=vol_runsum(:,fplevels(:,1));
Pvol=boxplot(vol_pareto,fplevels(:,2));
ylabel('Total Volume Seized (MT)','FontSize',12)
xlabel('Force Packages','FontSize',12)
set(gca,'FontSize',12)
saveas(h5,'pareto_totvol.png')

% total value
h6=figure;
set(h6,'Color','white')
val_runsum=reshape(sum(slsuc_val(:,2:TMAX,:),2)./1000000,MRUNS,BRUNS);
val_pareto=val_runsum(:,fplevels(:,1));
Pval=boxplot(val_pareto,fplevels(:,2));
ylabel('Total Value Seized ($Million)','FontSize',12)
xlabel('Force Packages','FontSize',12)
set(gca,'FontSize',12)
saveas(h6,'pareto_totval.png')

% total volume
h7=figure;
set(h7,'Color','white')
int_runsum=reshape(sum(slsuc_edges(:,2:TMAX,:),2),MRUNS,BRUNS);
int_pareto=int_runsum(:,fplevels(:,1));
Pint=boxplot(int_pareto,fplevels(:,2));
ylabel('Number of Successful Interdiction Events','FontSize',12)
xlabel('Force Packages','FontSize',12)
set(gca,'FontSize',12)
saveas(h7,'pareto_int.png')

% number of nodes
h8=figure;
set(h8,'Color','white')
node_runmdn=reshape(median(ntwk_nodes(:,2:TMAX,:),2),MRUNS,BRUNS);
node_pareto=node_runmdn(:,fplevels(:,1));
Pnode=boxplot(node_pareto,fplevels(:,2));
ylabel('Mdeian Number of Active Nodes','FontSize',12)
xlabel('Force Packages','FontSize',12)
set(gca,'FontSize',12)
saveas(h8,'pareto_nodes.png')

%%%%% Times series diagnostic %%%%%%%%
h30_1=figure;
set(h30_1,'Color','white')
subplot(6,1,1)
plot(3:TMAX,(slsuc_vol(mrun,3:TMAX,brun)./1000)./slsuc_edges(mrun,3:TMAX,brun),'-')
title('Volume (MT) per Event')
ylabel('S&L Volume (MT)')
subplot(6,1,2)
plot(3:TMAX,(slsuc_val(mrun,3:TMAX,brun)./1000000)./slsuc_edges(mrun,3:TMAX,brun),'-')
title('S&L Value ($M) per Event')
ylabel('Value ($M)')
subplot(6,1,3)
plot(3:TMAX,ntwk_edges(mrun,3:TMAX,brun),'-')
title('Number of Active Edges')
ylabel('#')
subplot(6,1,4)
plot(3:TMAX,ntwk_nodes(mrun,3:TMAX,brun),'-')
title('Number of Active Nodes')
ylabel('#')
subplot(6,1,5)
plot(3:TMAX,flowdiv(mrun,3:TMAX,brun),'-')
title('Flow Diversity')
ylabel('Flow Diversity')
subplot(6,1,6)
plot(3:TMAX,ntwk_degree(mrun,3:TMAX,brun),'-')
title('Average Node Degree')
ylabel('Degree')
xlabel('Time Step')
saveas(h30_1,'TimeSeries_line_9_18.png')

%%
%%%%%%%%%%%%%%% Distribution tests %%%%%%%%%%%%%%%%%%%%%%
fd=zeros(MRUNS*length(3:TMAX),BRUNS);
dgrs=zeros(MRUNS*length(3:TMAX),BRUNS);
nnds=zeros(MRUNS*length(3:TMAX),BRUNS);
vols=zeros(MRUNS*length(3:TMAX),BRUNS);
vals=zeros(MRUNS*length(3:TMAX),BRUNS);
sucints=zeros(MRUNS*length(3:TMAX),BRUNS);
for b=1:BRUNS
    fd(:,b)=reshape(flowdiv(:,3:TMAX,b)-flowdiv(:,2:TMAX-1,b),MRUNS*length(3:TMAX),1);
    dgrs(:,b)=reshape(ntwk_degree(:,3:TMAX,b)-ntwk_degree(:,2:TMAX-1,b),MRUNS*length(3:TMAX),1);
    nnds(:,b)=reshape(ntwk_nodes(:,3:TMAX,b)-ntwk_nodes(:,2:TMAX-1,b),MRUNS*length(3:TMAX),1);
    vols(:,b)=reshape(slsuc_vol(:,3:TMAX,b)./1000,MRUNS*length(3:TMAX),1);
    vals(:,b)=reshape(slsuc_val(:,3:TMAX,b)./1000000,MRUNS*length(3:TMAX),1);
    sucints(:,b)=reshape(slsuc_edges(:,3:TMAX,b),MRUNS*length(3:TMAX),1);
end
vol_edges=[0 quantile(reshape(vols,size(vols,1)*size(vols,2),1),[0.1 0.25 0.5 0.75 0.9]) volmax];
val_edges=[0 quantile(reshape(vals,size(vals,1)*size(vals,2),1),[0.1 0.25 0.5 0.75 0.9]) valmax];
int_edges=[0 quantile(reshape(sucints,size(sucints,1)*size(sucints,2),1),[0.1 0.25 0.5 0.75 0.9]) intmax];
fd_edges=[divchangemin quantile(reshape(fd,size(fd,1)*size(fd,2),1),[0.1 0.25 0.5 0.75 0.9]) divchangemax];
dgr_edges=[dgrchangemin quantile(reshape(dgrs,size(dgrs,1)*size(dgrs,2),1),[0.1 0.25 0.5 0.75 0.9]) dgrchangemax];
nnd_edges=[nndchangemin quantile(reshape(nnds,size(nnds,1)*size(nnds,2),1),[0.1 0.25 0.5 0.75 0.9]) nndchangemax];

%%% Volume, FD
[N_vol_div,~,~,binX_vol_div,binY_vol_div] = histcounts2(vols(:,1),fd(:,1),vol_edges,fd_edges);
bins_vol_div=1:size(N_vol_div,1)*size(N_vol_div,2);
expCounts_vol_div=reshape(N_vol_div,size(N_vol_div,1)*size(N_vol_div,2),1);

obsCounts_vol_div=zeros(size(expCounts_vol_div,1),BRUNS);
Cpatch_vol=zeros(size(expCounts_vol_div,1),BRUNS);
testmat_vol_div=zeros(size(N_vol_div,1),size(N_vol_div,2),BRUNS);
pmat_vol_div=zeros(size(N_vol_div,1),size(N_vol_div,2),BRUNS);
statmat_vol_div=zeros(size(N_vol_div,1),size(N_vol_div,2),BRUNS);

%%% Volume, Degree
[N_vol_dgr,~,~,binX_vol_dgr,binY_vol_dgr] = histcounts2(vols(:,1),dgrs(:,1),vol_edges,dgr_edges);
bins_vol_dgr=1:size(N_vol_dgr,1)*size(N_vol_dgr,2);
expCounts_vol_dgr=reshape(N_vol_dgr,size(N_vol_dgr,1)*size(N_vol_dgr,2),1);

obsCounts_vol_dgr=zeros(size(expCounts_vol_dgr,1),BRUNS);
testmat_vol_dgr=zeros(size(N_vol_dgr,1),size(N_vol_dgr,2),BRUNS);
pmat_vol_dgr=zeros(size(N_vol_dgr,1),size(N_vol_dgr,2),BRUNS);
statmat_vol_dgr=zeros(size(N_vol_dgr,1),size(N_vol_dgr,2),BRUNS);

%%% Value, FD
[N_val_div,Xedges,Yedges,binX_val_div,binY_val_div] = histcounts2(vals(:,1),fd(:,1),val_edges,fd_edges);
bins_val_div=1:size(N_val_div,1)*size(N_val_div,2);
expCounts_val_div=reshape(N_val_div,size(N_val_div,1)*size(N_val_div,2),1);

obsCounts_val_div=zeros(size(expCounts_val_div,1),BRUNS);
Cpatch_val=zeros(size(expCounts_val_div,1),BRUNS);
testmat_val_div=zeros(size(N_val_div,1),size(N_val_div,2),BRUNS);
pmat_val_div=zeros(size(N_val_div,1),size(N_val_div,2),BRUNS);
statmat_val_div=zeros(size(N_val_div,1),size(N_val_div,2),BRUNS);


for j=2:BRUNS
    [C_vol_div,~,~] = histcounts2(vols(:,j),fd(:,j),vol_edges,fd_edges);
    obsCounts_vol_div(:,j)=reshape(C_vol_div,size(C_vol_div,1)*size(C_vol_div,2),1);
    
    [C_vol_dgr,~,~] = histcounts2(vols(:,j),dgrs(:,j),vol_edges,dgr_edges);
    obsCounts_vol_dgr(:,j)=reshape(C_vol_dgr,size(C_vol_dgr,1)*size(C_vol_dgr,2),1);
    
    [C_val_div,~,~] = histcounts2(vals(:,j),fd(:,j),val_edges,fd_edges);
    obsCounts_val_div(:,j)=reshape(C_val_div,size(C_val_div,1)*size(C_val_div,2),1);
%     [h,p,st] = chi2gof(volbins,'Frequency',obsCounts(:,j),'Expected',expCounts);
    
    %%% Wilcoxon sign-ranked
    for k=1:length(bins_vol_div)
        [row,col]=ind2sub(size(N_vol_div),k);
        [p,h,stats]=signrank(fd(binX_vol_div == row & binY_vol_div == col,j),fd(binX_vol_div == row & binY_vol_div == col,1));
        testmat_vol_div(row,col,j)=h;
        pmat_vol_div(row,col,j)=p;
        statmat_vol_div(row,col,j)=stats.zval;
        
        [row,col]=ind2sub(size(N_vol_dgr),k);
        [p,h,stats]=signrank(dgrs(binX_vol_dgr == row & binY_vol_dgr == col,j),dgrs(binX_vol_dgr == row & binY_vol_dgr == col,1));
        testmat_vol_dgr(row,col,j)=h;
        pmat_vol_dgr(row,col,j)=p;
        statmat_vol_dgr(row,col,j)=stats.zval;
        
        [row,col]=ind2sub(size(N_val_div),k);
        [p,h,stats]=signrank(fd(binX_val_div == row & binY_val_div == col,j),fd(binX_val_div == row & binY_val_div == col,1));
        testmat_val_div(row,col,j)=h;
        pmat_val_div(row,col,j)=p;
        statmat_val_div(row,col,j)=stats.zval;
    end
    Cpatch_vol(:,j)=100*(obsCounts_vol_div(:,j)-expCounts_vol_div)./expCounts_vol_div;
    Cpatch_val(:,j)=100*(obsCounts_val_div(:,j)-expCounts_val_div)./expCounts_val_div;
end
% Create colormaps
[~,clredges_vol]=histcounts(Cpatch_vol,15);
clrstep_vol=clredges_vol(2)-clredges_vol(1);
zeropt=find(clredges_vol < 0,1,'last');
clrmap_vol=[ones(zeropt-1,1) linspace(0,0.8,zeropt-1)' linspace(0,0.8,zeropt-1)';
    1 1 1;
    linspace(0.8,0,length(zeropt+1:length(clredges_vol)))' ones(length(zeropt+1:length(clredges_vol)),1) linspace(0.8,0,length(zeropt+1:length(clredges_vol)))'];
% clrmap_vol=[ones(zeropt,1) linspace(0,0.8,zeropt)' linspace(0,0.8,zeropt)';...
%     linspace(0.8,0,length(zeropt+2:length(clredges_vol)))' ones(length(zeropt+2:length(clredges_vol)),1) linspace(0.8,0,length(zeropt+2:length(clredges_vol)))'];

[~,clredges_val]=histcounts(Cpatch_val,15);
clrstep_val=clredges_val(2)-clredges_val(1);
zeropt=find(clredges_val < 0,1,'last');
clrmap_val=[ones(zeropt-1,1) linspace(0,0.8,zeropt-1)' linspace(0,0.8,zeropt-1)';
    1 1 1;
    linspace(0.8,0,length(zeropt+1:length(clredges_val)))' ones(length(zeropt+1:length(clredges_val)),1) linspace(0.8,0,length(zeropt+1:length(clredges_val)))'];

for j=2:BRUNS
    %%%%%% Create patches %%%%%%
    %%% Volume, FD
    ycoords_vol_div=repmat(flipud(fd_edges'),1,length(vol_edges));
    xcoords_vol_div=repmat(vol_edges,length(fd_edges),1);
    % logxcoords=log10(xcoords+1);
    X_vol_div=zeros(4,size(N_vol_div,1)*size(N_vol_div,2));
    Y_vol_div=zeros(4,size(N_vol_div,1)*size(N_vol_div,2));
    for xx=1:size(N_vol_div,1)
        for yy=1:size(N_vol_div,2)
            ind=size(xcoords_vol_div,1)*(xx-1)+yy;
            polyind=size(N_vol_div,1)*(xx-1)+yy;
            X_vol_div(1,polyind)=xcoords_vol_div(ind);
            Y_vol_div(1,polyind)=ycoords_vol_div(ind);
            X_vol_div(2,polyind)=xcoords_vol_div(ind+1);
            Y_vol_div(2,polyind)=ycoords_vol_div(ind+1);
            X_vol_div(3,polyind)=xcoords_vol_div(ind+size(xcoords_vol_div,1)+1);
            Y_vol_div(3,polyind)=ycoords_vol_div(ind+size(xcoords_vol_div,1)+1);
            X_vol_div(4,polyind)=xcoords_vol_div(ind+size(xcoords_vol_div,1));
            Y_vol_div(4,polyind)=ycoords_vol_div(ind+size(xcoords_vol_div,1));
        end
    end
    cmap=colormap(winter(3));
    h100=figure;
    set(h100,'Color','white','Visible','off')
%     set(h100,'Color','white')
    scatter(reshape(slsuc_vol(:,2:TMAX-1,j)./1000,1,MRUNS*(TMAX-2)),...
        reshape(flowdiv(:,3:TMAX,j)-flowdiv(:,2:TMAX-1,j),1,MRUNS*(TMAX-2)),'.k')
    ylabel('Change in Flow Diversity, t+1','FontSize',12)
    xlabel('Interdiciton Volume (MT), t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(j,1),forcepkgs(j,2)))
    hold on
    % Cpatch=100*(obsCounts_vol_div(:,j)-expCounts_vol_div)./expCounts_vol_div;
%     colormap(h100,clrmap_vol)
    colormap(h100,cmap)
    patch(X_vol_div,Y_vol_div,flipud(Cpatch_vol(:,j)),'FaceAlpha',0.7,'EdgeColor','none')
    xlim([0 max(vol_edges)])
    ylim([min(fd_edges) max(fd_edges)])
    clb=colorbar;
    clb.Limits=[-50 300];
%     clb=colorbar('Ticks',[-50 -10 0 10 50 100],'TickLabels',{'< -50','-10','0','10','50','>100'});
%     clb.Ticks=[-50 -10 0 10 50 100];
%     clb.Ticks=linspace(-96,288,16);
    clb.Label.String='% Difference from Pac. 1, Carib. 1';
    plot(linspace(0,max(vol_edges),100),zeros(1,100),'-k')
    saveas(h100,sprintf('vol_div_diffcnt_pac%d_carib%d.png',forcepkgs(j,1),forcepkgs(j,2)))
    
    % %%% Volume, Degree
    % ycoords_vol_dgr=repmat(flipud(dgr_edges'),1,length(vol_edges));
    % xcoords_vol_dgr=repmat(vol_edges,length(dgr_edges),1);
    % % logxcoords=log10(xcoords+1);
    % X_vol_dgr=zeros(4,size(N_vol_dgr,1)*size(N_vol_dgr,2));
    % Y_vol_dgr=zeros(4,size(N_vol_dgr,1)*size(N_vol_dgr,2));
    % for xx=1:size(N_vol_dgr,1)
    %     for yy=1:size(N_vol_dgr,2)
    %         ind=size(xcoords_vol_dgr,1)*(xx-1)+yy;
    %         polyind=size(N_vol_dgr,1)*(xx-1)+yy;
    %         X_vol_dgr(1,polyind)=xcoords_vol_dgr(ind);
    %         Y_vol_dgr(1,polyind)=ycoords_vol_dgr(ind);
    %         X_vol_dgr(2,polyind)=xcoords_vol_dgr(ind+1);
    %         Y_vol_dgr(2,polyind)=ycoords_vol_dgr(ind+1);
    %         X_vol_dgr(3,polyind)=xcoords_vol_dgr(ind+size(xcoords_vol_dgr,1)+1);
    %         Y_vol_dgr(3,polyind)=ycoords_vol_dgr(ind+size(xcoords_vol_dgr,1)+1);
    %         X_vol_dgr(4,polyind)=xcoords_vol_dgr(ind+size(xcoords_vol_dgr,1));
    %         Y_vol_dgr(4,polyind)=ycoords_vol_dgr(ind+size(xcoords_vol_dgr,1));
    %     end
    % end
    % h101=figure;
    % set(h101,'Color','white')
    % scatter(reshape(slsuc_vol(:,2:TMAX-1,j)./1000,1,MRUNS*(TMAX-2)),...
    %     reshape(ntwk_degree(:,3:TMAX,j)-ntwk_degree(:,2:TMAX-1,j),1,MRUNS*(TMAX-2)),'.k')
    % ylabel('Change Average Node Degree, t+1','FontSize',12)
    % xlabel('Interdiciton Volume (MT), t','FontSize',12)
    % set(gca,'FontSize',12)
    % title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(j,1),forcepkgs(j,2)))
    % hold on
    % % Cpatch=100*(obsCounts_vol_dgr(:,j)-expCounts_vol_dgr)./expCounts_vol_dgr;
    % patch(X_vol_dgr,Y_vol_dgr,flipud(Cpatch),'FaceAlpha',0.7,'EdgeColor','none')
    % xlim([0 max(vol_edges)])
    % ylim([min(dgr_edges) max(dgr_edges)])
    % clb=colorbar;
    % clb.TickLabels=clredges_vol([2 4 6 8 10 12 14 16]);
    % clb.Label.String='% Difference from Pac. 1, Carib. 1';
    % saveas(h101,sprintf('vol_dgr_diffcnt_pac%d_carib%d.png',forcepkgs(j,1),forcepkgs(j,2)))
    
    %%% Value, FD
    ycoords_val_div=repmat(flipud(fd_edges'),1,length(val_edges));
    xcoords_val_div=repmat(val_edges,length(fd_edges),1);
    % logxcoords=log10(xcoords+1);
    X_val_div=zeros(4,size(N_val_div,1)*size(N_val_div,2));
    Y_val_div=zeros(4,size(N_val_div,1)*size(N_val_div,2));
    for xx=1:size(N_val_div,1)
        for yy=1:size(N_val_div,2)
            ind=size(xcoords_val_div,1)*(xx-1)+yy;
            polyind=size(N_val_div,1)*(xx-1)+yy;
            X_val_div(1,polyind)=xcoords_val_div(ind);
            Y_val_div(1,polyind)=ycoords_val_div(ind);
            X_val_div(2,polyind)=xcoords_val_div(ind+1);
            Y_val_div(2,polyind)=ycoords_val_div(ind+1);
            X_val_div(3,polyind)=xcoords_val_div(ind+size(xcoords_val_div,1)+1);
            Y_val_div(3,polyind)=ycoords_val_div(ind+size(xcoords_val_div,1)+1);
            X_val_div(4,polyind)=xcoords_val_div(ind+size(xcoords_val_div,1));
            Y_val_div(4,polyind)=ycoords_val_div(ind+size(xcoords_val_div,1));
        end
    end
    h102=figure;
    set(h102,'Color','white','Visible','off')
%     set(h102,'Color','white')
    scatter(reshape(slsuc_val(:,2:TMAX-1,j)./1000,1,MRUNS*(TMAX-2)),...
        reshape(flowdiv(:,3:TMAX,j)-flowdiv(:,2:TMAX-1,j),1,MRUNS*(TMAX-2)),'.k')
    ylabel('Change in Flow Diversity, t+1','FontSize',12)
    xlabel('Interdiciton value ($M), t','FontSize',12)
    set(gca,'FontSize',12)
    title(sprintf('Force Packages: Pac. %d, Carib. %d',forcepkgs(j,1),forcepkgs(j,2)))
    hold on
    % Cpatch=100*(obsCounts_val_div(:,j)-expCounts_val_div)./expCounts_val_div;
%     colormap(h102,clrmap_val)
    colormap(h102,cmap)
    patch(X_val_div,Y_val_div,flipud(Cpatch_val(:,j)),'FaceAlpha',0.7,'EdgeColor','none')
    xlim([0 max(val_edges)])
    ylim([min(fd_edges) max(fd_edges)])
    clb=colorbar;
    clb.Limits=[-100 200];
%     clb.Ticks=-50:50:300;
    clb.Label.String='% Difference from Pac. 1, Carib. 1';
    plot(linspace(0,max(val_edges),100),zeros(1,100),'-k')
    saveas(h102,sprintf('val_div_diffcnt_pac%d_carib%d.png',forcepkgs(j,1),forcepkgs(j,2)))
end
% % warning on
% % test for normality within bin
% normcheck=zeros(BRUNS,10);
% for bb=1:BRUNS
%     for i=1:10
%         normcheck(bb,i) = chi2gof(fd(bin_vol(:,bb) == i,bb),'cdf',...
%             {@normcdf,mean(fd(bin_vol(:,bb) == i,bb)),std(fd(bin_vol(:,bb) == i,bb))});
%     end
% end

% line(repmat(vol_edges,10,1),repmat(flipud(linspace(min(fd_edges),max(fd_edges),10)'),1,length(vol_edges)),'Color','k')
% hold on
% line(repmat(linspace(min(vol_edges),max(vol_edges),10)',1,length(fd_edges)),repmat(fd_edges,10,1),'Color','k')
% xlim([0 max(vol_edges)])
% ylim([min(fd_edges) max(fd_edges)])