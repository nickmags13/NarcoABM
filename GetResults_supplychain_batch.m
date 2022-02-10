%%%%%%%%%%%%% Get Results %%%%%%%%%%%%%%%%%%%%

% cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\model_results\SupplyChain_optint_071620
cd D:\Narcoscapes\Batch_ABM_Interdiction_Results\batch_pint
% cd D:\Narcoscapes\Batch_ABM_Interdiction_Results\batch_expndmax\SupplyChain_optint_batch_9
load \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM\GLOBALDIST.mat
load \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM\dept_nei_binning.mat

deptnei=readtable('\\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\My Documents\MATLAB\NarcoLogic\NarcoABM\dept_neighbors.txt');
fnames=dir;
fnamescell=struct2cell(fnames);
h=strncmp('supplychain_results_',fnamescell(1,:),20);
hind=find(h==1);

[CAadm0,CAattr0]=shaperead('D:\CentralAmerica\GADM\g2015_2014_0\CAadm0.shp',...
            'UseGeoCoords',true);
[CAadm1,CAattr1]=shaperead('D:\CentralAmerica\GADM\g2015_2014_1\CAadm1.shp',...
            'UseGeoCoords',true);

load(fnamescell{1,hind(1)})
NNODES=length(NodeTable.ID);
NodeTable.DeptCode([157 160])=4; %distinguish EPAC(3) from CARIB(4)
endnodeset=[156 161 162 163];
deptid=unique(cat(1,CAattr1.ADM1_CODE));
deptid=[1; deptid(~ismember(deptid,[1429 1432])); 3; 4]; %include producer node in deptid list
ineilist=1:length(deptid);
% deptid=[1; deptid(~ismember(deptid,[1429 1432]))]; %include producer node in deptid list
% ineilist=1:length(deptid);
% neinum=zeros(length(deptid));
% neicount=0;
% for n=1:length(deptid)
%     if deptid(n) == 1
%         ifirstnei=93670;
%     else
%         ifirstnei=deptnei.nbr_ADM1_CODE(deptnei.src_ADM1_CODE == deptid(n));
%     end
%     neinum(n,ismember(deptid,ifirstnei))=1;
%     unassigned=neinum(n,~ismember(ineilist,[1 n]));
%     while isempty(find(unassigned == 0,1)) == 0
%         neicount=neicount+1;
%         inextnei=deptnei.nbr_ADM1_CODE(ismember(deptnei.src_ADM1_CODE,deptid(neinum(n,:)==neicount)));
%         inextnei=unique(inextnei(~ismember(inextnei,[deptid(n); deptid(neinum(n,:)>0)])));
%         neinum(n,ismember(deptid,inextnei))=neicount+1;
%         unassigned=neinum(n,~ismember(ineilist,[1 n]));
%     end
%     neicount=0;
% end
% maxneidist=max(max(neinum));
% neinum=[neinum zeros(length(neinum(:,1)),2)];   %add dimension for EPAC and CARIB nodes
% neinum=[neinum; zeros(2,length(neinum(1,:)))];
% neinum(length(deptid)+1,:)=[(maxneidist+1)*ones(1,length(deptid)) 0 maxneidist+1]; %CARIB nodes
% neinum(:,length(deptid)+1)=[(maxneidist+1)*ones(1,length(deptid)) 0 maxneidist+1]';
% neinum(length(deptid)+2,:)=[(maxneidist+2)*ones(1,length(deptid)+1) 0]; %CARIB nodes
% neinum(:,length(deptid)+2)=[(maxneidist+2)*ones(1,length(deptid)+1) 0]';


TSTART=1;
TMAX=180;
MRUNS=30;
ERUNS=11;
ndto=2;
% batchind=[reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1) ...
%     reshape(repmat(1:MRUNS,1,ERUNS),MRUNS*ERUNS,1)];

BRUNS=24;   %batch runs of force package levels [batch, mrun, erun]
subbatchind=[reshape(repmat(1:MRUNS,1,ERUNS),MRUNS*ERUNS,1) ...
    reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1)];
batchind=[reshape(repmat(1:BRUNS,ERUNS*MRUNS,1),BRUNS*MRUNS*ERUNS,1) ...
    repmat(subbatchind,BRUNS,1)];
% sl_max=[75 100 125 150 175];
sl_max=125;
sl_min=ceil(sl_max/6);

% NACTNODES=cell(ndto,MRUNS*ERUNS);
% SLPEREVENT=cell(ndto,MRUNS*ERUNS);
NACTNODES=zeros(length(batchind(:,1)),TMAX,ndto);
SLPEREVENT=zeros(length(batchind(:,1)),TMAX,ndto);
VALPEREVENT=zeros(length(batchind(:,1)),TMAX,ndto);
dtoBDGT=cell(length(batchind(:,1)),1,ndto);
TMOV=cell(length(batchind(:,1)),1);
DEPTRPREM=cell(7,length(hind));
CNTRYRPREM=cell(7,length(hind));
DEPTTCOST=cell(7,length(hind));
CNTRYTCOST=cell(7,length(hind));
DEPTNODEPRICE=cell(7,length(hind));
CNTRYNODEPRICE=cell(7,length(hind));
SENDFLOW=zeros(NNODES,TMAX,length(batchind(:,1)));
SLVOL=zeros(NNODES,TMAX,length(batchind(:,1)));
INTRATE=zeros(NNODES,TMAX,length(batchind(:,1)));
DEPTFLOWS=cell(7,length(hind));
CNTRYFLOWS=cell(7,length(hind));
ACTROUTES=cell(length(hind),TMAX);
DISTROUTES=cell(length(hind),TMAX);
NEIBIN=zeros(max(max(neinum))+1,TMAX,length(hind));   %bin 1 is neighbor = 0, move within dept
DLVR=zeros(length(batchind(:,1)),TMAX);
ACTCHECK=zeros(1,length(batchind(:,1)));
% PRMYMV=zeros(1,length(batchind(:,1)));

edgecompare=cell(1,length(hind));

% for mr=1:length(hind)% MRUNS*EXPTRUNS
for mr=4802:length(hind) %batch9, run1, erun6   
    h=strcmp(sprintf('supplychain_results_optint_batch_%d_%d_%d.mat',...
        batchind(mr,1),batchind(mr,2),batchind(mr,3)),fnamescell(1,:));
%     h=strcmp(sprintf('supplychain_results_optint_batch_9_%d_%d.mat',...
%         batchind(mr,2),batchind(mr,3)),fnamescell(1,:));
%     h=strcmp(sprintf('supplychain_results_013018_%d_%d.mat',...
%         batchind(mr,1),batchind(mr,2)),fnamescell(1,:));
    filename=fnamescell{1,h};
    obj=matfile(filename);
    OUTFLOW=obj.OUTFLOW;
    slsuccess=obj.slsuccess;
    FLOW=obj.FLOW;
    nactnodes=obj.nactnodes;
    slperevent=obj.slperevent;
    slval=obj.slval;
    t_firstmov=obj.t_firstmov;
    
%     load(filename)
    
    % %%% Locate Department statistics
    cntrynames={'Guatemala','Panama','Panama','Costa Rica','Honduras','Panama','Nicaragua','Panama','Nicaragua','Guatemala','Honduras'};
    deptnames={'Pet','Dari','Ember','Puntarenas','Grac','Col','Atlantico Norte','Veraguas','Atlantico Sur','Izabal','Col'};
    cntryind=zeros(length(CAattr1),length(deptnames));
    deptind=zeros(length(CAattr1),length(deptnames));
    dind=zeros(1,length(deptnames));
    chind=cell(1,length(deptnames));
    slctnodes=cell(1,length(deptnames));
    cslctnodes=cell(1,length(deptnames));
    deptflows=zeros(length(deptnames),TMAX);
    deptflows_ts=zeros(length(deptnames),TMAX/12);
    cntryflows=zeros(length(deptnames),TMAX);
    cntryflows_ts=zeros(length(deptnames),TMAX/12);
    deptslvol=zeros(length(deptnames),TMAX);
    deptslvol_ts=zeros(length(deptnames),TMAX/12);
    deptintrt=zeros(length(deptnames),TMAX);
    deptintrt_ts=zeros(length(deptnames),TMAX/12);
%     deptriskprem=zeros(length(deptnames),TMAX);
%     deptriskprem_ts=zeros(length(deptnames),TMAX/12);
%     cntryriskprem=zeros(length(deptnames),TMAX);
%     cntryriskprem_ts=zeros(length(deptnames),TMAX/12);
%     deptctrans=zeros(length(deptnames),TMAX);
%     deptctrans_ts=zeros(length(deptnames),TMAX/12);
%     cntryctrans=zeros(length(deptnames),TMAX);
%     cntryctrans_ts=zeros(length(deptnames),TMAX/12);
    tsind=ceil((1:TMAX)./12);
    for jj=1:length(cntrynames)
        for ii=1:length(CAattr1)
            cntryind(ii,jj)=strncmp(cntrynames(jj),CAattr1(ii).ADM0_NAME,length(cntrynames{jj}));
            deptind(ii,jj)=strncmp(deptnames(jj),CAattr1(ii).ADM1_NAME,length(deptnames{jj}));
        end
        dind(jj)=find(cntryind(:,jj) == 1 & deptind(:,jj) == 1);
        chind(jj)=mat2cell(find(cntryind(:,jj) == 1),length(find(cntryind(:,jj) == 1)),1);
        nodeids=NodeTable.ID(NodeTable.DeptCode == CAattr1(dind(jj)).ADM1_CODE);
        cnodeids=NodeTable.ID(ismember(NodeTable.DeptCode,cat(1,CAattr1(chind{jj}).ADM1_CODE)));
        slctnodes(jj)=mat2cell(nodeids,length(nodeids),1);
        cslctnodes(jj)=mat2cell(cnodeids,length(cnodeids),1);
        deptflows(jj,:)=sum(OUTFLOW(slctnodes{jj},:),1);    %single run code
        cntryflows(jj,:)=sum(OUTFLOW(cslctnodes{jj},:),1);
%         RISKPREM(RISKPREM == 0)=NaN;
%         deptriskprem(jj,:)=mean(reshape(nanmean(RISKPREM(:,slctnodes{jj},1:TMAX),1),...
%             length(slctnodes{jj}),TMAX),1);
%         cntryriskprem(jj,:)=mean(reshape(nanmean(RISKPREM(:,cslctnodes{jj},1:TMAX),1),...
%             length(cslctnodes{jj}),TMAX),1);
%         CTRANS(CTRANS == 0)=NaN;
%         deptctrans(jj,:)=mean(reshape(nanmean(CTRANS(:,slctnodes{jj},1:TMAX),1),...
%             length(slctnodes{jj}),TMAX),1);
%         cntryctrans(jj,:)=mean(reshape(nanmean(CTRANS(:,cslctnodes{jj},1:TMAX),1),...
%             length(cslctnodes{jj}),TMAX),1);
        test=reshape(sum(slsuccess,1),NNODES,TMAX); %single run code
        deptslvol(jj,:)=sum(test(slctnodes{jj},:),1);  %single run code
        %     deptflows(jj,:)=sum(mediandlvr(slctnodes{jj},:),1);
        %     deptslvol(jj,:)=sum(medianslvol(slctnodes{jj},:),1);
        deptintrt(jj,:)=deptslvol(jj,:)./(deptslvol(jj,:)+deptflows(jj,:));
        for ts=1:TMAX/12
            deptflows_ts(jj,ts)=sum(deptflows(jj,tsind==ts));
            cntryflows_ts(jj,ts)=sum(cntryflows(jj,tsind==ts));
%             deptriskprem_ts(jj,ts)=mean(deptriskprem(jj,tsind==ts));
%             cntryriskprem_ts(jj,ts)=mean(cntryriskprem(jj,tsind==ts));
%             deptctrans_ts(jj,ts)=mean(deptctrans(jj,tsind==ts));
%             cntryctrans_ts(jj,ts)=mean(cntryctrans(jj,tsind==ts));
            deptslvol_ts(jj,ts)=sum(deptslvol(jj,tsind==ts));
            deptintrt_ts(jj,:)=deptslvol_ts(jj,:)./(deptslvol_ts(jj,:)+deptflows_ts(jj,:));
        end
        DEPTFLOWS(jj,mr)=mat2cell(deptflows_ts(jj,:),1,ts);
        CNTRYFLOWS(jj,mr)=mat2cell(cntryflows_ts(jj,:),1,ts);
%         DEPTRPREM(jj,mr)=mat2cell(deptriskprem_ts(jj,:),1,ts);
%         CNTRYRPREM(jj,mr)=mat2cell(cntryriskprem_ts(jj,:),1,ts);
%         DEPTTCOST(jj,mr)=mat2cell(deptctrans_ts(jj,:),1,ts);
%         CNTRYTCOST(jj,mr)=mat2cell(cntryctrans_ts(jj,:),1,ts);
    end
    
%     mrun_routes=cell(1,TMAX);
%     for t=1:TMAX
%         if isempty(cat(1,activeroute{:,t})) == 0
%             %            rtlength=length(cat(1,activeroute{:,t}));
% %             primary=cat(1,activeroute{1,t});
%             primary=int16(cat(1,activeroute{1,t}));
%             routes=cell([],1);
%             %             subroutes=mat2cell([1 primary(p)],1,2);
%             for p=1:length(primary)
%                 baseroute=mat2cell([1 primary(p)],1);
% %                 nextmov1=cat(1,activeroute{primary(p),t});
%                 nextmov1=int16(cat(1,activeroute{primary(p),t}));
%                 container=cell([],1);
%                 for k=1:length(nextmov1)
%                     container(size(container,1)+1,1)=mat2cell([baseroute{:} ...
%                         nextmov1(k)],ones(length(nextmov1(k)),1));
%                     if isempty(cat(1,activeroute{nextmov1(k),t})) == 0
% %                         nextmov2=cat(1,activeroute{nextmov1(k),t});
%                         nextmov2=int16(cat(1,activeroute{nextmov1(k),t}));
%                         if length(nextmov2)==1
%                             container(size(container,1),1)=mat2cell(...
%                                 [container{size(container,1),1} nextmov2],1);
%                         else
%                             edges=nextmov2;
%                             while length(edges) > 1
%                                 addind=size(container,1):size(container,1)+length(edges)-1;
%                                 container(size(container,1):size(container,1)+length(edges)-1,1)=...
%                                     mat2cell([repmat(container{size(container,1),1},length(edges),1) ...
%                                     edges],ones(length(edges),1));
%                                 
% %                                 actrt=cell([],1);
%                                 newrt=zeros([],1);
%                                 actedges=zeros([],1);
%                                 newnodes=zeros([],1);
%                                 nowrt=zeros(length(nextmov2),1);
%                                 clear inewedges iendedges
%                                 for g=1:length(nextmov2)
%                                     newmov=activeroute{nextmov2(g),t};
%                                     if isempty(find(newmov,1)) == 1 %no more nodes, no new routes
%                                         newrt(g,1)=0;
%                                         actedges(g,1)=0;
%                                         newnodes(g,1)=0;
%                                     elseif find(ismember(newmov,endnodeset),1) == 1 & ...
%                                             length(newmov) == 1 %new node to add to an existing route; send to 'nextmov3'
%                                         newrt(g,1)=0;
%                                         actedges(g,1)=1;
%                                         newnodes(g,1)=newmov;
%                                         nowrt(g,1)=1;
%                                     elseif find(ismember(newmov,endnodeset),1) == 1 & ...
%                                             length(newmov) > 1  
%                                         %new nodes to add to existing route
%                                         %and at least one new route; keep
%                                         %in this 'while loop'
%                                         newrt(g:g+length(newmov)-1,1)=ones(length(newmov),1);
%                                         actedges(g:g+length(newmov)-1,1)=...
%                                             ismember(newmov,endnodeset).*ones(length(newmov),1);
%                                         newnodes(g:g+length(newmov)-1,1)=newmov;
%                                         nowrt(g:g+length(newmov)-1,1)=ones(length(newmov),1);
%                                     end
% %                                     actrt(g,1)=activeroute(nextmov2(g),t)
%                                     inewedges=newrt == 1;
%                                     iendedges=actedges == 1;
%                                     inowrt=nowrt == 1;
%                                 end
% %                                 edges=cat(1,actrt{:});
%                                 edges=newnodes(inewedges);
% %                                 nextmov2=edges;
%                                 nextmov2=int16(edges);
%                             end
% %                             nextmov3=newnodes(inowrt);
%                             nextmov3=int16(newnodes(inowrt));
%                             container(addind(inowrt),1)=mat2cell(...
%                                     [cat(1,container{addind(inowrt),1}) nextmov3],ones(length(addind(inowrt)),1));
% %                             if length(nextmov3)==1
% % %                                 container(addind(iactedges),1)=mat2cell(...
% % %                                     [container{addind(iactedges),1} nextmov3],1);
% %                                 container(addind(inowrt),1)=mat2cell(...
% %                                     [container{addind(inowrt),1} nextmov3],1);
% %                             elseif length(nextmov3) > 1
% %                                 keyboard
% %                             end
%                         end
%                     end
%                 end
%                 routes(p)={container};
%             end
%             
%             
%                 
% %                     moves{size(moves,1)+1}=nextmov;
% %                 end
% %                 subroutes=mat2cell(repmat([1 primary(p)],length(moves),1),ones(length(moves),1));
% %                 
% %                 
% %                 moves=mat2cell([1; primary(p)],ones(2,1));
% %                 secondary=cat(1,activeroute{primary(p),t});
% % %                 moves{size(moves,1)+1}=secondary;
% %                 nextmov=secondary;
% %                 maxlength=0;
% %                 while isempty(nextmov) == 0
% %                     moves{size(moves,1)+1}=nextmov;
% %                     maxlength=max(maxlength,length(nextmov));
% %                     nextmov=cat(1,activeroute{nextmov,t});
% %                 end
% %                 subroutes=zeros(maxlength,size(moves,1));
% %                 subroutes(:,1:2)=repmat([1 primary(p)],maxlength,1);
% %                 for ii=1:maxlength
% %                     for jj=3:size(moves,1)
% %                     subroutes(:,ii)=
% %                 routes(size(routes,1)+1:size(routes,1)+size(moves,1),:)=
% %                     end
% %                 
% %                 end
% %             end
% %             
% %             for p=1:length(primary)
% %                 secondary=cat(1,activeroute{primary(p),t});
% %                 routes(size(routes,1)+1:size(routes,1)+length(secondary),:)=...
% %                     mat2cell(repmat([1 primary(p)],length(secondary),1),ones(length(secondary),1));
% % %                 subroutes=[routes{p} zeros(size(routes{p},1),1)];
% %                 subroutes=routes(size(routes,1)-length(secondary)+1:size(routes,1),1);
% %                 for s=1:length(secondary)
% %                     secflag=secondary(s);
% % %                     while isempty(secondary(s)) == 0
% % %                     secsteps=size(subroutes{s},2)+1;
% %                     while isempty(secflag) == 0
% %                     subroutes(s)=mat2cell([subroutes{s} secflag],1,length(subroutes{s})+1);
% % %                     secflag=cat(1,activeroute{secondary(s),t});
% %                     secflag=cat(1,activeroute{secflag,t});
% % %                     routes(p)=mat2cell([subroutes(s,:) secondary(s)],1,length(subroutes)+length(secondary(s)));
% % %                     secondary=cat(1,activeroute{secondary,t});
% % %                     secsteps=secsteps+1;
% %                     end
% %                 end
% %                 routes(size(routes,1)-length(secondary)+1:size(routes,1),1)=subroutes;
% %             end
%             mrun_routes(t)={routes};
%         end
%     end
%     ACTROUTES(mr)={mrun_routes};

    for t=1:TMAX
        flow_t=FLOW(:,:,t);
        if isempty(find(flow_t > 0,1)) == 1
            continue
        end
        actnetwork=flow_t > 0;
        ACTROUTES(mr,t)=mat2cell(actnetwork,NNODES,NNODES);
        [actsell,actbuy]=find(actnetwork == 1);
        deptsell=NodeTable.DeptCode(actsell);
        deptbuy=NodeTable.DeptCode(actbuy);
        rtdist=[];
        for j=1:length(actsell)
            rtdist=[rtdist; GLOBDIST(actsell(j),actbuy(j))];
            isell=find(deptid == NodeTable.DeptCode(actsell(j))); % reference neilist
            ibuy=find(deptid == NodeTable.DeptCode(actbuy(j)));
            NEIBIN(neinum(isell,ibuy)+1,t,mr)=NEIBIN(neinum(isell,ibuy)+1,t,mr)+1;
        end
        if isempty(find(rtdist,1)) == 1
            continue
        end
        DISTROUTES(mr,t)=mat2cell(rtdist(rtdist~=0),length(rtdist(rtdist~=0)),1);
    end
    
    %     ACTIVERTS{mr}=activeroute;
    SENDFLOW(:,:,mr)=OUTFLOW;
    SLVOL(:,:,mr)=reshape(sum(slsuccess,1),NNODES,TMAX);
    INTRATE(:,:,mr)=SLVOL(:,:,mr)./(SLVOL(:,:,mr)+SENDFLOW(:,:,mr));
    %     NACTNODES(:,mr)=mat2cell(nactnodes,ndto,TMAX);
    %     SLPEREVENT(:,mr)=mat2cell(slperevent,ndto,TMAX);
    for idt=1:ndto
        NACTNODES(mr,:,idt)=nactnodes(idt,:);
        SLPEREVENT(mr,:,idt)=slperevent(idt,:);
        VALPEREVENT(mr,:,idt)=slval(idt,:);
    end
%     dtoBDGT(mr,1,1)=mat2cell(DTOBDGT(1,:),1,TMAX);
%     dtoBDGT(mr,1,2)=mat2cell(DTOBDGT(2,:),1,TMAX);
    TMOV(mr,1)=mat2cell(t_firstmov',1,length(t_firstmov));
    
    edgecompare{mr}=EdgeTable.EndNodes;
    
    DLVR(mr,:)=SENDFLOW(1,:,mr)-sum(SLPEREVENT(mr,:,:),3);
    ACTCHECK(mr)=sum(sum(SENDFLOW(:,:,mr),2) > 0);
%     PRMYMV(mr)=length(unique(cat(1,activeroute{1,:})));

    clear nactnodes slperevent
end


% slevents=cell2mat(SLPEREVENT);
% actrtes=cell2mat(NACTNODES);
slevents=SLPEREVENT;
actrtes=NACTNODES;
totactrtes=sum(NACTNODES,3);
totslevents=sum(SLPEREVENT,3);
totvalevents=sum(VALPEREVENT./max(SLPEREVENT,1),3);
% dtobdgt_1=cell2mat(dtoBDGT(:,:,1));
% dtobdgt_2=cell2mat(dtoBDGT(:,:,2));
tmov=cell2mat(TMOV);

meanrtes=zeros(ndto,TMAX,ERUNS);
medianrtes=zeros(ndto,TMAX,ERUNS);
minrtes=zeros(ndto,TMAX,ERUNS);
maxrtes=zeros(ndto,TMAX,ERUNS);
cvrtes=zeros(ndto,TMAX,ERUNS);
meansl=zeros(ndto,TMAX,ERUNS);
mediansl=zeros(ndto,TMAX,ERUNS);
minsl=zeros(ndto,TMAX,ERUNS);
maxsl=zeros(ndto,TMAX,ERUNS);
cvsl=zeros(ndto,TMAX,ERUNS);
maxdto=zeros(length(hind),2);
meanfirstmov=zeros(ERUNS,NNODES);
maxflows=zeros(NNODES,ERUNS);
maxyear=zeros(NNODES,ERUNS);

medianslvol=zeros(NNODES,TMAX,ERUNS);
meanslvol=zeros(NNODES,TMAX,ERUNS);
mediandlvr=zeros(NNODES,TMAX,ERUNS);
meandlvr=zeros(NNODES,TMAX,ERUNS);
medianintrt=zeros(NNODES,TMAX,ERUNS);
meanintrt=zeros(NNODES,TMAX,ERUNS);
slvolQnt=cell(7,ERUNS);    % quantiles [0.025 0.05 0.25 0.5 0.75 0.95 0.975] for NNODES x TMAX
dlvrQnt=cell(7,ERUNS);
intrtQnt=cell(7,ERUNS);

%%% BRUN storage
brun_stats=cell(BRUNS,23);

for brun=1:BRUNS
% for brun=9
    subSLVOL=SLVOL(:,:,batchind(:,1)==brun);
    subSENDFLOW=SENDFLOW(:,:,batchind(:,1)==brun);
    subINTRATE=INTRATE(:,:,batchind(:,1)==brun);
    sub_actrtes=actrtes(batchind(:,1)==brun,:,:);
    sub_slevents=slevents(batchind(:,1)==brun,:,:);
    sub_tmov=tmov(batchind(:,1)==brun,:);
    
    for erun=1:ERUNS
%     for erun=6
        medianslvol(:,:,erun)=median(subSLVOL(:,:,subbatchind(:,2)==erun),3); %volume seized
        meanslvol(:,:,erun)=mean(subSLVOL(:,:,subbatchind(:,2)==erun),3); %volume seized
        mediandlvr(:,:,erun)=median(subSENDFLOW(:,:,subbatchind(:,2)==erun),3); %volume delivered
        meandlvr(:,:,erun)=mean(subSENDFLOW(:,:,subbatchind(:,2)==erun),3); %volume delivered
        medianintrt(:,:,erun)=median(subINTRATE(:,:,subbatchind(:,2)==erun),3); %interdiction rate
        meanintrt(:,:,erun)=mean(subINTRATE(:,:,subbatchind(:,2)==erun),3); %interdiction rate
        
        slvol_qnt=quantile(subSLVOL(:,:,subbatchind(:,2)==erun),[0.025 0.05 0.25 0.5 0.75 0.95 0.975],3);
        dlvr_qnt=quantile(subSENDFLOW(:,:,subbatchind(:,2)==erun),[0.025 0.05 0.25 0.5 0.75 0.95 0.975],3);
        intrt_qnt=quantile(subINTRATE(:,:,subbatchind(:,2)==erun),[0.025 0.05 0.25 0.5 0.75 0.95 0.975],3);
        
        for k=1:7
            slvolQnt(k,erun)=mat2cell(slvol_qnt(:,:,k),NNODES,TMAX);
            dlvrQnt(k,erun)=mat2cell(slvol_qnt(:,:,k),NNODES,TMAX);
            intrtQnt(k,erun)=mat2cell(slvol_qnt(:,:,k),NNODES,TMAX);
        end
        
        for idt=1:ndto
            %         istart=find(slevents(mr,:,idt)>0,1,'first');
            %         meanrtes(mr,:,idt)=mean(actrtes(mr,istart:TMAX,idt));
            %         medianrtes(idt,mr)=median(actrtes(mr,istart:TMAX,idt));
            %         minrtes(idt,mr)=min(actrtes(mr,istart:TMAX,idt));
            %         maxrtes(idt,mr)=max(actrtes(mr,istart:TMAX,idt));
            %         cvrtes(idt,mr)=var(actrtes(mr,istart:TMAX,idt))/mean(actrtes(mr,istart:TMAX,idt));
            %         meansl(idt,mr)=mean(slevents(mr,istart:TMAX,idt));
            %         mediansl(idt,mr)=median(slevents(mr,istart:TMAX,idt));
            %         minsl(idt,mr)=min(slevents(mr,istart:TMAX,idt));
            %         maxsl(idt,mr)=max(slevents(mr,istart:TMAX,idt));
            %         cvsl(idt,mr)=var(slevents(mr,istart:TMAX,idt))/mean(slevents(mr,istart:TMAX,idt));
            
            meanrtes(idt,:,erun)=mean(sub_actrtes(subbatchind(:,2)==erun,:,idt),1);
            medianrtes(idt,:,erun)=median(sub_actrtes(subbatchind(:,2)==erun,:,idt));
            minrtes(idt,:,erun)=min(sub_actrtes(subbatchind(:,2)==erun,:,idt));
            maxrtes(idt,:,erun)=max(sub_actrtes(subbatchind(:,2)==erun,:,idt));
            cvrtes(idt,:,erun)=var(sub_actrtes(subbatchind(:,2)==erun,:,idt))./...
                mean(sub_actrtes(subbatchind(:,2)==erun,:,idt));
            meansl(idt,:,erun)=mean(sub_slevents(subbatchind(:,2)==erun,:,idt));
            mediansl(idt,:,erun)=median(sub_slevents(subbatchind(:,2)==erun,:,idt));
            minsl(idt,:,erun)=min(sub_slevents(subbatchind(:,2)==erun,:,idt));
            maxsl(idt,:,erun)=max(sub_slevents(subbatchind(:,2)==erun,:,idt));
            cvsl(idt,:,erun)=var(sub_slevents(subbatchind(:,2)==erun,:,idt))./...
                mean(sub_slevents(subbatchind(:,2)==erun,:,idt));
        end
        
        for j=1:length(tmov(1,:))
            isubmov=sub_tmov(subbatchind(:,2)==erun,j) ~= 0;
            meanfirstmov(erun,j)=mean(sub_tmov(isubmov,j),1);
        end
    end

    summroutes=cell(length(sl_max),5);
    summsl=cell(length(sl_max),5);
    for g=1:length(sl_max)
        ind=batchind(:,1) == g;
        
        [rtemu,rtesigma,rtemuci,~]=normfit(sum(medianrtes(:,ind),1));
        rtequant=quantile(sum(medianrtes(:,ind),1),[0.025 0.25 0.5 0.75 0.975]);
        summroutes{g,1}=rtemu;
        summroutes{g,2}=rtesigma;
        summroutes(g,3)=mat2cell(rtemuci,2,1);
        summroutes{g,4}=rtequant(3);
        summroutes(g,5)=mat2cell([rtequant(2) rtequant(4)]',2,1);
        
        [slmu,slsigma,slmuci,~]=normfit(sum(mediansl(:,ind),1));
        slquant=quantile(sum(mediansl(:,ind),1),[0.025 0.25 0.5 0.75 0.975]);
        summsl{g,1}=slmu;
        summsl{g,2}=slsigma;
        summsl(g,3)=mat2cell(slmuci,2,1);
        summsl{g,4}=slquant(3);
        summsl(g,5)=mat2cell([slquant(2) slquant(4)]',2,1);
    end
    brun_stats(brun,1)=mat2cell(meanrtes,size(meanrtes,1),size(meanrtes,2),size(meanrtes,3));
    brun_stats(brun,2)=mat2cell(medianrtes,size(medianrtes,1),size(medianrtes,2),size(medianrtes,3));
    brun_stats(brun,3)=mat2cell(minrtes,size(minrtes,1),size(minrtes,2),size(minrtes,3));
    brun_stats(brun,4)=mat2cell(maxrtes,size(maxrtes,1),size(maxrtes,2),size(maxrtes,3));
    brun_stats(brun,5)=mat2cell(cvrtes,size(cvrtes,1),size(cvrtes,2),size(cvrtes,3));
    brun_stats(brun,6)=mat2cell(meansl,size(meansl,1),size(meansl,2),size(meansl,3));
    brun_stats(brun,7)=mat2cell(mediansl,size(mediansl,1),size(mediansl,2),size(mediansl,3));
    brun_stats(brun,8)=mat2cell(minsl,size(minsl,1),size(minsl,2),size(minsl,3));
    brun_stats(brun,9)=mat2cell(maxsl,size(maxsl,1),size(maxsl,2),size(maxsl,3));
    brun_stats(brun,10)=mat2cell(cvsl,size(cvsl,1),size(cvsl,2),size(cvsl,3));
    brun_stats(brun,11)=mat2cell(maxdto,size(maxdto,1),size(maxdto,2));
    brun_stats(brun,12)=mat2cell(meanfirstmov,size(meanfirstmov,1),size(meanfirstmov,2));
    brun_stats(brun,13)=mat2cell(maxflows,size(maxflows,1),size(maxflows,2));
    brun_stats(brun,14)=mat2cell(maxyear,size(maxyear,1),size(maxyear,2));
    
    brun_stats(brun,15)=mat2cell(medianslvol,size(medianslvol,1),size(medianslvol,2),size(medianslvol,3));
    brun_stats(brun,16)=mat2cell(meanslvol,size(meanslvol,1),size(meanslvol,2),size(meanslvol,3));
    brun_stats(brun,17)=mat2cell(mediandlvr,size(mediandlvr,1),size(mediandlvr,2),size(mediandlvr,3));
    brun_stats(brun,18)=mat2cell(meandlvr,size(meandlvr,1),size(meandlvr,2),size(meandlvr,3));
    brun_stats(brun,19)=mat2cell(medianintrt,size(medianintrt,1),size(medianintrt,2),size(medianintrt,3));
    brun_stats(brun,20)=mat2cell(meanintrt,size(meanintrt,1),size(meanintrt,2),size(meanintrt,3));
    brun_stats(brun,21)={slvolQnt};
    brun_stats(brun,22)={dlvrQnt};
    brun_stats(brun,23)={intrtQnt};
end

% % cd D:\Narcoscapes\Batch_ABM_Interdiction_Results\batch_rtcap
cd D:\Narcoscapes\Batch_ABM_Interdiction_Results\batch_pint
save('batch_pint_datastore.mat','ACTCHECK','brun_stats','CNTRYFLOWS','DEPTFLOWS','DLVR',...
    'INTRATE','SENDFLOW','SLPEREVENT','SLVOL','TMOV','totslevents',...
    'totvalevents','VALPEREVENT','ACTROUTES','DISTROUTES','NEIBIN','-v7.3')
% save('batch_pint_actroutes.mat','ACTROUTES','-v7.3')



