%%%%%%%%% Initialize NarcoLogic for execution from Python %%%%%%%%%%%%%%%%%
% %%% Notes for executing:

cd C:\Users\pcbmi\Box\NSF_D-ISN\Code\NarcoLogic

batchrun=9;
MRUNS=30;
ERUNS=11;
TSTART=1;
TMAX=180;   % 15 years at monthly time steps

rng default
load savedrngstate.mat

testflag=1;
erun=4;
mrun=1;

%%% Start initialization, set random number generator state for
%%% repeatability
rng(mrun)

%%% load experimental parameters file
[sl_max,sl_min,baserisk,riskmltplr,startstock,sl_learn,rt_learn,...
    losslim,prodgrow,targetseize,intcpctymodel,profitmodel,endstock,...
    growthmdl,timewght,locthink,expandmax,empSLflag,optSLflag,...
    suitflag,extnetflag,rtcap,basecap,p_sucintcpt]=load_expmntl_parms(ERUNS);

%%% Load landscape files
load coast_dist
load landsuit_file_default

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@ Agent Attributes @@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%%% Interdiction Agent %%%
delta_sl=sl_learn(erun);      % reinforcement learning rate for S&L vents (i.e., weight on new information)

%%% Network Agent %%%
ndto=2;         %initial number of DTOs
dtocutflag=zeros(ndto,1);
DTOBDGT=zeros(ndto,TMAX);
losstol=losslim(erun);        % tolerance threshold for loss due to S&L, triggers route fragmentation
stock_0=startstock(erun);     %initial cocaine stock at producer node
stock_max=endstock(erun);
startvalue=4500; %producer price, $385/kg: Zoe's numbers 4,500 in Panama
deltavalue=4.46;   %added value for distance traveled $8/kilo/km: Zoe's numbers $4.46
nodeloss=0;     % amount of cocaine that is normally lost (i.e., non-interdiction) at each node
ctrans_inland=371;  % transportation costs (kg/km) over-ground (3.5), includes
ctrans_coast=160;     % transportation costs (kg/km) via plane or boat (1.5)
ctrans_air=3486;
delta_rt=rt_learn(erun);       % reinforcement learning rate for network agent
%(i.e., weight on new information for successful routes)
% perceived risk model
alpharisk=2;
betarisk=0.5;
timewght_0=timewght(erun);
%slprob_0=alpharisk/(1+alpharisk+betarisk);     % baseline probability of seisure and loss event
slprob_0=1/(sum(timewght_0.^(0:12))+betarisk);
% cntrycpcty=[0.1 0.1 0.1 0.1 0.1 0.1 0.1];   %country-specific, per node trafficking capacity
bribepct=0.3;       % Annual proportion of gross profits from drug trafficking that go towards securing node territory
bribethresh=12;      % Maximum number of months a node can go without bribes to maintain control
rentcap=1-bribepct;     % proportion of value of shipments 'captured' by nodes
edgechange=expandmax(erun)*ones(ndto,1);       %Initial expansion of primary movement nodes, updated dynamically by route optimization

savedState=rng;
rng(thistate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Build trafficking network - NodeTable and EdgeTable   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load network_file_nodirect
EdgeTable.Capacity=rtcap(erun)*ones(height(EdgeTable),1);

nnodes=height(NodeTable);
mexnode=nnodes;
endnodeset=mexnode;
icoastdist=sub2ind(size(dcoast),NodeTable.Row,NodeTable.Col);
coastdist=dcoast(icoastdist);  %convert to km
NodeTable.CoastDist=coastdist;
NodeTable.CoastDist(1)=0;
NodeTable.CoastDist(nnodes)=0;
%%% Assign nodes to initials DTOs
for nn=2:nnodes-1
    westdir=NodeTable.Col(nn)-find(isnan(dcoast(NodeTable.Row(nn),1:...
        NodeTable.Col(nn)-1))==1,1,'last');
    eastdir=find(isnan(dcoast(NodeTable.Row(nn),NodeTable.Col(nn)+1:...
        size(LANDSUIT,2)))==1,1,'first');
    northdir=NodeTable.Row(nn)-find(isnan(dcoast(1:NodeTable.Row(nn)-1,...
        NodeTable.Col(nn)))==1,1,'last');
    southdir=find(isnan(dcoast(NodeTable.Row(nn)+1:size(LANDSUIT,1),...
        NodeTable.Col(nn)))==1,1,'first');
    [mindist,imindist]=min([westdir eastdir northdir southdir]);
    if westdir < 2.5*eastdir
        NodeTable.DTO(nn)=1;
    else
        NodeTable.DTO(nn)=2;
    end
end
load deptgrid
if extnetflag == 1
    [ext_NodeTable,ext_EdgeTable]=extend_network(nnodes,NodeTable,EdgeTable,Rdptgrid);
    NodeTable=ext_NodeTable;
    EdgeTable=ext_EdgeTable;
    nnodes=height(NodeTable);
    endnodeset=[mexnode 161:nnodes];
    EdgeTable.Capacity=basecap(erun)*ones(height(EdgeTable),1);
end

ADJ=zeros(nnodes);      % adjacency matrix for trafficking network
TRRTY=zeros(nnodes);    % control of nodes by each DTO
DIST=zeros(nnodes);     % geographic distance associated with edges
ADDVAL=zeros(nnodes);    % added value per edge in trafficking network
WGHT=ones(nnodes);     % dynamic weighting of edges
FLOW=zeros(nnodes,nnodes,TMAX);     % dynamic flows of cocaine between nodes
SLRISK=slprob_0*ones(nnodes);     % dynamic perceived risk of seisure and loss per edge by node agent
INTRISK=zeros(nnodes,TMAX); % dynamic perceived risk of interdiction at each node
CPCTY=zeros(nnodes);     % maximum flow possible between nodes
CTRANS=zeros(nnodes,nnodes,TMAX);   % transportation costs between nodes
RMTFAC=zeros(nnodes);   % landscape factor (remoteness) influencing S&L risk
COASTFAC=zeros(nnodes);   % landscape factor (distance to coast) influencing S&L risk
LATFAC=zeros(nnodes);   % decreased likelihood of S&L moving north to reflect greater DTO investment
BRDRFAC=zeros(nnodes);  % increased probability of S&L in department bordering an international border
SUITFAC=zeros(nnodes);

NEIHOOD=cell(nnodes,2);

STOCK=zeros(nnodes,TMAX);       %dynamic cocaine stock at each node
% DYNCAP=zeros(nnodes,TMAX);      %dynamic route capacity btw EPAC and CARB, based on JIATFS data
PRICE=zeros(nnodes,TMAX);       % $/kilo at each node
% RISKPREM=zeros(ndto,TMAX);      % risk premium after Caulkins et al. (1993)
RISKPREM=ones(nnodes,nnodes,TMAX);
INFLOW=zeros(nnodes,TMAX);       % dynamic stock of cocaine coming into at each node
OUTFLOW=zeros(nnodes,TMAX);       % dynamic stock of cocaine leaving from each node
TOTCPTL=zeros(nnodes,TMAX);     % total value of cocaine at each node
ICPTL=zeros(nnodes,TMAX);       % dynamic illicit capital accumulated at each node
LCPTL=zeros(nnodes,TMAX);       % dynamic legitimate capital accumulated at each node
BRIBE=zeros(nnodes,TMAX);       % annual bribe payments made at each node to maintain control
MARGIN=zeros(nnodes,TMAX);      % gross profit per node after purchasing, trafficking, and selling
RENTCAP=zeros(nnodes,TMAX);     % portion of MARGIN retained at node as profit
LEAK=zeros(nnodes,TMAX);        % dynamic amount of cocaine leaked at each node
activeroute=cell(nnodes,TMAX);  % track active routes
avgslrisk=cell(nnodes,TMAX);    % average S&L risk at each node given active routes
totslrisk=zeros(1,TMAX);        % network-wide average S&L risk
slcpcty=zeros(1,TMAX);

rng(savedState);
hitrngstate=rand(nnodes,1);

for k=1:nnodes
    % Create adjacency matrix (without graph toolbox)
    ADJ(k,EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==k,2))=1;
end

%%% Node Attributes
remotefac=[0; 1-NodeTable.PopSuit(2:nnodes)];
brdrfac=[0; NodeTable.DistBorderSuit(2:nnodes)];
suitfac=[0; NodeTable.LandSuit(2:nnodes)];
coastfac=[0; NodeTable.CoastDist(2:nnodes)./max(NodeTable.CoastDist(:))];
nwvec=sqrt(0.9.*NodeTable.Lat(2:nnodes).^2+0.1.*NodeTable.Lon(2:nnodes).^2);
latfac=[0; 1-nwvec./max(nwvec)];

% Create adjacency matrix (without graph toolbox)
iendnode=NodeTable.ID(NodeTable.DeptCode == 2);
ADJ(EdgeTable.EndNodes(EdgeTable.EndNodes(:,2)==iendnode,1),iendnode)=1;
iedge=find(ADJ == 1);
subneihood=zeros(size(LANDSUIT));
for j=1:nnodes
    %Create weight and capacity matrices
    WGHT(j,ADJ(j,:)==1)=EdgeTable.Weight(ADJ(j,:)==1);
    CPCTY(j,ADJ(j,:)==1)=EdgeTable.Capacity(ADJ(j,:)==1);
    
    % Create distance (in km) matrix
    latlon2=[NodeTable.Lat(ADJ(j,:)==1) NodeTable.Lon(ADJ(j,:)==1)];
    latlon1=repmat([NodeTable.Lat(j) NodeTable.Lon(j)],length(latlon2(:,1)),1);
    [d1km,d2km]=lldistkm(latlon1,latlon2);
    DIST(j,ADJ(j,:)==1)=d1km;
    if extnetflag == 1
        % add distance for extended network
        latlon1=repmat([NodeTable.Lat(1) NodeTable.Lon(1)],4,1);
        latlon2=[NodeTable.Lat([156; 161; 162; 163]) ...
            NodeTable.Lon([156; 161; 162; 163])];
        [d1km,d2km]=lldistkm(latlon1,latlon2);
        DIST(1,[156; 161; 162; 163])=d1km;
        
        latlon1=[NodeTable.Lat(157*ones(length(find(ADJ(157,:)==1)),1)) ...
            NodeTable.Lon(157*ones(length(find(ADJ(157,:)==1)),1))];
        latlon2=[NodeTable.Lat(find(ADJ(157,:)==1)') ...
            NodeTable.Lon(find(ADJ(157,:)==1)')];
        [d1km,d2km]=lldistkm(latlon1,latlon2);
        DIST(157,find(ADJ(157,:)==1)')=d1km;
        
        latlon1=[NodeTable.Lat(158*ones(length(find(ADJ(158,:)==1)),1)) ...
            NodeTable.Lon(158*ones(length(find(ADJ(158,:)==1)),1))];
        latlon2=[NodeTable.Lat(find(ADJ(158,:)==1)') ...
            NodeTable.Lon(find(ADJ(158,:)==1)')];
        [d1km,d2km]=lldistkm(latlon1,latlon2);
        DIST(158,find(ADJ(158,:)==1)')=d1km;
        
        latlon1=[NodeTable.Lat(159*ones(length(find(ADJ(159,:)==1)),1)) ...
            NodeTable.Lon(159*ones(length(find(ADJ(159,:)==1)),1))];
        latlon2=[NodeTable.Lat(find(ADJ(159,:)==1)') ...
            NodeTable.Lon(find(ADJ(159,:)==1)')];
        [d1km,d2km]=lldistkm(latlon1,latlon2);
        DIST(159,find(ADJ(159,:)==1)')=d1km;
        
        latlon1=[NodeTable.Lat(160*ones(length(find(ADJ(160,:)==1)),1)) ...
            NodeTable.Lon(160*ones(length(find(ADJ(160,:)==1)),1))];
        latlon2=[NodeTable.Lat(find(ADJ(160,:)==1)') ...
            NodeTable.Lon(find(ADJ(160,:)==1)')];
        [d1km,d2km]=lldistkm(latlon1,latlon2);
        DIST(160,find(ADJ(160,:)==1)')=d1km;
    end
    %Create added value matrix (USD) and price per node
    ADDVAL(j,ADJ(j,:)==1)=deltavalue.*DIST(j,ADJ(j,:)==1);
    if j == 1
        PRICE(j,TSTART)=startvalue;
    elseif ismember(j,endnodeset) == 1
        continue
    elseif ismember(j,157:160) == 1
        isender=EdgeTable.EndNodes(EdgeTable.EndNodes(:,2) == j,1);
        inextleg=EdgeTable.EndNodes(EdgeTable.EndNodes(:,1) == j,2);
        PRICE(j,TSTART)=PRICE(isender,TSTART)+ADDVAL(isender,j)+...
            PRICE(isender,TSTART)+mean(ADDVAL(j,inextleg));
        % even prices for long haul routes
        if j==160
            PRICE([157 160],TSTART)=min(PRICE([157 160],TSTART));
            PRICE([158 159],TSTART)=min(PRICE([158 159],TSTART));
        end
    else
        isender=EdgeTable.EndNodes(EdgeTable.EndNodes(:,2) == j,1);
        PRICE(j,TSTART)=mean(PRICE(isender,TSTART)+ADDVAL(isender,j));
    end
    for en=1:length(endnodeset)
        PRICE(endnodeset(en),TSTART)=max(PRICE(ADJ(:,endnodeset(en))==1,TSTART));
    end
    RMTFAC(j,ADJ(j,:)==1)=remotefac(ADJ(j,:)==1);
    COASTFAC(j,ADJ(j,:)==1)=coastfac(ADJ(j,:)==1);
    LATFAC(j,ADJ(j,:)==1)=latfac(ADJ(j,:)==1);
    BRDRFAC(j,ADJ(j,:)==1)=brdrfac(ADJ(j,:)==1);
    SUITFAC(j,ADJ(j,:)==1)=suitfac(ADJ(j,:)==1);
    
    % Transportation costs
    ireceiver=EdgeTable.EndNodes(EdgeTable.EndNodes(:,1) == j,2);
    idist_ground=(DIST(j,ireceiver)  >0 & DIST(j,ireceiver) <= 500);
    idist_air=(DIST(j,ireceiver) > 500);
    
    CTRANS(j,ireceiver(idist_ground),TSTART)=ctrans_inland.*...
        DIST(j,ireceiver(idist_ground))./DIST(1,mexnode);
    CTRANS(j,ireceiver(idist_air),TSTART)=ctrans_air.*...
        DIST(j,ireceiver(idist_air))./DIST(1,mexnode);

    if NodeTable.CoastDist(j) < 20 || ismember(j,157:159) == 1
        ireceiver=EdgeTable.EndNodes(EdgeTable.EndNodes(:,1) == j,2);
        idist_coast=(NodeTable.CoastDist(ireceiver) < 20);
        idist_inland=(NodeTable.CoastDist(ireceiver) >= 20);
        
        CTRANS(j,ireceiver(idist_coast),TSTART)=ctrans_coast.*...
            DIST(j,ireceiver(idist_coast))./DIST(1,mexnode);

        if ismember(j,157:159) == 1
            CTRANS(j,ireceiver(idist_coast),TSTART)=0;
            CTRANS(1,j,TSTART)=CTRANS(1,j,TSTART)+mean(ctrans_coast.*...
                DIST(j,ireceiver(idist_coast))./DIST(1,mexnode));
        end
    end
end

%%% Initialize Interdiction agent
% Create S&L probability layer
routepref=zeros(nnodes,nnodes,TMAX);   % weighting by network agent of successful routes
slevent=zeros(nnodes,nnodes,TMAX);  % occurrence of S&L event
intrdctobs=zeros(nnodes,nnodes,TMAX);
slnodes=cell(1,TMAX);
slsuccess=zeros(nnodes,nnodes,TMAX);    % volume of cocaine seized in S&L events
slvalue=zeros(nnodes,nnodes,TMAX);    % value of cocaine seized in S&L events
slcount_edges=zeros(1,TMAX);
slcount_vol=zeros(1,TMAX);
% intrdevent=zeros(nnodes,TMAX);
INTRDPROB=zeros(nnodes,TMAX);
SLPROB=zeros(nnodes,nnodes,TMAX);   % dynamic probability of S&L event per edge

if empSLflag(erun) == 1
    [empSLPROB,slctnodes] = build_SLemp(nnodes,TMAX,CAattr1,NodeTable,ADJ,ccdb);
    SLPROB=empSLPROB;
else
    facmat=LATFAC;
    facmat(:,:,2)=COASTFAC;
    facmat(:,:,3)=RMTFAC;
    facmat(:,:,4)=DIST./max(max(DIST));
    facmat(:,:,5)=BRDRFAC;
    facmat(:,:,6)=SUITFAC;
    SLPROB(:,:,TSTART)=mean(facmat(:,:,1:5),3);
    SLPROB(:,:,TSTART+1)=SLPROB(:,:,TSTART);
end
slmin=SLPROB(:,:,TSTART);
INTRDPROB(:,TSTART+1)=slprob_0*ones(nnodes,1); % dynamic probability of interdiction at nodes

%%% Initialize Node agents
STOCK(:,TSTART)=NodeTable.Stock(:);
% DYNCAP(NodeTable.DTO==1,TSTART)=rtcap(1,1)*basecap(erun);
% DYNCAP(NodeTable.DTO==2,TSTART)=rtcap(2,1)*basecap(erun);
TOTCPTL(:,TSTART)=NodeTable.Capital(:);
PRICE(:,TSTART+1)=PRICE(:,TSTART);
slcpcty_0=sl_min(erun);
slcpcty_max=sl_max(erun);
slcpcty(TSTART+1)=slcpcty_0;

% subjective risk perception with time distortion
twght=timewght_0*ones(nnodes,1);    % time weighting for dynamic, subjective perceived risk of interdiction event

%%% Set-up trafficking netowrk benefit-cost logic  %%%%%%%%%%%%
ltcoeff=locthink(erun)*ones(nnodes,1);
margval=zeros(nnodes,nnodes,TMAX);
for q=1:nnodes
    if isempty(find(ADJ(q,:)==1,1)) == 1
        continue
    end
    margval(q,q+1:nnodes,TSTART)=PRICE(q+1:nnodes,TSTART)-...
        PRICE(q,TSTART);
end
for nd=1:ndto
    idto=find(NodeTable.DTO == nd);
    margvalset=idto(~ismember(idto,endnodeset));
    routepref(1,idto,TSTART+1)=margval(1,idto)./max(margval(1,margvalset));
end
routepref(:,endnodeset,TSTART+1)=1;
totslrisk(TSTART+1)=1;

OWN=zeros(size(LANDSUIT));  % node agent land ownership
IOWN=cell(nnodes,TMAX);     % dynamic list of owned parcels
CTRANS(:,:,TSTART+1)=CTRANS(:,:,TSTART);
%%% Set-up figure for trafficking movie
MOV=zeros(nnodes,nnodes,TMAX);

% Output tables for flows(t) and interdiction prob(t-1)
t=TSTART +1;
if extnetflag(erun) == 1
    load init_flow_ext
else
    load init_flow
end
[rinit,cinit]=ind2sub([nnodes nnodes],find(FLOW(:,:,1) > 0));
for w=1:length(rinit)
    MOV(rinit(w),cinit(w),1)=FLOW(rinit(w),cinit(w),1);
end
[Tflow,Tintrd]=intrd_tables_batch(FLOW,slsuccess,SLPROB,NodeTable,EdgeTable,t,testflag,erun,mrun,batchrun);
clear t
%%% Write workspace .mat file to read by dynamic NarcoLogic function
wrksp_name='NarcoLogic_wrksp.mat';
save(wrksp_name)