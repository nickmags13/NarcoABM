%%%%%%%%%%%%%%%% Upload or build suitability layer %%%%%%%%%%%%%%%%%%%%%%%%
function [LANDSUIT]=load_suitability(buildsuit)

if buildsuit == 1       %load RAT suitability model
    
else
    % Tree cover
    [tcov,~]=geotiffread('D:\CentralAmerica\Model_inputs\clipped\treecov_clp.tif');
    treecov=tcov;
    %         tcovnodataval=255;
    treecov=double(treecov);
    treecov(cagrid_cntry==cntrynodataval)=NaN;
    itreecov=find(treecov > 0); %identify high forest cover areas
    treecovpct=quantile(treecov(itreecov),[0.025 0.25 0.50 0.75 0.975]);
    itreepick=find(treecov >= treecovpct(2) & treecov < treecovpct(3));
    avgtcov=zeros(length(dptcodes),1);    %average tree cover per county, weighting for generating trade nodes
    for cc=1:length(dptcodes)
        avgtcov(cc)=mean(treecov(ca_adm0 == dptcodes(cc)));
    end
    
    % Distance to coast and country borders
    [dcoast,~]=geotiffread('D:\CentralAmerica\Model_inputs\clipped\dcoast_clp.tif');
    dcoastnodataval=-9999;
    dcoast(cagrid_cntry==cntrynodataval)=NaN;
    dcoast(dcoast == dcoastnodataval)=NaN;
    dcoast_suit=1-dcoast./max(max(dcoast));
    
    % Population density as a proxy for remoteness
    [popden,~]=geotiffread('D:\CentralAmerica\Model_inputs\clipped\popden_clp.tif');
    popnodataval=-1;
    popden=double(popden);
    popden(cagrid_cntry==cntrynodataval)=NaN;
    popden(popden == popnodataval)=NaN;
    popq=quantile(reshape(popden,size(popden,1)*size(popden,2),1),...
        [0.25 0.5 0.75 0.95]);
    pop_suit=zeros(size(popden));
    pop_suit(popden > popq(3))=0;
    pop_suit(popden <= popq(3))=1-popden(popden <= popq(3))./popq(3);
    
    % Topography
    [slope,~]=geotiffread('D:\CentralAmerica\Model_inputs\ca_slope_250m.tif');
    slope(cagrid_cntry==cntrynodataval)=NaN;
    slopeclass=[8 16 30 31; 0 25 50 100]';   % GAEZ (see Magliocca et al., 2013, PLOS ONE)
    slp_suit=zeros(size(slope));
    slp_suit(slope < slopeclass(1,1))=1-slopeclass(1,2)/100;
    slp_suit(slope >= slopeclass(1,1) & slope < slopeclass(2,1))=1-slopeclass(2,2)/100;
    slp_suit(slope >= slopeclass(2,1) & slope < slopeclass(3,1))=1-slopeclass(3,2)/100;
    slp_suit(slope >= slopeclass(4,1))=1-slopeclass(4,2)/100;
    
    % Market Access
    [mktacc,~]=geotiffread('D:\CentralAmerica\Model_inputs\clipped\mktacc_clp.tif');
    manodataval=-9999;
    mktacc=double(mktacc);
    mktacc(cagrid_cntry==cntrynodataval)=NaN;
    mktacc(mktacc == manodataval)=NaN;
    mktacc_suit=mktacc;
    % submasuit=mktacc./median(mktacc(~isnan(mktacc)));
    
    %         % Maize Yield
    %         [mazyld,Rmazyld]=geotiffread('C:\Users\nrmagliocca\Box Sync\Data Drive\CentralAmerica\Model_inputs\clipped\mazyld_clp.tif');
    %         maznodataval=-9999;
    %         mazyld=double(mazyld);
    %         mazyld(cagrid_cntry==cntrynodataval)=NaN;
    %         mazyld(mazyld == maznodataval)=NaN;
    %
    %         % Oil Palm Yield
    %         [plmyld,Rplmyld]=geotiffread('C:\Users\nrmagliocca\Box Sync\Data Drive\CentralAmerica\Model_inputs\clipped\plmyld_clp.tif');
    %         plmnodataval=-9999;
    %         plmyld=double(plmyld);
    %         plmyld(cagrid_cntry==cntrynodataval)=NaN;
    %         plmyld(plmyld == plmnodataval)=NaN;
    %
    %         % Cattle density
    %         [ctlden,Rctlden]=geotiffread('C:\Users\nrmagliocca\Box Sync\Data Drive\CentralAmerica\Model_inputs\clipped\ctlden_clp.tif');
    %         plmnodataval=-9999;
    %         ctlden=double(ctlden);
    %         ctlden(cagrid_cntry==cntrynodataval)=NaN;
    %         ctlden(ctlden == plmnodataval)=NaN;
    
    % Initial land use
    % Land cover classes:
    % 1. built-up
    % 2. cropland — row crop agriculture
    % 3. shrubs
    % 4. trees
    % 5. pastureland/grassland — grazing land or natural grassland
    % 6. bare
    % 7. plantation — citrus, vineyard, coffee, etc.
    % 8. water
    % 9. plantation tree — eucalyptus, pine, etc.
    [luint,~]=geotiffread('D:\CentralAmerica\Model_inputs\clipped\luint_clp.tif');
    lunodataval=0;
    luint(cagrid_cntry==cntrynodataval)=NaN;
    luint(luint == lunodataval)=NaN;
    lu_suit=zeros(size(luint));
    lu_suit(luint == 1 | luint == 6 | luint == 7 | luint == 8 | luint ==9)=0;
    lu_suit(luint == 2)=0.5;
    lu_suit(luint == 3 | luint == 4 | luint == 5)=1;
    
    % Protected Areas
    [protarea,~]=geotiffread('D:\CentralAmerica\Model_inputs\clipped\protarea_clp.tif');
    protnodataval=255;
    protarea=double(protarea);
    protarea(cagrid_cntry==cntrynodataval)=NaN;
    protarea(protarea == protnodataval)=NaN;
    protsuit=1-isnan(protarea);
    
    %%% Land-base investment value
    invst_suit=zeros(size(luint));
    invst_suit(luint == 1 | luint == 6 | luint == 7 | luint == 8 | luint ==9)=0;
    invst_suit(luint == 2)=0.5;
    invst_suit(luint == 3 | luint == 5)=0.75;
    invst_suit(luint == 4)=1;
    
    %%% Weight each landscape attribute
    tcwght=1;       % tree cover (narco = 1)
    brdwght=1;      % distance to country border (narco = 1)
    dcstwght=0;     % distance to coast (always 0)
    mktwght=1;      % market access (always 1)
    popwght=1;      % population density - proxy for remoteness (narco = 1)
    slpwght=1;      % slope-constrained land suitability (always 1)
    luwght=0;       % suitability based on initial land use (narco = 0)
    invstwght=0;    % investment potential of initial land use (always 0)
    protwght=1;     % protected area status (always 1)
    
    wghts=[tcwght brdwght dcstwght mktwght popwght slpwght luwght invstwght protwght]./...
        sum([tcwght brdwght dcstwght mktwght popwght slpwght luwght invstwght protwght]);
    
    % %         %%% Null Model
    %         LANDSUIT=wghts(1).*treecov./100+wghts(2).*dbrdr_suit+wghts(3).*dcoast_suit+...
    %             wghts(4).*mktacc_suit+wghts(5).*pop_suit+wghts(6).*slp_suit+wghts(6).*...
    %             lu_suit+wghts(7).*invst_suit+wghts(8)*(1-protsuit);  % land suitability based on biophysical and narco variable predictors
    
    %%% Full model
    LANDSUIT=wghts(1).*treecov./100+wghts(2).*dbrdr_suit+wghts(3).*dcoast_suit+...
        wghts(4).*(1-mktacc_suit)+wghts(5).*pop_suit+wghts(6).*slp_suit+wghts(6).*...
        lu_suit+wghts(7).*invst_suit+wghts(8)*protsuit;  % land suitability based on biophysical and narco variable predictors
    
end