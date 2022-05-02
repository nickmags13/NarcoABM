%%%%%%%%% Function for executing NarcoLogic dynamics from Python %%%%%%%%%%

function [Tflow]=NarcoLogic_dynamics_python(t)

load NarcoLogic_wrksp.mat

[intrdct_events,intrdct_nodes]=optimize_interdiction_batch(t,ADJ,testflag,erun,mrun,batchrun);
slevent(:,:,t)=intrdct_events;
slnodes(t)=mat2cell(intrdct_nodes,size(intrdct_nodes,1),size(intrdct_nodes,2));

MOV(:,1,t)=NodeTable.Stock(:);

%%% Iterate through trafficking nodes
for n=1:nnodes %skip end nodes
    if isempty(find(ADJ(n,:)==1,1)) == 1
        continue
    end
    if ismember(n,endnodeset) == 1
        continue
    end
    
    %%%%%  Route cocaine shipmments %%%%%
    STOCK(n,t)=STOCK(n,t-1)+STOCK(n,t);
    rtdto=NodeTable.DTO(ADJ(n,:)==1);
    if isempty(find(rtdto == 0,1)) == 0
        rtdto(rtdto == 0)=NodeTable.DTO(n);
    end
    CPCTY(n,ADJ(n,:)==1)=basecap(erun)*rtcap(rtdto,floor(t/12)+1);
    TOTCPTL(n,t)=TOTCPTL(n,t-1)+TOTCPTL(n,t);
    %       ICPTL(n,t)=ICPTL(n,t-1)+ICPTL(n,t);
    if STOCK(n,t) > 0
        if n > 1
            LEAK(n,t)=nodeloss*STOCK(n,t); %drugs 'leaked' at each node
            STOCK(n,t)=STOCK(n,t)-LEAK(n,t);
        end
        %             inei=zeros([],1);
        if n == 1
            inei=find(ADJ(n,:) == 1 & routepref(n,:,t) > 0);
            %                         for nd=1:length(unique(NodeTable.DTO(2:nnodes-1)))
            for nd=1:length(unique(NodeTable.DTO(2:nnodes)))
                if isempty(find(NodeTable.DTO(inei) == nd,1)) == 1
                    idtombr=(NodeTable.DTO == nd);
                    subinei=find(ADJ(n,idtombr) == 1 & routepref(n,idtombr,t) > 0);
                    if isempty(find(subinei,1)) == 1
                        subinei=find(ADJ(n,idtombr) == 1 & routepref(n,idtombr,t) == ...
                            max(routepref(n,idtombr,t)));
                    end
                    inei=[inei subinei];
                end
            end
        else
            inei=find(ADJ(n,:) == 1 & routepref(n,:,t) > 0);
            inei=inei(ismember(inei,[find(NodeTable.DTO == ...
                NodeTable.DTO(n)); endnodeset']));
            if isempty(find(inei,1)) == 1
                inei=find(ADJ(n,:) == 1 & routepref(n,:,t) == ...
                    max(routepref(n,:,t)));
                inei=inei(ismember(inei,[find(NodeTable.DTO == ...
                    NodeTable.DTO(n)); endnodeset']));
            end
        end
        %%% Procedure for selecting routes based on expected profit %%%
        c_trans=CTRANS(n,inei,t);
        p_sl=SLRISK(n,inei);
        y_node=PRICE(inei,t)-PRICE(n,t);
        q_node=min(STOCK(n,t)./length(inei),CPCTY(n,inei));
        lccf=ltcoeff(n);
        totstock=STOCK(n,t);
        totcpcty=CPCTY(n,inei);
        tslrisk=totslrisk(t);
        rtpref=routepref(n,inei,t);
        dtonei=NodeTable.DTO(inei);
        profmdl=profitmodel(erun);
        cutflag=dtocutflag(unique(dtonei(dtonei~=0)));
        
        [neipick,neivalue,valuex]=calc_neival(c_trans,p_sl,y_node,...
            q_node,lccf,rtpref,tslrisk,dtonei,profmdl,cutflag,totcpcty,totstock,edgechange);
        
        % With top-down route optimization
        inei=inei(neipick);
        
        %  % Just bottom-up route optimization
        %  expmax=expandmax(erun);
        %  neipick=neipick(1:min(expmax,length(neipick)));
        %  inei=inei(neipick);
        
        % weight according to salience value fuction
        if isempty(find(valuex <= 0,1)) == 0
            WGHT(n,inei)=(1-SLRISK(n,inei))./sum(1-SLRISK(n,inei));
        else
            WGHT(n,inei)=(max(valuex(neipick),0)./sum(max(valuex(neipick),0)))';
        end
        activeroute(n,t)=mat2cell(inei',length(inei),1);
        
        neiset=unique(NodeTable.DTO(inei));
        
        FLOW(n,inei,t)=min(WGHT(n,inei)./sum(WGHT(n,inei)).*...
            STOCK(n,t),CPCTY(n,inei));
        OUTFLOW(n,t)=sum(FLOW(n,inei,t));
        STOCK(n,t)=STOCK(n,t)-OUTFLOW(n,t);
        nodecosts=sum(FLOW(n,inei,t).*CTRANS(n,inei,t));
        % Check for S%L event
        if isempty(find(ismember(find(slevent(n,:,t)),inei),1)) == 0
            isl=find(slevent(n,inei,t)==1);
            intrdctobs(n,inei(isl),t)=1;
            intcpt=min(p_sucintcpt(erun)*NodeTable.pintcpt(inei(isl)),1);
            
            %%% interception probability
            p_int=rand(length(intcpt),1);
            for p=1:length(intcpt)
                if p_int(p) <= intcpt(p)
                    slsuccess(n,inei(isl(p)),t)=FLOW(n,inei(isl(p)),t);
                    slvalue(n,inei(isl(p)),t)=FLOW(n,inei(isl(p)),t).*PRICE(inei(isl(p)),t)';
                    FLOW(n,inei(isl(p)),t)=0;     % remove from trafficking route due to S&L event
                else
                    slsuccess(n,inei(isl(p)),t)=0;
                    slvalue(n,inei(isl(p)),t)=0;
                end
            end
            STOCK(inei,t)=STOCK(inei,t)+FLOW(n,inei,t)';
            noderevenue=sum(FLOW(n,inei,t).*PRICE(inei,t)');
            TOTCPTL(inei,t)=TOTCPTL(inei,t)-(FLOW(n,inei,t)'.*PRICE(inei,t));
            ICPTL(n,t)=rentcap*sum(FLOW(n,inei).*ADDVAL(n,inei));
            MARGIN(n,t)=noderevenue-nodecosts+min(TOTCPTL(n,t),0);
            if n > 1
                BRIBE(n,t)=max(bribepct*MARGIN(n,t),0);
                if MARGIN(n,t) > 0
                    RENTCAP(n,t)=MARGIN(n,t)-BRIBE(n,t);
                else
                    RENTCAP(n,t)=MARGIN(n,t);
                end
                TOTCPTL(n,t)=max(TOTCPTL(n,t),0)+RENTCAP(n,t);  % losses on top of debt capture by MARGIN
            else
                RENTCAP(n,t)=MARGIN(n,t);
                TOTCPTL(n,t)=TOTCPTL(n,t)+RENTCAP(n,t);
            end
        else
            STOCK(inei,t)=STOCK(inei,t)+FLOW(n,inei,t)';
            nodecosts=sum(FLOW(n,inei,t).*CTRANS(n,inei,t));
            noderevenue=sum(FLOW(n,inei,t).*PRICE(inei,t)');
            TOTCPTL(inei,t)=TOTCPTL(inei,t)-(FLOW(n,inei,t)'.*PRICE(inei,t));
            ICPTL(n,t)=rentcap*sum(FLOW(n,inei).*ADDVAL(n,inei));
            MARGIN(n,t)=noderevenue-nodecosts+min(TOTCPTL(n,t),0);
            if n > 1
                BRIBE(n,t)=max(bribepct*MARGIN(n,t),0);
                if MARGIN(n,t) > 0
                    RENTCAP(n,t)=MARGIN(n,t)-BRIBE(n,t);
                else
                    RENTCAP(n,t)=MARGIN(n,t);
                end
                TOTCPTL(n,t)=max(TOTCPTL(n,t),0)+RENTCAP(n,t);
            else
                RENTCAP(n,t)=MARGIN(n,t);
                TOTCPTL(n,t)=TOTCPTL(n,t)+RENTCAP(n,t);
            end
        end
        %%%% Update perceived risk in response to S&L and Interdiction events
        timeweight=twght(n);
        % identify neighbors in network (without network toolbox)
        %       bcknei=EdgeTable.EndNodes(EdgeTable.EndNodes(:,2) == n,1)';
        fwdnei=inei;
        t_eff=0:12;
        if t == TSTART+1
            % Risk perception only updated when successful
            % interdiction takes place
            sloccur=[zeros(12,length(fwdnei)); (slsuccess(n,fwdnei,TSTART+1)>0)];
        elseif t > TSTART+1 && length(fwdnei) == 1
            sloccur=[zeros(13-length(max(TSTART+1,t-12):t),1); ...
                squeeze(slsuccess(n,fwdnei,max(TSTART+1,t-12):t)>0)];
        else
            sloccur=[zeros(13-length(max(TSTART+1,t-12):t),length(fwdnei)); ...
                squeeze(slsuccess(n,fwdnei,max(TSTART+1,t-12):t)>0)'];
        end
        [sl_risk,slevnt,tmevnt]=calc_intrisk(sloccur,...
            t_eff,alpharisk,betarisk,timeweight);
        SLRISK(n,fwdnei)=sl_risk;
        if isempty(find(sl_risk,1)) == 0
            avgslrisk(n,t)=mat2cell(SLRISK(n,activeroute{n,t}),1,...
                length(activeroute{n,t}));
        end
        % ICPTL(n,t)=ICPTL(n,t)-OUTFLOW(n,t)*VALUE(  %account for value retained at node
        
        NodeTable.Stock(:)=STOCK(:,t);
        NodeTable.Capital(:)=TOTCPTL(:,t);
    end
    RISKPREM(:,:,t)=max((1-delta_rt).*RISKPREM(:,:,t-1)+...
        delta_rt.*((SLRISK./baserisk(erun)).^riskmltplr(erun)),1);
    
    %%% Make trafficking movie
    MOV(:,n,t)=STOCK(:,t);      % Capture stock data after each node iteration
end

%%% Risk premium on cost of doing business (transport costs)
CTRANS(:,:,t+1)=CTRANS(:,:,t).*RISKPREM(:,:,t);

totslrisk(t+1)=mean(cat(2,avgslrisk{:,t}));
if empSLflag == 0
    %%% Updating interdiction event probability
    subslsuc=slsuccess(:,:,t);
    subslval=slvalue(:,:,t);
    subslprob=SLPROB(:,:,t);
    islcheck=(slevent(:,:,t) == 1);
    subslprob(islcheck)=(1-delta_sl).*subslprob(islcheck)+delta_sl.*...
        (subslsuc(islcheck) > 0);
    SLPROB(:,:,t+1)=subslprob;
end

% Reinforcement learning for successful routes
iactivenode=find(OUTFLOW(2:nnodes,t) > 0);
avgflow=STOCK(iendnode,t)/length(iactivenode);

activenodes=unique(cat(1,activeroute{:,t}));
actedge=activeroute(:,t);

% Calcuate updated marginal profit
for q=1:nnodes
    if isempty(find(ADJ(q,:)==1,1))==1
        continue
    end
    margval(q,q+1:nnodes,t)=PRICE(q+1:nnodes,t)-PRICE(q,t);
end
%%%%%%%%% Route Optimization %%%%%%%%%%%
for dt=1:ndto
    idto=find(NodeTable.DTO == dt);
    DTOBDGT(dt,t)=sum(STOCK(endnodeset,t).*PRICE(endnodeset,t)); %total DTO funds for expansion/viability
    dtorefvec=[1; idto; mexnode];
    subnnodes=length(idto);
    subroutepref=routepref(dtorefvec,dtorefvec,t);
    subactivenodes=activenodes(ismember(activenodes,idto));
    subactedges=cat(1,actedge{dtorefvec});
    ikeep=(NodeTable.DTO(subactedges)==dt);
    dtoACTEDGES=subactedges(ikeep);
    idtoactedges=find(ismember(dtorefvec,dtoACTEDGES)==1);
    subflow=FLOW(dtorefvec,dtorefvec,t);
    dtoslsuc=slsuccess(dtorefvec,dtorefvec,t);
    allflows=subflow+dtoslsuc;
    
    % locate active edges
    [irow,icol]=ind2sub(size(allflows),find(allflows > 0));
    sendedge=ismember(EdgeTable.EndNodes(:,1),dtorefvec);
    dtoEdgeTable=EdgeTable(sendedge,:);
    dtoEdgeTable=dtoEdgeTable(ismember(dtoEdgeTable.EndNodes(:,2),dtorefvec),:);
    dtoSLRISK=SLRISK(dtorefvec,dtorefvec);
    dtoADDVAL=margval(dtorefvec,dtorefvec,t);
    dtoCTRANS=CTRANS(dtorefvec,dtorefvec,t);
    %%% calculate losses from S&L events
    % volume-based - does not matter where in supply chain
    %     supplyfit=STOCK(iendnode,t)/stock_0;
    %     losstolval=losstol*stock_0;
    
    % value-based - price varies with location in supply chain
    ipossl=find(dtoslsuc > 0);
    [nrow,ncol]=ind2sub(size(dtoslsuc),ipossl);
    
    flowvalues=allflows(allflows > 0).*((PRICE(dtorefvec(icol),t)-...
        PRICE(dtorefvec(irow),t))-dtoCTRANS(allflows > 0));

    supplyfit=sum(dtoslsuc(ipossl).*((PRICE(dtorefvec(ncol),t)-...
        PRICE(dtorefvec(nrow),t))-dtoCTRANS(ipossl)));

    losstolval=losstol*max(flowvalues);
    
    if isempty(find(supplyfit ~= 0,1)) == 1 && isempty(find(losstolval ~= 0,1)) == 1
        supplyfit=0.1;
    end
    
    %%% Route capacity constrains flow volumes, need to expand
    %%% routes
    idtonet=dtorefvec(~ismember(dtorefvec,endnodeset));
    if sum(STOCK(idtonet,t)) >= max(dtoEdgeTable.Capacity)
        supplyfit=max(supplyfit,losstolval*sum(STOCK(idtonet,t))/rtcap(erun));   %triggers expansion of one route
    end
    
    %%call top-down route optimization
    expmax=expandmax(erun);
    [newroutepref,newedgechange]=optimizeroute_multidto(dtorefvec,allflows,subflow,supplyfit,expmax,...
        subroutepref,dtoEdgeTable,dtoSLRISK,dtoADDVAL,dtoCTRANS,losstolval,dtoslsuc);
    edgechange(dt)=newedgechange;
    % % Bottom-up route optimization
    % newroutepref=ADJ(dtorefvec,dtorefvec);
    
    routepref(dtorefvec,dtorefvec,t+1)=newroutepref;
    
%     if isempty(find(routepref(1,dtorefvec(1:length(dtorefvec)-1),t+1),1)) ==1
%         display('check route optimization for each dto')
%         keyboard
%     end
%     if isempty(find(newroutepref(1,:),1)) == 1
%         display('lost primary movement')
%     end
end

PRICE(:,t+1)=PRICE(:,t);
if growthmdl(erun) == 1
    STOCK(1,t+1)=stock_0+(prodgrow(erun)*ceil((t-TSTART)/12));    %additional production to enter network next time step
elseif growthmdl(erun) == 2
    STOCK(1,t+1)=(stock_max*stock_0*exp(prodgrow(erun)*floor(t/12)))/...
        (stock_max+stock_0*(exp(prodgrow(erun)*floor(t/12))-1));
end

STOCK(endnodeset,t+1)=0;    %remove stock at end node for next time step
NodeTable.Stock(1)=STOCK(1,t+1);
NodeTable.Stock(endnodeset)=0;

slcount_edges(t)=length(find(slsuccess(:,:,t) > 0));
h_slsuccess=slsuccess(:,:,t);
slcount_vol(t)=sum(h_slsuccess(h_slsuccess > 0));

%%%Output tables for flows(t) and interdiction prob(t-1)
[Tflow,Tintrd]=intrd_tables_batch(FLOW,slsuccess,SLPROB,NodeTable,EdgeTable,t,testflag,erun,mrun,batchrun);
Tflow=table2struct(Tflow,'ToScalar',true);
clear t
wrksp_name='NarcoLogic_wrksp.mat';
save(wrksp_name)