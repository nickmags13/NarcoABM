%%%%%%% Top-down supply chain optimization %%%%%%%%
% function newroutepref=optimizeroute(nnodes,subflow,supplyfit,activenodes,...
%     subroutepref,EdgeTable,SLRISK,ADDVAL,CTRANS,losstolval)
function newroutepref=optimizeroute_multidto(dtorefvec,subflow,supplyfit,subactivenodes,...
            subroutepref,dtoEdgeTable,dtoSLRISK,dtoADDVAL,dtoCTRANS,losstolval,dtoslsuc)

% allnodes=2:subnnodes;
iactiveedges=find(subflow > 0 | dtoslsuc > 0);
[actrow,actcol]=ind2sub(size(subflow),iactiveedges);
edgeparms=[subflow(iactiveedges) dtoSLRISK(iactiveedges) iactiveedges actrow actcol];
if supplyfit <= losstolval  %need to consolidate supply chain
    edgesort=sortrows(edgeparms,-2); %refernce this array when removing edges
    %    edgecut=1:min(length(iactiveedges)-ceil(length(iactiveedges)*...
    %        ((losstolval-supplyfit)/losstolval)),length(iactiveedges)-1);   %calc how many edges need to be removed
    iprimary=find(edgesort(:,4) == 1 & edgesort(:,5) ~= length(dtorefvec));   %primary movement
    edgecut=1:min(round(length(iactiveedges)*(supplyfit/...
        (supplyfit+losstolval))),length(iactiveedges)-1);
    %    iedgecut=edgecut(edgesort(edgecut,2)> min(edgeparms(:,2)));
    %    if length(iedgecut) < length(edgecut)    % are there less high risk edges than need to be removed?
    %        %remove as many high risk edges that are greater than minimum risk
    %        %edges, and then remove lowest volume edges
    %        if isempty(find(iedgecut,1)) == 1
    %            secondsort=sortrows(edgeparms,1);
    %        else
    %            secondsort=sortrows(edgeparms(~ismember(edgeparms(:,3),...
    %                edgesort(edgecut(iedgecut),3)),:),1);    %sort based on volume (low to high)
    %        end
    %        inewcuts=find(ismember(edgesort(:,3),secondsort(1:length(edgecut)-...
    %            length(iedgecut),3)) == 1);    %match based on unique iactiveedges
    %        edgecut=[edgecut(iedgecut)'; inewcuts];
    %    end
    
    %%% Preserve at least one primary movement
    
    minrisk_primary=min(edgesort(iprimary,2));
    ikeep_primary=find(edgesort(iprimary,2) == minrisk_primary);
    if length(ikeep_primary) == 1
        edgecut=edgecut(~ismember(edgecut,[iprimary(ikeep_primary); ...
            find(edgesort(edgecut,4)==edgesort(iprimary(ikeep_primary),5))]));
    else
        maxprofit_primary=max(edgesort(iprimary(ikeep_primary),1));
        ikeep_primary=ikeep_primary(edgesort(iprimary(ikeep_primary),1) == ...
            maxprofit_primary);
        if length(ikeep_primary) == 1
            edgecut=edgecut(~ismember(edgecut,[iprimary(ikeep_primary); ...
                find(edgesort(edgecut,4)==edgesort(iprimary(ikeep_primary),5))]));
        else
            ikeep_primary=ikeep_primary(ceil(length(ikeep_primary)*rand(1)));
            edgecut=edgecut(~ismember(edgecut,[iprimary(ikeep_primary); ...
                find(edgesort(edgecut,4)==edgesort(iprimary(ikeep_primary),5))]));
        end
    end
    
    % remove highest risk edges
    for j=1:length(edgecut)
       icheckroute=find(subflow(edgesort(edgecut(j),4),...
           ismember(dtorefvec,dtoEdgeTable.EndNodes(dtoEdgeTable.EndNodes(:,1)==...
           dtorefvec(edgesort(edgecut(j),4)),2))) > 0);
%        icheckroute=find(subflow(edgesort(edgecut(j),4),ismember(dtorefvec,...
%            dtoEdgeTable.EndNodes(ismember(dtoEdgeTable.EndNodes(:,1),...
%            cuttable(:,1)),2))) > 0); % check if affected node has more than one active recipient
       actroutes=dtoEdgeTable.EndNodes(dtoEdgeTable.EndNodes(:,1)==...
           dtorefvec(edgesort(edgecut(j),4)),1);
%        icheckroute=find(subflow(edgesort(edgecut(j),4),...
%            dtoEdgeTable.EndNodes(dtoEdgeTable.EndNodes(:,1)==...
%            edgesort(edgecut(j),4),2)) > 0);
%        actroutes=dtoEdgeTable.EndNodes(dtoEdgeTable.EndNodes(:,1)==edgesort(edgecut(j),4),1);
       checknoderoutes=(length(actroutes(icheckroute)) == length(find(edgesort(edgecut,4) == ...
           edgesort(edgecut(j),4))));  % check if all acctive edges will be removed
       if checknoderoutes == 1
           cutsenders=find(ismember(dtorefvec,dtoEdgeTable.EndNodes(...
               dtoEdgeTable.EndNodes(:,2)==dtorefvec(edgesort(edgecut(j),4)),1)) == 1);
%            cutsenders=find(ismember(dtorefvec,dtoEdgeTable.EndNodes(...
%                ismember(dtoEdgeTable.EndNodes(:,2),cuttable(:,2)),1)) == 1);
%            cutsenders=dtoEdgeTable.EndNodes(dtoEdgeTable.EndNodes(:,2)==...
%                edgesort(edgecut(j),4),1);
           cutind=sub2ind(size(dtoSLRISK),cutsenders,edgesort(edgecut(j),4)*...
               ones(length(cutsenders),1));
           subroutepref(cutind)=0;
       end
    
       if length(icheckroute) == 1
           subroutepref(edgesort(edgecut(j),3))=0; %remove high risk edge (only edge) from node
           irmvsender=(edgesort(:,5) == edgesort(edgecut(j),4));    %find nodes sending to affected node
           subroutepref(edgesort(irmvsender,3))=0; %remove sending edges
       else
           subroutepref(edgesort(edgecut(j),3))=0;
       end
   end
elseif supplyfit > losstolval    %need to expand supply chain
%     potnodes=allnodes(~ismember(allnodes,subactivenodes));
    potnodes=dtorefvec(~ismember(dtorefvec,[1; subactivenodes; ...
        dtorefvec(length(dtorefvec))]));
%     edgeadd=1:min(max(ceil(length(subactivenodes)*(1+(supplyfit-losstolval)/...
%         supplyfit))-length(subactivenodes),1),length(potnodes));
    edgeadd=1:min(max(ceil((supplyfit-losstolval)/supplyfit),1),length(potnodes));

    if isempty(find(potnodes,1)) == 1
%         display('No more nodes to expand')
%         subroutepref(iactiveedges)=1-SLRISK(iactiveedges);
        % could introduce new edges by modifying ADJ?
    else
        newedgeparms=[];
        for k=1:length(potnodes)
            potsenders=unique(dtoEdgeTable.EndNodes(ismember(dtoEdgeTable.EndNodes(:,2),...
                potnodes(k)),1));
            potsenders=potsenders(ismember(potsenders,[1; subactivenodes]));
            ipotsenders=find(ismember(potsenders,dtorefvec) == 1);
%             ipotedge=sub2ind(size(dtoSLRISK),potsenders,...
%                 potnodes(k)*ones(length(potsenders),1));
%             newedgeparms=[newedgeparms; (dtoADDVAL(potsenders,potnodes(k))-...
%                dtoCTRANS(potsenders,potnodes(k))).*(1-dtoSLRISK(potsenders,potnodes(k))) ...
%                 dtoSLRISK(potsenders,potnodes(k)) ipotedge k*ones(length(potsenders),1)];
            ipotedge=sub2ind(size(dtoSLRISK),find(ismember(dtorefvec,...
                potsenders(ipotsenders)) == 1),...
                find(dtorefvec == potnodes(k))*ones(length(potsenders),1));
            newedgeparms=[newedgeparms; (dtoADDVAL(ipotsenders,k)-...
               dtoCTRANS(ipotsenders,k)).*(1-dtoSLRISK(ipotsenders,k)) ...
                dtoSLRISK(ipotsenders,k) ipotedge k*ones(length(ipotsenders),1)];
        end
%         edgesort=sortrows(newedgeparms,2); %refernce this array when removing edges
        edgesort=sortrows(newedgeparms,-1);
        subroutepref(edgesort(edgeadd,3))=1;
        ireceivers=dtoEdgeTable.EndNodes(ismember(dtoEdgeTable.EndNodes(:,1),...
                unique(potnodes(edgesort(edgeadd,4)))),:);
%         addind=sub2ind(size(dtoSLRISK),ireceivers(:,1),ireceivers(:,2));
        send_row=[];
        rec_col=[];
        for jj=1:length(ireceivers(:,1))
            send_row=[send_row; find(ismember(dtorefvec,ireceivers(jj,1))==1)];
            rec_col=[rec_col; find(ismember(dtorefvec,ireceivers(jj,2))==1)];
        end
        addind=sub2ind(size(dtoSLRISK),send_row,rec_col);
        subroutepref(addind)=1;
        subroutepref(rec_col,length(dtorefvec))=1;
%         if length(unique(edgesort(edgeadd,4))) == length(potnodes)
%             % if edge is selected to send, then at least one of the sending
%             % edges of the new node must be activated
%             
%         else
%             subroutepref(edgesort(edgeadd,3))=1;
%         end
    end    
end
newroutepref=subroutepref;