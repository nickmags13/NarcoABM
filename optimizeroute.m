%%%%%%% Top-down supply chain optimization %%%%%%%%
function newroutepref=optimizeroute(nnodes,subflow,supplyfit,activenodes,...
    subroutepref,EdgeTable,SLRISK,ADDVAL,CTRANS,losstolval)

allnodes=2:nnodes;
iactiveedges=find(subflow > 0);
[actrow,actcol]=ind2sub(size(subflow),iactiveedges);
edgeparms=[subflow(iactiveedges) SLRISK(iactiveedges) iactiveedges actrow actcol];
if supplyfit >= losstolval   %need to consolidate supply chain
   edgesort=sortrows(edgeparms,-2); %refernce this array when removing edges
   edgecut=1:length(iactiveedges)-ceil(length(iactiveedges)*...
       ((supplyfit-losstolval)/losstolval));   %calc how many edges need to be removed
   iedgecut=edgecut(edgesort(edgecut,2)> min(edgeparms(:,2)));
   if length(iedgecut) < length(edgecut)    % are there less high risk edges than need to be removed?
       %remove as many high risk edges that are greater than minimum risk
       %edges, and then remove lowest volume edges
       if isempty(find(iedgecut,1)) == 1
           secondsort=sortrows(edgeparms,1);
       else
           secondsort=sortrows(edgeparms(~ismember(edgeparms(:,3),...
               edgesort(edgecut(iedgecut),3)),:),1);    %sort based on volume (low to high)
       end
       inewcuts=find(ismember(edgesort(:,3),secondsort(1:length(edgecut)-...
           length(iedgecut),3)) == 1);    %match based on unique iactiveedges
       edgecut=[edgecut(iedgecut)'; inewcuts]; 
   end
   % remove highest risk edges
   for j=1:length(edgecut)
       icheckroute=find(subflow(edgesort(edgecut(j),4),EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==...
           edgesort(edgecut(j),4),2)) > 0); % check if affected node has more than one active recipient
       actroutes=EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==edgesort(edgecut(j),4),1);
       checknoderoutes=(length(actroutes(icheckroute)) == length(find(edgesort(edgecut,4) == ...
           edgesort(edgecut(j),4))));  % check if all acctive edges will be removed
       if checknoderoutes == 1
           cutsenders=EdgeTable.EndNodes(EdgeTable.EndNodes(:,2)==...
               edgesort(edgecut(j),4),1);
           cutind=sub2ind(size(SLRISK),cutsenders,edgesort(edgecut(j),4)*...
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
elseif supplyfit < losstolval    %need to expand supply chain
    potnodes=allnodes(~ismember(allnodes,activenodes));
    edgeadd=1:min(max(floor(length(activenodes)*(1+(losstolval-supplyfit)/...
        losstolval))-length(activenodes),1),length(potnodes));
%     if length(potnodes) < length(edgeadd) || isempty(find(potnodes,1)) == 1
    if isempty(find(potnodes,1)) == 1
        display('No more nodes to expand')
        % could introduce new edges by modifying ADJ?
    else
        newedgeparms=[];
        for k=1:length(potnodes)
            potsenders=unique(EdgeTable.EndNodes(ismember(EdgeTable.EndNodes(:,2),...
                potnodes(k)),1));
%             potreceivers=unique(EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==...
%                 potnodes(k),2));
            potsenders=potsenders(ismember(potsenders,[1; activenodes]));
%             potreceivers=potreceivers(ismember(potreceivers,activenodes));
            ipotedge=sub2ind(size(SLRISK),potsenders,...
                potnodes(k)*ones(length(potsenders),1));
            newedgeparms=[newedgeparms; (ADDVAL(potsenders,potnodes(k))-...
                CTRANS(potsenders,potnodes(k))).*(1-SLRISK(potsenders,potnodes(k))) ...
                SLRISK(potsenders,potnodes(k)) ipotedge k*ones(length(potsenders),1)];
%             ipotedge=sub2ind(size(SLRISK),[potsenders; potreceivers],...
%                 potnodes(k)*ones(length([potsenders; potreceivers]),1));
%             newedgeparms=[newedgeparms; [ADDVAL(potsenders,potnodes(k)); ...
%                 ADDVAL(potnodes(k),potreceivers)'] [SLRISK(potsenders,potnodes(k));...
%                 SLRISK(potnodes(k),potreceivers)']  ipotedge ...
%                 k*ones(length([potsenders; potreceivers]),1)];
        end
%         edgesort=sortrows(newedgeparms,2); %refernce this array when removing edges
        edgesort=sortrows(newedgeparms,-1);
        subroutepref(edgesort(edgeadd,3))=1;
        ireceivers=EdgeTable.EndNodes(ismember(EdgeTable.EndNodes(:,1),...
                unique(potnodes(edgesort(edgeadd,4)))),:);
        addind=sub2ind(size(SLRISK),ireceivers(:,1),ireceivers(:,2));
        subroutepref(addind)=1;
        subroutepref(ireceivers(:,2),nnodes)=1;
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