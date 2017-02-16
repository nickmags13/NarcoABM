%%%%%%% Top-down supply chain optimization %%%%%%%%
function newroutepref=optimizeroute(nnodes,subflow,supplyfit,activenodes,...
    subroutepref,EdgeTable,SLRISK,ADDVAL,losstol)

allnodes=2:nnodes;
iactiveedges=find(subflow > 0);
[actrow,actcol]=ind2sub(size(subflow),iactiveedges);
edgeparms=[subflow(iactiveedges) SLRISK(iactiveedges) iactiveedges actrow actcol];
if supplyfit >= losstol   %need to consolidate supply chain
   edgesort=sortrows(edgeparms,-2); %refernce this array when removing edges
   edgecut=1:length(iactiveedges)-ceil(length(iactiveedges)*supplyfit);   %calc how many edges need to be removed
   iedgecut=edgecut(edgesort(edgecut,2)> min(edgeparms(:,2)));
   if length(iedgecut) < length(edgecut)    % are there less high risk nodes than need to be removed?
       %remove as many high risk edges that are greater than minimum risk
       %edges, and then remove lowest volume edges
       if isempty(find(iedgecut,1)) == 1
           secondsort=sortrows(edgeparms,1);
       else
           secondsort=sortrows(edgeparms(edgeparms(:,3)~=...
               edgesort(edgecut(iedgecut),3),:),1);    %sort based on volume (low to high)
       end
       inewcuts=find(ismember(edgesort(:,3),secondsort(1:length(edgecut)-...
           length(iedgecut),3)) == 1);    %match based on unique iactiveedges
       edgecut=[edgecut(iedgecut); inewcuts]; 
   end
   % remove highest risk edges
   for j=1:length(edgecut)
       checkroute=find(edgesort(:,4) == edgesort(edgecut(j),4)); % check if sending node has more than one recipient
       if length(checkroute) == 1
           subroutepref(edgesort(edgecut(j),3))=0; %remove high risk edge (only edge) from node
           irmvsender=(edgesort(:,5) == edgesort(edgecut(j),4));    %find nodes sending to affected node
           subroutepref(edgesort(irmvsender,3))=0; %remove sending edges
       else
           subroutepref(edgesort(edgecut(j),3))=0;
       end
   end
elseif supplyfit < losstol    %need to expand supply chain
    edgeadd=1:floor(length(activenodes)*(1+(losstol-supplyfit)))-...
        length(activenodes);
    potnodes=allnodes(~ismember(allnodes,activenodes));
    if length(potnodes) < length(edgeadd) || isempty(find(potnodes,1)) == 1
        display('No more nodes to expand')
        % could introduce new edges by modifying ADJ?
    else
        newedgeparms=[];
        for k=1:length(potnodes)
            potsenders=unique(EdgeTable.EndNodes(ismember(EdgeTable.EndNodes(:,2),...
                potnodes(k)),1));
            potsenders=potsenders(ismember(potsenders,[1; activenodes]));
            ipotedge=sub2ind(size(SLRISK),potsenders,potnodes(k)*...
                ones(length(potsenders),1));
            newedgeparms=[newedgeparms; ADDVAL(potsenders,potnodes(k)) ...
                SLRISK(potsenders,potnodes(k)) ipotedge];
        end
        edgesort=sortrows(newedgeparms,2); %refernce this array when removing edges
        subroutepref(edgesort(edgeadd,3))=1;
    end    
end
newroutepref=subroutepref;