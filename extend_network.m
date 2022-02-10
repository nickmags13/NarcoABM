function [ext_NodeTable,ext_EdgeTable]=extend_network(nnodes,NodeTable,EdgeTable,Rdptgrid)
warning('off')
ext_NodeTable=NodeTable;
ext_EdgeTable=EdgeTable;
%%%add intermediate nodes
%Caribbean
inn1=nnodes+1;
ext_NodeTable(inn1,1)={inn1};
ext_NodeTable.Lat(inn1)=15.36;
ext_NodeTable.Lon(inn1)=-78.9189;
[noderow,nodecol]=latlon2pix(Rdptgrid,ext_NodeTable.Lat(inn1),ext_NodeTable.Lon(inn1));
ext_NodeTable.Row(inn1)=round(noderow);
ext_NodeTable.Col(inn1)=round(nodecol);
ext_NodeTable.DeptCode(inn1)=3;
ext_NodeTable(inn1,7:width(ext_NodeTable))=ext_NodeTable(nnodes,7:width(ext_NodeTable));
ext_NodeTable.DTO(inn1)=2;
ext_NodeTable.CoastDist(inn1)=0;
ext_NodeTable.pintcpt(inn1)=0.01;

%north of Galapagos
inn2=nnodes+2;
ext_NodeTable(inn2,1)={inn2};
ext_NodeTable.Lat(inn2)=2.442;
ext_NodeTable.Lon(inn2)=-88.1542;
%             [noderow,nodecol]=latlon2pix(Rdptgrid,ext_NodeTable.Lat(inn2),ext_NodeTable.Lon(inn2));
ext_NodeTable.Row(inn2)=4800;
ext_NodeTable.Col(inn2)=100;
ext_NodeTable.DeptCode(inn2)=3;
ext_NodeTable(inn2,7:width(ext_NodeTable))=ext_NodeTable(nnodes,7:width(ext_NodeTable));
ext_NodeTable.DTO(inn2)=1;
ext_NodeTable.CoastDist(inn2)=0;
ext_NodeTable.pintcpt(inn2)=0.01;

%south of Galapagos
inn3=nnodes+3;
ext_NodeTable(inn3,1)={inn3};
ext_NodeTable.Lat(inn3)=-3.0039;
ext_NodeTable.Lon(inn3)=-93.1533;
%             [noderow,nodecol]=latlon2pix(Rdptgrid,ext_NodeTable.Lat(inn3),ext_NodeTable.Lon(inn3));
ext_NodeTable.Row(inn3)=4914;
ext_NodeTable.Col(inn3)=1;
ext_NodeTable.DeptCode(inn3)=3;
ext_NodeTable(inn3,7:width(ext_NodeTable))=ext_NodeTable(nnodes,7:width(ext_NodeTable));
ext_NodeTable.DTO(inn3)=1;
ext_NodeTable.CoastDist(inn3)=0;
ext_NodeTable.pintcpt(inn3)=0.01;

%%% add end nodes
%Dominican Rep.
inn4=nnodes+4;
ext_NodeTable(inn4,1)={inn4};
ext_NodeTable.Lat(inn4)=18.3383;
ext_NodeTable.Lon(inn4)=-70.3933;
[noderow,~]=latlon2pix(Rdptgrid,ext_NodeTable.Lat(inn4),ext_NodeTable.Lon(inn4));
ext_NodeTable.Row(inn4)=round(noderow);
ext_NodeTable.Col(inn4)=6537;
ext_NodeTable.DeptCode(inn4)=3;
ext_NodeTable(inn4,7:width(ext_NodeTable))=ext_NodeTable(nnodes,7:width(ext_NodeTable));
ext_NodeTable.DTO(inn4)=2;
ext_NodeTable.CoastDist(inn4)=0;
ext_NodeTable.pintcpt(inn4)=0.01;

%Yucatan
inn5=nnodes+5;
ext_NodeTable(inn5,1)={inn5};
ext_NodeTable.Lat(inn5)=19.0292;
ext_NodeTable.Lon(inn5)=-87.8444;
[~,nodecol]=latlon2pix(Rdptgrid,ext_NodeTable.Lat(inn5),ext_NodeTable.Lon(inn5));
ext_NodeTable.Row(inn5)=1;
ext_NodeTable.Col(inn5)=round(nodecol);
ext_NodeTable.DeptCode(inn5)=3;
ext_NodeTable(inn5,7:width(ext_NodeTable))=ext_NodeTable(nnodes,7:width(ext_NodeTable));
ext_NodeTable.DTO(inn5)=2;
ext_NodeTable.CoastDist(inn5)=0;
ext_NodeTable.pintcpt(inn5)=0;

%Chiapas
inn6=nnodes+6;
ext_NodeTable(inn6,1)={inn6};
ext_NodeTable.Lat(inn6)=16.0167;
ext_NodeTable.Lon(inn6)=-93.5794;
[noderow,~]=latlon2pix(Rdptgrid,ext_NodeTable.Lat(inn6),ext_NodeTable.Lon(inn6));
ext_NodeTable.Row(inn6)=round(noderow);
ext_NodeTable.Col(inn6)=50;
ext_NodeTable.DeptCode(inn6)=3;
ext_NodeTable(inn6,7:width(ext_NodeTable))=ext_NodeTable(nnodes,7:width(ext_NodeTable));
ext_NodeTable.DTO(inn6)=1;
ext_NodeTable.CoastDist(inn6)=0;
ext_NodeTable.pintcpt(inn6)=0;

%Oaxaca
inn7=nnodes+7;
ext_NodeTable(inn7,1)={inn7};
ext_NodeTable.Lat(inn7)=15.8114;
ext_NodeTable.Lon(inn7)=-96.5508;
[noderow,~]=latlon2pix(Rdptgrid,ext_NodeTable.Lat(inn7),ext_NodeTable.Lon(inn7));
ext_NodeTable.Row(inn7)=round(noderow);
ext_NodeTable.Col(inn7)=1;
ext_NodeTable.DeptCode(inn7)=3;
ext_NodeTable(inn7,7:width(ext_NodeTable))=ext_NodeTable(nnodes,7:width(ext_NodeTable));
ext_NodeTable.DTO(inn7)=1;
ext_NodeTable.CoastDist(inn7)=0;
ext_NodeTable.pintcpt(inn7)=0;

%%% Add custom edges
%source to Carib
ine1=height(EdgeTable)+1;
ext_EdgeTable.EndNodes(ine1,:)=[1 ext_NodeTable.ID(inn1)];
ext_EdgeTable(ine1,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%Carib to default Mexico(156)
ine2=height(EdgeTable)+2;
ext_EdgeTable.EndNodes(ine2,:)=[ext_NodeTable.ID(inn1) 156];
ext_EdgeTable(ine2,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%Carib to Yucatan
ine3=height(EdgeTable)+3;
ext_EdgeTable.EndNodes(ine3,:)=[ext_NodeTable.ID(inn1) ext_NodeTable.ID(inn5)];
ext_EdgeTable(ine3,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%source to Domincan Republic
ine4=height(EdgeTable)+4;
ext_EdgeTable.EndNodes(ine4,:)=[1 ext_NodeTable.ID(inn4)];
ext_EdgeTable(ine4,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%source to north of Galapagos
ine5=height(EdgeTable)+5;
ext_EdgeTable.EndNodes(ine5,:)=[1 ext_NodeTable.ID(inn2)];
ext_EdgeTable(ine5,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%north of Galapagos to Chiapas
ine6=height(EdgeTable)+6;
ext_EdgeTable.EndNodes(ine6,:)=[ext_NodeTable.ID(inn2) ext_NodeTable.ID(inn6)];
ext_EdgeTable(ine6,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%north of Galapagos to Oaxaca
ine7=height(EdgeTable)+7;
ext_EdgeTable.EndNodes(ine7,:)=[ext_NodeTable.ID(inn2) ext_NodeTable.ID(inn7)];
ext_EdgeTable(ine7,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%source to south of Galapagos
ine8=height(EdgeTable)+8;
ext_EdgeTable.EndNodes(ine8,:)=[1 ext_NodeTable.ID(inn3)];
ext_EdgeTable(ine8,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%south of Galapagos to Chiapas
ine9=height(EdgeTable)+9;
ext_EdgeTable.EndNodes(ine9,:)=[ext_NodeTable.ID(inn3) ext_NodeTable.ID(inn6)];
ext_EdgeTable(ine9,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%south of Galapagos to Oaxaca
ine10=height(EdgeTable)+10;
ext_EdgeTable.EndNodes(ine10,:)=[ext_NodeTable.ID(inn3) ext_NodeTable.ID(inn7)];
ext_EdgeTable(ine10,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

%Dom. Rep. to US/Mexico
ine11=height(EdgeTable)+11;
ext_EdgeTable.EndNodes(ine11,:)=[ext_NodeTable.ID(inn4) 156];
ext_EdgeTable(ine11,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));

% % Add Carib to HND and GT nodes
% deptlist=caadmid1(adm1_0 == 111 | adm1_0 == 103);   %GT and HND
% iadd=find(ismember(NodeTable.DeptCode,deptlist)==1);
% ine12=height(EdgeTable)+12:height(EdgeTable)+12+length(iadd);
% ext_EdgeTable.EndNodes(ine12,:)=[ones(length(iadd),1)*ext_NodeTable.ID(inn1) NodeTable.ID(iadd)];
% ext_EdgeTable(ine12,2:width(ext_EdgeTable))=ext_EdgeTable(156,2:width(ext_EdgeTable));
warning('on')