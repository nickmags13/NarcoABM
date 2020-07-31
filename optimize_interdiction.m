%%%%%%%%%%   Interdiction events from optimization model   %%%%%%%%%%%%%%%%
function [intrdct_events,intrdct_nodes]=optimize_interdiction(t,ADJ,testflag)

if testflag  == 1
    cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\ABM_test_files\MCI_3_061020
else
    cd C:\Users\nrmagliocca\'Google Drive'\NSF_EAGER_Models\INT_Nodes
end

readflag=0;
% fnames=dir;
% fnamescell=struct2cell(fnames);
t1=double(t-1 >= 100);
if t>= 101
    t2=floor((t-101)/10);
else
    t2=floor((t-1)/10);
end
t3=mod(t-1,10);
% t=t-2;
% t1=double(t >= 100);
% if t>= 100
%     t2=floor((t-100)/10);
% else
%     t2=floor(t/10);
% end
% t3=mod(t,10);
trgtfile=sprintf('%d%d%d_INT_3_MCI.txt',t1,t2,t3);
sprintf('Looking for %s',trgtfile)
while readflag == 0
   fnames=dir;
   fnamescell=struct2cell(fnames);
   h=strcmp(trgtfile,fnamescell(1,:));
   hind=find(h == 1);
   readflag=~isempty(find(hind,1));
end
sprintf('Interdiction Input File Success, t=%d',t)
Tintevent=readmatrix(trgtfile);
intrdct_events=zeros(size(ADJ));
intrdct_nodes=Tintevent;
for j=1:length(Tintevent)
    iupstream=(ADJ(:,Tintevent(j))==1);
    intrdct_events(iupstream,Tintevent(j))=1;
end
cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM