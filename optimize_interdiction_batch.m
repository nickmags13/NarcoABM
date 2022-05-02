%%%%%%%%%%   Interdiction events from optimization model   %%%%%%%%%%%%%%%%
function [intrdct_events,intrdct_nodes]=optimize_interdiction_batch(t,ADJ,testflag,erun,mrun,batchrun)

if testflag  == 1
%     cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\ABM_test_files\MCI_3_061020
    cd C:\Users\pcbmi\Box\NSF_D-ISN\Code\NarcoLogic\MTMCI_IntNodes
else
%     cd C:\Users\nrmagliocca\'Google Drive'\NSF_EAGER_Models\INT_Nodes
%     cd D:\NSF_EAGER_Models\INT_Nodes
    if batchrun == 1
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_1
    elseif batchrun == 2
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_2
    elseif batchrun == 3
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_3
    elseif batchrun == 4
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_4
    elseif batchrun == 5
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_5
    elseif batchrun == 6
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_6
    elseif batchrun == 7
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_7
    elseif batchrun == 8
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_8
    elseif batchrun == 9
        cd C:\Users\Penelope\Box\NSF_D-ISN\Data\IntegratedModels\INT_Nodes
    elseif batchrun == 10
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_10
    elseif batchrun == 11
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_11
    elseif batchrun == 12
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_12
    elseif batchrun == 13
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_13
    elseif batchrun == 14
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_14
    elseif batchrun == 15
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_15
    elseif batchrun == 16
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_16
    elseif batchrun == 17
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_17
    elseif batchrun == 18
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_18
    elseif batchrun == 19
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_19
    elseif batchrun == 20
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_20
    elseif batchrun == 21
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_21
    elseif batchrun == 22
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_22
    elseif batchrun == 23
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_23
    elseif batchrun == 24
        cd D:\NSF_EAGER_Models\INT_Nodes\INT_Nodes_24
    end
end

readflag=0;
% % fnames=dir;
% % fnamescell=struct2cell(fnames);
% t1=double(t-1 >= 100);
% if t>= 101
%     t2=floor((t-101)/10);
% else
%     t2=floor((t-1)/10);
% end
% t3=mod(t-1,10);
% mrun_t1=floor(mrun/10);
% mrun_t2=mod(mrun,10);
% erun_t1=floor(erun/100);
% erun_t2=floor(erun/10);
% erun_t3=mod(erun,10);
% trgtfile=sprintf('%d%d%d_INT_MTMCI_RUN_%d%d%d%d%d.txt',t1,t2,t3,mrun_t1,mrun_t2,erun_t1,erun_t2,erun_t3);
% trgtfile='C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\MTMCI_IntNodes\MTMCI_IntNodes.txt';
trgtfile='MTMCI_IntNodes.txt';
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
% cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM
cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic