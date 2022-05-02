%%%%%% Interdiction Initialization %%%%%%%%%%
function [Tflow,Tintrd]=intrd_tables_batch(FLOW,slsuccess,SLPROB,NodeTable,EdgeTable,t,testflag,erun,mrun,batchrun)
%function [TflowTest,Tintrd]=intrd_tables_batch(FLOW,slsuccess,SLPROB,NodeTable,EdgeTable,t,testflag,erun,mrun,batchrun)

varTypes={'double','double','double','double'};
Tflow=table('Size',[height(EdgeTable) 4],'VariableTypes',varTypes,...
    'VariableNames',{'End_Node','Start_Node','IntitFlow','DTO'});
startFLOW=FLOW(:,:,t)+slsuccess(:,:,t);
% intcptFLOW=FLOW(:,:,t).*repmat(NodeTable.pintcpt',height(NodeTable),1);
% startFLOW=intcptFLOW(:,:,t)+slsuccess(:,:,t);

Tintrd=table('Size',[height(EdgeTable) 3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'End_Node','Start_Node','IntitProb'});
if t == 1
    startSLPROB=SLPROB(:,:,1);
else
    startSLPROB=SLPROB(:,:,t-1);
end
sumprob=sum(sum(startSLPROB));
for g=1:height(EdgeTable)
    edge=table2array(EdgeTable(g,1));
    Tflow(g,1)={edge(2)};
    Tflow(g,2)={edge(1)};
    Tflow(g,3)={startFLOW(edge(1),edge(2))};
    Tflow(g,4)={NodeTable.DTO(edge(2))};
    
    Tintrd(g,1)={edge(2)};
    Tintrd(g,2)={edge(1)};
    Tintrd(g,3)={startSLPROB(edge(1),edge(2))/sumprob};
%     iedge=find(startSLPROB > 0);
%     tindex=height(Tintrd(:,1))+1:height(Tintrd(:,1))+length(iedge);
%     Tintrd(tindex,1)=array2table(iedge');
%     Tintrd(tindex,2)=array2table(ones(length(iedge),1)*i);
%     Tintrd(tindex,3)=array2table(SLPROB(i,iedge)');
end
t1=double(t >= 100);
if t>= 100
    t2=floor((t-100)/10);
else
    t2=floor(t/10);
end
t3=mod(t,10);

mrun_t1=floor(mrun/10);
mrun_t2=mod(mrun,10);
erun_t1=floor(erun/100);
erun_t2=floor(erun/10);
erun_t3=mod(erun,10);

if testflag == 0
    %     cd C:\Users\Penelope\Box\NSF_D-ISN\Data\IntegratedModels\ABM_Results
    if batchrun == 1
        cd C:\Users\pcbmi\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\ABM_Results_1
    elseif batchrun == 2
        cd C:\Users\pcbmi\Box\NSF_D-ISN\Data\IntegratedModels
    elseif batchrun == 3
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_3
    elseif batchrun == 4
        cd C:\Users\Penelope\Box\NSF_D-ISN\Data\IntegratedModels\ABM_Results
    elseif batchrun == 5
        cd C:\Users\Penelope\Box\NSF_D-ISN\Data\IntegratedModels\ABM_Results
    elseif batchrun == 6
        cd C:\Users\Penelope\Box\NSF_D-ISN\Data\IntegratedModels\ABM_Results
    elseif batchrun == 7
        cd C:\Users\Penelope\Box\NSF_D-ISN\Data\IntegratedModels\ABM_Results
    elseif batchrun == 8
        cd C:\Users\Penelope\Box\NSF_D-ISN\Data\IntegratedModels\ABM_Results
    elseif batchrun == 9
        cd C:\Users\Penelope\Box\NSF_D-ISN\Data\IntegratedModels\ABM_Results
    elseif batchrun == 10
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_10
    elseif batchrun == 11
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_11
    elseif batchrun == 12
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_12
    elseif batchrun == 13
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_13
    elseif batchrun == 14
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_14
    elseif batchrun == 15
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_15
    elseif batchrun == 16
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_16
    elseif batchrun == 17
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_17
    elseif batchrun == 18
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_18
    elseif batchrun == 19
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_19
    elseif batchrun == 20
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_20
    elseif batchrun == 21
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_21
    elseif batchrun == 22
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_22
    elseif batchrun == 23
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_23
    elseif batchrun == 24
        cd D:\NSF_EAGER_Models\ABM_Results\ABM_Results_24
    end
    
    writetable(Tflow,sprintf('%d%d%d_SL_MTMCI_RUN_%d%d%d%d%d.txt',t1,t2,t3,mrun_t1,mrun_t2,erun_t1,erun_t2,erun_t3))
%     writetable(Tintrd,sprintf('%d%d%d_SLPROB_3_MCI.txt',t1,t2,t3)) %apply to next time step's interdiction actions
end
% cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM
cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic