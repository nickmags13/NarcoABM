function parsave_illicit_supplychain_batch(savefname,EdgeTable,NodeTable,MOV,...
    FLOW,OUTFLOW,CTRANS,TOTCPTL,DTOBDGT,slsuccess,activeroute,STOCK,RISKPREM,...
    slperevent,slval,nactnodes,sltot,t_firstmov,PRICE,slcount_edges,slcount_vol,slnodes,batchrun)

% cd C:\Users\nrmagliocca\'Box Sync'\'Data Drive'\model_results\SupplyChain_full_021618
% cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\model_results\SupplyChain_optint_081420
% cd C:\Users\'HEIMA Lab'\Documents\MATLAB\NarcoLogic\model_results\SupplyChain_optint_testruns
if batchrun == 1
    cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\SupplyChain_optint_batch_1
elseif batchrun == 2
    cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\SupplyChain_optint_batch_2
elseif batchrun == 3
    cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\SupplyChain_optint_batch_3
elseif batchrun == 4
    cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\SupplyChain_optint_batch_4
elseif batchrun == 5
    cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\SupplyChain_optint_batch_5
elseif batchrun == 6
    cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\SupplyChain_optint_batch_6
elseif batchrun == 7
    cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\SupplyChain_optint_batch_7
elseif batchrun == 8
    cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\SupplyChain_optint_batch_8
elseif batchrun == 9
    cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic\TrialResults\SupplyChain_optint_batch_9
elseif batchrun == 10
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_10
elseif batchrun == 11
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_11
elseif batchrun == 12
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_12
elseif batchrun == 13
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_13
elseif batchrun == 14
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_14
elseif batchrun == 15
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_15
elseif batchrun == 16
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_16
elseif batchrun == 17
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_17
elseif batchrun == 18
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_18
elseif batchrun == 19
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_19
elseif batchrun == 20
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_20
elseif batchrun == 21
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_21
elseif batchrun == 22
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_22
elseif batchrun == 23
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_23
elseif batchrun == 24
    cd C:\Users\'HEIMA Lab 2'\Documents\MATLAB\NarcoLogic\model_results\batch_pint\SupplyChain_optint_batch_24
end

save(savefname,'EdgeTable','NodeTable','MOV','FLOW','OUTFLOW','CTRANS',...
    'TOTCPTL','DTOBDGT','slsuccess','activeroute','STOCK','RISKPREM',...
    'slperevent','slval','nactnodes','sltot','t_firstmov','PRICE',...
    'slcount_edges','slcount_vol','slnodes','-v7.3')
% cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM
cd C:\Users\Penelope\Box\NSF_D-ISN\Code\NarcoLogic


end