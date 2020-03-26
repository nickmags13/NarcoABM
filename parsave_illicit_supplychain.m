function parsave_illicit_supplychain(savefname,EdgeTable,NodeTable,MOV,...
    FLOW,OUTFLOW,CTRANS,TOTCPTL,DTOBDGT,slsuccess,activeroute,STOCK,RISKPREM,...
    slperevent,slval,nactnodes,sltot,t_firstmov,PRICE)

% cd C:\Users\nrmagliocca\'Box Sync'\'Data Drive'\model_results\SupplyChain_full_021618
cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\model_results\SupplyChain_optint_032520

save(savefname,'EdgeTable','NodeTable','MOV','FLOW','OUTFLOW','CTRANS',...
    'TOTCPTL','DTOBDGT','slsuccess','activeroute','STOCK','RISKPREM',...
    'slperevent','slval','nactnodes','sltot','t_firstmov','PRICE','-v7.3')
cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM


end