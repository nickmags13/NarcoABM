%%%%%%%% Master file to run experiments %%%%%%%%%%%%
tic
MRUNS=30;
ERUNS=5;

% sl_min=zeros(1,ERUNS);
% sl_max=zeros(1,ERUNS);
sl_max=[160 240 320 400 480];
sl_min=ceil(sl_max/6);

rng default

for erun=1:ERUNS
    for mrun=1:MRUNS
        
        illicit_supply_chain_v2
        
        filename=sprintf('supplychain_results_091817_%d_%d',erun,mrun);
        cd X:\model_results\SupplyChain_091817
        save(filename,'EdgeTable','NodeTable','MOV','FLOW',...
            'TOTCPTL','DTOBDGT','slsuccess','activeroute','STOCK',...
            'slperevent','nactnodes','sltot','-v7.3')
        cd C:\Users\nmagliocca\Documents\Matlab_code\NarcoLogic
    end
end
toc