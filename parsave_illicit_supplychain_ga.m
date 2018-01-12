function parsave_illicit_supplychain_ga(savefname,t_firstmov,deptflows_ts,cntryflows_ts)
tfirstmov=t_firstmov;
deptflowsts=deptflows_ts;
cntryflowsts=cntryflows_ts;
save(savefname,'tfirstmov','deptflowsts','cntryflowsts','-v7.3')

end