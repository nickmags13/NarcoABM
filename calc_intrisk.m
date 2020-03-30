function [sl_risk,slevnt,tmevnt]=calc_intrisk(sloccur,...
    t_eff,alpharisk,betarisk,timeweight)

slevnt=sum(sloccur.*repmat((timeweight.^(t_eff))',1,length(sloccur(1,:))),1);
% intrdevnt=sum(intrdoccur.*repmat(timeweight.^(t-(t:-1:max(TSTART+1,t-12))),length(intrdoccur(:,1)),1),2);
% intrdevnt=sum(intrdoccur.*repmat((timeweight.^(t_eff))',1,length(intrdoccur(1,:))),1);
% tmevnt=sum(timeweight.^(t-(t:-1:max(TSTART+1,t-12))));
tmevnt=sum(timeweight.^(t_eff));
sl_risk=(slevnt+alpharisk)./(tmevnt+alpharisk+betarisk);
% intrd_risk=(intrdevnt+alpharisk)./(tmevnt+alpharisk+betarisk);
end
