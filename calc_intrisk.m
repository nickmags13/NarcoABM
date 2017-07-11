function [sl_risk,intrd_risk,slevnt,intrdevnt,tmevnt]=calc_intrisk(sloccur,...
    intrdoccur,t,TSTART,alpharisk,betarisk,timeweight)

% slevnt=sum(sloccur.*repmat((timeweight.^(t-(t:-1:TSTART+1)))',1,length(sloccur(1,:))),1);
% intrdevnt=sum(intrdoccur.*repmat(timeweight.^(t-(t:-1:TSTART+1)),length(intrdoccur(:,1)),1),2);
% tmevnt=sum(timeweight.^(t-(t:-1:TSTART+1)));
slevnt=sum(sloccur.*repmat((timeweight.^(t-(t:-1:max(TSTART+1,t-12))))',1,length(sloccur(1,:))),1);
intrdevnt=sum(intrdoccur.*repmat(timeweight.^(t-(t:-1:max(TSTART+1,t-12))),length(intrdoccur(:,1)),1),2);
tmevnt=sum(timeweight.^(t-(t:-1:max(TSTART+1,t-12))));
sl_risk=(slevnt+alpharisk)./(tmevnt+alpharisk+betarisk);
intrd_risk=(intrdevnt+alpharisk)./(tmevnt+alpharisk+betarisk);
end
