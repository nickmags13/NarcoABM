function [sl_risk,intrd_risk,slevnt,intrdevnt,tmevnt]=calc_intrisk(sloccur,...
    intrdoccur,t,TSTART,alpharisk,betarisk,timeweight)

% slevnt=sum(sum(sloccur,1).*timeweight.^(t-1:-1:TSTART));
% intrdevnt=sum(sum(intrdoccur,2)'.*timeweight.^(t-1:-1:TSTART));
slevnt=sum(sloccur,1).*timeweight.^(t-1:-1:TSTART);
intrdevnt=sum(intrdoccur,2)'.*timeweight.^(t-1:-1:TSTART);
tmevnt=sum(timeweight.^(t-1:-1:TSTART));
sl_risk=(slevnt+alpharisk)./(tmevnt+alpharisk+betarisk);
intrd_risk=(intrdevnt+alpharisk)./(tmevnt+alpharisk+betarisk);
end
