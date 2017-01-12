function [sl_risk,intrd_risk,slevnt,intrdevnt,tmevnt]=calc_intrisk(sloccur,...
    intrdoccur,t,TSTART,alpharisk,betarisk,timeweight)

% slevnt=sum(sum(sloccur,1).*timeweight.^(t-1:-1:TSTART));
% intrdevnt=sum(sum(intrdoccur,2)'.*timeweight.^(t-1:-1:TSTART));
slevnt=sum(sloccur.*repmat((timeweight.^(t-1:-1:TSTART))',1,length(sloccur(1,:))),1);
intrdevnt=sum(intrdoccur.*repmat(timeweight.^(t-1:-1:TSTART),length(intrdoccur(:,1)),1),2);
tmevnt=sum(timeweight.^(t-1:-1:TSTART));
sl_risk=(slevnt+alpharisk)./(tmevnt+alpharisk+betarisk);
intrd_risk=(intrdevnt+alpharisk)./(tmevnt+alpharisk+betarisk);
end
