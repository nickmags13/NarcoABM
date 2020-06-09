% cd C:\Users\nrmagliocca\'Box Sync'\'Data Drive'\model_results\SupplyChain_full_021618
cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\model_results\SupplyChain_optint_032520
% load supplychain_results_021618_1_6.mat
load supplychain_results_optint_060420_1_1.mat

% [CAadm0,CAattr0]=shaperead('D:\CentralAmerica\GADM\g2015_2014_0\CAadm0.shp',...
%             'UseGeoCoords',true);
% [CAfull,CAfullattr]=shaperead('D:\CentralAmerica\GADM\CA_full_theater_0.shp','UseGeoCoords',...
%             true);
[CAfull,CAfullattr]=shaperead('D:\CentralAmerica\GADM\CA_full_theater_0_simple_01.shp','UseGeoCoords',...
            true);

[CAadm1,CAattr1]=shaperead('D:\CentralAmerica\GADM\g2015_2014_1\CAadm1.shp',...
            'UseGeoCoords',true); 
% lon_west=min(lon_list);
% lon_east=max(lon_list);
% lat_north=max(lat_list);
% lat_south=min(lat_list);
lon_west=min(NodeTable.Lon)-5;
lon_east=max(NodeTable.Lon)+5;
lat_south=min(NodeTable.Lat)-5;
lat_north=max(NodeTable.Lat)+5;
latlimit=[lat_south lat_north];
lonlimit=[lon_west lon_east];
TMAX=180;
NNODES=size(activeroute,1);
nnodes=NNODES;

% cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM
%%% Trafficking movie
writerObj = VideoWriter('intopt_MCI_3_long_fast_2.mp4','MPEG-4');
writerObj.FrameRate=10;
open(writerObj);

hmov=figure;
set(hmov,'color','white')
set(hmov,'Position',[100 10 860 680])
% H=geoshow(CAadm0,'FaceColor',[1 1 1],'EdgeColor',[0.7 0.7 0.7]);
H=geoshow(CAfull,'FaceColor',[1 1 1],'EdgeColor',[0.7 0.7 0.7]);
xlim([-98 -68])
ylim([-5 20])
% worldmap(latlimit,lonlimit);
% load coastlines
% geoshow(coastlat,coastlon);
hold on
idto1=(NodeTable.DTO == 1);
idto2=(NodeTable.DTO == 2);
h_nodes1=plot(NodeTable.Lon(idto1),NodeTable.Lat(idto1),'o','MarkerSize',3,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 0.5 1]);
h_nodes2=plot(NodeTable.Lon(idto2),NodeTable.Lat(idto2),'o','MarkerSize',3,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 0.8 0.5]);
h_producer=plot(NodeTable.Lon(1),NodeTable.Lat(1),'o','MarkerSize',3,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 0.5 0.5]);
h_consumer=plot(NodeTable.Lon(nnodes),NodeTable.Lat(nnodes),'o','MarkerSize',3,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 0 0]);
xlabel('Longitude')
ylabel('Latitude')
% set(H.Parent,'Visible','off')
% h_tempo=plot(NodeTable.Lon(1),NodeTable.Lat(1),'o','MarkerSize',4,'MarkerEdgeColor','red','MarkerFaceColor','red');
% h_tempx=plot(NodeTable.Lon(1),NodeTable.Lat(1),'xk','MarkerSize',20);
% 
% lgd=legend('DTO 1','DTO 2','Producer Node','Receiving Node','Cocaine Shipment','Interdiction Event');
% lgd.FontSize=14;

% for k=1:height(EdgeTable(:,1))
%     h_network(k)=plot([NodeTable.Lon(EdgeTable.EndNodes(k,1)) NodeTable.Lon(EdgeTable.EndNodes(k,2))],...
%                     [NodeTable.Lat(EdgeTable.EndNodes(k,1)) NodeTable.Lat(EdgeTable.EndNodes(k,2))],...
%                     '-','Color',[0.5 0.5 0.5]);
% end
for tt=1:TMAX
    h_title=title(sprintf('Month = %d',tt));
    h_title.FontSize = 14;
    if tt == 1
      startval=OUTFLOW(1,tt+1);
      mrksize=max(ceil(startval/25),4);
      h_start=plot(NodeTable.Lon(1),NodeTable.Lat(1),'o','MarkerSize',mrksize,...
          'MarkerEdgeColor','red','MarkerFaceColor','red');
      frame = getframe(hmov);
      writeVideo(writerObj,frame);
      delete(h_start)
    else
        % plot network at this time step
        [r,c]=ind2sub(size(FLOW(:,:,1)),find(FLOW(:,:,tt)>0));
        for i=1:length(r)
                h_edge(i)=plot([NodeTable.Lon(r(i)) NodeTable.Lon(c(i))],...
                    [NodeTable.Lat(r(i)) NodeTable.Lat(c(i))],'-','Color',[0.5 0.5 0.5]);
%             %%%% Color-code routes
%             if NodeTable.DTO(r(i)) == 1 || NodeTable.DTO(c(i)) == 1
%                 h_edge(i)=plot([NodeTable.Lon(r(i)) NodeTable.Lon(c(i))],...
%                     [NodeTable.Lat(r(i)) NodeTable.Lat(c(i))],'-','Color',[0.5 0.5 1]);
%             elseif NodeTable.DTO(r(i)) == 2 || NodeTable.DTO(c(i)) == 2
%                 h_edge(i)=plot([NodeTable.Lon(r(i)) NodeTable.Lon(c(i))],...
%                     [NodeTable.Lat(r(i)) NodeTable.Lat(c(i))],'-','Color',[0.5 0.8 0.5]);
%             end
        end
        
        % plot interdiction
        if isempty(find(slsuccess(:,:,tt)>0,1)) == 0
            [ir,ic]=ind2sub(size(slsuccess(:,:,1)),find(slsuccess(:,:,tt)>0));
            for ii=1:length(ir)
                h_intedge(ii)=plot([NodeTable.Lon(ir(ii)) NodeTable.Lon(ic(ii))],...
                    [NodeTable.Lat(ir(ii)) NodeTable.Lat(ic(ii))],'-m');
                h_intflow(ii)=plot(NodeTable.Lon(ic(ii)),NodeTable.Lat(ic(ii)),...
                    'xk','MarkerSize',20);
                %                 %%%% Color-code routes
                %                 if NodeTable.DTO(ir(ii)) == 1 || NodeTable.DTO(ic(ii)) == 1
                %                 h_intedge(ii)=plot([NodeTable.Lon(ir(ii)) NodeTable.Lon(ic(ii))],...
                %                     [NodeTable.Lat(ir(ii)) NodeTable.Lat(ic(ii))],'-','Color',[0.5 0.5 1]);
                %                 h_intflow(ii)=plot(NodeTable.Lon(ic(ii)),NodeTable.Lat(ic(ii)),...
                %                     'xk','MarkerSize',20);
                %                 elseif NodeTable.DTO(ir(ii)) == 2 || NodeTable.DTO(ic(ii)) == 2
                %                     h_intedge(ii)=plot([NodeTable.Lon(ir(ii)) NodeTable.Lon(ic(ii))],...
                %                     [NodeTable.Lat(ir(ii)) NodeTable.Lat(ic(ii))],'-','Color',[0.5 0.8 0.5]);
                %                 h_intflow(ii)=plot(NodeTable.Lon(ic(ii)),NodeTable.Lat(ic(ii)),...
                %                     'xk','MarkerSize',20);
                %                 end
                
            end
            frame = getframe(hmov);
            writeVideo(writerObj,frame);
            delete(h_intedge)
            delete(h_intflow)
        end
%         usenders=unique(r);
%         utakers=unique(c);
%         for j=1:length(usenders)
%             islct=find(r == usenders(j));
%             irec=find(r == usenders(utakers == usenders(j)));
%             snodes=r(islct);
%             tnodes=c(islct);
%             for s=1:length(islct)
%                 if j > 1
% %                     delete(h_flow(utakers == r(islct(s))))
%                 delete(h_flow(utakers(irec) == r(islct(s))))
%                 end
%                 flowval=max(ceil(MOV(c(islct(s)),r(islct(s)),tt)/10),4);
%                 h_flow(s)=plot(NodeTable.Lon(c(islct(s))),NodeTable.Lat(c(islct(s))),...
%                     'o','MarkerSize',flowval,'MarkerEdgeColor','red','MarkerFaceColor','red');
%             end
%             frame = getframe(hmov);
%             writeVideo(writerObj,frame);
%             delete(h_flow(s))
%         end
        utakers=find(sum(MOV(:,:,tt),2) > 0);
%         for n=1:nnodes-1
%             if n == 1
                isend=find(MOV(:,1,tt)> 0);
                for s=1:length(isend)
%                     flowval=max(ceil(MOV(isend(s),1,tt)/10),4);
                    flowval=max(ceil(MOV(isend(s),1,tt)/3661),4);
                    h_flow(s)=plot(NodeTable.Lon(isend(s)),NodeTable.Lat(isend(s)),...
                        'o','MarkerSize',flowval,'MarkerEdgeColor','red','MarkerFaceColor','red');
                end
                frame = getframe(hmov);
                writeVideo(writerObj,frame);
                delete(h_flow)
%             else
                for j=1:length(utakers)
                    isend=find(MOV(:,utakers(j),tt) > 0);
                    for s=1:length(isend)
%                     flowval=max(ceil(MOV(isend(s),utakers(j),tt)/10),4);
                    flowval=max(ceil(MOV(isend(s),utakers(j),tt)/3661),4);
                    h_flow(s)=plot(NodeTable.Lon(isend(s)),NodeTable.Lat(isend(s)),...
                        'o','MarkerSize',flowval,'MarkerEdgeColor','red','MarkerFaceColor','red');
                    end
                    frame = getframe(hmov);
                writeVideo(writerObj,frame);
                delete(h_flow)
                end
%             end
%         end

%         for n=1:nnodes-1
% %             isend=find(MOV(:,n,tt)>0);
%             isend=find(FLOW(n,:,tt)>0);
%             if isempty(find(isend>0,1)) == 1
%                 continue
% %             elseif n == nnodes
% %                 flowval=max(ceil(sum(FLOW(:,n,tt))/25),4);
% %                 h_flow=plot(NodeTable.Lon(n),NodeTable.Lat(n),'o','MarkerSize',...
% %                     flowval,'MarkerEdgeColor','red','MarkerFaceColor','red');
% %                 frame = getframe(hmov);
% %                 writeVideo(writerObj,frame);
% %                 delete(h_flow)
%             else
%                 %             lons=[NodeTable.Lon(n*ones(1,length(isend))) NodeTable.Lon(isend)];
%                 %             lats=[NodeTable.Lat(n*ones(1,length(isend))) NodeTable.Lat(isend)];
%                 %             h_edge=plot([NodeTable.Lon(n*ones(1,length(isend))) NodeTable.Lon(isend)],...
%                 %                 [NodeTable.Lat(n*ones(1,length(isend))) NodeTable.Lat(isend)],'-k');
%                 for s=1:length(isend)
%                     %                     if slsuccess(n,isend(s),tt) == 1
%                     %                         h_flow(s)=plot(NodeTable.Lon(isend(s)),NodeTable.Lat(isend(s)),...
%                     %                             'xk','MarkerSize',20);
%                     %                     else
%                     flowval=max(ceil(MOV(isend(s),n,tt)/10),4);
% %                     flowval=max(ceil(FLOW(n,isend(s),tt)/25),4);
%                     h_flow(s)=plot(NodeTable.Lon(isend(s)),NodeTable.Lat(isend(s)),...
%                         'o','MarkerSize',flowval,'MarkerEdgeColor','red','MarkerFaceColor','red');
%                     %                     h_flow=plot(NodeTable.Lon(n*ones(length(isend),1)),...
%                     %                         NodeTable.Lat(isend'),'o','MarkerSize',flowval,...
%                     %                         'MarkerEdgeColor','k','MarkerFaceColor','red');
%                     %                     end
%                 end
%                 keyboard
%                 frame = getframe(hmov);
%                 writeVideo(writerObj,frame);
%                 delete(h_flow)
%             end
%         end
        delete(h_edge)
%         frame = getframe(hmov);
%         writeVideo(writerObj,frame);
%         if exist('h_edge') == 1
%             delete(h_edge)
%         end
%         if exist('h_intedge') == 1
%             delete(h_intedge)
%         end
%         if exist('h_intflow') == 1
%             delete(h_intflow)
%         end
    end
end
close(writerObj);