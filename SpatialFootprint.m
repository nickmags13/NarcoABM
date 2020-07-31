% cd C:\Users\nrmagliocca\'Box Sync'\'Data Drive'\model_results\SupplyChain_full_021618
cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\model_results\SupplyChain_optint_072020
% load supplychain_results_021618_1_6.mat
load supplychain_results_optint_072020_1_1.mat

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
writerObj = VideoWriter('intopt_MCI_3_long_fast_180.mp4','MPEG-4');
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
idto1=(NodeTable.DTO == 1); %Pacific
idto2=(NodeTable.DTO == 2); %Caribbean
% h_nodes1=plot(NodeTable.Lon(idto1),NodeTable.Lat(idto1),'o','MarkerSize',3,'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[0.5 0.5 1]);
% h_nodes2=plot(NodeTable.Lon(idto2),NodeTable.Lat(idto2),'o','MarkerSize',3,'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[0.5 0.8 0.5]);
h_nodes1=plot(NodeTable.Lon(idto1),NodeTable.Lat(idto1),'o','MarkerSize',3,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 0.5 0.5]);
h_nodes2=plot(NodeTable.Lon(idto2),NodeTable.Lat(idto2),'o','MarkerSize',3,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 0.5 0.5]);
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
% for tt=1:60
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
        h_slnodes=plot(NodeTable.Lon(slnodes{tt}),NodeTable.Lat(slnodes{tt}),...
            '^','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','b');
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
            %             delete(h_intedge)
            %             delete(h_intflow)
        end
        utakers=find(sum(MOV(:,:,tt),2) > 0);
        isend=find(MOV(:,1,tt)> 0);
        for s=1:length(isend)
%             flowval=max(ceil(MOV(isend(s),1,tt)/100),4);
            flowval=max(ceil(MOV(isend(s),1,tt)/10000),4);
            h_flow(s)=plot(NodeTable.Lon(isend(s)),NodeTable.Lat(isend(s)),...
                'o','MarkerSize',flowval,'MarkerEdgeColor','red','MarkerFaceColor','red');
        end
%         legend([h_nodes1(1) h_edge(1) h_slnodes(1) h_intedge(1) h_intflow(1) h_flow(1)],...
%         'Trafficking Network Nodes','Trafficking Network Links',...
%             'Location of Interdiction Assets','Interdicted Network Link',...
%             'Interdicted Network Node','Cocaine Shipments',...
%             'Location','southeast')
        frame = getframe(hmov);
        writeVideo(writerObj,frame);
        delete(h_flow)
        for j=1:length(utakers)
            isend=find(MOV(:,utakers(j),tt) > 0);
            for s=1:length(isend)
%                 flowval=max(ceil(MOV(isend(s),utakers(j),tt)/100),4);
                flowval=max(ceil(MOV(isend(s),utakers(j),tt)/10000),4);
                h_flow(s)=plot(NodeTable.Lon(isend(s)),NodeTable.Lat(isend(s)),...
                    'o','MarkerSize',flowval,'MarkerEdgeColor','red','MarkerFaceColor','red');
            end
            frame = getframe(hmov);
            writeVideo(writerObj,frame);
            delete(h_flow)
        end
        delete(h_edge)
        delete(h_slnodes)
%         delete(h_intedge)
%         delete(h_intflow)
        %         frame = getframe(hmov);
        %         writeVideo(writerObj,frame);
        %         if exist('h_edge') == 1
        %             delete(h_edge)
        %         end
        if exist('h_intedge','var') == 1
            delete(h_intedge)
        end
        if exist('h_intflow','var') == 1
            delete(h_intflow)
        end
    end
end
close(writerObj);