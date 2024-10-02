function gccp_binmap(cbo,eventdata,handles)
%GCCP_BINMAP Map of GCCP nodes
%   GCCP_BINMAP is a callback function that creates a regional map of the
%   GCCP nodes and the Gaussian radius for stacking.

%   Author: Kevin C. Eagar
%   Date Created: 04/24/2011
%   Last Updated: 06/15/2011
%   Edit: Rob Porritt
%   Last Updated: 04/06/2015
%   Reason: Remove dependence on mapping tool box

% Get the handles of the FuncLab GUI
%--------------------------------------------------------------------------
h = guidata(gcbf);

MapDataPath = regexprep(which('funclab'),'funclab.m','map_data/');
ProjectDirectory = evalin('base','ProjectDirectory');
ProjectFile = evalin('base','ProjectFile');
load([ProjectDirectory ProjectFile],'FuncLabPreferences')
SetManagementCell = evalin('base','SetManagementCell');
CurrentSubsetIndex = evalin('base','CurrentSubsetIndex');
% Sets the name of the set. Removes white spaces and capitalizes between
% words
if CurrentSubsetIndex ~= 1
    SetName = SetManagementCell{CurrentSubsetIndex,5};
    % Remove spaces because they make filesystems go a little haywire
    % sometimes
    SpaceIDX = isspace(SetName);
    for idx=2:length(SetName)
        if isspace(SetName(idx-1))
            SetName(idx) = upper(SetName(idx));
        end
    end
    SetName = SetName(~SpaceIDX);
else
    SetName = 'Base';
end
if CurrentSubsetIndex == 1
    GCCPMatFile = fullfile(ProjectDirectory, 'GCCPData', 'GCCPData.mat');
else
    GCCPMatFile = fullfile(ProjectDirectory, 'GCCPData', SetName, 'GCCPData.mat');
end

% Set initial values
%--------------------------------------------------------------------------
load(GCCPMatFile,'GCCPGrid');
LatMax = max(cat(1,GCCPGrid{:,1}));  % 06/15/2011: CHANGED COLUMN IN GCCPGrid FROM 2 TO 1
LatMin = min(cat(1,GCCPGrid{:,1}));  % 06/15/2011: CHANGED COLUMN IN GCCPGrid FROM 2 TO 1
LonMax = max(cat(1,GCCPGrid{:,2}));  % 06/15/2011: CHANGED COLUMN IN GCCPGrid FROM 3 TO 2
LonMin = min(cat(1,GCCPGrid{:,2}));  % 06/15/2011: CHANGED COLUMN IN GCCPGrid FROM 3 TO 2
LatRange = abs(LatMax-LatMin);
LonRange = abs(LonMax-LonMin);
NorthLat = LatMax+0.1*LatRange;
SouthLat = LatMin-0.1*LatRange;
EastLon = LonMax+0.1*LonRange;
WestLon = LonMin-0.1*LonRange;

% Set map projection
%--------------------------------------------------------------------------
MapProjection = FuncLabPreferences{10};

% Setup the regional map using low-level functions (geoshow is giving me
% problems for some reason - keeps crashing and logging me out!)
%--------------------------------------------------------------------------
figure('Units','characters','NumberTitle','off','Name','GCCP Node Map','Units','characters',...
    'Position',[0 0 165 55],'Color',[1 1 1],'CreateFcn',{@movegui_kce,'center'});
% This part is where the mapping toolbox is needed. RWP adding switch to
% allow for crude map
% Checks for existance of mapping toolkit
if license('test', 'MAP_Toolbox')

    axesm('MapProjection',MapProjection,'MapLatLimit',[SouthLat NorthLat],'MapLonLimit',[WestLon EastLon],...
        'Frame','on','Grid','on','MeridianLabel','on','MLineLocation',1,'ParallelLabel','on','PLineLocation',1);
    set(gca,'XColor',[1 1 1],'YColor',[1 1 1]);
    get(gca,'position');
    switch FuncLabPreferences{12}
        case 'Coastlines (low-res)'
            % Load in coastline data
            %------------------------------------------------------------------
            coasts = gshhs([MapDataPath 'gshhs_c.b'],[SouthLat NorthLat],[WestLon EastLon]);
            geoshow(coasts, 'FaceColor',[.6 .4 .3],'EdgeColor',[0 0 0])
            linem([coasts.Lat],[coasts.Lon],'LineWidth',1,'Color',[0 0 0])
        case 'Coastlines (medium-res)'
            % Load in coastline data
            %------------------------------------------------------------------
            coasts = gshhs([MapDataPath 'gshhs_i.b'],[SouthLat NorthLat],[WestLon EastLon]);
            geoshow(coasts, 'FaceColor',[.6 .4 .3],'EdgeColor',[0 0 0])
            linem([coasts.Lat],[coasts.Lon],'LineWidth',1,'Color',[0 0 0])
        case 'Coastlines (hi-res)'
            % Load in coastline data
            %------------------------------------------------------------------
            coasts = gshhs([MapDataPath 'gshhs_f.b'],[SouthLat NorthLat],[WestLon EastLon]);
            geoshow(coasts, 'FaceColor',[.6 .4 .3],'EdgeColor',[0 0 0])
            linem([coasts.Lat],[coasts.Lon],'LineWidth',1,'Color',[0 0 0])
        case 'Topography (low-res)'
            % Load in ETOPO2v2 Topography and Bathymatry data and plot
            %------------------------------------------------------------------
            try
                [TopoGrid,topo_refvec] = etopo([MapDataPath 'ETOPO2v2c_i2_MSB.bin'],10,[SouthLat NorthLat],[WestLon EastLon]); % Updated ETOPO
                disp('GCCP Node map made using updated ETOPO2v2c topography.')
            catch
                [TopoGrid,topo_refvec] = etopo([MapDataPath 'ETOPO2.raw.bin'],10,[SouthLat NorthLat],[WestLon EastLon]); % Old ETOPO w/ 1-cell error
                disp('GCCP Node map made using older ETOPO2-2001 topography.')
            end
            [topo_lat,topo_lon] = meshgrat([SouthLat NorthLat],[WestLon EastLon],size(TopoGrid));
            surfacem(topo_lat,topo_lon,TopoGrid);
            demcmap(TopoGrid)
        case 'Topography (medium-res)'
            [TopoGrid,topo_refvec] = etopo([MapDataPath 'etopo1_bed_c.flt'],10,[SouthLat NorthLat],[WestLon EastLon]);
            [topo_lat,topo_lon] = meshgrat([SouthLat NorthLat],[WestLon EastLon],size(TopoGrid));
            surfacem(topo_lat,topo_lon,TopoGrid);
            demcmap(TopoGrid)
        case 'Topography (hi-res)'
            globeFile = globedems([SouthLat NorthLat],[WestLon EastLon]);
            for n = 1:length(globeFile)
                [TopoGrid,topo_refvec] = globedem([MapDataPath '/GLOBE_DEM/' globeFile{n}],10,[SouthLat NorthLat],[WestLon EastLon]);
                [topo_lat,topo_lon] = meshgrat([SouthLat NorthLat],[WestLon EastLon],size(TopoGrid));
                surfacem(topo_lat,topo_lon,TopoGrid);
                demcmap(TopoGrid)
            end
    end
    NodeColor = [.3 .3 .3];
    linem(cat(1,GCCPGrid{:,1}),cat(1,GCCPGrid{:,2}),'LineStyle','none','Marker','o','MarkerSize',4,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',NodeColor)   % 06/15/2011: CHANGED COLUMNS IN GCCPGrid
    [lat,lon] = scircle1(str2double(FuncLabPreferences{17}),str2double(FuncLabPreferences{18}),km2deg(GCCPGrid{1,3}));  % 06/15/2011: CHANGED COLUMNS IN GCCPGrid
    linem(lat,lon,'Color',[1 0 0],'LineWidth',2)
else
    if strcmp(MapProjectionID,'Stereographic') || strcmp(MapProjectionID,'Orthographic') || strcmp(MapProjectionID,'Azimuthal Equal-area') || ...
        strcmp(MapProjectionID,'Azimuthal Equidistant') || strcmp(MapProjectionID,'Gnomonic') || strcmp(MapProjectionID,'Satellite')    
        [~, ~, xdeg] = distaz(SouthLat, WestLon, NorthLat, EastLon);
        m_proj(MapProjectionID,'latitude',(SouthLat + NorthLat)/2, 'longitude',(WestLon + EastLon) / 2,'radius',xdeg/2);
    else
        m_proj(MapProjectionID,'latitude',[SouthLat NorthLat],'longitude',[WestLon EastLon]);
    end
    set(gca,'XColor',[1 1 1],'YColor',[1 1 1]);
    get(gca,'position');
    switch FuncLabPreferences{12}
        case 'Coastlines (low-res)'
            if ~exist([MapDataPath 'gshhs_l.b']) && ~exist('gshhs_c.b')
                disp('No gshhs_c.b coastline data exists.  Please download at http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/version2.1.1/gshhs_2.1.1.zip')
                disp('m_map version does not support variable resolution coast lines or topography.')
                m_gshhs_l('patch',[.9 .9 .9],'edgecolor','black');
                m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',12,'box','on');
            else
                m_gshhs_l('patch',[.9 .9 .9],'edgecolor','black');
                m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',12,'box','on');
            end
        case 'Coastlines (medium-res)'
            disp('m_map version does not fully support variable resolution coast lines or topography.')
            %m_coast('patch',[.9 .9 .9],'edgecolor','black');
            m_gshhs_i('patch',[.9 .9 .9],'edgecolor','black');
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',12,'box','on');
        case 'Coastlines (hi-res)'
            disp('m_map version does not support variable resolution coast lines or topography.')
            %m_coast('patch',[.9 .9 .9],'edgecolor','black');
            m_gshhs_h('patch',[.9 .9 .9],'edgecolor','black');
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',12,'box','on');
        case 'Topography (low-res)'
            disp('m_map version does not fully support variable resolution coast lines or topography.')
            m_gshhs_l('patch',[.9 .9 .9],'edgecolor','black');
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',12,'box','on');
            cptcmap('topography_rwp.cpt')
            m_elev;
        case 'Topography (medium-res)'
            disp('m_map version does not fully support variable resolution coast lines or topography.')
            %m_coast('patch',[.9 .9 .9],'edgecolor','black');
            m_gshhs_i('patch',[.9 .9 .9],'edgecolor','black');
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',12,'box','on');
            cptcmap('topography_rwp.cpt')
            m_elev;
        case 'Topography (hi-res)'
            disp('m_map version does not fully support variable resolution coast lines or topography.')
            %m_coast('patch',[.9 .9 .9],'edgecolor','black');
            m_gshhs_h('patch',[.9 .9 .9],'edgecolor','black');
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',12,'box','on');
            cptcmap('topography_rwp.cpt')
            m_elev;
    end
    
    NodeColor = [.3 .3 .3];
%     m_line(cat(1,GCCPGrid{:,2}),cat(1,GCCPGrid{:,1}),'LineStyle','none','Marker','o','MarkerSize',4,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',NodeColor) % 06/12/2011: CHANGED COLUMNS IN CCPGrid
%     for n=1:length(GCCPGrid(:,1))    % 06/12/2011: CHANGED COLUMNS IN CCPGrid
%         m_range_ring(GCCPGrid{n,2},GCCPGrid{n,1},GCCPGrid{n,3},10);
%         %[lat,lon] = scircle1(CCPGrid{n,1},CCPGrid{n,2},km2deg(CCPGrid{n,3}));   % 06/12/2011: CHANGED COLUMNS IN CCPGrid
%         %linem(lat,lon,'Color',[1 0 0])
%     end

    
    
    
    
%    % Version without mapping toolbox
%    disp('Making rough map due to lack of mapping toolbox')
%    load FL_coasts
%    plot(ncst(:,1),ncst(:,2),'k-')
%    axis([WestLon EastLon SouthLat NorthLat])
%    hold on
%    NodeColor = [.3 .3 .3];
%    binpts = [cat(1,GCCPGrid{:,1}),cat(1,GCCPGrid{:,2})];
%    plot(binpts(:,2), binpts(:,1),'LineStyle','none','Marker','o','MarkerSize',4,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',NodeColor)    
end