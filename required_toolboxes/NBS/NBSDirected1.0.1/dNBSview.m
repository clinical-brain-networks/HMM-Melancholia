function dNBSview(varargin)
%dNBSview  View networks in an NBS structure.
%   dNBSview will launch the viewer and prompt for an NBS structure to be
%   loaded.
%
%   dNBSview(NBS) will display the networks in NBS.mat.
%
%   An NBS structure contains the following fields:
%       NBS.n:            Number of networks in the NBS structure
%       NBS.con_mat{i}:   N x N upper-triangular, binary adjacency matrix
%                         specifying the ith network, i=1:n
%       NBS.node_coor:    N x 3 array of MNI node coordinates
%       NBS.node_label:   N x 1 cell array of node labels [optional]
%       NBS.test_stat:    N x N array of test-statistic values for each 
%                         edge [optional]
%       NBS.pval(i):      p-value for the ith network, i=1:n [optional]
%
%   Example NBS structure containing one network:
%       NBS.n=1;
%       NBS.con_mat=triu([0 1 0 0 0 
%                         1 0 0 1 1
%                         0 0 0 1 0
%                         0 1 1 0 0 
%                         0 1 0 0 0],1);
%       NBS.node_coor=[-20 31 43
%                      -32  52 -11
%                       40  -6  50
%                      -26   0 -18
%                       16 -68  -5];
%       NBS.node_label={'A';'B';'C';'D';'E'};
%       NBS.test_stat=[  0 2.3 3.7 0.1 3.0
%                      2.3   0 1.8 3.2 2.2
%                      3.7 1.8   0 2.6 2.6
%                      0.1 3.2 2.6   0 2.5
%                      3.0 2.2 2.6 2.5   0];
%       NBS.pval=0.02; 
%
%   azalesky@unimelb.edu.au



%Dependencies
%   MIP.mat: SPM background image used for maximum intensity projection
%   glicons.mat: Icons for toolbar taken from graphViz4Matlab
%
%Main Structures: 
%   S: structure containing most object handles
%   Net: structure describing display properties
%   Coor: structure specifying node coordinates
%
%Default network display properties are specified in function draw_net


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up background image 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPM maximum intensity projection in MNI space (mip96)
load(['icons',filesep,'MIP']);
img=mask_all+grid_all; %img=mip96+mask_all;
%Eliminate grid lines outside brain
ind_grid_outside_brain=find(img==2); 
img=img+mask_all;
img(ind_grid_outside_brain)=1; 
%Set background 
tmp=ones(size(img));
for i=26:198
    for j=22:157
        tmp(i,j)=0;
    end
    for j=212:336
        tmp(i,j)=0;
    end
end
for i=240:375
    for j=212:336
        tmp(i,j)=0; 
    end
end
%Set background color
%*0 = white
%*3 = black
img=img+tmp*0; 

%Colormap
C=[1 1 1;         %white: brain 
   0.7 0.7 0.7    %gray: grid lines in brain and non-brain
   0.4 0.4 0.4];  %dark gray: outside brain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up figure and figure sizing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Image aspect
Aspect=size(img);

%Screen aspect
sz=get(0,'ScreenSize');
ratio=sz(4)/sz(3); 

S.fh = figure('units','normalized',...
              'position',[0.055 0.055 1-0.2 1-0.2],...
              'menubar','none',...
              'name','dNBSview',...              
              'numbertitle','off',...
              'color','white',...
              'resize','on',...
              'colormap',C); 
          
offset=0; 
cntr_x=0.5;
cntr_y=0.5;

max_x=1-2*offset;
max_y=1-2*offset; 
width=Aspect(1);
height=Aspect(2); 

if height<width
    orig_x=cntr_x-max_x/2*ratio;
    orig_y=cntr_y-max_y/2*(height/width)*(max_x/max_y);  
    len_x=max_x*ratio; 
    len_y=max_y*(height/width)*(max_x/max_y);
else
    orig_y=cntr_y-max_y/2;
    orig_x=cntr_x-max_x/2*(width/height)*(max_y/max_x)*ratio; 
    len_x=max_x*(width/height)*(max_y/max_x)*ratio;
    len_y=max_y;
end
P.mip=[orig_x orig_y len_x len_y]; 

%Determined from mip96
P.list=[240 Aspect(2)-343]; 
P.list(1)=orig_x+P.list(1)/width*len_x; 
P.list(2)=orig_y+P.list(2)/height*len_y; 
x_end=orig_x+375/width*len_x;
y_end=orig_y+(Aspect(2)-208)/height*len_y;
P.list(3)=x_end-P.list(1);
P.list(4)=y_end-P.list(2);   

S.ax = axes('units','normalized',...
            'position',P.mip,...
            'box','off');
        
S.ls = uicontrol('style','list',...
                 'unit','normalized',...
                 'position',P.list,...
                 'min',0,'max',2,...
                 'fontsize',12,...
                 'BackgroundColor','white',...
                 'Value',[],...
                 'String','Load network...'); 



           
             
             
%Set up toolbar buttons
S.toolbar=uitoolbar(S.fh);            
load(['icons',filesep,'glicons']);
S.button(1)=uipushtool(S.toolbar,...
            'TooltipString','Decrease Node Size',...
            'CData',icons.downblue);
S.button(2)=uipushtool(S.toolbar,...
           'TooltipString','Increase Node Size',...
           'CData',icons.upblue);
S.button(3)=uipushtool(S.toolbar,...
           'TooltipString','Decrease Line Width',...
           'CData',icons.downdarkblue);
S.button(4)=uipushtool(S.toolbar,...
           'TooltipString','Increase Line Width',...
           'CData',icons.updarkblue);
S.button(5)=uipushtool(S.toolbar,...
           'TooltipString','Reset Sizing',...
           'CData',icons.flip);  
S.button(6)=uitoggletool(S.toolbar,...
           'TooltipString','Show All Nodes',...
           'CData',icons.tree);  

S.IH_sag = imagesc(rot90(img),'Parent',S.ax);  % Display the image.
axis image;
set(S.ax,'xtick',[],'ytick',[],'box','off','Visible','off');
hold on;

%Remove highlighting from listbox
set(S.ls,'Callback',{@bdfcn_ls,S});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display first network if NBS structure exists, otherwise open the figure
%without a network and prompt for NBS structure to be loaded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0
    S.menu_load=uimenu(S.fh,'Label','File');
    S.menu_load_sub(1)=uimenu(S.menu_load,'Label','Load',...
                              'Callback',{@load_nbs,S});
    S.menu_load_sub(2)=uimenu(S.menu_load,'Label','Close',...
                              'Callback',{@close_fig,S});
else 
    %Draw first network
    NBS=varargin{1};
    if isfield(NBS,'NBS')
        NBS=NBS.NBS; 
    end
    n=1;
    try draw_net(S,NBS,n);
        if NBS.n>1
            for n1=2:NBS.n
                draw_net(S,NBS,n1);
                
            end
        end
    catch; msg='Does not contain valid nbs structure';
           tmp=get(S.ls,'string');
           set(S.ls,'string',[{msg};tmp]);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load NBS sructure and draw network 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function load_nbs(varargin)
    S=varargin{3};
    [filename,pathname]=uigetfile('*.mat','Select the Matlab (*.mat) file');
    if ~isequal(filename,0)
        ok=0;
        try h=load([pathname,filename]); msg=sprintf('%s',filename); ok=1;
        catch;  msg=sprintf('%s could not be opened',filename); end
        if ok
            n=1;
            try NBS=h.nbs.NBS; draw_net(S,NBS,n); ok=1;
            catch; ok=0; end
            if ~ok
                try NBS=h.NBS; draw_net(S,NBS,n); ok=1;
                catch; ok=0; end
            end 
            if ~ok
                msg=sprintf('%s does not contain valid nbs structure',filename);
            end
        end
        if ok
            tmp=get(S.ls,'string');
            set(S.ls,'string',[{msg};tmp]);
        else
            tmp=get(S.ls,'string');
            set(S.ls,'string',[{msg};tmp]);
            %pause(5);
            set(S.ls,'string',tmp);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Close figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function close_fig(varargin)
    S=varargin{3};
    delete(S.fh); 
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adjust network display properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decrease_node_size(varargin)
    S=varargin{3}; 
    if isfield(S,'mark_axl')
       nw=get(S.mark_axl(1),'MarkerSize');
       if nw>3
            nw=nw-3;
            set(S.mark_axl(:),'MarkerSize',nw);
            set(S.mark_sag(:),'MarkerSize',nw);
            set(S.mark_cor(:),'MarkerSize',nw);
       end
    end
end
function increase_node_size(varargin)
    S=varargin{3}; 
    if isfield(S,'mark_axl')
       nw=get(S.mark_axl(1),'MarkerSize');
       nw=nw+3;
       set(S.mark_axl(:),'MarkerSize',nw);
       set(S.mark_sag(:),'MarkerSize',nw);
       set(S.mark_cor(:),'MarkerSize',nw);
    end
end
function decrease_line_width(varargin)
    S=varargin{3}; 
    if isfield(S,'line_axl')
       lw=get(S.line_axl(1),'LineWidth');
       if lw>0.5
            lw=lw-0.5;
            set(S.line_axl(:),'LineWidth',lw);
            set(S.line_sag(:),'LineWidth',lw);
            set(S.line_cor(:),'LineWidth',lw);
       end
    end
end
function increase_line_width(varargin)
    S=varargin{3}; 
    if isfield(S,'line_axl')
       lw=get(S.line_axl(1),'LineWidth');
            lw=lw+0.5;
            set(S.line_axl(:),'LineWidth',lw);
            set(S.line_sag(:),'LineWidth',lw);
            set(S.line_cor(:),'LineWidth',lw);
    end
end
function reset(varargin)
    S=varargin{3}; 
    Net=varargin{4}; 
    if isfield(S,'line_axl')
        set(S.line_axl(:),'LineWidth',Net.lw);
        set(S.line_sag(:),'LineWidth',Net.lw);
        set(S.line_cor(:),'LineWidth',Net.lw);
        set(S.mark_axl(:),'Color',Net.node_color);
        set(S.mark_sag(:),'Color',Net.node_color);
        set(S.mark_cor(:),'Color',Net.node_color);
    end   
    if isfield(S,'mark_axl')
       set(S.mark_axl(:),'MarkerSize',Net.nw);
       set(S.mark_sag(:),'MarkerSize',Net.nw);
       set(S.mark_cor(:),'MarkerSize',Net.nw);
       set(S.mark_axl(:),'Color',Net.node_color);
       set(S.mark_sag(:),'Color',Net.node_color);
       set(S.mark_cor(:),'Color',Net.node_color);
    end
end
function show_all_nodes(varargin)
    S=varargin{3};  
    if isfield(S,'mark_axl')
        for i=1:length(S.mark_axl)
            if strcmp('off',get(S.mark_axl(i),'Visible'))
                set(S.mark_axl(i),'Visible','on','UserData',1);
                set(S.mark_sag(i),'Visible','on','UserData',1);
                set(S.mark_cor(i),'Visible','on','UserData',1);
            end
        end
    end
end
function undo_show_all_nodes(varargin)
    S=varargin{3};  
    if isfield(S,'mark_axl')
        for i=1:length(S.mark_axl)
            if get(S.mark_axl(i),'UserData')
                set(S.mark_axl(i),'Visible','off','UserData',0);
                set(S.mark_sag(i),'Visible','off','UserData',0);
                set(S.mark_cor(i),'Visible','off','UserData',0);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Draw network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function draw_net(varargin)
     if length(varargin)>4
         S=varargin{3};
         NBS=varargin{4};
         n=varargin{5};
     else
         S=varargin{1};
         NBS=varargin{2};
         n=varargin{3};
     end
     %Delete any previous dropdown menu
     delete(get(S.ax,'UserData'));
     set(S.ax,'UserData',[]);
     
     if ~isempty(NBS.n) && NBS.n>0; 
         
         %Only create a dropdown menu if more than one signficant network
         %and dropdown menu not previously created
         if NBS.n>1 %&& isempty(get(S.ax,'UserData'))
             S.menu_view=uimenu(S.fh,'Label','View');
             for i=1:NBS.n
                 S.menu_view_net(i)=uimenu(S.menu_view,'Label',['Network ',num2str(i)],...
                     'Callback',{@draw_net,S,NBS,i});
                 
                 
             end
             %Use UserData of the axes object to indicate whether a menu has been
             %previously set up
             set(S.ax,'UserData',S.menu_view);
         end
         
         %Origins in MIP
         origin.axl=[126 275];
         origin.sag=[126 111];
         origin.cor=[308 111];
         for i=1:size(NBS.node_coor,1)
             Coor.x_axl(i)=origin.axl(1)+NBS.node_coor(i,2);
             Coor.y_axl(i)=origin.axl(2)+NBS.node_coor(i,1);
             Coor.x_sag(i)=origin.sag(1)+NBS.node_coor(i,2);
             Coor.y_sag(i)=origin.sag(2)-NBS.node_coor(i,3);
             Coor.x_cor(i)=origin.cor(1)+NBS.node_coor(i,1);
             Coor.y_cor(i)=origin.cor(2)-NBS.node_coor(i,3);
         end
         
         %Determine indices of edges
         if iscell(NBS.con_mat)
             [ind_i,ind_j]=find(NBS.con_mat{n});
         else
             [ind_i,ind_j]=find(NBS.con_mat);
         end
         %Line width
         lw=2;
         %Default line color
         line_color=[0 0 1];
         %Node width
         nw=30;
         %Default node color
         node_color=[0 0 0];
         %Make a structure
         Net.lw=lw; Net.line_color=line_color; Net.nw=nw; Net.node_color=node_color;
         if iscell(NBS.con_mat)
             Net.clicked_node=spalloc(length(NBS.con_mat{n}),1,1);
             Net.clicked_line=spalloc(length(NBS.con_mat{n}),length(NBS.con_mat{n}),1);
             Net.N=length(NBS.con_mat{n});
         else
             Net.clicked_node=spalloc(length(NBS.con_mat),1,1);
             Net.clicked_line=spalloc(length(NBS.con_mat),length(NBS.con_mat),1);
             Net.N=length(NBS.con_mat);
         end
         if isfield(NBS,'test_stat')
             Net.test_stat=NBS.test_stat;
         end
         Net.ind_i=ind_i;
         Net.ind_j=ind_j;
         
         
         
         
         
         %Delete previous network
         prev_handles=get(S.fh,'UserData');
         if ~isempty(prev_handles)
             %Delete handles
             delete(prev_handles(:));
         end
         
         %Reset toggle to off
         set(S.button(6),'OnCallback','');
         set(S.button(6),'OffCallback','');
         set(S.button(6),'State','off');
         
         %Put information in list box
         if isfield(NBS,'pval')
             pval=NBS.pval(n);
             if pval<0.001
                 s_pval='<0.001';
             else
                 s_pval=sprintf('%0.3f',pval);
             end
             s={['Displaying network ',num2str(n),' of ',num2str(NBS.n)];...
                 ['p-value: ',s_pval];...
                 ['Number of edges: ',num2str(length(Net.ind_i))];...
                 ['Number of nodes: ',num2str(length(unique([Net.ind_i;Net.ind_j])))];...
                 'Click on a node or edge...'};
         else
             s={['Displaying network ',num2str(n),' of ',num2str(length(NBS.n))];...
                 ['Number of edges: ',num2str(length(Net.ind_i))];...
                 ['Number of nodes: ',num2str(length(unique([Net.ind_i;Net.ind_j])))];...
                 'Click on a node or edge...'};
         end
         set(S.ls,'String',s);
         
         %Draw edges
         [S.line_axl,S.line_sag,S.line_cor]=draw_edges(Coor,Net);
         %Draw nodes
         [S.mark_axl,S.mark_sag,S.mark_cor]=draw_nodes(Coor,Net);
         %Assign callback to nodes
         for i=1:length(S.mark_axl)
             set(S.mark_axl(i),'buttondownfcn',{@bdfcn_node,S,i,Net,NBS});
             set(S.mark_sag(i),'buttondownfcn',{@bdfcn_node,S,i,Net,NBS});
             set(S.mark_cor(i),'buttondownfcn',{@bdfcn_node,S,i,Net,NBS});
         end
         %Assign callback to each line
         for i=1:length(S.line_axl)
             try
             set(S.line_axl(i),'buttondownfcn',{@bdfcn_line,S,i,Net,NBS});
             set(S.line_sag(i),'buttondownfcn',{@bdfcn_line,S,i,Net,NBS});
             set(S.line_cor(i),'buttondownfcn',{@bdfcn_line,S,i,Net,NBS});
             catch
                
             end
         end
         %Assign callback to toolbar buttons
         set(S.button(1),'ClickedCallback',{@decrease_node_size,S});
         set(S.button(2),'ClickedCallback',{@increase_node_size,S});
         set(S.button(3),'ClickedCallback',{@decrease_line_width,S});
         set(S.button(4),'ClickedCallback',{@increase_line_width,S});
         set(S.button(5),'ClickedCallback',{@reset,S,Net});
         set(S.button(6),'OnCallback',{@show_all_nodes,S});
         set(S.button(6),'OffCallback',{@undo_show_all_nodes,S});
         
         %Store line and node handles in UserData of the figure
         %This allows subsequnt deletion of the lines and nodes
         set(S.fh,'UserData',[S.mark_axl,S.mark_sag,S.mark_cor,...
             S.line_axl,S.line_sag,S.line_cor]);
         
         
         
         fname=strcat('NBS_network_',int2str(n),...
            '_of_',int2str(NBS.n));
        fname=fullfile(NBS.ResultPath,fname);
        set(gcf,'PaperPositionMode','auto');
        % print('-dpng',fname);
        % print('-djpeg',fname);
        saveas(gcf,fname,'jpg');
         
        pause(1); 
     else
         set(S.ls,'String','No network to display');
     end

 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Draw edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [line_axl,line_sag,line_cor]=draw_edges(Coor,Net)
    %Combine both arrays for checking Bidirectional Connections
    tempNet = [Net.ind_j Net.ind_i];
    %DEfine Color for directed Connections
%     minT =Net.test_stat(Net.ind_i(1):Net.ind_j(1));
%     maxT =Net.test_stat(Net.ind_i(1):Net.ind_j(1));
%     for i=1:length(Net.ind_i)
%         if (Net.test_stat(Net.ind_i(i):Net.ind_j(i))<minT)
%             minT =Net.test_stat(Net.ind_i(i):Net.ind_j(i));
%         elseif (Net.test_stat(Net.ind_i(i):Net.ind_j(i))>maxT)
%             maxT =Net.test_stat(Net.ind_i(i):Net.ind_j(i));
%         end
%     end
    n = 2;
    cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
    %Color define end
      
    
     j = 1;
     for i=1:length(Net.ind_i)
            
          tempEntry= [Net.ind_i(i) Net.ind_j(i)];
          [Result,LocResult] = ismember(tempEntry,tempNet,'rows');
          ab = abs(Net.test_stat(Net.ind_i(i):Net.ind_j(i)));
%           resDoubleDraw=0;
%           if Result == 1
%               try
%               EntryInverse = [Net.ind_j(i) Net.ind_i(i)];
%               %[resDoubleDraw,LocResult] = ismember(EntryInverse,Bidir,'rows');
%               catch
%                   resDoubleDraw = 0;
%               end
%               
%               try
%                    Bidir = [Bidir; tempEntry];
%               catch
%                   Bidir = tempEntry;
%               end    
%           end
%        if resDoubleDraw == 0 %if line not allready drawed           
           if (Result == false )%Draw Directed one sided  

                line_axl(j)=line([Coor.x_axl(Net.ind_i(i)),Coor.x_axl(Net.ind_j(i))],...
                    [Coor.y_axl(Net.ind_i(i)),Coor.y_axl(Net.ind_j(i))],...
                    'LineWidth',(Net.lw));%,...


                line_sag(j)=line([Coor.x_sag(Net.ind_i(i)),Coor.x_sag(Net.ind_j(i))],...
                    [Coor.y_sag(Net.ind_i(i)),Coor.y_sag(Net.ind_j(i))],...
                    'LineWidth',Net.lw);%,...,...

                line_cor(j)=line([Coor.x_cor(Net.ind_i(i)),Coor.x_cor(Net.ind_j(i))],...
                    [Coor.y_cor(Net.ind_i(i)),Coor.y_cor(Net.ind_j(i))],...
                    'LineWidth',Net.lw);%,...,...

                drawnow
                set(line_sag(j).Edge, 'ColorBinding','interpolated', 'ColorData',cd)
                %pause(0.03)
                set(line_axl(j).Edge, 'ColorBinding','interpolated', 'ColorData',cd)
                %pause(0.03)
                set(line_cor(j).Edge, 'ColorBinding','interpolated', 'ColorData',cd)

                %pause(0.03)
                j=j+1;
           else%Draw directed Both sided
               line_axl(j)=line([Coor.x_axl(Net.ind_i(i)),Coor.x_axl(Net.ind_j(i))],...
                [Coor.y_axl(Net.ind_i(i)),Coor.y_axl(Net.ind_j(i))],...
                'LineWidth',Net.lw+1,...
                'Color','g');%Net.line_color);

                line_sag(j)=line([Coor.x_sag(Net.ind_i(i)),Coor.x_sag(Net.ind_j(i))],...
                    [Coor.y_sag(Net.ind_i(i)),Coor.y_sag(Net.ind_j(i))],...
                    'LineWidth',Net.lw,...
                    'Color','g');

                line_cor(j)=line([Coor.x_cor(Net.ind_i(i)),Coor.x_cor(Net.ind_j(i))],...
                    [Coor.y_cor(Net.ind_i(i)),Coor.y_cor(Net.ind_j(i))],...
                    'LineWidth',Net.lw,...
                    'Color','g');
                 j=j+1;
           end
%        end
     end
    strFrom = {'from'};
    strTo = {'to'};
    line_Dir=line([270,340],[190,190],...
                'LineWidth',2);
    drawnow
    set(line_Dir.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    text(275,185,strFrom)
    text(325,185,strTo)
    
    strbi = {'bidirectional'};
    line_Bi=line([270,340],[170,170],...
                'LineWidth',2,...
                'Color','g');%Net.line_color);
    text(280,165,strbi)        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mark_axl,mark_sag,mark_cor]=draw_nodes(Coor,Net)
    for i=1:Net.N
        %Only draw node if a line orignates from it
        
        
        
        
        if ismember(i,[Net.ind_i;Net.ind_j])
            mark_axl(i)=line(Coor.x_axl(i),Coor.y_axl(i),...
                'Marker','.','MarkerSize',Net.nw,'Color',Net.node_color);
            mark_sag(i)=line(Coor.x_sag(i),Coor.y_sag(i),...
                'Marker','.','MarkerSize',Net.nw,'Color',Net.node_color);
            mark_cor(i)=line(Coor.x_cor(i),Coor.y_cor(i),...
                'Marker','.','MarkerSize',Net.nw,'Color',Net.node_color);
        else
            mark_axl(i)=line(Coor.x_axl(i),Coor.y_axl(i),...
                'Marker','.','MarkerSize',Net.nw,...
                'Visible','off','Color',Net.node_color);
            mark_sag(i)=line(Coor.x_sag(i),Coor.y_sag(i),...
                'Marker','.','MarkerSize',Net.nw,...
                'Visible','off','Color',Net.node_color);
            mark_cor(i)=line(Coor.x_cor(i),Coor.y_cor(i),...
                'Marker','.','MarkerSize',Net.nw,...
                'Visible','off','Color',Net.node_color);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Disable highlight in listbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bdfcn_ls(varargin)
    S = varargin{3};
    set(S.ls,'Value',[]);
 
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Buttondown for node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=bdfcn_node(varargin)
    % src - the object that is the source of the event
    % evnt - empty for this property
    [S,i,Net,NBS] = varargin{[3,4,5,6]};
    sel_typ = get(gcbf,'SelectionType');
    switch sel_typ
        case 'normal'
            %disp('User clicked left-mouse button')
            i_node=find(Net.clicked_node);
            [i_edge,j_edge]=find(Net.clicked_line);
            set(S.mark_axl(:),'Color',Net.node_color);
            set(S.mark_sag(:),'Color',Net.node_color);
            set(S.mark_cor(:),'Color',Net.node_color);
            set(S.line_axl(:),'LineWidth',Net.lw);%'Color',Net.line_color);
            set(S.line_sag(:),'LineWidth',Net.lw);%'Color',Net.line_color);
            set(S.line_cor(:),'LineWidth',Net.lw);%'Color',Net.line_color);
            if isempty(i_node) && isempty(i_edge)
                %Nothing previously clicked
                set(S.mark_axl(i),'Color','red');
                set(S.mark_sag(i),'Color','red');
                set(S.mark_cor(i),'Color','red');
                if isfield(NBS,'node_label')
                    ls=get(S.ls,'String');
                    set(S.ls,'String',[ls;{NBS.node_label{i}}]);
                else
                    ls=get(S.ls,'String');
                    set(S.ls,'String',[ls;{' '}]);
                end
                Net.clicked_node(i)=1;
            else
                %Node or edge previously clicked
                ls=get(S.ls,'String');
                [x,y] =size(ls);
               
                if x == 7
                    set(S.ls,'String',ls(1:end-2,:));
                elseif x == 6 
                    set(S.ls,'String',ls(1:end-1,:));
                end
               
                if i_node==i
                    %Same node clicked
                else
                    set(S.mark_axl(i),'Color','red');
                    set(S.mark_sag(i),'Color','red');
                    set(S.mark_cor(i),'Color','red');
                    if isfield(NBS,'node_label')
                        ls=get(S.ls,'String');
                        set(S.ls,'String',[ls;{NBS.node_label{i}}]);
                    else
                        ls=get(S.ls,'String');
                        set(S.ls,'String',[ls;{' '}]);
                    end
                    Net.clicked_node(i)=1;
                end
            end
            %Update Net
            Net.clicked_line(i_edge,j_edge)=0;
            Net.clicked_node(i_node)=0;
            for i=1:length(S.mark_axl)
                set(S.mark_axl(i),'buttondownfcn',{@bdfcn_node,S,i,Net,NBS});
                set(S.mark_sag(i),'buttondownfcn',{@bdfcn_node,S,i,Net,NBS});
                set(S.mark_cor(i),'buttondownfcn',{@bdfcn_node,S,i,Net,NBS});
            end
            for i=1:length(S.line_axl)
                set(S.line_axl(i),'buttondownfcn',{@bdfcn_line,S,i,Net,NBS});
                set(S.line_sag(i),'buttondownfcn',{@bdfcn_line,S,i,Net,NBS});
                set(S.line_cor(i),'buttondownfcn',{@bdfcn_line,S,i,Net,NBS});
            end
        case 'extend'
            %disp('User did a shift-click');
            %set(src,'Selected','on')
        case 'alt'
            %disp('User did a control-click');
            %set(src,'Selected','on')
            %set(src,'SelectionHighlight','off')
    end
end            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Buttondown for edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=bdfcn_line(varargin)
    % src - the object that is the source of the event
    % evnt - empty for this property
    [S,i,Net,NBS] = varargin{[3,4,5,6]};
     %Combine both arrays for checking Bidirectional Connections
    tempNet = [Net.ind_j Net.ind_i];%TmpNed is inverse
    %DEfine Color for directed Connections
    n=2;
    cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
    tempEntry= [Net.ind_i(i) Net.ind_j(i)];
    [ResultBidirectional,LocResult] = ismember(tempEntry,tempNet,'rows');
    
    sel_typ = get(gcbf,'SelectionType');
    switch sel_typ
        case 'normal'
            %disp('User clicked left-mouse button')
            if isfield(Net,'test_stat') && isfield(NBS,'node_label')
                tstat=sprintf(' %0.2f',Net.test_stat(Net.ind_i(i),Net.ind_j(i)));
                str=[NBS.node_label{Net.ind_i(i)},' - ',NBS.node_label{Net.ind_j(i)},'; test stat:',tstat];
                if ResultBidirectional == 1
                    tstat=sprintf(' %0.2f',Net.test_stat(Net.ind_i(LocResult),Net.ind_j(LocResult)));
                    str2= [NBS.node_label{Net.ind_j(i)},' - ',NBS.node_label{Net.ind_i(i)},'; test stat:',tstat];
                    str = append(str , newline, str2);
                end
            elseif isfield(NBS,'node_label')
                str=[NBS.node_label{Net.ind_i(i)},' - ',NBS.node_label{Net.ind_j(i)}];
            elseif isfield(Net,'test_stat')
                str=sprintf('test stat: %0.2f',Net.test_stat(Net.ind_i(i),Net.ind_j(i)));
            else
                str=' ';
            end
            drawnow
            i_node=find(Net.clicked_node);
            [i_edge,j_edge]=find(Net.clicked_line);
            set(S.mark_axl(:),'Color',Net.node_color);
            set(S.mark_sag(:),'Color',Net.node_color);
            set(S.mark_cor(:),'Color',Net.node_color);
            set(S.line_axl(:),'LineWidth',Net.lw);%'Color',Net.line_color);
            set(S.line_sag(:),'LineWidth',Net.lw);%'Color',Net.line_color);
            set(S.line_cor(:),'LineWidth',Net.lw);%'Color',Net.line_color);
            if isempty(i_node) && isempty(i_edge)
                %Nothing previously clicked
                drawnow
                set(S.line_axl(i),'LineWidth',7);%'Color','red');
                set(S.line_sag(i),'LineWidth',7);%'Color','red');
                set(S.line_cor(i),'LineWidth',7);%'Color','red');
                drawnow
                ls=get(S.ls,'String');
                set(S.ls,'String',[ls;{str}]);
                Net.clicked_line(Net.ind_i(i),Net.ind_j(i))=1;
            else
                %Node or edge previously clicked
                ls=get(S.ls,'String');
                [x,y] =size(ls);
               
                if x == 7
                    set(S.ls,'String',ls(1:end-2,:));
                elseif x == 6 
                    set(S.ls,'String',ls(1:end-1,:));
                end
                
                if i_edge==Net.ind_i(i) & j_edge==Net.ind_j(i)
                    %Same node clicked
                else
                    drawnow
                    set(S.line_axl(i),'LineWidth',7);%'Color','red');                   
                    set(S.line_sag(i),'LineWidth',7);%'Color','red');
                    set(S.line_cor(i),'LineWidth',7);%'Color','red');
                    ls=get(S.ls,'String');
                    set(S.ls,'String',[ls;{str}]);
                    Net.clicked_line(Net.ind_i(i),Net.ind_j(i))=1;
                end
            end
            drawnow
            %Update Net
            Net.clicked_line(i_edge,j_edge)=0;
            Net.clicked_node(i_node)=0;
            for i=1:length(S.mark_axl)
                set(S.mark_axl(i),'buttondownfcn',{@bdfcn_node,S,i,Net,NBS});
                set(S.mark_sag(i),'buttondownfcn',{@bdfcn_node,S,i,Net,NBS});
                set(S.mark_cor(i),'buttondownfcn',{@bdfcn_node,S,i,Net,NBS});
            end
            for i=1:length(S.line_axl)
                set(S.line_axl(i),'buttondownfcn',{@bdfcn_line,S,i,Net,NBS});
                set(S.line_sag(i),'buttondownfcn',{@bdfcn_line,S,i,Net,NBS});
                set(S.line_cor(i),'buttondownfcn',{@bdfcn_line,S,i,Net,NBS});
            end
        case 'extend'
            %disp('User did a shift-click');
            %set(src,'Selected','on')
        case 'alt'
            %disp('User did a control-click');
            %set(src,'Selected','on')
            %set(src,'SelectionHighlight','off')
    end
          

end

