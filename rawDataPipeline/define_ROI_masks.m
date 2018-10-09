function outputPos= define_ROI_masks(im,inputPos)
% Creats gui for defining and saving ROIs.
%
% INPUTS: imPath - path for individual registered movie (string)
%         inputPos - cell array of previously defined ROIs
%
% OUTPUTS: outputPos - array of ROIs, defined as image masks
%
% Amy LeMessurier 2014


if nargin==0;
    [fname,imPath] = uigetfile('.tif');
    imPath=strcat(imPath,fname);
    [im, meta] = LoadTIFF_SI5(imPath);
    inputPos=zeros(512,512);
end
  xpix=size(im,1);
  ypix=size(im,2);
  
set(0,'DefaultLineMarkerSize', 2)

[f,curr_im]=stackDispGui(im);
zoom on
hold on

% if isempty(inputPos)
%     inputPos=zeros(512,512);
% end

green = cat(3, zeros(size(curr_im)),ones(size(curr_im)), zeros(size(curr_im)));
            h=imshow(green);
            set(h,'AlphaData',0)
showButton=uicontrol(f,'Style','checkbox','Units','normalized','Position', [0.4 0.9 0.2 0.1], 'String', 'show ROIs', ...
    'Callback', @showROIs);
addButton=uicontrol(f,'Style','pushbutton','Units','normalized','Position', [0.1 0.9 0.2 0.1], 'String', 'add ROI', ...
    'Callback', @appendROI);%'UserData',1,'String','done!')

shiftU= uicontrol(f,'Style','pushbutton','Units','normalized','Position', [0.9 0.2 0.05 0.05], 'String','U',...
    'Callback', @shift_up);
shiftD= uicontrol(f,'Style','pushbutton','Units','normalized','Position', [0.9 0.1 0.05 0.05], 'String', 'D', ...
    'Callback', @shift_down);
shiftR= uicontrol(f,'Style','pushbutton','Units','normalized','Position', [0.95 0.15 0.05 0.05], 'String', 'R', ...
    'Callback', @shift_right);
shiftL= uicontrol(f,'Style','pushbutton','Units','normalized','Position', [0.85 0.15 0.05 0.05], 'String', 'L', ...
    'Callback', @shift_left);

rotateCC=uicontrol(f,'Style','pushbutton','Units','normalized','Position', [0.85 0.75 0.06 0.05], 'String', 'CC', ...
    'Callback', @rotate_cc);
rotateC=uicontrol(f,'Style','pushbutton','Units','normalized','Position', [0.92 0.75 0.06 0.05], 'String', 'C', ...
    'Callback', @rotate_c);


doneButton= uicontrol(f, 'Units', 'normalized', 'Position',[0.7 0.9 0.2 0.1], 'String', 'Done', 'Callback', @done);
% Static text


    function showROIs(hObject,EventData,Handles)
        
        if get(hObject,'Value')==get(hObject,'Max')
           
            masks=sum(inputPos,3);
            masks(masks>0)=1;
%             imshow(curr_im)

            hold on
            uistack(h,'top')
            set(h,'AlphaData',masks*0.5)
           
%         else
%             if ~isempty(r)
%                 cellfun(@(x)set(x,'Visible','off'),r,'Uni',0);
% %                 delete(r{:})
%             else
%             end
        end
    end


    function appendROI(hObject,EventData,Handles)
        ROI=imellipse;
        wait(ROI)
        thisROI=getVertices(ROI);
        thisROImask=poly2mask(thisROI(:,1),thisROI(:,2),xpix,ypix);
        zoom on
        inputPos=cat(3,inputPos, thisROImask);
        setColor(ROI,'c');
        %rectangle('Position',thisROI,'Curvature',[1 1],'FaceColor','c')
        hold on
    end
   function shift_up(hObject,EventData,Handles)
       inputPos=inputPos(2:end,:,:);
            inputPos=cat(1,inputPos,zeros(1,size(inputPos,2),size(inputPos,3))); 
       masks=sum(inputPos,3);
            masks(masks>0)=1;
            uistack(h,'top')
            set(h,'AlphaData',masks*0.5)
            
    end

    function shift_down(hObject,EventData,Handles)
        inputPos=inputPos(1:end-1,:,:);
            inputPos=cat(1,zeros(1,size(inputPos,2),size(inputPos,3)),inputPos); 
       masks=sum(inputPos,3);
            masks(masks>0)=1;
            uistack(h,'top')
            set(h,'AlphaData',masks*0.5)
    end

    function shift_right(hObject,EventData,Handles)
       inputPos=inputPos(:,1:end-1,:);
            inputPos=cat(2,zeros(size(inputPos,1),1,size(inputPos,3)),inputPos); 
       masks=sum(inputPos,3);
            masks(masks>0)=1;
            uistack(h,'top')
            set(h,'AlphaData',masks*0.5)
    end

    function shift_left(hObject,EventData,Handles)
         inputPos=inputPos(:,2:end,:);
            inputPos=cat(2,inputPos,zeros(size(inputPos,1),1,size(inputPos,3))); 
       masks=sum(inputPos,3);
            masks(masks>0)=1;
            uistack(h,'top')
            set(h,'AlphaData',masks*0.5);
    end

 function rotate_c(hObject,EventData,Handles)
     inputPos=imrotate(inputPos,-1,'nearest','crop');
     masks=sum(inputPos,3);
            masks(masks>0)=1;
            uistack(h,'top')
            set(h,'AlphaData',masks*0.5);
 end

function rotate_cc(hObject,EventData,Handles)
 inputPos=imrotate(inputPos,1,'nearest','crop');
     masks=sum(inputPos,3);
            masks(masks>0)=1;
            uistack(h,'top')
            set(h,'AlphaData',masks*0.5); 
end
    
    function done(hObject,eventdata,Handles)
        % Assign Output
%         emptyInds=sum(sum(inputPos));
%         emptyInds=emptyInds(:)>0;
        outputPos = inputPos;%(:,:,emptyInds);
        % Close figure
        delete(f); % close GUI
    end

% Pause until figure is closed ---------------------------------------%
waitfor(f);
end