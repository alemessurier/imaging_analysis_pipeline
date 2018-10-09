function [f,im] = stackDispGui( imstack )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%IMDISPGUI opens image with slider for interactively changing brightness
%  INPUTs: im - images converted to double (not scaled)
%
% AML 2014

f = figure('Visible','off');
ax = axes('Units','pixels');
%     normstack=zeros(size(imstack));
%     blstackTmp=imstack(:,:,1:19);
%     maxpixtmp=max(blstackTmp(:));
%     for j=1:19;
%         thisFrame=imstack(:,:,j);
%         normstack(:,:,j)=thisFrame/maxpixtmp;
%     end
%
%     for k=20:size(imstack,3)
%         thisFrame=imstack(:,:,k);
%         blstack=imstack(:,:,(k-19):k);
%         maxpix=max(blstack(:));
%         normstack(:,:,k)=thisFrame/maxpix;
%     end
%         imstack=uint8(imstack);

%     im=im/max(max(im));
maxLim=0.5;
im=imstack(:,:,1)
avg_im=mean(imstack,3);
avg_im=avg_im/max(max(avg_im));
avg_im=imadjust(avg_im);
imshow(im);
axis square


% Create slider
sldFrame = uicontrol(f,'Style', 'slider',...
    'Min',1,'Max',size(imstack,3),'Value',1,...
    'SliderStep',[1/size(imstack,3) 10/size(imstack,3)],...
    'Position', [100 10 300 20],...
    'Callback', @scrollFrame);

% Add a text uicontrol to label the slider.
txt = uicontrol(f,'Style','text',...
    'Position',[200 30 100 20 ],...
    'String','Frame num');

sldBright=uicontrol(f,'Style', 'slider',...
    'Min',0,'Max',1,'Value',.5,...
    'Units','pixels',...
    'Position', [40 10 20 300],...
    'Callback', @adjustLim);

% Add a text uicontrol to label the slider.
txt2 = uicontrol(f,'Style','text',...
    'Position',[40 300 100 20 ],...
    'String','Max pixel value');

avgButton=uicontrol(f,'Style','checkbox','Units','normalized','Position', [0.9 0.4 0.1 0.1], 'String', 'avg', ...
    'Callback', @showAvg);%'UserData',1,'String','done!')

%     Handles=guihandles(f)
% Make figure visble after adding all components
set(f,'Visible','on')


    function im=scrollFrame(hObject,EventData,Handles)
        i = ceil(sldFrame.Value);
        im=imstack(:,:,i);
        maxLim=sldBright.Value;
        im=imadjust(im,[0 maxLim],[]);
        imshow(im)
    end

    function maxLim=adjustLim(hObject,EventData,Handles)
        maxLim = get(hObject,'Value');
        imNew=imadjust(im,[0 maxLim],[]);
        imshow(imNew)
    end

    function showAvg(hObject,EventData,Handles)
        if get(hObject,'Value')==get(hObject,'Max')
            im=avg_im;
            imshow(im)
        else
            i = ceil(sldFrame.Value);
            im=imstack(:,:,i);
            maxLim=sldBright.Value;
            im=imadjust(im,[0 maxLim],[]);
            imshow(im)
        end
    end

end

