function  [rawTimeSeries_new,movies_exclude]=check_rawFluor( rawTimeSeries )


%% plot all traces
fns=fieldnames(rawTimeSeries);
spacing =100;
cellNames=fieldnames(rawTimeSeries.(fns{1}));
inds_ex=zeros(1,length(fns));
j=1;
g= figure;
% Create push button
keep_btn = uicontrol('Style', 'pushbutton', 'String', 'keep',...
    'Units','normalized','Position', [0.55 0.005 0.045 0.05],...
    'Callback', @keep_movie);

ex_btn = uicontrol('Style', 'pushbutton', 'String', 'exclude',...
    'Units','normalized','Position', [0.45 0.005 0.045 0.05],...
    'Callback', @exclude_movie);






fn=fns{j};
for i=1:length(cellNames)
    
    cn = cellNames{i};
    
    %         samp_rate=sampRate(j);
    %         time = [1:length(rawTimeSeries.(fn).(cn))]/samp_rate;
    plot(1:length(rawTimeSeries.(fn).(cn)), rawTimeSeries.(fn).(cn)+((i-1)*spacing),'k');
    hold on
end
%     ylabel('Fluorescence (AU)'); xlabel('Time (seconds)');
% axis([1 length(rawTimeSeries.(fn).(cn)) 0 (length(fns)*spacing)]);
%      scrollplot(200, 1:length(deltaF.(fn).(cn)), deltaF.(fn).(cn));


set(gca,'YTick',[0:spacing:((length(cellNames)*spacing))])
set(gca,'YTickLabel',cellNames)
title([num2str(j),'/',num2str(length(fns))])
% Pause until figure is closed ---------------------------------------%
waitfor(g);
inds_ex=logical(inds_ex);
inds_ex=inds_ex(1:length(fns));
movies_exclude=fns(inds_ex);
rawTimeSeries_new=rmfield(rawTimeSeries,movies_exclude);

    function keep_movie(source,event)
        j=j+1;
        if j<=length(fns)
            j<=length(fns)
            fn=fns{j};
            figure(g); hold off;
            for i=1:length(cellNames)
                
                cn = cellNames{i};
                
                %         samp_rate=sampRate(j);
                %         time = [1:length(rawTimeSeries.(fn).(cn))]/samp_rate;
                plot(1:length(rawTimeSeries.(fn).(cn)), rawTimeSeries.(fn).(cn)+((i-1)*spacing),'k');
                hold on
            end
            title([num2str(j),'/',num2str(length(fns))])

        end
    end

    function exclude_movie(source,event)
        inds_ex(j)=1;
        j=j+1;
        if j<=length(fns)
            fn=fns{j};
            figure(g); hold off;
            for i=1:length(cellNames)
                
                cn = cellNames{i};
                
                %         samp_rate=sampRate(j);
                %         time = [1:length(rawTimeSeries.(fn).(cn))]/samp_rate;
                plot(1:length(rawTimeSeries.(fn).(cn)), rawTimeSeries.(fn).(cn)+((i-1)*spacing),'k');
                hold on
            end
            title([num2str(j),'/',num2str(length(fns))])

        end
    end

end

