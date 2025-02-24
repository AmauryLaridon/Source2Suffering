function mf_XMinorTick(n_minortick,color,enable)
% Syntax: mf_XMinorTick
%
% Input:    
% {opt}     n_minortick  - number of steps of XMinorTick [values (numeric)]
%                          five steps by default
% 
% {opt}     color        - color used in XMinorTick, black by default
%
% {opt}     enable       - disable or enable XMinorTick, enbale='on' by 
%                          default
%
% Example 1:  Without input arguments
%    mf_XMinorTick
%
% Example 2:  Without Color input argument
% 
%    mf_XMinorTick(5)
%
% Example 3:  Set steps of XMinorTick and color XMinorTick
%
%    mf_XMinorTick(3,'r')
%
% Example 4:  Set disable  XMinorTick
%
%    mf_XMinorTick([],[],'off')
%
%
% mf_XMinorTick is part of  Flow uMet-LBA Qa/Qc system Toolbox 
% INPA - National Institute of Amazonian Research 
% LBA  - Large Scale Biosphere-Atmosphere Experiment in Amazonia
%
% Development : Paulo R. Teixeira and Team Micromet.
% FeedBack : paulo.ricardo.teixeira at gmail.com
%
% Note: Fix error - XMinorTick above XTick (solve: with NaN and change
% layer)
% Note: Add function - XMinorTick disable option

if strcmp(get(gca,'XminorTick'),'off')==0
    warning('MATLAB:mf_XMinorTick',...
        'XMinorTick Active! \n XMinorTick will be disabled')
    set(gca,'XMinorTick','off');
end
 if nargin<3;enable={'on'};end
if nargin<2;color= get(gca,'XColor');end
if nargin<1;n_minortick=5;color= get(gca,'XColor');end

if (ischar(enable) && strcmpi(enable,'off'));...
        try
        h = findobj(gca,'Type','line');
        delete(h(1:end-1,1))    
        return;
        catch err
            error([mfilename ':ErrorDelXMinorTick'],...
                'Could not disable XMinorTick');
        end
else
    
end

try
    x_value=get(gca,'Xtick');
    incremet=x_value(2)/n_minortick;
    new_xminortick=[min(x_value)+incremet:incremet:max(x_value)-incremet]';
    
    for j=1:size(x_value,2)
        id=new_xminortick==x_value(1,j);
        new_xminortick(id)=NaN;
        set(gca, 'Layer','top')
    end
    
    new_xminortick(:,2)=new_xminortick(:,1);
    
catch err
    error([mfilename ':ErrorRunFile'],'Could not run file');
end

y_value=get(gca,'YTick');
% height=[zeros(size(new_xminortick,1),1) ...
%     ones(size(new_xminortick,1),1).*min(get(gca,'Ytick'))+...
%     (y_value(2)-y_value(1))/10];
height=[zeros(size(new_xminortick,1),1) ...
    ones(size(new_xminortick,1),1).*min(get(gca,'Ytick'))+...
    (y_value(2)-y_value(1))/10] ./ 1.5;  % hack by Wim Thiery

for i=1:size(new_xminortick,1)
    line(new_xminortick(i,:),height(i,:),'Color',color)
end

end