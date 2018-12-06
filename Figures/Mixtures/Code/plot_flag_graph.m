%Better looking graph
function h = plot_flag_graph(Z, sigma, width, options)

options = fill_default_options(options);

res = size(Z,1)-1;
y1 = linspace(-width,width,res);
y2 = linspace(-width,width,res);
y1 = [y1,y1(res)+y1(res) - y1(res-1)];  %*awful hack so we display final row and column of data in plot
y2 = [y2,y2(res)+y2(res) - y2(res-1)];  %*
[Y1,Y2] = meshgrid(y1,y2);

h = pcolor(Y1,Y2,Z);
shading flat
axis equal tight
% cmap = gray(3);
% cmap = parula(3);
% cmap = summer(3);
% cmap = [0.3 0.3 0.3;0.2 0.4 0.3;0.3 0.3 0.3];
cmap = [255 189 76; 180 114 0;103 65 0]/255;    %Light to Dark Highlight
cmap = [103 65 0;180 114 0;255 189 76]/255;    %Dark to Light Highlight
cmap = [0 51 104;38 82 127;76 164 255]/255; %Darkt to Light Primary
% cmap = [76 164 255;38 82 127;0 51 104]/255; %Light to Dark Primary
colormap(cmap);
hcb = colorbar;

% set(hcb,'YTick',[1.333,2,2.667],'YTickLabel',[1,2,3],'Location','EastOutside','position',[0.73 .11 .05 .3],'FontSize',20)
set(hcb,'YTick',[1.333,2,2.667],'YTickLabel',[1,2,3],'Location','EastOutside','FontSize',20)

% Labels
if options.n == 3
    xlabel('x_2','FontSize',20)
    ylabel('x_3','FontSize',20)
end
if options.n == 2
    xlabel('x_1','FontSize',20)
    ylabel('x_2','FontSize',20)
end
% title('Number of probability masses required (x_1 = 0)','FontSize',20)

mid = 0 + width/(res-1);

%Axis lines
if options.axis_lines
    line([mid mid], [min(y1) max(y1)],'Color','k');  %y-axis
    line([min(y2) max(y2)], [mid mid],'Color','k');  %x-axis
end

% Ellipse
if options.ellipse
    colour = 'r';
    colour = [255 189 76]/255;
%     options.line_width = 2;
    
    t = linspace(0,2*pi,100);
    b = sqrt(3)*sigma;
    a = sqrt(3)*b; %What we would expect if this ellipse was the intersection
    %of the plane with a cyclinder along the axis x = y = z
    ellipse_x = (a*cos(t) - b*sin(t))/sqrt(2)+mid;
    ellipse_y = (a*cos(t) + b*sin(t))/sqrt(2)+mid;
    hold on
    plot(ellipse_x,ellipse_y,'Color',colour,'LineWidth',options.line_width)
    hold off
end

%Distance 2*sigma lines
if options.sigma_lines
%     colour = 'k';
%     colour = 'w';
    colour = [0.8,0.8,0.8];
%     options.line_width = 2;
    
%     dist = 2*sigma;
        dist = 2*sigma;   %For normal
%     dist = 2*sqrt(1/3)*sigma;   %For Cauchy
    
    hold on
    
    line([mid+dist mid+dist], [min(y1) max(y1)],'Color',colour,'LineWidth',options.line_width);  
    line([mid-dist mid-dist], [min(y1) max(y1)],'Color',colour,'LineWidth',options.line_width);  

    line([min(y1) max(y1)],[mid-dist mid-dist],'Color',colour,'LineWidth',options.line_width); 
    line([min(y1) max(y1)],[mid+dist mid+dist],'Color',colour,'LineWidth',options.line_width);

    shift = dist;
    x_line_1 = linspace(min(y1),max(y1)-shift,100);
    x_line_2 = linspace(min(y1)+shift,max(y1),100);
    y_line_1 = x_line_1 + shift;
    y_line_2 = x_line_2 - shift;
    plot(x_line_1,y_line_1,'Color',colour,'LineWidth',options.line_width)
    plot(x_line_2,y_line_2,'Color',colour,'LineWidth',options.line_width)
    hold off
end

if options.m2lines
%     colour = 'k';
%     colour = 'w';
    colour = [0.8,0.8,0.8];
    colour = [255 189 76]/255;
%     options.line_width = 2;
    dist = sqrt(2)*sqrt(3)*sigma;   %For normal
    
    hold on
    
    line([mid+dist mid+dist], [min(y1) max(y1)],'Color',colour,'LineWidth',options.line_width);  
    line([mid-dist mid-dist], [min(y1) max(y1)],'Color',colour,'LineWidth',options.line_width);  

    line([min(y1) max(y1)],[mid-dist mid-dist],'Color',colour,'LineWidth',options.line_width); 
    line([min(y1) max(y1)],[mid+dist mid+dist],'Color',colour,'LineWidth',options.line_width);

    shift = dist;
    x_line_1 = linspace(min(y1),max(y1)-shift,100);
    x_line_2 = linspace(min(y1)+shift,max(y1),100);
    y_line_1 = x_line_1 + shift;
    y_line_2 = x_line_2 - shift;
    plot(x_line_1,y_line_1,'Color',colour,'LineWidth',options.line_width)
    plot(x_line_2,y_line_2,'Color',colour,'LineWidth',options.line_width)
    hold off
end

% Linsday1983b bound
if options.lindsay_bound
    colour = 'w';
    colour = [255 217 152]/255;
    colour = [255 189 76]/255;
%     options.line_width = 2;
    
    dist = 2*sigma;   %For normal
%     dist = 2*sqrt(1/3)*sigma;   %For Cauchy
    
    line([dist+mid,mid],[dist+mid,dist+mid],'Color',colour,'LineWidth',options.line_width);
    line([dist+mid,dist+mid],[mid,dist+mid],'Color',colour,'LineWidth',options.line_width);

    line([-dist+mid,mid],[-dist+mid,-dist+mid],'Color',colour,'LineWidth',options.line_width);
    line([-dist+mid,-dist+mid],[mid,-dist+mid],'Color',colour,'LineWidth',options.line_width);

    line([mid,dist+mid],[-dist+mid,mid],'Color',colour,'LineWidth',options.line_width)
    line([-dist+mid,mid],[mid,dist+mid],'Color',colour,'LineWidth',options.line_width)
end

if options.save
    %Print to file... Make sure that the current folder is where you want to
    %save it!!!
    if strcmp(options.filename, '')
        filename = ['Sigma',num2str(sigma),'n3res',num2str(res),'width',num2str(width)];
        if options.ellipse
            filename = [filename,'ellipse'];
        end
        if options.lindsay_bound
            filename = [filename,'lindsaybound'];
        end
        filename = [filename,'-blue','.png'];
    else
        filename = options.filename;
    end

    % export_fig(filename,'-png','-transparent','-r600')
    % export_fig(filename,'-pdf')
    saveas(h, filename)
end
end

function options = fill_default_options(options)
    if ~isfield(options, 'ellipse')
        options.ellipse = false;
    end
    if ~isfield(options, 'lindsay_bound')
        options.lindsay_bound = false;
    end
    if ~isfield(options, 'sigma_lines')
        options.sigma_lines = false;
    end
    if ~isfield(options, 'save')
        options.save = false;
    end
    if ~isfield(options, 'axis_lines')
        options.axis_lines = true;
    end
    if ~isfield(options, 'line_width')
        options.line_width = 1;
    end
    if ~isfield(options, 'filename')
        options.filename = '';
    end
    if ~isfield(options, 'n')
        options.n = 3;
    end
    if ~isfield(options,'m2lines')
        options.m2lines = false;
    end
end