function [pt_point_positions, tikz_support_points, colorbar_limits] = surf3tikz(h_figure, export_name, cfg, debug)
%SURF3TIKZ converts a surf plot into a tikz file
% This function takes a figure containing a surf plot and creates a tikz file and an according png
% file. To achieve this, it automatically places data markers on points of the axes box around
% the plot and subsequently detects their pixel position.
% For views where only two of the three dimensions are visible a simpler approach is used which creates
% a 2D plot in TikZ, instead of a 3D one. This prevents unstable behaviour because of
% imperfect point positions that especially occur in 90° elevation (bird's eye view) cases. This
% behavior can be suppressed.
% Additionally the key values are available as output parameters to build your own tikz file from
% scratch.
%
% INPUT:
%   h_figure: Figure handle containing the figure plot. Ideally also contains a colorbar.
%   export_name: Base filename of the exported tikz and png file the function writes to
%   cfg: optional config struct
%      .export_dpi: dpi of exported png file, default: 300 dpi
%      .write_png: boolean to optionally suppress png output, default: true
%      .write_tikz: boolean to optionally suppress tikz file output, default: true
%      .write_fig: optional boolean to save the original figure as fig file, default: false
%      .screen_ppi: pixel per inch screen resolution, default: will try to determine system setting
%       (likely values are 96 or 72)
%      .inch2point_ratio: ratio of inch to point, default: 1/72 inch per point
%      .box_point_idc: Indices of the box points to use (1-8). Normally it should be fine to just
%       let the algorithm choose these points. But you can try other combinations if you want for more
%       information on how to choose these points see below in the code, default: []
%      .force_exact_values: per default, surf3tikz will round the pt_point_positions to the first
%       full pt, because otherwise the PGFPlot process seems to be stuck in certain cases. If you
%       don't want surf3tikz to do this, set this parameter to true, default: false
%      .force_3d: For views where one dimension is hidden to the viewer, a 2D approach is used by default.
%       To force the use of the 3D approach set this to true, default: false
%      .use_imagesc: For top down views, imagesc can be used as a way to save space for the image file
%       file, default: false, WARNING: EXPERIMENTAL AND USES THE FIRST SURFACE PLOT IT CAN FIND
%      .print_all: Per default, this function will write line plots into CSV files or writes them to the
%       TikZ file as coordinates. Overwrite this parameter, if instead a 1:1 graphics copy of your plots
%       should be created. Default: false
%      .ext_all: Per default, only line plots with more than ten points are written to CSV files.
%       You can force the function to export all line plots. Default: false
%      .inline_all: Per default, only line plots with fewer than ten points are written as inline
%       coordinates to the TikZ file. You can force the function to write all line plots to the TikZ
%       file instead. Default: false
%      .inline_limit: If you want to change the limit up to which amount of points, inline
%       coordinates should be used for line plots, you can set it with this option. Default: 10
%   debug: boolean, will switch on/off an additional debug plot with the set cursors and the found
%          pixel positions denoted by a white pixel ideally in the middle of every black cursor
%          rectangle, default: false
%
% OUTPUT:
%   pt_point_positions: the data point positions in the pt unit
%   act_data_points: the data point coordinates
%   colorbar_limits: cdata limits of the colorbar currently used
%
% NOTES:
% * If you try to run this script on a surf with a view of El=90°, tikz will throw a dimension error.
% * The choice of box points right now is still not fully tested. It seems to work for now, but
%   there is a chance that this has to be refined in the near future.
% * There is a point to be made to use imagesc instead of print in some situations. However, for now
%   it's not implemented in here yet.
%
% Written by Johannes Schlichenmaier, 2016
%

%% basic input argument handling
if nargin < 2
    error('You''ll need to specify a figure handle and an export base filename!');
end

if nargin < 3
    cfg = struct();
end

if ~isfield(cfg, 'export_dpi')
    cfg.export_dpi = 300;
end

if ~isfield(cfg, 'write_png')
    cfg.write_png = true;
end

if ~isfield(cfg, 'write_tikz')
    cfg.write_tikz = true;
end

if ~isfield(cfg, 'write_fig')
    cfg.write_fig = false;
end

if ~isfield(cfg, 'screen_ppi')
    cfg.screen_ppi = get(0,'ScreenPixelsPerInch');
end

if ~isfield(cfg, 'inch2point_ratio')
    cfg.inch2point_ratio = 1/72;
end

if ~isfield(cfg, 'box_point_idc')
    cfg.box_point_idc = [];
end

if ~isfield(cfg, 'force_exact_values')
    cfg.force_exact_values = false;
end

if ~isfield(cfg, 'force_3d')
    cfg.force_3d = false;
end

if ~isfield(cfg, 'use_imagesc')
    cfg.use_imagesc = false;
end

if ~isfield(cfg, 'print_all')
    cfg.print_all = false;
end

if ~isfield(cfg, 'ext_all')
    cfg.ext_all = false;
end

if ~isfield(cfg, 'inline_all')
    cfg.inline_all = false;
end

if ~isfield(cfg, 'inline_limit')
    cfg.inline_limit = 10;
end

if nargin < 4
    debug = false;
end

if cfg.print_all && cfg.use_imagesc
    warning('Setting both the print_all and use_imagesc options to true may result in undesired outcome!');
end

if cfg.print_all && (cfg.ext_all || cfg.inline_all)
    warning('Both the print_all and ext_all or inline_all options are set. Ignoring ext_all or inline_all!');
    cfg.ext_all = false;
    cfg.inline_all = false;
end

if cfg.ext_all && cfg.inline_all
    warning('Both the ext_all or inline_all options are set. Ignoring inline_all!');
    cfg.inline_all = false;
end


%% copy the current figure to work on it instead of the original
plot_handles.figure = copyobj(h_figure,0);
plot_handles.axes = plot_handles.figure.CurrentAxes;


%% prepare figure and gather some information
% hide everything that is not an axes object
plot_handles.axes.Visible = 'off';
title(plot_handles.axes,'')

colorbar_limits = nan;
colorbar_label = '';
i=1;
while true
    if isa(plot_handles.figure.Children(i), 'matlab.graphics.illustration.ColorBar')
        colorbar_limits = plot_handles.figure.Children(i).Limits;
        colorbar_label = plot_handles.figure.Children(i).Label.String;
    end
    if ~isa(plot_handles.figure.Children(i), 'matlab.graphics.axis.Axes')
        delete(plot_handles.figure.Children(i))
        i = i-1;
    end
    if i == numel(plot_handles.figure.Children)
        break;
    end
    i = i+1;
end
    

if iscellstr(plot_handles.axes.XLabel.String)
    horz_label_txt = plot_handles.axes.XLabel.String{1};
else
    horz_label_txt = plot_handles.axes.XLabel.String;
end

if iscellstr(plot_handles.axes.YLabel.String)
    vert_label_txt = plot_handles.axes.YLabel.String{1};
else
    vert_label_txt = plot_handles.axes.YLabel.String;
end

if iscellstr(plot_handles.axes.ZLabel.String)
    zlabel_txt = plot_handles.axes.ZLabel.String{1};
else
    zlabel_txt = plot_handles.axes.ZLabel.String;
end


%% decide if the three dimensional approach or the easy two dimensional one is neccessary
current_view_point = plot_handles.axes.View;

if ~cfg.force_3d && ~sum(mod(current_view_point,90))
    plot2d = true;
else
    plot2d = false;
end


%% find axes limits and explicitly set them

[axes_x, axes_y, axes_z, v_data_x, v_data_y, v_data_z] = get_plot_limits(plot_handles.axes);

plot_handles.axes.XLim = axes_x;
plot_handles.axes.YLim = axes_y;
plot_handles.axes.ZLim = axes_z;

if isnan(colorbar_limits)
    colorbar_limits = v_data_z;
end


%% gather plot data
current_view_point = plot_handles.axes.View;
use_imagesc = false;

% determine graphic objects to plot vs. to print
children_idc_plot_ext = [];
children_idc_plot_inline = [];
children_idc_print = [];
for i=1:numel(plot_handles.axes.Children)
    if isa(plot_handles.axes.Children(i), 'matlab.graphics.chart.primitive.Line')
        if numel(plot_handles.axes.Children(i).XData) > cfg.inline_limit
            children_idc_plot_ext(end+1) = i;
        else
            children_idc_plot_inline(end+1) = i;
        end
    else
        children_idc_print(end+1) = i;
    end
end

% overwriting depending on config options
if cfg.print_all
    children_idc_print = sort([children_idc_print, children_idc_plot_ext, children_idc_plot_inline]);
    children_idc_plot_ext = [];
    children_idc_plot_inline = [];
end
if cfg.ext_all
    children_idc_plot_ext = sort([children_idc_plot_ext, children_idc_plot_inline]);
    children_idc_plot_inline = [];
end
if cfg.inline_all
    children_idc_plot_inline = sort([children_idc_plot_inline, children_idc_plot_ext]);
    children_idc_plot_ext = [];
end

% gather plot information
if plot2d
    % basic idea:
    % for those special views (0,0), (180,0), (90,0), (-90,0),  (0,90), and (-90,90) no 3D
    % approach is necessary. A PNG file is saved and one has to specify the
    % coordinates of the lower left and the upper right point. This is done
    % in the following.
    % For an easy work flow the PNG file is trimed (the whitespace is
    % removed), hence the axis / xlim options and the real x,y,z values
    % must be taken into account for a correct axis label estimation
    
    con_axes = [axes_x; axes_y; axes_z];
    con_v_data = [v_data_x; v_data_y; v_data_z];

    [horz, vert] = viewpoint_to_img_axes( current_view_point );
    
    img_data_horz = con_v_data(abs(horz),:)';
    img_data_vert = con_v_data(abs(vert),:)';
    
    img_axes_horz = con_axes(abs(horz),:);
    img_axes_vert = con_axes(abs(vert),:);
    
    img_x_reversed = (horz < 0);
    img_y_reversed = (vert < 0);
    
    % point positions doesn't make sense
    pt_point_positions = [];
    
    % this is just xmin, xmax, ymin, and ymax
    tikz_support_points = [img_data_horz, img_data_vert];
    
    labels_txt = {horz_label_txt, vert_label_txt, zlabel_txt};
    horz_label_txt = labels_txt{abs(horz)};
    vert_label_txt = labels_txt{abs(vert)};
    
    if cfg.use_imagesc && abs(horz) < 3 && abs(vert) < 3
        % we are top down and can use imagesc if the user wants it
        use_imagesc = true;
    end
    
    [print_data_range_horz, print_data_range_vert] = get_print_data_range(plot_handles.axes.Children(children_idc_print), abs([horz, vert]));
    
    png_range_horz(1) = max(img_data_horz(1), print_data_range_horz(1));
    png_range_horz(2) = min(img_data_horz(2), print_data_range_horz(2));
    png_range_vert(1) = max(img_data_vert(1), print_data_range_vert(1));
    png_range_vert(2) = min(img_data_vert(2), print_data_range_vert(2));
    
else
    % determine plot support points for TikZ
    tikz_support_points = get_tikz_support_points(axes_x, axes_y, axes_z, current_view_point, cfg);
    
    % draw support markers and determine paper position
    [pt_point_positions, px_pos_order] = get_point_positions(tikz_support_points, plot_handles, cfg, debug);
end

%% print image data
if (cfg.write_png)
    if use_imagesc
        % select the first surf plot
        for surf_idx=1:numel(plot_handles.axes.Children)
            if isa(plot_handles.axes.Children(surf_idx), 'matlab.graphics.chart.primitive.Surface')
                break;
            elseif surf_idx == numel(plot_handles.axes.Children)
                error('Could not find a surface plot to convert to image!')
            end
        end
        
        xdata = unique(plot_handles.axes.Children(surf_idx).XData);
        ydata = unique(plot_handles.axes.Children(surf_idx).YData);
        cdata = plot_handles.axes.Children(surf_idx).CData;
        
        if (abs(horz(1)) == 1)
            png_range_x = png_range_horz;
            png_range_y = png_range_vert;
        else
            png_range_x = png_range_vert;
            png_range_y = png_range_horz;
        end
        
        [~, min_x_idx] = min(abs(png_range_x(1) - xdata));
        min_x_idx = max(1,min_x_idx-1);
        png_range_x(1) = xdata(min_x_idx);
        
        [~, max_x_idx] = min(abs(png_range_x(2) - xdata));
        max_x_idx = min(size(cdata,2),max_x_idx+1);
        png_range_x(2) = xdata(max_x_idx);
        
        [~, min_y_idx] = min(abs(png_range_y(1) - ydata));
        min_y_idx = max(1,min_y_idx-1);
        png_range_y(1) = ydata(min_y_idx);
        
        [~, max_y_idx] = min(abs(png_range_y(2) - ydata));
        max_y_idx = min(size(cdata,1),max_y_idx+1);
        png_range_y(2) = ydata(max_y_idx);
        
        if (abs(horz(1)) == 1)
            png_range_horz = png_range_x;
            png_range_vert = png_range_y;
        else
            png_range_horz = png_range_y;
            png_range_vert = png_range_x;
        end
        
        cdata = cdata(min_y_idx:max_y_idx, min_x_idx:max_x_idx);
        
        
        
        if abs(horz) ~= 1
            cdata = cdata';
        end
        if horz < 0
            cdata = cdata(:,end:-1:1);
        end
        if vert < 0
            cdata = cdata(end:-1:1,:);
        end
        cmap = plot_handles.figure.Colormap;
        im_scaled = round((cdata(end:-1:1,:)-min(cdata(:)))./(max(cdata(:))-min(cdata(:)))*size(cmap,1));
        imwrite(im_scaled, cmap, [export_name, '.png'], 'png')
    else
        if ~cfg.print_all
            for i=[children_idc_plot_ext, children_idc_plot_inline]
                plot_handles.axes.Children(i).Visible = 'off';
            end
        end
        print(plot_handles.figure, export_name, '-dpng', ['-r' num2str(cfg.export_dpi)]);
        system(['mogrify -transparent white ', export_name, '.png']);
        if plot2d
            system(['mogrify -trim ', export_name, '.png']);
        end
    end
end


%% write TikZ/PGFPlots data

[~, export_fname, ~] = fileparts(export_name);

if plot2d
    view_dims = [horz, vert];
else
    view_dims = [];
end

[plot_parameters_ext, plot_data_filenames] = process_plots_ext(plot_handles.axes.Children(children_idc_plot_ext), export_name, plot2d, abs(view_dims));
[plot_parameters_inline, plot_points_inline] = process_plots_inline(plot_handles.axes.Children(children_idc_plot_inline), plot2d, abs(view_dims));

if (cfg.write_tikz)
    tfile_h = fopen([export_name, '.tikz'], 'w');
    fprintf(tfile_h, '%% created with surf3tikz written by Johannes Schlichenmaier\n');
    fprintf(tfile_h, '%% based on surf2TikZ written by Fabian Roos\n');
    fprintf(tfile_h, '\\begin{tikzpicture}\n');
    fprintf(tfile_h, '\t \\begin{axis}[\n');
    fprintf(tfile_h, '\t \t grid,\n');
    fprintf(tfile_h, '\t \t enlargelimits = false,\n');
    
    if plot2d
        fprintf(tfile_h, '\t \t xmin = %f,\n', img_axes_horz(1));
        fprintf(tfile_h, '\t \t xmax = %f,\n', img_axes_horz(2));
        fprintf(tfile_h, '\t \t ymin = %f,\n', img_axes_vert(1));
        fprintf(tfile_h, '\t \t ymax = %f,\n', img_axes_vert(2));
        if img_x_reversed
            fprintf(tfile_h, '\t \t x dir = reverse,\n');
        end
        if img_y_reversed
            fprintf(tfile_h, '\t \t y dir = reverse,\n');
        end
        fprintf(tfile_h, '\t \t xlabel = {%s},\n', horz_label_txt);
        fprintf(tfile_h, '\t \t ylabel = {%s},\n', vert_label_txt);
    else
        
        fprintf(tfile_h, '\t \t xmin = %f,\n', axes_x(1));
        fprintf(tfile_h, '\t \t xmax = %f,\n', axes_x(2));
        fprintf(tfile_h, '\t \t ymin = %f,\n', axes_y(1));
        fprintf(tfile_h, '\t \t ymax = %f,\n', axes_y(2));
        fprintf(tfile_h, '\t \t zmin = %f,\n', axes_z(1));
        fprintf(tfile_h, '\t \t zmax = %f,\n', axes_z(2));
        fprintf(tfile_h, '\t \t xlabel = {%s},\n', horz_label_txt);
        fprintf(tfile_h, '\t \t ylabel = {%s},\n', vert_label_txt);
        fprintf(tfile_h, '\t \t zlabel = {%s},\n', zlabel_txt);
    end
    
    if (current_view_point(1) == 90 && current_view_point(2) == -90) || (current_view_point(1) == -90 && current_view_point(2) == 90)
    else
    end
    
    fprintf(tfile_h, '\t \t colorbar,\n');
    if isequal(plot_handles.figure.Colormap, parula)
        fprintf(tfile_h, '\t \t colormap name=parula,\n');
    elseif isequal(plot_handles.figure.Colormap(1,:), [0.26700401, 0.00487433, 0.32941519]) && ...
            isequal(plot_handles.figure.Colormap(end,:), [0.99324789, 0.90615657, 0.1439362])
        fprintf(tfile_h, '\t \t colormap/viridis,\n');
    end
    fprintf(tfile_h, '\t \t colorbar style = {%%,\n');
    fprintf(tfile_h, '\t \t \t ylabel = {%s},\n', colorbar_label);
    fprintf(tfile_h, '\t \t },\n');
    fprintf(tfile_h, '\t \t point meta min = %f,\n', colorbar_limits(1));
    fprintf(tfile_h, '\t \t point meta max = %f,\n', colorbar_limits(2));
    fprintf(tfile_h, '\t ]\n');
    if plot2d
        fprintf(tfile_h, '\t \t \\addplot graphics\n');
        fprintf(tfile_h, '\t \t \t [xmin=%f, xmax=%f, ymin=%f, ymax=%f]\n', ...
            png_range_horz(1), png_range_horz(2), png_range_vert(1), png_range_vert(2) );
        fprintf(tfile_h, '\t \t {%s.png};\n', export_fname);
    else
        fprintf(tfile_h, '\t \t \\addplot3 graphics[\n');
        fprintf(tfile_h, '\t \t \t points={%% important\n');
        for i=px_pos_order
            fprintf(tfile_h, '\t \t \t \t (%f,%f,%f) => (%f,%f)\n', tikz_support_points(i,1), ...
                tikz_support_points(i,2), ...
                tikz_support_points(i,3), ...
                pt_point_positions(i,1), ...
                pt_point_positions(i,2));
        end
        fprintf(tfile_h, '\t \t }]{%s.png};\n', export_fname);
    end
    
    fprintf(tfile_h, '\n');
    
    for i=sort([children_idc_plot_ext, children_idc_plot_inline])
        [~,ext_idx] = ismember(i, children_idc_plot_ext);
        [~,inline_idx] = ismember(i, children_idc_plot_inline);
        
        if plot2d
            fprintf(tfile_h, '\t \t \\addplot+[%%\n');
        else
            fprintf(tfile_h, '\t \t \\addplot3+[%%\n');
        end
        fprintf(tfile_h, '\t \t \t %% %s\n', 'mark=*, % x, +, o');
        fprintf(tfile_h, '\t \t \t %% %s\n', '% only marks');
        if ext_idx
            fprintf(tfile_h, '\t \t ]\n');
            fprintf(tfile_h, '\t \t table[col sep = comma]{%s};\n', plot_data_filenames{ext_idx});
        else
            fprintf(tfile_h, '\t \t ] plot coordinates {%%\n');
            if plot2d
                fprintf(tfile_h, '\t \t \t (%f,%f)\n', plot_points_inline{inline_idx});
            else
                fprintf(tfile_h, '\t \t \t (%f,%f,%f)\n', plot_points_inline{inline_idx});
            end
            fprintf(tfile_h, '\t \t };\n');
        end
        fprintf(tfile_h, '\n');
    end
    
    fprintf(tfile_h, '\t \\end{axis}\n');
    fprintf(tfile_h, '\\end{tikzpicture}');
    fclose(tfile_h);
end

close(plot_handles.figure)

%% save the original figure if desired
if cfg.write_fig
    savefig(h_figure, [export_name, '.fig'])
end
end

%% additional functions

function [ global_data_limit_x, global_data_limit_y, global_data_limit_z ] = get_data_limits( axes )
%GET_DATA_LIMITS get the max and min values of the data plotted in specific axes
%   This is not neccessarily the same as the plot ranges.

% [Inf, -Inf] is neccessary so e.g, lower data trumps the current min, starting from Inf
global_data_limit_x = [Inf, -Inf];
global_data_limit_y = [Inf, -Inf];
global_data_limit_z = [Inf, -Inf];

for i=1:numel(axes.Children)
    data_limits_x(1) = min(axes.Children(i).XData(:));
    data_limits_x(2) = max(axes.Children(i).XData(:));
    data_limits_y(1) = min(axes.Children(i).YData(:));
    data_limits_y(2) = max(axes.Children(i).YData(:));
    data_limits_z(1) = min(axes.Children(i).ZData(:));
    data_limits_z(2) = max(axes.Children(i).ZData(:));
    
    if data_limits_x(1) < global_data_limit_x(1)
        global_data_limit_x(1) = data_limits_x(1);
    end
    if data_limits_x(2) > global_data_limit_x(2)
        global_data_limit_x(2) = data_limits_x(2);
    end
    
    if data_limits_y(1) < global_data_limit_y(1)
        global_data_limit_y(1) = data_limits_y(1);
    end
    if data_limits_y(2) > global_data_limit_y(2)
        global_data_limit_y(2) = data_limits_y(2);
    end
    
    if data_limits_z(1) < global_data_limit_z(1)
        global_data_limit_z(1) = data_limits_z(1);
    end
    if data_limits_z(2) > global_data_limit_z(2)
        global_data_limit_z(2) = data_limits_z(2);
    end
end

end


function [ axes_x, axes_y, axes_z, v_data_x, v_data_y, v_data_z ] = get_plot_limits( axes )
%GET_PLOT_LIMITS returns the effective axes limits, as well as the limits of the visible data

axes_x = axes.XLim;
axes_y = axes.YLim;
axes_z = axes.ZLim;

[data_limit_x, data_limit_y, data_limit_z] = get_data_limits(axes);
%TODO: what to do when data is Inf and so are axes?

% use auto mode value for axes where neccessary
if isinf(axes_x(1))
    axes_x(1) = data_limit_x(1);
end
if isinf(axes_x(2))
    axes_x(2) = data_limit_x(2);
end

if isinf(axes_y(1))
    axes_y(1) = data_limit_y(1);
end
if isinf(axes_y(2))
    axes_y(2) = data_limit_y(2);
end

if isinf(axes_z(1))
    axes_z(1) = data_limit_z(1);
end
if isinf(axes_z(2))
    axes_z(2) = data_limit_z(2);
end

v_data_x = [max(axes_x(1), data_limit_x(1)), min(axes_x(2), data_limit_x(2))];
v_data_y = [max(axes_y(1), data_limit_y(1)), min(axes_y(2), data_limit_y(2))];
v_data_z = [max(axes_z(1), data_limit_z(1)), min(axes_z(2), data_limit_z(2))];

end


function [ horz, vert ] = viewpoint_to_img_axes( viewpoint )
%VIEPOINT_TO_IMG_AXES converts a figure viewport to a 2D image (X,Y) CS
% returns the horizontal and vertical dimension (x=1, y=2, z=3)

viewpoint = wrapTo180(viewpoint);

az = viewpoint(1);
el = viewpoint(2);

x = 1;
y = 2;
z = 3;

if sum(mod(viewpoint,90)) 
    error('No process yet for 2D and Az=%d, El=%d', current_view_point);
end

if abs(az) == abs(el) && abs(el) == 90
    horz = sign(az)*y;
    vert = -1*sign(az)*sign(el)*x;
else
    if abs(az) == 90
        horz = sign(az)*y;
    elseif abs(az) == 180
        horz = -x;
    else
        horz = x;
    end
    if abs(el) == 90
        vert = sign(horz)*sign(el)*y;
    elseif abs(el) == 180
        vert = -1*z;
    else
        vert = z;
    end
end

end


function [ support_points ] = get_tikz_support_points(limits_x, limits_y, limits_z, view_point, cfg)
%GET_TIKZ_SUPPORT_POINTS determines the support points for the PGFPlot image
% determine axes box outer points
box_points = nan(8,3);
box_points(1,:) = [limits_x(1),limits_y(1),limits_z(1)];
box_points(2,:) = [limits_x(1),limits_y(1),limits_z(2)];
box_points(3,:) = [limits_x(1),limits_y(2),limits_z(1)];
box_points(4,:) = [limits_x(1),limits_y(2),limits_z(2)];
box_points(5,:) = [limits_x(2),limits_y(1),limits_z(1)];
box_points(6,:) = [limits_x(2),limits_y(1),limits_z(2)];
box_points(7,:) = [limits_x(2),limits_y(2),limits_z(1)];
box_points(8,:) = [limits_x(2),limits_y(2),limits_z(2)];

% rotate box points and find out which points lie on top of each other

rot_matrix_Z = rotz(-view_point(1));
rot_matrix_X = rotx(view_point(2));
box_points_turned = (rot_matrix_X*(rot_matrix_Z*box_points'))';
[~, ~, ubp_rev_idc] = unique([box_points_turned(:,1), box_points_turned(:,3)], 'rows');

% As far as my testing showed, PGFPlots wants to have a set of points where you can find at least
% a pair of points, that stay the same in one dim and changes in the other two, for every dim.
% These pairs are precalculated. But to remind me and inform interested people, I'll leave this
% code in here anyways.
% drawer = nan(1,3,8);
% drawer(1,:,:) = box_points';
% box_matrix = repmat(box_points,1,1,8);
% sub_matrix = repmat(drawer, 8,1,1);
%
% diff_matrix = box_matrix - sub_matrix;
%
% X_pairs = squeeze(diff_matrix(:,1,:) == 0 & diff_matrix(:,2,:) ~= 0 & diff_matrix(:,3,:) ~= 0);
% Y_pairs = squeeze(diff_matrix(:,2,:) == 0 & diff_matrix(:,3,:) ~= 0 & diff_matrix(:,1,:) ~= 0);
% Z_pairs = squeeze(diff_matrix(:,3,:) == 0 & diff_matrix(:,1,:) ~= 0 & diff_matrix(:,2,:) ~= 0);
%
% [row, col] = find(X_pairs);
% X_pairs = unique(sort([row, col], 2), 'rows')
% [row, col] = find(Y_pairs);
% Y_pairs = unique(sort([row, col], 2), 'rows')
% [row, col] = find(Z_pairs);
% Z_pairs = unique(sort([row, col], 2), 'rows')
%
% X_pairs = [1,4;2,3;5,8;6,7];
% Y_pairs = [1,6;2,5;3,8;4,7];
% Z_pairs = [1,7;2,8;3,5;4,6];

% If you don't want to fulfill the aforementioned constraints with only three points and add
% an additional point, there are only two sets of points that work (I have the feeling that
% here is a high potential for bugs...):
box_point_sets = [1,4,6,7;2,3,5,8];

if ~isempty(cfg.box_point_idc);
    support_points = box_points(cfg.box_point_idc,:);
else
    if numel(unique(ubp_rev_idc(box_point_sets(1,:)))) == 4
        support_points = box_points(box_point_sets(1,:),:);
    elseif numel(unique(ubp_rev_idc(box_point_sets(2,:)))) == 4
        support_points = box_points(box_point_sets(2,:),:);
    else
        error('Neither of the proposed point trains seem to work, could not find 4 viable points...');
    end
end
end


function [ pt_point_positions, px_pos_order ] = get_point_positions( tikz_support_points, plot_handles, cfg, debug )
%GET_POINT_POSITIONS determines the positions of the support plots in points
% first, hide all other plots
for i=1:numel(plot_handles.axes.Children)
    plot_handles.axes.Children(i).Visible = 'off';
end

% prepare marker handle for paper position determination
hold(plot_handles.axes, 'on');
plot_handles.support_plot = plot3(plot_handles.axes, nan, nan, nan, ...
    'Marker', 'square', ...
    'MarkerFaceColor', [0,0,0], ...
    'MarkerEdgeColor', [1,1,1], ...
    'LineStyle', 'none', ...
    'Tag', 'BoxPoints', ...
    'MarkerSize', 10);
hold(plot_handles.axes, 'off');

% determining the paper positions of the selected box points
mid_row_pixel = zeros(4,1);
mid_col_pixel = zeros(4,1);
for i = 1:4
    % setting cursor position
    plot_handles.support_plot.XData = tikz_support_points(i,1);
    plot_handles.support_plot.YData = tikz_support_points(i,2);
    plot_handles.support_plot.ZData = tikz_support_points(i,3);
    % create 2D image
    frame_data = getframe(plot_handles.figure);
    [image_data,~] = frame2im(frame_data);
    flat_image_data = image_data(:,:,1)+image_data(:,:,2)+image_data(:,:,3);
    
    % find black rectangle
    min_idc = find(flat_image_data(:) == 0);
    [min_row, min_col] = ind2sub(size(flat_image_data), min_idc);
    
    % extract center point of black rectangle
    unique_min_row = unique(min_row);
    unique_min_col = unique(min_col);
    mid_row_pixel(i) = round((max(unique_min_row) + min(unique_min_row))/2);
    mid_col_pixel(i) = round((max(unique_min_col) + min(unique_min_col))/2);
end

% remap to origin down left
[ysize, ~] = size(flat_image_data);
pgf_pixel_points = [mid_col_pixel, repmat(ysize,4,1)-mid_row_pixel];

% calculate pts values
pt_point_positions = pgf_pixel_points./cfg.screen_ppi./cfg.inch2point_ratio;
if ~cfg.force_exact_values
    pt_point_positions = round(pt_point_positions);
end

% Apparently PGFPlot is a bit picky about the first two points. It apparently works best if both
% pixel positions are different in the first two points. So let's look for a solution.
drawer = nan(1,2,4);
drawer(1,:,:) = pt_point_positions';
box_matrix = repmat(pt_point_positions,1,1,4);
sub_matrix = repmat(drawer,4,1,1);
diff_matrix = box_matrix - sub_matrix;
possible_starts = squeeze(diff_matrix(:,1,:) ~= 0 & diff_matrix(:,2,:) ~= 0);

px_pos_order = [0,0,0,0];
for i=1:4
    poss_next_point = find(possible_starts(:,1));
    if ~isempty(poss_next_point)
        px_pos_order(1) = i;
        px_pos_order(2) = poss_next_point(1);
        px_pos_order(3:4) = setdiff(1:4, px_pos_order);
        break;
    elseif i==4
        warning('Could not find a good pixel sorting. Finishing up regardless!');
        px_pos_order = 1:4;
    end
end

for i=1:numel(plot_handles.axes.Children)
    plot_handles.axes.Children(i).Visible = 'on';
end
plot_handles.support_plot.Visible = 'off';

if debug
    plot_handles.cursor_manager.UpdateFcn = [];
    
    plot_handles.support_plot.XData = tikz_support_points(:,1);
    plot_handles.support_plot.YData = tikz_support_points(:,2);
    plot_handles.support_plot.ZData = tikz_support_points(:,3);
    plot_handles.support_plot.Visible = 'on';
    
    plot_handles.cursor_manager = datacursormode(plot_handles.figure);
    for i = 1:4
        debug_cursor = createDatatip(plot_handles.cursor_manager, plot_handles.support_plot);
        debug_cursor.Position = tikz_support_points(i,:);
    end
    
    debug_frame_data = getframe(plot_handles.figure);
    [debug_image_data,~] = frame2im(debug_frame_data);
    
    plot_handles.cursor_manager.removeAllDataCursors
    
    for i = 1:4
        debug_image_data(mid_row_pixel(i), mid_col_pixel(i), :) = 255*[1,1,1];
    end
    
    figure
    image(debug_image_data);
    figure(plot_handles.figure)
end

delete(plot_handles.support_plot);

end


function [plot_parameters, file_names ] = process_plots_ext(g_objects, write_name, plot2d, view_dims)
%PROCESS_PLOTS saves line plots into CSV files
num_objects = numel(g_objects);
[file_path, file_base_name, ~] = fileparts(write_name);

if isempty(file_path)
    file_path = '.';
end

file_names = cell(num_objects);
plot_parameters = cell(num_objects);

for i=1:num_objects
    XYZ_Data = [g_objects(i).XData;g_objects(i).YData;g_objects(i).ZData];
    file_names{i} = sprintf([file_base_name, '-plot-%d.csv'], i);
    f_handle = fopen([file_path filesep file_names{i}], 'w+');
    if plot2d
        XY_Data = XYZ_Data(view_dims, :);
        fprintf(f_handle, '%f,%f\n', XY_Data);
    else
        fprintf(f_handle, '%f,%f,%f\n', XYZ_Data);
    end
    fclose(f_handle);
end

end


function [plot_parameters, plot_points ] = process_plots_inline(g_objects, plot2d, view_dims)
%PROCESS_PLOTS saves line plots into CSV files
num_objects = numel(g_objects);

plot_points = cell(num_objects);
plot_parameters = cell(num_objects);

for i=1:num_objects
    XYZ_Data = [g_objects(i).XData;g_objects(i).YData;g_objects(i).ZData];
    if plot2d
        XY_Data = XYZ_Data(view_dims, :);
        plot_points{i} = XY_Data; % in Tikz_File sprintf('%f,%f\n', XY_Data);
    else
        plot_points{i} = XYZ_Data; % in Tikz_File sprintf('%f,%f,%f\n', XYZ_Data);
    end
end

end


function [ range_horz, range_vert] = get_print_data_range( g_objects, view_dims )
%GET_PRINT_DATA_RANGE returns the horizontal and vertical limits of the printed image data (2D
%approach only)

range_horz = [Inf, -Inf];
range_vert = [Inf, -Inf];

for i=1:numel(g_objects)
    XYZ_Data_c = {g_objects(i).XData(:)';g_objects(i).YData(:)';g_objects(i).ZData(:)'};
    XY_Data_c = XYZ_Data_c(view_dims);
    data_limits_horz(1) = min(XY_Data_c{1});
    data_limits_horz(2) = max(XY_Data_c{1});
    data_limits_vert(1) = min(XY_Data_c{2});
    data_limits_vert(2) = max(XY_Data_c{2});
    
    if data_limits_horz(1) < range_horz(1)
        range_horz(1) = data_limits_horz(1);
    end
    if data_limits_horz(2) > range_horz(2)
        range_horz(2) = data_limits_horz(2);
    end
    
    if data_limits_vert(1) < range_vert(1)
        range_vert(1) = data_limits_vert(1);
    end
    if data_limits_vert(2) > range_vert(2)
        range_vert(2) = data_limits_vert(2);
    end
end

end