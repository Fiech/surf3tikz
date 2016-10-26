function [pt_point_positions, tikz_support_points, colorbar_limits] = surf3tikz(h_figure, export_name, cfg, debug)
%SURF3TIKZ converts a surf plot into a tikz file
% This function takes a figure containing a surf plot and creates a tikz file and an according png
% file. To achieve this, it automatically places data markers on points of the axes box around
% the plot and subsequently detects their pixel position.
% For views (0,0), (180,0), (90,0), and (-90,0) a simpler approach is used which creates
% a 2D plot in TikZ, instead of a 3D one. This is done to prevent unstable behaviour because of
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
%      .force_3d: For special views (0,0), (180,0), (90,0), and (-90,0) instead of the 3D approach a
%       2D approach is used. To force the use of the 3D approach set this to true, default: false
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

if nargin < 4
    debug = false;
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

xlabel_txt = plot_handles.axes.XLabel.String;
ylabel_txt = plot_handles.axes.YLabel.String;


%% find axes limits and explicitly set them
axes_limit_x = plot_handles.axes.XLim;
axes_limit_y = plot_handles.axes.YLim;
axes_limit_z = plot_handles.axes.ZLim;

% find limits of all the plot objects of the current axes

[data_limits_x, data_limits_y, data_limits_z] = get_data_limits(plot_handles.axes);

for i=1:numel(plot_handles.axes.Children)
    plot_handles.axes.Children.Visible = 'off';
end

% use auto mode value for axes where neccessary
if isinf(axes_limit_x(1))
    axes_limit_x(1) = data_limits_x(1);
    plot_handles.axes.XLim = axes_limit_x;
end
if isinf(axes_limit_x(2))
    axes_limit_x(2) = data_limits_x(2);
    plot_handles.axes.XLim = axes_limit_x;
end

if isinf(axes_limit_y(1))
    axes_limit_y(1) = data_limits_y(1);
    plot_handles.axes.YLim = axes_limit_y;
end
if isinf(axes_limit_y(2))
    axes_limit_y(2) = data_limits_y(2);
    plot_handles.axes.YLim = axes_limit_y;
end

if isinf(axes_limit_z(1))
    axes_limit_z(1) = data_limits_z(1);
    plot_handles.axes.ZLim = axes_limit_z;
end
if isinf(axes_limit_z(2))
    axes_limit_z(2) = data_limits_z(2);
    plot_handles.axes.ZLim = axes_limit_z;
end

plot_handles.axes.XLim = axes_limit_x;
plot_handles.axes.YLim = axes_limit_y;
plot_handles.axes.ZLim = axes_limit_z;


%% decide if the three dimensional approach or the easy two dimensional one is neccessary
current_view_point = plot_handles.axes.View;

if ~cfg.force_3d && ~sum(mod(current_view_point,90))
    plot2d = true;
    % basic idea:
    % for those special views (0,0), (180,0), (90,0), and (-90,0) no 3D
    % approach is necessary. A PNG file is saved and one has to specify the
    % coordinates of the lower left and the upper right point. This is done
    % in the following.
    % For an easy work flow the PNG file is trimed (the whitespace is
    % removed), hence the axis / xlim options and the real x,y,z values
    % must be taken into account for a correct axis label estimation
    if (current_view_point(1) == 0 || current_view_point(1) == 180) && current_view_point(2) == 0
        % x-z view selected
        img_x(1) = max(axes_limit_x(1), data_limits_x(1));
        img_x(2) = min(axes_limit_x(2), data_limits_x(2));
        img_x(3:4) = axes_limit_x;
        img_y(1) = max(axes_limit_z(1), data_limits_z(1));
        img_y(2) = min(axes_limit_z(2), data_limits_z(2));
        img_y(3:4) = axes_limit_z;
        img_rev = false;
        if current_view_point(1) == 180
            img_rev = true;
        end
    elseif (current_view_point(1) == 90 || current_view_point(1) == -90) && current_view_point(2) == 0
        % y-z view selected
        img_x(1) = max(axes_limit_y(1), data_limits_y(1));
        img_x(2) = min(axes_limit_y(2), data_limits_y(2));
        img_x(3:4) = axes_limit_y;
        img_y(1) = max(axes_limit_z(1), data_limits_z(1));
        img_y(2) = min(axes_limit_z(2), data_limits_z(2));
        img_y(3:4) = axes_limit_z;
        img_rev = false;
        if current_view_point(1) == -90
            img_rev = true;
        end
    else
        % only a quick-and-dirty solution, solve this with imagesc way
        if current_view_point(1) == 0 && current_view_point(2) == 90
            img_x(1) = max(axes_limit_x(1), data_limits_x(1));
            img_x(2) = min(axes_limit_x(2), data_limits_x(2));
            img_x(3:4) = axes_limit_x;
            img_y(1) = max(axes_limit_y(1), data_limits_y(1));
            img_y(2) = min(axes_limit_y(2), data_limits_y(2));
            img_y(3:4) = axes_limit_y;
            img_rev = false;
        elseif (current_view_point(1) == 90 && current_view_point(2) == -90) || (current_view_point(1) == -90 && current_view_point(2) == 90)
            % y-x view selected and -y-x view selected (reverse axis)
            img_x(1) = max(axes_limit_y(1), data_limits_y(1));
            img_x(2) = min(axes_limit_y(2), data_limits_y(2));
            img_x(3:4) = axes_limit_y;
            img_y(1) = max(axes_limit_x(1), data_limits_x(1));
            img_y(2) = min(axes_limit_x(2), data_limits_x(2));
            img_y(3:4) = axes_limit_x;
            
            if current_view_point(1) == 90 && current_view_point(2) == -90
                img_rev = false;
            else
                % -y-x view selected
                img_rev = true;
            end
        else
            error('imagesc not yet implemented!')
        end
    end
    
    % point positions doesn't make sense
    pt_point_positions = [];
    
    % this is just xmin, xmax, ymin, and ymax
    tikz_support_points = [img_x, img_y];
else
    plot2d = false;
    
    
    %% choosing tikz support points
    % determine axes box outer points
    box_points = nan(8,3);
    box_points(1,:) = [axes_limit_x(1),axes_limit_y(1),axes_limit_z(1)];
    box_points(2,:) = [axes_limit_x(1),axes_limit_y(1),axes_limit_z(2)];
    box_points(3,:) = [axes_limit_x(1),axes_limit_y(2),axes_limit_z(1)];
    box_points(4,:) = [axes_limit_x(1),axes_limit_y(2),axes_limit_z(2)];
    box_points(5,:) = [axes_limit_x(2),axes_limit_y(1),axes_limit_z(1)];
    box_points(6,:) = [axes_limit_x(2),axes_limit_y(1),axes_limit_z(2)];
    box_points(7,:) = [axes_limit_x(2),axes_limit_y(2),axes_limit_z(1)];
    box_points(8,:) = [axes_limit_x(2),axes_limit_y(2),axes_limit_z(2)];
    
    % rotate box points and find out which points lie on top of each other
    
    rot_matrix_Z = rotz(-current_view_point(1));
    rot_matrix_X = rotx(current_view_point(2));
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
        tikz_support_points = box_points(cfg.box_point_idc,:);
    else
        if numel(unique(ubp_rev_idc(box_point_sets(1,:)))) == 4
            tikz_support_points = box_points(box_point_sets(1,:),:);
        elseif numel(unique(ubp_rev_idc(box_point_sets(2,:)))) == 4
            tikz_support_points = box_points(box_point_sets(2,:),:);
        else
            error('Neither of the proposed point trains seem to work, could not find 4 viable points...');
        end
    end
    
    
    %% draw support markers and determine paper position
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
    
end

%% tidy up
for i=1:numel(plot_handles.axes.Children)
    plot_handles.axes.Children(i).Visible = 'on';
end
plot_handles.support_plot.Visible = 'off';

if ~cfg.force_3d && ~sum(mod(current_view_point,90))
    if debug
        plot_handles.cursor_manager.UpdateFcn = [];
        
        plot_handles.support_plot.XData = box_points(:,1);
        plot_handles.support_plot.YData = box_points(:,2);
        plot_handles.support_plot.ZData = box_points(:,3);
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
        plot_handles.support_plot.Visible = 'off';
    end
    
end

%% write output files

% print png and make transparent
if (cfg.write_png)
    print(plot_handles.figure, export_name, '-dpng', ['-r' num2str(cfg.export_dpi)]);
    system(['mogrify -transparent white ', export_name, '.png']);
    if plot2d
        system(['mogrify -trim ', export_name, '.png']);
    end
end

if isnan(colorbar_limits)
    colorbar_limits = data_limits_z;
end

[~, export_fname, ~] = fileparts(export_name);

% write to TikZ file
if (cfg.write_tikz)
    tfile_h = fopen([export_name, '.tikz'], 'w');
    fprintf(tfile_h, '%% created with surf3tikz written by Johannes Schlichenmaier\n');
    fprintf(tfile_h, '%% based on surf2TikZ written by Fabian Roos\n');
    fprintf(tfile_h, '\\begin{tikzpicture}\n');
    fprintf(tfile_h, '\t \\begin{axis}[\n');
    fprintf(tfile_h, '\t \t grid,\n');
    fprintf(tfile_h, '\t \t enlargelimits = false,\n');
    
    if plot2d
        fprintf(tfile_h, '\t \t xmin = %f,\n', img_x(3));
        fprintf(tfile_h, '\t \t xmax = %f,\n', img_x(4));
        fprintf(tfile_h, '\t \t ymin = %f,\n', img_y(3));
        fprintf(tfile_h, '\t \t ymax = %f,\n', img_y(4));
        if img_rev
            fprintf(tfile_h, '\t \t x dir = reverse,\n');
        end
    else
        
        fprintf(tfile_h, '\t \t xmin = %f,\n', axes_limit_x(1));
        fprintf(tfile_h, '\t \t xmax = %f,\n', axes_limit_x(2));
        fprintf(tfile_h, '\t \t ymin = %f,\n', axes_limit_y(1));
        fprintf(tfile_h, '\t \t ymax = %f,\n', axes_limit_y(2));
        fprintf(tfile_h, '\t \t zmin = %f,\n', axes_limit_z(1));
        fprintf(tfile_h, '\t \t zmax = %f,\n', axes_limit_z(2));
    end
    
    if (current_view_point(1) == 90 && current_view_point(2) == -90) || (current_view_point(1) == -90 && current_view_point(2) == 90)
        fprintf(tfile_h, '\t \t xlabel = {%s},\n', ylabel_txt);
        fprintf(tfile_h, '\t \t ylabel = {%s},\n', xlabel_txt);
    else
        fprintf(tfile_h, '\t \t xlabel = {%s},\n', xlabel_txt);
        fprintf(tfile_h, '\t \t ylabel = {%s},\n', ylabel_txt);
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
            img_x(1), img_x(2), img_y(1), img_y(2) );
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
    
    fprintf(tfile_h, '\t \\end{axis}\n');
    fprintf(tfile_h, '\\end{tikzpicture}');
    fclose(tfile_h);
end

close(plot_handles.figure)

% optional, save the original figure if desired
if cfg.write_fig
    savefig(h_figure, [export_name, '.fig'])
end
end

function [ global_data_limits_x, global_data_limits_y, global_data_limits_z ] = get_data_limits( axes )
%GET_DATA_LIMITS get the max and min values of the data plotted in specific axes
%   This is not neccessarily the same as the plot ranges.

global_data_limits_x = [Inf, -Inf];
global_data_limits_y = [Inf, -Inf];
global_data_limits_z = [Inf, -Inf];

for i=1:numel(axes.Children)
    data_limits_x(1) = min(axes.Children(i).XData(:));
    data_limits_x(2) = max(axes.Children(i).XData(:));
    data_limits_y(1) = min(axes.Children(i).YData(:));
    data_limits_y(2) = max(axes.Children(i).YData(:));
    data_limits_z(1) = min(axes.Children(i).ZData(:));
    data_limits_z(2) = max(axes.Children(i).ZData(:));
    
    if data_limits_x(1) < global_data_limits_x(1)
        global_data_limits_x(1) = data_limits_x(1);
    end
    if data_limits_x(2) > global_data_limits_x(2)
        global_data_limits_x(2) = data_limits_x(2);
    end
    
    if data_limits_y(1) < global_data_limits_y(1)
        global_data_limits_y(1) = data_limits_y(1);
    end
    if data_limits_y(2) > global_data_limits_y(2)
        global_data_limits_y(2) = data_limits_y(2);
    end
    
    if data_limits_z(1) < global_data_limits_z(1)
        global_data_limits_z(1) = data_limits_z(1);
    end
    if data_limits_z(2) > global_data_limits_z(2)
        global_data_limits_z(2) = data_limits_z(2);
    end
end

end