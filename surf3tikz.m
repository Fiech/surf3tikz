function [pt_point_positions, tikz_support_points, colorbar_limits] = surf3tikz(h_figure, export_name, cfg, debug)
%SURF3TIKZ converts a surf plot into a tikz file
% This function takes a figure containing a surf plot and creates a tikz file and an according png
% file. To achieve this, it automatically places data cursors on extreme points of the box around the plot and
% subsequently detects their pixel position.
% Additionally the key values are available as output parameters to build your own tikz file from scratch.
%
% INPUT:
%   h_figure: Figure handle containing the figure plot. Ideally also contains a colorbar.
%   export_name: Base filename of the exported tikz and png file the function writes to
%   cfg: optional config struct
%      .export_dpi: dpi of exported png file, default: 300 dpi
%      .write_png: boolean to optionally suppress png output, default: true
%      .write_tikz: boolean to optionally suppress tikz file output, default: true
%      .screen_ppi: pixel per inch screen resolution, default: will try to determine system setting (normally 96)
%      .inch2point_ratio: ratio of inch to point, default: 1/72 inch per point
%   debug: boolean, will switch on/off an additional debug plot with the set cursors and the found pixel positions
%          denoted by a white pixel ideally in the middle of every black cursor rectangle, default: false
%
% OUTPUT:
%   pt_point_positions: the data point positions in the pt unit
%   act_data_points: the data point coordinates
%   colorbar_limits: cdata limits of the colorbar currently used
%

if nargin < 2 
	error('You''ll need to specify a figure handle and an export base filename');
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

if ~isfield(cfg, 'screen_ppi')
	cfg.screen_ppi = get(0,'ScreenPixelsPerInch');
end

if ~isfield(cfg, 'inch2point_ratio')
	cfg.inch2point_ratio = 1/72;
end

if nargin < 4
	debug = false;
end


plot_handles.figure = copyobj(h_figure,0);
plot_handles.axes = plot_handles.figure.CurrentAxes;

colorbar_limits = nan;
colorbar_label = '';
xlabel_txt = plot_handles.axes.XLabel.String;
ylabel_txt = plot_handles.axes.YLabel.String;

% remove everything that is not an axes object
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

global_data_limits_x = [Inf, -Inf];
global_data_limits_y = [Inf, -Inf];
global_data_limits_z = [Inf, -Inf];

% USERBREAK
plot_handles.axes.Visible = 'off';
title(plot_handles.axes,'')
% keyboard

% find limits of all the plot objects of the current axes
axes_limit_x = plot_handles.axes.XLim;
axes_limit_y = plot_handles.axes.YLim;
axes_limit_z = plot_handles.axes.ZLim;

for i=1:numel(plot_handles.axes.Children)
	data_limits_x(1) = min(plot_handles.axes.Children(i).XData(:));
	data_limits_x(2) = max(plot_handles.axes.Children(i).XData(:));
	data_limits_y(1) = min(plot_handles.axes.Children(i).YData(:));
	data_limits_y(2) = max(plot_handles.axes.Children(i).YData(:));
	data_limits_z(1) = min(plot_handles.axes.Children(i).ZData(:));
	data_limits_z(2) = max(plot_handles.axes.Children(i).ZData(:));
	
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
	
	plot_handles.axes.Children(i).Visible = 'off';
end

current_view_point = plot_handles.axes.View;

% use auto mode value for axes where neccessary
if isinf(axes_limit_x(1))
	axes_limit_x(1) = global_data_limits_x(1);
	plot_handles.axes.XLim = axes_limit_x;
end
if isinf(axes_limit_x(2))
	axes_limit_x(2) = global_data_limits_x(2);
	plot_handles.axes.XLim = axes_limit_x;
end

if isinf(axes_limit_y(1))
	axes_limit_y(1) = global_data_limits_y(1);
	plot_handles.axes.YLim = axes_limit_y;
end
if isinf(axes_limit_y(2))
	axes_limit_y(2) = global_data_limits_y(2);
	plot_handles.axes.YLim = axes_limit_y;
end

if isinf(axes_limit_z(1))
	axes_limit_z(1) = global_data_limits_z(1);
	plot_handles.axes.ZLim = axes_limit_z;
end
if isinf(axes_limit_z(2))
	axes_limit_z(2) = global_data_limits_z(2);
	plot_handles.axes.ZLim = axes_limit_z;
end

plot_handles.axes.XLim = axes_limit_x;
plot_handles.axes.YLim = axes_limit_y;
plot_handles.axes.ZLim = axes_limit_z;

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

rot_matrix_Z = rotz(-current_view_point(1));
rot_matrix_X = rotx(current_view_point(2));
box_points_turned = (rot_matrix_X*(rot_matrix_Z*box_points'))';

[~, ubpidc, ~] = unique([box_points_turned(:,1), box_points_turned(:,3)], 'rows');

box_points_X_min_idc = ubpidc(box_points_turned(ubpidc,1) == min(box_points_turned(ubpidc,1)));
box_points_X_max_idc = ubpidc(box_points_turned(ubpidc,1) == max(box_points_turned(ubpidc,1)));
box_points_Z_min_idc = ubpidc(box_points_turned(ubpidc,3) == min(box_points_turned(ubpidc,3)));
box_points_Z_max_idc = ubpidc(box_points_turned(ubpidc,3) == max(box_points_turned(ubpidc,3)));

box_points_idc = {box_points_X_min_idc, box_points_X_max_idc, box_points_Z_min_idc, box_points_Z_max_idc};

numel_X_min = numel(box_points_X_min_idc);
numel_X_max = numel(box_points_X_max_idc);
numel_Z_min = numel(box_points_Z_min_idc);
numel_Z_max = numel(box_points_Z_max_idc);

numels = [numel_X_min, numel_X_max, numel_Z_min, numel_Z_max];

tikz_support_points = nan(4,3);

if min(numels) == 1
	% first scenario: azimuth and elevation rotation
	% easy pick, no set dependecy, easy peasy: just take random (first) element from every set
	tikz_support_points(1,:) = box_points(box_points_idc{1}(1),:);
	tikz_support_points(2,:) = box_points(box_points_idc{2}(1),:);
	tikz_support_points(3,:) = box_points(box_points_idc{3}(1),:);
	tikz_support_points(4,:) = box_points(box_points_idc{4}(1),:);
elseif max(numels) == 4
	% second scenario: only azimuth rotation, but no coincided y axis
	% a bit more complex picking, semi dependend sets: start picking from sets in
	% ascending order of element number and don't use elements twice
	used_idc = [0,0,0,0];
	[~, sort_ind] = sort(numels);
	for i=1:4
		remaining_idc = setdiff(box_points_idc{sort_ind(i)},used_idc);
		tikz_support_points(i,:) = box_points(remaining_idc(1),:);
		used_idc(i) = remaining_idc(1);
	end	
elseif min(numels) == max(numels)
	% third scenario: coinciding x and y axis
	% complex picking, fully dependend sets: pick in every stage determines the net sext to pick from
	used_idc = [0,0,0,0];
	curr_set = 1;
	for i=1:4
		remaining_idc = setdiff(box_points_idc{curr_set},used_idc);
		tikz_support_points(i,:) = box_points(remaining_idc(1),:);
		used_idc(i) = remaining_idc(1);
		curr_set = setdiff(find(cellfun(@(x) sum(x == remaining_idc(1)), box_points_idc)),curr_set);
	end
end


hold(plot_handles.axes, 'on');
plot_handles.support_plot = plot3(plot_handles.axes, nan, nan, nan, ...
	'Marker', 'square', ...
	'MarkerFaceColor', [0,0,0], ...
	'MarkerEdgeColor', [1,1,1], ...
	'LineStyle', 'none', ...
	'Tag', 'BoxPoints', ...
	'MarkerSize', 10);
hold(plot_handles.axes, 'off');

mid_row_pixel = zeros(10,1);
mid_col_pixel = zeros(10,1);

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

%tidy up
for i=1:numel(plot_handles.axes.Children)	
	plot_handles.axes.Children(i).Visible = 'on';
end
plot_handles.support_plot.Visible = 'off';

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
end	

% remap to origin down left
[ysize, ~] = size(flat_image_data);
pgf_pixel_points = [mid_col_pixel, repmat(ysize,10,1)-mid_row_pixel];

% calculate pts values
pt_point_positions = pgf_pixel_points./cfg.screen_ppi./cfg.inch2point_ratio;

% print png and make transparent
if (cfg.write_png)
	print(plot_handles.figure, export_name, '-dpng', ['-r' num2str(cfg.export_dpi)]);
	system(['mogrify -transparent white ', export_name, '.png']);
end

if isnan(colorbar_limits)
	colorbar_limits = global_data_limits_z;
end

[~, export_fname, ~] = fileparts('../../LaTeX/wuff');

% write to TikZ file
if (cfg.write_tikz)
	tfile_h = fopen([export_name, '.tikz'], 'w');
	fprintf(tfile_h, '%% created with surf3tikz written by Johannes Schlichenmaier\n');
	fprintf(tfile_h, '%% based on surf2TikZ written by Fabian Roos\n');
	fprintf(tfile_h, '\\begin{tikzpicture}\n');
	fprintf(tfile_h, '\t \\begin{axis}[\n');
	fprintf(tfile_h, '\t \t grid,\n');
	fprintf(tfile_h, '\t \t enlargelimits = false,\n');
	fprintf(tfile_h, '\t \t zmin = %f,\n', global_data_limits_z(1));
	fprintf(tfile_h, '\t \t zmax = %f,\n', global_data_limits_z(2));
	fprintf(tfile_h, '\t \t xmin = %f,\n', global_data_limits_x(1));
	fprintf(tfile_h, '\t \t xmax = %f,\n', global_data_limits_x(2));
	fprintf(tfile_h, '\t \t ymin = %f,\n', global_data_limits_y(1));
	fprintf(tfile_h, '\t \t ymax = %f,\n', global_data_limits_y(2));
	fprintf(tfile_h, '\t \t xlabel = {%s},\n', xlabel_txt);
	fprintf(tfile_h, '\t \t ylabel = {%s},\n', ylabel_txt);
	fprintf(tfile_h, '\t \t colorbar,\n');
	fprintf(tfile_h, '\t \t colorbar style = {%%,\n');
	fprintf(tfile_h, '\t \t \t ylabel = {%s},\n', colorbar_label);
	fprintf(tfile_h, '\t \t },\n');
	fprintf(tfile_h, '\t \t point meta min = %f,\n', colorbar_limits(1));
	fprintf(tfile_h, '\t \t point meta max = %f,\n', colorbar_limits(2));
	fprintf(tfile_h, '\t ]\n');
	fprintf(tfile_h, '\t \t \\addplot3 graphics[\n');
	fprintf(tfile_h, '\t \t \t points={%% important\n');
	for i=1:size(tikz_support_points,1)
		fprintf(tfile_h, '\t \t \t \t (%f,%f,%f) => (%f,%f)\n', tikz_support_points(i,1), ...
			tikz_support_points(i,2), ...
			tikz_support_points(i,3), ...
			pt_point_positions(i,1), ...
			pt_point_positions(i,2));
	end
	fprintf(tfile_h, '\t \t }]{%s.png};\n', export_fname);
	fprintf(tfile_h, '\t \\end{axis}\n');
	fprintf(tfile_h, '\\end{tikzpicture}');
	fclose(tfile_h);
end

close(plot_handles.figure)

% edge detection
% no! noooooooooooooooooooooooooooooooooooooooooo
end

