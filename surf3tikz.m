function [pt_point_positions, act_data_points, colorbar_limits] = surf3tikz(h_figure, export_name, cfg, debug)
%SURF3TIKZ converts a surf plot into a tikz file
% This function takes a figure containing a surf plot and creates a tikz file and an according png
% file. To achieve this, it automatically places data cursors on extreme points of the surf plot and
% subsequently detects their pixel position
%
% INPUT:
%   h_figure: Figure handle containing the figure plot. Ideally also contains a colorbar.
%   export_name: Base filename of the exported tikz and png file
%   cfg: optional config struct
%      

ppi = get(0,'ScreenPixelsPerInch');
inch2point_ratio = 1/72;
export_dpi = 300;

% if nargin < 2 
% 	error('');

plot_handles.figure = copyobj(h_figure,0);
plot_handles.axes = plot_handles.figure.CurrentAxes;

colorbar_limits = [0,100];
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

% determine surf object
for i=1:numel(plot_handles.axes.Children)
	if isa(plot_handles.axes.Children(i), 'matlab.graphics.chart.primitive.Surface')
		plot_handles.surf = plot_handles.axes.Children(i);
	end
end

axis off
title('')

% keyboard

x_limits(1) = min(plot_handles.surf.XData(:));
x_limits(2) = max(plot_handles.surf.XData(:));
y_limits(1) = min(plot_handles.surf.YData(:));
y_limits(2) = max(plot_handles.surf.YData(:));
[x_max_idx,y_max_idx] = size(plot_handles.surf.XData);
[~, maxidx] = max(plot_handles.surf.ZData(:));
[~, minidx] = min(plot_handles.surf.ZData(:));


data_points = zeros(10,3);
act_data_points = zeros(10,3);
mid_row_pixel = zeros(10,1);
mid_col_pixel = zeros(10,1);

data_point_idc(1) = sub2ind([x_max_idx,y_max_idx], 1, 1);
data_point_idc(2) = sub2ind([x_max_idx,y_max_idx], 1, round(y_max_idx/2));
data_point_idc(3) = sub2ind([x_max_idx,y_max_idx], 1, y_max_idx);
data_point_idc(4) = sub2ind([x_max_idx,y_max_idx], round(x_max_idx/2), 1);
data_point_idc(5) = sub2ind([x_max_idx,y_max_idx], x_max_idx, 1);
data_point_idc(6) = sub2ind([x_max_idx,y_max_idx], round(x_max_idx/2), round(y_max_idx/2));
data_point_idc(7) = sub2ind([x_max_idx,y_max_idx], x_max_idx, round(y_max_idx/2));
data_point_idc(8) = sub2ind([x_max_idx,y_max_idx], x_max_idx, y_max_idx);
data_point_idc(9) = maxidx;
data_point_idc(10) = minidx;

for i = 1:numel(data_point_idc)
	data_points(i,:) = [plot_handles.surf.XData(data_point_idc(i)), ...
	plot_handles.surf.YData(data_point_idc(i)), ...
	plot_handles.surf.ZData(data_point_idc(i))];
end

% tmp_point = cursor_info.Position

plot_handles.cursor_manager = datacursormode(plot_handles.figure);
plot_handles.cursor_manager.removeAllDataCursors
plot_handles.cursor_manager.UpdateFcn = 'datatip_dummy_function';
plot_handles.cursor = createDatatip(plot_handles.cursor_manager, ...
	plot_handles.surf);

plot_handles.cursor.MarkerEdgeColor = [1,1,1];
plot_handles.cursor.MarkerFaceColor = [0,0,0];

for i = 1:numel(data_point_idc)
	% setting cursor position
	plot_handles.cursor.Position = data_points(i,:);
	act_data_points(i,:) = plot_handles.cursor.Position;
	
	% create 2D image
	frame_data = getframe(plot_handles.figure);
	[image_data,~] = frame2im(frame_data);
	flat_image_data = image_data(:,:,1)+image_data(:,:,2)+image_data(:,:,3);
	
	% find black rectangle
% 	min_idc = find(flat_image_data(:) == min(flat_image_data(:)));
	min_idc = find(flat_image_data(:) == 0);
	[min_row, min_col] = ind2sub(size(flat_image_data), min_idc);
	
	% extract center point of black rectangle
	
	unique_min_row = unique(min_row);
	unique_min_col = unique(min_col);
	
	mid_row_pixel(i) = round((max(unique_min_row) + min(unique_min_row))/2);
	mid_col_pixel(i) = round((max(unique_min_col) + min(unique_min_col))/2);
	
end

%tidy up
plot_handles.cursor_manager.removeAllDataCursors

if debug
	plot_handles.cursor_manager.UpdateFcn = [];
	
	for i = 1:size(data_points,1)
		debug_cursor = createDatatip(plot_handles.cursor_manager, plot_handles.surf);
		debug_cursor.Position = data_points(i,:);
	end
	
	debug_frame_data = getframe(plot_handles.figure);
	[debug_image_data,~] = frame2im(debug_frame_data);
	
	plot_handles.cursor_manager.removeAllDataCursors
	
	for i = 1:size(data_points,1)
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
pt_point_positions = pgf_pixel_points./ppi./inch2point_ratio;

% print png and make transparent
print(plot_handles.figure, export_name, '-dpng', ['-r' num2str(export_dpi)]);
system(['mogrify -transparent white ', export_name, '.png']);

tfile_h = fopen([export_name, '.tikz'], 'w');

fprintf(tfile_h, '%% created with surf3tikz written by Johannes Schlichenmaier\n');
fprintf(tfile_h, '%% based on surf2TikZ written by Fabian Roos\n');
fprintf(tfile_h, '\\begin{tikzpicture}\n');
fprintf(tfile_h, '\t \\begin{axis}[\n');
fprintf(tfile_h, '\t \t grid,\n');
fprintf(tfile_h, '\t \t enlargelimits = false,\n');
fprintf(tfile_h, '\t \t zmin = %f,\n', colorbar_limits(1));
fprintf(tfile_h, '\t \t zmax = %f,\n', colorbar_limits(2));
fprintf(tfile_h, '\t \t xmin = %f,\n', x_limits(1));
fprintf(tfile_h, '\t \t xmax = %f,\n', x_limits(2));
fprintf(tfile_h, '\t \t ymin = %f,\n', y_limits(1));
fprintf(tfile_h, '\t \t ymax = %f,\n', y_limits(2));
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
for i=1:size(data_points,1)
	fprintf(tfile_h, '\t \t \t \t (%f,%f,%f) => (%f,%f)\n', act_data_points(i,1), ...
		act_data_points(i,2), ...
		act_data_points(i,3), ...
		pt_point_positions(i,1), ...
		pt_point_positions(i,2));
end
fprintf(tfile_h, '\t \t }]{%s.png};\n', export_name);
fprintf(tfile_h, '\t \\end{axis}\n');
fprintf(tfile_h, '\\end{tikzpicture}');

fclose(tfile_h);

close(plot_handles.figure)

% edge detection
% no! noooooooooooooooooooooooooooooooooooooooooo
end


function [ output_text ] = datatip_dummy_function( ~, event_obj )
%DATATIP_DUMMY_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
output_text = '';
end

