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

