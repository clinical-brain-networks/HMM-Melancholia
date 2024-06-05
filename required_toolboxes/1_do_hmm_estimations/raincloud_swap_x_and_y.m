function raincloud_swap_x_and_y(in)

for i_object = 1:numel(in)
    
    this_obj = in{i_object};
        switch get(this_obj,'type')
            
            
            case {'area', 'scatter', 'line'}
                
                old_xdata = get(this_obj,'XData');
                old_ydata = get(this_obj,'YData');
                
                new_xdata = old_ydata;
                new_ydata = old_xdata;
                
                set(this_obj, 'XData', new_xdata);
                set(this_obj, 'YData', new_ydata);
             
                
            case 'rectangle'
                
                old_pos = get(this_obj,'position');
                new_pos = [old_pos(2) old_pos(1) old_pos(4) old_pos(3)];
                set(this_obj,'position', new_pos);
                
                
        end
end

parent_axes_obj = get(in{1},'parent');
old_parent_xlim = get(parent_axes_obj,'Xlim');
old_parent_ylim = get(parent_axes_obj,'Ylim');
set(parent_axes_obj,'Xlim', old_parent_ylim);
set(parent_axes_obj,'Ylim', old_parent_xlim);
                
                
                