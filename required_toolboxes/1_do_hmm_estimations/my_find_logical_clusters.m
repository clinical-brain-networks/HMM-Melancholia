function out = my_find_logical_clusters(v)

info = struct();

i=1;
if v(i) == 0
    started = 0;
else
    started = 1;
    this_b = 1;
    this_e = 1;
    this_count = 1;
end


for i=2:numel(v)
    if v(i)==1 && started
        try
        this_count = this_count + 1;
        this_e = this_e + 1;
        catch
        keyboard;
        end
        
    end
    
    if v(i) == 0 && started
        % keyboard;
        started = 0;
        info(end+1).b = this_b;
        info(end).e = this_e;
        info(end).count = this_count;
    end
    
    if v(i) == 1 && ~started
        started = 1;
        this_b = i;
        this_e = i;
        this_count = 1;
    end
    
end

out=info;
        