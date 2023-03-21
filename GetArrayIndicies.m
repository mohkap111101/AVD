function [index, index_1] = GetArrayIndicies(num, a)

    for i = 2:length(a)
        if num <= a(i)
            index = i-1;
            index_1 = i;
            
            return
        end
    end

end