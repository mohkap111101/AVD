function output = isa_func(param,height)

if (param == "a")
    [~,output,~,~] = atmosisa(distdim(height,'ft','m'));

elseif (param == "rho")
    [~,~,~,output] = atmosisa(distdim(height,'ft','m'));

elseif (param == "t")
    [output,~,~,~] = atmosisa(distdim(height,'ft','m'));

elseif (param == "p")
    [~,~,output,~] = atmosisa(distdim(height,'ft','m'));

else 
    disp("Error with ISA function")
end

end