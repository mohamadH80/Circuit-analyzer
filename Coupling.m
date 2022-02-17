classdef Coupling < handle
    %coupling between inductors with ids 1&2 by coupling Coefficient k
    properties
        id1 
        id2
        k
    end
    
    methods
        function obj = Coupling(id1,id2,k)
            obj.k=k;
            obj.id1=id1;
            obj.id2=id2;
        end
    end
end

