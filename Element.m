classdef Element < handle
    
    properties
        id
        posnode
        negnode
        type
        value
        volres
        curres
        dep %it's for dependent source;which elements it depend to.
    end
    
    methods
        function obj = Element(id,posnode,negnode,type,value,de)
            obj.id=id;
            obj.posnode=posnode;
            obj.negnode=negnode;
            obj.type=type;
            obj.value=value;    
            obj.dep=de;
        end
        
    end
end

