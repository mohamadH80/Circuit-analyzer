classdef Input < handle
    
    properties
        circuit
    end
    
    methods
        function obj = Input()
            obj.circuit = CircuitSimulator();
            obj.getData();
             
        end
        
        function getData(obj)
            obj.circuit.n = input('input nodes number:');
            b = input('input elements number:');
            obj.circuit.elements=[];
            for i=1:b
                de=[];
                id=input('id:');
                posnode=input('posnode:');
                negnode=input('negnode:');
                type=input('type:');
                if (type=='Iv' | type=='Ii' | type=='Vi' | type=='Vv')
                    de=input('inter id of dependents element:');
                end
                value=input('value:');
                e=Element(id,posnode,negnode,type,value,de);
                obj.circuit.elements=[obj.circuit.elements e];
            end   
            getCouplingData(obj);
            
        end
        
        function getCouplingData(obj)
            cn=input('numbers of coupling in circuit:');
            for i=1:cn
                id1=input('id1:');
                id2=input('id2:');
                k=input('coupling Coefficient:');
                c=Coupling(id1,id2,k);
                obj.circuit.couplings=[obj.circuit.couplings c];
            end
        end
        
    end
      
end

