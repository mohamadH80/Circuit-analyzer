classdef CircuitSimulator < handle
    
    properties
        elements  %array of elements
        n         %nodes number
        delements %array of deleted elements
        couplings %coupling object
    end
    
    methods
        %solve circuit for one frequency
        function [element] = simulatorr(obj,id,omega)
            voltageSourceTransform(obj,omega);
            A=aMatrix(obj);
            Ybo=ybMatrix(obj,omega);
            Vso=vsMatrix(obj);
            Jso=jsMatrix(obj);
            [Yb,Js,Vs]=IdependI(obj,Ybo,Jso,Vso);
            Yn=A*Yb*A';
            Is=A*Yb*Vs-A*Js;
            E=Yn\Is;
            V=A'*E;
            J=Yb*V+Js-Yb*Vs;
            for k=1:length(obj.elements)
                obj.elements(k).volres=[obj.elements(k).volres V(k)];
                obj.elements(k).curres=[obj.elements(k).curres J(k)];
                if obj.elements(k).id==id
                    element=obj.elements(k);
                end
            end
            for m=1:length(obj.delements)
                if obj.delements(m).id==id
                    element=obj.delements(m);
                end
            end                      
        end
        %solve circuit for some frequency
        function [element] = simulator(obj,id,omegas)
            for w=omegas
                if w==0
                    e=simulatorr(obj,id,1e-10);
                    estimate(obj);
                else
                    e=simulatorr(obj,id,w);
                end
            end
            element=e;
        end
        %find A matrix
        function A = aMatrix(obj)
            b=length(obj.elements);
            a = zeros(obj.n-1,b);
            for i=1:obj.n-1
                for j=1:length(obj.elements)
                    if obj.elements(j).posnode==i
                        a(i,j)=1;
                    elseif obj.elements(j).negnode==i
                        a(i,j)=-1;
                    end
                end
            end
            A=a;
        end
        
        %conveert series elements with source voltage to super elements
        function voltageSourceTransform(obj,om)
            indxs=zeros(1,length(obj.elements));
            indxs=boolean(indxs);
            for ee=1:length(obj.elements)
                e=obj.elements(ee);
                if e.type=='V'
                    nr=e.posnode;%the node that must be removed                 
                    Adm=0;
                    for w=1:length(obj.elements)
                        if (obj.elements(w).posnode==e.posnode && obj.elements(w).negnode==e.negnode)
                            indxs(w)=true;
                            obj.elements(w).volres=[obj.elements(w).volres e.value];
                            adm=0;
                            if obj.elements(w).type=='R'
                                adm=(obj.elements(w).value)^-1;
                            elseif obj.elements(w).type=='L'
                                adm=(obj.elements(w).value*1i*om)^-1;
                            elseif obj.elements(w).type=='C'
                                adm=obj.elements(w).value*1i*om;
                            end
                            obj.elements(w).curres=[obj.elements(w).curres e.value*adm];
                        elseif (obj.elements(w).posnode==e.negnode && obj.elements(w).negnode==e.posnode)
                            indxs(w)=true;
                            obj.elements(w).volres=[obj.elements(w).volres -e.value];
                            adm=0;
                            if obj.elements(w).type=='R'
                                adm=(obj.elements(w).value)^-1;
                            elseif obj.elements(w).type=='L'
                                adm=(obj.elements(w).value*1i*om)^-1;
                            elseif obj.elements(w).type=='C'
                                adm=obj.elements(w).value*1i*om;
                            end
                            obj.elements(w).curres=[obj.elements(w).curres -e.value*adm];
                        end
                        Adm=Adm+adm;
                    end
                    
                    obj.elements(ee).curres=[obj.elements(ee).curres e.value*Adm];
                    
                    for ww=1:length(obj.elements)
                        w=obj.elements(ww);
                        if w.id~=e.id && (w.posnode==nr || w.negnode==nr)
                            w.type=[w.type 'T'];
                            w.value=[w.value e.value];
                            if w.posnode==nr
                                w.posnode=e.negnode;
                            else
                                w.negnode=e.negnode;
                            end
                        end
                        obj.elements(ww)=w;
                        
                    end
   
                    indxs(ee)=true;                    
                end
                
            end
            
            indxs=boolean(indxs);
            obj.delements=obj.elements(indxs);
            obj.elements(indxs)=[];
            
        end
        %find Yb matrix
        function [Yb] = ybMatrix(obj,omega)
            b=length(obj.elements);
            y=zeros(b,b);
            for i=1:b
                e=obj.elements(i);
                if e.type(1)=='R'
                    y(i,i)=(e.value(1))^-1;
                elseif e.type(1)=='L'
                    k=findK(obj,e.id);
                    y(i,i)=(e.value(1)*(1-k^2)*omega*1i)^-1;
                elseif e.type(1)=='C'
                    y(i,i)=e.value(1)*omega*1i;
                elseif e.type=='Iv'
                    j=findanElementinElements(obj,e.dep);
                    y(i,j)=e.value;
                end
            end
            Yb=y;
        end
        %find Js & Vs matrix
        function [Js] = jsMatrix(obj)
            b=length(obj.elements);
            m=zeros(b,1);
            for k=1:b
                t=obj.elements(k).type;
                if t(1)=='I' && t~="Ii" && t~="Iv"
                    m(k)=obj.elements(k).value(1);
                end
            end
            Js=m;
        end
        function [Vs] = vsMatrix(obj)
            b=length(obj.elements);
            m=zeros(b,1);
            for k=1:b
                if length(obj.elements(k).type)==2
                    if obj.elements(k).type(2)=='T'
                        m(k)=obj.elements(k).value(2);
                    end
                end
            end
            Vs=m;
        end
        
        function estimate(obj)
            for k=1:length(obj.elements)
                for kv=1:length(obj.elements(k).volres)
                    if real(obj.elements(k).volres(kv))<=1e-7
                        obj.elements(k).volres(kv)=complex(0,imag(obj.elements(k).volres(kv)));
                    end
                    if imag(obj.elements(k).volres(kv))<=1e-7
                        obj.elements(k).volres(kv)=real(obj.elements(k).volres(kv));
                    end
                end
                for kj=1:length(obj.elements(k).curres)
                    if real(obj.elements(k).curres(kj))<=1e-7
                        obj.elements(k).curres(kj)=complex(0,imag(obj.elements(k).curres(kj)));
                    end
                    if imag(obj.elements(k).curres(kj))<=1e-7
                        obj.elements(k).curres(kj)=real(obj.elements(k).curres(kj));
                    end
                end
            end
                        
            for k=1:length(obj.delements)
                for kv=1:length(obj.delements(k).volres)
                    if real(obj.delements(k).volres(kv))<=1e-7
                        obj.delements(k).volres(kv)=complex(0,imag(obj.delements(k).volres(kv)));
                    end
                    if imag(obj.delements(k).volres(kv))<=1e-7
                        obj.delements(k).volres(kv)=real(obj.delements(k).volres(kv));
                    end
                end
                for kj=1:length(obj.delements(k).curres)
                    if real(obj.delements(k).curres(kj))<=1e-7
                        obj.delements(k).curres(kj)=complex(0,imag(obj.delements(k).curres(kj)));
                    end
                    if imag(obj.delements(k).curres(kj))<=1e-7
                        obj.delements(k).curres(kj)=real(obj.delements(k).curres(kj));
                    end
                end
            end
        end

        function [Yb,Js,Vs]=IdependI(obj,Y,js,vs)
            for i=1:length(obj.elements)
                if obj.elements(i).type=='Ii'
                    j=findanElementinElements(obj,obj.elements(i).dep);
                    r=obj.elements(i).value;
                    Y(i,:)=Y(j,:);
                    js(i)=r*js(j);
                    vs(i)=r*vs(j);
                end
            end
            Yb=Y;
            Js=js;
            Vs=vs;
        end
        
        function [j]=findanElementinElements(obj,id)
            for i=1:length(obj.elements)
                if obj.elements(i).id==id
                    j=i;
                end
            end
        end
        
        function k=findK(obj,id)
            k=0;
            for i=obj.couplings
                if i.id1==id || i.id2==id
                    k=i.k;
                end
            end
        end
        
    end
end
