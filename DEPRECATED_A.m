classdef A
    properties
        prop
    end
    methods
        function obj = A(k)
            obj.prop = obj.assign_prop(k);
        end
        function prop = assign_prop(~, k)
            prop = k;
        end
    end
end