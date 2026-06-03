classdef (Abstract) Base_optimization < handle
    %BASE_OPTIMIZATION 此处显示有关此类的摘要
    %   此处显示详细说明

    properties
        optParaNum % number of parameters to be optimized
        optParaRecord % record of history of parameters to be optimized: best fitness and the corresponding parameters, and the iteration when it occured
        core % intrinsic parameters of the optimization method
    end

    methods (Abstract)
        update_core(obj) % update optPara and hyperPara
    end

    methods
        function fixedPara = fix_parameters(obj, para)
            % Name: fix_parameters
            % Description: clip parameters so that their values are within
            %   reasonable ranges.  Bounds are defined in the single-source
            %   muscle_reflex_param_defs().
            defs = muscle_reflex_param_defs();
            fixedPara = para;
            for i = 1 : numel(defs)
                d = defs(i);
                fixedPara(d.idx) = obj.min_max(para(d.idx), d.lb, d.ub);
            end
        end

        function x_ = min_max(~, x, lowerTh, upperTh)
            % Name: min_max
            % Description: clip parameters within [lowerTh, upperTh].
            x_ = x;
            x_(x_<lowerTh) = lowerTh;
            x_(x_>upperTh) = upperTh;
        end
    end
end

