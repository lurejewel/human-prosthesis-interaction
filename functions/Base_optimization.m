classdef (Abstract) Base_optimization
    %BASE_OPTIMIZATION 此处显示有关此类的摘要
    %   此处显示详细说明

    properties
        optPara % parameters to be optimized
        optParaNum % number of parameters to be optimized
        optParaRecord % record of history of parameters to be optimized: best fitness and the corresponding parameters, and the iteration when it occured
        hyperPara % parameters for optimization
        intrinsicPara % intrinsic parameters of the optimization method
    end

    methods (Abstract)
        update(obj) % update optPara and hyperPara
    end

    methods
        function fixedPara = fix_parameters(~, para)
            % Name: fix_parameters
            % Description: clip parameters so that their values are within
            % reasonable ranges.
            fixedPara(1) = min_max(para(1), -10, 10); % para.tib.KL
            fixedPara(2) = min_max(para(2), 0.5, 0.8); % para.tib.L0
            fixedPara(3) = min_max(para(3), -10, 10); % para.tib_sol.KF
            fixedPara(4) = min_max(para(4), -10, 10); % para.sol.KF
            fixedPara(5) = min_max(para(5), -10, 10); % para.gas.KF
            fixedPara(6) = min_max(para(6), -10, 10); % para.ili_pelvis_tilt.KP
            fixedPara(7) = min_max(para(7), -10, 10); % para.ili_pelvis_tilt.KV
            % fixedPara(8) = min_max(para(8), -0.5, 0.5); % para.pelvis_tilt.P0
            fixedPara(8) = min_max(para(8), -1, 1); % para.ili_pelvis_tilt.C0
            fixedPara(9) = min_max(para(9), -1, 1); % para.ili.C0
            fixedPara(10) = min_max(para(10), -10, 10); % % para.ili.KL
            fixedPara(11) = min_max(para(11), 0, 2); % para.ili.L0
            fixedPara(12) = min_max(para(12), -1, 1); % para.ili_pelvis_tilt.P02
            fixedPara(13) = min_max(para(13), -10, 10); % para.ili_pelvis_tilt.KP2
            fixedPara(14) = min_max(para(14), -10, 10); % para.ili_pelvis_tilt.KV2
            fixedPara(15) = min_max(para(15), -10, 10); % para.ili_ham.KL
            fixedPara(16) = min_max(para(16), 0, 2); % para.ili_ham.L0
            fixedPara(17) = min_max(para(17), -10, 10); % para.ham_pelvis_tilt.KP
            fixedPara(18) = min_max(para(18), -10, 10); % para.ham_pelvis_tilt.KV
            fixedPara(19) = min_max(para(19), -1, 1); % para.ham_pelvis_tilt.C0
            fixedPara(20) = min_max(para(20), 0, 10); % para.ham_glu.KF
            fixedPara(21) = min_max(para(21), -10, 10); % para.glu_pelvis_tilt.KP
            fixedPara(22) = min_max(para(22), -10, 10); % para.glu_pelvis_tilt.KV
            fixedPara(23) = min_max(para(23), -1, 1); % para.glu_pelvis_tilt.C0
            fixedPara(24) = min_max(para(24), -10, 10); % para.glu.KF
            % fixedPara(25) = min_max(para(25), -1, 1); % para.glu.C0
            fixedPara(25) = min_max(para(25), -10, 10); % para.vas.KF1
            fixedPara(26) = min_max(para(26), -10, 10); % para.vas.KF2
            fixedPara(27) = min_max(para(27), -1, 1); % para.vas.C0
            fixedPara(28) = min_max(para(28), -1, 0); % para.vas_knee.pos_max
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

