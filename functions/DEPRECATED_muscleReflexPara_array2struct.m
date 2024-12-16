% function muscleReflexParaStruct = muscleReflexPara_array2struct(muscleReflexParaArray)
% % Name: muscleReflexPara_array2struct
% % Description: convert the muscle-reflex parameters from array form to
% % struct form.
% 
% muscleReflexParaStruct.tib.KL                     = muscleReflexParaArray(1);
% muscleReflexParaStruct.tib.L0                     = muscleReflexParaArray(2);
% muscleReflexParaStruct.tib_sol.KF                 = muscleReflexParaArray(3);
% muscleReflexParaStruct.sol.KF                     = muscleReflexParaArray(4);
% muscleReflexParaStruct.gas.KF                     = muscleReflexParaArray(5);
% muscleReflexParaStruct.ili_pelvis_tilt.KP         = muscleReflexParaArray(6);
% muscleReflexParaStruct.ili_pelvis_tilt.KV         = muscleReflexParaArray(7);
% muscleReflexParaStruct.ili_pelvis_tilt.C0         = muscleReflexParaArray(8);
% muscleReflexParaStruct.ili.C0                     = muscleReflexParaArray(9);
% muscleReflexParaStruct.ili.KL                     = muscleReflexParaArray(10);
% muscleReflexParaStruct.ili.L0                     = muscleReflexParaArray(11);
% muscleReflexParaStruct.ili_pelvis_tilt.P02        = muscleReflexParaArray(12);
% muscleReflexParaStruct.ili_pelvis_tilt.KP2        = muscleReflexParaArray(13);
% muscleReflexParaStruct.ili_pelvis_tilt.KV2        = muscleReflexParaArray(14);
% muscleReflexParaStruct.ili_ham.KL                 = muscleReflexParaArray(15);
% muscleReflexParaStruct.ili_ham.L0                 = muscleReflexParaArray(16);
% muscleReflexParaStruct.ham_pelvis_tilt.KP         = muscleReflexParaArray(17);
% muscleReflexParaStruct.ham_pelvis_tilt.KV         = muscleReflexParaArray(18);
% muscleReflexParaStruct.ham_pelvis_tilt.C0         = muscleReflexParaArray(19);
% % muscleReflexParaStruct.ham.KF                     = muscleReflexParaArray(20);
% muscleReflexParaStruct.ham_glu.KF                 = muscleReflexParaArray(20);
% muscleReflexParaStruct.glu_pelvis_tilt.KP         = muscleReflexParaArray(21);
% muscleReflexParaStruct.glu_pelvis_tilt.KV         = muscleReflexParaArray(22);
% muscleReflexParaStruct.glu_pelvis_tilt.C0         = muscleReflexParaArray(23);
% muscleReflexParaStruct.glu.KF                     = muscleReflexParaArray(24);
% % muscleReflexParaStruct.glu.C0                     = muscleReflexParaArray(25);
% muscleReflexParaStruct.vas.KF1                    = muscleReflexParaArray(25); % 26 -> 25
% muscleReflexParaStruct.vas.KF2                    = muscleReflexParaArray(26);
% muscleReflexParaStruct.vas.C0                     = muscleReflexParaArray(27);
% muscleReflexParaStruct.vas_knee.pos_max           = muscleReflexParaArray(28);
% 
% end