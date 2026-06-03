function defs = muscle_reflex_param_defs()
% Name: muscle_reflex_param_defs
% Description: Single source of truth for the 28-element muscle-reflex
%   parameter vector layout.  Returns a struct array with fields:
%     .idx   — position in the 28×1 parameter vector
%     .path  — cell array of field names for struct assignment
%              e.g. {'tib','KL'} → para.tib.KL
%     .lb    — lower clip bound (for CMA-ES)
%     .ub    — upper clip bound (for CMA-ES)
%
% Usage:
%   defs = muscle_reflex_param_defs();
%   % clip:  fixedPara(d.idx) = min(max(para(d.idx), d.lb), d.ub);
%   % assign: para.(d.path{1}).(d.path{2}) = arx(d.idx);

defs = [ ...
    struct('idx', 1,  'path', {{'tib',              'KL'}},      'lb', -10, 'ub', 10,  'note', 'tibialis anterior: length gain')
    struct('idx', 2,  'path', {{'tib',              'L0'}},      'lb', 0.5, 'ub', 0.8, 'note', 'tibialis anterior: length offset')
    struct('idx', 3,  'path', {{'tib_sol',          'KF'}},      'lb', -10, 'ub', 10,  'note', 'tibialis anterior: soleus force coupling')
    struct('idx', 4,  'path', {{'sol',              'KF'}},      'lb', -10, 'ub', 10,  'note', 'soleus: force gain')
    struct('idx', 5,  'path', {{'gas',              'KF'}},      'lb', -10, 'ub', 10,  'note', 'gastrocnemius: force gain')
    struct('idx', 6,  'path', {{'ili_pelvis_tilt',  'KP'}},      'lb', -10, 'ub', 10,  'note', 'iliopsoas←pelvis tilt: position gain')
    struct('idx', 7,  'path', {{'ili_pelvis_tilt',  'KV'}},      'lb', -10, 'ub', 10,  'note', 'iliopsoas←pelvis tilt: velocity gain')
    struct('idx', 8,  'path', {{'ili_pelvis_tilt',  'C0'}},      'lb', -1,  'ub', 1,   'note', 'iliopsoas←pelvis tilt: constant offset')
    struct('idx', 9,  'path', {{'ili',              'C0'}},      'lb', -1,  'ub', 1,   'note', 'iliopsoas: constant stimulation')
    struct('idx', 10, 'path', {{'ili',              'KL'}},      'lb', -10, 'ub', 10,  'note', 'iliopsoas: length gain')
    struct('idx', 11, 'path', {{'ili',              'L0'}},      'lb', 0,   'ub', 2,   'note', 'iliopsoas: length offset')
    struct('idx', 12, 'path', {{'ili_pelvis_tilt',  'P02'}},     'lb', -1,  'ub', 1,   'note', 'iliopsoas←pelvis tilt: offset #2')
    struct('idx', 13, 'path', {{'ili_pelvis_tilt',  'KP2'}},     'lb', -10, 'ub', 10,  'note', 'iliopsoas←pelvis tilt: position gain #2')
    struct('idx', 14, 'path', {{'ili_pelvis_tilt',  'KV2'}},     'lb', -10, 'ub', 10,  'note', 'iliopsoas←pelvis tilt: velocity gain #2')
    struct('idx', 15, 'path', {{'ili_ham',          'KL'}},      'lb', -10, 'ub', 10,  'note', 'iliopsoas: hamstring length coupling')
    struct('idx', 16, 'path', {{'ili_ham',          'L0'}},      'lb', 0,   'ub', 2,   'note', 'iliopsoas←hamstring: length offset')
    struct('idx', 17, 'path', {{'ham_pelvis_tilt',  'KP'}},      'lb', -10, 'ub', 10,  'note', 'hamstrings←pelvis tilt: position gain')
    struct('idx', 18, 'path', {{'ham_pelvis_tilt',  'KV'}},      'lb', -10, 'ub', 10,  'note', 'hamstrings←pelvis tilt: velocity gain')
    struct('idx', 19, 'path', {{'ham_pelvis_tilt',  'C0'}},      'lb', -1,  'ub', 1,   'note', 'hamstrings←pelvis tilt: constant offset')
    struct('idx', 20, 'path', {{'ham_glu',          'KF'}},      'lb', 0,   'ub', 10,  'note', 'hamstrings: gluteus force coupling')
    struct('idx', 21, 'path', {{'glu_pelvis_tilt',  'KP'}},      'lb', -10, 'ub', 10,  'note', 'gluteus←pelvis tilt: position gain')
    struct('idx', 22, 'path', {{'glu_pelvis_tilt',  'KV'}},      'lb', -10, 'ub', 10,  'note', 'gluteus←pelvis tilt: velocity gain')
    struct('idx', 23, 'path', {{'glu_pelvis_tilt',  'C0'}},      'lb', -1,  'ub', 1,   'note', 'gluteus←pelvis tilt: constant offset')
    struct('idx', 24, 'path', {{'glu',              'KF'}},      'lb', -10, 'ub', 10,  'note', 'gluteus: force gain')
    struct('idx', 25, 'path', {{'vas',              'KF1'}},     'lb', -10, 'ub', 10,  'note', 'vasti: force gain (early stance)')
    struct('idx', 26, 'path', {{'vas',              'KF2'}},     'lb', -10, 'ub', 10,  'note', 'vasti: force gain (late stance)')
    struct('idx', 27, 'path', {{'vas',              'C0'}},      'lb', -1,  'ub', 1,   'note', 'vasti: constant stimulation')
    struct('idx', 28, 'path', {{'vas_knee',         'pos_max'}}, 'lb', -1,  'ub', 0,   'note', 'vasti: knee angle threshold')
];

end
