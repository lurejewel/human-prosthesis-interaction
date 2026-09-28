function def = foot_contact_defs()
% Name: foot_contact_defs
% Description: single source of truth for the foot-ground contact model of
%   human0918.osim.
%
%   In this model each foot's ground contact is NOT a single lumped force
%   anymore: it is split into three independent HuntCrossleyForce components
%   (heel / toe / tip), each acting on its own ContactSphere:
%
%       force  'ground_heel_r'  <->  ContactSphere 'heel_r'
%       force  'ground_toe_r'   <->  ContactSphere 'toe_r'
%       force  'ground_tip_r'   <->  ContactSphere 'tip_r'
%       ... and the same three for the left foot ('_l')
%
%   The resultant ground reaction force is the sum of the three components
%   and the centre of pressure is their normal-force-weighted application
%   point (see cal_grf.m).  Keeping the names here means a future change to
%   the contact model only has to be made in one place.
%
% Output struct:
%   def.sides       {'r','l'}                     - foot side order
%   def.points      {'heel','toe','tip'}          - contact point order
%   def.nPoints     number of contact points per foot
%   def.geomNames   {nSides x nPoints}            - ContactSphere names
%   def.forceNames  {nSides x nPoints}            - HuntCrossleyForce names

def.sides   = {'r', 'l'};
def.points  = {'heel', 'toe', 'tip'};
def.nPoints = numel(def.points);

def.geomNames  = cell(numel(def.sides), def.nPoints);
def.forceNames = cell(numel(def.sides), def.nPoints);
for i = 1 : numel(def.sides)
    for j = 1 : def.nPoints
        def.geomNames{i, j}  = sprintf('%s_%s', def.points{j}, def.sides{i});
        def.forceNames{i, j} = ['ground_' def.geomNames{i, j}];
    end
end

end
