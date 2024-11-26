% -------------------------------------------------------------------------
% model_construction.m
%
% AUTHOR: Wei Jin, Peking Univ., wjin24@stu.pku.edu.cn
% LAST UPDATE: 24-10-18
%
% FUNCTION: define and assemble components of the coupled human-prosthesis
% musculoskeletal model.
% 
% OUTPUT: a musculoskeletal model (.osim file), which can be visualized
% either through MATLAB or in OpenSim GUI
% 
% DETIALED COMPONENTS:
% [gravity]
% [ground]
% [bodies]
% - pelvis: pelvis
% - thigh, both sides: femur_r/l
% - shank, both sides: tibia_r/l
% - foot, both sides: calcn_r/l
% - HAT: torso
% [joints]
% - ground ~ pelvis: ground_pelvis, with 
% - pelvis ~ HAT: back
% - pelvis ~ thigh, both sides: hip_r/l
% - thigh ~ shank, both sides: knee_r/l
% - shank ~ foot, both sides: ankle_r/l
% [muscles]
% - hamstrings: hamstrings_r/l
% - glutues maximum: glut_max_r/l
% - ilipsoas: ilipsoas_r/l
% - vastus: vasti_r/l
% - gastrocnemius: gastroc_r/l
% - soleus: soleus_r/l
% - tibialis anterior: tibia_r/l
% [contact geometries]
% 
% [forces](other than muscles)
% 
% [markers]
% - sternum: sternum
% - acromium: acromium_r/l
% - top head: top_head
% - anterior superior iliac spine: asis_r/l
% - sacrum: sacrum_v
% - upper/front/rear thigh: thigh_upper/front/rear_r/l
% - lateral/medial knee: knee_lat/med_r/l
% - upper/front/rear shank: shank_upper/front/rear_r/l
% - lateral/medial ankle: ankle_lat/med_r/l
% - heel: heel_r/l
% - superior/lateral midfoot: midfoot_sup/lat_r/l
% - lateral/medial/tip toe: toe_lat/med/tip_r/l
%
%
% HISTORY:
% 24-10-11: basic file structure construction; [bodies] finished
% 24-10-12: [joints], [muscles] finished
% 24-10-14: [markers], [contact] finished
% 24-10-18: [markers] removed (seperately defined in
% assets/scale_markerSet.xml; [joints] hip DOF refined
% 24-10-25: [joints] ground_pelvis DOF refined
% TODO: rra, so, cmc trial; muscle reflex; cma-es optimization
% 
% -------------------------------------------------------------------------
addpath('assets\', 'model\','model\Geometry\');
import org.opensim.modeling.* % import opensim libraries

%% Initialization
model = Model(); % init a model object
model.setName('coupled_human-prosthesis_model');
model.setUseVisualizer(true); % visualizer open during simulation

% --- gravity --- %
model.setGravity(Vec3(0, -9.80665, 0));

% --- ground --- %
ground = model.getGround();

%% Bodies
% 目前先实现gait0914（SCONE模型）的复杂度
% pelvis
bodies.pelvis = Body('pelvis', 11.777, Vec3(-0.0707,0,0), Inertia(0.1028,0.0871,0.0579)); % construction function, para: <name>, <mass>, <COM>, <inertia>
bodies.pelvis.attachGeometry(Mesh('sacrum.vtp')); % add visual mesh files to the body
bodies.pelvis.attachGeometry(Mesh('pelvis.vtp'));
bodies.pelvis.attachGeometry(Mesh('l_pelvis.vtp'));
model.addBody(bodies.pelvis);

% right thigh: femur_r
bodies.femur_r = Body('femur_r', 9.3014, Vec3(0,-0.17,0), Inertia(0.1339,0.0351,0.1412));
bodies.femur_r.attachGeometry(Mesh('femur.vtp'));
model.addBody(bodies.femur_r);

% right shank: tibia_r
bodies.tibia_r = Body('tibia_r', 3.7075, Vec3(0,-0.1867,0), Inertia(0.0504,0.0051,0.0511));
bodies.tibia_r.attachGeometry(Mesh('tibia_r.vtp'));
bodies.tibia_r.attachGeometry(Mesh('fibula.vtp'));
model.addBody(bodies.tibia_r);

% right foot: calcn_r
bodies.calcn_r = Body('calcn_r', 1.25, Vec3(0.1,0.03,0), Inertia(0.0014,0.0039,0.0041));
frames.calcn_r1 = PhysicalOffsetFrame(); % manually create an offset frame, only for attaching visual mesh files (talus.vtp in this case)
frames.calcn_r1.setName('calcn_r_geom_frame_2');
frames.calcn_r1.set_translation(Vec3(0.04877,0.04195,-0.00792)); % translational offset w.r.t parent frame
frames.calcn_r1.connectSocket_parent(bodies.calcn_r); % specify parent frame (reference frame)
frames.calcn_r1.attachGeometry(Mesh('talus.vtp'));
bodies.calcn_r.attachGeometry(Mesh('foot.vtp')); % "body" itself can be a frame too
frames.calcn_r2 = PhysicalOffsetFrame();
frames.calcn_r2.setName('calcn_r_geom_frame_3');
frames.calcn_r2.set_translation(Vec3(0.1788,-0.002,0.00108));
frames.calcn_r2.connectSocket_parent(bodies.calcn_r)
frames.calcn_r2.attachGeometry(Mesh('bofoot.vtp'));
bodies.calcn_r.addComponent(frames.calcn_r1);
bodies.calcn_r.addComponent(frames.calcn_r2);
model.addBody(bodies.calcn_r);

% left thigh: femur_l
bodies.femur_l = Body('femur_l', 9.3014, Vec3(0,-0.17,0), Inertia(0.1339,0.0351,0.1412));
bodies.femur_l.attachGeometry(Mesh('l_femur.vtp'));
model.addBody(bodies.femur_l);

% left shank: tibia_l
bodies.tibia_l = Body('tibia_l', 3.7075, Vec3(0,-0.1867,0), Inertia(0.0504,0.0051,0.0511));
bodies.tibia_l.attachGeometry(Mesh('tibia_l.vtp'));
bodies.tibia_l.attachGeometry(Mesh('l_fibula.vtp'));
model.addBody(bodies.tibia_l);

% left foot: calcn_l
bodies.calcn_l = Body('calcn_l', 1.25, Vec3(0.1,0.03,0), Inertia(0.0014,0.0039,0.0041));
frames.calcn_l1 = PhysicalOffsetFrame(); % manually create an offset frame, only for attaching visual mesh files (talus.vtp in this case)
frames.calcn_l1.setName('calcn_l_geom_frame_2');
frames.calcn_l1.set_translation(Vec3(0.04877,0.04195,0.00792));
frames.calcn_l1.connectSocket_parent(bodies.calcn_l);
frames.calcn_l1.attachGeometry(Mesh('l_talus.vtp'));
bodies.calcn_l.attachGeometry(Mesh('l_foot.vtp'));
frames.calcn_l2 = PhysicalOffsetFrame();
frames.calcn_l2.setName('calcn_l_geom_frame_3');
frames.calcn_l2.set_translation(Vec3(0.1788,-0.002,-0.00108));
frames.calcn_l2.connectSocket_parent(bodies.calcn_l);
frames.calcn_l2.attachGeometry(Mesh('l_bofoot.vtp'));
bodies.calcn_l.addComponent(frames.calcn_l1);
bodies.calcn_l.addComponent(frames.calcn_l2);
model.addBody(bodies.calcn_l);

% HAT: torso
bodies.torso = Body('torso', 34.2366, Vec3(-0.03,0.32,0), Inertia(1.4745,0.7555,1.4314));
bodies.torso.attachGeometry(Mesh('hat_spine.vtp'));
bodies.torso.attachGeometry(Mesh('hat_jaw.vtp'));
bodies.torso.attachGeometry(Mesh('hat_skull.vtp'));
bodies.torso.attachGeometry(Mesh('hat_ribs.vtp'));
model.addBody(bodies.torso);

%% Joints
% ground ~ pelvis: ground_pelvis
spatialTransform = SpatialTransform(); % define DOF of the joint here
spatialTransform.updTransformAxis(0).setCoordinateNames(ArrayStr('pelvis_list',1));
spatialTransform.updTransformAxis(0).setFunction(LinearFunction(1,0)); % rx: rotate as it should be (0~5: rx, ry, rz, tx, ty, tz; x-forward, y-upward; z-rightward)
spatialTransform.updTransformAxis(1).setCoordinateNames(ArrayStr('pelvis_rotation',1));
spatialTransform.updTransformAxis(1).setFunction(LinearFunction(1,0));
spatialTransform.updTransformAxis(2).setCoordinateNames(ArrayStr('pelvis_tilt',1)); % funny thing is that the coordinate will NOT be created unless you name it
spatialTransform.updTransformAxis(2).setFunction(LinearFunction(1,0)); % rz: rotate as it should be
spatialTransform.updTransformAxis(3).setCoordinateNames(ArrayStr('pelvis_tx',1));
spatialTransform.updTransformAxis(3).setFunction(LinearFunction(1,0));
spatialTransform.updTransformAxis(4).setCoordinateNames(ArrayStr('pelvis_ty',1));
spatialTransform.updTransformAxis(4).setFunction(LinearFunction(1,0));
spatialTransform.updTransformAxis(5).setCoordinateNames(ArrayStr('pelvis_tz',1));
spatialTransform.updTransformAxis(5).setFunction(LinearFunction(1,0));
joints.ground_pelvis = CustomJoint('ground_pelvis', ground, Vec3(0), Vec3(0), bodies.pelvis, Vec3(0), Vec3(0), spatialTransform); % construction function of CustomJoint, para: <name>, <parent body/frame>, <location in parent>, <orientation in parent>, <child body/frame>, <location in child>, <orientation in parent>, <spatial transform>

model.addJoint(joints.ground_pelvis);
coords.pelvis_list = model.updCoordinateSet().get('pelvis_list'); % stretch coordinate set from the model. `upd` means writable; any changes in coords.pelvis_tilt will be synchronized in the model
coords.pelvis_list.setRange([deg2rad(-90), deg2rad(90)]);
coords.pelvis_rotation = model.updCoordinateSet().get('pelvis_rotation');
coords.pelvis_rotation.setRange([deg2rad(-90), deg2rad(90)]);
coords.pelvis_tilt = model.updCoordinateSet().get('pelvis_tilt');
coords.pelvis_tilt.setRange([deg2rad(-90), deg2rad(90)]);
coords.pelvis_tx = model.updCoordinateSet().get('pelvis_tx');
coords.pelvis_tx.setRange([-5 5]);
coords.pelvis_ty = model.updCoordinateSet().get('pelvis_ty');
coords.pelvis_ty.setDefaultValue(0.95);
coords.pelvis_ty.setRange([-1 2]);
coords.pelvis_tz = model.updCoordinateSet().get('pelvis_tz');
coords.pelvis_tz.setRange([-3 3]);

% pelvis ~ right thigh: hip_r
spatialTransform = SpatialTransform();
spatialTransform.updTransformAxis(0).setCoordinateNames(ArrayStr('hip_adduction_r',1));
spatialTransform.updTransformAxis(0).setFunction(LinearFunction(1,0));
spatialTransform.updTransformAxis(1).setCoordinateNames(ArrayStr('hip_rotation_r',1));
spatialTransform.updTransformAxis(1).setFunction(LinearFunction(1,0));
spatialTransform.updTransformAxis(2).setCoordinateNames(ArrayStr('hip_flexion_r',1));
spatialTransform.updTransformAxis(2).setFunction(LinearFunction(1,0));
spatialTransform.updTransformAxis(3).setFunction(Constant(0));
spatialTransform.updTransformAxis(4).setFunction(Constant(0));
spatialTransform.updTransformAxis(5).setFunction(Constant(0));
joints.hip_r = CustomJoint('hip_r', bodies.pelvis, Vec3(-0.0707,-0.0661,0.0835), Vec3(0), bodies.femur_r, Vec3(0), Vec3(0), spatialTransform);

model.addJoint(joints.hip_r);
coords.hip_adduction_r = model.updCoordinateSet().get('hip_adduction_r');
coords.hip_adduction_r.setRange([deg2rad(-120), deg2rad(120)]);
coords.hip_rotation_r = model.updCoordinateSet().get('hip_rotation_r');
coords.hip_rotation_r.setRange([deg2rad(-120), deg2rad(120)]);
coords.hip_flexion_r = model.updCoordinateSet().get('hip_flexion_r');
coords.hip_flexion_r.setRange([deg2rad(-120), deg2rad(120)]);

% pelvis ~ right thigh: hip_l
spatialTransform = SpatialTransform();
spatialTransform.updTransformAxis(0).setCoordinateNames(ArrayStr('hip_adduction_l',1));
spatialTransform.updTransformAxis(0).setFunction(LinearFunction(1,0));
spatialTransform.updTransformAxis(1).setCoordinateNames(ArrayStr('hip_rotation_l',1));
spatialTransform.updTransformAxis(1).setFunction(LinearFunction(1,0));
spatialTransform.updTransformAxis(2).setCoordinateNames(ArrayStr('hip_flexion_l',1));
spatialTransform.updTransformAxis(2).setFunction(LinearFunction(1,0));
spatialTransform.updTransformAxis(3).setFunction(Constant(0));
spatialTransform.updTransformAxis(4).setFunction(Constant(0));
spatialTransform.updTransformAxis(5).setFunction(Constant(0));
joints.hip_l = CustomJoint('hip_l', bodies.pelvis, Vec3(-0.0707,-0.0661,-0.0835), Vec3(0), bodies.femur_l, Vec3(0), Vec3(0), spatialTransform);

model.addJoint(joints.hip_l);
coords.hip_adduction_l = model.updCoordinateSet().get('hip_adduction_l');
coords.hip_adduction_l.setRange([deg2rad(-120), deg2rad(120)]);
coords.hip_rotation_l = model.updCoordinateSet().get('hip_rotation_l');
coords.hip_rotation_l.setRange([deg2rad(-120), deg2rad(120)]);
coords.hip_flexion_l = model.updCoordinateSet().get('hip_flexion_l');
coords.hip_flexion_l.setRange([deg2rad(-120), deg2rad(120)]);

% right thigh ~ right shank: knee_r
joints.knee_r = PinJoint('knee_r', bodies.femur_r, Vec3(0,-0.396,0), Vec3(0), bodies.tibia_r, Vec3(0), Vec3(0));
joints.knee_r.upd_coordinates(0).setRange([deg2rad(-120), deg2rad(10)])
joints.knee_r.upd_coordinates(0).setName('knee_flexion_r');
model.addJoint(joints.knee_r);

% left thigh ~ left shank: knee_l
joints.knee_l = PinJoint('knee_l', bodies.femur_l, Vec3(0,-0.396,0), Vec3(0), bodies.tibia_l, Vec3(0), Vec3(0));
joints.knee_l.upd_coordinates(0).setRange([deg2rad(-120), deg2rad(10)])
joints.knee_l.upd_coordinates(0).setName('knee_flexion_l');
model.addJoint(joints.knee_l);

% right shank ~ right foot: ankle_r
joints.ankle_r = PinJoint('ankle_r', bodies.tibia_r, Vec3(0,-0.43,0), Vec3(0), bodies.calcn_r, Vec3(0.04877,0.04195,-0.00792), Vec3(0));
joints.ankle_r.upd_coordinates(0).setRange([deg2rad(-60),deg2rad(60)]);
joints.ankle_r.upd_coordinates(0).setName('ankle_dorsiflexion_r');
model.addJoint(joints.ankle_r);

% left shank ~ left foot: ankle_l
joints.ankle_l = PinJoint('ankle_l', bodies.tibia_l, Vec3(0,-0.43,0), Vec3(0), bodies.calcn_l, Vec3(0.04877,0.04195,0.00792), Vec3(0));
joints.ankle_l.upd_coordinates(0).setRange([deg2rad(-60),deg2rad(60)]);
joints.ankle_l.upd_coordinates(0).setName('ankle_dorsiflexion_l');
model.addJoint(joints.ankle_l);

% pelvis ~ HAT: back
joints.back = WeldJoint('back', bodies.pelvis, Vec3(-0.1007,0.0815,0), Vec3(0), bodies.torso, Vec3(0), Vec3(0)); % weld jonit is a type of "fixed" joint with no coordinate/DOF
model.addJoint(joints.back);

clear spatialTransform

%% Muscles
% hamstrings_r
muscles.hamstrings_r = Millard2012EquilibriumMuscle('hamstrings_r', 2594, 0.109, 0.31, 0); % construction fuction, para: <name>, <max isometric force>, <optimal fiber length>, <tendon slack length>, <pennation angle at optimal>
muscles.hamstrings_r.addNewPathPoint('bifemlh_r-P1', bodies.pelvis, Vec3(-0.12596,-0.10257,0.06944)); % add muscle path point
muscles.hamstrings_r.addNewPathPoint('bifemlh_r-P2', bodies.tibia_r, Vec3(-0.028,-0.02,0.02943));
muscles.hamstrings_r.addNewPathPoint('bifemlh_r-P3', bodies.tibia_r, Vec3(-0.021,-0.04,0.0343));
model.addForce(muscles.hamstrings_r);

% glut_max_r
muscles.glut_max_r = Millard2012EquilibriumMuscle('glut_max_r', 1944, 0.147, 0.127, 0);
muscles.glut_max_r.addNewPathPoint('glut_max2_r-P1', bodies.pelvis, Vec3(-0.1349,0.0176,0.0563));
muscles.glut_max_r.addNewPathPoint('glut_max2_r-P2', bodies.pelvis, Vec3(-0.1376,-0.052,0.0914));
muscles.glut_max_r.addNewPathPoint('glut_max2_r-P3', bodies.femur_r, Vec3(-0.0426,-0.053,0.0293));
muscles.glut_max_r.addNewPathPoint('glut_max2_r-P4', bodies.femur_r, Vec3(-0.0156,-0.1016,0.0419));
model.addForce(muscles.glut_max_r);

% iliopsoas_r
muscles.ilipsoas_r = Millard2012EquilibriumMuscle('ilipsoas_r', 2342, 0.1, 0.163, deg2rad(8));
muscles.ilipsoas_r.addNewPathPoint('psoas_r-P1', bodies.pelvis, Vec3(-0.0647,0.0887,0.0289));
muscles.ilipsoas_r.addNewPathPoint('psoas_r-P2', bodies.pelvis, Vec3(-0.03,-0.01,0.076));
muscles.ilipsoas_r.addNewPathPoint('psoas_r-P4', bodies.femur_r, Vec3(0.033,-0.035,0.0038));
muscles.ilipsoas_r.addNewPathPoint('psoas_r-P5', bodies.femur_r, Vec3(-0.0188,-0.0597,0.0104));
model.addForce(muscles.ilipsoas_r);

% vasti_r
muscles.vasti_r = Millard2012EquilibriumMuscle('vasti_r', 4530, 0.087, 0.136, deg2rad(30));
muscles.vasti_r.addNewPathPoint('vas_int_r-P1', bodies.femur_r, Vec3(0.029,-0.1924,0.031));
muscles.vasti_r.addNewPathPoint('vas_int_r-P2', bodies.femur_r, Vec3(0.0335,-0.2084,0.0285));
muscles.vasti_r.addNewPathPoint('default', bodies.tibia_r, Vec3(0.04,0.025,0.0018)); % 为什么OpenSim中取名叫default？
model.addForce(muscles.vasti_r);

% gastroc_r
muscles.gastroc_r = Millard2012EquilibriumMuscle('gastroc_r', 2241, 0.06, 0.39, deg2rad(17));
muscles.gastroc_r.addNewPathPoint('med_gas_r-P1', bodies.femur_r, Vec3(-0.02,-0.386,-0.024));
muscles.gastroc_r.addNewPathPoint('med_gas_r-P2', bodies.calcn_r, Vec3(0,0.031,-0.0053));
model.addForce(muscles.gastroc_r);

% soleus_r
muscles.soleus_r = Millard2012EquilibriumMuscle('soleus_r', 3549, 0.05, 0.25, deg2rad(25));
muscles.soleus_r.addNewPathPoint('soleus_r-P1', bodies.tibia_r, Vec3(-0.0024,-0.1533,0.0071));
muscles.soleus_r.addNewPathPoint('soleus_r-P2', bodies.calcn_r, Vec3(0,0.031,-0.0053));
model.addForce(muscles.soleus_r);

% tibia_r
muscles.tibia_r = Millard2012EquilibriumMuscle('tibia_r', 1759, 0.098, 0.223, deg2rad(5));
muscles.tibia_r.addNewPathPoint('tib_ant_r-P1', bodies.tibia_r, Vec3(0.0179,-0.1624,0.0115));
muscles.tibia_r.addNewPathPoint('tib_ant_r-P2', bodies.tibia_r, Vec3(0.0329,-0.3951,-0.0177));
muscles.tibia_r.addNewPathPoint('tib_ant_r-P3', bodies.calcn_r, Vec3(0.1166,0.0178,-0.0305));
model.addForce(muscles.tibia_r);

% hamstrings_l
muscles.hamstrings_l = Millard2012EquilibriumMuscle('hamstrings_l', 2594, 0.109, 0.31, 0);
muscles.hamstrings_l.addNewPathPoint('bifemlh_l-P1', bodies.pelvis, Vec3(-0.12596,-0.10257,-0.06944));
muscles.hamstrings_l.addNewPathPoint('bifemlh_l-P2', bodies.tibia_l, Vec3(-0.028,-0.02,-0.02943));
muscles.hamstrings_l.addNewPathPoint('bifemlh_l-P3', bodies.tibia_l, Vec3(-0.021,-0.04,-0.0343));
model.addForce(muscles.hamstrings_l);

% glut_max_l
muscles.glut_max_l = Millard2012EquilibriumMuscle('glut_max_l', 1944, 0.147, 0.127, 0);
muscles.glut_max_l.addNewPathPoint('glut_max2_l-P1', bodies.pelvis, Vec3(-0.1349,0.0176,-0.0563));
muscles.glut_max_l.addNewPathPoint('glut_max2_l-P2', bodies.pelvis, Vec3(-0.1376,-0.052,-0.0914));
muscles.glut_max_l.addNewPathPoint('glut_max2_l-P3', bodies.femur_l, Vec3(-0.0426,-0.053,-0.0293));
muscles.glut_max_l.addNewPathPoint('glut_max2_l-P4', bodies.femur_l, Vec3(-0.0156,-0.1016,-0.0419));
model.addForce(muscles.glut_max_l);

% iliopsoas_l
muscles.ilipsoas_l = Millard2012EquilibriumMuscle('ilipsoas_l', 2342, 0.1, 0.163, deg2rad(8));
muscles.ilipsoas_l.addNewPathPoint('psoas_l-P1', bodies.pelvis, Vec3(-0.0647,0.0887,-0.0289));
muscles.ilipsoas_l.addNewPathPoint('psoas_l-P2', bodies.pelvis, Vec3(-0.03,-0.01,-0.076));
muscles.ilipsoas_l.addNewPathPoint('psoas_l-P4', bodies.femur_l, Vec3(0.033,-0.035,-0.0038));
muscles.ilipsoas_l.addNewPathPoint('psoas_l-P5', bodies.femur_l, Vec3(-0.0188,-0.0597,-0.0104));
model.addForce(muscles.ilipsoas_l);

% vasti_l
muscles.vasti_l = Millard2012EquilibriumMuscle('vasti_l', 4530, 0.087, 0.136, deg2rad(30));
muscles.vasti_l.addNewPathPoint('vas_int_l-P1', bodies.femur_l, Vec3(0.029,-0.1924,-0.031));
muscles.vasti_l.addNewPathPoint('vas_int_l-P2', bodies.femur_l, Vec3(0.0335,-0.2084,-0.0285));
muscles.vasti_l.addNewPathPoint('default', bodies.tibia_l, Vec3(0.04,0.025,-0.0018)); % 为什么OpenSim中取名叫default？
model.addForce(muscles.vasti_l);

% gastroc_l
muscles.gastroc_l = Millard2012EquilibriumMuscle('gastroc_l', 2241, 0.06, 0.39, deg2rad(17));
muscles.gastroc_l.addNewPathPoint('med_gas_l-P1', bodies.femur_l, Vec3(-0.02,-0.386,0.024));
muscles.gastroc_l.addNewPathPoint('med_gas_l-P2', bodies.calcn_l, Vec3(0,0.031,0.0053));
model.addForce(muscles.gastroc_l);

% soleus_l
muscles.soleus_l = Millard2012EquilibriumMuscle('soleus_l', 3549, 0.05, 0.25, deg2rad(25));
muscles.soleus_l.addNewPathPoint('soleus_l-P1', bodies.tibia_l, Vec3(-0.0024,-0.1533,-0.0071));
muscles.soleus_l.addNewPathPoint('soleus_l-P2', bodies.calcn_l, Vec3(0,0.031,0.0053));
model.addForce(muscles.soleus_l);

% tibia_l
muscles.tibia_l = Millard2012EquilibriumMuscle('tibia_l', 1759, 0.098, 0.223, deg2rad(5));
muscles.tibia_l.addNewPathPoint('tib_ant_l-P1', bodies.tibia_l, Vec3(0.0179,-0.1624,-0.0115));
muscles.tibia_l.addNewPathPoint('tib_ant_l-P2', bodies.tibia_l, Vec3(0.0329,-0.3951,0.0177));
muscles.tibia_l.addNewPathPoint('tib_ant_l-P3', bodies.calcn_l, Vec3(0.1166,0.0178,0.0305));
model.addForce(muscles.tibia_l);

%% Contact Geometry
% ground
contactGeo.ground = ContactHalfSpace(Vec3(0,0,0), Vec3(0,0,deg2rad(-90)), model.getGround(), 'ground_contact_geo'); % construction func, para: <contact location>, <contact orientation>, <attached body>, <name>
model.addContactGeometry(contactGeo.ground);
% heel
contactGeo.heel_r = ContactSphere(0.03, Vec3(0.015,0.015,-0.005), bodies.calcn_r, 'heelR_contact_geo'); % construction func, para: <radius>, <contact location>, <attached body>, <name>
model.addContactGeometry(contactGeo.heel_r);
contactGeo.heel_l = ContactSphere(0.03, Vec3(0.015,0.015,0.005), bodies.calcn_l, 'heelL_contact_geo');
model.addContactGeometry(contactGeo.heel_l);
% toe
contactGeo.toe_r = ContactSphere(0.03, Vec3(0.185,0.015,0), bodies.calcn_r, 'toeR_contact_geo');
model.addContactGeometry(contactGeo.toe_r);
contactGeo.toe_l = ContactSphere(0.03, Vec3(0.185,0.015,0), bodies.calcn_l, 'toeL_contact_geo');
model.addContactGeometry(contactGeo.toe_l);

%% Forces
% foot-ground contact force
stiffness =2e6;
dissipation = 1;
staticFriction = 0.9;
dynamicFriction = 0.6;
viscousFriction = 0.6;
transitionVelocity = 0.15;

forces.ground_heel_r = HuntCrossleyForce();
forces.ground_heel_r.setName('heelR_ground_contact_force');
forces.ground_heel_r.addGeometry('heelR_contact_geo'); % every contact force is related to two contact geometries
forces.ground_heel_r.addGeometry('ground_contact_geo');
forces.ground_heel_r.setStiffness(stiffness);
forces.ground_heel_r.setDissipation(dissipation);
forces.ground_heel_r.setStaticFriction(staticFriction);
forces.ground_heel_r.setDynamicFriction(dynamicFriction);
forces.ground_heel_r.setViscousFriction(viscousFriction);
forces.ground_heel_r.setTransitionVelocity(transitionVelocity);
model.addForce(forces.ground_heel_r);

forces.ground_heel_l = HuntCrossleyForce();
forces.ground_heel_l.setName('heelL_ground_contact_force');
forces.ground_heel_l.addGeometry('heelL_contact_geo');
forces.ground_heel_l.addGeometry('ground_contact_geo');
forces.ground_heel_l.setStiffness(stiffness);
forces.ground_heel_l.setDissipation(dissipation);
forces.ground_heel_l.setStaticFriction(staticFriction);
forces.ground_heel_l.setDynamicFriction(dynamicFriction);
forces.ground_heel_l.setViscousFriction(viscousFriction);
forces.ground_heel_l.setTransitionVelocity(transitionVelocity);
model.addForce(forces.ground_heel_l);

forces.ground_toe_r = HuntCrossleyForce();
forces.ground_toe_r.setName('toeR_ground_contact_force');
forces.ground_toe_r.addGeometry('toeR_contact_geo');
forces.ground_toe_r.addGeometry('ground_contact_geo');
forces.ground_toe_r.setStiffness(stiffness);
forces.ground_toe_r.setDissipation(dissipation);
forces.ground_toe_r.setStaticFriction(staticFriction);
forces.ground_toe_r.setDynamicFriction(dynamicFriction);
forces.ground_toe_r.setViscousFriction(viscousFriction);
forces.ground_toe_r.setTransitionVelocity(transitionVelocity);
model.addForce(forces.ground_toe_r);

forces.ground_toe_l = HuntCrossleyForce();
forces.ground_toe_l.setName('toeL_ground_contact_force');
forces.ground_toe_l.addGeometry('toeL_contact_geo');
forces.ground_toe_l.addGeometry('ground_contact_geo');
forces.ground_toe_l.setStiffness(stiffness);
forces.ground_toe_l.setDissipation(dissipation);
forces.ground_toe_l.setStaticFriction(staticFriction);
forces.ground_toe_l.setDynamicFriction(dynamicFriction);
forces.ground_toe_l.setViscousFriction(viscousFriction);
forces.ground_toe_l.setTransitionVelocity(transitionVelocity);
model.addForce(forces.ground_toe_l);

clear stiffness dissipation staticFriction dynamicFriction viscousFriction transitionVelocity

% coordinate limit force: string + 6 double + bool
stiffness = 2;
upperLim = 0;
lowerLim = -120;
damping = 0.2;
transition = 0.01;

forces.knee_limit_r = CoordinateLimitForce('knee_flexion_r', upperLim, stiffness, lowerLim, stiffness, damping, transition, true);
forces.knee_limit_l = CoordinateLimitForce('knee_flexion_l', upperLim, stiffness, lowerLim, stiffness, damping, transition, true);
model.addForce(forces.knee_limit_r);
model.addForce(forces.knee_limit_l);

clear stiffness upperLim lowerLim damping transition
%% Markers
% markers.sternum = Marker('Sternum', bodies.torso, Vec3(0.07,0.3,0)); % construction func, para: <name>, <attached body>, <location>
% model.addMarker(markers.sternum);
% markers.acromium_r = Marker('R.Acromium', bodies.torso, Vec3(-0.03,0.44,0.15));
% model.addMarker(markers.acromium_r);
% markers.acromium_l = Marker('L.Acromium', bodies.torso, Vec3(-0.03,0.44,-0.15));
% model.addMarker(markers.acromium_l);
% markers.top_head = Marker('Top.Head', bodies.torso, Vec3(0.00084,0.657,0));
% model.addMarker(markers.top_head);
% markers.asis_r = Marker('R.ASIS', bodies.pelvis, Vec3(0.02,0.03,0.128));
% model.addMarker(markers.asis_r);
% markers.asis_l = Marker('L.ASIS', bodies.pelvis, Vec3(0.02,0.03,-0.128));
% model.addMarker(markers.asis_l);
% markers.v_sacral = Marker('V.Sacral', bodies.pelvis, Vec3(-0.16,0.04,0));
% model.addMarker(markers.v_sacral);
% markers.thigh_upper_r = Marker('R.Thigh.Upper', bodies.femur_r, Vec3(0.018,-0.2,0.064));
% model.addMarker(markers.thigh_upper_r);
% markers.thigh_front_r = Marker('R.Thigh.Front', bodies.femur_r, Vec3(0.08,-0.25,0.0047));
% model.addMarker(markers.thigh_front_r);
% markers.thigh_rear_r = Marker('R.Thigh.Rear', bodies.femur_r, Vec3(0.01,-0.3,0.06));
% model.addMarker(markers.thigh_rear_r);
% markers.knee_lat_r = Marker('R.Knee.Lat', bodies.femur_r, Vec3(0,-0.404,0.05));
% model.addMarker(markers.knee_lat_r);
% markers.knee_med_r = Marker('R.Knee.Med', bodies.femur_r, Vec3(0,-0.404,-0.05));
% model.addMarker(markers.knee_med_r);
% markers.shank_upper_r = Marker('R.Shank.Upper', bodies.tibia_r, Vec3(0.005,-0.065,0.05));
% model.addMarker(markers.shank_upper_r);
% markers.shank_front_r = Marker('R.Shank.Front', bodies.tibia_r, Vec3(0.05,-0.08,0));
% model.addMarker(markers.shank_front_r);
% markers.shank_rear_r = Marker('R.Shank.Rear', bodies.tibia_r, Vec3(0.02,-0.13,0.05));
% model.addMarker(markers.shank_rear_r);
% markers.ankle_lat_r = Marker('R.Ankle.Lat', bodies.tibia_r, Vec3(-0.005,-0.41,0.053));
% model.addMarker(markers.ankle_lat_r);
% markers.ankle_med_r = Marker('R.Ankle.Med', bodies.tibia_r, Vec3(0.006,-0.3888,-0.038));
% model.addMarker(markers.ankle_med_r);
% markers.heel_r = Marker('R.Heel', bodies.calcn_r, Vec3(-0.02,0.02,0));
% model.addMarker(markers.heel_r);
% markers.midfoot_sup_r = Marker('R.Midfoot.Sup', bodies.calcn_r, Vec3(0.13,0.03,-0.03));
% model.addMarker(markers.midfoot_sup_r);
% markers.midfoot_lat_r = Marker('R.Midfoot.Lat', bodies.calcn_r, Vec3(0.1,0.02,0.04));
% model.addMarker(markers.midfoot_lat_r);
% markers.toe_lat_r = Marker('R.Toe.Lat', bodies.calcn_r, Vec3(0.19,0,0.065));
% model.addMarker(markers.toe_lat_r);
% markers.toe_mid_r = Marker('R.Toe.Mid', bodies.calcn_r, Vec3(0.19,0.005,-0.04));
% model.addMarker(markers.toe_mid_r);
% markers.toe_tip_r = Marker('R.Toe.Tip', bodies.calcn_r, Vec3(0.26,0.005,0));
% model.addMarker(markers.toe_tip_r);
% 
% markers.thigh_upper_l = Marker('L.Thigh.Upper', bodies.femur_l, Vec3(0.018,-0.2,-0.064));
% model.addMarker(markers.thigh_upper_l);
% markers.thigh_front_l = Marker('L.Thigh.Front', bodies.femur_l, Vec3(0.08,-0.25,-0.0047));
% model.addMarker(markers.thigh_front_l);
% markers.thigh_rear_l = Marker('L.Thigh.Rear', bodies.femur_l, Vec3(0.01,-0.3,-0.06));
% model.addMarker(markers.thigh_rear_l);
% markers.knee_lat_l = Marker('L.Knee.Lat', bodies.femur_l, Vec3(0,-0.404,-0.05));
% model.addMarker(markers.knee_lat_l);
% markers.knee_med_l = Marker('L.Knee.Med', bodies.femur_l, Vec3(0,-0.404,0.05));
% model.addMarker(markers.knee_med_l);
% markers.shank_upper_l = Marker('L.Shank.Upper', bodies.tibia_l, Vec3(0.005,-0.065,-0.05));
% model.addMarker(markers.shank_upper_l);
% markers.shank_front_l = Marker('L.Shank.Front', bodies.tibia_l, Vec3(0.05,-0.08,0));
% model.addMarker(markers.shank_front_l);
% markers.shank_rear_l = Marker('L.Shank.Rear', bodies.tibia_l, Vec3(0.02,-0.13,-0.05));
% model.addMarker(markers.shank_rear_l);
% markers.ankle_lat_l = Marker('L.Ankle.Lat', bodies.tibia_l, Vec3(-0.005,-0.41,-0.053));
% model.addMarker(markers.ankle_lat_l);
% markers.ankle_med_l = Marker('L.Ankle.Med', bodies.tibia_l, Vec3(0.006,-0.3888,0.038));
% model.addMarker(markers.ankle_med_l);
% markers.heel_l = Marker('L.Heel', bodies.calcn_l, Vec3(-0.02,0.02,0));
% model.addMarker(markers.heel_l);
% markers.midfoot_sup_l = Marker('L.Midfoot.Sup', bodies.calcn_l, Vec3(0.13,0.03,0.03));
% model.addMarker(markers.midfoot_sup_l);
% markers.midfoot_lat_l = Marker('L.Midfoot.Lat', bodies.calcn_l, Vec3(0.1,0.02,-0.04));
% model.addMarker(markers.midfoot_lat_l);
% markers.toe_lat_l = Marker('L.Toe.Lat', bodies.calcn_l, Vec3(0.19,0,-0.065));
% model.addMarker(markers.toe_lat_l);
% markers.toe_mid_l = Marker('L.Toe.Mid', bodies.calcn_l, Vec3(0.19,0.005,0.04));
% model.addMarker(markers.toe_mid_l);
% markers.toe_tip_l = Marker('L.Toe.Tip', bodies.calcn_l, Vec3(0.26,0.005,0));
% model.addMarker(markers.toe_tip_l);

%% Finalize
model.finalizeConnections(); 
model.print('model\coupled_human-prosthesis_model.osim');

%% Display
model = Model('model\coupled_human-prosthesis_model.osim');
model.setUseVisualizer(true);
state = model.initSystem();
model.updVisualizer().show(state);