% -------------------------------------------------------------------------
% Name: read_trc_file
% Description: visualize the .trc file (marker movement from motion capture
% data).
% Note 1: OpenSim API does not provide direct tool related with the .trc
% file, so visualization is done via MATLAB plotting tools.
% Note 2: .trc file is not related to model file.
% Note 3: the xyz definition in MATLAB plotting tool is not consistent with
% that in .trc file, so adjustment is needed.
% -------------------------------------------------------------------------
addpath('assets\', 'model\');
import org.opensim.modeling.*

%% load trc file
trcPath = 'assets\subject01_walk1.trc';
trc = MarkerData(trcPath); % read .trc file
rate    = trc.getCameraRate(); % camera rate; time between the rate
markerNum = trc.getNumMarkers(); % number of markers

%% show marker data in every frame
% init matlab plotting tool
figure, hold on, grid on, axis equal, view(45,30)
xlabel('X/mm'), ylabel('Z/mm'), zlabel('Y/mm');
title(['motion data read from ' trcPath], 'Interpreter', 'none');
markerPlot = scatter3(nan(markerNum,1),nan(markerNum,1),nan(markerNum,1),'filled');
% flow: get marker number -> get marker position -> visualize
for frame = 0 : trc.getNumFrames() - 1 % for each frame
    currentFrame = trc.getFrame(frame); % get the frame
    currentPos   = zeros(markerNum, 3); % init position vector for all markers at the frame
    for i = 0 : markerNum-1 % for each marker
        tempPos = currentFrame.getMarker(i); % get position (vec3)
        currentPos(i+1,:) = [tempPos.get(0), tempPos.get(2), tempPos.get(1)]; % turn vec3 to matrix, and adjust the axes' order so that it matches those in motion capture system
    end
    % assign the markers' position data to the scatter handle
    markerPlot.XData = currentPos(:,1);
    markerPlot.YData = currentPos(:,2);
    markerPlot.ZData = currentPos(:,3);
    drawnow; % fresh the figure
    pause(1/rate); % set appropriate pause time to simulate the camera rate
end