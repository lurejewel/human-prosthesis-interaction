function plot_simulation_results(model, modelInfo, dataCategory)
% Name: plot_simulation_results
% Description
% TODO：后续走的时间长了，需要先把数据对步态周期进行归一化

muscleNum = height(modelInfo.muscleExcitations);
time = modelInfo.time;
switch dataCategory
    case 'excitation'
        for muscleIndex = 0 : muscleNum-1
            figure, plot(time, modelInfo.muscleExcitations(muscleIndex+1, 1:end-modelInfo.muscleReflexDelay), 'LineWidth',2);
            grid on, title(['excitation of ' char(model.getMuscles.get(muscleIndex).getName)], 'Interpreter','none')
            xlabel('time/s')
        end

    case 'activation'
        for muscleIndex = 0 : muscleNum-1
            figure, plot(time, modelInfo.muscleActivations(muscleIndex+1, :), 'LineWidth',2);
            grid on, title(['activation of ' char(model.getMuscles.get(muscleIndex).getName)], 'Interpreter','none')
            xlabel('time/s')
        end

    case 'excitation&activation'
        for muscleIndex = 0 : muscleNum-1
            figure, plot(time, modelInfo.muscleExcitations(muscleIndex+1, 1:end-modelInfo.muscleReflexDelay), 'LineWidth',2);
            hold on, plot(time, modelInfo.muscleActivations(muscleIndex+1, :), '--', 'LineWidth',2);
            grid on, title(['excitation & activation of ' char(model.getMuscles.get(muscleIndex).getName)], 'Interpreter','none')
            xlabel('time/s'), legend('excitation', 'activation')
        end

    case 'muscleForce'
        for muscleIndex = 0 : muscleNum-1
            figure, plot(time, modelInfo.muscleForces(muscleIndex+1, :), 'LineWidth',2);
            grid on, title(['force of ' char(model.getMuscles.get(muscleIndex).getName)], 'Interpreter','none')
            xlabel('time/s')
        end

    case 'phase'
        figure, plot(time, modelInfo.phaseR, 'LineWidth',2);
        xlabel('time/s'), grid on, title('gait phase of the right leg');
        figure, plot(time, modelInfo.phaseL, 'LineWidth',2);
        xlabel('time/s'), grid on, title('gait phase of the left leg');

    case 'grf'
        figure, plot(time, modelInfo.grf.normalR, 'LineWidth',2);
        xlabel('time/s'), grid on, title('normal component of the ground reaction force exerted on the right leg');
        figure, plot(time, modelInfo.grf.normalL, 'LineWidth',2);
        xlabel('time/s'), grid on, title('normal component of the ground reaction force exerted on the left leg');
        figure, plot(time, modelInfo.grf.frictionR, 'LineWidth',2);
        xlabel('time/s'), grid on, title('frictional component of the ground reaction force exerted on the right leg');
        figure, plot(time, modelInfo.grf.frictionL, 'LineWidth',2);
        xlabel('time/s'), grid on, title('frictional component of the ground reaction force exerted on the left leg');

end

