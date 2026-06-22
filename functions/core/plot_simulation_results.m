function plot_simulation_results(model, modelInfo, dataCategory)
% Name: plot_simulation_results
% Description
% TODO：后续走的时间长了，需要先把数据对步态周期进行归一化

nMuscles = numel(modelInfo.st.muscle.names);
time = modelInfo.st.simInfo.timeSeries;
switch dataCategory
    case 'excitation'
        for i = 0 : nMuscles-1
            figure, plot(time, modelInfo.dy.muscle.exc(i+1, 1:end-modelInfo.st.muscle.delay), 'LineWidth',2);
            grid on, title(['excitation of ' char(model.getMuscles.get(i).getName)], 'Interpreter','none')
            xlabel('time/s')
        end

    % case 'activation'    % DEPRECATED: act field removed
    %     for i = 0 : nMuscles-1
    %         figure, plot(time, modelInfo.dy.muscle.act(i+1, :), 'LineWidth',2);
    %         grid on, title(['activation of ' char(model.getMuscles.get(i).getName)], 'Interpreter','none')
    %         xlabel('time/s')
    %     end

    % case 'excitation&activation'    % DEPRECATED: act field removed
    %     for i = 0 : nMuscles-1
    %         figure, plot(time, modelInfo.dy.muscle.exc(i+1, 1:end-modelInfo.st.muscle.delay), 'LineWidth',2);
    %         hold on, plot(time, modelInfo.dy.muscle.act(i+1, :), '--', 'LineWidth',2);
    %         grid on, title(['excitation & activation of ' char(model.getMuscles.get(i).getName)], 'Interpreter','none')
    %         xlabel('time/s'), legend('excitation', 'activation')
    %     end

    % case 'muscleForce'    % DEPRECATED: fMTU field removed
    %     for i = 0 : nMuscles-1
    %         figure, plot(time, modelInfo.dy.muscle.fMTU(i+1, :), 'LineWidth',2);
    %         grid on, title(['force of ' char(model.getMuscles.get(i).getName)], 'Interpreter','none')
    %         xlabel('time/s')
    %     end

    case 'phase'
        figure, plot(time, modelInfo.dy.phase.r, 'LineWidth',2);
        xlabel('time/s'), grid on, title('gait phase of the right leg');
        ylabel('0-ES; 1-LS; 2-LO; 3-SW; 4-LD')
        figure, plot(time, modelInfo.dy.phase.l, 'LineWidth',2);
        xlabel('time/s'), grid on, title('gait phase of the left leg');
        ylabel('0-ES; 1-LS; 2-LO; 3-SW; 4-LD')

    case 'grf'
        figure, plot(time, modelInfo.grf.fyr, 'LineWidth',2);
        xlabel('time/s'), grid on, title('normal component of grf exerted on the right leg');
        figure, plot(time, modelInfo.grf.fyl, 'LineWidth',2);
        xlabel('time/s'), grid on, title('normal component of grf exerted on the left leg');
        figure, plot(time, modelInfo.grf.fxr, 'LineWidth',2);
        xlabel('time/s'), grid on, title('frictional component of grf exerted on the right leg');
        figure, plot(time, modelInfo.grf.fxl, 'LineWidth',2);
        xlabel('time/s'), grid on, title('frictional component of grf exerted on the left leg');

end

