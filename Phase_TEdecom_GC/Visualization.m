%     %% 6. plot redundant information change in TE Decomposition
%     figure
%     plot(g_red,'-*k','LineWidth', 2);hold on
%     plot(2:cc_num+1,g_red_surr(:,condition_index),'-or');hold on
%     plot(p_x (sig_idx_red), g_red(sig_idx_red), 'g*', 'MarkerSize', 8, 'LineWidth', 2);
%     % xlim([0 3])
%     for i = 1:length(sig_idx_red)
%         x_val_red = p_x(sig_idx_red(i));
%         y_val_red = g_red(sig_idx_red(i));
%         text(x_val_red, y_val_red, sprintf('%d', Channels_idx(drivers_red(condition_index(sig_idx_red(i)-1)))), ...
%             'VerticalAlignment', 'top', ...
%             'HorizontalAlignment', 'center', ...
%             'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
%     end
%     xlabel('Number of Conditioned Variables')
%     ylabel(sprintf('TE values with order %d and lag %d', best_p, best_d))
%     title(['Redundant information flow test at time window', num2str(window_index)])
%     legend('Raw signal','Surrogate Signal')
%     %% plot synergistic information change in TE Decomposition
%     figure
%     plot(g_syn,'-*k','LineWidth', 2);hold on
%     plot(2:cc_num+1,g_syn_surr(:,condition_index),'-or');
%     plot(p_x (sig_idx_syn), g_syn(sig_idx_syn), 'g*', 'MarkerSize', 8, 'LineWidth', 2);
%     % xlim([0 3])
%     for i = 1:length(sig_idx_syn)
%         x_val_syn = p_x(sig_idx_syn(i));
%         y_val_syn = g_syn(sig_idx_syn(i));  
%         text(x_val_syn , y_val_syn, sprintf('%d', Channels_idx(drivers_syn(condition_index(sig_idx_syn(i)-1)))), ...
%             'VerticalAlignment', 'top', ...
%             'HorizontalAlignment', 'center', ...
%             'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
%     end
%     xlabel('Number of Conditioned Variables')
%     ylabel(sprintf('TE values with order %d and lag %d', best_p, best_d))
%     title(['Synergistic information flow test at time window', num2str(window_index)])
%     legend('Raw signal','Surrogate Signal')
% end
% 
