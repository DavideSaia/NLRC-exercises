%Solution of Algebric Riccati Equation 
function P_end = ARE(A,B,Q,R,S)

max_iter = 200; 
dt = 1; %delta t 
tf = 1; %final time
threshold = 1e-6; %threshold for early stopping 
n = size(A,1);%dimension of P vector 

%Plot variables
P_history = [];
time_points = []; 

    %Loop to converge to an ARE solution
    for iter = 1:max_iter
    
        tspan = [tf,0];%time vector
        [~,P_v]=ode45(@(t,P) DRE_vec(P,A,B,Q,R),tspan,S(:)); %calculation of P for tf instant
        %Extraction at time t=0 of last P 
        P_end = reshape(P_v(end,:),n,n);
        
        % Save P and the current time for the plot
        P_history = [P_history; P_v(end, :)];
        time_points = [time_points; tf];

        %Check if norm is greater than a threshold 
        if iter>1 && norm(P_end-P_old,'fro') < threshold
            break
        end
        
        P_old = P_end;%save the last P calculated
        tf = tf + dt; %increase time 
    
    end


    iterations = iter; %Number of iteration for convergence 
    
    %Plot data
    P_history_matrix = reshape(P_history, [], n, n);
    
    %Figure with custom size 
    figure('Position', [100 100 800 600]);
    line_styles = {'-', '--', ':', '-.'};
    colors = lines(n^2); % Genera colori distinti
    
    % Plot each elements of P  
    hold on;
    legend_entries = cell(n^2, 1);
    counter = 1;
    for i = 1:n
        for j = 1:n
            plot(time_points, squeeze(P_history_matrix(:, i, j)), ...
                 'LineStyle', line_styles{mod(counter-1,4)+1}, ...
                 'Color', colors(counter,:), ...
                 'LineWidth', 2);
            legend_entries{counter} = sprintf('P_{%d%d}', i, j);
            counter = counter + 1;
        end
    end
    
    
    legend(legend_entries, 'Location', 'bestoutside');
    title_str = sprintf('Evoluzione della matrice P\nSoluzione finale P = \n%s', ...
                        mat2str(P_end, 4));
    title(title_str, 'FontSize', 12);
    
    xlabel('Tempo finale tf');
    ylabel('Valore degli elementi');
    grid on;
    box on;
    hold off;
    
    % Numbers of iteration for the convergence 
    fprintf('Convergenza raggiunta in %d iterazioni\n', iterations);
    fprintf('Soluzione finale P:\n');
    disp(P_end);
end

