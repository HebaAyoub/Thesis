clear
clc
% gamma= 25;
m_THD = [10^-5,10^-6,10^-7];

% 
results = [];
for th=1:length(m_THD)

    for i=1:20
        [h(i),s(i)] = get_handover_rate(m_THD(th));
        results = [results; m_THD(th), h(i), s(i)];
    end
end
writematrix(results, 'handover_data.csv');
% 
function [handoverRate,snr_avg]= get_handover_rate(m_THD)
% Gamma_THD = 25; % SNR threshold (𝛤𝑇𝐻𝐷)
% m_THD = -10^-5; % 

% Parameters for handover
n = 10; % Moving average window size
PI = 0.5;
handoverCount = 0;
% Parameters for LTE-A Network
numLTE_BS = 7; % Number of LTE-A macrocells (based on hexagonal grid)
numUsers = 100; % Number of users
numTimeSteps = 100; % Number of time steps for user movement simulation

% LTE-A Parameters
LTE_power_dBm = 46; % dBm
LTE_power_watt = 10.^(LTE_power_dBm/10);
LTE_radius_km = 3; % Cell radius in km
LTE_bandwidth_MHz = 100; % Bandwidth in MHz

% DSRC Parameters
numDSRC_RSU = 5; % Number of DSRC RSUs
DSRC_power_dBm = 13; % dBm
DSRC_power_watt = 10.^(DSRC_power_dBm/10);
DSRC_radius_km = 0.3; % DSRC RSU radius in km
DSRC_bandwidth_MHz = 20; % Bandwidth in MHz

% Convert LTE-A and DSRC cell radius from km to meters
LTE_radius_m = LTE_radius_km * 1000;
DSRC_radius_m = DSRC_radius_km * 1000;

% Area Length for user positions and DSRC/LTE-A layout (in meters)
areaLength = 3 * LTE_radius_m; % Define a large area to accommodate hexagonal grid


% Generate LTE-A BS positions (3-row by 3-column hexagonal grid)
LTE_BS_positions = hexagonal_grid(3, 3, LTE_radius_m);

spacing = areaLength/numDSRC_RSU; % Distance between the RSUs

% Generate positions for RSUs in one line (along x-axis)
x_positions = (0:numDSRC_RSU-1) * spacing; % x-coordinates with constant spacing
y_positions = 4500*ones(1, numDSRC_RSU); % y-coordinates remain the same (e.g., y = 0)

% Combine x and y coordinates into position matrix
DSRC_RSU_positions = [x_positions(:), y_positions(:)];

% Ensure they fit within the area
if max(x_positions) > areaLength
    error('RSUs exceed the area length, increase areaLength or reduce spacing.');
end

% Generate Random Initial Positions for Users
user_positions = areaLength * rand(numUsers, 2);

% Create Figure for Animation
figure;
hold on;
grid on;
axis equal;
axis([0 areaLength 0 areaLength]);
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Heterogeneous Network Layout: LTE-A Hexagonal Grid & DSRC RSUs (User Movement)');

% Plot LTE-A Base Stations and their coverage
for i = 1:size(LTE_BS_positions, 1)
    plot(LTE_BS_positions(i,1), LTE_BS_positions(i,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    viscircles(LTE_BS_positions(i,:), LTE_radius_m, 'LineStyle', '--', 'Color', 'r');
end

% Plot DSRC RSUs and their coverage
for i = 1:size(DSRC_RSU_positions, 1)
    plot(DSRC_RSU_positions(i,1), DSRC_RSU_positions(i,2), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    viscircles(DSRC_RSU_positions(i,:), DSRC_radius_m, 'LineStyle', '--', 'Color', 'b');
end

% Plot Initial User Positions (to be updated dynamically)
user_plot = plot(user_positions(:,1), user_positions(:,2), 'k*', 'MarkerSize', 5);

% SINR and User Association Initialization
noisePower_dBm = -90; % Noise power in dBm

% Mobility Parameters
user_speed_mps = 30; % User speed in meters per second
timeStep = 1; % Time step in seconds

% Array to track the current serving cell for each user
% Initially, users are not connected to any cell (indicated by 0)
currentServingCell = zeros(numUsers, 1); 

% Cell entry log (storing the time steps when users enter a new cell)
cellEntryLog = cell(numUsers, 1); % A cell array where each user will have a log of time steps

% Simulation time 
sim_time = 50; % Total simulation time in seconds
time_steps = sim_time / PI; % Number of measurement steps

% Initialize variables
RSS_t = zeros(1, time_steps); % RSS measurements for current MI
S_t = zeros(1, time_steps); % Moving average of RSS
m = zeros(1, time_steps); % Slope of the moving average
retain_MI = 0; % Indicator for whether to retain the current MI
% MI_selected = 1; % Initially connected to MI 1

 for u = 1:numUsers
        
       % Calculate distances to LTE-A BS
            distances_LTE = sqrt(sum((LTE_BS_positions - user_positions(u,:)).^2, 2));
        
            pl_lte =103.8+20.9*log(distances_LTE.*10^(-3));
        
            % Calculate received powers from LTE-A BS
%             receivedPower_LTE = LTE_power_dBm - pl_lte ;
            receivedPower_LTE_watt = LTE_power_watt*10.^(-pl_lte/10);
        
            % Calculate distances to DSRC RSUs
%             s= (rand * 12) - 6;
            distances_DSRC = sqrt(sum((DSRC_RSU_positions - user_positions(u,:)).^2, 2));
            pl_dsrc=47.9+18*log(distances_DSRC.*10^(-3));
    
            % Calculate received powers from DSRC RSUs
%             receivedPower_DSRC = DSRC_power_dBm - pl_dsrc;
            receivedPower_DSRC_watt = DSRC_power_watt*10.^(-pl_dsrc/10);

        % Choose the best network (LTE-A or DSRC)
        [servingPower_LTE, bestLTE] = max(receivedPower_LTE_watt);
        [servingPower_DSRC, bestDSRC] = max(receivedPower_DSRC_watt);
       
%             
        % Determine which network serves the user (LTE-A or DSRC)
        if servingPower_LTE > servingPower_DSRC;
            % User is served by LTE-A (cell number is LTE BS index)
            newServingCell = bestLTE;
            rss = servingPower_LTE;
        else
            % User is served by DSRC (cell number is negative DSRC RSU index)
            newServingCell = -bestDSRC;
            rss = servingPower_DSRC;
        end
        currentServingCell(u) = newServingCell;
        previousServingCell(u) = currentServingCell(u);
 end
 
 % Main algorithm loop
 rss_log = [];
for t = 2:1:time_steps
%     User Mobility: Random movement
    random_angles = 2 * pi * rand(numUsers, 1); % Random angles of movement
    user_positions(:,1) = user_positions(:,1) + user_speed_mps * timeStep * cos(random_angles);
    user_positions(:,2) = user_positions(:,2) + user_speed_mps * timeStep * sin(random_angles);
    for u = 1:numUsers
        if currentServingCell(u)>0
            currentServingCell_edit = abs(currentServingCell(u));
            distance = sqrt((LTE_BS_positions(currentServingCell_edit,1) - user_positions(u,1))^2 + (LTE_BS_positions(currentServingCell_edit,2) - user_positions(u,2))^2);
            pl_lte =103.8+20.9*log(distance*10^(-3));
            % Calculate received powers from LTE-A BS
%             RSS_t(t) = LTE_power_dBm - pl_lte ;
             RSS_t_watt(t) = LTE_power_watt*10^(-pl_lte/10);
             power_db(t) = 10 * log10(RSS_t_watt(t));

%             RSS_t_watt(t) = 10.^((RSS_t(t))/10);
            Gamma_t(t) = power_db(t)+90;
%             Gamma_t(t) = 10*log10(RSS_t_watt(t)/(LTE_bandwidth_MHz* 10^-9));
            
        else
            currentServingCell_edit = abs(currentServingCell(u));
            distance = sqrt((LTE_BS_positions(currentServingCell_edit,1) - user_positions(u,1))^2 + (LTE_BS_positions(currentServingCell_edit,2) - user_positions(u,2))^2);
            pl_dsrc=47.9+18*log(distance*10^(-3));
            % Calculate received powers from DSRC RSUs
%             RSS_t(t) = DSRC_power_dBm - pl_dsrc;
            RSS_t_watt(t) = DSRC_power_watt*10^(-pl_dsrc/10);
            power_db(t) = 10 * log10(RSS_t_watt(t));

            Gamma_t(t) = power_db(t)+90;

        end

        
        % Step 4: Calculate SNR (Gamma_t)
    if RSS_t_watt(t) > m_THD
        retain_MI = 1; % Maintain connection with current MI

    else
        retain_MI = 0; % Consider handover
    end
   
     % Step 16: Handover decision
    if retain_MI == 0
        handoverCount = handoverCount + 1 % Increment handover count
        all_RSU_distance = sqrt(sum((DSRC_RSU_positions - user_positions(u,:)).^2, 2));
        all_lte_distance = sqrt(sum((LTE_BS_positions - user_positions(u,:)).^2, 2));
        pl_lte =103.8+20.9*log(all_RSU_distance.*10^(-3));
        pl_dsrc = 47.9+18*log(all_RSU_distance.*10^(-3));
%             RSS_t_lte = LTE_power_dBm-pl_lte;
        RSS_t_lte_watt = LTE_power_watt*10.^(-pl_lte/10);
%             RSS_t_rsu = DSRC_power_dBm-pl_dsrc;
        RSS_t_rsu_watt = DSRC_power_watt*10.^(-pl_dsrc/10);

        rss_lte_watt = max(RSS_t_lte_watt);
        rss_rsu_watt = max(RSS_t_rsu_watt);
        if rss_rsu_watt > rss_lte_watt
            currentServingCell(u) = find(RSS_t_rsu_watt == rss_rsu_watt);
            RSS_t_watt(t) = rss_rsu_watt;
            power_db(t) = 10 * log10(RSS_t_watt(t));
            Gamma_t(t) = power_db(t)+90;

%              Gamma_t(t) = 10*log10(RSS_t_watt(t)/(DSRC_bandwidth_MHz*10^6 * 10^-9));

        else
            currentServingCell(u) = find(RSS_t_lte == rss_lte_watt);
            RSS_t_watt(t) = rss_lte_watt;
            power_db(t) = 10 * log10(RSS_t_watt(t));
            Gamma_t(t) = power_db(t)+90;
        end
         % Detect handover
%         if previousServingCell(u) ~= currentServingCell(u)
%             handoverCount = handoverCount + 1 % Increment handover count
% %             disp(['Time ', num2str(t * PI), 's: Handover to MI ', num2str(currentServingCell(u))]);
%         end
%         
        
        % Update previousServingCell after each step
        previousServingCell(u) = currentServingCell(u);
        

    end
       
    end
    
end


% Calculate total handover rate
totalHandovers = sum(handoverCount);
handoverRate = totalHandovers / (numUsers * time_steps)
snr_avg = mean(Gamma_t);
end

function positions = hexagonal_grid(numRows, numCols, radius)
    positions = [];
    for row = 0:numRows-1
        for col = 0:numCols-1
            % Calculate position based on hexagonal grid formula
            x = radius * sqrt(3) * col + mod(row, 2) * (radius * sqrt(3) / 2);
            y = radius * 1.5 * row;
            positions = [positions; x, y];
        end
    end
end