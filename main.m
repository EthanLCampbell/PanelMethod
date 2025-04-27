%%-----------------------------------------------------------------------%%
% MATLAB Single-Vortex Constant-Strength Source Panel Method
% Author: Ethan Labianca-Campbell
%%-----------------------------------------------------------------------%%

clear; clc; close all;

%% ASSUMPTIONS
% - inviscid flow
% - icompressible flow (const density, low mach number)
% - irrotational flow (except for single vortex)
% - No 3D effects

%% CORE METHOD
% 1. Represent airfoil surface with a series of straight line segments
% 2. Place a source element of const strength (lambda or sigma) on each
%    panel
% 3. Ass a single vortex of strength Gamma at a fixed point (quarter chord,
%    center) to model overall circulation needed for lift
% 4. Apply boundary condition (flow tangency to airfoil surface)
% 5. Apply kutta condition by requiring tangential velocities above and
%    below trailing edge to be equal and opposite
% 6. Solve system of linear eqns to find unknown source strengths and
%    vortex strength
% 7. Compute velocities and pressure coeffs on surface
% 8. Integrate pressures to find lift (or use kutta-joukowski thrm)

%% STEP 0: Setup & Airfoil Geometry

% --- Flow Conditions ---
Vinf = 1.0;       % Freestream velocity (normalize to 1 for coefficients)
alpha_deg = 15;  % Angle of attack in degrees
alpha_rad = deg2rad(alpha_deg); % Angle of attack in radians7

% --- Airfoil Geometry ---
% Example: Load NACA 0012 coordinates (assuming a file 'naca0012.dat')
% The file should have x in the first column, y in the second,
% ordered from TE over top surface to LE, then bottom surface to TE.

[x_af, y_af, success, message] = readAirfoilData('naca0012.txt');

%{
% Option B: 4-Panel Thin Diamond (Better)
% TE -> UpperMid -> LE -> LowerMid -> TE (Counter-clockwise)
chord = 1.0;
thick = 0.01; % Small thickness
x_af = [chord; chord/2; 0; chord/2; chord];
y_af = [0; thick/2; 0; -thick/2; 0];
fprintf('Using 4-Panel Diamond Airfoil.\n');
success = true; % Mark as success for the rest of the script
%}
N_points = length(x_af);
N_panels = N_points - 1; % since N-1 surfaces between N points

fprintf('Number of points: %d\n', N_points);
fprintf('Number of panels: %d\n', N_panels);

% Plot geometry
%%{
figure;
plot(x_af, y_af, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
axis equal; grid on;
title('Airfoil Geometry and Panels');
xlabel('x/c'); ylabel('y/c');
%}

%% STEP 1: Discretize Geometry (panel generation)

% --- Panel Geometry Calculations ---
X_start = x_af(1:N_panels);     % x-coord of starting point of panels
Y_start = y_af(1:N_panels);     % y-coord of starting point of panels
X_end   = x_af(2:N_points);     % x-coord of ending point of panels
Y_end   = y_af(2:N_points);     % y-coord of ending point of panels

% Panel midpoints (control points)
X_ctrl = (X_start + X_end) / 2;
Y_ctrl = (Y_start + Y_end) / 2;

% Panel lengths
S = sqrt((X_end - X_start).^2 + (Y_end - Y_start).^2);

% Panel angles (angle of panel tangent vector w.r.t positive x-axis)
phi = atan2(Y_end - Y_start, X_end - X_start); % In radians

% Angle of panel normal vector w.r.t positive x-axis
% Normal points outwards from the airfoil
beta = phi - pi/2;
% Ensure beta is within [-pi, pi] or [0, 2pi] if needed (atan2 does this)
beta = atan2(sin(beta), cos(beta)); % Normalize angle

% Tangent and normal vectors (optional, angles phi and beta are often sufficient)
tx = cos(phi);
ty = sin(phi);
nx = cos(beta);
ny = sin(beta);

% --- Sanity Check Normal Vector Direction ---
try % Use try-catch in case centroid calculation fails for weird shapes
    % Calculate approximate centroid
    cent_x = mean(x_af(1:N_panels));
    cent_y = mean(y_af(1:N_panels));

    % Check normal direction for a few panels
    check_indices = [1, round(N_panels/4), round(N_panels/2), round(N_panels*3/4), N_panels];
    %fprintf("Normal Vector Check (Dot product > 0 means likely outward):\n");
    for k = check_indices
        % Vector from centroid to control point
        vec_cx = X_ctrl(k) - cent_x;
        vec_cy = Y_ctrl(k) - cent_y;
        % Dot product with normal vector
        dot_prod = vec_cx * nx(k) + vec_cy * ny(k);
        %fprintf("  Panel %d: Dot Product = %f\n", k, dot_prod);
        if dot_prod <= 0
            warning('Normal vector for panel %d might be pointing inward!', k);
        end
    end
catch ME_centroid
    warning('Could not perform centroid check');
end
% --- End Sanity Check ---

% Add control points to plot
%%{
hold on;
plot(X_ctrl, Y_ctrl, 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
quiver(X_ctrl, Y_ctrl, nx, ny, 0.1, 'r'); % Plot normal vectors (optional)
legend('Airfoil Points', 'Control Points', 'Location', 'northwest');
hold off;
%}

%% STEP 2: Compute Influence Coefficients
% we need velocity induced at control point i by a unit strength source
% panel j (lambda_j = 1) & we need velocity induced by a unit strength
% vortex (gamma=1) placed at a specific loc (origin, 0,0)

% Initialize influence matrices
I = zeros(N_panels, N_panels); % Normal influence of sources
J = zeros(N_panels, N_panels); % Tangential influence of sources
K_n = zeros(N_panels, 1);     % Normal influence of unit vortex
K_t = zeros(N_panels, 1);     % Tangential influence of unit vortex

% Location of the single vortex (at quarter chord)
xv = 0.2501;
yv = 0.0;

for i = 1:N_panels  % Loop through control points (where velocity is evaluated)
    for j = 1:N_panels % Loop through source panels (which induce velocity)

        % Transform control point i into panel j's local coordinate system
        xt = X_ctrl(i) - X_start(j); % x-coord relative to panel j start
        yt = Y_ctrl(i) - Y_start(j); % y-coord relative to panel j start

        % Rotate into panel j's tangent/normal frame
        x_local = xt * cos(-phi(j)) - yt * sin(-phi(j));
        y_local = xt * sin(-phi(j)) + yt * cos(-phi(j));

        % Use a function to get induced velocity (safer)
        [u, v] = calculate_source_velocity(S(j), x_local, y_local);
        if i == 2
            %fprintf('DEBUG i=%d, j=%d: Helper returned u=%.4f, v=%.4f (x_local=%.4f, y_local=%.4f)\n', i, j, u, v, x_local, y_local); % DEBUG LINE
        end
        % Special handling for self-induced term (i == j)
        if i == j
            u = 0.0; % Tangential velocity from panel itself at midpoint
            v = -0.5; % Normal velocity lambda/2 (CHECK SIGN BASED ON NORMAL DEF)
        end
        % Transform induced velocity (u, v) back to global frame
        phi_j = phi(j);
        if i == 2
            %fprintf('DEBUG i=%d, j=%d: Transforming u=%.4f, v=%.4f with phi(j)=%.4f\n', i, j, u, v, phi_j); % DEBUG LINE
        end
        Vx_global = u * cos(phi(j)) - v * sin(phi(j));
        Vy_global = u * sin(phi(j)) + v * cos(phi(j));
        if i == 2
            %fprintf('DEBUG i=%d, j=%d: Result Vx_global=%.4f, Vy_global=%.4f\n', i, j, Vx_global, Vy_global); % DEBUG LINE
        end
        % Project global velocity onto panel i's normal and tangent directions
        % Normal velocity component at panel i due to source panel j
        I(i, j) = Vx_global * nx(i) + Vy_global * ny(i);
        % Tangential velocity component at panel i due to source panel j
        J(i, j) = Vx_global * tx(i) + Vy_global * ty(i);

        % --- START DEBUG PRINT FOR J ---
        %{
        % Print details only for combinations relevant to Kutta condition (i=1 or N, j=1 or N)
        % For N=4, N_panels=4. Check i=1,4 and j=1,4.
        if (i == 1 || i == N_panels) && (j == 1 || j == N_panels)
            fprintf('\n--- Checking J Calculation ---\n');
            fprintf('i=%d, j=%d:\n', i, j);
            fprintf('  Inputs: Vx_global=%.6f, Vy_global=%.6f, tx(%d)=%.6f, ty(%d)=%.6f\n', Vx_global, Vy_global, i, tx(i), i, ty(i));
            J_recalc = Vx_global * tx(i) + Vy_global * ty(i); % Recalculate for verification
            fprintf('  Result: J_recalc=%.6f, Stored J(%d,%d)=%.6f\n', J_recalc, i, j, J(i, j));
        end
        %}
        % --- END DEBUG PRINT FOR J ---
    end % End loop panel j
    
    % --- Calculate influence of unit vortex (Gamma = 1) at (xv, yv) ---
    dx_vortex = X_ctrl(i) - xv;
    dy_vortex = Y_ctrl(i) - yv;
    r_sq_vortex = dx_vortex^2 + dy_vortex^2;

    % --- ADD CHECK FOR ZERO DISTANCE ---
    if r_sq_vortex < 1e-12 % Check if distance is essentially zero
        Vx_vortex = 0;
        Vy_vortex = 0;
    else
        % Velocity induced by unit vortex (Gamma=1) in global frame
        Vx_vortex = +(1 / (2 * pi)) * dy_vortex / r_sq_vortex;
        Vy_vortex = -(1 / (2 * pi)) * dx_vortex / r_sq_vortex;
    end % Add closing 'end' for the if statement

    % Project vortex-induced velocity onto panel i's normal and tangent
    K_n(i) = Vx_vortex * nx(i) + Vy_vortex * ny(i);
    K_t(i) = Vx_vortex * tx(i) + Vy_vortex * ty(i);
    
    % --- START DEBUG PRINT FOR K_t ---
    %{
    % Print details only for panels relevant to Kutta condition (i=1 or N)
    if (i == 1 || i == N_panels)
         fprintf('\n--- Checking K_t Calculation ---\n');
         fprintf('i=%d:\n', i);
         fprintf('  Inputs: Vx_vortex=%.6f, Vy_vortex=%.6f, tx(%d)=%.6f, ty(%d)=%.6f\n', Vx_vortex, Vy_vortex, i, tx(i), i, ty(i));
         Kt_recalc = Vx_vortex * tx(i) + Vy_vortex * ty(i); % Recalculate for verification
         fprintf('  Result: Kt_recalc=%.6f, Stored K_t(%d)=%.6f\n', Kt_recalc, i, K_t(i));
    end
    %}
    % --- END DEBUG PRINT FOR K_t ---
end % End loop panel i


%% STEP 3: Build system of eqns 
% we need N_panels + 1 equations for the N_panels source strengths and the
% 1 vortex strength

% --- Build Linear System Ax = b ---
A = zeros(N_panels + 1, N_panels + 1);
b = zeros(N_panels + 1, 1);

% Fill boundary condition rows (N_panels equations)
% Sum(lambda_j * I(i,j)) + Gamma * K_n(i) = -Vinf * NormalComponent
for i = 1:N_panels
    A(i, 1:N_panels) = I(i, :);         % Source influences on normal velocity
    A(i, N_panels + 1) = K_n(i);         % Vortex influence on normal velocity
    b(i) = -Vinf * cos(beta(i) - alpha_rad); % RHS: Normal component of freestream
end

% Fill Kutta condition row (1 equation)
% Sum(lambda_j * (J(1,j)+J(N,j))) + Gamma * (K_t(1)+K_t(N)) = -Vinf * TangentComponentsSum
A(N_panels + 1, 1:N_panels) = J(1, :) + J(N_panels, :); % Source influences on tangential velocity sum
A(N_panels + 1, N_panels + 1) = K_t(1) + K_t(N_panels); % Vortex influence on tangential velocity sum
b(N_panels + 1) = -Vinf * (cos(phi(1) - alpha_rad) + cos(phi(N_panels) - alpha_rad)); % RHS: Tangential freestream components sum

% --- CHECK MATRIX STEP 3 ---
%{
fprintf('\n--- Verifying Kutta Row Construction (A(N+1,:)) ---\n');
N = N_panels; % Use N for brevity
fprintf('Expected A(N+1, 1) = J(1,1)+J(N,1) = %.6f + %.6f = %.6f\n', J(1,1), J(N,1), J(1,1)+J(N,1));
fprintf('Actual A(N+1, 1)   = %.6f\n', A(N+1, 1));
fprintf('Expected A(N+1, N) = J(1,N)+J(N,N) = %.6f + %.6f = %.6f\n', J(1,N), J(N,N), J(1,N)+J(N,N));
fprintf('Actual A(N+1, N)   = %.6f\n', A(N+1, N));
fprintf('Expected A(N+1, N+1) = K_t(1)+K_t(N) = %.6f + %.6f = %.6f\n', K_t(1), K_t(N), K_t(1)+K_t(N));
fprintf('Actual A(N+1, N+1) = %.6f\n', A(N+1, N+1));
fprintf('Actual b(N+1)      = %.6f\n', b(N+1));
%}
% --- End Verification Section ---

%% STEP 4: Solve the System

% --- Solve the Linear System ---
x_sol = A \ b;

% Extract source strengths and vortex strength
lambda = x_sol(1:N_panels);
Gamma = x_sol(N_panels + 1);

fprintf('Calculated Gamma (Circulation): %f\n', Gamma);

%% STEP 5: Compute surface velocities & pressure coefficients
% --- Calculate Tangential Velocities and Pressure Coefficients ---
Vt = zeros(N_panels, 1);
Cp = zeros(N_panels, 1);

for i = 1:N_panels
    % Ensure lambda is column vector for multiplication
    if size(lambda, 1) ~= N_panels; error('Lambda is wrong size!'); end
    if size(J, 1) ~= N_panels || size(J, 2) ~= N_panels; error('J is wrong size!'); end

    term1 = J(i, :) * lambda;
    term2 = Gamma * K_t(i);
    term3 = Vinf * cos(phi(i) - alpha_rad); % Note: for alpha=0, alpha_rad=0

    Vt_calc_in_loop = term1 + term2 + term3;
    Cp_calc_in_loop = 1.0 - (Vt_calc_in_loop / Vinf)^2;

    % Store calculated values
    Vt(i) = Vt_calc_in_loop;
    Cp(i) = Cp_calc_in_loop;

    % --- DEBUG PRINTING (for selected panels) ---
    %{
    if i==1 || i==round(N_panels/4) || i==round(N_panels/2) || i==round(N_panels*3/4) || i==N_panels
         fprintf('DEBUG Step 5: i=%3d, term1=%.4f, term2=%.4f, term3=%.4f -> Vt(i)=%.4f -> Cp(i)=%.4f\n', ...
                 i, term1, term2, term3, Vt(i), Cp(i));
    end
    %}
    % --- END DEBUG PRINTING ---
end

%% STEP 6: Calculate Aerodynamic Coefficients

% --- Calculate Lift Coefficient ---
% Method 1: Integrating pressure forces (Normal Force Coefficient approx Cl)
CN = -sum(Cp .* S .* sin(phi)); % Normal force coefficient
CA = -sum(Cp .* S .* cos(phi)); % Axial force coefficient (less accurate in potential flow)

Cl_alt_int = sum(-Cp .* S .* (nx * (-sin(alpha_rad)) + ny * cos(alpha_rad)));
Cd_alt_int = sum(-Cp .* S .* (nx * cos(alpha_rad) + ny * sin(alpha_rad)));

% Method 2: Kutta-Joukowski Theorem (Generally more accurate for Cl)
Cl_KJ = 2 * Gamma / (Vinf * 1.0); % Assuming chord c=1.0 for coefficient

fprintf('Lift Coefficient Cl (Pressure Integration): %f\n', Cl_alt_int);
fprintf('Lift Coefficient Cl (Kutta-Joukowski):    %f\n', Cl_KJ);
fprintf('Drag Coefficient Cd (Pressure Integration): %f (a measure of numerical inaccuracy - should be 0)\n', Cd_alt_int); % Expect near zero

% --- Calculate Moment Coefficient (about quarter-chord, c/4) ---
xc_ref = 0.25; % Quarter chord reference point (x-coordinate)
yc_ref = 0.0;  % Quarter chord reference point (y-coordinate, assumes symmetric section)

Cm_c4 = sum(Cp .* S .* ( -(X_ctrl - xc_ref).*ny + (Y_ctrl - yc_ref).*nx ));
fprintf('Moment Coefficient Cm_c/4 (Pressure Integration): %f\n', Cm_c4);

%% STEP 7: Plot results

% --- Plot Cp Distribution ---
figure;
% Plot Cp vs x, distinguishing upper and lower surfaces
mid_point_index = round(N_panels / 2); % Approximate index separating upper/lower

plot(X_ctrl(1:mid_point_index), Cp(1:mid_point_index), 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Lower Surface');
hold on;
plot(X_ctrl(mid_point_index+1:N_panels), Cp(mid_point_index+1:N_panels), 'rs-', 'LineWidth', 1.5, 'DisplayName', 'Upper Surface');

set(gca, 'YDir','reverse'); % Standard Cp plot convention
grid on;
xlabel('x/c');
ylabel('C_p');
title(['Pressure Coefficient Distribution, \alpha = ', num2str(alpha_deg), '^{\circ}']);
legend('Location', 'best');
hold off;

% --- Plot Vt Distribution ---
figure;
upper_indices = find(Y_ctrl > 0);
lower_indices = find(Y_ctrl <= 0);
[~, upper_sort_idx] = sort(X_ctrl(upper_indices));
[~, lower_sort_idx] = sort(X_ctrl(lower_indices));

plot(X_ctrl(lower_indices(lower_sort_idx)), Vt(lower_indices(lower_sort_idx))/Vinf, 'bo-', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', 'Lower Surface');
hold on;
plot(X_ctrl(upper_indices(upper_sort_idx)), Vt(upper_indices(upper_sort_idx))/Vinf, 'rs-', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', 'Upper Surface');
grid on;
xlabel('x/c');
ylabel('V_t / V_{\infty}');
title(['Tangential Velocity Distribution, \alpha = ', num2str(alpha_deg), '^{\circ}, N = ', num2str(N_panels)]);
legend('Location', 'best');
hold off;


%% HELPER FUNCTIONS (because coding is hard) 

function [u, v] = calculate_source_velocity(PanelLength, x, y)
    % Calculates velocity (u,v) induced by a constant strength
    % source panel (lambda=1) of length PanelLength at local point (x,y).
    % Panel lies on x-axis from x=0 to x=PanelLength.

    epsilon = 1e-8; % Small number to avoid singularity

    if abs(y) < epsilon
        % Point is on or very near the panel's axis
        % Handle carefully based on x position relative to panel ends (0, PanelLength)
        if x > epsilon && x < PanelLength - epsilon % Point is ON the panel
            u = 0; % Tangential velocity is tricky, often taken as 0 for midpoint approx.
            v = 0.5 * sign(y + epsilon); % +/- 0.5 depending on which side y approaches from
                                        % Using sign(y+epsilon) avoids sign(0) issues
                                        % and defaults to +0.5 if y is exactly 0
        elseif x <= epsilon % Point is off the start end
            theta1 = atan2(y, x);
            theta2 = atan2(y, x - PanelLength);
            % log term: Be careful if term1 or term2 is zero
            term1_sq = max(x^2 + y^2, epsilon^2); % Avoid log(0)
            term2_sq = max((x - PanelLength)^2 + y^2, epsilon^2);
             % REMOVED extra /2 factor here
            u = (1 / (4 * pi)) * log(term1_sq / term2_sq);
            v = (1 / (2 * pi)) * (theta2 - theta1);
        else % Point is off the end end (x >= PanelLength - epsilon)
            theta1 = atan2(y, x);
            theta2 = atan2(y, x - PanelLength);
             % REMOVED extra /2 factor here
            term1_sq = max(x^2 + y^2, epsilon^2);
            term2_sq = max((x - PanelLength)^2 + y^2, epsilon^2);
            u = (1 / (4 * pi)) * log(term1_sq / term2_sq);
            v = (1 / (2 * pi)) * (theta2 - theta1);
        end
    else
        % Standard case: Point is off the panel axis
        theta1 = atan2(y, x);
        theta2 = atan2(y, x - PanelLength);
         % REMOVED extra /2 factor here
        term1_sq = x^2 + y^2;
        term2_sq = (x - PanelLength)^2 + y^2;
        u = (1 / (4 * pi)) * log(term1_sq / term2_sq);
        v = (1 / (2 * pi)) * (theta2 - theta1);
    end
end


function [x_af, y_af, success, message] = readAirfoilData(filename)
    %readAirfoilData Reads airfoil coordinates from a specific text format.
    %   [x_af, y_af, success, message] = readAirfoilData(filename) reads the
    %   airfoil data from the specified file. The file is expected to have:
    %   - Line 1: Header text (ignored)
    %   - Line 2: Two numbers (Np_upper, Np_lower) - points per surface
    %   - Line 3 onwards: Np_upper lines of (x,y) for upper surface (LE to TE)
    %   - Next Np_lower lines of (x,y) for lower surface (LE to TE)
    %
    %   Outputs:
    %   x_af      - Column vector of x-coordinates ordered counter-clockwise
    %               from a closed trailing edge. Empty if failed.
    %   y_af      - Column vector of y-coordinates ordered counter-clockwise
    %               from a closed trailing edge. Empty if failed.
    %   success   - Logical true if successful, false otherwise.
    %   message   - String containing status or error information.
    %
    %   The function reorganizes the data, reverses the upper surface order,
    %   averages the trailing edge y-coordinates to create a closed TE, and
    %   ensures counter-clockwise ordering suitable for panel methods.
    
    % Initialize outputs
    x_af = [];
    y_af = [];
    success = false;
    message = '';
    
    % --- Input Validation ---
    if nargin < 1 || ~ischar(filename) && ~isstring(filename) || isempty(filename)
        message = 'Error: Invalid or empty filename provided.';
        return;
    end
    if ~isfile(filename)
         message = sprintf('Error: File not found: %s', filename);
         return;
    end
    
    % --- File Reading and Processing ---
    fid = -1; % Initialize file ID
    try
        fid = fopen(filename, 'r');
        if fid == -1
            % This condition is technically redundant due to isfile check above,
            % but kept for robustness in case of permission issues etc.
            error('readAirfoilData:FileOpenError', 'Cannot open file: %s', filename);
        end
    
        % Read header line 1 (string) and ignore
        fgetl(fid);
    
        % Read header line 2 (number of points)
        header2 = textscan(fid, '%f %f', 1);
        if isempty(header2{1}) || isempty(header2{2})
             error('readAirfoilData:HeaderError', 'Could not read point counts from header line 2.');
        end
        Np_upper = header2{1};
        Np_lower = header2{2};
        if Np_upper < 2 || Np_lower < 2
            error('readAirfoilData:HeaderError', 'Invalid point counts in header: N_upper=%d, N_lower=%d.', Np_upper, Np_lower);
        end
        fprintf('Expecting %d upper and %d lower surface points from file.\n', Np_upper, Np_lower);
    
        % Read the coordinate data based on header counts
        data = cell2mat(textscan(fid, '%f %f', Np_upper + Np_lower));
        fclose(fid); % Close file now that data is read
    
        % Validate data read
        if size(data, 1) ~= (Np_upper + Np_lower) || size(data, 2) ~= 2
            error('readAirfoilData:DataReadError', 'Expected %d total points, but read %d.', Np_upper + Np_lower, size(data, 1));
        end
    
        % --- Separate and Reshape Data ---
        % Extract upper surface data (as read, likely LE to TE)
        xu_read = data(1:Np_upper, 1);
        yu_read = data(1:Np_upper, 2);
    
        % Extract lower surface data (as read, likely LE to TE)
        xl_read = data(Np_upper+1 : end, 1);
        yl_read = data(Np_upper+1 : end, 2);
    
        % Reverse the upper surface data to go from TE to LE
        xu_rev = flipud(xu_read);
        yu_rev = flipud(yu_read);
    
        % Create a closed trailing edge by averaging the y-coordinates
        % TE is the first point of reversed upper data & last point of lower data
        x_TE = xu_rev(1); % Should be ~1.0
        if abs(x_TE - 1.0) > 0.01 || abs(xl_read(end) - 1.0) > 0.01
            warning('Trailing edge x-coordinate is not close to 1.0.');
        end
        y_TE_avg = (yu_rev(1) + yl_read(end)) / 2.0;
    
        % Assemble the final coordinate arrays in counter-clockwise order
        % Format: [TE_avg; Upper (TEex -> LEinc); Lower (LEex -> TEex); TE_avg]
        % xu_rev(2:end) -> Point after TE down to LE (inclusive)
        % xl_read(2:end-1) -> Point after LE down to point before TE
        x_af = [x_TE; xu_rev(2:end); xl_read(2:end-1); x_TE];
        y_af = [y_TE_avg; yu_rev(2:end); yl_read(2:end-1); y_TE_avg];
    
        % Ensure output is column vector
        x_af = x_af(:);
        y_af = y_af(:);
    
        success = true;
        message = sprintf('Successfully processed %s', filename);
        fprintf('%s\n', message);
    
    catch ME % Catch errors during file operations or processing
        if fid ~= -1
            fclose(fid); % Ensure file is closed even if error occurred
        end
        message = sprintf('Error processing file %s: %s (Line %d)', filename, ME.message, ME.stack(1).line);
        fprintf('%s\n', message);
        % Return empty arrays and success=false
        x_af = [];
        y_af = [];
        success = false;
    end
    
end % End of function readAirfoilData

