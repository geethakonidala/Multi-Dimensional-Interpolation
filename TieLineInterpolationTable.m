classdef TieLineInterpolationTable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        base_pressures;
%         base_coords;
%         F
        composition_interpolants;
    end
    
    methods
        function obj = TieLineInterpolationTable(table_filename)
            %table_filename = 'D:\mkhait\projects\ADGPRS\ADGPRS_original\TestSuite\Denis\OUT_acc_0.txt';
            % interpolate given values
            
            %%%%%%%%%%  TABLE LOAD AND INIT %%%%%%%%%%%%%%
            fid = fopen(table_filename);    % Opening the file

            i = 1;
            %A = fscanf(fid, '%f');
            A = readtable(table_filename, 'Delimiter',' ','ReadVariableNames',false);
            A = table2array(A);
            fclose(fid);
            
            obj.base_pressures = unique(A(:,1));
            %%
            % Interpolating the operator values only at base pressures
%             obj.base_coords = A(:,1:end - 1); % Excluding Pressure,Operator
%             v = A(:, end);                    % Operator Values
%             obj.F = scatteredInterpolant(obj.base_coords, v, 'linear');
%%
            % Creating the separate interpolant function for every base pressure
            % Each interpolant is based on all compositions belonging
            % correspondent pressure value
            %obj.composition_interpolants = scatteredInterpolant.empty(length(obj.base_pressures));
            for i = 1:length(obj.base_pressures)
                filtered_A = A(A(:,1)==obj.base_pressures(i),:); % pick only those rows of A which have obj.base_pressures(i) pressure
                filtered_coords = filtered_A(:,2:end - 1); % only compositions
                filtered_vals = filtered_A(:,end);
                obj.composition_interpolants{i} = scatteredInterpolant(filtered_coords, filtered_vals,'linear');  
            end
     end
            
        function res = interpolate(obj, values)
            dx = 0.0000001;
            %%%%%%%%%%%%% RESULTS %%%%%%%%%%%%%

            res(1) = obj.F(values);
            for a = 1 : obj.n_axes
                values_p = values;
                values_p(a) = values(a) + dx;
                plus = obj.F(values_p);
                values_p(a) = values(a) - dx;
                minus = obj.F(values_p);
                res(a + 1) = (plus - minus) / (2 * dx);
            end
        end
        
        function res = interpolate_v(obj, param_space_coords)    % Values here is [P, Z1, Z2]
            res = zeros(length(param_space_coords),1);
            
            for i = 1:length(obj.base_pressures) - 1
                p_low = obj.base_pressures(i);
                p_high = obj.base_pressures(i + 1);
                bigger = param_space_coords(:,1) > p_low;
                smaller = param_space_coords(:,1) <= p_high;
                filter = logical(bigger .* smaller);
                z_filtered = param_space_coords(filter,2:end); % only compositions
                p_filtered = param_space_coords(filter,1); % only pressure
                
                interpolated_vals_low = obj.composition_interpolants{i}(z_filtered);
                interpolated_vals_high = obj.composition_interpolants{i + 1}(z_filtered);
                weights_high  = (p_filtered - p_low) / (p_high-p_low);
                weights_low  = (p_high - p_filtered) / (p_high-p_low);
                res(filter) = interpolated_vals_low .* weights_low + interpolated_vals_high .* weights_high;
            end
            
        end
    end
    
end

