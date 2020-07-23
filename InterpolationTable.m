classdef InterpolationTable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        B;
        axes_npoints;
        axes_min;
        axes_max;
        axes_step;
        axes_step_inv;
        n_axes;
        n_vertex;
        rev_axes_idx;
        wrkspc_array;
        low_high;
    end
    
    methods
        function obj = InterpolationTable(table_filename)
            %table_filename = 'D:\mkhait\projects\ADGPRS\ADGPRS_original\TestSuite\Denis\OUT_acc_0.txt';
            % interpolate given values
            
            %%%%%%%%%%  TABLE LOAD AND INIT %%%%%%%%%%%%%%
            fid = fopen(table_filename);    % Opening the file

            i = 1;
            axis_names = cell(10,1);
            axes_params = [];
            table_dimens = [];
            while 1 
                tline = fgets(fid);
                %axis[i] = sscanf(tline, '%d %f %f');
                [str, count] = sscanf(tline, '%s %d %f %f');
                if (count == 1)
                    break;
                end
                s = char(str(1:end-3));
                axis_names{i} = reshape (s, [1, length(s)]);
                axes_params(i, 1) = str(end - 2);
                axes_params(i, 2) = str(length(str) - 1);
                axes_params(i, 3) = str(length(str));
                table_dimens(i) = axes_params(i, 1);
                i = i + 1;
            end

            A = fscanf(fid, '%f');
            fclose(fid);

            table_dimens = fliplr(table_dimens);
            obj.B = reshape(A, table_dimens);

            obj.axes_npoints = axes_params(:,1).';   % Transposing 
            obj.axes_min = axes_params(:,2).';
            obj.axes_max = axes_params(:,3).';
            obj.axes_step = (obj.axes_max - obj.axes_min) ./ (obj.axes_npoints - 1);
            obj.axes_step_inv = 1 ./ obj.axes_step; % to make 100% same result as C++ code
            obj.n_axes = size(obj.axes_step, 2);
            obj.n_vertex = 2^obj.n_axes;
            obj.rev_axes_idx = fliplr (1:obj.n_axes);
            obj.wrkspc_array = zeros (obj.n_axes + 1, obj.n_vertex);
            obj.low_high = zeros(obj.n_vertex, obj.n_axes);
            for p = 1 : obj.n_vertex
                obj.low_high(p, :) = rem(floor((p-1) ./ 2.^(obj.rev_axes_idx -1)),2) + 1;
            end
        end
        
        function res = interpolate(obj, values)
            %values and same table

            axis_indexes(:, 1) = ceil ((values - obj.axes_min) ./ obj.axes_step);
            axis_indexes(:, 2) = ceil ((values - obj.axes_min) ./ obj.axes_step + 1);
            axis_low = (axis_indexes(:,1).' - 1) .* obj.axes_step + obj.axes_min;
            mult = (values - axis_low) .* obj.axes_step_inv;

            for p = 1 : obj.n_vertex
                low_high_p = obj.low_high(p, :); 
                idx2 = axis_indexes(obj.n_axes,low_high_p(obj.n_axes));
                mu = obj.axes_npoints(obj.n_axes);
                for a = obj.n_axes - 1 : -1: 1
                    idx2 = idx2 + (axis_indexes(a,low_high_p(a)) - 1) * mu;
                    mu = mu * obj.axes_npoints (a);
                end    
                obj.wrkspc_array(1, p) = obj.B(idx2);
                
                %idx = diag(axis_indexes(:,low_high)).';
                %idx = fliplr(idx); % numeration in B is in reversed order (from last, fastest axis, to first)
                %subCell = num2cell(idx);
                %obj.wrkspc_array(1, p) = obj.B(subCell{:});
                %if obj.B(subCell{:}) ~= obj.B(idx2)
                %    print ('error')
                %end    
                    
            end    


            %calculate interpolation
            pwr = obj.n_vertex / 2; % distance between high and low values
            for i = 1 : obj.n_axes
                % calc own derivative
                obj.wrkspc_array(i + 1, 1 : pwr) = (obj.wrkspc_array(1, pwr+1:2*pwr) - obj.wrkspc_array(1, 1:pwr)) .* obj.axes_step_inv (i);
                % update upper derivatives
                obj.wrkspc_array(2 : i, 1 : pwr) = obj.wrkspc_array(2 : i, 1 : pwr) + (obj.wrkspc_array(2 : i, pwr+1:2*pwr) - obj.wrkspc_array(2 : i, 1:pwr)) .* mult(i);
                % calc value
                obj.wrkspc_array(1, 1 : pwr) = obj.wrkspc_array(1, 1 : pwr) + (values(i) - axis_low(i)) .* obj.wrkspc_array(i + 1, 1 : pwr);
                pwr = pwr / 2;
            end    


            %%%%%%%%%%%%% RESULTS %%%%%%%%%%%%%

            res(1) = obj.wrkspc_array(1,1);
            res(2:obj.n_axes+1) = obj.wrkspc_array(2:obj.n_axes+1,1);
        end
        
        function res = interpolate_v(obj, values)
            %values and same table
            nval = size(values, 1);
            axis_mins = repmat(obj.axes_min, nval, 1);
            axis_steps = repmat(obj.axes_step, nval, 1);
            axis_steps_inv = repmat(obj.axes_step_inv, nval, 1);
            axis_indexes(:,:, 1) = ceil ((values - axis_mins) ./ axis_steps); % Rounding to the nearest integer > the value
            axis_indexes(:,:, 2) = ceil ((values - axis_mins) ./ axis_steps + 1);
            axis_low = (axis_indexes(:,:,1) - 1) .* axis_steps + axis_mins;
            mult = (values - axis_low) .* axis_steps_inv;
            
            obj.wrkspc_array = zeros (nval, obj.n_axes + 1, obj.n_vertex);

            for p = 1 : obj.n_vertex
                low_high_p = obj.low_high(p, :); 
                idx2 = axis_indexes(:,obj.n_axes,low_high_p(obj.n_axes));
                %idx2 = repmat (idx2, nval, 1);
                mu = obj.axes_npoints(obj.n_axes);
                for a = obj.n_axes - 1 : -1: 1
                    idx2 = idx2 + (axis_indexes(:,a,low_high_p(a)) - 1) .* mu;
                    mu = mu * obj.axes_npoints (a);
                end    
                obj.wrkspc_array(:, 1, p) = obj.B(idx2);
                
                %idx = diag(axis_indexes(:,low_high)).';
                %idx = fliplr(idx); % numeration in B is in reversed order (from last, fastest axis, to first)
                %subCell = num2cell(idx);
                %obj.wrkspc_array(1, p) = obj.B(subCell{:});
                %if obj.B(subCell{:}) ~= obj.B(idx2)
                %    print ('error')
                %end    
                    
            end    


            %calculate interpolation
            pwr = obj.n_vertex / 2; % distance between high and low values
            for i = 1 : obj.n_axes
                mults_i = repmat (mult(:,i), 1, i - 1, pwr);
                % calc own derivative
                obj.wrkspc_array(:,i + 1, 1 : pwr) = (obj.wrkspc_array(:,1, pwr+1:2*pwr) - obj.wrkspc_array(:,1, 1:pwr)) .* obj.axes_step_inv (i);
                % update upper derivatives
                obj.wrkspc_array(:,2 : i, 1 : pwr) = obj.wrkspc_array(:,2 : i, 1 : pwr) + (obj.wrkspc_array(:,2 : i, pwr+1:2*pwr) - obj.wrkspc_array(:,2 : i, 1:pwr)) .* mults_i;
                % calc value
                difs = values(:,i) - axis_low(:,i);
                difs = repmat(difs, 1, 1, pwr);
                obj.wrkspc_array(:,1, 1 : pwr) = obj.wrkspc_array(:,1, 1 : pwr) + difs .* obj.wrkspc_array(:,i + 1, 1 : pwr);
                pwr = pwr / 2;
            end    


            %%%%%%%%%%%%% RESULTS %%%%%%%%%%%%%

            res(:,1) = obj.wrkspc_array(:,1,1);
            res(:,2:obj.n_axes+1) = obj.wrkspc_array(:,2:obj.n_axes+1,1);
        end
    end
    
end

