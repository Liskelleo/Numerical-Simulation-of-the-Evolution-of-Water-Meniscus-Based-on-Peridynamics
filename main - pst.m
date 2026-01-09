% 加载数据
data = load('C:\Users\ASUS\Desktop\soil_board.dat');
indices = data(:, 1); % 第一列是索引
x_0 = data(:, 2); % 第二列是x坐标
y_0 = data(:, 3); % 第三列是y坐标
labels = data(:, 4); % 第四列是标签

data1 = load('C:\Users\ASUS\Desktop\air_backround.dat');
x1 = data1(:, 2);
y1 = data1(:, 3);

% 设置基本参数
deltax = 1e-2; % 离散步长0.00125
r = 3.015 * deltax; % 搜索半径3.6*deltax
dt = 5e-7; % 时间步长1e-5
num_steps = 30000; % 时间步数3e4
F_B = [0, 0]; % 忽略体力，如需考虑要换成函数去调用
v_wall = [0, 0]; % 固体墙速度矢量
a_wall = [0, 0]; % 固体墙加速度矢量
v_fmax = 1; % 最大流体速度，这里取为1；若对于溃坝问题，等于sqrt(gH)，g为重力加速度
c_f =12 * v_fmax; % 流体声速，与最大流体速度有关，这里之前取的是10，貌似这个参数影响很大
gamma_f = 1; % 流体压强常数
density_f0 = 1; % 初始流体密度
reference_pressure = density_f0 * (c_f)^2 / gamma_f; % 参考压强
theta_eq = 5 * pi / 180; % 用户定义的theta_eq，示例值为5度
beta = 1; % 表面张力系数
miu = 0.2; % 动力粘度系数
C_pst = 1e-4;
% 初始化物理量
neighbors = cell(length(x_0), 1);
density = zeros(length(x_0), 1);
density(labels == 1) = density_f0; % 初始化流体粒子的密度
pressure = zeros(length(x_0), 1);
pressure(labels == 1) = reference_pressure * ((density(labels == 1) / density_f0).^gamma_f - 1); % 初始化流体粒子的压强
displacement = zeros(length(x_0), 2);pst_displacement = zeros(length(x_0), 2); % [新增] 初始化PST位移数组
velocity = zeros(length(x_0), 2);
velocity_magnitude = zeros(length(x_0), 1);
acceleration_old = zeros(length(x_0), 2);
acceleration_now = zeros(length(x_0), 2);


% 预设改进的高斯权函数，平滑长度h=1.2*deltax
w_0 = @(D, deltax) (exp(-(D ./ (1.2 * deltax)).^2) - exp(-9)) / (1.44 * pi * (1 - 10 * exp(-9))); 
% 匿名函数形式的 vfrac_s
calculate_vfrac_s = @(D, deltax) (ones(size(D)) );

% calculate_vfrac_s = @(D, deltax) (...
%     (3.015 * deltax - D <= 0.5 * deltax & 3.015 * deltax - D >= 0) .* (0.5 + (3.015 * deltax - D) ./ deltax) + ...
%     (3.015 * deltax - D >= 0.5 * deltax) .* 1 + ...
%     (3.015 * deltax - D < 0) .* 0);

% 以下逻辑有待讨论
% calculate_vfrac_s = @(D, deltax) (...
%     (abs(D - 3.015 * deltax) < 0.5 * deltax) .* (0.5 + (3.015 * deltax - D) ./ deltax) + ...
%     (3.015 * deltax - D >= 0.5 * deltax) .* 1 + ...
%     (abs(D - 3.015 * deltax) >= 0.5 * deltax & 3.015 * deltax - D < 0.5 * deltax) .* 0);

% 主循环
for step = 1:num_steps

    disp(['Time step ', num2str(step), ';']);
    body_force = zeros(length(x_0), 2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % 加载项初始化 % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        
        if step ~=1
            % 加载数据
            data = load('C:\Users\ASUS\Desktop\soil_board_now.dat');
            indices = data(:, 1); % 第一列是索引
            x = data(:, 2); % 第二列是x坐标
            y = data(:, 3); % 第三列是y坐标
            labels = data(:, 4); % 第四列是标签
        elseif step == 1
            x = x_0; y=y_0;
        end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % 寻找近场粒子 % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        
        % 调用子程序计算邻近点集，但是为了节约时间，不选择在for循环里插入外部函数
        % neighbors = element_form_continum(x, y, indices, labels, r);
        
        % 使用 rangesearch 函数获取邻近点集
        [neighborIndices, ~] = rangesearch([x, y], [x, y], r);
        % 使用索引存储邻近点集
        for i = 1:length(x)
            current_indices = neighborIndices{i};
            current_labels = labels(current_indices);

            if labels(i) == 2
                % 剔除所有索引为2的邻近粒子
                current_indices = current_indices(current_labels ~= 2);
                current_labels = current_labels(current_labels ~= 2);
            end

            if isempty(current_indices) || all(current_labels == 2)
                continue; % 如果周围全是索引为2的邻近粒子，跳过这个中心粒子
            end

            neighbors{i} = indices(current_indices); % 存储索引而不是坐标
        end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % 设置固体边界 % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        % 利用固体边界条件，求解虚拟层速度矢量
        for i = 1:length(x)
            if labels(i) == 2 % 只对固体虚拟层粒子进行求解
                current_neighbors = neighbors{i};
                if isempty(current_neighbors)
                    continue; 
                end
%                 disp(length(current_neighbors));
                current_neighbors = current_neighbors(current_neighbors ~= indices(i)); % 排除自身点
%                 disp(length(current_neighbors));
%                 disp('******************************************');

                dx = x(current_neighbors) - x(i);
                dy = y(current_neighbors) - y(i);
                D = sqrt(dx.^2 + dy.^2);
 
                % 计算 vfrac_s 和 w_0
                vfrac_s = calculate_vfrac_s(D, deltax);
%                 disp(vfrac_s);
%                 if abs(D - 3.015 * deltax) < 0.5 * deltax
%                     vfrac_s = 0.5 + (3.015 * deltax - D) ./ deltax;
%                 elseif 3.015 * deltax - D >= 0.5 * deltax
%                     vfrac_s = ones(size(D));
%                 else
%                     vfrac_s = zeros(size(D));
%                 end
                weight = w_0(D, deltax);

                % 计算邻近点速度的加权和
                v_j = velocity(current_neighbors, :);
%                 disp(size(D));disp(size(weight));disp(size(v_j));
                weighted_sum_velocity = sum([vfrac_s, vfrac_s] .* [weight, weight] .* v_j, 1);

                % 更新速度矢量
                velocity(i, :) = 2 * v_wall - weighted_sum_velocity;
                
                % 计算邻近点密度和压强的加权和
                density_j = density(current_neighbors);
                pressure_j = pressure(current_neighbors);
                F_B_a_wall = repmat(F_B,length(current_neighbors),1) - [density_j, density_j] .* repmat(a_wall,length(current_neighbors),1);
                weighted_sum_pressure = sum(vfrac_s .* (pressure_j + (F_B_a_wall(:, 1) .* dx + F_B_a_wall(:, 2) .* dy))  .* weight, 1);
%                 disp(size(weighted_sum_pressure));

                % 更新压强
                pressure(i) = weighted_sum_pressure / sum(vfrac_s .* weight);
                
                % 更新密度
                density(i) = density_f0 * (pressure(i) / reference_pressure)^(1 / gamma_f);
            end
        end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % 计算表面张力 % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

            x_filtered = x(labels == 1);
            y_filtered = y(labels == 1);
            indices_filtered = indices(labels == 1);
%             k = boundary(x_filtered, y_filtered,0.9);
            % 找到边界 (使用 Alpha Shape)
            shp = alphaShape(x, y, 1.5*deltax); % 适当调整 alpha 值
            k = boundaryFacets(shp);
            k = k(:,1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%             if step ==1
%                 % 初始化空数组
%                 xA = [];yA = [];
%                 xB = [];yB = [];
%                 
%                 % 遍历所有点
%                 for i = 1:length(k)
%                         % 将该点的横坐标和纵坐标分别添加到 xA 和 yA 中
%                         xA(end + 1) = x(k(i));;  % 添加到 xA
%                         yA(end + 1) = y(k(i));  % 添加到 yA
%                 end
%                 
%                 % 遍历所有点
%                 for i = 1:length(x)
%                     % 如果当前点的标签为2
%                     if labels(i) == 2
%                         % 将该点的横坐标和纵坐标分别添加到 xB 和 yB 中
%                         xB(end + 1) = x(i);  % 添加到 xB
%                         yB(end + 1) = y(i);  % 添加到 yB
%                     end
%                 end
%                 
%                 % 组合点集A和点集B的坐标
%                 A = [xA', yA'];  % 将xA和yA组合成二维点集矩阵
%                 B = [xB', yB'];  % 将xB和yB组合成二维点集矩阵
%                 
%                 % 使用 ismember 判断点集A中的点是否属于点集B
%                 isMember = ismember(A, B, 'rows');
%                 
%                 % 提取点集A中不属于点集B的点
%                 notInB = A(~isMember, :);
%                 
%                 % 提取这些点的横纵坐标
%                 xA_notInB = notInB(:, 1);  % 不属于B的点的横坐标
%                 yA_notInB = notInB(:, 2);  % 不属于B的点的纵坐标
%                 
% %                 % 显示结果
% %                 disp('点集A中不属于点集B的点：');
% %                 disp(notInB);
% %                 disp('横坐标：');
% %                 disp(xA_notInB);
% %                 disp('纵坐标：');
% %                 disp(yA_notInB);
%                 
%                 % 根据 y 坐标分为两部分
%                 y_greater_than_0_5 = notInB(notInB(:, 2) > 0.5, :);  % y > 0.5 的点
%                 y_less_than_or_equal_0_5 = notInB(notInB(:, 2) <= 0.5, :);  % y ≤ 0.5 的点
%                 
%                 % 保存为 .dat 文件
%                 filename1 = 'points_y_greater_than_0_5.dat';
%                 filename2 = 'points_y_less_than_or_equal_0_5.dat';
%                 
%                 % 写入文件
%                 dlmwrite(filename1, y_greater_than_0_5, 'delimiter', '\t', 'precision', 6);
%                 dlmwrite(filename2, y_less_than_or_equal_0_5, 'delimiter', '\t', 'precision', 6);
%                 
%                 % 输出文件名
%                 disp(['文件保存完成：', filename1]);
%                 disp(['文件保存完成：', filename2]);
% 
%             end
%                         break



%             if step == 110
%                 % 初始化空数组
%                 A = [];  % 用于存储点集A的坐标和对应的第五列值
%                 xB = []; yB = [];
% 
%                 % 遍历所有点，并将坐标和对应的data第五列值一起存储到A中
%                 for i = 1:length(k)
%                     % 提取点的坐标和对应的第五列值
%                     x_val = x(k(i));
%                     y_val = y(k(i));
%                     data5_val = data(k(i), 5);  % 提取data的第五列值
% 
%                     % 将坐标和第五列值作为一行加入到A中
%                     A(end + 1, :) = [x_val, y_val, data5_val];
%                 end
% 
%                 % 遍历所有点，提取标签为2的点的坐标到xB和yB
%                 for i = 1:length(x)
%                     if labels(i) == 2
%                         xB(end + 1) = x(i);
%                         yB(end + 1) = y(i);
%                     end
%                 end
% 
%                 % 组合点集B的坐标
%                 B = [xB', yB'];  % 将xB和yB组合成二维点集矩阵
% 
%                 % 使用 ismember 判断点集A中的点是否属于点集B
%                 % 注意：这里只比较前两列（x和y坐标）
%                 isMember = ismember(A(:, 1:2), B, 'rows');
% 
%                 % 提取点集A中不属于点集B的点（包括坐标和对应的第五列值）
%                 notInB = A(~isMember, :);
% 
%                 % 根据 y 坐标分为两部分
%                 y_greater_than_0_5 = notInB(notInB(:, 2) > 0.5, :);  % y > 0.5 的点
%                 y_less_than_or_equal_0_5 = notInB(notInB(:, 2) <= 0.5, :);  % y ≤ 0.5 的点
% 
% %             % 显示结果
% %             disp('更新后的 y_greater_than_0_5:');
% %             disp(y_greater_than_0_5);
% %             disp('更新后的 y_less_than_or_equal_0_5:');
% %             disp(y_less_than_or_equal_0_5);
% 
%                 % 保存为 .dat 文件
%                 filename1 = 'points_y_greater_than_0_5.dat';
%                 filename2 = 'points_y_less_than_or_equal_0_5.dat';
% 
%                 % 写入文件
%                 dlmwrite(filename1, y_greater_than_0_5, 'delimiter', '\t', 'precision', 6);
%                 dlmwrite(filename2, y_less_than_or_equal_0_5, 'delimiter', '\t', 'precision', 6);
% 
%                 % 输出文件名
%                 disp(['文件保存完成：', filename1]);
%                 disp(['文件保存完成：', filename2]);
%             end
% %                         break



%             if step == 50 || step == 75 || step== 120 || step== 150
%                 figure;
%                 scatter(x(labels == 2), y(labels == 2), 'r', 'filled');hold on;
%                 scatter(x(labels == 1), y(labels == 1), 'b');hold on;
%                 xlabel('x');
%                 ylabel('y');axis equal
%                 xlim([- 0.1, 0.6]);
%                 ylim([0.5, 1.05]);
%             end
            
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

            % 设置缓冲区距离
            buffer_distance =- 0.25 * deltax;

            % 计算边界点的缓冲区
            buffer_x = [];
            buffer_y = [];
            for i = 1:length(k)
                p1 = [x(k(i)), y(k(i))];
                if i == length(k)
                    p2 = [x(k(1)), y(k(1))];  % 闭合边界
                else
                    p2 = [x(k(i+1)), y(k(i+1))];
                end
                direction = (p2 - p1) / norm(p2 - p1);
                normal = [-direction(2), direction(1)];  % 法向量
                buffer_x = [buffer_x; p1(1) + buffer_distance * normal(1)];
                buffer_y = [buffer_y; p1(2) + buffer_distance * normal(2)];
            end

            % 确保缓冲区多边形闭合
            buffer_x = [buffer_x; buffer_x(1)];
            buffer_y = [buffer_y; buffer_y(1)];
            
            [mask, ~] = inpolygon(x1, y1, buffer_x, buffer_y);
%             mask = in & ~on; % 保留在多边形内部的点，去除边缘上的点
%             hold on; axis equal;
%             plot(x_filtered(k), y_filtered(k));
%             mask = inpolygon(x1, y1, x_filtered(k), y_filtered(k));
%             hold on; scatter(x1(~mask), y1(~mask));

            % 筛选后的 data1 (空气粒子保留)
            data1_filtered = data1(~mask, :);

            % 合并 data 和 data1_filtered，保留完整的第四参数
            % 如果记事本文件需要追加多个物理量，在下一行修改，原本是[data; data1_filtered]容易报错
            data_combined = [data(:,1:4); data1_filtered]; 
            x_combined = data_combined(:, 2);
            y_combined = data_combined(:, 3);
            labels_combined = data_combined(:, 4);
            indices_combined = data_combined(:, 1);
            
            % 搜索样本为(x_filtered, y_filtered)，邻近半径为4e-2
            [idx, ~] = rangesearch([x_filtered, y_filtered], [x(k), y(k)], 5e-2);

            % 整合第一层筛选后的邻近点索引
            all_neighbors_idx = unique([idx{:}]);
            all_neighbors_mask = ismember(indices_combined, indices_filtered(all_neighbors_idx));
%             scatter(x_combined(all_neighbors_mask), y_combined(all_neighbors_mask), 'g');
%             hold on;

            % 第二次筛选，搜索样本为合并的(x_combined, y_combined)，邻近半径为3.015e-2
            [second_idx, second_dist] = rangesearch([x_combined, y_combined], [x_combined(all_neighbors_mask), y_combined(all_neighbors_mask)], 4e-2);
%             for i = 1:length(second_idx) disp(['Point ', num2str(i), ': ', num2str(second_idx{i})]); end

            % 筛选符合条件的点
            second_mask = cellfun(@(neighbors_idx) any(labels_combined(neighbors_idx) == 0), second_idx);
            valid_points_indices = find(second_mask);

            filtered_points_x = x_combined(all_neighbors_mask);
            filtered_points_y = y_combined(all_neighbors_mask);
            filtered_points_idx = indices_combined(all_neighbors_mask);

            filtered_points_x = filtered_points_x(valid_points_indices);
            filtered_points_y = filtered_points_y(valid_points_indices);
            filtered_points_idx = filtered_points_idx(valid_points_indices);
%             if step == 70
%             scatter(filtered_points_x, filtered_points_y, 'r', 'filled'); % 显示筛选后的点
%             hold on;
%             end
            % 初始化特征向量模长和单位特征向量数组
            feature_magnitudes = zeros(size(filtered_points_x));
            feature_vectors = zeros(length(filtered_points_x), 2);
            unit_feature_vectors = zeros(length(filtered_points_x), 2);
            
            % 计算特征向量和单位特征向量
            for i = 1:length(valid_points_indices)
                x_point = filtered_points_x(i);
                y_point = filtered_points_y(i);
                neighbors_idx = second_idx{valid_points_indices(i)}; % neighbors_idx是x_combined为背景的索引信息
                D = second_dist{valid_points_indices(i)}(:); % 确保 D 是列向量

                % 排除自身粒子
                valid_neighbors_idx = neighbors_idx(D ~= 0);
                valid_D = D(D ~= 0);
                dx = x_combined(valid_neighbors_idx) - x_point;
                dy = y_combined(valid_neighbors_idx) - y_point;

                % 计算 C_ij 和 vfrac_s
                C_ij = 1 * (labels_combined(valid_neighbors_idx) == 0); % 该值设置为正，保证法向量的方向由自由表面指向空气
                vfrac_s = calculate_vfrac_s(valid_D, deltax);

                % 计算 sum_j
                sum_j_x = sum(vfrac_s .* C_ij .* (dx ./ (valid_D .^ 2)));
                sum_j_y = sum(vfrac_s .* C_ij .* (dy ./ (valid_D .^ 2)));

                feature_vector = 2 / (pi * 3.015 * 3.015) * [sum_j_x, sum_j_y];
                feature_magnitude = norm(feature_vector);
                if feature_magnitude < 1e-2 * deltax % 设置截断误差
                    unit_feature_vector = [0, 0];
                else 
                    unit_feature_vector = feature_vector / feature_magnitude;
                end

                % 考虑附壁效应
                if abs(x_point) <= deltax * 1.5
                    unit_feature_vector = -[0, unit_feature_vector(2)/abs(unit_feature_vector(2))] .* sin(theta_eq) + [-1, 0] * cos(theta_eq);
%                     feature_vector = (feature_vector+unit_feature_vector*feature_magnitude)/2;
                    feature_vector = -unit_feature_vector*feature_magnitude;
                elseif abs(x_point - 0.5) <= deltax * 1.5
                    unit_feature_vector = -[0, unit_feature_vector(2)/abs(unit_feature_vector(2))] .* sin(theta_eq) + [1, 0] * cos(theta_eq);
%                     feature_vector = (feature_vector+unit_feature_vector*feature_magnitude)/2;
                    feature_vector = -unit_feature_vector*feature_magnitude;
                end

                % 存储特征向量模长和单位特征向量
                feature_vectors(i, :) = feature_vector;
                feature_magnitudes(i) = feature_magnitude;
                unit_feature_vectors(i, :) = unit_feature_vector;
            end
%             hold on;scatter(filtered_points_x, filtered_points_y, 'r', 'filled');
% quiver(filtered_points_x, filtered_points_y, unit_feature_vectors(:,1), unit_feature_vectors(:,2), 'AutoScale', 'on', 'AutoScaleFactor', mean_magnitude, 'Color', 'black');
% hold off;

            % 逐点计算曲率 kappa 和表面张力
            for i = 1:length(valid_points_indices)
                x_point = filtered_points_x(i);
                y_point = filtered_points_y(i);
                original_index = filtered_points_idx(i); % 获取原始数据的索引

                % 创建 meshgrid 进行坐标差值计算
                dx = filtered_points_x - x_point;
                dy = filtered_points_y - y_point;
                D = sqrt(dx.^2 + dy.^2);

                % 排除空粒子
                valid_mask = D ~= 0;
                dx = dx(valid_mask);
                dy = dy(valid_mask);
                D = D(valid_mask);
                % 这里不好改是因为ufv_x是以filtered_points_x作为背景信息，不是以x_combined作为背景信息，
                % 所以不能用上一条循环的逻辑来读取
                ufv_x = unit_feature_vectors(valid_mask, 1);
                ufv_y = unit_feature_vectors(valid_mask, 2);
                % 计算 vfrac_s，可以不用又去筛选每一点的近场特征向量，是一个很好的简化代码的思路
                vfrac_s = calculate_vfrac_s(D, deltax);
%                 vfrac_s = 0.5150;

                % 边界粒子表面校正
                weight = w_0(D, deltax);
                zeta = sum(weight) / sum((ufv_x .* ufv_x + ufv_y .* ufv_y) .* weight);
                
                % 计算曲率 kappa
                kappa = -2 * zeta / (pi * 3.015 * 3.015) * sum(vfrac_s .* ((ufv_x - unit_feature_vectors(i, 1)) .* dx + (ufv_y - unit_feature_vectors(i, 2)) .* dy) ./ (D .^ 2));
                
                % 计算表面张力
                surface_force = beta * kappa * feature_vectors(i,:); %* unit_feature_vectors(i, :) * feature_magnitudes(i);
                if norm(surface_force) < 1e-2 * deltax     % 应用截断误差
                    body_force(original_index, :) = [0 0];
                end
                body_force(original_index, :) = surface_force; % 使用原始数据的索引存储表面张力
            end
%             disp(body_force(364,:));
% if step ==50
%             quiver(x, y, body_force(:, 1), body_force(:, 2), 'b');hold on; % 显示表面张力矢量
%             break
% end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % 计算剩余体力 % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

         % 开始流体计算
         for i = 1:length(x)
            if labels(i) == 1 % 只对固体虚拟层粒子进行求解
                current_neighbors = neighbors{i};
                current_neighbors = current_neighbors(current_neighbors ~= indices(i)); % 排除自身点
                if isempty(current_neighbors)
                    continue;
                end
                
                dx = x(current_neighbors) - x(i);
                dy = y(current_neighbors) - y(i);
                D = sqrt(dx.^2 + dy.^2);
                U = sum((mean(D(:))^2) .* [dx./(D.^3), dy./(D.^3)], 1);
                delta_r_pst = C_pst * abs(max(velocity_magnitude)) * dt * U; 
                
                % 设置截断
                delta_r_pst(delta_r_pst >= 0.2 * deltax) = 0.2 * deltax;
                
                % 存入数组供最后更新位置使用
                pst_displacement(i, :) = delta_r_pst;

 
                % 计算 vfrac_s 
                vfrac_s = calculate_vfrac_s(D, deltax);
                
                % 计算压强梯度力
                pressure_j = pressure(current_neighbors);
                pressure_i = repmat(pressure(i),length(current_neighbors),1);
                body_force(i, 1) = body_force(i, 1)  - 2 / (pi * 3.015 * 3.015) * sum(vfrac_s .* (pressure_j - pressure_i) .* dx ./ (D .^ 2), 1);
                body_force(i, 2) = body_force(i, 2)  - 2 / (pi * 3.015 * 3.015) * sum(vfrac_s .* (pressure_j - pressure_i) .* dy ./ (D .^ 2), 1);
                
                % 计算粘性阻力
                velocity_j = velocity(current_neighbors, :);
                velocity_i = repmat(velocity(i, :),length(current_neighbors),1);
                body_force(i, :) = body_force(i, :) + 6 * miu / (pi * 3.015 * 3.015 * 3.015 * deltax) * sum([vfrac_s, vfrac_s] .* (velocity_j - velocity_i) ./ [D, D], 1);
                dv = velocity_j - velocity_i;
                
                % 更新粒子密度
                density(i) = density(i) - density(i) * dt * 2 / (pi * 3.015 * 3.015) * sum(vfrac_s .* (dv(:, 1) .* dx + dv(:, 2) .* dy)  ./ (D .^ 2), 1);
                
                % 更新粒子压强
                pressure(i) = reference_pressure * ((density(i)/density_f0) ^ gamma_f - 1);

                % 更新粒子加速度
                acceleration_now(i, :) = (body_force(i, :) + F_B) / density(i);
                
                % 更新粒子速度
                velocity(i, :) = velocity(i, :) + dt * (acceleration_old(i, :) + acceleration_now(i, :)) / 2;

                if any(pst_displacement(i, :))
                    % 预提取位移分量，让后面的式子稍微短一点（MATLAB也会自动优化）
                    du_x = pst_displacement(i, 1);
                    du_y = pst_displacement(i, 2);
                
                    velocity(i, 1) = velocity(i, 1) + ...
                        ... % --- [一次项] ---
                        (2 / (pi * 3.015^2)) * sum(vfrac_s .* (velocity_j(1) - velocity_i(1)) .* (dx .* du_x + dy .* du_y) ./ (D.^2), 1) + ...
                        ... % --- [二次项] ---
                        (3 * deltax^2 / (2 * pi * (3.015 * deltax)^3)) * sum(vfrac_s .* (velocity_j(1) - velocity_i(1)) .* ...
                        ((3 * dx.^2 - dy.^2) * du_x^2 + 8 * dx .* dy * du_x * du_y + (-dx.^2 + 3 * dy.^2) * du_y^2) ./ (D.^3), 1);
                    
                    velocity(i, 2) = velocity(i, 2) + ...
                        ... % --- [一次项] ---
                        (2 / (pi * 3.015^2)) * sum(vfrac_s .* (velocity_j(2) - velocity_i(2)) .* (dx .* du_x + dy .* du_y) ./ (D.^2), 1) + ...
                        ... % --- [二次项] ---
                        (3 * deltax^2 / (2 * pi * (3.015 * deltax)^3)) * sum(vfrac_s .* (velocity_j(2) - velocity_i(2)) .* ...
                        ((3 * dx.^2 - dy.^2) * du_x^2 + 8 * dx .* dy * du_x * du_y + (-dx.^2 + 3 * dy.^2) * du_y^2) ./ (D.^3), 1);
                    
                    pressure(i) = pressure(i) + ...
                        ... % --- [一次项] ---
                        (2 / (pi * 3.015^2)) * sum(vfrac_s .* (pressure_j - pressure_i) .* (dx .* du_x + dy .* du_y) ./ (D.^2), 1) + ...
                        ... % --- [二次项] ---
                        (3 * deltax^2 / (2 * pi * (3.015 * deltax)^3)) * sum(vfrac_s .* (pressure_j - pressure_i) .* ...
                        ((3 * dx.^2 - dy.^2) * du_x^2 + 8 * dx .* dy * du_x * du_y + (-dx.^2 + 3 * dy.^2) * du_y^2) ./ (D.^3), 1);
                end

%                 if x(i) < 0.5 * deltax
%                      if velocity(i, 1) < 0
%                          disp(velocity(i, :));
%                         velocity(i, :) = velocity(i, :) - 2 * velocity(i, 1) * [1, 0];
%                         disp(velocity(i, :));
%                      end
%                 elseif x(i) > 0.5 - 0.5 * deltax
%                     if velocity(i, 1) > 0
%                         velocity(i, :) = velocity(i, :) - 2 * velocity(i, 1) * [1, 0];
%                     end
%                 end
                
                % 更新粒子位移
                displacement(i, :) = displacement(i, :) + velocity(i, :) * dt + acceleration_old(i, :) * dt * dt / 2;
                if any(pst_displacement(i, :))
                    displacement(i, :) = displacement(i, :) + pst_displacement(i, :);
                end

%                 disp(acceleration_now(i, :) );
                % 更新粒子位置
                x(i) = x_0(i) + displacement(i, 1);
                y(i) = y_0(i) + displacement(i, 2);
                
                % 更新粒子速度
%                 velocity(i, :) = velocity(i, :) + dt * (acceleration_old(i, :) + acceleration_now(i, :)) / 2;
                if x(i) < 0.75 * deltax
                     if velocity(i, 1) < 0
                        displacement(i, :) = displacement(i, :) - velocity(i, :) * dt;
                        velocity(i, :) = velocity(i, :) - 2 * velocity(i, 1) * [1, 0];
                        displacement(i, :) = displacement(i, :) + velocity(i, :) * dt;
                        x(i) = x_0(i) + displacement(i, 1);
                        y(i) = y_0(i) + displacement(i, 2);
                     end
                elseif x(i) > 0.5 - 0.75 * deltax
                    if velocity(i, 1) > 0
                        displacement(i, :) = displacement(i, :) - velocity(i, :) * dt;
                        velocity(i, :) = velocity(i, :) - 2 * velocity(i, 1) * [1, 0];
                        displacement(i, :) = displacement(i, :) + velocity(i, :) * dt;
                        x(i) = x_0(i) + displacement(i, 1);
                        y(i) = y_0(i) + displacement(i, 2);
                    end
                end

%                 % 更新粒子速度
%                 velocity(i, :) = velocity(i, :) + dt * (acceleration_old(i, :) + acceleration_now(i, :)) / 2;
%                 if x(i) < 0.5 * deltax || x(i) > 0.5 - 0.5 * deltax
%                     disp(velocity(i, :));
%                      velocity(i, :) = velocity(i, :) - 2 * velocity(i, 1) .* [1, 0];
%                 end

                % 求速度绝对值
                velocity_magnitude(i, :) = sqrt(velocity(i, 1)^2 + velocity(i, 2)^2);
                
            end
            
                % 将当前时间步的加速度存储为上一个时间步的加速度 
                acceleration_old = acceleration_now;
                
         end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % 储存相关信息 % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        fidnew=fopen('C:\Users\ASUS\Desktop\soil_board_now.dat','w+'); % 打开文件用于写入数据
        
        for i = 1:length(indices)
            fprintf(fidnew, '%d\t', indices(i));          % 写入每个粒子的索引编号
            fprintf(fidnew, '%15.5e\t', x(i));            % 写入每个粒子的位置坐标 x
            fprintf(fidnew, '%15.5e\t', y(i));            % 写入每个粒子的位置坐标 y
            fprintf(fidnew, '%d\t', labels(i));           % 写入每个粒子的边界信息
            fprintf(fidnew, '%15.5e\t', velocity_magnitude(i));      % 写入每个粒子的密度
            fprintf(fidnew, '\n');                        % 换行
        end

        fclose(fidnew); % 关闭文件
        
        
% %         figure;
%         scatter(x(labels == 2), y(labels == 2), 'r', 'filled'); hold on; % 绘制固体边界
%         scatter(x(labels == 1), y(labels == 1), 'b'); hold on; % 绘制流体粒子
%         xlabel('x');
%         ylabel('y');
%         axis equal;
%         xlim([min(filtered_points_x) - 0.1, max(filtered_points_x) + 0.1]);
%         ylim([min(filtered_points_y) - 0.2, max(filtered_points_y) + 0.2]);
% %         title(['Time Step: ', num2str(step)]);
%         filename = sprintf('C:\\Users\\ASUS\\Desktop\\video\\frame_%04d.png', step); % 保存为PNG文件
%         saveas(gcf, filename); % 保存图像
%         close(gcf); % 关闭当前图形窗口

end

% 展示结果
final_data = load('C:\Users\ASUS\Desktop\soil_board_now.dat');
x_new = final_data(:, 2); % 第二列是x坐标
y_new = final_data(:, 3); % 第三列是y坐标
labels = final_data(:, 4); % 第四列是标签
density_new = final_data(:, 5); % 第五列是密度
scatter(x_new(labels == 2), y_new(labels == 2), 'r', 'filled');hold on;
scatter(x_new(labels == 1), y_new(labels == 1), 'b');hold on;
scatter(x_new(labels == 1), y_new(labels == 1), 36, density_new(labels == 1), 'filled'); hold on;
% scatter(xA_notInB, yA_notInB);
c = colorbar; % 添加颜色条以显示密度值
c.Label.String = '速度 (cm/s)';
xlabel('x');
ylabel('y');
axis equal;
xlim([min(filtered_points_x) - 0.1, max(filtered_points_x) + 0.1]);
ylim([min(filtered_points_y) - 0.2, max(filtered_points_y) + 0.2]);