clear all
close all
clc

addpath('../IMUDataSets');
addpath('../Madgwick/quaternion_library');

data = load('flat_motion2.txt');
conv = diag([-1 1 -1]);
euler_true = data(:, 1 : 3) * conv;
acc = data(:, 4 : 6) * conv;
gyro = data(:, 7 : 9) * conv; 

num = 5;
Dr = zeros(3, num);
Db = zeros(3, num);
weight = abs(randn(num, 1));
weight = weight ./ sum(weight);

len = length(data(:, 1));
dt = 1 / 500;
gain = 0.01;
chi = 0.1;
q = [1; 0; 0; 0];
qy = q;
time = dt * (1 : len);

quaternion = zeros(len, 4);
quaternion_meas = zeros(len, 4);

for i = 1 : len
    C_true = euler2rotMat( euler_true(i, 1) * pi / 180, ...
                           euler_true(i, 2) * pi / 180, ...
                           euler_true(i, 3) * pi / 180)';
    for j = 1 : num
        Dr(:, j) = randn(3, 1);
        Dr(:, j) = Dr(:, j) ./ norm(Dr(:, j));
        
        Db(:, j) = C_true * Dr(:, j) + randn(3, 1) * 0.1;
        Db(:, j) = Db(:, j) ./ norm(Db(:, j));
    end
    
    
    [q, qy] = GCF(gyro(i, :), dt, gain, chi, q, qy, Db, Dr, num, weight);
    
    quaternion(i, :) = q';
    quaternion_meas(i, :) = qy';
end

euler = quatern2euler(quaternConj(quaternion)) * 180 / pi;
euler_meas = quatern2euler(quaternConj(quaternion_meas)) * 180 / pi;


figure(1);
subplot(3, 1, 1);
plot(time, euler(:, 1), 'LineWidth', 1); hold on;
plot(time, euler_true(:, 1), 'LineWidth', 1); hold off;
legend('GCF', 'True');
xlabel('Time (s)');
ylabel('Roll (deg)');

subplot(3, 1, 2);
plot(time, euler(:, 2), 'LineWidth', 1); hold on;
plot(time, euler_true(:, 2), 'LineWidth', 1); hold off;
legend('GCF', 'True');
xlabel('Time (s)');
ylabel('Pitch (deg)');

subplot(3, 1, 3);
plot(time, euler(:, 3), 'LineWidth', 1); hold on;
plot(time, euler_true(:, 3), 'LineWidth', 1); hold off;
legend('GCF', 'True');
xlabel('Time (s)');
ylabel('Yaw (deg)');



figure(2);
subplot(3, 1, 1);
plot(time, euler_meas(:, 1), 'LineWidth', 1); hold on;
plot(time, euler_true(:, 1), 'LineWidth', 1); hold off;
legend('GCF Meas', 'True');
xlabel('Time (s)');
ylabel('Roll (deg)');

subplot(3, 1, 2);
plot(time, euler_meas(:, 2), 'LineWidth', 1); hold on;
plot(time, euler_true(:, 2), 'LineWidth', 1); hold off;
legend('GCF Meas', 'True');
xlabel('Time (s)');
ylabel('Pitch (deg)');

subplot(3, 1, 3);
plot(time, euler_meas(:, 3), 'LineWidth', 1); hold on;
plot(time, euler_true(:, 3), 'LineWidth', 1); hold off;
legend('GCF Meas', 'True');
xlabel('Time (s)');
ylabel('Yaw (deg)');
