% The Generalized Complementary Filter (GCF) proposed by Jin Wu
% Author: Jin Wu
% e-mail: jin_wu_uestc@hotmail.com


function [q, qy] = GCF(omega, dt, gain, chi, q_, qy_, Db, Dr, num, weight)

     wx = omega(1);     wy = omega(2);     wz = omega(3);

     total = 0;

     for i = 1 : num
         M1 = M1_matrix(Db(:, i));
         M2 = M2_matrix(Db(:, i));
         M3 = M3_matrix(Db(:, i));
             
         total = total + ...
                 weight(i) * (eye(4) - Dr(1, i) * M1 - Dr(2, i) * M2 - Dr(3, i) * M3);   
     end

     qy = (eye(4) - 2 * chi * total) * qy_;
     qy = qy ./ norm(qy);

     q = q_ + 0.5 * dt * quaternProd(q_', [0, wx, wy, wz])';
     q = q ./ norm(q);
     q = (1 - gain) * q + gain * qy;
     q = q ./ norm(q);
end