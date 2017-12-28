% The Generalized Complementary Filter (GCF) proposed by Jin Wu
% Author: Jin Wu
% e-mail: jin_wu_uestc@hotmail.com


function [q, qy] = GCF(omega, dt, gain, chi, q_, Db, Dr, valid, weight)

     wx = omega(1);     wy = omega(2);     wz = omega(3);

     Omg = [   0, - wx, - wy, - wz;
              wx,    0,   wz, - wy;
              wy, - wz,    0,   wx;
              wz,   wy, - wx,    0];

     total = 0;

     for i = 1 : length(valid)
         index = valid(i);
    
         if(index ~= 0)
             M1 = M1_matrix(Db(index, :));
             M2 = M2_matrix(Db(index, :));
             M3 = M3_matrix(Db(index, :));
             
             total = total + ...
                 weight(index) * (eye(4) - Dr(index, 1) * M1 - Dr(index, 2) * M2 - Dr(index, 3) * M3);
         end    
     end

     qy = (eye(4) - 2 * chi * total) * q_;
     qy = qy ./ norm(qy);

     q = q_ + (eye(4) - gain) ./ 2 .* dt * Omg * q_;
     q = (1 - gain) * q + gain * qy;
     q = q ./ norm(q);
end