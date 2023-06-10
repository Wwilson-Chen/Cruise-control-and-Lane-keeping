function delta = str_ang(er)

 e1 = er(1);
 e1_d = er(2);
 e1_i = er(10);
 e2 = er(3);
 e2_d = er(4);
 e2_i = er(11);
 Kp1 = 5000; % 1000
 Kd1 = 500; % 100
 Ki1 = 10; % 0
 Kp2 = 500; % 0
 Kd2 = 50; % 0
 Ki2 = 0;

 delta = - Kp1*e1 - Kd1*e1_d - 50*Ki1*e1_i - Kp2*e2 - Kd2*e2_d - Ki2*e2_i;
end