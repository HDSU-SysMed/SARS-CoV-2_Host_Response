function dx = sim_SARS_CoV2_model_ODEs(t, x)

global k

dx = zeros(4, 1);
dx( 1) = -k(2)*x(1)+(k(1)*x(1)/(1+x(4)/k(9))); % V
dx( 2) = k(7)*x(1)/(k(11)+x(1))-(k(8)*x(2)); % P
dx( 3) = k(3)* (x(1)/(x(1)+(k(12)+k(10))))/(1+x(2)/k(4))-(k(5)*x(3)); % M
dx( 4) = k(6)*x(3)/(k(11)+x(3))-(k(6)*x(4)); % A
