clc
clear vars

%% constants
red = 52.3813;
green = 23.6906;
blue = 57.4895;

l78 = 38.0938;
l67 = 44.8596;

rod_dia = 3/32*25.4; % inches*mm/1inch
Area = rod_dia^2*pi; % mm^
E = 125e9;

F_load = 9.81*1; % N

theta_1_12 = asind(green/blue);
theta_6_7 = asind(green/l67);

%% external forces
F2_y = F_load*(5*red + 1*l78)/red;

F1_y = F_load - F2_y;

%% internal stresses - solving
T1_12 = -F1_y/sind(theta_1_12);
T1_2 = -T1_12*cosd(theta_1_12);

T6_7 = -F_load/sind(theta_6_7);
T7_8 = -T6_7*cosd(theta_6_7);

T6_8 = -T6_7*sind(theta_6_7);
T5_6 = T6_7*cosd(theta_6_7);

T5_8 = -T6_8/sind(theta_1_12);
T8_9 = T7_8 - T5_8*cosd(theta_1_12);

T5_9 = 0;
T9_10 = T8_9;

T5_10 = -T5_8;
T4_5 = -T5_10*cosd(theta_1_12) + T5_8*cosd(theta_1_12) + T5_6;

T3_4 = T4_5;
T4_10 = 0;

T3_10 = -T5_10;
T10_11 = T9_10 - T3_10*cosd(theta_1_12) + T5_10*cosd(theta_1_12);

T11_12 = T10_11;
T3_11 = 0;

T3_12 = -T3_10;
T2_3 = T3_4 + T3_10*cosd(theta_1_12) - T3_12*cosd(theta_1_12);

T2_12 = -F2_y;
T1_2 = T2_3;

T12_13 = -T1_12*cosd(theta_1_12) + T11_12 + T3_12*cosd(theta_1_12);
T1_13 = 0;
T12_13 = 0;

link = ["T1_2"; "T1_12"; "T1_13";"T2_3"; "T2_12"; "T3_4"; "T3_10"; "T3_11"; "T3_12";
	"T4_5"; "T4_10";"T5_6"; "T5_8"; "T5_9";"T5_10"; "T6_7"; "T6_8"; "T7_8";
	"T8_9"; "T9_10";"T10_11"; "T11_12"; "T12_13"];

lengths = [red; blue; green; red; green; red; blue; green; blue;
	red; green; red; blue; green; blue; l67; green; l78;
	red; red; red; red; red];

stresses = [T1_2; T1_12; T1_13; T2_3; T2_12; T3_4; T3_10; T3_11; T3_12;
	T4_5; T4_10; T5_6; T5_8; T5_9; T5_10; T6_7; T6_8; T7_8;
	T8_9; T9_10; T10_11; T11_12; T12_13];

deflection = stresses.*lengths./(Area*E);

buckling = stresses.*(stresses<0)

table = [link, lengths, stresses, deflection]


%% internal stresses - system of equations
% x = ["T1_2"; "T1_12"; "T1_13";
% 	"T2_3"; "T2_12"; "T3_12"; "T3_4"; "T3_10"; "T3_11";
% 	"T4_5"; "T4_10";
% 	"T5_6"; "T5_8"; "T5_9";
% 	"T5_10"; "T6_7"; "T6_8"; "T7_8"; "T8_9"; "T9_10";
% 	"T10_11"; "T11_12"; "T12_13"]';
% 
% A1 = 1*(x=="T12_13");% T12_13 = 0;
% A2 = 1*(x=="T1_13");% T1_13 = 0;
% A3 = 1*(x=="T5_9");% T5_9 = 0;
% A4 = 1*(x=="T4_10");% T4_10 = 0;
% A5 = 1*(x=="T3_11");% T3_11 = 0;
% 
% 
% A6 = -1*(x=="T7_8") -1*(x=="T6_7")/cosd(theta_1_12); % -T7_8 -T6_7*sind(theta_1_12);
% A7 = -1*(x=="T6_8") -sind(theta_1_12)*(x=="T6_7"); %-T5_6 + T6_7*cosd(theta_1_12);
% A8 = -1*(x=="T_6") + cosd(theta_1_12)*(x=="T6_7");
% A9 = -1*(x=="T5_6") -1*(x=="T6_8")/sind(theta_1_12);% -T5_8 -T6_8/sind(theta_1_12);
% A10 = -1*(x=="T8_9") + (x=="T7_8") -1*(x=="T5_8")*cosd(theta_1_12);% -T8_9 + T7_8-T5_8*cosd(theta_1_12);
% 
% A11 = -1*(x=="T9_10") + (x=="T8_9");% -T9_10 + T8_9;
% A12 = -1*(x=="T5_10") -1*(x=="T5_8");% -T5_10 -T5_8;
% % -T4_5 + -T5_10*cosd(theta_1_12) + T5_8*cosd(theta_1_12) + T5_6;
% A13 = -1*(x=="T4_5") + -1*(x=="T5_10")*cosd(theta_1_12) + (x=="T5_8")*cosd(theta_1_12) + (x=="T5_6");
% A14 = -1*(x=="T3_4") + (x=="T4_5");% -T3_4 + T4_5;
% A15 = -1*(x=="T3_10") -1* (x=="T5_10");% -T3_10 -T5_10;
% 
% % -T10_11 + T9_10 - T3_10*cosd(theta_1_12) + T5_10*cosd(theta_1_12);
% A16 = -1*(x=="T10_11") + (x=="T9_10") -1*(x=="T3_10")*cosd(theta_1_12) + (x=="T5_10")*cosd(theta_1_12);
% A17 = -1*(x=="T11_12") + (x=="T10_11");% -T11_12 + T10_11;
% A18 = -1*(x=="T3_12") -1*(x=="T3_10");% -T3_12 -T3_10;
% % -T2_3 + T3_4+T3_10*cosd(theta_1_12) - T3_12*cosd(theta_1_12);
% A19 = -1*(x=="T2_3") + (x=="T3_4")+(x=="T3_10")*cosd(theta_1_12) - (x=="T3_12")*cosd(theta_1_12);
% A20 = -1*(x=="T1_2") -1*(x=="T1_12")/cosd(theta_1_12);% -T1_2 -T1_12/cosd(theta_1_12);
% 
% A21 = -1*(x=="T1_2") + (x=="T2_3");% -T1_2 + T2_3;
% % -T1_12 + (-T2_12 - T3_12*sind(theta_1_12))/sind(theta_1_12);
% A22 = -1*(x=="T1_12") + (-1*(x=="T2_12")-1*(x=="T3_12")*sind(theta_1_12))/sind(theta_1_12);
% % -T12_13 -T1_12*cosd(theta_1_12) + T11_12 + T3_12*cosd(theta_1_12);
% A23 = -1*(x=="T12_13") -1*(x=="T1_12")*cosd(theta_1_12) + (x=="T11_12") + (x=="T3_12")*cosd(theta_1_12);
% A24 = -1*(x=="T6_7");
% A25 = -1*(x=="T1_12");% -T1_12 -F1_y/sind(theta_1_12);
% A26 = -1*(x=="T2_12");% -T2_12 -F2_y;
% 
% A = [A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;
% 	A11;A12;A13;A14;A15;A16;A17;A18;A19;A20;
% 	A21;A22;A23;A24;A25;A26];
% 
% b = [zeros(23,1);
% -F_load/sind(theta_1_12);
% -1*F1_y/sind(theta_1_12);
% -F2_y;
% ];
% 
% x = -inv(A)*b
% 
