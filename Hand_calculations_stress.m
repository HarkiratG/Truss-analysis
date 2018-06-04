clc
clear vars

%% constants
green = 23.6906/1000;
red = 52.3813/1000;
blue = 57.4895/1000;

l78 = 38.0938/1000;
l67 = 44.8596/1000;

rod_dia = 3/32*25.4/1000; % inches*mm/1inch*1m/1000m
Area = rod_dia^2*pi; % m^2
E = 100e9; 

F_load = 9.81*1; % N

theta_1_12 = asind(green/blue);
theta_6_7 = asind(green/l67);

%% p0 system


% external forces
F2_y_p0 = F_load*(5*red + 1*l78)/red;

F1_y_p0 = F_load - F2_y_p0;

% internal stresses - solving
T1_12_p0 = -F1_y_p0/sind(theta_1_12);
T1_2_p0 = -T1_12_p0*cosd(theta_1_12);

T6_7_p0 = -F_load/sind(theta_6_7);
T7_8_p0 = -T6_7_p0*cosd(theta_6_7);

T6_8_p0 = -T6_7_p0*sind(theta_6_7);
T5_6_p0 = T6_7_p0*cosd(theta_6_7);

T5_8_p0 = -T6_8_p0/sind(theta_1_12);
T8_9_p0 = T7_8_p0 - T5_8_p0*cosd(theta_1_12);

T5_9_p0 = 0;
T9_10_p0 = T8_9_p0;

T5_10_p0 = -T5_8_p0;
T4_5_p0 = -T5_10_p0*cosd(theta_1_12) + T5_8_p0*cosd(theta_1_12) + T5_6_p0;

T3_4_p0 = T4_5_p0;
T4_10_p0 = 0;

T3_10_p0 = -T5_10_p0;
T10_11_p0 = T9_10_p0 - T3_10_p0*cosd(theta_1_12) + T5_10_p0*cosd(theta_1_12);

T11_12_p0 = T10_11_p0;
T3_11_p0 = 0;

T3_12_p0 = -T3_10_p0;
T2_3_p0 = T3_4_p0 + T3_10_p0*cosd(theta_1_12) - T3_12_p0*cosd(theta_1_12);

T2_12_p0 = -F2_y_p0;
T1_2_p0 = T2_3_p0;

T12_13_p0 = -T1_12_p0*cosd(theta_1_12) + T11_12_p0 + T3_12_p0*cosd(theta_1_12);
T1_13_p0 = 0;
T12_13_p0 = 0;

link = ["T1_2"; "T1_12"; "T1_13";"T2_3"; "T2_12"; "T3_4"; "T3_10"; "T3_11"; "T3_12";
	"T4_5"; "T4_10";"T5_6"; "T5_8"; "T5_9";"T5_10"; "T6_7"; "T6_8"; "T7_8";
	"T8_9"; "T9_10";"T10_11"; "T11_12"; "T12_13"];

lengths = [red; blue; green; red; green; red; blue; green; blue;
	red; green; red; blue; green; blue; l67; green; l78;
	red; red; red; red; red];

inertia_diagonal = 1;
inertia_vertical = 2;
inertia_horizontal = 3;
inertia_67 = 4;
inertia_78 = 5;

Inertia = [(lengths==green).*inertia_vertical + ...
	(lengths==red).*inertia_horizontal + ...
	(lengths==blue).*inertia_diagonal + ...
	(lengths==l67).*inertia_67 + ...
	(lengths==l78).*inertia_78];

tension_p0 = [T1_2_p0; T1_12_p0; T1_13_p0; T2_3_p0; T2_12_p0; T3_4_p0; T3_10_p0; T3_11_p0; T3_12_p0;
	T4_5_p0; T4_10_p0; T5_6_p0; T5_8_p0; T5_9_p0; T5_10_p0; T6_7_p0; T6_8_p0; T7_8_p0;
	T8_9_p0; T9_10_p0; T10_11_p0; T11_12_p0; T12_13_p0];

stresses = tension_p0./Area;

deflection0 = tension_p0.*lengths./(Area*E);

T_v1 = [(link=="T10_11")*(-red/blue) + (link=="T3_11")*(-green/blue) + ...
	(link=="T3_10")*1 + (link=="T3_4")*(-red/blue) + (link=="T3_10")*(-green/blue)];

W1 = T_v1.*deflection0;

work_sum_p0 = sum(W1);

buckling = (tension_p0<0)*pi^2*E.*Inertia./(lengths.^2);

p0_system = ["Links" "Lengths" "Member Forces" "Deflections" "Virtual Internal Force" "Virtual internal work";
	link, lengths, tension_p0, deflection0, T_v1, W1]

%% p1 system
deflection1 = T_v1.*lengths./(Area*E);
W2 = T_v1.*deflection1;
p1_system = ["Links" "Lengths" "Member Forces" "Deflections" "Virtual Internal Force" "Virtual internal work";
	link, lengths, T_v1+"P", deflection1+"P", T_v1, W1+"P"]

work_sum_p1 = sum(W2);

P = -work_sum_p0/work_sum_p1

%% p0 + p1 system
% external forces
F_load = 1; %1 N
F2_y_p01 = F_load*(5*red + 1*l78)/red;

F1_y_p01 = F_load - F2_y_p01;

% internal stresses - solving
T1_12_p01 = -F1_y_p01/sind(theta_1_12);
T1_2_p01 = -T1_12_p01*cosd(theta_1_12);

T6_7_p01 = -F_load/sind(theta_6_7);
T7_8_p01 = -T6_7_p01*cosd(theta_6_7);

T6_8_p01 = -T6_7_p01*sind(theta_6_7);
T5_6_p01 = T6_7_p01*cosd(theta_6_7);

T5_8_p01 = -T6_8_p01/sind(theta_1_12);
T8_9_p01 = T7_8_p01 - T5_8_p01*cosd(theta_1_12);

T5_9_p01 = 0;
T9_10_p01 = T8_9_p01;

T5_10_p01 = -T5_8_p01;
T4_5_p01 = -T5_10_p01*cosd(theta_1_12) + T5_8_p01*cosd(theta_1_12) + T5_6_p01;

T3_4_p01 = T4_5_p01;
T4_10_p01 = 0;

T3_10_p01 = -T5_10_p01;
T10_11_p01 = T9_10_p01 - T3_10_p01*cosd(theta_1_12) + T5_10_p01*cosd(theta_1_12);

T11_12_p01 = T10_11_p01;
T3_11_p01 = 0;

T3_12_p01 = -T3_10_p01;
T2_3_p01 = T3_4_p01 + T3_10_p01*cosd(theta_1_12) - T3_12_p01*cosd(theta_1_12);

T2_12_p01 = -F2_y_p01;
T1_2_p01 = T2_3_p01;

T12_13_p01 = -T1_12_p01*cosd(theta_1_12) + T11_12_p01 + T3_12_p01*cosd(theta_1_12);
T1_13_p01 = 0;
T12_13_p01 = 0;

tension_p0_1 = tension_p0 + T_v1.*P;

deflection_01 = tension_p0_1.*lengths./(Area*E);

T_v0_1 = [T1_2_p01; T1_12_p01; T1_13_p01; T2_3_p01; T2_12_p01; T3_4_p01; T3_10_p01; T3_11_p01; T3_12_p01;
	T4_5_p01; T4_10_p01; T5_6_p01; T5_8_p01; T5_9_p01; T5_10_p01; T6_7_p01; T6_8_p01; T7_8_p01;
	T8_9_p01; T9_10_p01; T10_11_p01; T11_12_p01; T12_13_p01];

W01 = T_v0_1.*deflection_01;
p0_1_system = ["Links" "Lengths" "Member Forces" "Deflections" "Virtual Internal Force" "Virtual internal work";
	link, lengths, tension_p0_1, deflection_01, T_v0_1, W01]

final_deflection = sum(W01)
%%

