%%   CS 543 Homework 3 
%%   Chongxin Luo
%%   March. 19, 2017
%%   Running file for Question 1: Single-View Metrology

%% Part (a): Find three vanish point in the image 
clc; clear all; close all;
im = im2double(imread('kyoto_street.JPG'));
disp('Indicate lines for right vanish point');
p1 = getVanishingPoint_shell(im);
disp(p1);
disp('Indicate lines for left vanish point');
p2 = getVanishingPoint_shell(im);
disp(p2);
disp('Indiate lines for top vanish point');
p3 = getVanishingPoint_shell(im);
disp(p3);
imagesc(im);
figure(1), hold off, imagesc(im)
hold on 
plot(p1(1)/p1(3), p1(2)/p1(3), '*r');
plot(p2(1)/p2(3), p2(2)/p2(3), '*r');
plot([p1(1)/p1(3) p2(1)/p2(3)],[p1(2)/p1(3) p2(2)/p2(3)],'r');
axis image
horizon = cross(p1,p2)
horizon(1) = horizon(1)/horizon(3);
horizon(2) = horizon(2)/horizon(3);
horizon(3) = horizon(3)/horizon(3);
norm = horizon(1)*horizon(1) + horizon(2)*horizon(2);
horizon(1) = horizon(1)/norm; horizon(2) = horizon(2)/norm;
horizon
%% part(a) Results
% Right vanish point:
%    1.0e+09 *
% 
%    -1.3618
%    -0.2830
%    -0.0002
% Left vanish point:
%    1.0e+08 *
% 
%     0.4947
%    -1.1220
%    -0.0007
% Top vanish point:
%    1.0e+08 *
% 
%    -0.4996
%     6.1746
%    -0.0004
% horizon =
% 
%    1.0e+17 *
% 
%    -0.0000
%    -0.0011
%     1.6679
% horizon =
% 
%    1.0e+03 *
% 
%    -0.0242
%    -1.5616
%     0.0010
%%
clc;
line = [-24.2, -1561.6, 1];
n = norm(line);
normline = line/n;
disp(normline);
disp(normline(1)*normline(1) + normline(2)*normline(2));


%% Part (b): Calculate focal length (f) and optimal center (u0,v0)
% initialize vanish point and matrix k 
syms f u v
vp1 = [-49960000; 617460000; -40000];       % top vanish
vp2 = [-1361800000; -283000000;-200000];    % right vanish
vp3 = [49470000; -112200000; -70000];       % left vanish

vp1 = vp1/vp1(3); vp2 = vp2/vp2(3); vp3 = vp3/vp3(3);

k = [f 0 u; 0 f v; 0 0 1];
kinv = [1/f 0 -u/f; 0 1/f -v/f; 0 0 1];

% setup system of equation and solve system of equation 
eq1 = transpose(vp1) * transpose(kinv) * kinv * vp2 == 0;
eq2 = transpose(vp1) * transpose(kinv) * kinv * vp3 == 0;
eq3 = transpose(vp2) * transpose(kinv) * kinv * vp3 == 0;

result = solve([eq1, eq2, eq3], [f;u;v]);
u0 = double(result.u(1));
v0 = double(result.v(1));
f = double(result.f(2));
disp(double(result.u));
disp(double(result.v));
disp(double(result.f));

%% Part (c): Compute the camera's rotation matrix R

kinv = [1/f 0 -u0/f; 0 1/f -v0/f; 0 0 1];
r1 = kinv * vp2; 
r3 = kinv * vp3;
r2 = cross(r3,r1);
R = [r1,r2,r3];
R = normc(R);
disp(R);
disp(R*R');

%% Part(d): Find horizon line, estimate heights 
clc; clear all; close all;
im = im2double(imread('CIMG6476.JPG'));
disp('Indicate lines for right vanish point');
p1 = getVanishingPoint_shell(im);
disp(p1);
disp('Indicate lines for left vanish point');
p2 = getVanishingPoint_shell(im);
disp(p2);
disp('Indiate lines for top vanish point');

imagesc(im);
figure(1), hold off, imagesc(im)
hold on 
plot(p1(1)/p1(3), p1(2)/p1(3), '*r');
plot(p2(1)/p2(3), p2(2)/p2(3), '*r');
plot([p1(1)/p1(3) p2(1)/p2(3)],[p1(2)/p1(3) p2(2)/p2(3)],'r', 'Linewidth', 1);
horizon = cross(p1,p2);
disp(p1); disp(p2);
disp(horizon);
% horizon(1) = horizon(1)/horizon(3);
% horizon(2) = horizon(2)/horizon(3);
% horizon(3) = horizon(3)/horizon(3);
% norm = horizon(1)*horizon(1) + horizon(2)*horizon(2);
% horizon(1) = horizon(1)/norm; horizon(2) = horizon(2)/norm;
% disp(horizon);
% horizon
%    1.0e+17 *
% 
%     0.0001
%    -0.0075
%     9.1653
% right vanish
%    1.0e+09 *
% 
%     1.4834
%     0.3735
%     0.0003
% left vanish
%    1.0e+08 *
% 
%    -4.9666
%     4.9281
%     0.0041
%%
clc; close all; clear all;
im = im2double(imread('CIMG6476.JPG'));
imagesc(im);
hold on
rightvp = [1.4834 0.3735 0.0003]*1.0e+09;
leftvp = [-4.9666 4.9281 0.0041]*1.0e+08;
plot([rightvp(1)/rightvp(3) leftvp(1)/leftvp(3)],[rightvp(2)/rightvp(3) leftvp(2)/leftvp(3)],'r', 'Linewidth', 1);
% six points for three heights : (bottom, top)
% sign: (989, 1435) (986, 1353)
% tractor: (1500, 1510) (1500, 1391)
% building: (915, 1317) (918, 693)

horizon = [0.0001; -0.0075; 9.1653]*1.0e+17;
sign1 = [1037; 2191; 1]; sign2 = [1036; 2065; 1];
tractor1 = [1761; 2313; 1]; tractor2 = [1762; 2123; 1];
building1 = [933; 1990; 1]; building2 = [930; 923; 1];

hold on 
plot([sign1(1) sign2(1)], [sign1(2) sign2(2)], 'g', 'Linewidth', 1);
plot([tractor1(1) tractor2(1)], [tractor1(2) tractor2(2)], 'g', 'Linewidth', 1);
plot([building1(1) building2(1)], [building1(2) building2(2)], 'g', 'Linewidth', 1);

% find building height
bottomline = cross(sign1, building1);
intersect = cross(bottomline, horizon);
intersect = intersect/intersect(3);
plot([intersect(1) sign1(1)], [intersect(2) sign1(2)], 'b', 'Linewidth', 1);
plot(intersect(1), intersect(2), '*y');
topline = cross(sign2, intersect);
plot([intersect(1) sign2(1)], [intersect(2) sign2(2)], 'b', 'Linewidth', 2);
buildingline = cross(building1, building2);
buildingcross = cross(topline, buildingline);
buildingcross = buildingcross/buildingcross(3);
plot(buildingcross(1), buildingcross(2), '*y');
buildingheight = (abs(building1(2)-building2(2))*1.65) / abs(buildingcross(2)-building1(2));
disp(buildingheight);

% find tractor height
bottomline = cross(sign1, tractor1);
intersect = cross(bottomline, horizon);
intersect = intersect/intersect(3);
plot([intersect(1) tractor1(1)], [intersect(2) tractor1(2)], 'b', 'Linewidth', 1);
plot(intersect(1), intersect(2), '*y');
topline = cross(sign2, intersect);
tractorline = cross(tractor1, tractor2);
tractorcross = cross(topline, tractorline);
tractorcross = tractorcross/tractorcross(3);
plot([intersect(1) tractorcross(1)], [intersect(2) tractorcross(2)], 'b', 'Linewidth', 2);
plot(tractorcross(1), tractorcross(2), '*y');
tractorheight = (abs(tractor1(2)-tractor2(2))*1.65) / abs(tractorcross(2)-tractor1(2));
disp(tractorheight);


