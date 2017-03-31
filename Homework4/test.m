%% Convert int to string test 
clc; close all; clear all;
for ii = 0:10;
    disp(num2str(ii+1));
end 

%% Matlab Map test
clc; close all; clear all;
table = containers.Map;

table('foo') = 1;
table
% for ii = 0:10
%     key = num2str(ii) + num2str(ii);
%     disp(key);
%     table(key) = 1;
% end
% table

%% for loop increment by 10 test
for ii = 1:10:100
    disp(ii);
end 

%% matrix test
m = [2,2;
     4,2;
     6,2;
     8,2];
n = [2;2;2;2];
 disp(m);
m = m ./ n;
 disp(m)