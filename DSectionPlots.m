% D section fitting

Input = load('DSectionRawData_a_b.mat');
One = polyfit(x1,y1,5);

OnePointFive = polyfit(x1_5,y1_5,5);

Two = polyfit(x2,y2,5);

Three = polyfit(x3,y3,5);

% Plots to chek validity
x_coord = 0:0.01:10;
y1_coord = polyval(One,x_coord)
y1_5_coord = polyval(OnePointFive,x_coord)
y2_coord = polyval(Two,x_coord)
y3_coord = polyval(Three,x_coord)

hold on
plot(x_coord,y1_coord)
plot(x_coord,y1_5_coord)
plot(x_coord,y2_coord)
plot(x_coord,y3_coord)
hold off

%%

Input = load('DSectionRawData_b_a.mat');
bOne = polyfit(x_b_a_1,y_b_a_1,5);

bOnePointFive = polyfit(x_b_a_1_5,y_b_a_1_5,5);

bFive = polyfit(x_b_a_5,y_b_a_5,5);

% Plots to chek validity
bx_coord = 0:0.01:10;
by1_coord = polyval(bOne,bx_coord)
by1_5_coord = polyval(bOnePointFive,bx_coord)
by5_coord = polyval(bFive,bx_coord)

hold on
plot(bx_coord,by1_coord)
plot(bx_coord,by1_5_coord)
plot(bx_coord,by5_coord)
hold off


save("DSectionFnb.mat","bOne","bOnePointFive","bFive");