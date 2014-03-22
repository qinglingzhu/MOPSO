%输入pareto端面的值，画出该端面。2目标或3目标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Frontshow(Epa,Marker)
%--------------------------------------------------------------------------
if size(Epa,2)==2 %如果Epa的列数为2  即目标函数的个数有两个
    %%plot pareto fronts in 2-D
    f1=Epa(:,1);%第一个目标函数的值赋给f1
    f2=Epa(:,2);%第二个目标函数的值赋给f2
    plot(f1,f2,Marker,'MarkerSize',4); %以红色的星星来显示目标函数值的分布情况
    grid on; %显示网格
    xlabel('Function 1');%x坐标轴的标记为Function 1
    ylabel('Function 2');%y坐标轴的标记为Function 2
    hold on   %在一张图上同时显示多组数据
elseif size(Epa,2)>=3%如果Epa的列数大于3  即目标函数值的个数有三个以上
    %%plot pareto fronts in 3-D
    f1=Epa(:,1);f2=Epa(:,2);f3=Epa(:,3);
    plot3(f1,f2,f3,'kd'); %以黑色菱形来显示目标函数值的分布情况
    grid on; %显示网格
    xlabel('Function 1');%x坐标轴的标记为Function 1
    ylabel('Function 2');%y坐标轴的标记为Function 2
    zlabel('Function 3');%z坐标轴的标记为Function 3
    hold on     %在一张图上同时显示多组数据
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%