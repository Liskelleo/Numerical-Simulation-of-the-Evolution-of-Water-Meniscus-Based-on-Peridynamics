% 本脚本用于单液桥简化模型的几何建模
% 作者: Liskelleo
clc; clear all;
fid=fopen('soil_board.dat','w+'); % 打开文件用于写入数据
fida=fopen('air_backround.dat','w+'); % 打开文件用于写入数据
lat_c=0.01; % 网格单元的尺寸
nx=55; % x方向的网格点数
ny=101; % y方向的网格点数
x0=-0.02; % x方向的起始坐标
y0=0; % y方向的起始坐标
coord=zeros(nx*ny,4); % 初始化坐标数组
coord2=zeros(nx*ny,4); % 初始化背景数组
ii=0; ai=0; % 初始化索引计数器
for i=1:nx
    for j=1:ny
        xi=x0+(i-1)*lat_c; % 计算x坐标
        yi=y0+(j-1)*lat_c; % 计算y坐标
        if xi<=0 ||  xi>=0.5 % 检查坐标是否为固体
            ii=ii+1;    
            coord(ii,:)=[ii,xi,yi,2]; % 存储索引和坐标，固体虚拟层标记为2
        else
            ai=ai+1;    
            coord2(ai, :)=[ai,xi yi,0]; % 存储索引和坐标，空气背景层标记为0
            if yi >=0.2 && yi<=0.8
                ii=ii+1;
                coord(ii,:)=[ii,xi,yi,1]; % 存储索引和坐标，流体粒子标记为1
            end
        end
    end
end 
npart=ii; apart=ai; % 根据ii,ai的最终值设置总粒子个数npart,apart
 
%初始化 zuobiao 数组，用于存储流固粒子的索引和坐标
zuobiao=zeros(npart,4);
zuobiao(:,1)=1:npart;
zuobiao(1:npart,2:4)=coord(1:ii,2:4);
%fprintf(fid,'%d\t ',npart); % 写入粒子数量
%fprintf(fid,'%d\t ',lat_c);  % 写入粒子尺寸
%fprintf(fid,'\n');

for i=1:npart
    b=zuobiao(i,1);
    c=[zuobiao(i,2) zuobiao(i,3)];
    d=zuobiao(i,4);
    fprintf(fid,'%d\t',b);         % 写入每个粒子的索引编号
    fprintf(fid,'%15.5e\t',c);   % 写入每个粒子的位置坐标
    fprintf(fid,'%d\t',d);         % 写入每个粒子的边界信息
    fprintf(fid,'\n');
end

%初始化 azuobiao 数组，用于存储空气粒子的索引和坐标
azuobiao=zeros(apart,4);
azuobiao(:,1)=1:apart;
azuobiao(1:apart,2:4)=coord2(1:ai,2:4);
%fprintf(fida,'%d\t ',apart); % 写入粒子数量
%fprintf(fida,'%d\t ',lat_c);  % 写入粒子尺寸
%fprintf(fida,'\n');

for i=1:apart
    ab=azuobiao(i,1);
    ac=[azuobiao(i,2) azuobiao(i,3)];
    ad=azuobiao(i,4);
    fprintf(fida,'%d\t',ab+npart);         % 写入每个粒子的索引编号
    fprintf(fida,'%15.5e\t',ac);   % 写入每个粒子的位置坐标
    fprintf(fida,'%d\t',ad);         % 写入每个粒子的边界信息
    fprintf(fida,'\n');
end

% 根据粒子位置坐标寻找边界粒子
point=find((zuobiao(:,4)==2));
mpart=length(point);
bc=zeros(mpart,4);
bc(:,1)=zuobiao(point,1);
bc(1:mpart,2:3)=zuobiao(point,2:3);

point2=find((zuobiao(:,4)==1));
mpart2=length(point2);
bcf=zeros(mpart2,4);
bcf(:,1)=zuobiao(point2,1);
bcf(1:mpart2,2:3)=zuobiao(point2,2:3);

x=azuobiao(:,2);
y=azuobiao(:,3);
scatter(x,y,'g','filled');
hold on
xf=bcf(:,2);
yf=bcf(:,3);
scatter(xf,yf,'b');
hold on
x1=bc(:,2);
x2=bc(:,3);
scatter(x1,x2,'r','filled');
axis equal %daspect([1 1 1]);
hold off
