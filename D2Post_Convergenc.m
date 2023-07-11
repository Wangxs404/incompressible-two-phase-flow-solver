%% 绘制收敛曲线
close all

X=0.4./[128,64,32,16];
Y=[ 0.00031766,0.00070658,0.0020242,0.01812 ];

loglog(X, Y,'b-s','LineWidth',1.5,'MarkerSize',10);
hold on

tr1 = [1e-2,0.01812]; tr2 = 0.003125; ratio = 2;
triangle(tr1,tr2,ratio,'1','2',1.2);

tr1 = [0.0075,0.00031766]; tr2 = 0.025; ratio = 1;
triangle(tr1,tr2,ratio,'1','1',1.2);

grid on;
axis square;

xlim([1e-3,1e-1]);
ylim([1e-4,1e-1]);
xlabel('D/h')
ylabel('L∞')

function triangle(tr1,tr2,ratio,text1,text2,pos_adjust)
    x0 = tr1(1); y0 = tr1(2); x1 = tr2; 
    y1 = y0*(x1/x0)^ratio;
    loglog([x0,x1],[y0,y0],'-k','HandleVisibility','off','Linewidth',1);
    loglog([x1,x1],[y0,y1],'-k','HandleVisibility','off','Linewidth',1);
    loglog([x1,x0],[y1,y0],'-k','HandleVisibility','off','Linewidth',1);
    text(sqrt(x0*x1),y0*pos_adjust,text1,'fontsize',10);
    text(x1*pos_adjust,sqrt(y0*y1),text2,'fontsize',10);
end

%  三角形参数说明
%         tr1：水平边和斜边交点坐标
%         tr2：水平边另一个端点坐标
%         ratio：斜边斜率
%         text1：水平边文字
%         text2：垂直边文字
%         pos adjust：文字坐标的微调
