function [fitresult, gof] = createFit3(identity_com, tmscore_com)
%CREATEFIT3(IDENTITY_COM,TMSCORE_COM)
%  创建一个拟合。
%
%  要进行 '无标题拟合 1' 拟合的数据:
%      X 输入: identity_com
%      Y 输出: tmscore_com
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 05-May-2023 14:01:51 自动生成


%% 拟合: '无标题拟合 1'。
[xData, yData] = prepareCurveData( identity_com, tmscore_com );

% 设置 fittype 和选项。
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.987276742091642 0.128137317156022 -0.000565109537689303];

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% 绘制数据拟合图。
figure( 'Name', '无标题拟合 1' );
h = plot( fitresult, xData, yData );
legend( h, 'tmscore_com vs. identity_com', '无标题拟合 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 为坐标区加标签
xlabel( 'Sequence identity', 'Interpreter', 'none' );
ylabel( 'TM-score', 'Interpreter', 'none' );
grid on


