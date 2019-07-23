hfig = figure();
haxes1 = subplot(1,2,1, 'Parent', hfig);
haxes2 = subplot(1,2,2, 'Parent', hfig);

x = 0:0.01:4*pi;
hplot1a = plot(haxes1,x,sin(x));
hold(haxes1,'on');
hplot2a = plot(haxes1,x,cos(x));
hold(haxes1,'on');

hplot1b = plot(haxes2,x,sin(x+pi/4));
hold(haxes2,'on');
hplot2b = plot(haxes2,x,cos(x+pi/4));
hold(haxes2,'on');

surf3tikz(hfig,'test_2d_multi_axes');