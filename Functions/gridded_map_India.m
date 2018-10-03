function [pl] = gridded_map_India(X,Y,inc,cm,pl);
%X and Y are the center coordinates of gridcells, cm is a nx3
%matrix that specifies the color for each pair X,Y. inc is the
%width of the box, pl is the subplot handle

deg=struct();
for ct1=1:length(X);
    x=X(ct1);
    y=Y(ct1);
    deg(ct1).X=[x-inc x-inc x-inc x x+inc x+inc x+inc x x-inc];
    deg(ct1).Y=[y-inc y y+inc y+inc y+inc y y-inc y-inc y-inc];
end

%call the axes
axes(pl);
hold on;
for ct1=1:length(deg)
    Cs=cm(ct1,:);
    lb=patch(deg(ct1).X, deg(ct1).Y,Cs,'EdgeColor','none');
end 
load coast;
%plot(long,lat,'linewidth',2,'Color',[0.8 0.8 0.8]);
grid on; box on;

