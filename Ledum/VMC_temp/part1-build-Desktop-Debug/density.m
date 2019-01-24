%% 
clear all;
close all;
clc;


%% Load in binary data file
x0 = fopen('x0.bin');
x1 = fopen('x1.bin');
x2 = fopen('x2.bin');
x3 = fopen('x3.bin');
x4 = fopen('x4.bin');
x5 = fopen('x5.bin');
x6 = fopen('x6.bin');
x7 = fopen('x7.bin');
x8 = fopen('x8.bin');
x9 = fopen('x9.bin');

x0 = fread(x0, 'double');
x1 = fread(x1, 'double');
x2 = fread(x2, 'double');
x3 = fread(x3, 'double');
x4 = fread(x4, 'double');
x5 = fread(x5, 'double');
x6 = fread(x6, 'double');
x7 = fread(x7, 'double');
x8 = fread(x8, 'double');
x9 = fread(x9, 'double');

y0 = fopen('y0.bin');
y1 = fopen('y1.bin');
y2 = fopen('y2.bin');
y3 = fopen('y3.bin');
y4 = fopen('y4.bin');
y5 = fopen('y5.bin');
y6 = fopen('y6.bin');
y7 = fopen('y7.bin');
y8 = fopen('y8.bin');
y9 = fopen('y9.bin');

y0 = fread(y0, 'double');
y1 = fread(y1, 'double');
y2 = fread(y2, 'double');
y3 = fread(y3, 'double');
y4 = fread(y4, 'double');
y5 = fread(y5, 'double');
y6 = fread(y6, 'double');
y7 = fread(y7, 'double');
y8 = fread(y8, 'double');
y9 = fread(y9, 'double');

randx = [x0; x1; x2; x3; x4; x5; x6; x7; x8; x9];
randxx = [x0(length(x0)/10:end); 
          x1(length(x0)/10:end); 
          x2(length(x0)/10:end); 
          x3(length(x0)/10:end); 
          x4(length(x0)/10:end); 
          x5(length(x0)/10:end); 
          x6(length(x0)/10:end); 
          x7(length(x0)/10:end);
          x8(length(x0)/10:end);
          x9(length(x0)/10:end)];
randy = [y0; y1; y2; y3; y4; y5; y6; y7; y8; y9];
randyy = [y0(length(x0)/10:end); 
         y1(length(x0)/10:end); 
         y2(length(x0)/10:end); 
         y3(length(x0)/10:end); 
         y4(length(x0)/10:end); 
         y5(length(x0)/10:end); 
         y6(length(x0)/10:end); 
         y7(length(x0)/10:end);
         y8(length(x0)/10:end);
         y9(length(x0)/10:end)];

minx = min(randx);
maxx = max(randx);
miny = min(randy);
maxy = max(randy);

randxx = [x0(length(x0)/10:end); 
          x1(length(x0)/10:end)];
randyy = [y0(length(x0)/10:end); 
         y1(length(x0)/10:end)];

colormap('gray');
% randx =randn(100000,1);
% randy =randn(100000,1);
[z,N] = hist3([randyy randxx],[800 800]);
% surf(x)

xb = linspace(min(randxx),max(randxx),size(z,1)+1);
yb = linspace(min(randyy),max(randyy),size(z,1)+1);

imagesc(xb,yb,max(max(z))-z);
xlabel('x');
ylabel('y');

minx = min(randx);
maxx = max(randx);
miny = min(randy);
maxy = max(randy);

%%
    


randx = [x0(length(x0)/10:end) x1(length(x0)/10:end)];
randy = [y0(length(x0)/10:end) y1(length(x0)/10:end)];

minx = min(min(randx));
maxx = max(max(randx));
miny = min(min(randy));
maxy = max(max(randy));

j = [1 2];
for k=1:2
    i = j(k);
    h = figure(k);
%     subplot(1,2,k);
    colormap('gray');
    
    randxx = [randx(:,i); minx; maxx];
    randyy = [randy(:,i); miny; maxy];
    
    [z,N] = hist3([randyy randxx], [400 400]);
    % surf(x)

    xb = linspace(min(randxx),max(randxx),size(z,1)+1);
    yb = linspace(min(randyy),max(randyy),size(z,1)+1);

    imagesc(xb, yb, max(max(z))-z); %max(max(z))-
    
    axis equal
    
    % Adjusts white space around figures for better exporting to LaTeX
    % documents.

    ti = get(gca,'TightInset');                                           
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);         
    set(gca,'units','centimeters')                                         
    pos = get(gca,'Position');                                             
    ti = get(gca,'TightInset');                                            
    set(gcf, 'PaperUnits','centimeters');                                  
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);        
    set(gcf, 'PaperPositionMode', 'manual');                               
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]); 
    %     axis equal
    h = figure(k);
    str = sprintf('~/Comp2/Figurer/HeAtomWITHJASTROWParticle%d.pdf', k);
    saveas(gcf, str,'pdf')
    
end


%%

n = 1;

X = sqrt((x0(n:end) - x1(n:end)).^2+(y0(n:end) - y1(n:end)).^2+(y2(n:end)-x2(n:end)).^2);
mean(X)

% He
% - uten jastrow <r_ij> = 1.1973 (1.1973 hvis vi kapper av 10%).
% - med jastrow  <r_ij> = 1.3876 (1.3856 hvis vi kapper av 10%).
