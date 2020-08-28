clear
close all
clc

%%
cd output/


%%
s=what;
MatFiles = s.mat;
OutMatFiles = MatFiles(contains(MatFiles,'microenvironment0'));
OutMatFiles(1) = [];
OutMatFiles(1) = [];
CMicEnvOutMatFiles = MatFiles(contains(MatFiles,'microenvironment1'));
names={'Oxygen','Glucose','Lactate'};
%%
filename1="MicroEnv.gif";
filename2="CMicroEnv.gif";
%for i = 1:length(OutMatFiles)
for i = 1:45

    out = read_microenvironment(strcat(OutMatFiles{i}));
    h = figure(1);
    set(gcf, 'Position',  [0, 0, 5000, 4000])
    plot_microenvironment(out,names)
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename1,'gif','WriteMode','append');
    end

    
    load(CMicEnvOutMatFiles{i,1});
    j=figure(2);
    set(gcf, 'Position',  [0, 0, 1800, 1000])
    subplot(2,2,1);
    plot(multiscale_microenvironment(5,:))
    title('oxygen')
    subplot(2,2,2);
    plot(multiscale_microenvironment(6,:))
    title('glucose')
    subplot(2,2,3);
    plot(multiscale_microenvironment(7,:))
    title('lactate')
    frame = getframe(j);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename2,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename2,'gif','WriteMode','append');
    end
    
    
end 
%%
%for i = 1:length(CMicEnvOutMatFiles)