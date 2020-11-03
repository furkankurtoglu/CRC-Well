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
names={'Oxygen','Glucose','Glutamine,','Lactate'};
Save_MicEnv = 'Y';
Save_CMicEnv = 'N';

%%
filename1="MicroEnv.gif";
filename2="CMicroEnv.gif";
for i = 1:length(OutMatFiles)+1
%for i = 1:1

    if (Save_MicEnv == 'Y')
    out = read_microenvironment(strcat(OutMatFiles{i}));
    h = figure(1);
    set(gcf, 'Position',  [0, 0, 1344, 756])
    plot_microenvironment(out,names)
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
        print('hey')
    else
        imwrite(imind,cm,filename1,'gif','WriteMode','append');
    end
    end
    
    if (Save_CMicEnv == 'Y')
    load(CMicEnvOutMatFiles{i,1});
    centers=(multiscale_microenvironment(1,:));
    j=figure(2);
    set(gcf, 'Position',  [0, 0, 1344, 756])
    subplot(2,2,1);
    plot(centers,multiscale_microenvironment(5,:))
    xlim([-500 10200])
    title('oxygen')
    subplot(2,2,2);
    plot(centers,multiscale_microenvironment(6,:))
    xlim([-500 10200])
    title('glucose')
    subplot(2,2,3);
    plot(centers,multiscale_microenvironment(7,:))
    xlim([-500 10200])
    title('glutamine')
    subplot(2,2,4);
    plot(centers,multiscale_microenvironment(8,:))
    SubGroupTitle=strcat('Time =',num2str(i), ' min');
    suptitle(SubGroupTitle)
    xlim([-500 10200])
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
    
end 
%%
cd ..
