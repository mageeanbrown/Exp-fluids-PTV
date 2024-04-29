clear all

% Specify the path to the video file
% videoFilePath = "C:\Users\rossh\Documents\Brown\Classes\2024Spring\ENGN2912T - Experiemental Fluid Mechanics\Labs\HW05 PTV and Microrheology\Videos\group 3 1 um 25 hz fluid 1.avi";
videoFilePath = "C:\Users\rossh\Documents\Brown\Classes\2024Spring\ENGN2912T - Experiemental Fluid Mechanics\Labs\HW05 PTV and Microrheology\Videos\group 4 1 um 25 hz fluid 2.avi";
% videoFilePath = "C:\Users\rossh\Documents\Brown\Classes\2024Spring\ENGN2912T - Experiemental Fluid Mechanics\Labs\HW05 PTV and Microrheology\Videos\group 3 200 nm 25 hz fluid 1.avi";

% Create a VideoReader object to read the video file
videoReader = VideoReader(videoFilePath);

% Get video properties
numFrames = videoReader.NumberOfFrames;
frameHeight = videoReader.Height;
frameWidth = videoReader.Width;

% Initialize the 3D matrix to store the video frames
a = zeros(frameHeight, frameWidth, numFrames, 'uint8');

% Read frames and store them in the 3D matrix
for idx = 1:numFrames
    frame = readFrame(videoReader);
    a(:,:,idx) = rgb2gray(frame);
end

% Display a message
disp('Video file has been converted to a 3D matrix successfully.');
%% Use Videos: group 3 1 um 25 hz fluid 1.avi, group 4 1 um 25 hz fluid 2.avi
% These were the parameters for the videos I found worked best...

lobjectList = [29,16];
extentList = [25,19];
massList = [200000,115000];

% For vid group 3... set z = 1
% z = 1;
% For vid group 4... set z = 2
z = 2;
%%
% Input Parameters:
n = 1;
a1 = a(:,:,n);  % The n-th image of the movie

lnoise=1;       % The high freq cutoff 


lobject=lobjectList(z);      % The low freq cutoff (try diameter of particles in pixels)

% The output is a filtered version of the input image ... 
b = bpass(a1,lnoise,lobject); 

% To check the filtered image, display it ...
colormap('gray'); imagesc(b); drawnow
% It should have nice bright blobs (the particles) on a black background. 
%%
% Input Parameters:
InputImage=b; % Finds features in a single image. The input image is the output image of bpass.
extent = extentList(z);   % Slightly larger than the size of the blobs in bpass  
f = findfeatures(InputImage,extent,masscut=massList(z));
%,minpix=150);
% See header for more information about input parameters.  The 'min' parameter sets a minimum 
% peak brightness needed to identify a blob as a feature.  
% disp(f(1:size(f,1),:))
%%
% Display the image (a1) and the found features in that image (f).  
fo=fover2d(a1, f); drawnow

% Using some keywords can make it easier to see ...
figure(10)
fo=fover2d(a1, f, circle='y', radius=40); drawnow
%%
pt_all=epretrack(a,bplo=lnoise,bphi=lobject,extent=extent,mass=massList(z));
disp(size(pt_all))
%%
% titleBeg = "Group 3: 1 $\mu$m, 25 Hz, Fluid 1: ";
titleBeg = "Group 4: 1 $\mu$m, 25 Hz, Fluid 2: ";
figure(2)
plot(pt_all(:,1),pt_all(:,2),'.'); drawnow
figure(3)
histogram(pt_all(:,3))
title(titleBeg + "Brightness", "Interpreter","latex","FontSize",16); drawnow
xlabel("Brightness (summed pixel brightness)")
ylabel("Number of Found Particles")
figure(4)
histogram(pt_all(:,4))
title(titleBeg + "Size", "Interpreter","latex","FontSize",16); drawnow
xlabel("Size [pixels]")
ylabel("Number of Found Particles")
figure(5)
histogram(pt_all(:,5))
title(titleBeg + "Eccentricity", "Interpreter","latex","FontSize",16); drawnow
xlabel("Eccentricity")
ylabel("Number of Found Particles")


%% Examing particle properties across time series
figure(6)
plot(pt_all(:,3),pt_all(:,4),'.'); %drawnow 
xlabel('Brightness')
ylabel('Size')

figure(7)
plot(pt_all(:,3),pt_all(:,5),'.'); %drawnow 
xlabel('Brightness')
ylabel('Eccentricity')

figure(8)
plot(pt_all(:,4),pt_all(:,5),'.'); %drawnow 
xlabel('Size')
ylabel('Eccentricity')

%%Everything below this is still a WPI

% %%
% w=find(pt_all(:,5) > 0.1); % index of particles with an eccentricity > 0.01 to check for non-cicular particles
% fo=fover2d(a(:,:,1:100),pt_all(w,:),circle='y',radius=40,big='y'); % Make a movie of those particles 
%%
w=find(pt_all(:,5) < 0.15); % index of particles with an eccentricity > 0.1
pt=pt_all(w,:); % New array of 'good' particles
% *** CHECK ***
% Look at the size of pt to see how many particles were found ... 
plot(pt(:,1),pt(:,2),'.'); drawnow
%%
figure()
histogram(mod(pt(:,1:2),1),20); drawnow
title(titleBeg + "Pixel Bias", "Interpreter","latex","FontSize",16); drawnow
xlabel("Interpolated Fractional Part of Center of Particle")
ylabel("Number of Found Particles")

%%

t=track(pt_all,6,dim=2,memory=2,goodenough=10);

figure(6)
plottr(t);

figure(7)
plottr(t, ID=5);

figure(9)
% m=msd(t);

m=msd(t,timestep=1/25,minN=200);


figure(10)
res=linearfit(m(:,1),m(:,6));

%%

msdArray = zeros(numFrames,4); % 1st is deltaT, 2nd Column is count, 3rd is msd, 4th std,
msdArray(:,1) = 1:numFrames;

%% msd -var

indParticles = unique(t(:,7)); % Numbers of all the particles that are found
particleTracks = cell(1,length(indParticles));
for z1 = 1:length(indParticles)
    logicalFound = t(:,7) == indParticles(z1);
    particleTracks{z1} = t(logicalFound,:);
end

%%
for r = 1:length(particleTracks)
    for q = 1:length(particleTracks{r})-1
        for p = q+1:length(particleTracks{r})
            t_1 = particleTracks{r}(q,6); t_2 = particleTracks{r}(p,6); 
            msdCur = (particleTracks{r}(q,1) - particleTracks{r}(p,1))^2 + (particleTracks{r}(q,2) - particleTracks{r}(p,2))^2;
            deltaT = t_2-t_1;
            msdArray(deltaT,2) = msdArray(deltaT,2) + 1;
            msdArray(deltaT,3) = msdArray(deltaT,3) + msdCur;
        end
    end
end
msdArray(:,3) = msdArray(:,3)./msdArray(:,2);

%%
for r = 1:length(particleTracks)
    for q = 1:length(particleTracks{r})-1
        for p = q+1:length(particleTracks{r})
            t_1 = particleTracks{r}(q,6); t_2 = particleTracks{r}(p,6); 
            stdPartCur = (particleTracks{r}(q,1) - particleTracks{r}(p,1))^2 + (particleTracks{r}(q,2) - particleTracks{r}(p,2))^2;
            deltaT = t_2-t_1;
            % msdArray(deltaT,2) = msdArray(deltaT,2) + 1; % Do not want to
            % double count...
            msdArray(deltaT,4) = msdArray(deltaT,4) + (msdArray(deltaT,3) - stdPartCur)^2;
        end
    end
end
msdArray(:,4) = msdArray(:,4)./msdArray(:,2);
msdArray(:,4) = sqrt(msdArray(:,4));


%% Remove places with no count
msdArrayFin = msdArray(~msdArray(:,2) == 0, :);
%%
frameRate = 25; % Hz
msdArrayFin(:,1) = msdArrayFin(:,1)./frameRate;
group3_1um_pixConv = (1/(19.75*10^6))^2;
group4_1um_pixConv = (1/(17.5939*10^6))^2;

msdArrayFin(:,3:4) = msdArrayFin(:,3:4) * group3_1um_pixConv;

%%
hex = ["#4c0e59","#47c9c1"];
if length(msdArrayFin) > 475
    poi = 1:400;
else
    poi = 1:length(msdArrayFin)-3;
end
poi = 1:355;
% poi = 1:30; % Set to the values you want to display (good to limit as it gets weird after a certain point)
t_s = msdArrayFin(poi,1);
msdPlot = msdArrayFin(poi,3);
lowerStd = msdArrayFin(poi,3)-msdArrayFin(poi,4)/2;
upperStd = msdArrayFin(poi,3)+msdArrayFin(poi,4)/2;
%% Fit linear slope to log-spaced data

polyLin = polyfit(log10(t_s), log10(msdPlot), 1);
% polyLin = polyfit(t_s, msdPlot, 1);
YpolyLin = polyval(polyLin, log10(t_s));
foundSlope = polyLin(1);
foundInter = polyLin(2);

%% Main Figure Plot, MSD standart dev and fitted linear slope
figure()
loglog(t_s,msdPlot,"LineWidth",2.5,"Color",hex(1))
hold on 
% loglog(t_s,lowerStd,"LineWidth",0.0025,"Color",hex(1))
% loglog(t_s,upperStd,"LineWidth",0.0025,"Color",hex(1))
patch([t_s; flipud(t_s)],[lowerStd; flipud(upperStd)], 'm', 'FaceAlpha',0.2, 'EdgeColor','none');
% patch([t_s; flipud(t_s)],[lowerStd; flipud(upperStd)], "FaceColor", 'blue', 'FaceAlpha',0.2, 'EdgeColor','none');
semilogx(t_s, 10.^(YpolyLin),"LineWidth",2.5,'LineStyle',":","Color",hex(2))
title("MSD($\Delta t$)","Interpreter","latex","FontSize",16)

subtitle([titleBeg,"Slope: " + num2str(foundSlope)], "Interpreter","latex","FontSize",16)            % Change SubTitle for Specific Videos

ylabel("MSD($\Delta t$) [m$^2$]","Interpreter","latex","FontSize",16)
xlabel("$\Delta t$ [s]","Interpreter","latex","FontSize",16)
xlim([0,inf])
ylim([min(lowerStd),inf])
legend(["MSD","Standard Deviation","Fitted Line"],"Location","best")

%% Calculate Viscosity of Fluid based on Stokes-Einstein Eq
diffusionCoef = 10^foundInter/2;
% diffusionCoef = foundSlope/2;
display(titleBeg,"Interp")
display("Diffusion Coeffiencent, D: " + num2str(diffusionCoef))
boltz = 1.380649 * 10^(-23); % Boltzman Constant
temp = 298.15; % Absolute Temperature (Assuming 25 Celcius room temp)
rad = ( 1*10^(-6) )/2; % Radius of particle - (change depending on video)

visc = (boltz*temp)/(6*pi*diffusionCoef*rad);
display("Viscosity, u: " + num2str(visc))





