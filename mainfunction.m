clear all

% Specify the path to the video file
videoFilePath = '/Users/adampoche/Downloads/group 4 200 nm 25 hz fluid 2.avi';

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
%%
% Input Parameters:
n = 1;
a1 = a(:,:,n);  % The n-th image of the movie

lnoise=1;       % The high freq cutoff 


lobject=35;      % The low freq cutoff (try diameter of particles in pixels)

% The output is a filtered version of the input image ... 
b = bpass(a1,lnoise,lobject); 

% To check the filtered image, display it ...
colormap('gray'); imagesc(b); drawnow
% It should have nice bright blobs (the particles) on a black background. 
%%
% Input Parameters:
InputImage=b; % Finds features in a single image. The input image is the output image of bpass.
extent = 35;   % Slightly larger than the size of the blobs in bpass  
f = findfeatures(InputImage,extent,masscut=90000);
%,minpix=150);
% See header for more information about input parameters.  The 'min' parameter sets a minimum 
% peak brightness needed to identify a blob as a feature.  
disp(f(1:size(f,1),:))
%%
% Display the image (a1) and the found features in that image (f).  
fo=fover2d(a1, f); drawnow

% Using some keywords can make it easier to see ...
figure(1)
fo=fover2d(a1, f, circle='y', radius=40); drawnow
%%
pt_all=epretrack(a,bplo=lnoise,bphi=lobject,extent=15,mass=90000);
disp(size(pt_all))
%%
figure(2)
plot(pt_all(:,1),pt_all(:,2),'.'); drawnow
figure(3)
histogram(pt_all(:,3))
title('Brightness'); drawnow
figure(4)
histogram(pt_all(:,4))
title('Size'); drawnow
figure(5)
histogram(pt_all(:,5))
title('Eccentricity'); drawnow


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

% %%
% w=find(pt_all(:,5) > 0.1); % index of particles with an eccentricity > 0.01 to check for non-cicular particles
% fo=fover2d(a(:,:,1:100),pt_all(w,:),circle='y',radius=40,big='y'); % Make a movie of those particles 
%%
w=find(pt_all(:,5) < 0.06); % index of particles with an eccentricity > 0.1
pt=pt_all(w,:); % New array of 'good' particles
% *** CHECK ***
% Look at the size of pt to see how many particles were found ... 
plot(pt(:,1),pt(:,2),'.'); drawnow
%%
histogram(mod(pt(:,1:2),1),20); drawnow

%%

t=track(pt_all,6,dim=2,memory=3,goodenough=10);

figure(6)
plottr(t);

figure(7)
plottr(t, ID=5);

figure(9)
m=msd(t);

figure(10)
res=linearfit(m(:,1),m(:,6));