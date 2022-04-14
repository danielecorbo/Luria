%%%%%%%this script, created by Daniele Corbo, automatically calculates the kinematics of a motor task recorded with a camera.  
%%%%%%%%%%%%%%%%%
clc
clear
%define the directory
ro_dir='/*';
files=dir(['*.avi']) ;

for t=1:size(files)
fname=files(t).name;
info= VideoReader([fname]); 
sx = info.Width; 
sy = info.Height;
first_frame = 1;
last_frame = info.NumberOfFrames  %info.NumFrames
new_sx = round(sx/1); 
new_sy = round(sy/1);
flow = zeros([new_sy new_sx 2 last_frame]); 
n_frames = 5; %can be changed
n_scales = 6; % can be changed
thres_lin = 1; % can be changed
nc_min = 4; % can be changed
filtname = 'Gt11B0.0833f0.25.mat'; 
stable = 0; % can be changed

tic;
%keyboard
tmp1 = read(info); 
for start_frame = (first_frame-1):(last_frame-n_frames) 
  II = zeros([ new_sy new_sx n_frames ]);
  for frame = 1:n_frames  
    
    tmp= tmp1(:,:,:,start_frame+frame);
   
    if (size(tmp,3)==3)
      tmp = double(rgb2gray(tmp)) ./ 255; % rgb to gray
    else
      tmp = double(tmp) ./ 255;
    end
   tmp = imresize_old(tmp,1,'nearest',0);  

    II(:,:,frame) = tmp(end:-1:1,:); 
  
  end
  
 disp(start_frame);
  O = (ctf_optic_flow(II,n_scales,thres_lin, ...   
                     nc_min,filtname,stable,[],[]));
  O(isnan(O)) = 0; 
  O = cat(3,medfilt2(O(:,:,1),[11 11]),medfilt2(O(:,:,2),[11 11])); 
  tmp = single(complex(O(:,:,1),O(:,:,2))); 
  flow(:,:,:,start_frame+3) = O; 
  fprintf('frame %04d/%04d - %3.0f s\r', ...
          start_frame+3,last_frame-n_frames+3,toc);
  tic;
end
flow_file=[fname(1:length(fname)-4) '_flow'];
fprintf('\n');
save(flow_file,'flow');
clear flow tmp II
end

%%%%%calculate motion/degree from flow
clear;clc
%%%%define directory
ro_dir='/*';
files=dir(['*.mat']) ;
for r=1:size(files)
fname=files(r).name;
load(fname)
t= size(flow);
%%%%%%%%%%%resize flow%%%%%%%%%%%%%%%%%%
W=imresize(flow,[240 400]);
t= size(W);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_x(:,:,:)= W(:,:,1,:);
A_y(:,:,:)= W(:,:,2,:);
a_x=mean(A_x,3);
a_y=mean(A_y,3);
A_xx= flipdim(a_x,1);
A_yy=flipdim(a_x,1);

%% 
for i=1:t(1)
    for u= 1:t(2)
        M(i,u)= sqrt(A_xx(i,u)^2+A_yy(i,u)^2);
    end
end
figure; imagesc(M);
moto=mean(mean(M));
%% 
moto_file=[fname(1:length(fname)-14) '_motion'];
    eval(['fname' num2str(r) '=motion;'])
clear a_x A_x A_xx a_y A_yy i M A_y 
save kinematics_.txt fname* -ascii -double -tabs
end




