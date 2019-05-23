% script to load tow raster files, trim all areas not overlapping, then
% laign in cell aray

%% input params
clear
files=cellstr(ls([dir_in, '*.tif']));
file_in=[dir_in, files{n}];

%% load files
fprintf('Processing file:\n\t%s\n', files{n})

    % choose input file A
[bw, bw_R]=geotiffread(file_in);
disp('Importing file...')

    % Choose input file B
test_name=['WC',files{n}(2:end)];
test_in=[im_dir_in, test_name];

[test, test_R]=geotiffread(test_in);
test=test==1;
    % clip extents to match
im={~bw, test};
    %% prep to pad smaller file to align pixels in geographic locationshelp
    % find out which file is most NW (assuming cols start from N and rows
    % start from W
R=[bw_R; test_R];
    % find which image to pad on each side
[top,north]=min(vertcat(R(:).YWorldLimits)); northFile=north(1); % file with most N upper extent
    top=max(top);
[left,west]=max(vertcat(R(:).XWorldLimits)); westFile=west(1); % file with most W left extent
    left=min(left);
[btm,south]=max(vertcat(R(:).YWorldLimits)); southFile=south(2); % file with most W left extent
    btm=min(btm);
[right,east]=min(vertcat(R(:).XWorldLimits)); eastFile=east(2); % file with most W left extent
    right=max(right);

    % compute pad
for i=1:length(im)
    [NW_intrins(i,1), NW_intrins(i,2)] = worldToIntrinsic(R(i),left, top); % outside bounds. order: x,y
    [SE_intrins(i,1), SE_intrins(i,2)] = worldToIntrinsic(R(i),right, btm);
    buf(i,1)=NW_intrins(i,2)-R(i).YIntrinsicLimits(1); % top buffer
    buf(i,2)=NW_intrins(i,1)-R(i).XIntrinsicLimits(1); % left buffer
    buf(i,3)=-SE_intrins(i,2)+R(i).YIntrinsicLimits(2); % btm buffer
    buf(i,4)=-SE_intrins(i,1)+R(i).XIntrinsicLimits(2); % right buffer
end
    % pad images and place in cell array
    buf=round(buf);
for i=1:length(im)
    im{i}=im{i}(1+buf(i,1):end-buf(i,3), 1+buf(i,2):end-buf(i,4));
end