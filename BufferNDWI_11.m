function labeled=BufferNDWI_11(cir, InitialBuffer);
% Next finalize the hole filtering and start validating!
% Playing around with structering element size and filtering inital guess
% Next, rewrite all calcs into structure, or remove
% GoF metric for Otsu thresh --part 2
% Code takes loaded image, divides into NDWI-delineated regions, labels...
% them, creates external buffer.
% THIS ONE uses global thresh 
% Initial buffer is optional and specifies whether to buffer and recacl
% thresh before loop

%%Load and display

%%%%%%%%%%%%%%%%%%%%%%%% non-function paramters
% clear
% bw_in='D:\ArcGIS\FromMatlab\Adapt_Clip2.tif';
% img_in='D:\ArcGIS\FromMatlab\Clip_Square_2.tif'
% img_in='D:\ArcGIS\out\20170816 Ortho Segment_3-Pond25-26-clip'
% img_in='D:\AboveData\Ortho 20170806 Canada Albers.tif';
% cir=imread(img_in);
% [bw, bw_R]=geotiffread(bw_in);
% [cir, cir_R]=geotiffread(img_in);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' -Loading data-')
if~exist('InitialBuffer')
    InitialBuffer=false;
end
% cir=cir(1:5000, 1:5000, :); %<---------- REMOVE
% tile=123-16; %tile to import % used to be 124
% cir=importchunks(img_in, 16, 16, tile);
% cir=cir.tile;

%% Params %%%%%%%%%%%%%%%%%%%%%%%%%
MinSize=25; %Min region size to classify as open water
Len=5000; %Size to split river reaches
Num=10; %Number of dilations per river reach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check for empty tile and compute ndwi
% cir=im2double(cir);
NoValues=or(cir(:,:,1)==192 & cir(:,:,2)==192 & cir(:,:,3)==192,...
    cir(:,:,1)==255 & cir(:,:,2)==255 & cir(:,:,3)==255);
[a,b]=size(cir(:,:,1));
Switch=true;
Black=0;
if length(unique(cir))<=2
    Black=1;
    disp(upper('This tile is outside image region (black pixels)'))
    Switch=false;
end
cir=im2single(cir);
cir_ndwi=(cir(:,:,3)-cir(:,:,1))./(cir(:,:,3)+cir(:,:,1));
bw=imbinarize(cir_ndwi, 'adaptive', 'Sensitivity', 0.75);
clear cir
disp(' -Data loaded-')

%% Initial buffer of all lakes 
% Relabel
if Switch & InitialBuffer
    [L, j]=bwlabel(bw, 4);
    stats = regionprops(L, 'PixelList', 'Centroid',...
        'Area', 'BoundingBox', 'Perimeter');
    per=sum([stats.Perimeter]);
    ar=sum([stats.Area]);
    buf=1.5*ar/per;
    if buf<=0 | buf >1000 | isinf(buf) | isnan(buf)
        buf=1;
    end
    SE = strel('diamond',round(buf));
    disp('Running initial dilation/threshold on guess...')
    bw_dil_lakeN=logical(imdilate(L, SE));
    [level, EM] = graythresh((cir_ndwi(bw_dil_lakeN)));
    bw=cir_ndwi>=level;
end
[L, j]=bwlabel(bw, 4);
stats = regionprops(L, 'PixelList', 'Centroid',...
    'Area', 'BoundingBox', 'Perimeter');
%%%%%%%%%%%%%%shortcut
% labeled=imread('D:\ArcGIS\FromMatlab\LocalThreshSplitClip2.tif');
%%%%%%%%%%%%%%%%%%%%%%

%% prelims before loop

% Geometry
CircPerim=2*pi*sqrt([stats.Area]/pi) ;
PerSin=CircPerim./[stats.Perimeter];

%% Loop
disp('iterating...')
labeled=zeros(a,b);
% fileID = fopen('log.txt','a');

tic
j=1; %counter
SE = strel('diamond',10);

% Recorded vars: buf, ratios, buf area,

% for i=632:634
i=1; %counts all regions, including split parts
% vector of small regions to keep that are actually river segments
keep=logical(zeros(1,length(stats))); 
while Switch
    % filtering out small regions (prob shadows)
    % don't proceed until size condition is met
%     if stats(i).Area>MinSize
%         break
%     end
    if i>=length(stats) | i>=length(keep) % condition to end loop
        break
    end
    while stats(i).Area<MinSize & ~keep(i) % condition to filter out smalls
        i=i+1;
%         L(L==i)=0;
%         disp(i)
        if i>=length(stats) | i>=length(keep) % condition to end loop
            break
        end
    end
    if mod(j, 25)==0
        disp(['Frame ', num2str(j), '     ', datestr(now)])
        
    end
    if i>=length(stats)
            break
        end
    if PerSin(i)<0.09 && stats(i).Area > 4e5 %River condition
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        fprintf('region %d is probably a river.  Splitting...\n', i)
        per=stats(i).Perimeter;
        ar=stats(i).Area;
        box=stats(i).BoundingBox;
        buf=1*ar/per;
        box(1)=max(box(1)-buf, 1); box(2)=max(box(2)-buf, 1);
        box(3)=min(box(3)+2*buf, b-box(1)); box(4)=min(box(4)+2*buf, a-box(2));
        box=round(box);
        tic
        roi=L(box(2):box(2)+box(4),box(1):box(1)+box(3)); % subset of image for splitting
        Lp_roi=splitregion_5_2(roi, i, Len, Num); %good Len:Num ratio =40
        Lp=zeros(a,b); Lp(box(2):box(2)+box(4), box(1):box(1)+box(3))=Lp_roi;
        if ~isequal(size(L), size(Lp))
            error('EK: Caution: Size of L is not equal to size Lp')
        end
        if ~isequal(size(L(L==i)), size(L(Lp>0)))
            warning('EK: Caution: Regions ins Lp are not same size as regions in L==i')
            size(L(L==i)) %bigger.  why?
            size(L(Lp>0))
        end
        offset=max(max(Lp));
        L(L>i) = L(L>i) + offset-1;
        L(L==i) = L(L==i)+Lp(L==i)-1; % I fixed this line
        keep=[keep, zeros(1,offset-1)];
        keep(i:i+offset-1)=1;
        clear stats
        disp('    intermediate regionprops...')
        stats = regionprops(L,'Centroid', 'Area', 'BoundingBox',...
             'Perimeter');
        CircPerim=2*pi*sqrt([stats.Area]/pi) ;
        PerSin=CircPerim./[stats.Perimeter];
        disp('    intermediate regionprops finished.')
        toc
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    end
    per(j)=stats(i).Perimeter+1;
    ar(j)=stats(i).Area;
%         display(['Area of initial guess: ', num2str(stats(i).Area)])
%     disp('Structuring element...')
%         buf(j)=(-rm-rM+sqrt((rm+rM)^2+4*rm*rM))/2;
    buf(j)=ar(j)/per(j);
%         if ~isreal(buf)
%             disp('IMAGINARY BUFFER LENGTH.  EK')
%         end
    SE = strel('diamond',round(buf(j)));
    bw_dil_lakeN=imdilate(L==i, SE);
    [level(j), EM(j)] = graythresh((cir_ndwi(bw_dil_lakeN==1)));
%             if EM
    labeled=labeled | (bw_dil_lakeN.*(cir_ndwi>=level(j)));

    bar(j)=sum(sum(bw_dil_lakeN));
    ratios(j)=stats(i).Area/bar(j);
    PerSinrec(j)=PerSin(i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Plotting
%     subplot(221)
%     imshow(labeled)  
%     title(['Region ', num2str(i)])
%     box=stats(i).BoundingBox;
%     [rm, c]=find(L==i); %select region (Lake C) and label it
%     hold on
%     plot(c,rm, 'r.', 'MarkerSize', 15)
%     hold off
%     subplot(222)
%     imshow(cir)
%     axis([box(1), box(1)+box(3), box(2), box(2)+box(4)])
%     subplot(2,2,[3,4])
%     histogram((cir_ndwi(bw_dil_lakeN==1)))
%     hax=get(gca, 'YLim');
%     line([level(j), level(j)], hax)
%     G(j)=getframe(gcf);
%     pause(.3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    indx(j)=i; %record region number
    j=j+1;
    if i>=length(stats)
        Switch=false;
    end
    i=i+1;

        
end

%% Final filter on region size
L=bwlabel(labeled);
stats=regionprops('table', L, 'Area', 'Perimeter');
allowableAreaIndexes = (stats.Area > MinSize);
keeperIndexes = find(allowableAreaIndexes); 
newlabeled = ismember(L, keeperIndexes); 
clear L
stats=regionprops('table', newlabeled, 'Area', 'Perimeter', 'Centroid');


% fclose(FileID)
% CircPerim=2*pi*avg([[stats(1:j-1).MajorAxisLength]',...
%         [stats(1:j-1).MinorAxisLength]'])';
% PerSin=CircPerim./[stats(1:j-1).Perimeter];

% if ~exist('indx')
%     [indx, ar, bar, per, buf, ratios, EM, level,...
%         PerSinrec]=deal(1,1,1,1,1,1,1,1,1);
% end

%% Record stats and Extra stats
global COUNTER
if isempty(COUNTER)
    COUNTER=-9999;
end

if exist('indx')
    tabl=table([1:j-1]', indx', ar', bar', per', buf', ratios', EM', level',...
        PerSinrec', 'VariableNames',...
        {'Lake_number', 'Region_Number', 'Area', 'Buffered_Area', 'Perimeter', 'Buffer_Size',...
        'Area_ratio', 'EM', 'Thresh','PerSinuosity'})
    writetable(tabl, ['D:\CodeData\A_P_', date,'-',num2str(COUNTER), '.csv'], 'WriteVariableNames', 1)
    writetable(stats, ['D:\CodeData\FinalStats_', date,'-',num2str(COUNTER), '.csv'], 'WriteVariableNames', 1)
    disp('Log tables Written:')
    disp(['     D:\CodeData\A_P_', date,'-',num2str(COUNTER), '.csv'])
    disp(['     D:\CodeData\FinalStats_', date,'-',num2str(COUNTER), '.csv'])
%     fprintf('Average ratio of region to buffered region is %.2f\n',...
%         avg(ratios))
end
if Black
    labeled=zeros(a,b);
end
labeled=logical(labeled);
% clear ar bar per buf ratios EM PerSin indx

% imagesc(labeled)
%% Save classifcation and video

% Video

% viddir='D:\AboveData\pics\vid\';
% disp(['Saving to directory: ', viddir])
% v= VideoWriter([viddir, 'localthreshclassification5.mp4'], 'MPEG-4');
% v.FrameRate=2;
% open(v)
% writeVideo(v, G);
% close(v)

% Image

% gtinfo=geotiffinfo(img_in);
% disp('Saving to directory: D:\ArcGIS\FromMatlab\')
% geotiffwrite('D:\ArcGIS\FromMatlab\20170816 Ortho Segment_3-Pond25-26-clip_localThresh.tif',...
%     labeled,...
%     cir_R, 'GeoKeyDirectoryTag',gtinfo.GeoTIFFTags.GeoKeyDirectoryTag);
toc