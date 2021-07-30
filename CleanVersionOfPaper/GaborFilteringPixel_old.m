function [eimg] = GaborFilteringPixel(img,mask,fimg,oimg,opts)
% global gaborfilter sze

if isa(img,'uint8')
    img = double(img);
end

gaborfilter = opts.gaborfilter;
sze = opts.sze;
    
% load('GaborFilters2','gaborfilter','sze');
% load GaborFilters
mask = double(mask);
MaxCount = 3;

[h,w] = size(mask);


eimg = ones(size(mask));


if nargin<4   
    [gx,gy] = gradient(img);
    Gxx = gx.*gx;
    Gxy = gx.*gy;
    Gyy = gy.*gy;
    
    gh = fspecial('gaussian',10,5);
    Gxx = imfilter(Gxx,gh);
    Gxy = imfilter(Gxy,gh);
    Gyy = imfilter(Gyy,gh);
    
    oimg = 0.5*atan2(2*Gxy,Gxx - Gyy)+pi/2;
end

if size(oimg,1)~=h | size(oimg,2)~=w | size(fimg,1)~=h | size(fimg,2)~=w
    error('the mask, oimg and fimg should be of the same size.');
end
% % generate filters for ridge generation
% fimg = round(fimg*100)/100;
% 
% 
% [validr,validc] = find(fimg > 0.03 & fimg < 0.3);  % find where there is valid frequency data.
% ind = sub2ind([h,w], validr, validc);
% 
% % Round the array of frequencies to the nearest 0.01 to reduce the
% % number of distinct frequencies we have to deal with.
% % freq(ind) = round(freq(ind)*100)/100;
% 
% % Generate an array of the distinct frequencies present in the array
% % freq
% unfreq = unique(fimg(ind));
% 
% % Generate a table, given the frequency value multiplied by 100 to obtain
% % an integer index, returns the index within the unfreq array that it
% % corresponds to
% freqindex = ones(100,1);
% for k = 1:length(unfreq)
%     freqindex(round(unfreq(k)*100)) = k;
% end

% Generate filters corresponding to these distinct frequencies and
% orientations in 'angleInc' increments.
angleInc = 3;  % Fixed angle increment between filter orientations in
% degrees. This should divide evenly into 180


% filter = cell(length(unfreq),180/angleInc);
% sze = zeros(length(unfreq),1);


% kx =0.4;
% ky =0.4;
% for k = 1:length(unfreq)
%     sigmax = 1/unfreq(k)*kx;
%     sigmay = 1/unfreq(k)*ky;
%     %     sgimax = 5;sigmay = 5;
%     %     sigmax =4;
%     %     sigmay = 4;
%     sze(k) = round(3*max(sigmax,sigmay));
%     [x,y] = meshgrid(-sze(k):sze(k));
%     reffilter = exp(-(x.^2/sigmax^2 + y.^2/sigmay^2)/2)...
%         .*(cos(2*pi*unfreq(k)*x) );%+ sqrt(-1) * sin(2*pi*unfreq(k)*x) );
%     
%     % Generate rotated versions of the filter.  Note orientation
%     % image provides orientation *along* the ridges, hence +90
%     % degrees, and imrotate requires angles +ve anticlockwise, hence
%     % the minus sign.
%     for o = 1:180/angleInc
%         filter{k,o} = imrotate(reffilter,-(o*angleInc+90),'bilinear','crop');
%         filter{k,o}  = filter{k,o}   - mean(filter{k,o} (:));
%         %         filter{k,o}  = filter{k,o} /(sum(sum(abs(filter{k,o}) )));
%     end
% end




% if showfilter % Display largest scale filter for inspection
%     figure(7), imshow(filter{1,end},[]); title('filter');
% end

% Find indices of matrix points greater than maxsze from the image
% boundary
maxsze = sze(1);


% Convert orientation matrix values from radians to an index value
% that corresponds to round(degrees/angleInc)
maxorientindex = round(180/angleInc);
orientindex = round(oimg/pi*180/angleInc);
% keyboard
i = find(orientindex < 1);   orientindex(i) = orientindex(i)+maxorientindex;
i = find(orientindex > maxorientindex);
orientindex(i) = orientindex(i)-maxorientindex;
% Finally do the filtering

fimg(fimg>0.3) = 0.3;
fimg(fimg<0.05) = 0.05;
% fimg(is())
for i=1:h
    for j = 1:w
        if mask(i,j)==0
             continue;
         end
        % find filter corresponding to freq(r,c)
     
        sfreq = round(fimg(i,j)*100) ;
%         if(sfreq== 0)
%             continue;
%         end
%         filterindex = freqindex(sfreq);
%         sfreq
        s = sze(sfreq);
        
        if i<=s | j<=s | i>h-s | j>w-s
            continue;
        end
        
        temp = zeros(2*s+1,2*s+1);
        y1 = i - s;  fy1=1; y2 = i+s;   fy2=2*s+1;
        x1 = j - s;  fx1=1; x2 = j + s; fx2=2*s+1;
%         if(y1<1)
%             fy1 = 1 + 1 - y1;
%             y1=1;
%         end
%         
%         if x1<1
%             fx1 = 1 + 1 - x1;
%             x1 = 1;
%         end
%         if y2>h
%             fy2 = fy2+  h - y2;
%             y2 = h;
%         end
%         if x2>w
%             fx2 =fx2 + w - x2;
%             x2 = w;
%         end
        blk = img(y1:y2, x1:x2);
%         blk = blk-mean(blk(:));
%         blk = blk./std(blk(:));
        temp(fy1:fy2,fx1:fx2) = blk;%img(y1:y2, x1:x2);
        %             if find(isnan(temp))
        %                 keyboard
        %             end
        singlefilter = gaborfilter{sfreq,orientindex(i,j)};
        blk = real(temp.*singlefilter);
        eimg(i,j) = sum(blk(:));
        
    end
end
%     img = mat2gray(enhImg);
% %      keyboard
% end
% keyboard
