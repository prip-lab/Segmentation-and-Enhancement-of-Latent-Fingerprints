function [response, enhImage ] = GaborRespose(img, blkOriImg,freq,blkSize)

border = 0;
[blkH,blkW] = size(blkOriImg);
sinOriImg = imresize(sin(2*blkOriImg), blkSize,'bilinear' );
cosOriImg = imresize(cos(2*blkOriImg), blkSize,'bilinear' );
oimg = atan2(sinOriImg,cosOriImg)*0.5;

% mask = imresize(blkMask, blkSize,'bilinear' );

[h,w]=size(img);


angleInc = 3;  % Fixed angle increment between filter orientations in
% degrees. This should divide evenly into 180


[rows, cols] = size(img);
enhImage = zeros(rows,cols);
response = zeros(rows,cols);
 freq = round(freq*100)/100;
 
 
[validr,validc] = find(freq > 0.03 & freq < 0.3);  % find where there is valid frequency data.
ind = sub2ind([blkH,blkW], validr, validc);

% Round the array of frequencies to the nearest 0.01 to reduce the
% number of distinct frequencies we have to deal with.
% freq(ind) = round(freq(ind)*100)/100;

% Generate an array of the distinct frequencies present in the array
% freq
unfreq = unique(freq(ind));

% Generate a table, given the frequency value multiplied by 100 to obtain
% an integer index, returns the index within the unfreq array that it
% corresponds to
freqindex = ones(100,1);
for k = 1:length(unfreq)
    freqindex(round(unfreq(k)*100)) = k;
end

% Generate filters corresponding to these distinct frequencies and
% orientations in 'angleInc' increments.
filter = cell(length(unfreq),180/angleInc);
sze = zeros(length(unfreq),1);


kx =0.4;
ky =0.4;
for k = 1:length(unfreq)
    sigmax = 1/unfreq(k)*kx;
    sigmay = 1/unfreq(k)*ky;
%     sgimax = 5;sigmay = 5;
%     sigmax =4;
%     sigmay = 4;
    sze(k) = round(3*max(sigmax,sigmay));
    [x,y] = meshgrid(-sze(k):sze(k));
    reffilter = exp(-(x.^2/sigmax^2 + y.^2/sigmay^2)/2)...
        .*(cos(2*pi*unfreq(k)*x) );%+ sqrt(-1) * sin(2*pi*unfreq(k)*x) );
    
    % Generate rotated versions of the filter.  Note orientation
    % image provides orientation *along* the ridges, hence +90
    % degrees, and imrotate requires angles +ve anticlockwise, hence
    % the minus sign.
    for o = 1:180/angleInc
        filter{k,o} = imrotate(reffilter,-(o*angleInc+90),'bilinear','crop');
        filter{k,o}  = filter{k,o}   - mean(filter{k,o} (:));
%         filter{k,o}  = filter{k,o} /(sum(sum(abs(filter{k,o}) )));
    end
end

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
freq = imresize(freq, blkSize,'nearest' );
% Finally do the filtering
for i=border+1: blkH*blkSize + border
    for j = border+1:blkW*blkSize + border
        % find filter corresponding to freq(r,c)
        r = round(i/blkSize)+1;  
        if r >blkH 
            r = blkH;
        end
        c = round(j/blkSize)+1;
          if c >blkW 
            c = blkW;
          end
          sfreq = round(freq(i,j)*100) ;
        if(sfreq== 0)
            continue;
        end
        filterindex = freqindex(sfreq);
        s = sze(filterindex);
        
        
        temp = zeros(2*s+1,2*s+1);
        y1 = i - s;  fy1=1; y2 = i+s;   fy2=2*s+1;
        x1 = j - s;  fx1=1; x2 = j + s; fx2=2*s+1;
        if(y1<1)
            fy1 = 1 + 1 - y1;
            y1=1;
        end
        
        if x1<1
            fx1 = 1 + 1 - x1;
            x1 = 1;
        end
        if y2>h
            fy2 = fy2+  h - y2;
            y2 = h;
        end
        if x2>w
            fx2 =fx2 + w - x2;
            x2 = w;
        end
        temp(fy1:fy2,fx1:fx2) = img(y1:y2, x1:x2);
        
        singlefilter = filter{filterindex,orientindex(i,j)};
        blk = real(temp.*singlefilter);
        enhImage(i,j) = sum(blk(:));
        response(i,j) =sum(abs(blk(:)));% abs(sum(blk(:)));
    end
end
    
    