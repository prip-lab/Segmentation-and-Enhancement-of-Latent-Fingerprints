function [val, indr, indc] = findpeaks2D(X, minormax)
% find peaks
% 
% input: 
%          X           --- m x n data
%          minormax    --- 'min': find local min, or
%                          'max': find local max
% 
% output:
%          val         --- peaks val
%          indr        --- row ind of peaks
%          indc        --- colomn ind of peaks
% 
% Fanglin Chen, 2010-05-20

val = [];
indr = [];
indc = [];
if ndims(X) ~= 2
    disp('Error in findpeak2D: the dim of input X is not equal to 2!!');
    return;
end

if ~(strcmp(minormax, 'min') || strcmp(minormax, 'max'))
    disp('Error in findpeaks2D: input minormax should be ''min'' of ''max''');
    return;
end

[m, n] = size(X);
minX = min(min(X));
maxX = max(max(X));
Y1 = X;
if strcmp(minormax, 'max') == 1
    for ii = 2 : m - 1
        for jj = 2 : n - 1
            if Y1(ii, jj) < X(ii, jj+1) || Y1(ii, jj) < X(ii+1, jj-1) ...
                    || Y1(ii, jj) < X(ii+1, jj) || Y1(ii, jj) < X(ii+1, jj+1)
                Y1(ii, jj) = minX;
            end
        end
    end
    
    tmp = diff(Y1(1, :));
    Y1(1, tmp > 0) = minX;
    
    tmp = diff(Y1(m, :));
    Y1(m, tmp > 0) = minX;
    
    tmp = diff(Y1(:, 1));
    Y1(tmp > 0, 1) = minX;
    
    tmp = diff(Y1(:, n));
    Y1(tmp > 0, n) = minX;
elseif strcmp(minormax, 'min') == 1
    for ii = 2 : m - 1
        for jj = 2 : n - 1
            if Y1(ii, jj) > X(ii, jj+1) || Y1(ii, jj) > X(ii+1, jj-1) ...
                    || Y1(ii, jj) > X(ii+1, jj) || Y1(ii, jj) > X(ii+1, jj+1)
                Y1(ii, jj) = maxX;
            end
        end
    end
    
    tmp = diff(Y1(1, :));
    Y1(1, tmp < 0) = maxX;
    
    tmp = diff(Y1(m, :));
    Y1(m, tmp < 0) = maxX;
    
    tmp = diff(Y1(:, 1));
    Y1(tmp < 0, 1) = maxX;
    
    tmp = diff(Y1(:, n));
    Y1(tmp < 0, n) = maxX;
end
        
Y2 = Y1;
if strcmp(minormax, 'max') == 1
    for ii = m - 1 : -1 : 2
        for jj = n - 1 : -1 : 2
            if Y2(ii, jj) <= Y1(ii, jj-1) || Y2(ii, jj) <= Y1(ii-1, jj+1) ...
                    || Y2(ii, jj) <= Y1(ii-1, jj) || Y2(ii, jj) <= Y1(ii-1, jj-1)
                Y2(ii, jj) = minX;
            end
        end
    end
    
    tmp = diff(Y1(1, :));    
    Y2(1, find(tmp <= 0) + 1) = minX;
    
    tmp = diff(Y1(m, :));
    Y2(m, find(tmp <= 0) + 1) = minX;
    
    tmp = diff(Y1(:, 1));
    Y2(find(tmp <= 0) + 1, 1) = minX;
    
    tmp = diff(Y1(:, n));
    Y2(find(tmp <= 0) + 1, n) = minX;
    
    ind = find(Y2 ~= minX);
elseif strcmp(minormax, 'min') == 1
    for ii = m - 1 : -1 : 2
        for jj = n - 1 : -1 : 2
            if Y2(ii, jj) >= Y1(ii, jj-1) || Y2(ii, jj) >= Y1(ii-1, jj+1) ...
                    || Y2(ii, jj) >= Y1(ii-1, jj) || Y2(ii, jj) >= Y1(ii-1, jj-1)
                Y2(ii, jj) = maxX;
            end
        end
    end
    
    tmp = diff(Y1(1, :));
    Y2(1, find(tmp >= 0) + 1) = maxX;
    
    tmp = diff(Y1(m, :));
    Y2(m, find(tmp >= 0) + 1) = maxX;
    
    tmp = diff(Y1(:, 1));
    Y2(find(tmp >= 0) + 1, 1) = maxX;
    
    tmp = diff(Y1(:, n));
    Y2(find(tmp >= 0) + 1, n) = maxX;
    
    ind = find(Y2 ~= maxX);
end

if isempty(ind)
    disp('Warning in findpeaks2D: did not find!!');
    return;
end

val = X(ind);
[indr, indc] = ind2sub([m, n], ind);