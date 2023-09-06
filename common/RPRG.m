function [index, intersct] = RPRG(index, thre)
%
% Ridge path regrouping (RPRG) 
%
% -------------- Input ------------------
%  iniindex: initial time/frequency indexes of ridge curves, each curve indexes should be listed in one row
%  thre: threshold for finding the intersection intervals
%
% -------------- Output -----------------
%  index: final curve indexes after regrouping
%  intersect: detected intersetction intervals for ridge curves
%
% Author: Yuan JIANG
% Time: 2023-09-06

%% Finding all intersection intervals
[M, N] = size(index);    % M curves with the length of N
pairs = nchoosek(1: M, 2);     % curve pairs when randomly selecting 2 from M curves
Npair = size(pairs, 1);

num = 1;
intersct = {};
for i = 1: Npair
    tempindex = find( abs( index(pairs(i, 1), :) - index(pairs(i, 2), :) ) <= thre );
    if ~isempty(tempindex)
        tempindex2 = find(diff(tempindex) >= N/70); % find if there are multiple intersections for two overlapped ridge curves
        if ~isempty(tempindex2)
            for j = 1: length(tempindex2) + 1
                if j == 1
                    intersct(num, 1) = {[pairs(i,1), pairs(i,2)]};   % ridge numbers of two intersected curves
                    intersct(num, 2) = {[tempindex(1), tempindex(tempindex2(1))]};   % starting and ending time for intersected intervals
                    num = num + 1;
                else
                    if j == length(tempindex2) + 1
                        intersct(num, 1) = {[pairs(i,1), pairs(i,2)]};
                        intersct(num, 2) = {[tempindex(tempindex2(end)+1), tempindex(end)]};
                        num = num + 1;
                    else
                        intersct(num, 1) = {[pairs(i,1), pairs(i,2)]};
                        intersct(num, 2) = {[tempindex(tempindex2(j-1)+1), tempindex(tempindex2(j))]};
                        num = num + 1;
                    end
                end
            end
        else
            intersct(num, 1) = {[pairs(i,1), pairs(i,2)]};
            intersct(num, 2) = {[tempindex(1), tempindex(end)]};
            num = num + 1;
        end
    end
end

%% Merge close intersection intervals
if isempty(intersct)
    return;
end

[Ninter, ~] = size(intersct);
if Ninter > 2
    pairs = nchoosek(1: Ninter, 2); % intersection pairs when randomly selecting 2
    Npair = size(pairs, 1);
    
    for i = 1: Npair
        % interval 1
        comp1 = intersct{pairs(i,1), 1};
        interv1 = intersct{pairs(i,1), 2};
        locy1 = mean(index(comp1, interv1(1)));
        locx1 = interv1(1);
        % interval 2
        comp2 = intersct{pairs(i,2), 1};
        interv2 = intersct{pairs(i,2), 2};
        locy2 = mean(index(comp2, interv2(1)));
        locx2 = interv2(1);
        
        if ((locx1 - locx2)^2 + (locy1 - locy2)^2) <= 4*thre^2  % the distance of 2 intersection initervals
            ucomp = union(comp1, comp2);        % merge ridge numbers
            intersct(pairs(i,1), 1) = {ucomp};
            uset = union(interv1, interv2);     % merge time instants
            intersct(pairs(i,1), 2) = {[uset(1), uset(end)]};
            intersct(pairs(i,2), :) = intersct(pairs(i,1), :);  % assign merging results to two intersection intervals
        end
    end
    
    % delete the same intersection intervals
    douinter = zeros(Ninter, 100);
    for i = 1: Ninter
        ddd = cat(2, intersct{i, :});
        douinter(i, 1:length(ddd)) = ddd;
    end
    [~, m, ~] = unique(douinter, 'rows');   % merge rows with same elements
    dd = setdiff(1: Ninter, m);
    intersct(dd, :) = [];   % delete
    
    % sort intersection intervals in ascending order of ending time of each interval
    [Ninter, ~] = size(intersct);
    rightend = zeros(1, Ninter);
    for i = 1: Ninter
        rightend(i) = intersct{i,2}(2);
    end
    [~, sorin] = sort(rightend);
    for i = 1: Ninter
        sortset(i,:) = intersct(sorin(i), :);
    end
    intersct = sortset;
    
end

%% regroup ridge curves
len = floor(N/70);
for i = 1: Ninter
    comp = intersct{i, 1};
    Ncomp = length(comp);
    left = intersct{i,2}(1);
    right = intersct{i,2}(2);
    
    if (left - len) < 1
        slopr = (index(comp, right+len) - index(comp, right)) / len;
        for ii = 1: Ncomp
            index(comp(i), 1:right) = round(index(comp(i),right) + slopr(i)*((1:right)-right)); % linear prediction
        end
    end
    
    if (right + len) > N
        slopl = (index(comp, left) - index(comp, left-len)) / len;
        for ii = 1: Ncomp
            index(comp(i), left:N) = round(index(comp(i), left) + slopl(i)*((left:N)-left));    % linear prediction
        end
    end
    
    if (left-len) >= 1 && (right+len) <= N
        slopl = index(comp, left) - index(comp, left-len);  % slopes
        slopr = index(comp, right+len) - index(comp, right);    
        slopl = slopl(:);
        slopr = slopr(:);
        deltaslop = abs(bsxfun(@minus, slopr, slopl')); % connection matrix
        linearinterpol = zeros(Ncomp, right);
        for ii = 1: Ncomp
            [jj, kk] = find(deltaslop == min(deltaslop(:)));
            jj = jj(1);
            kk = kk(1);
            deltaslop(:, kk) = inf;
            deltaslop(jj, :) = inf;
            linearinterpol(jj, :) = [index(comp(kk), 1:(left-1)), round(linspace(index(comp(kk),left), index(comp(jj),right), right-left+1))];  % linear interpolation
        end
        for jj = 1: Ncomp
            index(comp(jj), 1:right) = linearinterpol(jj, :);
        end
    end
end
