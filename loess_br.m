function [loess_br_fit] = loess_br(x,y,NF,cutoff_slope)
%[loess_br] = loess_br_trough(x,y,NF,cutoff_slope)
%
% x is a vector with regularly spaced x-coordinates and should be positive in
% upstream direction. 
%
% y is a vector containing the z-coordinates, corresponding to x. z is positive in upward
% direction. 
%
% NF is the correlation scale. A higher values of NF results in a smoother LOESS
% curve. 
%
% cutoff_slope is the slope of the primary lee side slope at which a break
% needs to be introduced. 
%
%
%  J Zomer, 2021-01-28

ydiff_cutoff = cutoff_slope; %default = 0.03

% Parameters set by user. 
window = 15; %m; 
distance=0.005; % For minima to 'identify as peak', the vertical distance to 
       % neighbouring extremes (to trough to neighbouring peaks) should be larger than 0.5cm. 
distance_v = 0.05; % Vertical distance between minimum and next 'trough'. 
window_trough = 1; % m
  
%%%% STEP 1: Unmodified LOESS curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lowpass = loess_interp(x,y,0,x,NF,2);

%%%% STEP 2: Approximate break locations based on unmodified LOESS curve %%
[troughs(:,2),troughs(:,1)] = findpeaks(-lowpass,x);troughs(:,2)=-troughs(:,2);
[crests(:,2),crests(:,1)] = findpeaks(lowpass,x);

breakpoints=[];
ydiff = diff(lowpass)./diff(x); % Local slope of LOESS curve

% Find leesides
for i=1:size(troughs(:,1))
 lss = find(crests(:,1) > troughs(i,1)); 
 if ~isempty(lss) 
 lss = lss(1); % Find next trough
 ls = find((x >= troughs(i,1)) & (x <= crests(lss,1))); % Select section between trough and next crest
    if max(ydiff(ls))>ydiff_cutoff  % If max lee slope is higher than ydiff_cutoff, breakpoint at trough location
     breakpoints(i)=troughs(i,1);
    end
 end
end

% If there are no breakpoints, apply LOESS, else, continue
if isempty(breakpoints) 
    lpassc_2 = loess_interp(x,y,0,x,NF,2);
else breakpoints(breakpoints==0)=[]; 
    
%%%% STEP 3: Update break location based on bed elevation profiles %%%%%%%%
% Find minimum y (unfiltered data) within specified distance ('window') from 
% breakpoint
    for i=1:length(breakpoints)
      idx = (x>breakpoints(i)-window & x<breakpoints(i)+window);
      xsel = x(idx); ysel=y(idx);
      idx = ysel == min(ysel);
      breakpoints(i) = xsel(idx); % breakpoints is updated
    end
    
%%%% STEP 4: Update break location based on local minima within window %%%%
    % Find all peaks
    [peaks, locs] = findpeaks(-y,x);peaks=-peaks;
    [tops, locs_tops] = findpeaks(y,x);
    [locs_all,ind_sort] = sort([locs; locs_tops]);
    peaks_all = [peaks; tops]; peaks_all = peaks_all(ind_sort);
    
    % All peaks are removed that have a vertical distance to 
    % up and downstreams crest < 5mm. 
    for n=1:size(peaks,1)  
        idx = find(locs_all==locs(n));
        if  ~(idx==1) && abs(peaks(n)-peaks_all(idx-1))<distance   
            peaks(n)=nan;
        elseif idx<size(peaks_all,1)
            if abs(peaks(n)-peaks_all(idx+1))<distance
            peaks(n)=nan; 
            end
        end
    end
    locs(isnan(peaks))=[];
    peaks(isnan(peaks))=[];

    % Find peaks within window    
    for i=1:length(breakpoints)
    idx=find(locs>breakpoints(i) & locs < breakpoints(i)+window);         
       locs_sel=locs(idx);
       peaks_sel=peaks(idx);

    idxx = find(abs(peaks_sel-y(x==breakpoints(i)))<distance_v);         % Absolute difference should be less than 5 cm. 
        if ~isempty(idxx);breakpoints(i) = max(locs_sel(idxx));end       % Take most upstream location of peaks. 
    clear xx yy
    end
            
%%%% STEP 5: Move breakpoint upstream if slope y is very low %%%%%%%%%%%%%%
        for i=1:length(breakpoints)
            idx = (x >= breakpoints(i) & x < breakpoints(i)+window_trough);
            xx = x(idx); yy=y(idx);
            for id=2:length(xx)
                if (yy(id)-yy(id-1)) < 0.005 % 
                    breakpoints(i) = xx(id);
                else; break
                end
            end
        end
    clear xx yy
    
%%%% STEP 6:Break up the transect based on x-values in breakpoints and 
%%%% apply LOESS stepwise
    breakpointsc = [min(x), breakpoints,max(x)]; % Add 0 and max x to breakpoints
    xxc=[];yyc=[];lpassc=[];
    for i=1:length(breakpointsc)-1
      xx{i} = x(x >= breakpointsc(i) & x <= breakpointsc(i+1));
      yy{i} = y(x >= breakpointsc(i) & x <= breakpointsc(i+1));
      [idx,~]=find(xx{i}==breakpointsc(2:end-1));
      xxl{i} = sort([xx{i};repmat(xx{i}(idx),50,1)]);
      yy{i} = interp1(xx{i}, yy{i}, xxl{i}); xx{i}=xxl{i};
      lpass{i} = loess_interp(xx{i},yy{i},0,xx{i},NF,2);
      [~,idx] = unique(xx{i});idx=[find(xx{i}==0);idx];
      xx{i} = xx{i}(idx);
      yy{i} = yy{i}(idx);
      lpass{i} = lpass{i}(idx);
      xxc = [xxc;xx{i}];
      yyc = [yyc;yy{i}];
      lpassc = [lpassc;lpass{i}];
    end
    lpassc_2 = loess_interp(xxc,lpassc,0,x,NF/10,2); % Apply another loess to get smooth sharp transitions at trough
end

loess_br_fit=lpassc_2;
end


