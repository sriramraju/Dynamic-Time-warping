% Modified from original script by Quan Wang (2013)
% This version returns warp path P and corrects norm to sum of squares

% dynamic time warping of two signals s and t with windowing constraint w

function [d, P] = dtw(s,t,w)
% s: signal 1, size is ns*k, row for time, column for channel
% t: signal 2, size is nt*k, row for time, column for channel
% w: window parameter
%      if s(i) is matched with t(j) then |i-j|<=w
% d: resulting distance
% P: resulting warp path matrix

if nargin<3
    w=Inf;
end

ns=size(s,1);
nt=size(t,1);
if size(s,2)~=size(t,2)
    error('Error in dtw(): the dimensions of the two input signals do not match.');
end
w=max(w, abs(ns-nt)); % adapt window size

% initialization

D=zeros(ns+1,nt+1)+Inf; % cache matrix
D(1,1)=0;

% begin dynamic programming

for i=1:ns
    for j=max(i-w,1):min(i+w,nt)
        %oost=norm(s(i,:)-t(j,:));
        %oost=sumsqr(s(i,:)-t(j,:));
        vec = s(i,:)-t(j,:);
        oost=vec*vec';
        D(i+1,j+1)=oost+min([D(i,j+1), D(i+1,j), D(i,j)]);     
    end
end

% recover warp path if requested

if nargout > 1

    i = ns;
    j = nt;

    P = zeros(ns,nt);

    while i>1 && j>1

       P(i,j) = 1;
       [~,m] = min([D(i,j), D(i+1,j), D(i,j+1)]);

       switch m
           case 1
               i = i-1;
               j = j-1;
           case 2
               j = j-1;
           case 3
               i = i-1;
       end
    end

    % Deal with edges

    while i>1
        P(i,j) = 1;
        i = i-1;   
    end

    while j>1
        P(i,j) = 1;
        j = j-1;
    end

    % Add the first point

    P(1,1) = 1;
    
end

% Return the distance

d=D(ns+1,nt+1);

end
