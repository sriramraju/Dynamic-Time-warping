function DTW_dist = alignDTW(in_1,in_2)

% Author: Sriram Raju Dandu
% Alignment of cycles and obtaining DTW distance

% in_1 and in_2 are 3D accelerometer values of similiar length

MAGin_1 = sqrt(in_1(:,1).^2 + in_1(:,2).^2 + in_1(:,3).^2);
MAGin_2 = sqrt(in_2(:,1).^2 + in_2(:,2).^2 + in_2(:,3).^2);

[corrVals lags] = xcorr(MAGin_1,MAGin_2);
shiftPT = lags(find(corrVals == max(corrVals)));

shiftIN_2 = circshift(in_2,shiftPT);

[DTW_dist ~] = dtw(in_1,shiftIN_2);

