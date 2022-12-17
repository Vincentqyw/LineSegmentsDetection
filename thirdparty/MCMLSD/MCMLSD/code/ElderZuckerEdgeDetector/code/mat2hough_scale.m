function [ r, th] = mat2hough_scale(m, n, r_max, r_res, th_res)
% MAT2HOUGH converts matrix coordinates to Hough map values
%
%   Inputs:
%               r_max:      Maximum value r can take
%               m:          vertical matrix index
%               n:          horizontal matrix index
%
%   Output:
%               r:          Hough domain r value equivalent to n
%               th:         Hough domain \theta value equivalent to m
%
%   Author:
%               Ron Tal
%
%   Date:
%               February 16, 2009
%
%   Description:
%               This function converts the coordinates of a matlab matrix
%               representing a Hough map, to the equivalent Hough domain
%               values
%
    th = round2frac((m - 1) * th_res, th_res);
    r = round2frac((n - 1) * r_res - r_max, r_res);
end