function [ options ] = fuseset( varargin )
%FUSESET Summary of this function goes here
%   Detailed explanation goes here

options = [];
for i = 1:2:length(varargin)
    options = setfield(options,varargin{i},varargin{i+1});
end

end

