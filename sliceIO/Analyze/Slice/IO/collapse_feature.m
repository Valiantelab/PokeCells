function [all_features] = collapse_feature(F, features, norm_type)
% USAGE: [all_features] = collapse_feature(F, varargin)
%
% Takes  a feature cell array and generates a matrix of spike features
%
%Input:
%   F: feature cell array
%   varargin: the features that one wants:
%       
all_features = [];
for i=1:numel(F)
    af2 = [];
    for j= 1:numel(F{i})
        if ~isempty(F{i}{j})
            count = 0;
            af1 = [];
            for k=1:numel(F{i}{j})
                if are_fields(F{i}{j}{k}, features) 
                    count = count + 1;
                    for v=1:numel(features)
                        af1(count,v) = getfield(F{i}{j}{k},features{v});
                    end
                    af1(count,v+1) = k;
                else
                    error('A specific field does not exist');
                end
            end
            if norm_type == 1;
             
                % Normalize within each epoch
                if min(size(af1)) > 1
                    m =  max(af1);
                    m(numel(m)) = 1;
                    af1 = af1./repmat(m, size(af1,1), 1);
                else
                    af1 = af1./af1;
                end
            end
            af2 = [af2' af1']';
        end
    end
    if norm_type == 2;
        % Normalize over all epochs
        if min(size(af1)) > 1
            m =  max(af2);
            m(numel(m)) = 1;
            af2 = af2./repmat(m, size(af2,1), 1);
        end
    end
    all_features = [all_features' af2']';
end

function [okay] = are_fields(f, args)

okay = true;
for i=1:numel(args)
    okay = okay & isfield(f, args{1});
end