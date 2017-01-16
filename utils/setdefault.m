function [ opts ] = setdefault( opts,deft )
%SETDEFAULT set default values to options by means of traversing the fields
%of the default structure
%   SETDEFAULT(opts,deft) sets the default opts with all the fields in
%   deft.
names = fieldnames(deft);
for i = 1:length(names)
    if ~isfield(opts,names{i})
        opts.(names{i}) = deft.(names{i}); % [default] initial value of penalty factor mu
    end
end

