function [elf] = assert_elf_model(model)
% ASSERT_ELF_MODEL  Assert that a structure is an ELF model.
%
%   [ELF] = ASSERT_ELF_MODEL(MODEL)
%
%   Checks that the structure MODEL is an ELF model.  If not, converts
%   MODEL to an ELF model and warns that this procedure is not efficient
%   for repeated calls to the parent function.
%

fields = {'A','param'};

if ~all(isfield(model,fields))
    % convert model
    fprintf('This model is not an ELF model.  It will be automatically\n');
    fprintf('converted.  For repeated calls to this function, it is\n');
    fprintf('more efficient to convert beforehand.\n');
    
    elf = make_elf_model(model);
else
    elf = model;
end
