function [mdef_type, mdef_delta] = parse_mass_definition(mdef)
    if strcmpi(mdef, 'vir')
        mdef_type  = 'vir';
        mdef_delta = [];
    elseif ~isempty(regexp(mdef, '^\d+c$', 'once'))
        mdef_type  = 'c';
        mdef_delta = str2double(mdef(1:end-1));
    elseif ~isempty(regexp(mdef, '^\d+m$', 'once'))
        mdef_type  = 'm';
        mdef_delta = str2double(mdef(1:end-1));
    else
        error('parse_mass_definition: invalid mdef ''%s''.', mdef);
    end
end