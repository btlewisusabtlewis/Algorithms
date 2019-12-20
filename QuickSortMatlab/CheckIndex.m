function CheckIndex(l, r, index)
    % Verify that index is between indices l and r.
    if (index < l) || (index > r)
        error('** ERROR: index %i is not in [%i..%i]\n', index, l, r);
    end;
end

