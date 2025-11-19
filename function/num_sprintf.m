function table_output = num_sprintf(table_input, format)
%NUM_SPRINTF Convert table with doubles to sprintf formats
%   table_input: input table
%   format: string format passed to the format function
    arguments
        table_input (:, :) table
        format (:, :) cell = {'%0.4f'}
    end
    n_var = size(table_input, 2);
    name_var = table_input.Properties.VariableNames;
    if (numel(format) ~= n_var) && (numel(format) ~= 1)
        error('The length of format must be one or the number of table variables')
    end
    for i_var = 1:n_var
        if isequal(class(table_input{:, i_var}), 'double')
            if numel(format) == 1
                format_i = format{1};
            else
                format_i = format{i_var};
            end
            table_input.(name_var{i_var}) = cellfun(@(x) sprintf(format_i, x), table2cell(table_input(:, i_var)), 'UniformOutput', false);
        end
    end
    table_output = table_input;
end