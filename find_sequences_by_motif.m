function [motif_found_sequences] = find_sequences_by_motif(sequences_to_search, motif_pattern_input)
    % Finds sequences containing a specific protein motif pattern.
    % Supports standard characters and ambiguity, e.g., 'X' for any residue.
    % Also supports PROSITE-like syntax such as 'X(2)' for repeats.

    if isempty(sequences_to_search), error('Target sequence set is empty.'); end
    if isempty(motif_pattern_input), error('Motif pattern cannot be empty.'); end

    fprintf('Searching for motif: "%s"...\n', motif_pattern_input);

    % Convert PROSITE-like syntax to a standard regular expression
    % e.g., 'P-X(2)-G' -> 'PXXG' -> 'P..G'
    % 1. Expand X(n) notation
    motif_pattern = regexprep(motif_pattern_input, 'X\((\d+)\)', '${repmat(''X'',1,str2double($1))}');
    % 2. Remove hyphens and convert 'X' to '.' (any character)
    regex_pattern = strrep(motif_pattern, '-', '');
    regex_pattern = strrep(regex_pattern, 'X', '.');

    fprintf('Converted to RegEx pattern: "%s"\n', regex_pattern);

    found_indices = false(1, length(sequences_to_search));
    for i = 1:length(sequences_to_search)
        if ~isempty(regexp(sequences_to_search(i).Sequence, regex_pattern, 'once'))
            found_indices(i) = true;
        end
    end
    
    motif_found_sequences = sequences_to_search(found_indices);

    if ~isempty(motif_found_sequences)
        fprintf('Successfully found %d sequence(s) containing the motif.\n', length(motif_found_sequences));
    else
        disp('No sequences were found containing the specified motif.');
    end
end