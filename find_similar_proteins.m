function [sorted_results_table] = find_similar_proteins(query_sequence, target_sequences, output_filename, top_n)
    % Finds proteins similar to a query sequence using Smith-Waterman local alignment.
    % Results are ranked by percent identity and saved to a FASTA file.

    % --- Alignment Parameters ---
    SCORING_MATRIX = blosum62;
    GAP_OPEN_PENALTY = 10;
    GAP_EXTEND_PENALTY = 1;

    % --- Input Validation ---
    empty_table = table('Size', [0, 4], 'VariableTypes', {'string', 'string', 'double', 'double'}, ...
                        'VariableNames', {'Header', 'Sequence', 'Identity', 'Score'});
    if isempty(query_sequence), disp('Error: Query sequence is empty.'); sorted_results_table = empty_table; return; end
    if isempty(target_sequences), disp('Error: Target sequence set is empty.'); sorted_results_table = empty_table; return; end

    num_targets = length(target_sequences);
    results = cell(num_targets, 4); % Pre-allocate for performance

    fprintf('Performing local alignment against %d target sequences...\n', num_targets);
    fprintf('Parameters -> Matrix: BLOSUM62, Gap Open: %d, Gap Extend: %d\n', GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY);
    
    h_wait = waitbar(0, 'Calculating sequence similarities...');

    for i = 1:num_targets
        target_seq_header = target_sequences(i).Header;
        target_seq = target_sequences(i).Sequence;

        if isempty(target_seq), continue; end

        [score, alignment] = swalign(query_sequence, target_seq, ...
                                     'ScoringMatrix', SCORING_MATRIX, ...
                                     'GapOpen', GAP_OPEN_PENALTY, ...
                                     'ExtendGap', GAP_EXTEND_PENALTY);
        
        % --- Percent Identity Calculation ---
        % Definition: (Number of identical matches) / (Length of alignment excluding gaps).
        % This is a standard method that evaluates identity only over the aligned blocks.
        aligned_query = alignment(1,:);
        aligned_target = alignment(3,:);
        
        is_match = (aligned_query == aligned_target);
        is_not_gap_in_either = (aligned_query ~= '-') & (aligned_target ~= '-');
        
        num_identical = sum(is_match & is_not_gap_in_either);
        effective_length = sum(is_not_gap_in_either);

        identity_percent = 0;
        if effective_length > 0
            identity_percent = (num_identical / effective_length) * 100;
        end

        results(i,:) = {string(target_seq_header), string(target_seq), identity_percent, score};
        
        if mod(i, 20) == 0, waitbar(i/num_targets, h_wait); end
    end
    close(h_wait);

    results_table = cell2table(results, 'VariableNames', {'Header', 'Sequence', 'Identity', 'Score'});
    results_table = rmmissing(results_table); % Remove rows for any skipped empty sequences

    % Sort by Identity (primary) and Score (secondary)
    sorted_results_table = sortrows(results_table, {'Identity', 'Score'}, {'descend', 'descend'});

    num_to_display = min(top_n, height(sorted_results_table));
    if num_to_display == 0 || sorted_results_table.Identity(1) == 0
        disp('No sequences with significant similarity were found.');
        sorted_results_table = empty_table;
        return;
    end

    fprintf('\n--- Top %d Similar Sequences ---\n', num_to_display);
    for i = 1:num_to_display
        fprintf('%d. Header: %s | Identity: %.2f%% | Score: %.2f\n', ...
                i, sorted_results_table.Header(i), sorted_results_table.Identity(i), sorted_results_table.Score(i));
    end

    % Save top N sequences to a new FASTA file
    top_sequences = struct('Header', {}, 'Sequence', {});
    for i = 1:num_to_display
        top_sequences(i).Header = sprintf('%s | Identity: %.2f%% | Score: %.2f', ...
                                          sorted_results_table.Header(i), ...
                                          sorted_results_table.Identity(i), ...
                                          sorted_results_table.Score(i));
        top_sequences(i).Sequence = sorted_results_table.Sequence{i};
    end

    try
        fastawrite(output_filename, top_sequences);
        fprintf('\nTop %d similar sequences saved to: %s\n', num_to_display, output_filename);
    catch ME
        warning('Failed to write output file: %s', ME.message);
    end
end