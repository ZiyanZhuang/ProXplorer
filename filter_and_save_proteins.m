function [filtered_sequences, filtered_lengths] = filter_and_save_proteins(sequences, output_filename, min_len, max_len)
    % Filters protein sequences based on a specified length range and saves them to a new FASTA file.
    % For very large datasets, pre-allocation of 'filtered_sequences' could offer a minor performance gain.

    if ~isstruct(sequences) || ~all(isfield(sequences, {'Header', 'Sequence'}))
        error('Input must be a valid FASTA structure array with Header and Sequence fields.');
    end
    
    % Use logical indexing for efficient filtering
    lengths = arrayfun(@(s) length(s.Sequence), sequences);
    keep_indices = lengths >= min_len & lengths <= max_len;
    
    filtered_sequences = sequences(keep_indices);
    filtered_lengths = lengths(keep_indices);

    if ~isempty(filtered_sequences)
        try
            fastawrite(output_filename, filtered_sequences);
        catch ME
            warning('Failed to write to FASTA file (写入文件失败): %s', ME.message);
        end
    else
        disp('No sequences matched the specified length criteria; output file was not created.');
    end
end