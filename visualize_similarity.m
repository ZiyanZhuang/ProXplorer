function visualize_similarity(query_sequence_header, sorted_results_table)
    % Visualizes similarity search results as a scatter plot of
    % percent identity versus sequence length, with alignment score as color.

    if isempty(sorted_results_table)
        disp('No data available for visualization.');
        return;
    end

    lengths = cellfun(@length, sorted_results_table.Sequence);
    identities = sorted_results_table.Identity;
    scores = sorted_results_table.Score;
    headers = sorted_results_table.Header;

    figure('Name', ['Similarity Results for ', query_sequence_header]);
    scatter(lengths, identities, 50, scores, 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2]);
    
    % Use 'parula' colormap (perceptually uniform) instead of 'jet'
    colormap(parula);
    cb = colorbar;
    ylabel(cb, 'Alignment Score');

    title(['Sequences Similar to "', query_sequence_header, '"']);
    xlabel('Sequence Length (amino acids)');
    ylabel('Percent Identity (%)');
    grid on;
    set(gca, 'FontSize', 12);

    % Interactive data tips
    dcm = datacursormode(gcf);
    set(dcm, 'UpdateFcn', @(src, event) createDataTipText(event, headers));
    disp('Visualization created. Click on points to see sequence details.');
end

function output_txt = createDataTipText(event_obj, headers)
    % Custom data tip update function
    pos = get(event_obj, 'Position');
    idx = get(event_obj, 'DataIndex');
    output_txt = {
        sprintf('Header: %s', headers{idx}), ...
        sprintf('Length: %d', pos(1)), ...
        sprintf('Identity: %.2f%%', pos(2)), ...
        sprintf('Score: %.2f', event_obj.Target.CData(idx))
    };
end