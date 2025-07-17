function build_phylogenetic_tree(sequences, num_bootstraps)
    % Constructs a phylogenetic tree from a set of protein sequences.
    % Workflow:
    % 1. Multiple Sequence Alignment (MSA) using MUSCLE or ClustalW algorithm.
    % 2. Calculate pairwise evolutionary distances using the JTT model.
    % 3. Construct a tree using the Neighbor-Joining (NJ) method.
    % 4. Perform bootstrap analysis to assess node confidence.
    % 5. Visualize the final tree with bootstrap support values.

    if nargin < 2
        num_bootstraps = 100; % Default bootstrap replicates
    end
    if isempty(sequences) || length(sequences) < 2
        errordlg('At least two sequences are required to build a tree. (至少需要两条序列)', 'Input Error');
        return;
    end

    headers = {sequences.Header};
    seqs_cell = {sequences.Sequence};

    % --- 1. Multiple Sequence Alignment (MSA) ---
    disp('Performing Multiple Sequence Alignment (MSA)... This may take some time.');
    try
        % multialign uses a ClustalW-like progressive alignment algorithm
        msa = multialign(seqs_cell, 'ScoringMatrix', 'blosum62', 'GapOpen', 10, 'GapExtend', 0.2);
    catch ME
        errordlg(['Multiple Sequence Alignment failed: ', ME.message], 'MSA Error');
        return;
    end
    disp('MSA complete.');
    
    % --- 2. Calculate Pairwise Distance Matrix ---
    disp('Calculating pairwise distance matrix using JTT model...');
    try
        % Use 'JTT' (Jones-Taylor-Thornton), a standard model for protein evolution.
        dist_matrix = seqpdist(msa, 'Method', 'JTT', 'Alpha', 'pairwise');
    catch ME
        errordlg(['Failed to calculate distance matrix: ', ME.message], 'Distance Matrix Error');
        return;
    end
    disp('Distance matrix calculation complete.');
    
    % --- 3. Construct Neighbor-Joining Tree ---
    disp('Constructing Neighbor-Joining tree...');
    try
        % *** CRITICAL FIX ***
        % Correctly specify the method ('NJ') and pass names as a parameter pair.
        tree = seqneighjoin(dist_matrix, 'Method', 'NJ', 'Names', headers);
    catch ME
        errordlg(['Failed to build phylogenetic tree: ', ME.message], 'Tree Construction Error');
        return;
    end
    disp('Initial tree construction complete.');

    % --- 4. Bootstrap Analysis ---
    if num_bootstraps > 0
        disp(['Performing bootstrap analysis with ', num2str(num_bootstraps), ' replicates... This is computationally intensive.']);
        try
            % Generate bootstrap replicates by resampling MSA columns
            boot_msa_samples = seqbootstrp(msa, 'NBoot', num_bootstraps);
            
            % Create a function handle for parallel processing
            dist_fun = @(s) seqpdist(s, 'Method', 'JTT', 'Alpha', 'pairwise');
            
            % Calculate distances for all bootstrap samples
            boot_dists = cellfun(dist_fun, boot_msa_samples, 'UniformOutput', false);

            % Build a tree for each bootstrap sample
            boot_trees = cell(num_bootstraps, 1);
            for k = 1:num_bootstraps
                boot_trees{k} = seqneighjoin(boot_dists{k}, 'Method', 'NJ', 'Names', headers);
            end
            
            % Get bootstrap support values and attach them to the original tree
            bootstrap_values = getbybranches(tree, boot_trees);
            
            % Format values as percentages for labels
            branch_labels = arrayfun(@(x) sprintf('%.0f', x*100), bootstrap_values, 'UniformOutput', false);
            tree.BranchNames = branch_labels;

            disp('Bootstrap analysis complete. Support values added to branches.');
        catch ME
            % *** CORRECTED LINE ***
            % The warning message is now formatted correctly using '%s' to prevent misinterpretation of ME.message content.
            warning('Bootstrap analysis failed: %s. Displaying tree without support values.', ME.message);
        end
    end

    % --- 5. Visualize the Tree ---
    disp('Visualizing phylogenetic tree...');
    h_tree = plot(tree, 'Orient', 'left');
    set(h_tree.BranchLabels, 'FontSize', 9, 'Color', [0.1 0.1 0.8]);
    set(h_tree.LeafLabels, 'FontSize', 10);
    title(sprintf('Phylogenetic Tree (Neighbor-Joining, JTT model, %d bootstraps)', num_bootstraps));
    xlabel('Evolutionary Distance');
end