function main_protein_analysis()
    % Main function to orchestrate a protein sequence analysis workflow.
    % Provides functionalities for sequence filtering, similarity search (local alignment),
    % motif discovery, and phylogenetic analysis.
    %
    % The workflow guides the user from data loading through various analysis options.
    % All user interactions are designed to be clear and bilingual.

    % --- 1. Load Primary FASTA Dataset ---
    [filename, pathname] = uigetfile(...
        {'*.fasta;*.fa;*.fna;*.faa', 'FASTA Files (*.fasta, *.fa, *.fna, *.faa)'; '*.*', 'All Files (*.*)'}, ...
        'Select your source FASTA file (选择您的原始FASTA文件)');
    if isequal(filename, 0)
        disp('Operation cancelled by user (用户取消操作).');
        return;
    end
    fasta_filepath = fullfile(pathname, filename);
    disp(['Reading file (正在读取文件): ', fasta_filepath]);

    try
        original_sequences = fastaread(fasta_filepath);
    catch ME
        errordlg(['Failed to read FASTA file (读取FASTA文件失败): ', ME.message], 'File Read Error (文件读取错误)');
        return;
    end

    if isempty(original_sequences)
        warndlg('No sequences found in the FASTA file (FASTA文件中没有序列).', 'Empty File (空文件)');
        return;
    end

    % --- 2. Initial Data Characterization ---
    original_lengths = arrayfun(@(s) length(s.Sequence), original_sequences);
    fprintf('\n--- Initial Dataset Summary (原始序列统计) ---\n');
    fprintf('Total number of sequences (总序列数): %d\n', length(original_lengths));
    fprintf('Min length (最小长度): %d\n', min(original_lengths));
    fprintf('Max length (最大长度): %d\n', max(original_lengths));
    fprintf('Mean length (平均长度): %.2f\n', mean(original_lengths));
    fprintf('Median length (中位长度): %d\n', median(original_lengths));

    figure('Name', 'Length Distribution of Original Sequences (原始序列长度分布)');
    histogram(original_lengths, 'BinMethod', 'auto');
    title('Length Distribution of Original Protein Sequences');
    xlabel('Sequence Length (amino acids)');
    ylabel('Frequency (序列数量)');
    grid on;

    % --- 3. Optional Sequence Filtering by Length ---
    sequences_to_analyze = original_sequences; % Default to original set
    
    filter_choice = questdlg('Would you like to filter sequences by length? (是否按长度筛选序列?)', ...
        'Sequence Filtering Option (序列筛选选项)', 'Yes (是)', 'No (否)', 'No (否)');

    if strcmp(filter_choice, 'Yes (是)')
        prompt = {'Enter minimum sequence length (输入最小序列长度):', ...
                  'Enter maximum sequence length (输入最大序列长度):', ...
                  'Enter output filename for filtered sequences (输入筛选后文件名):'};
        dlg_title = 'Filtering Parameters (筛选参数)';
        default_ans = {'100', '1000', 'filtered_proteins.fasta'};
        answer = inputdlg(prompt, dlg_title, [1, 60], default_ans);

        if isempty(answer)
            disp('Filtering cancelled. Proceeding with original dataset. (用户取消筛选，使用原始序列分析)');
        else
            min_len = str2double(answer{1});
            max_len = str2double(answer{2});
            output_filename = answer{3};

            if isnan(min_len) || isnan(max_len) || min_len < 0 || min_len > max_len
                warndlg('Invalid length range. Proceeding with original dataset. (无效长度范围，使用原始序列分析)', 'Input Error (输入错误)');
            else
                fprintf('Filtering sequences between %d and %d aa...\n', min_len, max_len);
                [filtered_seqs, filtered_lengths] = filter_and_save_proteins(original_sequences, output_filename, min_len, max_len);

                if ~isempty(filtered_seqs)
                    sequences_to_analyze = filtered_seqs; % Update the working set of sequences
                    fprintf('\n--- Filtered Dataset Summary (筛选后序列统计) ---\n');
                    fprintf('Number of sequences matching criteria (符合条件的序列数): %d\n', length(filtered_lengths));
                    fprintf('Filtered sequences saved to (已保存到): %s\n', output_filename);

                    figure('Name', 'Length Distribution of Filtered Sequences (筛选后序列长度分布)');
                    histogram(filtered_lengths, 'BinMethod', 'auto');
                    title('Length Distribution of Filtered Protein Sequences');
                    xlabel('Sequence Length (amino acids)');
                    ylabel('Frequency (序列数量)');
                    grid on;
                else
                    disp('No sequences met the filtering criteria. Proceeding with original dataset. (无序列符合条件，使用原始序列分析)');
                end
            end
        end
    end

    % --- 4. Main Analysis Menu ---
    while true
        choice_options = {
            'Similarity Search (序列相似性查找)', ...
            'Motif Search (基序查找)', ...
            'Build Phylogenetic Tree (构建进化树)', ...
            'Exit (退出程序)'
        };
        
        analysis_choice_idx = menu('Select an Analysis to Perform (请选择要执行的分析):', choice_options);

        if analysis_choice_idx == 0 || analysis_choice_idx == length(choice_options) % User closed menu or chose Exit
            disp('Exiting program (程序已退出).');
            close all; % Close all figures
            return;
        end
        
        analysis_choice = choice_options{analysis_choice_idx};

        switch analysis_choice
            case 'Similarity Search (序列相似性查找)'
                similarity_search_wrapper(sequences_to_analyze);
            case 'Motif Search (基序查找)'
                motif_search_wrapper(sequences_to_analyze);
            case 'Build Phylogenetic Tree (构建进化树)'
                phylogenetic_tree_wrapper(sequences_to_analyze);
            otherwise
                disp('Invalid selection (无效选择).');
        end
        fprintf('\n--- Analysis complete. Please select your next action. (当前分析完成，请选择下一步操作) ---\n\n');
    end
end

%% --- WRAPPER AND HELPER FUNCTIONS ---

function similarity_search_wrapper(target_sequences)
    % Encapsulates the workflow for finding proteins similar to a query sequence.
    if isempty(target_sequences)
        warndlg('Target sequence set is empty. Cannot perform search. (目标序列集为空)', 'Empty Set');
        return;
    end

    % Get query sequence
    [q_filename, q_pathname] = uigetfile('*.fasta;*.fa', 'Select Query FASTA File (选择查询FASTA文件)');
    if isequal(q_filename, 0), disp('Operation cancelled.'), return; end
    
    try
        query_struct = fastaread(fullfile(q_pathname, q_filename));
        if isempty(query_struct), warndlg('Query file is empty. (查询文件为空)', 'File Error'), return; end
        if length(query_struct) > 1, warning('Query file contains multiple sequences. Using the first one. (查询文件包含多条序列，仅使用第一条)'); end
        query_data = query_struct(1);
    catch ME
        errordlg(['Failed to read query file: ', ME.message], 'File Read Error');
        return;
    end

    % Get parameters
    answer = inputdlg({'Enter number of top hits to save (输入保存的序列数量 N):'}, 'Top N Hits', [1 40], {'20'});
    if isempty(answer), disp('Operation cancelled.'), return; end
    top_n = str2double(answer{1});
    if isnan(top_n) || top_n <= 0, warndlg('Invalid number for Top N.', 'Input Error'), return; end

    % Define output filename
    clean_header = regexprep(query_data.Header, '[^a-zA-Z0-9_-]', '_');
    default_out_name = sprintf('top_%d_similar_to_%s.fasta', top_n, clean_header);
    out_answer = inputdlg({'Enter output filename:'}, 'Save Similar Sequences', [1 70], {default_out_name});
    if isempty(out_answer), disp('Operation cancelled.'), return; end
    output_filename = out_answer{1};

    % Run the analysis
    fprintf('Searching for sequences similar to "%s"...\n', query_data.Header);
    sorted_results_table = find_similar_proteins(query_data.Sequence, target_sequences, output_filename, top_n);

    if isempty(sorted_results_table)
        disp('No significant similarity found.');
        return;
    end
    
    % Optional visualization
    vis_choice = questdlg('Visualize similarity results (identity vs. length)? (是否可视化相似性结果?)', ...
        'Visualization Option', 'Yes (是)', 'No (否)', 'Yes (是)');
    if strcmp(vis_choice, 'Yes (是)')
        visualize_similarity(query_data.Header, sorted_results_table);
    end
end

function motif_search_wrapper(target_sequences)
    % Encapsulates the workflow for finding sequences with a specific motif.
    if isempty(target_sequences), warndlg('Target sequence set is empty.', 'Empty Set'), return; end

    prompt = {'Enter Motif Pattern (e.g., P-X(2)-G or P..G): (输入基序模式)'};
    answer = inputdlg(prompt, 'Motif Search (基序查找)', [1 60], {''});
    if isempty(answer) || isempty(answer{1}), disp('Operation cancelled or motif is empty.'), return; end
    motif_pattern = answer{1};

    found_sequences = find_sequences_by_motif(target_sequences, motif_pattern);

    if ~isempty(found_sequences)
        save_choice = questdlg('Save the found sequences to a new FASTA file? (是否保存找到的序列?)', ...
            'Save Motif Results', 'Yes (是)', 'No (否)', 'Yes (是)');
        if strcmp(save_choice, 'Yes (是)')
            out_answer = inputdlg({'Enter output filename:'}, 'Save Motif Results', [1 50], {'motif_matches.fasta'});
            if ~isempty(out_answer)
                try
                    fastawrite(out_answer{1}, found_sequences);
                    fprintf('Sequences with motif saved to: %s\n', out_answer{1});
                catch ME
                    errordlg(['Failed to save file: ', ME.message], 'File Write Error');
                end
            end
        end
    end
end

function phylogenetic_tree_wrapper(sequences_to_analyze)
    % Encapsulates the workflow for building a phylogenetic tree.
    if isempty(sequences_to_analyze)
        warndlg('No sequences available to build a tree. (无可用序列)', 'Empty Set');
        return;
    end
    
    fprintf('You have %d sequences available for tree building.\n', length(sequences_to_analyze));

    % Let user select sequences for the tree
    headers = {sequences_to_analyze.Header};
    [indices, ok] = listdlg('PromptString', {'Select Sequences for Tree (按住Ctrl多选):', '(At least 3 sequences recommended)'}, ...
                            'SelectionMode', 'multiple', ...
                            'ListString', headers, ...
                            'Name', 'Select Sequences for Tree (选择序列建树)');
    
    if ~ok || isempty(indices)
        disp('Tree construction cancelled.');
        return;
    end
    
    if length(indices) < 3
        warndlg('At least 3 sequences are recommended for a meaningful tree. (建议至少选择3条序列)', 'Warning');
        if length(indices) < 2
            disp('At least 2 sequences are required. Aborting.');
            return;
        end
    end
    
    selected_sequences = sequences_to_analyze(indices);

    % Get bootstrap parameter
    answer = inputdlg({'Enter number of bootstrap replicates (e.g., 1000): (输入Bootstrap重复次数)'}, ...
                      'Bootstrap Parameter', [1 50], {'1000'});
    if isempty(answer), disp('Operation cancelled.'), return; end
    num_bootstraps = str2double(answer{1});
    if isnan(num_bootstraps) || num_bootstraps <= 0
        num_bootstraps = 1000; % Default value
        warning('Invalid bootstrap number. Defaulting to 1000. (无效输入，默认1000次)');
    end

    % Build the tree
    build_phylogenetic_tree(selected_sequences, num_bootstraps);
end