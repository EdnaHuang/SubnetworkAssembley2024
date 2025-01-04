function mi_matrix = calculate_mi_matrix_trace(data, L)
%CALCULATE_MI_MATRIX Calculates the mutual information matrix for input data.
%https://www.mathworks.com/matlabcentral/fileexchange/13289-fast-mutual-information-of-two-images-or-signals
% Jose Delpiano (2024). Fast mutual information of two images or signals
% (https://www.mathworks.com/matlabcentral/fileexchange/13289-fast-mutual-information-of-two-images-or-signals)
% MATLAB Central File Exchange. Retrieved November 25, 2024.
% Edit the function to handel our data

%   mi_matrix = calculate_mi_matrix(data, L) calculates the mutual information
%   (MI) between all pairs of rows in the input matrix `data`. The MI values
%   are stored in a symmetric matrix.
%
%   Input:
%   - data: m x n matrix (m = rows, n = columns)
%   - L (optional): Number of bins for histogram (default: 256)
%
%   Output:
%   - mi_matrix: m x m matrix of mutual information values

    % Set default number of bins
    if nargin < 2
        L = 256; % Default number of bins
    end

    % Convert the input data to double for calculations
    data = double(data);

    % Get the number of rows in the input matrix
    [m, ~] = size(data);

    % Initialize the MI matrix
    mi_matrix = zeros(m, m);

    % Loop through each pair of rows to calculate MI
    for i = 1:m
        for j = i:m
            % Calculate mutual information for the pair of rows
            mi_value = compute_mi(data(i, :), data(j, :), L);

            % Store the MI value in the symmetric matrix
            mi_matrix(i, j) = mi_value;
            mi_matrix(j, i) = mi_value; % Since MI is symmetric
        end
    end
end

% -----------------------------------
% Helper function: Compute Mutual Information
function I = compute_mi(A, B, L)
    % Compute the mutual information between two signals A and B
    % A, B: Input signals (vectors)
    % L: Number of bins for the histogram

    % Calculate histograms for A and B
    na = histcounts(A(:), L, 'Normalization', 'probability'); % Probability of A
    nb = histcounts(B(:), L, 'Normalization', 'probability'); % Probability of B

    % Calculate the joint histogram
    n2 = compute_hist2(A, B, L);
    n2 = n2 / sum(n2(:)); % Normalize joint histogram

    % Calculate the mutual information
    I = sum(minf(n2, na' * nb));
end

% -----------------------------------
% Helper function: 2D Histogram
function n = compute_hist2(A, B, L)
    % Compute the joint histogram of two signals A and B using L bins

    % Determine the min and max values for A and B
    ma = min(A(:)); MA = max(A(:));
    mb = min(B(:)); MB = max(B(:));

    % Scale and round A and B to fit in the range {0, ..., L-1}
    A = round((A - ma) * (L - 1) / (MA - ma + eps));
    B = round((B - mb) * (L - 1) / (MB - mb + eps));

    % Initialize joint histogram
    n = zeros(L);

    % Compute the joint histogram using indexing
    x = 0:L-1;
    for i = 0:L-1
        counts = histc(B(A == i), x); % Count occurrences
        n(i+1, :) = counts(:)';      % Transpose to match row-wise assignment
    end
end

% -----------------------------------
% Helper function: Mutual Information Calculation
function y = minf(pab, papb)
    % Calculate the contribution to mutual information
    % pab: Joint probability distribution
    % papb: Product of marginal distributions

    % Find indices where both probabilities are non-zero
    I = find(papb(:) > 1e-12 & pab(:) > 1e-12);

    % Compute the contribution to mutual information
    y = pab(I) .* log2(pab(I) ./ papb(I));
end
