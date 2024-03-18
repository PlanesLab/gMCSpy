% Assuming you have a sparse matrix named 'sparseMatrix'
% Initialize a counter
counter = 0;

% Iterate through all rows of the matrix
for rowNumber = 1:size(G, 1)
    % Sum the columns of the current row
    rowSum = sum(G(rowNumber, :));
    
    % Check if the sum is greater than 5
    if rowSum == 1
        % If the sum is greater than 5, increment the counter
        counter = counter + 1;
    end
end

% Display the final counter value
fprintf('Total rows with : %d\n', counter);
