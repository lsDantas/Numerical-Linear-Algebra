///////////////////////////////////////////////////////////
// Gaussian Elimination
//
// Description: Solves a system of linear equations 
// (Ax = B) using the Gaussian elimination algorithm.
///////////////////////////////////////////////////////////
// Input:
//   A: an n x n square matrix
//   B: an n x 1 vector
///////////////////////////////////////////////////////////
// Output:
//   x:  the solution of the system 
//       (assuming such solution exists)
//   LU: LU decomposition of A in a single matrix
//         for i > j: LU(i,j) = L(i,j) 
//         for j >= i: LU(i,j) = U(i,j)
///////////////////////////////////////////////////////////

function [x, LU, permutation_matrix] = solve_gaussian_elimination(A, b)
    // Concatenate Matrices
    LU = [A,b];
    [n] = size(LU,1);

    // Permutation matrix for putting LU's rows back in order
    permutation_matrix = eye(n,n);
    
    // Gaussian Elimination Process
    for j=1:(n-1)
        // Select pivot with greatest absolute value (heuristic)
        [greatest_pivot, id_pivot] = max(abs(LU(j:n,j)));
        id_pivot = id_pivot + j - 1;

        // Switch rows to change pivot
        LU([j id_pivot], :) = LU([id_pivot j], :);
        permutation_matrix([j id_pivot], :) = permutation_matrix([id_pivot j], :);
        
        // Row operations
        for i=(j+1):n
            // Identify multiplicative constant
            LU(i,j) = LU(i,j)/LU(j,j);

            // Sum rows
            LU(i,j+1:n+1) = LU(i,j+1:n+1) - LU(i,j)*LU(j,j+1:n+1);
        end
    end


    // Find x through back substitution
    x=zeros(n,1);
    x(n) = LU(n,n+1)/LU(n,n);
    for i=n-1:-1:1
        x(i) = (LU(i,n+1) - LU(i,i:n)*x(i:n)) / LU(i,i);
    end

    // Removes last column
    LU=LU(1:n,1:n);
endfunction