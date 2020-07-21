///////////////////////////////////////////////////////////
// Solve with LU
//
// Description: Solves linear systems of the 
// form Ax = B using the LU decomposition of A.
///////////////////////////////////////////////////////////
// Input: 
//   LU: the LU decomposition of A in a single matrix
//         for i > j: LU(i,j) = L(i,j) 
//         for j >= i: LU(i,j) = U(i,j)
//   B:  a n x m matrix
//   P:  the permutation matrix produced 
//       during LU decomposition 
//       (identity matrix if none)
///////////////////////////////////////////////////////////
// Output: 
//   x: the solution of the system 
//      (assuming such solution exists)
///////////////////////////////////////////////////////////

function [x] = solve_LU(LU, B, P)
    // Obtain matrix dimensions and adjust row order
    [n]=size(LU,1);
    [m]=size(B, 2);
    B = P*B;

    // Split L and U from LU
    L = tril(LU) - eye(LU).*LU + eye(n,n);
    U = triu(LU);

    // Solve Ly = B, with y = Ux
    y = zeros(n,m);
    // Forward substitution
    for i=1:n
        // Sums products of known terms of y with elements of L 
        sum_terms = zeros(1,m);
        for j=1:i-1
            sum_terms = sum_terms + L(i,j)*y(j, :);
        end
        
        // Obtain yi through known elements
        y(i,:) = (B(i,:) - sum_terms)/L(i,i);
    end

    // Solve Ux = y
    x = zeros(n,m);
    // Back Substitution
    for i=n:-1:1
        // Sums product of known terms of x with elements of X
        sum_terms = zeros(1,m);
        for j=n:-1:i+1
           sum_terms = sum_terms + U(i,j)*x(j,:);
        end
        
        // Obtain xi through known elements
        x(i,:) = (y(i,:) - sum_terms)/U(i,i);
    end

endfunction