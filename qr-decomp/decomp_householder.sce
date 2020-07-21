///////////////////////////////////////////////////////////
// Householder Method
//
// Description: Implements the Householder decomposition
// method which can further be refined for QR 
// decomposition.
///////////////////////////////////////////////////////////
// Input:
//   A: a m x n full column rank matrix
///////////////////////////////////////////////////////////
// Output:
//   U: a residual matrix from which an orthogonal matrix
//      Q may later be obtained
//   R: an upper triangular matrix
///////////////////////////////////////////////////////////

function [U,R] = decomp_householder(A)
    // Determine matrix dimensions and initialize U
    [m n]=size(A);
    U = zeros(m,n);

    // Generate Householder reflection vectors
    for i=1:n
        // Select vector to be Householder reflected
        x = A(i:m, i);

        // Ensure obtuse angle with axis
        if x(1) < 0
            x(1) = x(1) - norm(x);
        else
            x(1) = x(1) + norm(x);
        end

        // Normalize vector and update U
        u = x/norm(x);
        U(i:m, i) = u;
        
        // Update A using Householder reflection
        A(i:m, i:n) = A(i:m, i:n) - 2 * u * (u'*A(i:m,i:n));
    end

    // Extract Upper Triangular Matrix
    R = triu(A);
endfunction