///////////////////////////////////////////////////////////
// Modified Gram-Schmidt Process
//
// Description: Orthonormalizes a set of vectors
// corresponding to the columns of a matrix A.
// Yields QR decomposition of A when applied to
// full column rank matrix. Modified implementation
// for greater numerical stability.
///////////////////////////////////////////////////////////
// Input:
//   Q: a m x n full column rank matrix
///////////////////////////////////////////////////////////
// Output:
//   Q: an orthogonal matrix
//   R: an upper triangular matrix
///////////////////////////////////////////////////////////

function [Q,R] = decomp_modified_gram_schmidt(Q)
    // Determine matrix dimensions
    [m n]=size(Q);

    // Initialize R
    R = zeros(m,m);

    // Find Orthonormal Vectors
    for i=1:n
        // Identify vector projections
        for j=1:(i-1)
            // Fill i-th column of R
            R(j,i) = Q(:,j)'*Q(:,j);

            // Subtract projection from vector
            Q(:,i) = Q(:,i) - R(j,i)*Q(:,j);
        end
        
        // Find vector norm
        R(i,i) = norm(Q(:,i));

        // Update Q
        Q(:,i) = Q(:,i)/R(i,i);
    end

endfunction