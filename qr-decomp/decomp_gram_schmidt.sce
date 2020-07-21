///////////////////////////////////////////////////////////
// Gram-Schmidt Process
//
// Description: Orthonormalizes a set of vectors
// corresponding to the columns of a matrix A.
// Yields QR decomposition of A when applied to
// full column rank matrix. Traditional 
// implementation.
///////////////////////////////////////////////////////////
// Input:
//   A: a m x n full column rank matrix
///////////////////////////////////////////////////////////
// Output:
//   Q: an orthogonal matrix
//   R: an upper triangular matrix
///////////////////////////////////////////////////////////

function [Q,R] = decomp_gram_schmidt(A)
    // Determine matrix dimensions
    [m n]=size(A);

    // Initialize Q and R
    Q = eye(m,m);
    R = zeros(m,n);

    // Find Orthonormal Vectors
    for i=1:n
        // Select current vector
        Q(:,i) = A(:,i);

        // Identify vector projections
        for j=1:(i-1)
            // Fill i-th column of R
            R(j,i) = A(:,i)'*Q(:,j);

            // Subtract projection from vector
            Q(:,i) = Q(:,i) - R(j,i)*Q(:,j);
        end

        // Find vector norm
        R(i,i) = norm(Q(:,i));

        // Update Q
        Q(:,i) = Q(:,i)/R(i,i);
    end
endfunction