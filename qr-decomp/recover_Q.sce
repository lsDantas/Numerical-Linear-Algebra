///////////////////////////////////////////////////////////
// Recover Q
//
// Description: Produces an orthogonal matrix Q
// from the U matrix produced by the Householder 
// method.
///////////////////////////////////////////////////////////
// Input:
//   U: a m x n matrix produced by the Householder
//      algorithm
///////////////////////////////////////////////////////////
// Output:
//   Q: an m x m orthogonal matrix
///////////////////////////////////////////////////////////

function [Q] = recover_Q(U)
    // Determine matrix dimensions and initialize Q
    [m n]=size(U);
    Q = eye(m,m);

    // Compute Q from U
    for i=n:-1:1
        // Build Householder matrix H_i
        H = eye(m,m); 
        H(i:m,i:m) = H(i:m,i:m) - 2 * U(i:m, i) * U(i:m, i)';
        
        // Update Q
        Q = H*Q;
    end

endfunction