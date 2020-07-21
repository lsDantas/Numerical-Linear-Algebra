///////////////////////////////////////////////////////////
// QR Spectrum
//
// Description: Finds eigenvalues of given matrix A
// through the QR algorithm.
///////////////////////////////////////////////////////////
// Input:
//   A: a m x n full column rank matrix
///////////////////////////////////////////////////////////
// Output:
//   S: a vector with the eigenvalues of A
///////////////////////////////////////////////////////////

function [S] = qr_spectrum(A)
    // Initialize iteration count
    counter = 0;
    max_iter = 1000;

    // QR algorithm
    while counter <= max_iter
        // Obtain QR decomposition of A
        [U R] = decomp_householder(A);
        Q = recover_Q(U);

        // Update A and counter
        A = R*Q;
        counter = counter + 1; 
    end
    
    // Obtain eigenvalues from diagonal
    S = diag(A);
endfunction