///////////////////////////////////////////////////////////
// SVD Image Compression
//
// Description: Compresses an image using singular 
// value decomposition. Implementated using native
// Scilab matrix operations.
///////////////////////////////////////////////////////////
// Input:
//   A: a black-and-white image
//   p: compression level expressed as a percentage
//      (a number between 0 and 1)
///////////////////////////////////////////////////////////
// Output:
//   C: a compressed version of A
///////////////////////////////////////////////////////////

function [C] = svd_compress_parallel(A,p)
    // Convert image to real matrix
    A = double(A);
    [m n] = size(A);

    // Obtain singular value decomposition
    [U, S, V, r] = svd(A);

    // Establish number of eigenvectors to be used
    s = max(1, floor(r*p));

    // Rebuild image using largest eigenvectors
    C = U(:,1:s)*S(1:s,1:s)*V'(1:s,:);

    // Convert C to integer matrix
    C = iconvert(C,11);
endfunction




