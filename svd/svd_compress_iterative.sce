///////////////////////////////////////////////////////////
// SVD Image Compression
//
// Description: Compresses an image using singular 
// value decomposition. Implementated using a for 
// loop.
///////////////////////////////////////////////////////////
// Input:
//   A: a black-and-white image
//   p: compression level expressed as a percentage
//      (a number between 0 and 1)
///////////////////////////////////////////////////////////
// Output:
//   C: a compressed version of A
///////////////////////////////////////////////////////////

function [C] = svd_compress_iterative(A,p)
    // Convert image to real matrix
    A = double(A);
    [m n] = size(A);

    // Obtain singular value decomposition
    [U, S, V, r] = svd(A);

    // Establish number of eigenvectors to be used
    s = max(1, floor(r*p));

    // Rebuild image using largest eigenvectors
    C = zeros(m,n);
    for k=1:s
        C = C + S(k,k)*U(:,k)*V'(k,:)
    end

    // Convert C to integer matrix
    C = iconvert(C,11);
endfunction




