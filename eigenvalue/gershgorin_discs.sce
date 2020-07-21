///////////////////////////////////////////////////////////
// Gershgorin discs
//
// Description: Finds Gershgorin discs associated
// with a matrix
///////////////////////////////////////////////////////////
// Input:
//   A: a n x n matrix with real eigenvalues
///////////////////////////////////////////////////////////
// Output:
//   c: vector with centers of discs
//   r: vector radii of discs
///////////////////////////////////////////////////////////

function [c,r] = gershgorin_discs(A)
    // Obtain matrix dimensions
    [n]=size(A,1);

    // Determine centers from diagonal
    c = diag(A);

    // Determine radii
    r = zeros(n,1);
    for i=1:n
        for j=1:n
            // Sum elements outside diagonal
            if i ~= j
                r(i,1) = r(i,1) + A(i,j);
            end
        end
    end
endfunction