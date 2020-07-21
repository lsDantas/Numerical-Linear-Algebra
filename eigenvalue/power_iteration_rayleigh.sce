///////////////////////////////////////////////////////////
// Power Iteration with Rayleigh Quotient
//
// Description: Finds approximations for the dominant 
// eigenvalue and eigenvector pair of a matrix A using 
// the power iteration method. Uses the Rayleigh 
// quotient for faster convergence.
///////////////////////////////////////////////////////////
// Input:
//   A:        a real n x n diagonalizable matrix 
//             with a dominant eigenvalue
//   x0:       initial approximation for dominant 
//             eigenvector
//   epsilon:  precision used as halt criterion
//   max_iter: maximum number of iterations
///////////////////////////////////////////////////////////
// Output:
//   lambda:      dominant eigenvalue of A
//                (approximation)
//   x1:          dominant eigenvector of A
//                (approximation)
//   counter:     number of iterations executed
//   norm_error:  Euclidean norm of difference
//                between subsequent iterations
///////////////////////////////////////////////////////////

function [lambda, x1, counter, norm_error] = power_iteration_rayleigh(A, x0, epsilon, max_iter)
    // First Iteration
    counter = 1;
    x0 = x0/norm(x0,2);
    x1 = A*x0;

    // Other Iterations
    while counter <= max_iter
        // Approximate lambda with Rayleight quotient
        lambda = x0'*x1;

        // Normalize x1
        x1 = x1/norm(x1,2);

        // Verify precision
        norm_error = norm(x1 - x0, 2);
        if norm_error < epsilon
            disp("Desired precision reached!");
            return;
        end

        // Update dominant eigenvector approximation
        x0 = x1;
        x1 = A*x0;

        // Update counter
        counter = counter + 1;
    end
endfunction