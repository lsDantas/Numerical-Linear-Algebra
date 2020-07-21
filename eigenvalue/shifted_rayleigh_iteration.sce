///////////////////////////////////////////////////////////
// Shifted Rayleigh Quotient Iteration
// 
// Description: Approximates closest eigenvalue to a
// particular value along with its corresponding 
// eigenvector. Uses the shifted rayleigh quotient
// iteration method.
///////////////////////////////////////////////////////////
// Input:
//   A:        a real n x n diagonalizable matrix 
//   x0:       initial approximation for eigenvector
//   epsilon:  precision used as halt criterion
//   alpha:    reference value for eigenvalue search
//   max_iter: maximum number of iterations
///////////////////////////////////////////////////////////
// Output:
//   lambda:      closest eigenvalue to alpha
//   x1:          eigenvector corresponding to 
//                lambda
//   counter:     number of iterations executed
//   norm_error:  Euclidean norm of difference
//                between subsequent iterations
///////////////////////////////////////////////////////////

function [lambda, x1, counter, norm_error] = shifted_rayleigh_iteration(A, x0, epsilon, alpha, max_iter)
    // Determine matrix dimensions and start counter
    [n]=size(A,1);
    counter = 1;

    // Take alpha as initial lambda
    lambda = alpha;

    // Normalize x0
    x0 = x0/norm(x0,2);

    // Method iterations
    while counter <= max_iter
        // Update eigenvector through linear system
        x1 = solve_gaussian_elimination(A - lambda * eye(n,n), x0);
        x1 = x1/norm(x1,2);

        // Approximate lambda
        lambda = x1'*(A*x1);

        // Verify precision
        norm_error = norm(x1 - x0, 2);
        if norm_error < epsilon
            disp("Desired precision reached!");
            return;
        end

        // Store iteration and update counter
        x0 = x1;
        counter = counter + 1;
    end
endfunction