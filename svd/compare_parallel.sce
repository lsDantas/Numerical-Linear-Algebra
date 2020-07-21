///////////////////////////////////////////////////////////
// Compressed Image Comparisons (Parallel)
//
// Description: Displays an image compressed to 
// different degrees. Uses parallel implementation.
///////////////////////////////////////////////////////////
// Input:
//   original_image:    string with path to image to be 
//                      compressed
//   compression_rates: row vector of compression rates
//                      expressed as numbers from 0 to 1
///////////////////////////////////////////////////////////

function compare_parallel(original_image, compression_rates)
    // Read image and convert to black-and-white
    A = imread(original_image);
    if size(size(A))(2) == 3
        A = rgb2gray(A);
    end

    // Sort compression rates in increasing order
    n = size(compression_rates,2);
    compression_rates = gsort(compression_rates, 'g', 'i');

    // Build display grid
    dim = ceil(sqrt(n + 1));
    subplot(dim, dim, 1);

    // Display original image
    imshow(A);
    title("Original Image");

    // Compress and display images
    for i=1:n
        // Compress
        C = svd_compress_parallel(A, compression_rates(i));

        // Display
        subplot(dim, dim, i+1);
        title("p = " + string(compression_rates(i)));
        imshow(C);
    end

endfunction