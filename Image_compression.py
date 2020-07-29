import numpy

import time

import cv2


# initialise matrix for discrete cosine transform
transform_matrix = numpy.zeros((8, 8))

# first row is all 1/sqrt(8)
transform_matrix[0, :] = numpy.ones((1, 8)) / numpy.sqrt(8)

# remaining entries are cosines
for row in range(1, 8):
    for column in range(0, 8):
        transform_matrix[row, column] = (1/2) * numpy.cos(
                                                          (2 * (column)
                                                           + 1
                                                           )
                                                          * (row)
                                                          * numpy.pi / 16
                                                          )

# initialise JPEG standard quantisation matrix
quantisation_matrix = numpy.array([[16, 11, 10, 16, 24, 40, 51, 61],
                                   [12, 12, 14, 19, 26, 58, 60, 55],
                                   [14, 13, 16, 24, 40, 57, 69, 56],
                                   [14, 17, 22, 29, 51, 87, 80, 62],
                                   [18, 22, 37, 56, 68, 109, 103, 77],
                                   [24, 35, 55, 64, 81, 104, 113, 92],
                                   [49, 64, 78, 87, 103, 121, 120, 101],
                                   [72, 92, 95, 98, 112, 100, 103, 99]
                                   ])

# import image
original_image = cv2.imread('TestImage.jpg')

# resize image (comment out if not wanted)
scale_percent = 50  # percent of original size
new_width = int(original_image.shape[1] * scale_percent / 100)
new_height = int(original_image.shape[0] * scale_percent / 100)
dim = (new_width, new_height)
# resize image
original_image = cv2.resize(original_image, dim, interpolation=cv2.INTER_AREA)
# cv2.imwrite('Resized_original.jpg', original_image)

# convert from RGB to YCbCr format
YCbCr_image = cv2.cvtColor(original_image, cv2.COLOR_BGR2YCR_CB)

# extract image dimensions
width = len(YCbCr_image[0, :, 0])
height = len(YCbCr_image[:, 0, 0])

# crop so dimensions are multiples of 8
width = width - (width % 8)
height = height - (height % 8)

YCbCr_image = YCbCr_image[0:height, 0:width, :].astype(numpy.uint8)

# set compression factors for Y, Cb, and Cr components, from 0.0001 to 50
Yq = 0.0001
Cbq = 0.01
Cbr = 0.01

# combine quality factors into array
quality = [Yq, Cbq, Cbr]

# initialise compressed image
compressed_image = numpy.zeros((height, width, 3))
non_zero = numpy.zeros((2, 3))  # to keep track of how many more zeros we get

start_time = time.time()  # for timing compression

# compression algorithm with Discrete Cosine Transform
for component in range(len(quality)):  # loop through 3 components of YCbCr
    image_components = YCbCr_image[:, :, component]  # isolate component in img

    non_zero[0, component] = numpy.count_nonzero(image_components)  # count
# non-zero entries in original for comparison with compressed version later

    image_components -= 128  # translate scale from 0-255 to -128-127 for DCT

    # initialise quantisation matrix for this component
    component_quantisation_matrix = quantisation_matrix * quality[component]

    # break up into 8x8 blocks and work with those
    for block_row in range(0, height, 8):
        for block_column in range(0, width, 8):
            current_block = image_components[
                                             block_row:(block_row + 8),
                                             block_column:(block_column + 8)
                                             ]

            transformed_block = numpy.dot(
                                          transform_matrix,
                                          numpy.dot(
                                                    current_block,
                                                    numpy.transpose(
                                                            transform_matrix)
                                                   )
                                         )

            quantised_block = transformed_block//component_quantisation_matrix

            # count non-zero entries in compressed image for comparison
            non_zero[1, component] += numpy.count_nonzero(quantised_block)

            dequantised_block = numpy.multiply(
                                               quantised_block,
                                               component_quantisation_matrix
                                               )

            detransformed_block = numpy.dot(
                                            numpy.transpose(transform_matrix),
                                            numpy.dot(
                                                      dequantised_block,
                                                      transform_matrix
                                                      )
                                            ).astype(numpy.uint8)

            compressed_image[
                             block_row:(block_row + 8),
                             block_column:(block_column + 8),
                             component
                             ] = detransformed_block

    compressed_image[:, :, component] += 128  # translate back to 0-255 scale

compressed_image = compressed_image.astype(numpy.uint8)

# convert back to RGB
RGB_compressed_image = cv2.cvtColor(compressed_image, cv2.COLOR_YCR_CB2BGR)

# display images
cv2.imshow('original', original_image)
cv2.imshow('compressed', RGB_compressed_image)
cv2.waitKey(0)
cv2.destroyAllWindows()
cv2.imwrite('CompressedImage.jpg', RGB_compressed_image)

print("Compression percentages by component (Y, Cr, Cb):")
number_of_entries = original_image[:, :, 0].size
for component in range(3):
    old_zeros = number_of_entries - non_zero[0, component]
    new_zeros = number_of_entries - non_zero[1, component]
    percentage_change = ((new_zeros - old_zeros) / number_of_entries) * 100
    print("%.0f" % percentage_change)

print("--- %s seconds ---" % (time.time() - start_time))
