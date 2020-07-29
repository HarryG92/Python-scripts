# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 20:07:37 2020

A jpeg compression script I wrote as an application of Fourier analysis for
a student I was tutoring

@author: gulli
"""

import numpy
import time
import cv2


class my_JPEG:
    def __init__(self, compression_vector, height, width):
        self.compression_vector = compression_vector
        self.height = height
        self.width = width
        self.image = numpy.zeros((height, width, 3))

    def set_block(self, block_row, block_column, component, block):
        self.image[
                   block_row:(block_row + 8),
                   block_column:(block_column + 8),
                   component
                   ] = block


class my_JPEG_compressor:
    def __init__(self, file_name):
        # import image
        self.original_image = cv2.imread(file_name)
        self.working_image = self.original_image

        # extract image dimensions
        self.true_width = len(self.original_image[0, :, 0])
        self.true_height = len(self.original_image[:, 0, 0])

        # crop so dimensions are multiples of 8
        self.width = self.true_width - (self.true_width % 8)
        self.height = self.true_height - (self.true_height % 8)
        self.working_image = self.working_image[
                                                0:self.height,
                                                0:self.width,
                                                :
                                                ].astype(numpy.uint8)

        # initialise matrix for discrete cosine transform
        self.transform_matrix = numpy.zeros((8, 8))

        # first row is all 1/sqrt(8)
        self.transform_matrix[0, :] = numpy.ones((1, 8)) / numpy.sqrt(8)

        # remaining entries are cosines
        for row in range(1, 8):
            for column in range(0, 8):
                self.transform_matrix[row, column] = (
                                                      (1/2)
                                                      * numpy.cos(
                                                                  2
                                                                  * column
                                                                  + 1
                                                                  )
                                                      * row
                                                      * numpy.pi
                                                      / 16
                                                      )

        # initialise JPEG standard quantisation matrix
        self.quantisation_matrix = numpy.array(
                [[16, 11, 10, 16, 24, 40, 51, 61],
                 [12, 12, 14, 19, 26, 58, 60, 55],
                 [14, 13, 16, 24, 40, 57, 69, 56],
                 [14, 17, 22, 29, 51, 87, 80, 62],
                 [18, 22, 37, 56, 68, 109, 103, 77],
                 [24, 35, 55, 64, 81, 104, 113, 92],
                 [49, 64, 78, 87, 103, 121, 120, 101],
                 [72, 92, 95, 98, 112, 100, 103, 99]
                 ])

    # display current working image
    def display(self):
        cv2.imshow('Current working image', self.working_image)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    # carry out discrete cosine transform on an 8x8 block
    def cosine_transform(self, matrix):
        return numpy.dot(
                         self.transform_matrix,
                         numpy.dot(
                                   matrix,
                                   numpy.transpose(self.transform_matrix)
                                   )
                         )

    # carry out inverse discrete cosine transform on an 8x8 block
    def inverse_cosine_transform(self, matrix):
        return numpy.dot(
                         numpy.transpose(self.transform_matrix),
                         numpy.dot(
                                   matrix,
                                   self.transform_matrix
                                   )
                         ).astype(numpy.uint8)

    # scale original image by a given percentage
    def resize(self, scale_percent):
        new_width = int(self.original_image.shape[1] * scale_percent / 100)
        new_height = int(self.original_image.shape[0] * scale_percent / 100)
        dim = (new_width, new_height)
        # resize image
        self.working_image = cv2.resize(
                                        self.original_image,
                                        dim,
                                        interpolation=cv2.INTER_AREA
                                        )

        # crop so dimensions are multiples of 8
        self.width = new_width - (new_width % 8)
        self.height = new_height - (new_height % 8)
        self.working_image = self.working_image[
                                                0:self.height,
                                                0:self.width,
                                                :
                                                ].astype(numpy.uint8)

    def convert_to_YCrCb(self, jpeg):
        return cv2.cvtColor(jpeg, cv2.COLOR_BGR2YCR_CB).astype(numpy.uint8)

    def convert_to_BGR(self, jpeg):
        return cv2.cvtColor(jpeg, cv2.COLOR_YCR_CB2BGR).astype(numpy.uint8)

    def compress(self, compression):
        # initialise compressed image
        compressed_image = my_JPEG(
                                   compression,
                                   self.height,
                                   self.width
                                   )

        # convert to YCrCb
        YCrCb_image = self.convert_to_YCrCb(self.working_image)

        # create array to keep track of how many more zeros we get
        non_zero = numpy.zeros((2, 3))

        start_time = time.time()  # for timing compression

        # compression algorithm with Discrete Cosine Transform
        for component in range(3):  # loop through 3 components of YCrCb
            # isolate component in image
            image_component = YCrCb_image[:, :, component]

            # count non-zero entries in original for comparison with compressed
            # version later
            non_zero[0, component] = numpy.count_nonzero(image_component)

            # translate scale from 0-255 to -128-127 for DCT
            image_component -= 128

            # initialise quantisation matrix for this component
            component_quantisation_matrix = (
                                             self.quantisation_matrix
                                             * compression[component]
                                             )

            # break up into 8x8 blocks and work with those
            for block_row in range(0, self.height, 8):
                for block_column in range(0, self.width, 8):
                    current_block = image_component[
                                             block_row:(block_row + 8),
                                             block_column:(block_column + 8)
                                             ]

                    # take the DCT of the current block
                    transformed_block = self.cosine_transform(current_block)

                    # quantise by dividing entrywise by the quantisation matrix
                    # and taking floor function
                    quantised_block = (
                                       transformed_block
                                       // component_quantisation_matrix
                                       )
                    # count non-zero entries in compressed image for comparison
                    non_zero[1, component] += numpy.count_nonzero(
                                                        quantised_block)

                    # paste 8x8 block in with others in compressed image
                    compressed_image.set_block(
                                               block_row,
                                               block_column,
                                               component,
                                               quantised_block
                                               )

        print("Compression completed in %s seconds" % (time.time()-start_time))
        print("Compression percentages by component (Y, Cr, Cb):")
        number_of_entries = self.original_image[:, :, 0].size
        for component in range(3):
            old_zeros = number_of_entries - non_zero[0, component]
            new_zeros = number_of_entries - non_zero[1, component]
            percentage_change = (
                                 (new_zeros - old_zeros)
                                 / number_of_entries
                                 * 100
                                 )
            print("%.0f" % percentage_change)
        return compressed_image

    def decompress(
                   self,
                   compressed_jpeg
                   ):

        # read off compression factors
        compression = compressed_jpeg.compression_vector

        # initialise decompressed image
        decompressed_image = my_JPEG(compression, self.height, self.width)

        start_time = time.time()  # for timing decompression

        # compression algorithm with Discrete Cosine Transform
        for component in range(3):  # loop through 3 components of YCrCb
            # isolate component in image
            image_component = compressed_jpeg.image[:, :, component]

            # translate scale from 0-255 to -128-127 for DCT
            image_component -= 128

            # initialise quantisation matrix for this component
            component_quantisation_matrix = (
                                             self.quantisation_matrix
                                             * compression[component]
                                             )

            # break up into 8x8 blocks and work with those
            for block_row in range(0, self.height, 8):
                for block_column in range(0, self.width, 8):
                    current_block = image_component[
                                             block_row:(block_row + 8),
                                             block_column:(block_column + 8)
                                             ]

                    # dequantise by multiplying by the quantisation matrix
                    dequantised_block = numpy.multiply(
                                               current_block,
                                               component_quantisation_matrix
                                               )

                    # inverse discrete cosine transform
                    detransformed_block = self.inverse_cosine_transform(
                            dequantised_block)

                    detransformed_block += 128  # revert to 0-255 scale

                    # paste 8x8 block in with others in decompressed image
                    decompressed_image.set_block(
                                                 block_row,
                                                 block_column,
                                                 component,
                                                 detransformed_block
                                                 )

        print(
              "Decompression completed in %s seconds"
              % (time.time() - start_time)
              )
        return decompressed_image

    def compress_and_compare(
                             self,
                             luminance_compression,
                             red_compression,
                             blue_compression
                             ):

        compression = [
                       luminance_compression,
                       red_compression,
                       blue_compression]

        compressed_image = self.compress(compression)
        decompressed_image = self.convert_to_YCrCb(
                self.decompress(compressed_image).image.astype(numpy.uint8))

        # display images
        cv2.imshow('original', self.working_image)
        cv2.imshow('after compression', decompressed_image)
        cv2.waitKey(0)
        cv2.destroyAllWindows()









image = my_JPEG_compressor("TestImage.jpg")
image.resize(10)
# image.display()

image.compress_and_compare(0.00001, 0.0001, 0.00001)
