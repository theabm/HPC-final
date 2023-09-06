# this is a python file that takes a certain number of test cases - 
# think unit tests - and outputs the corresponding pgm and the expected pgm for 
# 1 step of evolution (computed by hand)

# these generated files will be used for testing. The pgm will be passed to our 
# code, and its output will be compared with the expected pgm file

# cases:

# Rules:
# 1) a dead cell with three live neighbors comes to life.
# 2) a live cell with 2 or 3 live neighbors survives
# 3) a cell with less than 2 or more than 3 live neighbors becomes/stays dead

# the examples each try to only use one of the rules, to test them separately. 
# we start with a basic example, then we try a few examples on the borders 
# to test infinite boundaries.

rule1_1 = [
        0,0,0,0,
        1,1,0,0,
        0,0,0,0,
        0,1,0,0,
        ]
rule1_1_expected = [
        0,0,0,0,
        0,0,0,0,
        0,1,0,0,
        0,0,0,0,
        ]

rule1_2 = [
        0,0,0,0,
        0,0,1,0,
        0,0,0,0,
        1,0,0,1,
        ]
rule1_2_expected = [
        0,0,0,1,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        ]

rule1_3 = [
        1,0,1,0,
        0,0,0,0,
        0,1,0,0,
        0,0,0,0,
        ]
rule1_3_expected = [
        0,0,0,0,
        0,1,0,0,
        0,0,0,0,
        0,1,0,0,
        ]


# 2) a live cell with 2 or 3 live neighbors survives
rule2_1 = [
        0,0,0,0,
        0,1,0,0,
        0,1,1,0,
        0,0,0,0,
        ]
rule2_1_expected = [
        0,0,0,0,
        0,1,1,0,
        0,1,1,0,
        0,0,0,0,
        ]

rule2_2 = [
        0,0,0,0,
        0,1,0,0,
        0,1,0,0,
        0,1,0,0,
        ]
rule2_2_expected = [
        0,0,0,0,
        0,0,0,0,
        1,1,1,0,
        0,0,0,0,
        ]

rule2_3 = [
        0,0,0,0,
        0,1,1,0,
        0,1,1,0,
        0,0,0,0,
        ]
rule2_3_expected = [
        0,0,0,0,
        0,1,1,0,
        0,1,1,0,
        0,0,0,0,
        ]

rule2_4 = [
        1,0,0,0,
        1,0,0,1,
        0,0,0,0,
        0,0,0,0,
        ]
rule2_4_expected = [
        1,0,0,1,
        1,0,0,1,
        0,0,0,0,
        0,0,0,0,
        ]

rule2_5 = [
        1,0,0,1,
        1,0,0,1,
        0,0,0,0,
        0,0,0,0,
        ]
rule2_5_expected = [
        1,0,0,1,
        1,0,0,1,
        0,0,0,0,
        0,0,0,0,
        ]

rule2_6 = [
        0,1,1,0,
        0,0,0,0,
        0,0,0,0,
        0,1,1,0,
        ]
rule2_6_expected = [
        0,1,1,0,
        0,0,0,0,
        0,0,0,0,
        0,1,1,0,
        ]

rule2_7 = [
        1,0,0,1,
        0,0,0,0,
        0,0,0,0,
        1,0,0,1,
        ]
rule2_7_expected = [
        1,0,0,1,
        0,0,0,0,
        0,0,0,0,
        1,0,0,1,
        ]

# there is no need to test rule 3 as it is indirectly tested in all the other 
# examples

# now we test some more complex cases that are documented

glider = [
        0,0,0,0,0,
        0,0,1,0,0,
        0,0,0,1,0,
        0,1,1,1,0,
        0,0,0,0,0,
        ]

glider_1 = [
        0,0,0,0,0,
        0,0,0,0,0,
        0,1,0,1,0,
        0,0,1,1,0,
        0,0,1,0,0,
        ]

glider_2 = [
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,1,0,
        0,1,0,1,0,
        0,0,1,1,0,
        ]

glider_3 = [
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,1,0,0,
        0,0,0,1,1,
        0,0,1,1,0,
        ]

glider_4 = [
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,1,0,
        0,0,0,0,1,
        0,0,1,1,1,
        ]

glider_5 = [
    #   0,0,0,1,1,0
        0,0,0,1,0,#0
        0,0,0,0,0,#0
        0,0,0,0,0,#0
        0,0,1,0,1,#0
        0,0,0,1,1,#0
    #   0,0,0,1,0,0
        ]

glider_6 = [
        0,0,0,1,1,#0
        0,0,0,0,0,#0
        0,0,0,0,0,#0
        0,0,0,0,1,#0
        0,0,1,0,1,#0
    #   0,0,0,1,1,0,
        ]

glider_7 = [
        0,0,0,1,1,#0
        0,0,0,0,0,#0
        0,0,0,0,0,#0
        0,0,0,1,0,#0
        1,0,0,0,1,#1
    #   0,0,0,1,1,0,
        ]

glider_8 = [
        1,0,0,1,1,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,1,
        1,0,0,0,0,
        ]

glider_9 = [
        1,0,0,0,1,
        0,0,0,0,1,
        0,0,0,0,0,
        0,0,0,0,0,
        1,0,0,1,0,
        ]

glider_10 = [
        1,0,0,1,0,
        1,0,0,0,1,
        0,0,0,0,0,
        0,0,0,0,0,
        1,0,0,0,0,
        ]

# time step 15
glider_15 = [
        1,0,0,0,0,
        0,1,1,0,0,
        1,1,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        ]

def list_to_pgm(input_list, width, height, output_file):
    if len(input_list) != width * height:
        raise ValueError("Input list size doesn't match specified width and height")

    with open(output_file, "wb") as f:
        # Write PGM header
        f.write(b"P5\n")
        f.write(f"{width} {height}\n".encode())
        f.write(b"255\n")

        # Write pixel values
        for pixel_value in input_list:
            f.write(bytes([pixel_value]))

# Example usage:
custom_list = [0, 1, 0, 1, 1, 0, 0, 1, 0]  # Replace with your custom list of
# zeros and ones
width = 3  # Replace with the width of your image
height = 3  # Replace with the height of your image
output_file_name = "output.pgm"

list_to_pgm(custom_list, width, height, output_file_name)

