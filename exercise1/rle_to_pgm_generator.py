from PIL import Image
import numpy as np

# Define the size of the grid, generation number, and custom header
grid_size = 65  # Adjust this based on your desired grid size

# Define the Snark loop pattern for each generation (in RLE format)
pgm_patterns = {
    "snark_loop":"27b2o$27bobo$29bo4b2o$25b4ob2o2bo2bo$25bo2bo3bobob2o$28bobobobo$29b2obobo$33bo2$19b2o$20bo8bo$20bobo5b2o$21b2o$35bo$36bo$34b3o2$25bo$25b2o$24bobo4b2o22bo$31bo21b3o$32b3o17bo$34bo17b2o2$45bo$46b2o12b2o$45b2o14bo$3b2o56bob2o$4bo9b2o37bo5b3o2bo$2bo10bobo37b2o3bo3b2o$2b5o8bo5b2o35b2obo$7bo13bo22b2o15bo$4b3o12bobo21bobo12b3o$3bo15b2o22bo13bo$3bob2o35b2o5bo8b5o$b2o3bo3b2o37bobo10bo$o2b3o5bo37b2o9bo$2obo56b2o$3bo14b2o$3b2o12b2o$19bo2$11b2o17bo$12bo17b3o$9b3o21bo$9bo22b2o4bobo$38b2o$39bo2$28b3o$28bo$29bo$42b2o$35b2o5bobo$35bo8bo$44b2o2$31bo$30bobob2o$30bobobobo$27b2obobo3bo2bo$27bo2bo2b2ob4o$29b2o4bo$35bobo$36b2o!"
}

# Create a blank PGM image
img = Image.new("L", (grid_size, grid_size))
PAD = 4

# Function to convert RLE format to a grid of 0s and 1s
def rle_to_grid(rle):
    grid = []
    mod_grid_size = grid_size + 2*PAD

    for i in range(PAD):
        row = []
        row.extend([0]*mod_grid_size)
        grid.append(row)

    count_rows = len(rle.split('$'))

    for line in rle.split('$'):
        row = []
        repeat = ""
        count = grid_size
        row.extend([0]*PAD)
        for token in line:
            if token.isnumeric():
                repeat = repeat + token
                continue
            elif token == 'b':
                if len(repeat) == 0:
                    repeat = 1
                row.extend([0]*int(repeat))
                count -= int(repeat)
                repeat=""
                continue
            elif token == 'o':
                if len(repeat) == 0:
                    repeat = 1
                row.extend([1]*int(repeat))
                count -= int(repeat)
                repeat=""
                continue

        row.extend([0]*count)
        row.extend([0]*PAD)
        if len(row) != mod_grid_size:
            print("problem")

        grid.append(row)

    missing_rows = grid_size - count_rows

    for i in range(missing_rows):
        row = []
        row.extend([0]*grid_size)
        grid.append(row)

    for i in range(PAD):
        row = []
        row.extend([0]*mod_grid_size)
        grid.append(row)

    if np.shape(grid) != (mod_grid_size,mod_grid_size):
        print("problem")

    return grid

# Convert the RLE pattern for the specified generation to a grid
pattern = pgm_patterns["snark_loop"]
grid = rle_to_grid(pattern)

# Open the file for writing with a custom header
with open("snark_loop.pgm", "wb") as f:
    # Write the PGM header
    header = f"P5 {grid_size} {grid_size} 1\n"
    f.write(header.encode())

    # Write the PGM pixel data (0 for alive, 255 for dead)
    for y in range(min(grid_size, len(grid))):
        for x in range(min(grid_size, len(grid[y]))):
            if grid[y][x] == 0:
                f.write(bytes([0]))  # Alive
            else:
                f.write(bytes([1]))  # Dead
