import cv2
import numpy as np

# Replace "path/to/your/images" with the actual directory containing your images
image_paths = [f"c_code/figures/fullSystem_ID/PI_basicTest_{i}.png" for i in range(1,16)]

# Load the images
images = [cv2.imread(path) for path in image_paths]

# Get the dimensions of the first image
height, width, channels = images[0].shape

# Create a blank image for the grid
grid_image = np.zeros((height * 4, width * 4, channels), dtype=np.uint8)

# Arrange the images in the grid
for i in range(4):
    for j in range(4):
        print(i * 4 + j)
        if (i * 4 + j) < 15:
            grid_image[i * height:(i + 1) * height, j * width:(j + 1) * width] = images[i * 4 + j]

# Save the grid image
cv2.imwrite("grid_image.png", grid_image)