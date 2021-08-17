import numpy
from PIL import Image

def pyramid_denoise(noisy_image):
    input_image = noisy_image
    w, h = input_image.size
    circle = np.zeros(w, h)
    for x in range(0, w+1):
        for y in range(0, h+1):
            if numpy.sqrt((x-(w/2))^2 + (y-h/2)^2) > 180:
                circle(x,y) = 1
    numpy.where(input_image==1, 0, input_image)
    numpy.where(circle==0, 0, input_image)
    #bad_pixels = 