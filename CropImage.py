def crop_im():

    with Image.open("shape_{}.png".format(i)) as im:
        width, height = im.size
        print(width)
        print(height)
        # The crop method from the Image module takes four coordinates as input.
        # The right can also be represented as (left+width)
        # and lower can be represented as (upper+height).
        left = 64
        right = 700
        upper = 64
        lower = 700

        # Here the image "im" is cropped and assigned to new variable im_crop
        im_crop = im.crop((left, upper, right, lower))
        im_crop.save("shape_{}.png".format(i+1000), "png", quality=100)

for i in range(1,486):
    if __name__ == "__main__":
        crop_im()
    
    
