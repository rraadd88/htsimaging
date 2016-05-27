


@pims.pipeline
def average_z(image):
    return image.mean(axis=0)  # the same as image[0] + ... + image[4]
