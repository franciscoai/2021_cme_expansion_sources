import numpy as np

def PILOOP(datad):

    """
        Computes the PIL line characteristic from a magnetogram:
            input = INPUT ---> SINGLE AR LOS MAGNETOGRAM
            OUTPUT ---> DICT WITH COORDINATES X0,Y0, AND PIL ANGLE THETA0
                        theta0 = tilt angle
                        x0,y0 = point where the line passes
        Author: Mariano Poisson
    """
    #  datad[abs(datad)<thr]=0.0
    sz2, sz1 = (np.shape(datad))
    x1 = np.linspace(0, int((sz1 - 1)), int((sz1)))
    y1 = np.linspace(0, int((sz2 - 1)), int((sz2)))
    xx, yy = np.meshgrid(x1, y1)

    # COMPUTE X POSITION OF BARYCENTERS
    xp = np.sum((datad > 0) * xx * (datad)) / np.sum((datad > 0) * (datad))
    xn = np.sum((datad < 0) * xx * (datad)) / np.sum((datad < 0) * (datad))

    # COMPUTE SIGN OF LEADING POLARITY
    hsm = 1
    if xp > xn:
        hsm = -1

    # COMPUTE TOTAL FLUX WEIGHTED CENTER POSITION
    xc = np.sum(xx * np.abs(datad)) / np.sum(np.abs(datad))
    yc = np.sum(yy * np.abs(datad)) / np.sum(np.abs(datad))

    # DEFINITION OF PIL FUNCTION
    def PIL(xs=None, ys=None, theta=None, hem=None):

        pil1 = (xx - xc - xs) * np.cos(theta * np.pi / 180.) + (yy - yc - ys) * np.sin(
            theta * np.pi / 180.)

        mask1 = (pil1 > 0) * np.ones_like(pil1)
        mask2 = (pil1 < 0) * np.ones_like(pil1)

        return hem * (mask2 - mask1) * np.abs(datad[:, :])

    # FIRST RUN ---------------------------

    I = np.zeros(shape=(20, 20))
    for i in range(0, 20):
        for j in range(0, 20):
            I[j, i] = np.sum(np.abs(datad - PIL(xs=np.arange(-10, 10, 1)[i], ys=0,
                                                theta=np.arange(0, 180, 9)[j],
                                                hem=hsm)) ** 2)

    ind = np.unravel_index(np.argmin(I, axis=None), I.shape)

    # FIRST RUN PARAMETERS
    xs1 = np.arange(-10, 10, 1)[ind[1]]
    theta1 = np.arange(0, 180, 9)[ind[0]]

    # SECOND RUN (ZOOM IN)---------------------------

    I2 = np.zeros(shape=(20, 20))
    for i in range(0, 20):
        for j in range(0, 20):
            I2[j, i] = np.sum(np.abs(
                datad - PIL(xs=xs1 + np.arange(-1, 1, 0.1)[i], ys=0,
                            theta=theta1 + np.arange(-9, 9, .9)[j], hem=hsm)) ** 2)

    ind2 = np.unravel_index(np.argmin(I2, axis=None), I2.shape)

    # SECOND RUN PARAMETERS
    xs2 = np.arange(-1, 1, .1)[ind2[1]]
    theta2 = np.arange(-9, 9, .9)[ind2[0]]

    return {'x0': xc + xs1 + xs2, 'y0': yc, 'theta0': theta1 + theta2}


