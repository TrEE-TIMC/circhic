import sys
import numpy as np


def generate_borders(data, granularity=0.5, inner_radius=0.5, resolution=0,
                     inner_gdis=800000, outer_gdis=800000, origin=1,
                     chromosome_type='circular', frac_lin=0.7, rotate_lin=0,
                     thick_r=0.005, thick_extreme=0.002):

    if not resolution:
        raise ValueError("Data resolution must be set")
    if inner_radius >= 1:
        sys.exit('inner_radius is normalized wrt 1, i.e. it must be <= 1')

    Lg = resolution*len(data)  # genome length given the input resolution
    if origin > Lg:
        sys.exit('origin must be <= Lg')

    N, Nc = int(len(data)), int(len(data)/granularity)

    # converting indexes of the output matrix into real values (from -1 to 1)
    # such that the size of the matrix is 2 and, hence, the outer radius of
    # the strip is 1

    V = np.tile(2*np.arange(Nc)/Nc-1, (Nc, 1))

    # corresponding radii (R)
    V2 = V*V
    R = np.sqrt(V2+np.transpose(V2))

    # indexes of non-zero R entries just at the edges
    iR_bord = (
        ((R >= inner_radius) &
         (R <= inner_radius+thick_r)) | ((R >= 1-thick_r) & (R <= 1)))

    if chromosome_type == 'circular':
        C = np.empty((Nc, Nc)) * np.nan  # background mask
        for ic, jc in np.argwhere(iR_bord > 0):
            C[ic, jc] = 1
        return C
    elif chromosome_type == 'linear':

        iR = (R >= inner_radius) & (R <= 1)

        # corresponding CLOCKWISE polar angles (Theta); here, we need to
        # consider separately the angles given by acos and asin REM: genomic
        # coordinate goes clockwise
        Theta_cos, Theta_sin = np.zeros((Nc, Nc)), np.zeros((Nc, Nc))
        Theta_cos[iR] = np.arccos(V[iR]/R[iR])
        Theta_sin[iR] = np.arcsin(np.flip(np.transpose(V)[iR]/R[iR]))

        # default: the CLOCKWISE angle is given by the acos when the asin is
        # negative
        Theta = np.array(Theta_cos)

        # if asin is positive, the angle depends on the value of acos
        iT = (Theta_sin >= 0) & (Theta_cos > np.pi/2)
        Theta[iT & iR] = np.pi+Theta_sin[iT & iR]

        iT = (Theta_sin >= 0) & (Theta_cos <= np.pi/2)
        Theta[iT & iR] = 2*np.pi - Theta_sin[iT & iR]

        # setting the origin of the clockwise polar angle at twelve o'clock

        Theta[iR] = (
            (Theta[iR] + np.pi/2+frac_lin*np.pi-rotate_lin/360*2*np.pi) %
            (2*np.pi))

        # Bref is the bin "x" of the ms where we consider a contact between
        # x-s/2 and x+s/2
        Bref = np.zeros((Nc, Nc), dtype=int)

        # with this defition, angles larger than frac_lin*2*pi are irrelevant,
        # which is handled at the very end

        Bref[iR] = (((Theta[iR]/2/np.pi/frac_lin)*N).astype(int)) % N

        # Half_s is the corresponding "s/2" in x-s/2 and x+s/2
        Half_s = np.zeros((Nc,  Nc),  dtype=int)
        # Half_s[iR] = ((outer_gdis*N/Lg - (R[iR] - 1)/(inner_radius -
        # 1)*(-inner_gdis+outer_gdis)*N/Lg)/2).astype(int)
        if inner_gdis <= 0:
            r_mid = (
                (outer_gdis*inner_radius+np.abs(inner_gdis)*1) /
                (outer_gdis+np.abs(inner_gdis)))
            iS = (R <= r_mid)
            Half_s[iR & iS] = (
                ((inner_gdis*N/Lg+(R[iR & iS] - inner_radius) /
                 (inner_radius - r_mid)*(1+inner_gdis*N/Lg))/2).astype(int))
            iS = (R >= r_mid)
            Half_s[iR & iS] = (
                ((2+(R[iR & iS] - r_mid)/(1 - r_mid) *
                 (outer_gdis*N/Lg-2))/2).astype(int))
        else:
            Half_s[iR] = (
                (outer_gdis * N / Lg -
                 (R[iR] - 1) / (inner_radius - 1) *
                 (-inner_gdis+outer_gdis) * N / Lg) / 2).astype(int)

        # building bins x-s/2 and x+s/2
        Ih,  Jh = np.zeros((Nc, Nc), dtype=int),  np.zeros((Nc, Nc), dtype=int)
        Ih[iR] = (Bref[iR] - Half_s[iR] + N) % N
        Jh[iR] = (Bref[iR] + Half_s[iR] + N) % N

        # filling up the final matrix
        C = np.empty((Nc, Nc)) * np.nan  # background mask

        # handling angles larger than frac_lin*2*pi
        # handling non-periodic conditions
        iT = (Theta <= frac_lin*2*np.pi)
        which_indices = np.argwhere(
            (iR_bord > 0) & (iT > 0) &
            (Bref - np.abs(Half_s) > 0) &
            (Bref + np.abs(Half_s) < N))

        for ic, jc in which_indices:
            C[ic, jc] = 1

        thickness = 0.08 * thick_extreme * 2 * Nc / frac_lin

        which_indices = np.argwhere(
            (iR > 0) & (iT > 0) &
            (((Bref - np.abs(Half_s) > 0) &
              (Bref - np.abs(Half_s) <= thickness)) |
             ((Bref + np.abs(Half_s) < N) &
              (Bref + np.abs(Half_s) >= N - thickness))))

        for ic, jc in which_indices:
            C[ic, jc] = 1

        return C

    else:
        sys.exit('Unknown chromosome_type for generating circular data')


def generate_circular_map(data, granularity=0.5, inner_radius=0.5,
                          resolution=None, inner_gdis=800000,
                          outer_gdis=800000,
                          mode="reflect",
                          origin=1, chromosome_type='circular', frac_lin=0.7,
                          rotate_lin=0):
    """Generate circular Strip Data


    Parameters
    ----------
    data : ndarray
        Input data matrix to circularize with size (N, N)

    granularity : float <= 1
        granularityarity of display: the size of the output matrix is (Nc, Nc)
        where Nc=N/granularity; the smaller, the neater but also the longer

    inner_radius : float <= 1
        Inner radius of the strip, supposing that the outer radius is equal to
        1

    resolution : integer
        resolution of data (in bp)

    inner_gdis : integer
        Genomic distance corresponding to the inner circle

    outer_gdis: integer
        Genomic distance corresponding to the outer circle

    origin: integer
        Genomic position at the vertical top

    chromosome_type : string
        'circular': circular chromosome, default
        'lin': linear chromosome

    frac_lin : float <= 1
        fraction of strip filled by linear chromosome

    rotate_lin : float <= 360
        rotation angle (in degree) of the chromosome (the linear chromosome is
        symmetric with respect to the vertical axis by default)

    Returns
    -------
    (Nc, Nc) ndarray containing the count data projected on a circular strip.
    """
    if resolution is None:
        raise ValueError("Data resolution must be set")

    if inner_radius >= 1:
        raise ValueError("Inner radius should be <= 1")

    if mode not in ["reflect", "distant"]:
        raise ValueError(
            "mode %s is unknown. Possible values for mode are "
            "'reflect', and 'distant'.")

    Lg = resolution*len(data)  # genome length given the input resolution
    if origin > Lg:
        raise ValueError('origin must be <= Lg')

    N, Nc = int(len(data)), int(len(data)/granularity)

    # converting indexes of the output matrix into real values (from -1 to 1)
    # such that the size of the matrix is 2 and, hence, the outer radius of
    # the strip is 1
    V = np.tile(2*np.arange(Nc)/Nc-1, (Nc, 1))

    # corresponding radii (R)
    V2 = V*V
    R = np.sqrt(V2+np.transpose(V2))

    # indexes of non-zero R entries (to work only with the circular strip
    # entries)
    iR = (R > 0) & (R >= inner_radius) & (R <= 1)

    # corresponding CLOCKWISE polar angles (Theta); here, we need to consider
    # separately the angles given by acos and asin
    # REM: genomic coordinate goes clockwise
    Theta_cos, Theta_sin = np.zeros((Nc, Nc)), np.zeros((Nc, Nc))
    Theta_cos[iR] = np.arccos(V[iR]/R[iR])
    Theta_sin[iR] = np.arcsin(np.flip(np.transpose(V)[iR]/R[iR]))

    # default: the CLOCKWISE angle is given by the acos when the asin is
    # negative
    Theta = np.array(Theta_cos)

    # if asin is positive, the angle depends on the value of acos
    iT = (Theta_sin >= 0) & (Theta_cos > np.pi/2)
    Theta[iT & iR] = np.pi+Theta_sin[iT & iR]

    iT = (Theta_sin >= 0) & (Theta_cos <= np.pi/2)
    Theta[iT & iR] = 2*np.pi - Theta_sin[iT & iR]

    # setting the origin of the clockwise polar angle at the top vertical
    if chromosome_type == 'circular':
        Theta[iR] = (Theta[iR] + np.pi/2) % (2*np.pi)
    elif chromosome_type == 'linear':
        Theta[iR] = (
            (Theta[iR] + np.pi/2 + frac_lin*np.pi-rotate_lin/360*2*np.pi) %
            (2*np.pi))
    else:
        raise ValueError("chromosome_type should be 'circular' or 'linear'.")

    # Bref is the bin "x" of the ms where we consider a contact between x-s/2
    # and x+s/2
    Bref = np.zeros((Nc, Nc), dtype=int)
    if chromosome_type == 'circular':
        Bref[iR] = ((((origin-1)/Lg+Theta[iR]/2/np.pi)*N).astype(int)) % N
    else:
        # with this defition, angles larger than frac_lin*2*pi are irrelevant,
        # which is handled at the very end
        Bref[iR] = (((Theta[iR]/2/np.pi/frac_lin)*N).astype(int)) % N

    # Half_s is the corresponding "s/2" in x-s/2 and x+s/2
    Half_s = np.zeros((Nc, Nc), dtype=int)

    if mode == "reflect":
        r_mid = (
            (outer_gdis*inner_radius+inner_gdis) /
            (outer_gdis+inner_gdis))
        iS = (R <= r_mid)
        Half_s[iR & iS] = (
            (-inner_gdis*N/Lg+(R[iR & iS] - inner_radius) /
             (inner_radius - r_mid) *
             (1-inner_gdis*N/Lg))/2).astype(int)
        iS = (R >= r_mid)
        Half_s[iR & iS] = (
            (2+(R[iR & iS] - r_mid) /
             (1 - r_mid) *
             (outer_gdis*N/Lg-2))/2).astype(int)
    else:
        Half_s[iR] = (
            ((outer_gdis*N/Lg - (R[iR] - 1) /
             (inner_radius - 1)*(-inner_gdis+outer_gdis)*N/Lg)/2).astype(int))

    # building bins x-s/2 and x+s/2
    Ih, Jh = np.zeros((Nc, Nc), dtype=int), np.zeros((Nc, Nc), dtype=int)
    Ih[iR] = (Bref[iR] - Half_s[iR] + N) % N
    Jh[iR] = (Bref[iR] + Half_s[iR] + N) % N

    # filling up the final matrix
    C = np.empty((Nc, Nc)) * np.nan  # background mask

    if chromosome_type == 'circular':
        for ic, jc in np.argwhere(iR > 0):
            C[ic, jc] = data[Ih[ic, jc], Jh[ic, jc]]
    else:
        # handling angles larger than frac_lin*2*pi
        # handling non-periodic conditions
        which_indices = np.argwhere(
            (iR > 0) &
            (Theta <= frac_lin*2*np.pi) & (Bref - np.abs(Half_s) > 0) &
            (Bref + np.abs(Half_s) < N))
        for ic, jc in which_indices:
            C[ic, jc] = data[Ih[ic, jc], Jh[ic, jc]]

    return C


def _convert_from_x_to_theta(x, lengths, resolution=None):
    """
    Converts from genomic distance to theta.

    Parameters
    ----------
    gdis : array (n, )
        Array of genomic distances to convert

    lengths : array (l, )
        Lengths of the chromosomes

    resolution : integer, optional, default: None
        Resolution
    """
    if resolution is None:
        resolution = 1

    theta = x * 2 * np.pi / lengths.sum() * resolution
    return theta


def convert_xy_to_thetar(coordinates, lengths, resolution=None):
    """
    Convert from (x, y) coordinates to (theta, r)

    Parameters
    ----------
    coordinates : tuple of ndarray
        (x, y) where x and y are in genomics coordinates

    lengths : ndarray (l, )
        the lengths of the chromosome

    resolution : integer, optional, default: None
        Resolution

    Returns
    -------
    (theta, s) tuple of polar coordinates
    """
    if resolution is None:
        resolution = 1

    x, y = coordinates
    if x.shape[0] != y.shape[0]:
        raise ValueError("coordinates should be of the same length")

    if np.any(x > (lengths.sum() / resolution)):
        raise ValueError("value in x larger than expected")

    if np.any(y > (lengths.sum() / resolution)):
        raise ValueError("value in x larger than expected")

    s = y - x
    theta = (x + s / 2) * 2 * np.pi / lengths.sum() * resolution
    return (theta, s)
