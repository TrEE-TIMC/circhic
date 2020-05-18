import sys
import numpy as np


def generate_circular_map(data, bin_circ=0.5, inner_radius=0.5, res=0,
                          inner_gdis=800000, outer_gdis=800000, origin=1,
                          mode='circ', frac_lin=0.7, rotate_lin=0):

    """Generate Circular Strip Data


    Parameters
    ----------
    data : ndarray
        Input data matrix to circularize with size (N, N)

    bin_circ : float <= 1
        The size of the output matrix is (Nc, Nc) where Nc=N/bin_circ

    inner_radius : float <= 1
        Inner radius of the strip, supposing that the outer radius is equal to
        1

    res : integer
        Resolution of data (in bp)

    inner_gdis : integer
        Genomic distance corresponding to the inner circle

    outer_gdis: integer
        Genomic distance corresponding to the outer circle

    origin: integer
        Genomic position at the vertical top

    mode : string
        'circ': circular chromosome, default
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
    if not res:
        raise ValueError("Data resolution must be set")
    if inner_radius >= 1:
        raise ValueError("Inner radius should be <= 1")

    Lg = res*len(data)  # genome length given the input resolution
    if origin > Lg:
        sys.exit('origin must be <= Lg')

    N, Nc = int(len(data)), int(len(data)/bin_circ)

    # converting indexes of the output matrix into real values (from -1 to 1)
    # such that the size of the matrix is 2 and, hence, the outer radius of
    # the strip is 1
    V = np.tile(2*np.arange(Nc)/Nc-1, (Nc, 1))

    # corresponding radii (R)
    V2 = V*V
    R = np.sqrt(V2+np.transpose(V2))

    # indexes of non-zero R entries (to work only with the circular strip
    # entries)
    iR = (R >= inner_radius) & (R <= 1)

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
    if mode == 'circ':
        Theta[iR] = (Theta[iR] + np.pi/2) % (2*np.pi)
    elif mode == 'lin':
        Theta[iR] = (
            (Theta[iR] + np.pi/2 + frac_lin*np.pi-rotate_lin/360*2*np.pi) %
            (2*np.pi))
    else:
        raise ValueError("mode should be 'circ' or 'lin'.")

    # Bref is the bin "x" of the ms where we consider a contact between x-s/2
    # and x+s/2
    Bref = np.zeros((Nc, Nc), dtype=int)
    if mode == 'circ':
        Bref[iR] = ((((origin-1)/Lg+Theta[iR]/2/np.pi)*N).astype(int)) % N
    else:
        # with this defition, angles larger than frac_lin*2*pi are irrelevant,
        # which is handled at the very end
        Bref[iR] = (((Theta[iR]/2/np.pi/frac_lin)*N).astype(int)) % N

    # Half_s is the corresponding "s/2" in x-s/2 and x+s/2
    Half_s = np.zeros((Nc, Nc), dtype=int)

    Half_s[iR] = np.round(
        (outer_gdis*N/Lg - (R[iR] - 1) /
         (inner_radius - 1) * (inner_gdis+outer_gdis) *
         N / Lg) /
        2).astype(int)

    # building bins x-s/2 and x+s/2
    Ih, Jh = np.zeros((Nc, Nc), dtype=int), np.zeros((Nc, Nc), dtype=int)
    Ih[iR] = (Bref[iR] - Half_s[iR] + N) % N
    Jh[iR] = (Bref[iR] + Half_s[iR] + N) % N

    # filling up the final matrix
    C = np.empty((Nc, Nc)) * np.nan  # background mask

    if mode == 'circ':
        for ic, jc in np.argwhere(iR > 0):
            C[ic, jc] = data[Ih[ic, jc], Jh[ic, jc]]
    else:
        # handling angles larger than frac_lin*2*pi
        # handling non-periodic conditions
        which_indices = np.argwhere(
            (iR > 0) &
            (Theta <= frac_lin*2*np.pi) & (Bref - Half_s > 0) &
            (Bref + Half_s < N))
        for ic, jc in which_indices:
            C[ic, jc] = data[Ih[ic, jc], Jh[ic, jc]]

    return C
