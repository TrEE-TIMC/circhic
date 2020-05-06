import sys
import numpy as np
from scipy import sparse
# from IPython.core.display import clear_output

import matplotlib.pyplot as plt
from matplotlib import colors

#def display(disp = ''):
#    clear_output();
#    print(disp);
#    sys.stdout.flush()
    

def load_HiC(filename, n_bins):
    """Load HiC data
    
    Reads the Hi-C data from a matrix filename.
    
    Parameters
    ----------
    filename : string
        path to the file containing the Hi-C count data.
        
    n_bins : integer
        number of bins of the Hi-C matrix
        
    Returns
    -------
    (n_bins, n_bins) ndarray containing the count data.
    """

    counts = np.zeros((n_bins, n_bins))
    for line in open(filename).readlines():
        content = line.split()

        # ndarray index starts at 0 while the Hi-C pro
        # matrix format starts at 1.
        i, j = int(content[0])-1, int(content[1])-1 
        counts[i, j] = float(content[2])
        counts[j, i] = counts[i, j]

    return counts

def genCircData(HiC,bin_circ=0.5,r_in=0.5,res=0,s_in=800000,s_out=800000,pos0=1):
    """Generate Circular Strip Data
    
    
    Parameters
    ----------
    HiC : ndarray
        Input HiC matrix to circularize with size (N, N)
        
    bin_circ : float <= 1
        The size of the output matrix is (Nc, Nc) where Nc=N/bin_circ
        
    r_in : float <= 1
        Inner radius of the strip, supposing that the outer radius is equal to 1
        
    res : integer
        Resolution of HiC (in bp)
        
    s_in : integer
        Genomic distance corresponding to the inner circle
        
    s_out: integer
        Genomic distance corresponding to the outer circle
        
    pos0: integer
        Genomic position at the vertical top        
        
    Returns
    -------
    (Nc, Nc) ndarray containing the count data projected on a circular strip.
    """
            
    if not res: sys.exit('HiC resolution must be set')    
    if r_in>=1: sys.exit('r_in is normalized wrt 1, i.e. it must be <= 1')    
        
    Lg=res*len(HiC) #genome length given the input resolution
    if pos0>Lg: sys.exit('pos0 must be <= Lg')
    
    N,Nc=int(len(HiC)),int(len(HiC)/bin_circ)
    
    #converting indexes of the output matrix into real values (from 0 to 2) such that the size of the matrix is 2 and, hence, the outer radius of the strip is 1
    X= np.tile(2*np.arange(Nc)/Nc-1,(Nc,1))    

    #corresponding radii (R)
    X2=X*X
    R = np.sqrt(X2+np.transpose(X2))
        
    # indexes of non-zero R entries (to work only with the circular strip entries)
    iR = (R >= r_in) & (R <= 1)
        
    #corresponding polar angles (Theta)
    ThetaX, ThetaY = np.zeros((Nc,Nc)), np.zeros((Nc,Nc))
    ThetaX[iR] = np.arccos(X[iR]/R[iR])
    ThetaY[iR] = np.arcsin(np.transpose(X)[iR]/R[iR])
    
    Theta = np.array(ThetaX)
    
    iT = (ThetaY<=0) & (ThetaX>np.pi/2)
    Theta[iT & iR] = ThetaX[iT & iR] - 2*ThetaY[iT & iR]
    
    iT = (ThetaY<=0) & (ThetaX<=np.pi/2)
    Theta[iT & iR] = 2*np.pi + ThetaY[iT & iR]
        
    shift_theta=(pos0-1)/Lg*2*np.pi    
    Theta[iR] = (Theta[iR] + np.pi/2+shift_theta)%(2*np.pi)            

    #Bref is the bin "x" of the ms where we consider a contact between x-s/2 and x+s/2
    Bref = np.zeros((Nc,Nc),dtype=int)
    Bref[iR] = np.minimum((Theta[iR]/2/np.pi*N).astype(int),N-1) #should work without the minimum... (just a precaution)

    #Half_s is the corresponding "s/2" in x-s/2 and x+s/2
    Half_s = np.zeros((Nc,Nc),dtype=int)
    Half_s[iR] = ((s_out*N/Lg - (R[iR] - 1)/(r_in - 1)*(s_in+s_out)*N/Lg)/2).astype(int)
    
    #building bins x-s/2 and x+s/2
    Ih, Jh = np.zeros((Nc,Nc),dtype=int), np.zeros((Nc,Nc),dtype=int)
    Ih[iR] = (Bref[iR] - Half_s[iR] + N)%N
    Jh[iR] = (Bref[iR] + Half_s[iR] + N)%N
    
    #filling up the final matrix
    C = np.empty((Nc,Nc)) * np.nan #background mask
    for ic,jc in np.argwhere(iR>0): C[ic,jc] = HiC[Ih[ic,jc],Jh[ic,jc]]
    
    return C

def plotCirc_simple(Circ,H):
    """Plot a circular strip
        
    Parameters
    ----------
    Circ : ndarray
        Circular strip matrix
        
    H : ndarray
        corresponding original HiC matrix, to compute vmax and vmin
             
    Returns
    -------
    An image can be used with e.g. plt.subplot
    """
    
    vmax=np.max([H[i,(i+1)%len(H)] for i in range(len(H))])
    vmin=vmax/50
    min_non_zero=np.min(H[H>0])
    
    plt.imshow(Circ,norm=colors.SymLogNorm(min_non_zero),vmin=vmin,vmax=vmax,cmap='viridis')
    plt.xticks([]); plt.yticks([]);
    plt.axis('off')
    
    return

def plotCirc(H,r_in,s_in,s_out,genomic_data,marks={},bin_circ=0.5,res=10000,pos0=1):
    """Plot a circular strip with the option to add genomic data and marks
        
    The plot requires the original HiC matrix (to generate the circular strip matrix) because it makes radii handling more robust 
    
    Parameters
    ----------
    H : ndarray
        original HiC matrix
        
    r_in : float <= 1
        Inner radius of the strip, supposing that the outer radius is equal to 1
        
    s_in : integer
        Genomic distance corresponding to the inner circle
        
    s_out: integer
        Genomic distance corresponding to the outer circle
        
    genomic_data: dictionary
        'data' : input data (per bin of the HiC matrix)
        'r_min' : 'corresponding radius of the minimum value'
        'r_max' : 'corresponding radius of the maximum value'

    marks: dictionary (keys = mark names) of dictionary with:
        'bin' : localtion of the mark
        'marker' : marker type
        'ms' : marker size
        'color' : marker color

    bin_circ : float <= 1
        The size of the output matrix is (Nc, Nc) where Nc=N/bin_circ
        
    res : integer
        Resolution of HiC (in bp)
        
    pos0: integer
        Genomic position at the vertical top        
             
    Returns
    -------
    An image that cannot be used with e.g. plt.subplot
    """
    
    #circular strip
    vmax=np.max([H[i,(i+1)%len(H)] for i in range(len(H))])
    vmin=vmax/50
    min_non_zero=np.min(H[H>0])
    
    Circ=genCircData(H,bin_circ=bin_circ,r_in=r_in,res=res,pos0=pos0,s_in=s_in,s_out=s_out)
    
    #genomic data
    if 'data' in genomic_data:
        dataG, rmin_g, rmax_g = genomic_data['data'], genomic_data['r_min'], genomic_data['r_max']
    else:
        dataG, rmin_g, rmax_g = [], 0, 1
    
    
    #figure with three different axes to handle 1) the circular strip, 2) the genomic data and 3) the genomic marks
    fig = plt.figure()

    #circular strip
    axis_c=np.array([np.max(((1-1/rmax_g)/2,0)), np.max(((1-1/rmax_g)/2,0)), np.min((1/rmax_g,1)), np.min((1/rmax_g,1))])
    ax_c=fig.add_axes(axis_c)
    ax_c.imshow(Circ,norm=colors.SymLogNorm(min_non_zero),vmin=vmin,vmax=vmax,cmap='viridis')
    plt.xticks([]); plt.yticks([]);
    plt.axis('off')

    #genomic data
    if 'data' in genomic_data:
        axis_g=np.array([np.max(((1-rmax_g)/2,0)), np.max(((1-rmax_g)/2,0)), np.min((rmax_g,1)), np.min((rmax_g,1))])
        ax_g = fig.add_axes(axis_g, polar=True, frameon=False)
        theta=np.array([np.pi/2-i*2*np.pi/len(dataG) for i in range(len(dataG))])
        r=(dataG-np.min(dataG))/(np.max(dataG)-np.min(dataG))*(1-rmin_g/rmax_g)+rmin_g/rmax_g
        ax_g.plot(np.concatenate((theta,[theta[0]])),np.concatenate((r,[r[0]])),'r-')
        ax_g.set_rmax(1)
        plt.rgrids([]); plt.thetagrids([]); plt.axis('off');

    #genomic marks
    if marks:
        ax_m = fig.add_axes(axis_c,polar=True,frameon=False)
        for name,mark in marks.items():

            if 'color' not in mark: color='white'
            else: color=mark['color']

            if 'marker' not in mark: marker='o'
            else: marker=mark['marker']

            if 'ms' not in mark: ms=14
            else: ms=mark['ms']

            theta=[np.pi/2-mark['bin']*2*np.pi/len(H)]    
            r=[r_in + s_in*(1 - r_in)/(s_out + s_in)]    
            ax_m.plot(theta,r,marker,ms=ms,color=color)
            ax_m.set_rmax(1)        
        plt.rgrids([]); plt.thetagrids([]); plt.axis('off');

    return
