import numpy as np
import argparse

def get_parser():
    '''
    Get the command line argument parser
    '''

    parser = argparse.ArgumentParser(description="Get the triangles for chosen side binning")
    parser.add_argument("--output_file", type=str, help="Filename", default="triangles.dat")
    parser.add_argument("--rMin", type=float, help="Minimum side in Mpc/h", default=0.)
    parser.add_argument("--rMax", type=float, help="Minimum side in Mpc/h", default=160.)
    parser.add_argument("--binSize", type=float, help="Minimum side in Mpc/h", default=5.)
    parser.add_argument("--all_combinations", help="Return every combination", action="store_true", default=False)

    return parser


def get_triangles(rMin, rMax, deltaR, all_combinations=False):
    """
    Get all the triangles from a minimum to a maximum separation,
    for a chosen bin size. 

    The triangle sides are defined as the bin centers:

    r_i = rMin + (i + 0.5) * deltaR 

    for i from 0 to the total number of bins

    By default triangles are sorted so that:
       
       - r23 >= r13 >=r12
    
    Moreover, only sides combination that form a closed 
    triangles (sum of internal angles==np.pi) are selected.

    To get all combinations use
    all_combination = False
    
    Parameters
    - rMin: Minimum separation
    - rMax: Maximum separation
    - binSize: bin size
    - all_combinations: False, default, only return side combinations that
                          form a triangles in increasing order; 
                          True, return all combinations

    Return
    - the total number of triangles
    - Array containing the triangle sides
    """
    
    nBins = int((rMax-rMin)/deltaR)

    edges = np.linspace(rMin, rMax, nBins+1)
    bins = np.array([(edges[i+1]+edges[i])*0.5 for i in  range(nBins)])

    triangles = []
    for i in range(nBins):
        r12Min, r12Max, r12 = edges[i], edges[i+1], bins[i]
        for j in range(i, nBins):
            r13Min, r13Max, r13 = edges[j], edges[j+1], bins[j]
        
            r23Min, r23Max = np.max([0, r13Min-r12Max]), r12Max+r13Max
            nBins_r23 = int((r23Max-r23Min)/deltaR)
            for k in range(nBins_r23):
                r23 = r23Min+(k+0.5)*deltaR
                
                mu12 = (r12**2+r13**2-r23**2)/(2*r12*r13)
                mu13 = (r12**2+r23**2-r13**2)/(2*r12*r23)
                mu23 = (r23**2+r13**2-r12**2)/(2*r23*r13)

                isOK = np.all(np.abs([mu12, mu13, mu23]) <=1 )
                isIncreasing = ((r23>=r13) & (r23<=rMax))

                if (isIncreasing & isOK) or all_combinations:
                    triangles.append([r12, r13, r23])
    return len(triangles), np.array(triangles)

if __name__ == "__main__":

    # Get command line arguments
    args = get_parser().parse_args()

    # Minimum side in Mpc/h
    rMin = args.rMin
    
    # Maximum side in Mpc/h
    rMax = args.rMax

    # Bin size in Mpc/h
    binSize = args.binSize

    # Return only close triangles
    all_combinations = args.all_combinations

    # Get the output file
    output_file = args.output_file

    # Compute the triangles sides
    nT, triangles = get_triangles(rMin, rMax, binSize, all_combinations)
    
    # Saving file
    np.savetxt(output_file, triangles, header = "# r12[Mpc/h] r13[Mpc/h] r23[Mpc/h]", fmt=["%g"]*3)
    print("I wrote the file %s with %d triangles"%(output_file, nT))
