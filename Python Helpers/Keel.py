#%%

import Airfoil
import Bezier
import numpy as np

#%%
def genBulbSlices(c1,z1,z2,c3,z3, xx=12, nSlices=10, slope = False, error = 0.07):
    
    '''
    c0,z0 = 810,0
    c1,z1 = 990,20
    z2 = 180
    c3,z3 = 0,226
    I = 1.956e08
    '''
    c0, z0 = 836, 0

    M  = 181.44     #(kg) -> 400 lbs
    density = 11343 #(kg/m3) Lead
    naca_factor = 8221.0/12000.0 # = A/(t*c^2)

    h = z3

    v_intersection = naca_factor * (12/100) * (0.0417951312013333*h**3 + 296.025928*h**2 + 698896.0*h)
    v_total = (M/density) * 1e9 #Volume in mm3
    v_bulb = (v_total + v_intersection)*(1-error)

    I = v_bulb/(naca_factor*(xx/100.0))

    bCurve = Bezier.BezierCurve(
        c0,z0,
        c1,z1,
        z2,
        c3,z3,
        I
    )

    ccurve, zcurve = bCurve.getBezierCurve()

    x_z_func = lambda z : -0.354098*z  

    ztab, ctab = [],[]
    slices = []

    for z in Airfoil.arange(z0,z3,nSlices):
        z_find,i = Airfoil.find_closest(z,zcurve)
        c = ccurve[i]
        xyz = Airfoil.genNACAxyz(xx,c,z)
        if slope:
            dx = x_z_func(z)
            xyz = movexyz(xyz,dx=dx)
        slices.append(xyz)

    return slices

def genBulbFile(folderpath,c1,z1,z2,c3,z3,xx=12,nSlices=10,slope=False, error = 0.07):

    slices = genBulbSlices(c1,z1,z2,c3,z3,xx,nSlices,slope,error)

    csv = ''
    for n in range(len(slices[0][0])):
        row = ''
        for slice in slices:
            slice_csv = ';'.join([ str(slice[0][n]) , str(slice[1][n]) , str(slice[2][n]) ])
            row += slice_csv + '|'
        row = row[:-1]
        csv += row + '\n'

    filename = '\\' +  '_'.join([
        str(xx),
        str(c1),
        str(z1),
        str(z2),
        str(c3),
        str(z3)
    ]) 

    f = open(folderpath + filename + '.csv','w')
    f.write(csv)
    f.close()

def readSlices(filepath):

    f = open(filepath,'r')
    lines = f.readlines()
    nSlices = len(lines[0].split('|'))
    slices = [[] for n in range(nSlices)]
    for line in lines:
        slices_row = line.replace('\n','').split('|')
        for n in range(nSlices):
            xyz = [float(i) for i in slices_row[n].split(';')]
            slices[n].append(xyz)

    return slices

def movexyz(xyz_array,dx=0,dy=0,dz=0):
        x_array = [j + dx for j in xyz_array[0]]
        y_array = [j + dy for j in xyz_array[1]]
        z_array = [j + dz for j in xyz_array[2]]

        return (x_array,y_array,z_array)

def keelSlices():

    keel_sections = [
                 [838 , 0     ],
                 [940 , 305   ],
                 [1054, 305*2 ],
                 [1165, 305*3 ],
                 [1289, 305*4 ]
    ]

    keel_slices = []
    top_chord = keel_sections[-1][0]
    bottom_chord = keel_sections[0][0]
    for c,z in keel_sections:
        keel_naca = Airfoil.genNACAxyz(xx=12,c=c,z=z)
        keel_slices.append(movexyz(keel_naca,dx=(bottom_chord-c)))
    
    return keel_slices



# %%
