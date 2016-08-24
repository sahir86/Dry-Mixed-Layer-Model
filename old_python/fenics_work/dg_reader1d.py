import numpy as np
import pylab as pl

"""
Function to animate DG plots from xml files. Inputs are:
Xfile - a file containing the x coordinates of DG pts
Ffile - a file prefix (which will be extended to <Ffile>_<n>.xml)
p - degree of DG elements
pts - number of points to plot within each element (default=10)
"""
def animate(Xfile,Ffile,p,n_dumps,pts=10):
    Ffile_n = Ffile+"_0.xml"
    Xvals, Fvals = read(Xfile,Ffile_n,p,pts)
    pl.ion()
    lines = pl.plot(Xvals.T,Fvals.T,'k')
    pl.axis([0.,1.0,-2.,2.])
    for n in range(1,n_dumps):
        Ffile_n = Ffile+"_"+str(n)+".xml"
        Xvals, Fvals = read(Xfile,Ffile_n,p,pts)
        for i, line in enumerate(lines):                
            line.set_ydata(Fvals[i,:])
        pl.draw()
        pl.savefig('Ffile'+'_'+str(n)+'.png')

"""
Function to read in xml files (assumed to be DG) and return 2 ele x
(pts+1) array for plotting, one with the X vals and one with the
function vals. Both input files need to be from the same DG space.
"""
def read(Xname,Fname,p,pts=10):
    Xfile = file(Xname, 'r')
    Ffile = file(Fname, 'r')
    assert(Xfile.readline().strip()=="<?xml version=\"1.0\"?>")
    assert(Xfile.readline().strip()=="<dolfin xmlns:dolfin=\"http://fenicsproject.org\">")
    assert(Ffile.readline().strip()=="<?xml version=\"1.0\"?>")
    assert(Ffile.readline().strip()=="<dolfin xmlns:dolfin=\"http://fenicsproject.org\">")

    s = Xfile.readline().strip()
    l = len(s)
    assert(s[:21]=="<function_data size=\"")
    NDofs = int(s[21:l-2])
    assert(s[l-2:]=="\">")
    NEle = NDofs/(p+1)
    s = Ffile.readline().strip()
    l = len(s)
    assert(s[:21]=="<function_data size=\"")
    NDofs2 = int(s[21:l-2])
    assert(s[l-2:]=="\">")
    assert(NDofs==NDofs2)

    svals = np.arange(0.,1.0*(pts+1))/pts
    Xvals = np.zeros((NEle,np.size(svals)))
    Fvals = np.zeros((NEle,np.size(svals)))

    dof = -1
    for ele in range(NEle):
        Xele = np.zeros((p+1,))
        Fele = np.zeros((p+1,))
        for iloc in range(p+1):
            dof += 1
            ldof = len(str(dof))
            sX = Xfile.readline().strip().split()
            assert(len(sX)==6)
            assert(sX[0]=="<dof")
            indexs = sX[1]
            assert(indexs[:7]=="index=\"")
            indexs = indexs[7:]
            l = len(indexs)-1
            assert(indexs[l]=="\"")
            indexs = indexs[:l]
            assert(dof==int(indexs))
            l = len(sX[2])
            assert(sX[2][:7]=="value=\"")
            assert(sX[2][l-1]=="\"")
            val = float(sX[2][7:l-1])
            Xele[iloc]=val

            sF = Ffile.readline().strip().split()
            assert(len(sF)==6)
            assert(sF[0]=="<dof")
            indexs = sF[1]
            assert(indexs[:7]=="index=\"")
            indexs = indexs[7:]
            l = len(indexs)-1
            assert(indexs[l]=="\"")
            indexs = indexs[:l]
            assert(dof==int(indexs))
            l = len(sF[2])
            assert(sF[2][:7]=="value=\"")
            assert(sF[2][l-1]=="\"")
            val = float(sF[2][7:l-1])
            Fele[iloc]=val
        pX = np.poly1d(np.polyfit(np.arange(0.,1.0*(p+1))/p,Xele,p))
        pF = np.poly1d(np.polyfit(np.arange(0.,1.0*(p+1))/p,Fele,p))
        Xvals[ele,:] = pX(svals)
        Fvals[ele,:] = pF(svals)
            
    sX=Xfile.readline().strip()
    assert(sX=="</function_data>")
    sF=Ffile.readline().strip()
    assert(sF=="</function_data>")
    sX=Xfile.readline().strip()
    assert(sX=="</dolfin>")
    sF=Ffile.readline().strip()
    assert(sF=="</dolfin>")

    return Xvals, Fvals
