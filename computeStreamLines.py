#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import time
import sys
import msh
import multiprocessing

#Parallel functions
def fun(f, q_in, q_out):
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))
def parmap(f, X, nprocs=multiprocessing.cpu_count()):
    if nprocs>1:
        q_in = multiprocessing.Queue(1)
        q_out = multiprocessing.Queue()
        proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out))
                for _ in range(nprocs)]
        for p in proc:
            p.daemon = True
            p.start()
        sent = [q_in.put((i, x)) for i, x in enumerate(X)]
        [q_in.put((None, None)) for _ in range(nprocs)]
        res = [q_out.get() for _ in range(len(sent))]
        [p.join() for p in proc]
        return [x for i, x in sorted(res)]
    else:
        return [f(x) for x in X]

#Maths
def isPointInTetra(pt,pts):
    A=np.ones((4,4))
    A[:,:3]=pts
    mats = np.zeros((5,4,4))
    for i in range(5):
        if i == 0:
            mats[0] = np.copy(A)
        else:
            mats[i] = np.copy(A)
            mats[i,i-1,:3]=pt
    dets = np.linalg.det(mats)
    if np.sum(np.sign(dets))==5 or np.sum(np.sign(dets))==-5:
        return True
    else:
        return False
def interpolate(pt, pts, vectors):
    dist = np.linalg.norm(pts-pt,axis=1)
    s = np.sum(1./dist)
    vec = [ np.sum([1./d*v for d,v in zip(dist,vectors[:,j])])/s for j in range(3) ]
    return np.array(vec)
def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.
    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

#Octree class and methods
class Octree:
    def __init__(self, mesh, depth):
        self.time1 = 0
        self.time2 = 0
        self.time3 = 0
        self.depth = depth
        self.mesh = mesh
        self.getOctreeIndices()
        self.distributeTetras()
        self.makeTetraSet()
    def getOctreeIndices(self):
        tmpPts = np.copy(self.mesh.verts[:,:3])
        tmpPts-=[self.mesh.xmin, self.mesh.ymin, self.mesh.zmin]
        tmpPts/=self.mesh.dims
        tmpPts*=2**self.depth
        tmpPts[tmpPts==2**self.depth]=2**self.depth-1
        self.ptsIndices = tmpPts.astype(np.int16, copy=False)
    def distributeTetras(self):
        done = np.zeros((len(self.mesh.verts)),dtype="bool_")
        self.tetsIndices = []
        for i,tet in enumerate(self.mesh.tets):
            ptInds = self.ptsIndices[tet[:4]]
            mi, ma = np.min(ptInds,axis=0), np.max(ptInds,axis=0)
            diff = ma-mi
            nb = np.sum(diff)
            if nb>1:
                if nb==2:
                    ptInds=mi+cartesian([range(diff[0]+1),range(diff[1]+1),range(diff[2]+1)])
                else:
                    ptInds=mi+cartesian([[1,0],[1,0],[1,0]])
            self.tetsIndices.append(ptInds)
        self.tetras = [[[[] for k in range(2**self.depth)] for j in range(2**self.depth)] for j in range(2**self.depth)]
        for i, tet in enumerate(self.mesh.tets):
            for ptIndex in tet[:4]:
                inds = self.tetsIndices[i]
                for x in inds:
                    self.tetras[x[0]][x[1]][x[2]].append(i)
    def makeTetraSet(self):
        for i in range(2**self.depth):
            for j in range(2**self.depth):
                for k in range(2**self.depth):
                    self.tetras[i][j][k] = set(self.tetras[i][j][k])
    def getPtIndices(self,pt):
        ind = np.copy(pt)
        ind-=[self.mesh.xmin, self.mesh.ymin, self.mesh.zmin]
        ind/=self.mesh.dims
        ind*=2**self.depth
        ind-=1e-8
        ind = ind.astype(np.int16, copy=False)
        ind[ind==2**self.depth]=2**self.depth-1
        return ind
    def computeVelocity(self,pt):
        velocity = None
        ptIndex = self.getPtIndices(pt)
        for i in self.tetras[ptIndex[0]][ptIndex[1]][ptIndex[2]]:
            ptInds = self.mesh.tets[i,:4]
            pts = self.mesh.verts[ptInds][:,:3]
            if isPointInTetra(pt, pts):
                velocity = interpolate(pt,pts,self.mesh.vectors[ptInds])
                break
        if velocity is None:
            if pt[0]<self.mesh.xmin or pt[0]>self.mesh.xmax or pt[1]<self.mesh.ymin or pt[1]>self.mesh.ymax or pt[2]<self.mesh.zmin or pt[2]>self.mesh.zmax:
                return None
            else:
                print "Error: velocity is None at point",pt
                print "writing debug.mesh"
                self.saveSubdivision(ptIndex, "debug.mesh", point=pt, wholeDomain=True)
        return velocity
    def computeTrajectory(self, pt, step, maxIt):
        pos = []
        velocity = self.computeVelocity(pt)
        pos.append(pt)
        for it in range(maxIt):
            pt = pt + step*velocity
            pos.append(pt)
            velocity = self.computeVelocity(pt)
            if velocity is None:
                break
        return pos
    def saveSubdivision(self, inds, outfile, point=None, wholeDomain=None):
        #Retaining the tetra of the considered subdivision
        tetraMesh = msh.Mesh(self.mesh.path)
        tetraMesh.readSol()
        for t in self.tetras[inds[0]][inds[1]][inds[2]]:
            tetraMesh.tets[t,-1]=10
        for ind in range(10):
            tetraMesh.removeRef(ind)
        tetraMesh.tris = np.array([])
        tetraMesh.discardUnused()

        #Creating a cube of the subdivision boundaries
        div = self.mesh.dims/8
        mins = [self.mesh.xmin + inds[0]*div[0], self.mesh.ymin + inds[1]*div[1], self.mesh.zmin + inds[2]*div[2]]
        boundaryMesh=msh.Mesh(cube=[mins[0], mins[0]+div[0], mins[1], mins[1]+div[1], mins[2], mins[2]+div[2]])
        boundaryMesh.replaceRef(0,2)
        tetraMesh.fondre(boundaryMesh)

        #Creating a cube around the point
        if point is not None:
            off = 0.01
            point = msh.Mesh(cube=[point[0]-off, point[0]+off, point[1]-off, point[1]+off, point[2]-off, point[2]+off])
            point.replaceRef(0,1)
            tetraMesh.fondre(point)

        #Creating a box of the whole domain
        if wholeDomain is not None:
            domain = msh.Mesh(cube=[self.mesh.xmin, self.mesh.xmax, self.mesh.ymin, self.mesh.ymax, self.mesh.zmin, self.mesh.zmax])
            tetraMesh.fondre(domain)

        #Writing the file
        tetraMesh.write(outfile)
        tetraMesh.writeSol(outfile[:-5]+".sol")

if __name__ == "__main__":
    print "1 - Opening the .mesh file and the corresponding .sol file"
    t = time.time()
    mesh = msh.Mesh("demo/snailbox3d1.mesh")
    mesh.readSol()
    print "Mesh reading in", time.time() - t,"s."

    print "2 - Creating the octree structure"
    t = time.time()
    octree = Octree(mesh, 3)
    print "Octree structure created in", time.time() - t,"s."

    print "3 - Computing the streamlines"
    t = time.time()
    #Parameters
    nPoints = 10
    maxIt   = 100
    step    = 0.02
    #First points
    initialXZ = np.random.normal(loc=0.5,scale=0.15,size=(nPoints,2))
    initialPts = np.insert(initialXZ, 1, 0.02, axis=1)-0.5
    initialPts = initialPts*octree.mesh.dims + octree.mesh.center
    #Actual computing
    positions = parmap(lambda x: octree.computeTrajectory(x, step, maxIt), initialPts, nprocs=8)
    print "Streamlines computed in", time.time()-t
    print "Length of each streamline:", [len(p) for p in positions]
    print"N'accepte pas ce changement inutile!" 
