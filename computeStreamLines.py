#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import time
import sys
import msh

def isPointInTetra(pt,pts):
    A=np.ones((4,4))
    A[:,:3]=pts
    detA=np.linalg.det(A)
    for i in range(4):
        tmp = np.copy(A)
        tmp[i,:3] = pt
        det = np.linalg.det(tmp)
        if np.sign(det)!=np.sign(detA):
            return False
    return True
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

class Octree:
    def __init__(self, mesh, depth):
        self.depth = depth
        self.mesh = mesh
        self.getOctreeIndices()
        self.distributeTetras()
    def getOctreeIndices(self):
        tmpPts = np.copy(self.mesh.verts[:,:3])
        tmpPts-=[self.mesh.xmin, self.mesh.ymin, self.mesh.zmin]
        tmpPts/=self.mesh.dims
        tmpPts*=2**self.depth
        tmpPts[tmpPts==2**self.depth]=2**self.depth-1
        self.indices = tmpPts.astype(np.int16, copy=False)
    def distributeTetras(self):
        for indTetra, tetra in enumerate(self.mesh.tets):
            ptInds = self.indices[tetra[:4]]
            mi, ma = np.min(ptInds,axis=0), np.max(ptInds,axis=0)
            diff = ma-mi
            nb = np.sum(diff)
            if nb>1:
                if nb==2:
                    ptInds=mi+cartesian([range(diff[0]+1),range(diff[1]+1),range(diff[2]+1)])
                else:
                    ptInds=mi+cartesian([[1,0],[1,0],[1,0]])
        self.tetras = [[[[] for k in range(2**self.depth)] for j in range(2**self.depth)] for j in range(2**self.depth)]
        for indTetra, tetra in enumerate(self.mesh.tets):
            for ptIndex in tetra[:4]:
                i = self.indices[ptIndex]
                self.tetras[i[0]][i[1]][i[2]].append(indTetra)
    def getPtIndices(self,pt):
        ind = np.copy(pt)
        ind-=[self.mesh.xmin, self.mesh.ymin, self.mesh.zmin]
        ind/=self.mesh.dims
        ind*=2**self.depth
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
        return velocity

if __name__ == "__main__":
    print "1 - Opening the .mesh file and the corresponding .sol file"
    t = time.time()
    mesh = msh.Mesh("demo/snailbox3d1.mesh")
    mesh.readSol()
    print "Mesh reading in", time.time() - t,"s."

    print "2 - Creating the octree structure"
    t = time.time()
    octree = Octree(mesh, 4)
    print "Octree structure created in", time.time() - t,"s."

    print "3 - Computing the streamlines"
    t = time.time()
    nPoints = 10
    maxIterations = 1000
    step = 0.01
    initialXZ = np.random.random(size=(nPoints, 2))
    initialPts = np.insert(initialXZ, 1, 0.02, axis=1)-0.5
    initialPts = initialPts*octree.mesh.dims + octree.mesh.center
    positions  = []
    for traj, pt in enumerate(initialPts):
        print "Computing trajectory", traj
        pos = []
        pos.append(pt)
        velocity = octree.computeVelocity(pt)
        if velocity == None:
            print "First velocity is none at point",pt
            sys.exit()
        for it in range(maxIterations):
            pt = pt + step*velocity
            pos.append(pt)
            velocity = octree.computeVelocity(pt)
            if velocity == None:
                print "Velocity is none at point",pt
                break
            if it==maxIterations-1:
                print "Reached end of streamline at point",pt
        positions.append(pos)

    print "Streamlines computed in ", time.time() - t, "s."
    print [len(p) for p in positions]
    print positions[0]
