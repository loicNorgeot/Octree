#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import time
import sys
import msh

class aNode:
    def __init__(self, maxDepth, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1, depth=0):
        self.maxDepth = maxDepth
        self.xmin = xmin
        self.ymin = ymin
        self.zmin = zmin
        self.xmax = xmax
        self.ymax = ymax
        self.zmax = zmax
        self.depth = depth
        if self.depth<self.maxDepth:
            self.subdivide()
        else:
            self.nodes = []
        self.tetras = []
    def subdivide(self):
        xmin, xmax, ymin, ymax, zmin, zmax = self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax
        n1 = aNode(self.maxDepth, xmin, xmin+(xmax-xmin)/2., ymin, ymin+(ymax-ymin)/2., zmin, zmin+(zmax-zmin)/2., self.depth+1)
        n2 = aNode(self.maxDepth, xmin, xmin+(xmax-xmin)/2., ymin, ymin+(ymax-ymin)/2., zmin+(zmax-zmin)/2., zmax, self.depth+1)
        n3 = aNode(self.maxDepth, xmin, xmin+(xmax-xmin)/2., ymin+(ymax-ymin)/2., ymax, zmin, zmin+(zmax-zmin)/2., self.depth+1)
        n4 = aNode(self.maxDepth, xmin, xmin+(xmax-xmin)/2., ymin+(ymax-ymin)/2., ymax, zmin+(zmax-zmin)/2., zmax, self.depth+1)
        n5 = aNode(self.maxDepth, xmin+(xmax-xmin)/2., xmax, ymin, ymin+(ymax-ymin)/2., zmin, zmin+(zmax-zmin)/2., self.depth+1)
        n6 = aNode(self.maxDepth, xmin+(xmax-xmin)/2., xmax, ymin, ymin+(ymax-ymin)/2., zmin+(zmax-zmin)/2., zmax, self.depth+1)
        n7 = aNode(self.maxDepth, xmin+(xmax-xmin)/2., xmax, ymin+(ymax-ymin)/2., ymax, zmin, zmin+(zmax-zmin)/2., self.depth+1)
        n8 = aNode(self.maxDepth, xmin+(xmax-xmin)/2., xmax, ymin+(ymax-ymin)/2., ymax, zmin+(zmax-zmin)/2., zmax, self.depth+1)
        self.nodes = [n1, n2, n3, n4, n5, n6, n7, n8]
    def isInDomain(self, pt):
        if pt[0]>=self.xmin and pt[0]<=self.xmax and pt[1]>=self.ymin and pt[1]<=self.ymax and pt[2]>=self.zmin and pt[2]<=self.zmax:
            return True
        else:
            return False
    def attributePoint(self, pt):
        if self.depth != self.maxDepth:
            for n in self.nodes:
                if n.isInDomain(pt):
                    return n.attributePoint(pt)
        else:
            return self
    def caracterize(self):
        print self
        print self.depth
        print "[%.2f %.2f] [%.2f %.2f] [%.2f %.2f]" % (self.xmin,self.xmax, self.ymin, self.ymax, self.zmin, self.zmax)

class node:
    def __init__(self):
        self.tets = []
        self.pts  = []

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
def computeVelocity(pt, octree):
    division = octree.attributePoint(pt)
    velocity = None
    for i in division.tetras:
        inds = mesh.tets[i,:4]
        pts = mesh.verts[inds][:,:3]
        if isPointInTetra(pt, pts):
            velocity = interpolate(pt,pts,mesh.vectors[inds])
            break
    return velocity

def getOctreeIndices(mesh,depth):
    tmpPts = np.copy(mesh.verts[:,:3])
    tmpPts-=[mesh.xmin, mesh.ymin,mesh.zmin]
    tmpPts/=mesh.dims
    tmpPts*=2**depth
    tmpPts[tmpPts==2**depth]=2**depth-1
    tmpPts = tmpPts.astype(np.int16, copy=False)
    return tmpPts
def getPtIndices(pt,depth):
    ind = np.copy(pt)
    ind-=[mesh.xmin, mesh.ymin,mesh.zmin]
    ind/=mesh.dims
    ind*=2**depth
    ind = ind.astype(np.int16, copy=False)
    ind[ind==2**depth]=2**depth-1
    return ind


if __name__ == "__main__":
    print "1 - Opening the .mesh file and the corresponding .sol file"
    mesh = msh.Mesh("demo/snailbox3d1.mesh")
    mesh.readSol()

    print "2 - Creating the octree structure"
    depth = 2#512 div
    octree = aNode(depth, mesh.xmin,mesh.xmax,mesh.ymin,mesh.ymax,mesh.zmin,mesh.zmax)

    t = time.time()
    ptsNodes   = [octree.attributePoint(v) for v in mesh.verts]
    print "ptsNodes:", time.time() - t

    t = time.time()
    indices = getOctreeIndices(mesh,depth)
    print "Indices:", time.time() - t

    t = time.time()
    points = [[[[] for k in range(2**depth)] for j in range(2**depth)] for j in range(2**depth)]
    for i,pt in zip(indices,mesh.verts):
        points[i[0]][i[1]][i[2]].append(pt)
    print "Points:",time.time() - t

    t = time.time()
    tetras = [[[[] for k in range(2**depth)] for j in range(2**depth)] for j in range(2**depth)]
    ok = 0
    pasok = 0
    for indTetra, tetra in enumerate(mesh.tets):
        for ptIndex in tetra:
            i = indices[ptIndex]
            tetras[i[0]][i[1]][i[2]].append(indTetra)
        if np.any(indices[tetra[0]]-indices[tetra[1]]) or np.any(indices[tetra[0]]-indices[tetra[2]]) or np.any(indices[tetra[0]]-indices[tetra[3]]):
            pasok+=1
        else:
            ok+=1
    print "Tetrahedra:", time.time() - t
    print pasok, ok


    t = time.time()
    tetraNodes = [set([ptsNodes[i] for i in tet[:4]]) for tet in mesh.tets]
    print "TetraNodes:", time.time() - t

    t = time.time()
    for i,nodeSet in enumerate(tetraNodes):
        for n in nodeSet:
            n.tetras.append(i)
    print "Adding tets:", time.time() - t
    print len(set(ptsNodes)), len(tetraNodes),len([t for t in tetraNodes if len(t)>1])


    """
    print "3 - Computing the streamlines"
    nPoints = 10
    maxIterations = 1000
    step = 0.001
    initialXZ = np.random.random(size=(nPoints, 2))
    initialPts = np.insert(initialXZ, 1, 0.02, axis=1)-0.5
    initialPts = initialPts*mesh.dims + mesh.center
    positions  = []
    velocities = []
    t0 = time.time()

    for position in initialPts:
        pos = []
        vel = []
        #First information
        pos.append(position)
        velocity = computeVelocity(position, octree)
        if velocity is not None:
            vel.append(np.linalg.norm(velocity))
        else:
            print "Had to stop at", position
            division = octree.attributePoint(position)
            division.caracterize()
            break
        #Iterating along the streamlines
        for it in range(maxIterations):
            #New position
            position += step*velocity
            pos.append(position)
            #New velocity
            velocity = computeVelocity(position, octree)
            if velocity is not None:
                vel.append(np.linalg.norm(velocity))
            else:
                break
        positions.append(pos)
        velocities.append(vel)
    print "Resolution en ", time.time() - t0, "s."
    for p,v in zip(positions,velocities):
        print len(p), len(v)
    """
