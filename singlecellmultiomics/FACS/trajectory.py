#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy
import numpy as np
import scipy.stats


class trajectory():

	def __init__(self):

		self.trajectory={}
		self.trajectory['progression'] = np.concatenate([np.arange(0,1,1.0/300.0),[1.0]])
		self.trajectory['x'] = []
		self.trajectory['y'] = []


	def set(self,x,y,startpoint,endpoint,iterations=10,sampling=50,downsample=0.06,log10=True, nk=10, outlierNumber=10):
		self.orix = x
		self.oriy = y
		self.iterations = iterations
		self.sampling = sampling
		self.startpoint = startpoint
		self.endpoint = endpoint
		self.downsample = downsample
		self.log10 = log10
		if log10:
			self.startpoint=np.log10(self.startpoint)
			self.endpoint=np.log10(self.endpoint)
		self.nk = nk
		self.outliers=outlierNumber

		self.prepareData()


	def prepareData(self):
		#print('filtering data')
		#f = (self.orix>0)&(self.oriy>0)
		x = self.orix#=np.min(x[x>0])/2
		y = self.oriy#=np.min(y[y>0])/2
		if self.log10:
			x = np.log10(self.orix)
			y = np.log10(self.oriy)
		self.filteredx = x
		self.filteredy = y

		xy = np.concatenate([[x],[y]]).T
		# second downsample (triming)
		ix = np.random.permutation(np.arange(len(x)).astype('int'))[:int(len(x)/2)] #/5
		xy2 = np.concatenate([[x[ix]],[y[ix]]]).T
		for i in range(self.outliers):# len(xy)>int(termination-termination*0.1):
			meanDist = np.percentile(scipy.spatial.distance.cdist(xy,xy2),2.5,axis=1)
			ix = np.argsort(meanDist)
			xy = np.delete(xy, ix[-1],0)
		self.filteredx = xy[:,0]
		self.filteredy = xy[:,1]


	def sampleData(self):
		ix = np.random.permutation(np.arange(len(self.filteredx)).astype('int'))[:int(len(self.filteredx)/2)]#/4
		x=self.filteredx[ix]
		y=self.filteredy[ix]
		xy = np.concatenate([[x],[y]]).T
		return xy

	def buildnetwork(self,xy):
		termination = int(len(xy)*self.downsample*4)
		#print('runing first downsample to %d points'%termination)
		while len(xy)>termination:
			print(len(xy))
			meanDist = np.percentile(scipy.spatial.distance.cdist(xy,xy),1,axis=0)
			ix = np.argsort(meanDist)
			xy = np.delete(xy, np.random.choice(ix[0:10],2),0)
		# move points randomly
		xy=xy+0.000010*np.concatenate([[(np.random.rand(len(xy))-0.5)*np.mean(xy[:,0])],[(np.random.rand(len(xy))-0.5)*np.mean([xy[:,1]])]]).T
		# first downsample
		#print('runing first downsample to %d points'%termination)
		#while len(xy)>termination:
		#	print(len(xy))
		#	meanDist = np.percentile(scipy.spatial.distance.cdist(xy,xy),1,axis=0)
		#	ix = np.argsort(meanDist)
		#	xy = np.delete(xy, ix[0:2],0)
		# move points randomly
		#xy=xy+0.33*np.concatenate([[(np.random.rand(len(xy))-0.5)*np.mean(xy[:,0])],[(np.random.rand(len(xy))-0.5)*np.mean([xy[:,1]])]]).T
		self.xy_pre = xy
		#print('running second downsample')
		# second downsample (triming)
		#xy2 = np.concatenate([[self.orix],[self.oriy]]).T
		#for i in range(self.outilers):# len(xy)>int(termination-termination*0.1):
		#	meanDist = np.percentile(scipy.spatial.distance.cdist(xy,xy2),50,axis=1)
		#	ix = np.argsort(meanDist)
		#	xy = np.delete(xy, ix[-1],0)

		# final random selection
		ix = np.random.permutation(np.arange(len(xy))).astype('int')[:int(len(xy)/2)]
		xy=xy[ix,:]
		print('final size of data is %d'%len(ix))


		# add start and end points
		xy = np.concatenate([xy,[self.startpoint],[self.endpoint]],axis=0)
		print(xy[-3:,:])
		d = scipy.spatial.distance.cdist(xy,xy)
		neighbours = {}
		for i,p in enumerate(xy):
			if i not in neighbours.keys():
				neighbours[i]={}
				neighbours[i]['nk']=[]
			neighbours[i]['pos']=p
			neighbours[i]['type']='node'
			currentdist = d[i,:]
			currentdist[currentdist==0]=np.max(currentdist)
			ix=np.argsort(currentdist)
			ix = ix[ix!=i][:self.nk]
			neighbours[i]['nk']=np.unique(np.concatenate([neighbours[i]['nk'],ix.ravel()])).astype('int')
			for nk in neighbours[i]['nk']:
				if nk not in neighbours.keys():
					neighbours[nk]={}
					neighbours[nk]['nk']=[]
				neighbours[nk]['nk']=np.unique(np.concatenate([neighbours[nk]['nk'],[i]])).astype('int')
		return xy, neighbours

	def samplePaths(self,xy,neighbours):
		#print('start of sampling paths')
		startix=len(xy)-2
		endix=len(xy)-1
		sampledPaths=[]
		while len(sampledPaths)<self.sampling:
			path = [startix]
			i = 0
			while True:
				f = np.array([n not in path for n in neighbours[path[i]]['nk']])
				potential = neighbours[path[i]]['nk'][f.ravel()]
				f=[]
				for p in potential:
					#print neighbours[p]['nk'],path[i]
					if np.sum([n not in path for n in neighbours[p]['nk'] ]) >= (len(neighbours[p]['nk'])-2):
						f.append(True)
					else:
						f.append(False)
				potential=potential[np.array(f).ravel()]
				#print(np.sum(f))
				if not(np.any(potential)):
					break
				path.append(np.random.choice(potential.ravel(),1)[0])
				if endix in path:
					break
				i+=1
			if endix in path:
				#print('success')
				sampledPaths.append(path)
		#print('end of sampling paths')
		return sampledPaths


	def run(self):
		self.xy_history=[]

		for iteration in range(self.iterations):
			#print('Currently running iteration %d.'%(iteration+1))
			#print('Sampling data')

			#print('building neighbours')
			#xy = np.concatenate([[x],[y]]).T
			downsampledxy = self.sampleData()
			self.xy,self.neighbours = self.buildnetwork(downsampledxy)
			self.sampledPaths = self.samplePaths(self.xy,self.neighbours)
			self.xy_history.append(self.xy)


			for j in range(len(self.sampledPaths)):
				good=True
				f = np.random.rand(len(self.sampledPaths[j]))>0.6666
				f[0]=True
				f[-1]=True
				subpath = np.array(self.sampledPaths[j])[f]

				traj = {}
				traj['x']=[]
				traj['y']=[]
				for i,k in enumerate(subpath[:-1]):
					root = (100.0*np.sqrt(np.sum((self.neighbours[k]['pos']-self.neighbours[subpath[i+1]]['pos'])**2)))
					deltax = self.neighbours[subpath[i+1]]['pos'][0]-self.neighbours[k]['pos'][0]
					m = (self.neighbours[subpath[i+1]]['pos'][1]-self.neighbours[k]['pos'][1])/deltax
					b = m*self.neighbours[k]['pos'][0]-self.neighbours[k]['pos'][1]
					cx = np.arange(self.neighbours[k]['pos'][0],self.neighbours[subpath[i+1]]['pos'][0],deltax/root)
					if len(cx)==0:
						print('bad')
						good = False
						break
					cy = cx*m-b
					traj['x'] = np.concatenate([traj['x'],cx])
					traj['y'] = np.concatenate([traj['y'],cy])
				if good:

					progression = np.arange(len(traj['x'])).astype('float')/len(traj['x'])
					ix = np.argmin(scipy.spatial.distance.cdist(np.array([self.trajectory['progression']]).T,np.array([progression]).T),axis=1)
					self.trajectory['x'].append(traj['x'][ix])
					self.trajectory['y'].append(traj['y'][ix])
		self.trajectory['mean']={}
		self.trajectory['mean']['x']=np.mean(self.trajectory['x'],axis=0)
		self.trajectory['mean']['y']=np.mean(self.trajectory['y'],axis=0)
