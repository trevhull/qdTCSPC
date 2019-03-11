#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 13:16:51 2018

@author: Trevor Hull
www.trevorhull.com
github.com/trevhull

TO DO:
    o fittings for all lifetimes
    o figure out inheritance
    o add iteration script
    o finish importing modules
        o loop_CT
        o create_CPfile
        o plot_CP
    o write documentation for modules
    o fix some shitty modules
        o generalize digitize
        o generalize blinkingLifetime
        o add fitting to all lifetimes
    o write new modules
        o background subtraction
        o converting laser power
    o kwargs and optional shit
    o write meta modules
        o automatic loading, plotting
        x blink
        o BX
        o blink
        o on/off
    o save as .json or hdf5
    o GUI?
    o implement CPA


"""

#import os
import numpy as np
import matplotlib.pyplot as plt
from phconvert import pqreader
import pycorrelate as pyc
from scipy.optimize import curve_fit
from matplotlib.gridspec import GridSpec
from scipy.integrate import simps
#import PySimpleGUI as sg

def monoexfit(x, a, b, e):
    '''
    (a*(np.exp((-(x/b)))))+e
 
    '''
    return (a*(np.exp((-(x/b)))))+e
 
def biexfit(x, a, b, c, d, e):
     return (a*(np.exp((-(x/b)))))+(c*(np.exp((-(x/d)))))+e
 
def triexfit(x, a, b, c, d, e, f, g):
     return (a*(np.exp((-(x/b)))))+(c*(np.exp((-(x/d)))))+(e*(np.exp((-(x/f)))))+g

def truncpowerfunc(x, a, b, c, d):
    return (a*(x**(-b)))*(d*(np.exp(-x/c)))

def powerfunc(x, a, b):
    return a*(x**(-b))

class ptu:

    
    def __init__(self, path, file):
        
        self.path = path
        self.name = file
        

    #def loadptu(self):
        filename = self.path + self.name
        timestamps, detectors, nanotimes, self.meta = pqreader.load_ptu(filename)
        nanotimes = nanotimes[detectors !=127]
        timestamps = timestamps[detectors !=127]
        detectors = detectors[detectors !=127]

        self.cins = int(round(1/self.meta['nanotimes_unit']/self.meta['tags']['TTResult_SyncRate']['value']))
        
        
          #This is here because there's some problem in picoquant's record taking where you get
          #photon counts in your lifetime measurement that are impossibly long, i.e. they are longer
          #then the time between laser pulses. I'm not 100% sure why this happens but I think it
          #has something to do with the tcspc resolution. However, you should be careful about the
          #way python converts integers in cins because sometimes if cins = 3999.999 the int rounding
          #will make it 3999. so make sure cins has round() in it
        self.timestamps = timestamps[nanotimes<self.cins]
        self.detectors = detectors[nanotimes<self.cins]
        self.nanotimes = nanotimes[nanotimes<self.cins]
        
        # Need some units to get to truetime aka T2, timestamps_unit & nanotimes)unit also provided by pqreader
        self.truetime = (((self.timestamps*self.meta['timestamps_unit'])+(self.nanotimes*self.meta['nanotimes_unit'])))
        #longtime is gonna be out cutoff, last value of truetime. might not actually need it for this program ...
        self.longtime = self.truetime[-1]
        '''
        file_dict = {'longtime': longtime,
                'truetime': truetime,
                'timestamps':timestamps,
                'detectors':detectors,
                'nanotimes':nanotimes,
                'meta':meta}
        
        self.file_dict = file_dict
        '''
        return #file_dict
    
    def window(self, lowerBound, upperBound):
        '''
        Create a new ptu object that only contains photons during some duration
        of the experiment between lowerBound and upperBound
        '''
        newself = ptu(self.path, self.name)
        
        newself.timestamps = newself.timestamps[(newself.truetime > lowerBound) & (newself.truetime < upperBound)]
        newself.nanotimes = newself.nanotimes[(newself.truetime > lowerBound) & (newself.truetime < upperBound)]
        newself.detectors = newself.detectors[(newself.truetime > lowerBound) & (newself.truetime < upperBound)]
        newself.truetime = newself.truetime[(newself.truetime > lowerBound) & (newself.truetime < upperBound)]

        newself.longtime = newself.truetime[-1]
        
        print('data has been truncated between ' + str(lowerBound) + ' and ' + str(upperBound) +' seconds')
        
        return newself
    
    def lifetime(self, fitFunc):
        #meta = self.meta
        bins = int(1/self.meta['nanotimes_unit']/self.meta['tags']['TTResult_SyncRate']['value'])
        nanotimes = self.nanotimes
        lifeIntensity, lifeBin_ = np.histogram(nanotimes, bins)
        lifeBins = lifeBin_[:(len(lifeBin_))-1]*self.meta['nanotimes_unit']*1e6
        #fig = plt.figure(figsize=(3,2))
        #ax1 = plt.subplot(1,1,1)
        #ax1.set_yscale('log')
        #plt.plot(bubin, plint)
        #plt.xlabel('time ($\mu$s)')
        #plt.ylabel('intensity (counts)')
        #plt.show()
        
        
        self.lifeBins = lifeBins
        self.lifeIntensity = lifeIntensity
        
        return# plint, bubin
        
    def lifePlot(self, fitFunc, lifeBound):
        '''
        '''
        self.lifetime(fitFunc)
        plt.figure(figsize=(5,4))
        ax1 = plt.subplot(1,1,1)
        ax1.set_yscale('log')
        plt.plot(self.lifeBins, self.lifeIntensity, '.')
        
        if fitFunc == monoexfit:
            self.fitMono(lifeBound)
        elif fitFunc == biexfit:
            self.fitBi(lifeBound)
        elif fitFunc == triexfit:
            self.fitTri(lifeBound)
        else:
            print('not a valid fit, sorry')
        
        plt.plot(self.fitBins, fitFunc(self.fitBins, *self.lifepopt), 'r')
        
        
        
        plt.xlabel('time ($\mu$s)')
        plt.ylabel('intensity (counts)')
        
        
        
        
        
        plt.show()
        #print(popt)
        
        
        return
        
    def calcBlink(self, resolution):
        '''
        '''
        bins = int(self.longtime/resolution)
        blinkY,cins = np.histogram(self.truetime,bins)
        blinkX = cins[:(len(cins)-1)]
        
        freqX, dins = np.histogram(blinkY, 1000)
        freqY = dins[:(len(dins)-1)]
        
        self.blinkX = blinkX
        self.blinkY = blinkY
        self.freqX = freqX
        self.freqY = freqY
        self.blinkRes = resolution
        
                
        return #file_dict
    
    def plotBlink(self, resolution):
        '''
        plot some damn blinking
        '''
        self.calcBlink(resolution)
        plt.figure(figsize=(16,8))
        gs = GridSpec(2, 5)
        # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
        ax1 = plt.subplot(gs[0, :-1])
        ax1.tick_params(labelleft=True)
        plt.xlim(min(self.blinkX),self.longtime)
        plt.xlabel('Time (sec)', fontsize=12)
        plt.ylabel('Intensity (Counts)', fontsize=12)
        plt.plot(self.blinkX, self.blinkY)
        ax2 = plt.subplot(gs[0, -1])
        plt.xlabel('frequency (counts)', fontsize=12)
        plt.plot(self.freqX, self.freqY)
        ax2.tick_params(labelleft=False)
        plt.show()
        return #blinkX, blinkY, freqX, freqY
    '''
    def digitize(self, threshold):
        i = 0
        
        off = np.zeros(len(self.blinkY))
        for i in range(0,len(self.blinkY)):
            if self.blinkY[i] > threshold:
                off[i] = 1
            else:
                off[i] = 0
                
        self.digital = off
        self.threshold = threshold
        return #blink_dict
    '''
    def digitize(self, *threshold):
        i = 0
        
        off = np.empty(len(self.blinkY)) * np.NaN
        
        
        #print(threshold)
        for j in range(0,len(threshold)): 
            
            '''
            if type(threshold[j]) != tuple:
                print('this is not my beautiful tuple')
                fixThreshold = (min(self.blinkY), threshold[j]), (threshold[j], max(self.blinkY))
                for i in range(0,len(self.blinkY)): 
                        if (self.blinkY[i] > fixThreshold[j][0]) & (self.blinkY[i] < fixThreshold[j][1]):
                            off[i] = j
                            #print(j)
                        else:
                            pass

            if type(threshold[j]) != tuple:
                print('this is not my beautiful tuple')
                fixThreshold = (min(self.blinkY), threshold[j]), (threshold[j], max(self.blinkY))
            
            '''
            
            if type(threshold[j]) == tuple:
                print(threshold[j])
                #print(j)
                for i in range(0,len(self.blinkY)): 
                    if (self.blinkY[i] > threshold[j][0]) & (self.blinkY[i] < threshold[j][1]):
                        off[i] = j
                        #print(j)
                    else:
                        pass

            else:
                print('this is not my beautiful tuple')
                #fixThreshold = (min(self.blinkY), threshold[0]), (threshold[0], max(self.blinkY))
                
                for i in range(0,len(self.blinkY)): 
                    if self.blinkY[i] > threshold[j]:
                        off[i] = 1
                    else:
                        off[i] = 0
                
                threshold = ((min(self.blinkY), threshold[0]),(threshold[0], max(self.blinkY)))


        self.digital = off
        self.threshold = threshold

        return #blink_dict
    
    
    
    def countDigital(self, choice):
    #NOTE: YOU NEVER ACTUALLY TYPE IN THIS COMMAND IT ALL GOES TO digitallot
    #first get some info from our file dictionary        
    #then we're gonna make 'timebin' which will hold the on/off time as we're counting, this will end
    #up in an array later on
        
        if max(self.digital) == 1:
            timebin = 0
            n = len(self.blinkX)
            #generalized for on or off, depending on what you specify!
            for i in range(n-1):
                #remember you have to digitize before you can count! 1 is on, 0 is off
                if choice == 'off':
                    num = 0
                elif choice == 'on':
                    num = 1

                else:
                    print('bad choice! choose "on" or "off"')
                

                if (self.digital[i] == num):
                    timebin += 1 
                elif (timebin == 0) & (self.digital[i]!=num):
                    pass
                else:
                    yield timebin
                    #print(timebin*0.01)
                    timebin = 0
        
        elif max(self.digital) == 2:
            timebin = 0
            n = len(self.blinkX)
            #generalized for on or off, depending on what you specify!
            for i in range(n-1):
                #remember you have to digitize before you can count! 1 is on, 0 is off

                
                if choice == 'off':
                    num = 0
                elif choice == 'grey':
                    num = 1
                elif choice == 'on':
                    num = 2
                else:
                    print('bad choice! choose "on" or "off" or maybe "grey"')
                                     
                
                if (self.digital[i] == num):
                    timebin += 1 
                elif (timebin == 0) & (self.digital[i]!=num):
                    pass
                else:
                    yield timebin
                    #print(timebin*0.01)
                    timebin = 0
                    
                
            
        return #file_dict
    
    def digitalPlot(self, choice, fitFunc):
    
        #blinkX = file_dict['blinkX']
        #digital = file_dict['digital']
        #choice = str(choice)
        choiceTime = np.fromiter(self.countDigital(choice), dtype=int)
        correctTime = choiceTime*self.blinkRes
        #digBins =  int(len(/0.08863636363636365)
        digBins =  int(self.longtime/self.blinkRes)
        statY, tatX = np.histogram(correctTime, bins = digBins)
        statX = (tatX[:(len(tatX)-1)])
    
        #to properly weight things we need to remove the zeros
        fixStatY = statY[statY!=0]
        probX = statX[statY!=0]
    
        #let's make an array the right length to save memory for our iteration
        probY = np.zeros(len(fixStatY))
        #So we need to weight each value by the probablity or the time distance (dt)
        #between nearest neighbors a and b
        for i in range(0,len(fixStatY)):
            if i == 0:
                dt = abs(probX[i+1] - probX[i])
            elif probX[i] == probX[-1]:
                dt = abs(probX[i-1] - probX[i])
            else:
                a = abs(probX[i-1] - probX[i])
                b = abs(probX[i+1] - probX[i])
                dt = (a+b)/2
            #then take the value fixStatY divide by the total # of records sum(fixStatY) * 1 /dt
            probY[i]= ((fixStatY[i]/sum(fixStatY))*(1/dt))
            
        
        #let's try and fit the data
        try:
            popt, pcov = curve_fit(fitFunc, probX, probY)        
        except RuntimeError:
            print('fit error')
            pass
        
        ax = plt.subplot(1,1,1)
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.plot(probX, probY, 'o')
        try:
            plt.plot(probX, fitFunc(probX, *popt), 'r-')
        except UnboundLocalError:
            pass
        plt.show()
        
        if fitFunc == powerfunc:
            print('P(t) = ' + str(round(popt[0],2)) + '*tau^(-' + str(round(popt[1],2)) +')')
        elif fitFunc == truncpowerfunc:
            print('P(t) = ' + str(round(popt[0],2)) + 'xtau^(-' + str(round(popt[1],2)) + ')*' + str(round(popt[2],2)) + 'e^(-tau/' + str(round(popt[3],2)) + ')')
            
        if choice == 'on':
            self.onProbX = probX
            self.onProbY = probY
            self.onpop = popt
        else:
            self.offProbX = probX
            self.offProbY = probY
            self.offpop = popt
        
        return #popt, probX, probY, digBins
    
    def calcDig(self, onfitFunc, offfitFunc):
        print('ON')
        self.digitalPlot('on',onfitFunc)
        print('OFF')
        self.digitalPlot('off',offfitFunc)
        return
    
    def onFrac(self):
        self.onFrac = len(self.digital[self.digital==1])/len(self.digital)
        print(round(self.onFrac, 2))
        
        return
        
    def offFrac(self):
        self.offFrac = len(self.digital[self.digital==0])/len(self.digital)
        print(round(self.offFrac, 2))
        
        return
    
    def antibunching(self,samples,timegate):
        
        rep = self.meta['tags']['TTResult_SyncRate']['value']
        maxCorr = round(1/rep*1e6,1)+(round((1/rep*1e6)/2,1))
        l = -maxCorr*1e-6
        m = maxCorr*1e-6
        sp = samples
        p = (m-l)/sp
        lags = np.arange(l,m,p)
        self.antibunchX = (lags[:len(lags)-1])* 1e6
        
        correctNanotimes = self.nanotimes*self.meta['nanotimes_unit']
        a = self.truetime[(self.detectors==0)&(correctNanotimes>(timegate/1e9))]
        b = self.truetime[(self.detectors==1)&(correctNanotimes>(timegate/1e9))]
        
        self.G = pyc.pcorrelate(a, b, lags, 1)
        self.H = pyc.pcorrelate(b, a, lags, 1)
        self.antibunchY = (self.G + self.H)/2
        
    
        return self.antibunchY, self.antibunchX
    
    def plotAB(self, samples, timegate):
        self.antibunching(samples, timegate)
        plt.plot(self.antibunchX, self.antibunchY)
        plt.xlabel('delay time ($\mu$s)')
        plt.ylabel('G[t]')
        plt.show()
        
    def BXratio(self, bound):
        rep = round(1/self.meta['tags']['TTResult_SyncRate']['value']*1e6,1)
        lcbound = 0 - (bound/2)
        rcbound = 0 + (bound/2)
        lrbound = rep - (bound)
        rrbound = rep + (bound)
        
        cent = self.antibunchY[(self.antibunchX > lcbound) & (self.antibunchX < rcbound)]
        xcent = self.antibunchX[(self.antibunchX > lcbound) & (self.antibunchX < rcbound)]
    
        right = self.antibunchY[(self.antibunchX > lrbound) & (self.antibunchX < rrbound)]
        xright = self.antibunchX[(self.antibunchX > lrbound) & (self.antibunchX < rrbound)]
    
        area = simps(cent, xcent)
        rarea = simps(right, xright)
        self.bx = area/rarea
        #print(self.bx)
        return self.bx
        
    def printBX(self, bound):
        self.BXratio(bound)
        print(self.bx)
    
    
    
    
    #COME BACK AND WORK ON THIS
    def loopBX(self, maxGateTime, inc, samples):
        '''
        maxGateTime 
        
        
        
        '''
        i = 0
    
        rep = 1/(self.meta['tags']['TTResult_SyncRate']['value'])*1e6
        #maxCorr = rep+(rep/2)
        bound = rep/4
        
        tg = np.zeros((int(maxGateTime/inc)))
        bxarray = np.zeros((int(maxGateTime/inc)))
        for i in range(int(maxGateTime/inc)):
            I, plags = self.antibunching(samples, (i*inc))
            #bx = self.BXratio(bound)
            bxarray[i] = self.BXratio(bound)
            tg[i] = (i*inc)
            
        
        self.GDT = tg
        self.RTG = bxarray    
        #plt.plot(tg,bxarray)
        
        plt.plot(self.GDT, self.RTG, 'o')
        
        plt.show()
        return #file_dict

    def crosstalk(self,samples, timegate):
        self.antibunching(samples, timegate)
        
        #l,m = antibunching_G(file_dict, samples, timegate)
        #n, o = antibunching_H(file_dict, samples,timegate)
        p = np.append(self.H[self.antibunchX<0],self.G[self.antibunchX>0])
        plt.figure(figsize=(3,3))
        plt.plot(self.antibunchX,p)
        plt.xlabel('correlation time ($\mu$s)')
        plt.ylabel('coincidences')
        plt.show()
        
        return #l,m,n,o,p
    
    def calcFlid(self, fitFunc):
        #nanotimes = file_dict['nanotimes']
        #truetime = file_dict['truetime']
        #blinkX = file_dict['blinkX']
        lifetime_resolution = 1
        
        self.flid = np.zeros(len(self.blinkX))
        i=0
        #banotimes = nanotimes[nanotimes<2500]
        while (i+1)*lifetime_resolution < len(self.blinkX):
            dtrins = self.nanotimes[(self.blinkX[(i*lifetime_resolution)]< self.truetime) & (self.truetime < self.blinkX[(i +1)*lifetime_resolution])]
            #ax = plt.subplot(1,1,1)
            #ax.set_yscale('log')
            #print(i)
            #dbins = len(self.blinkX)
            dbins = int(round(len(self.blinkX)/10))
            #ax = plt.subplot(1,1,1)
            #ax.set_yscale('log')
            #plt.xlim
            #plt.show(plt.hist(dtrins, dbins, histtype = 'step'))
            a, bi = np.histogram(dtrins, dbins)
            bbins = bi[:(len(bi)-1)]*self.meta['nanotimes_unit']*1e6
    
    
            fa = a[bbins > 0.004]
            fbins = bbins[bbins > 0.004]
            try:
                popt, pcov = curve_fit(fitFunc, fbins, fa, bounds = (0,[100, 0.2, 100])) #4.44194979  0.14768228 39.53100168  0.05597095
    
            except RuntimeError:
                 #print("Found an error")
                    
                 pass
            if fitFunc == biexfit:
                self.flid[i] = ((popt[0]*popt[1]) + (popt[2]*popt[3]))/(popt[0]+popt[2])
                #self.flid[i] = tave
            else:
                #tave = popt[1]
                self.flid[i] = popt[1]
            #plt.hist(dtrins, dbins, histtype = 'step')
            #plt.plot(fbins, *popt(fbins))
            #plt.show()
            #self.flid = flid
            i+=1
        return #
    
    def heatFlid(self, xlim):
        #self.flid = file_dict['flid']
        #blink_y = file_dict['blink_y']
        extent = [min(self.flid),xlim,min(self.blinkY),max(self.blinkY)]
        plt.hexbin(self.flid,self.blinkY,extent=extent,gridsize=80,bins='log')
        plt.xlim(0,xlim)
        plt.ylim(min(self.blinkY),max(self.blinkY))
        plt.xlabel('lifetime ($\mu$s)', fontsize=12)
        plt.ylabel('Intensity (Counts)', fontsize=12)
        plt.show()
        return
    
    def blinkLifetime(self):
        '''
        calculate lifetime histograms in a specified intensity range.
        Requires blinking trace and digitization.
        
        TODO:
            o generalization
            o separate calculation from plotting
            o calculate bins instead of name it
        
        
        '''
        #check which bins are 'on' (==1) and 'off' (==0) and assign each photon
        #in those bins the on/off 1/0 of the bin so we can gather up all nanotimes
        #for the lifetime fitting
        self.nanodigital = np.zeros(len(self.nanotimes))
        i = 0
        j = 0
        while self.blinkX[i] < self.blinkX[-1]:
            if self.truetime[j] > self.blinkX[i]:
                    i+=1
            else:
                #if self.digital[i] == 0:
                    #print('low')
                self.nanodigital[j] = self.digital[i]
                j+=1
    
               # if self.digital[i] == 1:
                    #print('high')
                #    self.nanodigital[j] = 1
                 #   j+=1
    
        #Calculate lifetime histograms
        #"ON" histogram
        dbins = 1000
        
        if max(self.nanodigital) == 1:        

            self.onY, onX_needsTrim = np.histogram(self.nanotimes[self.nanodigital == 1], dbins)#='auto')
            self.onX = onX_needsTrim[:(len(onX_needsTrim))-1]*self.meta['nanotimes_unit']*1e6
            
            #"OFF" histogram
            self.offY, offX_needsTrim = np.histogram(self.nanotimes[self.nanodigital == 0], dbins)#='auto')
            self.offX = offX_needsTrim[:(len(offX_needsTrim))-1]*self.meta['nanotimes_unit']*1e6
    
        elif max(self.nanodigital) == 2:        

            self.onY, onX_needsTrim = np.histogram(self.nanotimes[self.nanodigital == 2], dbins)#='auto')
            self.onX = onX_needsTrim[:(len(onX_needsTrim))-1]*self.meta['nanotimes_unit']*1e6
            
            self.greyY, greyX_needsTrim = np.histogram(self.nanotimes[self.nanodigital == 1], dbins)#='auto')
            self.greyX = greyX_needsTrim[:(len(greyX_needsTrim))-1]*self.meta['nanotimes_unit']*1e6
            
            #"OFF" histogram
            self.offY, offX_needsTrim = np.histogram(self.nanotimes[self.nanodigital == 0], dbins)#='auto')
            self.offX = offX_needsTrim[:(len(offX_needsTrim))-1]*self.meta['nanotimes_unit']*1e6
        
        
        
        return #onY, onX_corrected, offY, offX_corrected
        
    def plotBlinkLifetime(self, choice):
        '''
        print the blinking, frequencies, and lifetimes you calculated using BlinkLife
        
        '''
        #plotting:
        self.blinkLifetime()
        
        
        
        plt.figure(figsize=(16,4))
        gs = GridSpec(1, 8)
                # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
        ax1 = plt.subplot(gs[0, :4])
        plt.xlim(0,self.longtime)
        plt.ylim(0,max(self.blinkY + (self.blinkY*0.05)))
        plt.xlabel('time (sec)', fontsize=12)
        plt.ylabel('Intensity (Counts)', fontsize=12)
        #off
        ax1.axhspan(self.threshold[0][0],self.threshold[0][1], color='r', alpha = 0.5)
        #on
        ax1.axhspan(self.threshold[-1][0],self.threshold[-1][1], color='b', alpha = 0.5)
        
        if max(self.nanodigital) == 2:
            ax1.axhspan(self.threshold[-2][0],self.threshold[-2][1], color='g', alpha = 0.5)
        
        plt.plot(self.blinkX, self.blinkY, 'black')
        #plt.plot(self.blinkX[self.blinkY < 60], self.blinkY[self.blinkY < 60], 'r.')
        #plt.axhline(y=60, color='r', linestyle='--')
    
    
        ax2 = plt.subplot(gs[0, 4:5])
        plt.ylim(0,max(self.blinkY + (self.blinkY*0.05)))
        #plt.xlim(0,500)
        plt.xlabel('frequency (counts)', fontsize=12)
        ax2.axhspan(self.threshold[0][0],self.threshold[0][1], color='r', alpha = 0.5)
        #on
        ax2.axhspan(self.threshold[-1][0],self.threshold[-1][1], color='b', alpha = 0.5)
        
        if max(self.nanodigital) == 2:
            ax2.axhspan(self.threshold[-2][0],self.threshold[-2][1], color='g', alpha = 0.5)
        plt.plot(self.freqX, self.freqY, 'black')
        ax2.tick_params(labelleft=False)
    
    
    
        ax3 = plt.subplot(gs[0, -2:])
        plt.xlabel('time ($\mu$s)', fontsize=12)
        plt.ylabel('log intensity (counts)', fontsize=12)
        ax3.set_yscale('log')   
        plt.plot(self.onX, (self.onY-min(self.onY))/max(self.onY), 'b.')#, histtype = 'step')
        plt.plot(self.offX, (self.offY-min(self.offY))/max(self.offY), 'r.')#, histtype = 'step')
        if max(self.nanodigital) == 2:
            plt.plot(self.greyX, (self.greyY-min(self.greyY))/max(self.greyY), 'g.')#, histtype = 'step')

        if choice == 1:
    
            try:
                fitonX = self.onX[self.onX > 0.008]
                fitonY = self.onY[self.onX > 0.008]
                popt, pcov = curve_fit(monoexfit, fitonX, fitonY, p0 =( 0.01, 0.03, 0.001))#, bounds = (0,[1, 0.02, 1]))
                print(popt)
            except RuntimeError:
                print('this is not my beautiful fit!')
    
            plt.plot(fitonX, monoexfit(fitonX, *popt), 'r--')
    
        else:
            pass
            
        plt.show()
        
        return
    
    def fitTri(self, lifeBound):
        try:
            self.fitIntensity = self.lifeIntensity[self.lifeBins > lifeBound]
            self.fitBins = self.lifeBins[self.lifeBins > lifeBound]
            popt, pcov = curve_fit(triexfit, self.fitBins, self.fitIntensity, p0 =( 0.7, 0.070, .3, 0.03, 0.1, 0.01, self.lifeIntensity[-1]))
            #plt.plot(fitBins, fitFunc(fitBins, *popt), '-')
            self.lifepopt = popt
        except RuntimeError:
            print('unable to fit lifetime, sorry')
        return
    
    def fitBi(self, lifeBound):
        try:
            self.fitIntensity = self.lifeIntensity[self.lifeBins > lifeBound]
            self.fitBins = self.lifeBins[self.lifeBins > lifeBound]
            popt, pcov = curve_fit(biexfit, self.fitBins, self.fitIntensity, p0 =( 0.7, 0.070, .3, 0.01, self.lifeIntensity[-1]))
            #plt.plot(fitBins, fitFunc(fitBins, *popt), '-')
            self.lifepopt = popt
        except RuntimeError:
            print('unable to fit lifetime, sorry')
        return
    
    def fitMono(self, lifeBound):
        try:
            self.fitIntensity = self.lifeIntensity[self.lifeBins > lifeBound]
            self.fitBins = self.lifeBins[self.lifeBins > lifeBound]
            popt, pcov = curve_fit(monoexfit, self.fitBins, self.fitIntensity, p0 =( 0.7, 0.070, .3))
            #plt.plot(fitBins, fitFunc(fitBins, *popt), '-')
            self.lifepopt = popt
        except RuntimeError:
            print('unable to fit lifetime, sorry')
        return
            
    def laser(self):
        au = self.meta['tags']['UsrPowerDiode']['value']
        opticalPower = 0.00391 * au
        print('laser power is ' + str(round(opticalPower,3)) + ' uW')
        self.opticalPower = opticalPower
        transmittance405 = 0.83
        fwhm = 250/1e9*100
        #powerDensity = 0.88*((opticalPower/1e6*transmittance405)/(fwhm**2))
        powerDensity = ((2*opticalPower/1e6*transmittance405)/(np.pi*(fwhm/1.18)**2))

        print('laser power density is ' + str(round(powerDensity,2)) + ' mW/cm^2')
        
        return
    

    
    
    
    
    
class blink(ptu):
    
    def __init__(self, path, file, resolution):
        ptu.__init__(self, path, file)
        self.path = path
        self.name = file
        #self = ptu(self.path,self.name)
        #self.loadptu()
        self.plotBlink(resolution)
        return


