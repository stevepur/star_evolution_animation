import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
import matplotlib.animation as animation
from isochrones.mist import MISTEvolutionTrackGrid # requires pip install isochrones


auInSolarRadii = 215.032
lyInAu = 63241.1
largeReferenceThreshold = 0.55

teffColors = pd.read_csv("teffColorTable.txt",
                         delim_whitespace=True, comment='#', dtype={"bitcode":object})
teffColors = teffColors[teffColors.type=="2deg"]


def age_text(age):
    return str("Age: " + f'{int(age):,}' + " Years")
   
def mass_text(mass):
    formatStr = "{:." + str(3) + "f}"    
    return str("Mass: " + formatStr.format(mass) + " Solar Masses")

def radius_text(radius):
    earthToSunRadius = 0.0091577
    formatStr = "{:." + str(3) + "f}"    
    rStr = "Radius: " + formatStr.format(radius) + " Solar Radii" 
    if radius < 0.08:
        rStr = rStr + " = " + formatStr.format(radius/earthToSunRadius) + " Earth Radii"
    return rStr

def star_color(logTeff):
    teff = 10**logTeff
    starColor = [np.interp(teff, teffColors.temp, teffColors.R),
                np.interp(teff, teffColors.temp, teffColors.G),
                np.interp(teff, teffColors.temp, teffColors.B)]
    return starColor

def earthInstellation(radius, logTeff):
    teff = (10**logTeff)/5778
    return radius**2 * teff**4
    
def instellation_text(radius, logTeff):
    instellation = earthInstellation(radius, logTeff)
    hzb = get_hz_boundaries(radius, logTeff)
    formatStr = "{:." + str(2) + "f}"    
    return str("Earth flux: " + formatStr.format(instellation) 
               + " hz: [" + formatStr.format(hzb[0][0]) + "," + formatStr.format(hzb[1][0]) + "] AU")

def get_hz_flux(teff, hzType = "conservative"):
    if np.isscalar(teff):
        teff = np.array([teff])
    Ts = teff - 5780
#     # erratum Kopperapu 2013
#     KoppHzOptIn = 1.7763 + 1.4335e-4*Ts + 3.3954e-9*Ts**2 + -7.6364e-12*Ts**3 + -1.1950e-15*Ts**4
#     KoppHzInRunGreen = 1.0385 + 1.2456e-4*Ts + 1.4612e-8*Ts**2 + -7.6345e-12*Ts**3 + -1.711e-15*Ts**4
#     KoppHzPessIn = 1.0146 + 8.1884e-5*Ts + 1.9394e-9*Ts**2 + -4.3618e-12*Ts**3 + -6.8260e-16*Ts**4
#     KoppHzOut = 0.3507 + 5.9578e-5*Ts + 1.6707e-9*Ts**2 + -3.0058e-12*Ts**3 + -5.1925e-16*Ts**4

    if hzType == "optimistic":
        hzIndices = [0, 3]
    elif hzType == "conservative":
        hzIndices = [1, 2]
    else:
        raise ValueError('Bad catalog name');
    
    
    hzLabels = ["Recent Venus",
                "Runaway Greenhouse",
                "Maximum Greenhouse",
                "Early Mars",
                "Runaway Greenhouse for 5 ME",
                "Runaway Greenhouse for 0.1 ME"]
    seffsun  = [1.776,1.107, 0.356, 0.320, 1.188, 0.99]
    a = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4]
    b = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8]
    c = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12]
    d = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15]

    hz = np.zeros((len(hzIndices), len(Ts)))
    for i in range(len(hzIndices)):
        hz[i,:] = seffsun[hzIndices[i]] + a[hzIndices[i]]*Ts + b[hzIndices[i]]*Ts**2 + c[hzIndices[i]]*Ts**3 + d[hzIndices[i]]*Ts**4
    
#     return KoppHzOptIn, KoppHzOut, KoppHzPessIn, KoppHzInRunGreen
    return hz

def get_hz_boundaries(radius, logTeff, hzType = "conservative"):
    teff = np.clip(10**logTeff, 2600, 7200)
    hzFlux = get_hz_flux(teff, hzType)
    return (radius * (teff/5778)**2)/np.sqrt(np.abs(hzFlux))

class starEvolutionSimulator:
    def __init__(self, starMass = 1, starMetallicity = 0, ssScale = 3):
        self.solarSystemScale = 3/ssScale # 3 is a good default, larger ssScale shows more of the solar system
        
        grid_tracks = MISTEvolutionTrackGrid()
        if grid_tracks.df.index.names[0] is None:
            grid_tracks.df.index.set_names(['initial_feh', 'initial_mass', 'EEP'], inplace=True)
        minLum = -3.4
        
        # stars from the Gaia catalog of nearby stars https://gucds.inaf.it/GCNS/Original/
        # all stars within 100pc
        # 20% subset of bright enough stars with RUWE < 1.2
        gcnsExtract = pd.read_csv("GCNS_subset.csv")
        gcnsLogTeff = gcnsExtract.logTeff
        gcnsLogLum = gcnsExtract.logLum
        
        masses = np.unique(grid_tracks.df.index.get_level_values('initial_mass').values)
        fehs = np.unique(grid_tracks.df.index.get_level_values('initial_feh').values)

        if (starMass < np.min(masses)) | (starMass > np.max(masses)):
            raise ValueError('starMass must be in the range [' + str(np.min(masses)) + ', ' + str(np.max(masses)) + ']')
        if starMass not in masses:
            oldMass = starMass
            starMass = masses[np.searchsorted(masses, starMass, side="left")]
            print(str(oldMass) + " is not an allowed mass, choosing " + str(starMass))
            
        if (starMetallicity < np.min(fehs)) | (starMetallicity > np.max(fehs)):
            raise ValueError('starMetallicity must be in the range [' + str(np.min(fehs)) + ', ' + str(np.max(fehs)) + ']')
        if starMetallicity not in fehs:
            oldMetallicity = starMetallicity
            starMetallicity = fehs[np.searchsorted(fehs, starMetallicity, side="left")]
            print(str(oldMetallicity) + " is not an allowed metallicity, choosing " + str(starMetallicity))
            
        # compute middle age main sequence
        msMassList = np.arange(0.1, 300, 0.1)
        msLogLum = []
        msLogTeff = []
        for m in range(len(msMassList)):
            try:
                star = grid_tracks.df.xs((0.0, np.round(msMassList[m], 1)))
                maxAge = np.max(star.star_age.values)
                minAge = np.min(star.star_age.values)
                midAge = minAge + (maxAge - minAge)/2
                midIdx = np.searchsorted(star.star_age.values, midAge)
                msLogLum.append(star.logL.values[midIdx])
                msLogTeff.append(star.logTeff.values[midIdx])
            except:
                continue
        
        self.star = grid_tracks.df.xs((starMetallicity, starMass))
        
        if starMass > 20:
            maxLum = 7.6
        elif starMass > 10:
            maxLum = 6.6
        else:
            maxLum = 5.6
            
        self.maxAge = np.max(self.star.star_age.values)
        self.minAge = np.min(self.star.star_age.values)
        
        # pre-compute the habitable zone boundaries
        self.hz = np.nan_to_num(get_hz_boundaries(self.star.radius.values, self.star.logTeff.values), nan=500) # conservative HZ
        self.ohz = np.nan_to_num(get_hz_boundaries(self.star.radius.values, self.star.logTeff.values, hzType = "optimistic"), nan=500) # optimistic HZ
        
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(6, 3)
        ax.set_ylim(minLum, maxLum)
        ax.set_facecolor('k')
        
        # plot background stars that show the main sequence
        ax.plot(gcnsLogTeff, gcnsLogLum, '.', color=[.2, .2, .2], markersize = 0.5, alpha=1, zorder=100)
        # plot middle age main sequence
        ax.plot(msLogTeff, msLogLum, color=[.3, .3, .3], alpha=1, zorder=101)
        
        # plot planet orbits
        self.orbitRad = [0.3870993, 0.72336, 1, 1.52371, 5.2029, 9.537, 19.189, 30.0699]
        orbitColor = ['silver', 'w', 'lightblue', 'orange', 'ivory','goldenrod', 'cyan', 'turquoise']
        orbitLw = [1, 1.3, 2.5, 1.3, 1.5, 1.5, 1, 1]
        orbitAlpha = [0.2, 0.3, 0.5, 0.3, 0.3, 0.3, 0.2, 0.2]
        self.orbits = []
        for i in range(len(self.orbitRad)):
            self.orbits.append(ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
                          s=(self.solarSystemScale*self.orbitRad[i]*auInSolarRadii)**2, lw=orbitLw[i], edgecolors=orbitColor[i],
                               facecolors='none', alpha=orbitAlpha[i], zorder=300))
                               
        # plot large scale references
        refInAu = np.array([100, 200, 300, 400, 500, 600, 700, 800])
        self.largeReferenceCircles = []
        self.largeReferenceRad = []
        self.largeReferenceText = []
        textFontSize = 9
        for i in range(len(refInAu)):
            r = refInAu[i]*auInSolarRadii
            self.largeReferenceRad.append(r)
            self.largeReferenceCircles.append(ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
                  s=(self.solarSystemScale*r)**2, lw=1, ls="-", edgecolors='c',
                       facecolors='none', alpha=0.5, zorder=300))
            # self.largeReferenceText.append(ax.text(self.star.logTeff.values[0]-self.solarSystemScale*r,
            #                                           self.star.logL.values[0], str(refInAu[i]) + " AU", fontsize=textFontSize, c='c',
            #                                           alpha=0.5, zorder=300))
        # refInLy = np.arange(0.1, 1, 0.1)
        # for i in range(len(refInLy)):
        #     r = refInLy[i]*lyInAu*auInSolarRadii
        #     self.largeReferenceRad.append(r)
        #     self.largeReferenceCircles.append(ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
        #           s=(self.solarSystemScale*r)**2, lw=1, ls="-", edgecolors='w',
        #                facecolors='none', alpha=0.5, zorder=300))
            # # self.largeReferenceText.append(ax.text(self.star.logTeff.values[0]-self.solarSystemScale*r,
            # #                                           self.star.logL.values[0], str(refInLy[i]) + " ly", fontsize=textFontSize, c='c',
            # #                                           alpha=0.5, zorder=300))
        # earthOrbit = ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
        #                   s=(self.solarSystemScale*auInSolarRadii)**2, lw=3, edgecolors='lightblue',
        #                        facecolors='none', zorder=300)
        
        self.innerHz = ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
                          s=(self.hz[0][0]*self.solarSystemScale*auInSolarRadii)**2, lw=2, edgecolors='none',
                               facecolors='k', alpha=1, zorder=11)
        self.innerHz2 = ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
                          s=(self.hz[0][0]*self.solarSystemScale*auInSolarRadii)**2, lw=2, edgecolors='none',
                               facecolors='lightgreen', alpha=0.2, zorder=12)
        self.outerHz = ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
                          s=(self.hz[1][0]*self.solarSystemScale*auInSolarRadii)**2, lw=2, edgecolors='none',
                               facecolors='yellow', alpha=0.1, zorder=10)
        self.innerOHz = ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
                          s=(self.hz[0][0]*self.solarSystemScale*auInSolarRadii)**2, lw=2, edgecolors='none',
                               facecolors='k', alpha=1, zorder=13)
        self.outerOHz = ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
                          s=(self.hz[1][0]*self.solarSystemScale*auInSolarRadii)**2, lw=2, edgecolors='none',
                               facecolors='lightgreen', alpha=0.2, zorder=9)
        
        self.starLine, = ax.plot(self.star.logTeff.values[0], self.star.logL.values[0], lw=0.5, c='w',
                          alpha=1, zorder=300)
        self.scat = ax.scatter(self.star.logTeff.values[0], self.star.logL.values[0],
                          s=(self.solarSystemScale*self.star.radius.values[0])**2, lw=0.5, edgecolors='k',
                          color=star_color(self.star.logTeff.values[0]), zorder=200)
        self.starLoc, = ax.plot(self.star.logTeff.values[0], self.star.logL.values[0], 'k+', alpha=0.2, zorder=400)
        self.xl = ax.get_xlim()
        self.yl = ax.get_ylim()
        ax.plot([np.log10(5778),np.log10(5778)], self.yl, 'r', lw=0.3, zorder=300)
        ax.plot(self.xl, [0,0], 'r', lw=0.3, zorder=300)
        textFontSize = 18
        self.ageText = ax.text(self.xl[0]-0.02, self.yl[0]+1.3, age_text(self.star.star_age.values[0]), 
                               fontsize=textFontSize, c='w', zorder=300)
        self.massText = ax.text(self.xl[0]-0.02, self.yl[0]+1.0, mass_text(self.star.mass.values[0]), 
                                fontsize=textFontSize, c='w', zorder=300)
        self.radiusText = ax.text(self.xl[0]-0.02, self.yl[0]+0.7, radius_text(self.star.radius.values[0]), 
                                  fontsize=textFontSize, c='w', zorder=300)
        self.instellationText = ax.text(self.xl[0]-0.02, self.yl[0]+0.4, 
                                        instellation_text(self.star.radius.values[0], self.star.logTeff.values[0]), 
                                        fontsize=textFontSize, c='w', zorder=300)
        
        hz = np.nan_to_num(get_hz_boundaries(self.star.radius.values, self.star.logTeff.values), nan=0.1)
        
        self.timeLum = self.yl[0]+0.3 # at bottom
        # self.timeLum = self.yl[1]-0.3 # at top
        ax.plot(self.xl, [self.timeLum, self.timeLum], 'w', lw=2, alpha=0.5, zorder=300)
        self.timePip, = plt.plot(self.xl[0], self.timeLum, 'yd',  markersize=5, zorder=300)
        self.maxAge = np.max(self.star.star_age.values)
        self.minAge = np.min(self.star.star_age.values)
        ax.text(self.xl[0]-0.02, self.timeLum-0.2, f'{int(self.minAge):,}' + " Years", fontsize=10, c='w', ha="left", zorder=300)
        ax.text(self.xl[1]+0.02, self.timeLum-0.2, f'{int(self.maxAge):,}' + " Years", fontsize=10, c='w', ha="right", zorder=300)
        
        lumLegendText = ["1/10 as bright", 
                         "Solar brightness",
                         "10 x brighter", 
                         "100 x brighter", 
                         "1,000 x brighter", 
                         "10,000 x brighter", 
                         "100,000 x brighter", 
                         "1,000,000 x brighter",
                         "10,000,000 x brighter"
                        ]
        lumLegendValues = [-1, 0, 1, 2, 3, 4, 5, 6, 7]
        for i in range(len(lumLegendValues)):
            ax.text(self.xl[0]-0.01, lumLegendValues[i]-0.05, lumLegendText[i], fontsize=10, c='w', zorder=300)
            ax.plot(self.xl, [lumLegendValues[i],lumLegendValues[i]], 'w', lw=0.3, alpha=0.5, zorder=300)
            
        for t in [1000, 10000, 100000]:
            for tp in [3, 5, 10]:
                tp = t*tp
                if tp == 1e6:
                    continue
                ax.text(np.log10(tp), self.yl[1]-0.2, f'{int(tp):,}' + "K", fontsize=10, c='w', ha='center', zorder=300)

        self.animation = animation.FuncAnimation(fig, self.animateStar, frames=self.star.shape[0], interval=100, repeat=False)

        self.paused = False
        fig.canvas.mpl_connect('button_press_event', self.toggle_pause)
        fig.canvas.mpl_connect('key_press_event', self.key_pressed)
        
    def toggle_pause(self, *args, **kwargs):
        if self.paused:
            self.animation.resume()
        else:
            self.animation.pause()
        self.paused = not self.paused
        
    def key_pressed(self, keyData):
        if not self.paused: # don't scale when the animation is paused and the user can't see the result
            if keyData.key == ',':
                self.solarSystemScale = 1.2*self.solarSystemScale
            elif keyData.key == '.':
                self.solarSystemScale = self.solarSystemScale/1.2


    def animateStar(self, frame_num):
        ii = frame_num % self.star.shape[0]
        
#        self.solarSystemScale = 0.75*(self.xl[1] - self.xl[0])/(self.hz[1][ii]*self.solarSystemScale*auInSolarRadii)
        
        # Update the scatter collection, with the new colors, sizes and positions.
        self.scat.set_facecolors(star_color(self.star.logTeff.values[ii]))
        self.scat.set_sizes([(self.solarSystemScale*self.star.radius.values[ii])**2])
        self.scat.set_offsets([self.star.logTeff.values[ii], self.star.logL.values[ii]])
        self.starLoc.set_data([self.star.logTeff.values[ii]], [self.star.logL.values[ii]])
    
        for i in range(len(self.orbits)):
            self.orbits[i].set_offsets([self.star.logTeff.values[ii], self.star.logL.values[ii]])
            self.orbits[i].set_sizes([(self.solarSystemScale*self.orbitRad[i]*auInSolarRadii)**2])
        # earthOrbit.set_offsets([self.star.logTeff.values[ii], self.star.logL.values[ii]])
        if self.solarSystemScale < largeReferenceThreshold:
            for i in range(len(self.largeReferenceCircles)):
                self.largeReferenceCircles[i].set_offsets([self.star.logTeff.values[ii], self.star.logL.values[ii]])
                self.largeReferenceCircles[i].set_sizes([(self.solarSystemScale*self.largeReferenceRad[i])**2])
                # self.largeReferenceText[i].set_x(self.star.logTeff.values[ii]+0.0015*self.solarSystemScale*self.largeReferenceRad[i])
                # # self.largeReferenceText[i].set_x(self.star.logTeff.values[ii])
                # self.largeReferenceText[i].set_y(self.star.logL.values[ii])
                    
    
        self.innerHz.set_offsets([self.star.logTeff.values[ii], self.star.logL.values[ii]])
        self.innerHz2.set_offsets([self.star.logTeff.values[ii], self.star.logL.values[ii]])
        self.outerHz.set_offsets([self.star.logTeff.values[ii], self.star.logL.values[ii]])
        self.innerHz.set_sizes([(self.hz[0][ii]*self.solarSystemScale*auInSolarRadii)**2])
        self.innerHz2.set_sizes([(self.hz[0][ii]*self.solarSystemScale*auInSolarRadii)**2])
        self.outerHz.set_sizes([(self.hz[1][ii]*self.solarSystemScale*auInSolarRadii)**2])    
    
        self.innerOHz.set_offsets([self.star.logTeff.values[ii], self.star.logL.values[ii]])
        self.outerOHz.set_offsets([self.star.logTeff.values[ii], self.star.logL.values[ii]])
        self.innerOHz.set_sizes([(self.ohz[0][ii]*self.solarSystemScale*auInSolarRadii)**2])
        self.outerOHz.set_sizes([(self.ohz[1][ii]*self.solarSystemScale*auInSolarRadii)**2])    
    
        age = self.star.star_age.values[ii]
        self.ageText.set_text(age_text(age))
        timeFraction = (age - self.minAge)/(self.maxAge - self.minAge)
        self.timePip.set_data([self.xl[0] + timeFraction*(self.xl[1] - self.xl[0])], [self.timeLum])
    
        self.massText.set_text(mass_text(self.star.mass.values[ii]))
        self.radiusText.set_text(radius_text(self.star.radius.values[ii]))
        self.instellationText.set_text(instellation_text(self.star.radius.values[ii], self.star.logTeff.values[ii]))
    
        self.starLine.set_data((self.star.logTeff.values[0:ii], self.star.logL.values[0:ii]))
        
        return  
        
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--starMass", type=float, default=1, help="star mass in Solar masses in range [0.1,300]")
parser.add_argument("-f", "--starMetallicity", type=float, default=0, help="star metallicity relative to the sun in dex in range [-4, 0.5]")
parser.add_argument("-s", "--ssScale", type=float, default=1, help="Solar System scale")
args = parser.parse_args()
sim = starEvolutionSimulator(starMass = args.starMass, starMetallicity = args.starMetallicity, ssScale = args.ssScale)

plt.show()
