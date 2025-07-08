1. [Introduction](#intrp)
2. [How to Install](#install)
3. [How to Run](#run)
4. [Interaction During the Animation](#interact)

<h3>Introduction<a id="intro"></a></h3>

An animation of stellar evolution, showing both motion on the HR diagram and changes in the star size and habitable zones.  This animation is based on [Tim Morton's isochrones package](https://isochrones.readthedocs.io/en/latest/) and [Kopperapu et al's definition of the habitable zone](https://iopscience.iop.org/article/10.1088/0004-637X/765/2/131)

Star_evolution_animation plots two kinds of things on the same figure:
<ul>
  <li>The star's ***position*** on the HR diagram with coordinates (effective temperature, intrinsic brightness).</li>
  <li>The star's ***radius*** in arbitrary units that have nothing to do with the HR diagram.  The radii of the habitable zones Solar System classical planet orbits, and large reference circles are shown on the same scale.  Below, this scale is called the Solar System Scale.</li>
</ul>

![startup image](./images/startup_annotated.png)

<h3>How to Install<a id="install"></a></h3>

Requirements: a good scientific python 3.10 or later environment that includes the astropy package.  Installing Anaconda is the easiest way to create this environment.

You need the isochrones package, which can be installed with <code>pip install isochrones</code>

<h3>How to Run<a id="run"></a></h3>

Clone or download this repository.  Then in the repository directory, run by typing the following on the command line:
<code>python star_evolution_animation.py</code>,
which simulates a 1 Solar mass star with Solar metallicity.

There are three startup options:

<code>python star_evolution_animation.py -m M</code> sets the initial stellar mass to M Solar masses (default M=1)

<code>python star_evolution_animation.py -f F</code> sets the initial stellar metallicity to F (default F=0, Solar metallicity)

<code>python star_evolution_animation.py -s S</code> sets the Solar System scale to S (default S=1)

So, for example, <code>python star_evolution_animation.py -m 2 -f 0.25 -s 100</code> simulates a 2-Solar-mass star with metallicity 0.25 zoomed out from the default scale by a factor of 100.

You can see these options by typing <code>python star_evolution_animation.py --help</code>, which prints       
<pre>usage: star_evolution_animation.py [-h] [-m STARMASS] [-f STARMETALLICITY] [-s SSSCALE]

options:
  -h, --help            show this help message and exit
  -m STARMASS, --starMass STARMASS
                        star mass in Solar masses in range [0.1,300]
  -f STARMETALLICITY, --starMetallicity STARMETALLICITY
                        star metallicity relative to the sun in dex in range [-4, 0.5]
  -s SSSCALE, --ssScale SSSCALE
                        Solar System scale
</pre>

Star_evolution_animation uses grids in stellar mass and metallicity from the isochrones package.  If you specify a mass or metallicity that is not on that grid, the mass or metallicity will be rounded up to the next value that is on the grid.

<h3>Interaction During the Animation<a id="interact"></h3>
While the animation is running, you can 
<ul>
  <li>Right-click in the animation window to pause and resume the animation.</li>
  <li>Press the ',' and '.' key to zoom the scale in and out.  This does not work when the animation is paused.</li>
</ul>
