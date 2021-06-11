# -*- coding: utf-8 -*-
# 
# University of Illinois Open Source License
# Copyright 2008-2018 Luthey-Schulten Group,
# All rights reserved.
# 
# Developed by: Luthey-Schulten Group
#                           University of Illinois at Urbana-Champaign
#                           http://www.scs.uiuc.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# - Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# - Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# - Neither the names of the Luthey-Schulten Group, University of Illinois at
# Urbana-Champaign, nor the names of its contributors may be used to endorse or
# promote products derived from this Software without specific prior written
# permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.
# 
# Author(s): Joseph R. Peterson
# 
#

import lm
import pyLM

def getReactionString(rct, prd, rate):
	rxnStr = ""
	if isinstance(rct,str):
		if rct == '':
			rct = '∅'
		rct = [rct]
	if isinstance(prd,str):
		if prd == '':
			prd = '∅'
		prd = [prd]
	rxnStr += " + ".join(rct) 
	rxnStr += " ⟶ "
	rxnStr += " + ".join(prd)
	units = "?"
	if len(rct) == 0:
		units = "molecules/s"
	elif len(rct) == 1:
		units = "s⁻¹"
	elif len(rct) == 2:
		units = "molecule⁻¹sec⁻¹"
	return rxnStr, rate, units


def writeTable(columnNames,rows):
	# Table definitions
	cols = len(columnNames)
	headerStyle = '<td style="text-align:%s"><b>%s</b></td>'
	rowStyle    = '<td style="text-align:%s">%s</td>'
	def align(idx):
		if idx == 0:
			return "left"
		return "center"
	s  = ''
	# Write header
	s += '<table>'
	s += '<tr>%s</tr>'%("".join([headerStyle%(align(i), x) for i,x in enumerate(columnNames)]))
	for row in rows:
		s += '<tr>%s</tr>'%("".join([rowStyle%(align(i), x) for i,x in enumerate(row)]))
	# Write footer
	s += '</table>'
	return s
	

def visualizeRDMEInitialConditions(sim):
    '''Visualize the RDME simulation volume inside Jupyter.
    
    Use the '%matplotlib notebook' magic to enable
    interactive picking of points.

    Args:
        sim (RDMESimulation): The simulation to visualize.

    '''
    # Error checking
    if not isinstance(sim, pyLM.RDME.RDMESimulation):
        raise TypeError("Can only visualize RDME simulations.")
    
    # Start widgets
    import ipywidgets as widgets
    from IPython import display
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    if matplotlib.get_backend() != 'nbAgg':
        raise RuntimeError("This widget requires the nbAgg backend. Restart the notebook with `%matplotlib notebook'")
    plt.ion()
    
    # Get a discretized version of the lattice and dimensions
    lattice = None
    try:
        lattice = sim.lattice
    except:
        lattice = sim.getLattice()
    siteLattice = lattice.getSiteLatticeView()
    particleLattice = lattice.getParticleLatticeView()
    nx, ny, nz = lattice.getXSize()-1, lattice.getYSize()-1, lattice.getZSize()-1
    ns = {"X":nx, "Y":ny, "Z":nz}
    
    # Selector and slider
    # Top Section
    label1 = widgets.HTML(value="<h2>Site Data</h2>")
    toggle = widgets.ToggleButtons(options=["X","Y","Z"], description="Slice:", disabled=False)
    slider = widgets.IntSlider(value=0, min=0, max=nx, description="X:", orientation='horizontal', readout=False)
    intField = widgets.IntText(value = 0)
    # Bottom Section
    label2 = widgets.HTML(value="<h2>Particle Data</h2>")
    particleSelector = widgets.SelectMultiple(options=sim.species_id, description="Particles")
    # Info Section
    label3 = widgets.HTML(value="<h2>Info</h2>")
    position = widgets.HTML(value="Position: N/A")
    siteType = widgets.HTML(value="Site Type: N/A")
    particles = widgets.HTML(value="Particles: N/A")
    # Advanced Settings
    siteCmapSelector = widgets.Dropdown(options = sorted(m for m in plt.cm.datad if not m.endswith("_r")),
                                        value = 'Accent',
                                        description = 'Site Type Colormap:')
    particleCmapSelector = widgets.Dropdown(options = sorted(m for m in plt.cm.datad if not m.endswith("_r")),
                                        value = 'Set2',
                                        description = 'Particle Colormap:')
    radiusPicker = widgets.IntText(value=0, description='Picker Radius (pixels):')
    settingsBox = widgets.VBox([radiusPicker,siteCmapSelector,particleCmapSelector])
    settings = widgets.Accordion(children=[settingsBox],selected_index=None)
    settings.set_title(0, 'Advanced Options')
    
    # Create plot and plot change function
    plt.figure(figsize=(12,5))
    ax = plt.gca()
    fig = plt.gcf()
    
    # Show site lattice
    sdata = siteLattice[:,:,slider.value].transpose()
    image = ax.imshow(sdata, 
                      interpolation='none', 
                      vmin=0, vmax=len(sim.siteTypes.items()),
                      cmap=plt.get_cmap(siteCmapSelector.value, len(sim.siteTypes.items()))) # Z, Y, X
    ax.set_xlabel("Z (nm)")
    ax.set_ylabel("Y (nm)")
    ax.set_xticklabels([i * sim.latticeSpacing/1e-9 for i in ax.get_xticks()])
    ax.set_yticklabels([i * sim.latticeSpacing/1e-9 for i in ax.get_yticks()])
    cbar = fig.colorbar(image, ticks=np.arange(len(sim.siteTypes))+0.5)
    cbar.ax.set_yticklabels([x[0] for x in sorted(sim.siteTypes.items(), key=lambda x:x[1])])
    
    # Show particle lattice
    pdata = np.zeros(shape=sdata.shape)
    if len(particleSelector.value) > 0:
        for p in particleSelector.value:
            pidx = sim.species_id.index(p)+1
            for i, j in np.argwhere(particleLattice[0,:,:,slider.value] == pidx)[:,0:2]:
                pdata[j,i] = pidx
    image2 = ax.imshow(np.ma.masked_where(pdata==0, pdata), 
                       interpolation='none', 
                       vmin=0, vmax=len(sim.species_id)+1, 
                       cmap=plt.get_cmap(particleCmapSelector.value, len(sim.species_id)+1))
    cbar2 = fig.colorbar(image2, ticks=np.arange(len(sim.species_id)+1)+0.5)
    cbar2.ax.set_yticklabels([''] + sim.species_id)

    # Update the plot
    def updatePlot():
        sdata = None
        if toggle.value == "X":
            sdata = siteLattice[:,::-1,slider.value].transpose()
            ax.set_xlabel("Z (nm)")
            ax.set_ylabel("Y (nm)")
            image.set_extent([0,nz,0,ny])
            image2.set_extent([0,nz,0,ny])
        elif toggle.value == "Y":
            sdata = siteLattice[:,slider.value,::-1].transpose()
            ax.set_xlabel("Z (nm)")
            ax.set_ylabel("X (nm)")
            image.set_extent([0,nz,0,nx])
            image2.set_extent([0,nz,0,nx])
        elif toggle.value == "Z":
            sdata = siteLattice[slider.value,::-1,:].transpose()
            ax.set_xlabel("X (nm)")
            ax.set_ylabel("Y (nm)")
            image.set_extent([0,nx,0,ny])
            image2.set_extent([0,nx,0,ny])
            
        # Extract particle data
        pdata = np.zeros(shape=sdata.shape)
        if len(particleSelector.value) > 0:
            for p in particleSelector.value:
                pidx = sim.species_id.index(p)+1
                if toggle.value == "X":
                    for i, j in np.argwhere(particleLattice[:,:,:,slider.value] == pidx)[:,1:3]:
                        pdata[ny-j,i] = pidx
                elif toggle.value == "Y":
                    for i, j in np.argwhere(particleLattice[:,:,slider.value,:] == pidx)[:,1:3]:
                        pdata[nx-j,i] = pidx
                elif toggle.value == "Z":
                    for i, j in np.argwhere(particleLattice[:,slider.value,:,:] == pidx)[:,1:3]:
                        pdata[j,nx-i] = pidx
                        
        ax.set_xticklabels([i * sim.latticeSpacing/1e-9 for i in ax.get_xticks()])
        ax.set_yticklabels([i * sim.latticeSpacing/1e-9 for i in ax.get_yticks()])
            
        image.set_data(sdata)            
        image2.set_data(np.ma.masked_where(pdata==0, pdata))
    
    # Define GUI callback functions
    def onDimChange(change):
        slider.description = change["new"] + ":"
        slider.max = ns[change["new"]]
        slider.value = 0
        # Update the matplotlib image
        updatePlot()
        
    def onSliderChange(change):
        # Update value in the int field
        intField.value = change["new"]
        # Update the matplotlib image
        updatePlot()
        
    def onParticleSelector(change):
        # Update the matplotlib image
        updatePlot()
        
    def onTextChange(change):
        # Update value in the int field
        newVal = change["new"]
        if newVal < slider.min:
            newVal = slider.min
        if newVal > slider.max:
            newVal = slider.max
        slider.value = newVal
        intField.value = newVal
        # Update the matplotlib image
        updatePlot()
        
    def changeColormap(change):
        # Update colormaps
        image.set_cmap(plt.get_cmap(siteCmapSelector.value, len(sim.siteTypes.items())))
        image2.set_cmap(plt.get_cmap(particleCmapSelector.value, len(sim.species_id)+1))
        cbar.ax.set_yticklabels([x[0] for x in sorted(sim.siteTypes.items(), key=lambda x:x[1])])
        cbar2.ax.set_yticklabels([''] + sim.species_id)
        # Update the matplot images
        updatePlot()
        
        
    # Register GUI interactions
    toggle.observe(onDimChange, names="value")
    slider.observe(onSliderChange, names="value")
    particleSelector.observe(onParticleSelector, names="value")
    intField.observe(onTextChange, names="value")
    siteCmapSelector.observe(changeColormap, names="value")
    particleCmapSelector.observe(changeColormap, names="value")
    
    # Create GUI elements
    topBox = widgets.VBox([label1,toggle,widgets.HBox([slider,intField])])
    particleBox = widgets.VBox([label2,particleSelector])
    infoBox = widgets.VBox([label3, position, siteType, particles])
    bottomBox = widgets.HBox([particleBox, infoBox])
    totalBox = widgets.VBox([topBox,bottomBox,settings])
    
    # Interaction with matplotlib interactive
    def onRelease(event):
        # Update the position that was clicked
        position.value = "Position: %0.1fnm, %0.1fnm (%d, %d)"%(event.xdata*sim.latticeSpacing/1e-9, event.ydata*sim.latticeSpacing/1e-9, event.xdata, event.ydata)
        X, Y = int(event.xdata), int(event.ydata)
        r = radiusPicker.value
        
        # Update the info for the site type that was clicked
        localsd = None
        thest = None
        if toggle.value == "X":
            localsd = siteLattice[:,:,slider.value].transpose()
            thest = [x[0] for x in sorted(sim.siteTypes.items(), key=lambda x:x[1])][localsd[Y,X]]
        elif toggle.value == "Y":
            localsd = siteLattice[:,slider.value,:].transpose()
            thest = [x[0] for x in sorted(sim.siteTypes.items(), key=lambda x:x[1])][localsd[Y,X]]
        elif toggle.value == "Z":
            localsd = siteLattice[slider.value,:,:].transpose()
            thest = [x[0] for x in sorted(sim.siteTypes.items(), key=lambda x:x[1])][localsd[X,Y]]
        siteType.value="Site Type: %s"%(thest)

        # Update the particles at the clicked location
        pmap = {}
        for partName, partID in sim.particleMap.items(): 
            particleLattice = sim.lattice.getParticleLatticeView() 
            if toggle.value == "X":
                val = np.sum(particleLattice[:,max(0,X-r):min(X+r,nz),max(0,Y-r):min(Y+r,ny),slider.value,:] == partID)
                if val > 0:
                    pmap[partName] = val
            elif toggle.value == "Y":
                val = np.sum(particleLattice[:,max(0,X-r):min(X+r,nz),slider.value,max(0,Y-r):min(Y+r,nx),:] == partID)
                if val > 0:
                    pmap[partName] = val
            elif toggle.value == "Z":
                val = np.sum(particleLattice[:,slider.value,max(0,Y-r):min(Y+r,ny),max(0,X-r):min(X+r,nx),:] == partID)
                if val > 0:
                    pmap[partName] = val
        if len(pmap) == 0:
            particles.value = "Particles: N/A"
        else:
            particles.value = "Particles: %s"%(",".join(["%dx%s"%(count, part) for part, count in pmap.items()]))
            
    # Register interactive commands
    cid_up = fig.canvas.mpl_connect('button_release_event', onRelease)
    
    # Execute widget
    display.display(totalBox)

