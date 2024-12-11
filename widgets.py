## QT Imports
from queue import Empty
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm

font = {'size'   : 2}
matplotlib.rc('font', **font)

import os
import numpy as np
from scintools.dynspec import Dynspec,BasicDyn
import scintools.ththmod as thth
from astropy.io import fits
import astropy.units as u

import pickle as pkl
import sys
from fitsio import *

   
def do_svd(ds):
    goodf = np.isfinite(np.nanstd(ds,1))
    goodt = np.isfinite(np.nanstd(ds,0))
    dsSVD = np.copy(ds)
    ds_f = dsSVD[goodf]
    ds_ft = ds_f[:,goodt]
    ds_ft/=thth.svd_model(ds_ft).real
    ds_f*=0
    dsSVD*=0
    ds_f[:,goodt]=ds_ft
    dsSVD[goodf]=ds_f
    return(dsSVD)

class MplCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=300):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(MplCanvas, self).__init__(self.fig)
        # policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        # policy.setHeightForWidth(True)
        # self.setSizePolicy(policy)
        self.aspect=float(height)/float(width)
    def heightForWidth(self, width):
        return width*self.aspect

    def save(self,directory='./',filename=''):
        plotFile, filter = QFileDialog.getSaveFileName(self,'Save File',os.path.join(directory,filename))
        self.fig.savefig(plotFile)

class MyTableWidget(QWidget):
    
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.dyn = DynWidget(self)
        self.layout = QVBoxLayout(self)
        self.lastIndex = 0
        
        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tab1 = GeneralTab(self)
        self.tab2 = ThThTab(self,self.tab1)
        self.tab3 = PhaseRetrievalTab(self)


        self.dyn.newLoaded.connect(self.on_dyn_load)
        self.tab1.emptyFolder.connect(self.on_folder_empty)
        
        # Add tabs
        self.tabs.addTab(self.tab1,"General")
        self.tabs.addTab(self.tab2,"Theta-Theta Prep")
        self.tabs.addTab(self.tab3,"Phase Retrieval")
        self.tabs.setTabEnabled(1,False)
        self.tabs.setTabEnabled(2,False)
        
        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        self.tabs.currentChanged.connect(self.changed)
        self.create_display_layout()
        self.tab1.displayLayout.addLayout(self.displayHolder)
        self.tab2.parametersSet.connect(self.thetatheta_parameters_set)

        self.tab1.set_folder(os.getcwd())
        self.tab2.saveParams.connect(self.tab1.save_push)

    def thetatheta_parameters_set(self):
        self.tabs.setTabEnabled(2,True)

    def on_dyn_load(self):
        self.tab1.on_dyn_load()
        self.tab2.on_dyn_load()
        self.tab3.on_dyn_load()
        self.clear_plots()
        self.update_plots()
        if self.tabs.currentWidget() == self.tab2:
            self.tab2.update_overlay()
        self.fdSlider.setValue(100)
        self.tauSlider.setValue(100)
        self.tabs.setTabEnabled(2,False)
        self.fdSlider.setEnabled(True)
        self.tauSlider.setEnabled(True)

    def create_display_layout(self):
        self.displayHolder = QVBoxLayout(self)
        self.displayDspec = MplCanvas()
        self.displaySspec = MplCanvas()
        self.dspecPlot = self.displayDspec.fig.add_subplot(111)
        self.sspecPlot = self.displaySspec.fig.add_subplot(111)
        self.dspecLines = []
        self.sspecLines = []
        splitter = QHBoxLayout()
        # splitter.setAlignment(Qt.AlignHCenter)
        splitter.addWidget(self.displayDspec)
        splitter.addWidget(self.displaySspec)
        self.displayHolder.addLayout(splitter)

        tauHolder = QHBoxLayout()
        self.tauSlider = QSlider(Qt.Orientation.Horizontal)
        self.tauSlider.setMinimum(1)
        self.tauSlider.setMaximum(100)
        self.tauSlider.setValue(100)
        tauLabel = QLabel()
        # tauLabel.setAlignment(Qt.AlignHCenter)
        tauLabel.setText('\u03C4 Extent')
        tauHolder.addWidget(tauLabel)
        tauHolder.addWidget(self.tauSlider)
        self.displayHolder.addLayout(tauHolder)

        fdHolder = QHBoxLayout()
        self.fdSlider = QSlider(Qt.Orientation.Horizontal)
        self.fdSlider.setMinimum(1)
        self.fdSlider.setMaximum(100)
        self.fdSlider.setValue(100)
        fdLabel = QLabel()
        # fdLabel.setAlignment(Qt.AlignHCenter)
        fdLabel.setText('fD Extent')
        fdHolder.addWidget(fdLabel)
        fdHolder.addWidget(self.fdSlider)
        self.displayHolder.addLayout(fdHolder)

        self.fdSlider.sliderMoved.connect(self.fd_scale)
        self.tauSlider.sliderMoved.connect(self.tau_scale)

        self.fdSlider.setEnabled(False)
        self.tauSlider.setEnabled(False)

    def fd_scale(self):
        fd = thth.fft_axis(self.dyn.time,u.mHz)
        extent = thth.ext_find(fd,fd)
        mx= extent[1]
        fraction = self.fdSlider.value()/100.
        self.sspecPlot.set_xlim((-fraction*mx,fraction*mx))
        self.displaySspec.draw()

    def tau_scale(self):
        tau = thth.fft_axis(self.dyn.freq,u.us)
        extent = thth.ext_find(tau,tau)
        mx= extent[1]
        fraction = self.tauSlider.value()/100.
        self.sspecPlot.set_ylim((0,fraction*mx))
        self.displaySspec.draw()

    def clear_plots(self):
        self.dspecPlot.cla()
        self.sspecPlot.cla()
        self.remove_overlay()

    def update_plots(self):
        dyn,dsExtent,ss,ssExtent = self.dyn.basic_plot()

        self.dspecPlot.imshow(dyn,origin='lower',aspect='auto',extent=dsExtent,interpolation='nearest',cmap='magma')
        self.dspecPlot.set_xlabel(r'$t~\left(\rm{s}\right)$')
        self.dspecPlot.set_ylabel(r'$\nu~\left(\rm{MHz}\right)$')
        
        self.sspecPlot.imshow(ss,origin='lower',aspect='auto',extent=ssExtent,norm=LogNorm(vmin=np.median(ss)),interpolation='nearest',cmap='magma')
        self.sspecPlot.set_xlabel(r'$f_D~\left(\rm{mHz}\right)$')
        self.sspecPlot.set_ylabel(r'$\tau~\left(\mu\rm{s}\right)$')
        self.sspecPlot.set_ylim((0,ssExtent[3]))
        # self.displayFigure.fig.suptitle(self.dyn.filename)
        self.displaySspec.fig.tight_layout()
        self.displayDspec.fig.tight_layout()
        try:
            self.displaySspec.draw()
        except:
            self.sspecPlot.cla()
            self.displaySspec.draw()
        try:
            self.displayDspec.draw()
        except:
            self.dspecPlot.cla()
            self.displayDspec.draw()

    def add_overlay(self,ncf,cwf,nct,cwt,edges_lim,tau_mask,eta_min,eta_max):
        upperCuts = self.dyn.dyn.freqs[np.linspace(1,ncf,ncf).astype(int)*cwf-1] + self.dyn.dyn.df/2
        rightCuts = self.dyn.dyn.times[np.linspace(1,nct,nct).astype(int)*cwt-1] + self.dyn.dyn.dt/2
        for cut in upperCuts:
            self.dspecLines.append(self.dspecPlot.axhline(cut,color='r',lw=1))
        for cut in rightCuts:
            self.dspecLines.append(self.dspecPlot.axvline(cut,color='r',lw=1))
        self.sspecLines.append(self.sspecPlot.axvline(edges_lim,color='r',lw=1))
        self.sspecLines.append(self.sspecPlot.axvline(-edges_lim,color='r',lw=1))
        if tau_mask>0:
            self.sspecLines.append(self.sspecPlot.axhline(tau_mask,color='r',lw=1))
        fd=thth.fft_axis(self.dyn.time,u.mHz)
        self.sspecLines.append(self.sspecPlot.plot(fd,eta_min*fd**2,lw=1,c='tab:blue')[0])
        self.sspecLines.append(self.sspecPlot.plot(fd,eta_max*fd**2,lw=1,c='tab:orange')[0])
        self.displaySspec.fig.tight_layout()
        self.displayDspec.fig.tight_layout()
        self.displaySspec.draw()
        self.displayDspec.draw()

    def remove_overlay(self):
        for line in self.dspecLines:
            line.remove()
        for line in self.sspecLines:
            line.remove()
        self.sspecLines = []
        self.dspecLines = []

    def changed(self,idx):
        if self.tabs.widget(idx) == self.tab2:
            self.tab2.pause = True
            if self.tab2.fileSelect.currentText != self.dyn.filename:
                if self.tab2.fileSelect.findText(self.dyn.filename) == -1:
                    self.tab2.filter_from_dyn()
                    self.tab2.fileSelect.setEditable(True)
                    self.tab2.fileSelect.setCurrentText(self.dyn.filename)
                    self.tab2.fileSelect.setEditable(False)
                else:
                    self.tab2.fileSelect.setCurrentText(self.dyn.filename)
            self.tab2.pause=False
            self.tab2.update_items()
            self.displayHolder.setParent(None)
            self.tab2.displayLayout.addLayout(self.displayHolder)
            self.tab2.update_overlay()
            self.displaySspec.fig.tight_layout()
            self.displayDspec.fig.tight_layout()
            self.displaySspec.draw()
            self.displayDspec.draw()
            self.tab2.on_open()
        if self.tabs.widget(idx) == self.tab1:
            if self.tab1.fileList[self.tab1.currentID] != self.tab2.fileSelect.currentText():
                self.tab1.currentID = np.argwhere(self.tab1.fileList == self.tab2.fileSelect.currentText()).min()
                self.tab1.set_file()
            self.displayHolder.setParent(None)    
            self.tab1.displayLayout.addLayout(self.displayHolder)
            self.displaySspec.fig.tight_layout()
            self.displayDspec.fig.tight_layout()
            self.remove_overlay()
        self.lastIndex = idx

    def on_folder_empty(self):
        self.tabs.setTabEnabled(1,False)
        self.tabs.setTabEnabled(2,False)
        self.dyn.set = False
        self.fdSlider.setEnabled(False)
        self.tauSlider.setEnabled(False)

class GeneralTab(QWidget):
    emptyFolder = pyqtSignal()

    def __init__(self,parent) -> None:
        super().__init__()
        self.layout = QVBoxLayout(parent)
        self.create_layout()
        self.parent = parent
        self.currentID = 0
        
    def create_layout(self):
        self.loadButton = QPushButton("Select Directory")
        self.loadButton.clicked.connect(self.loadpush)
        self.layout.addWidget(self.loadButton)
        self.directoryLabel = QLabel()
        self.directoryLabel.setAlignment(Qt.AlignHCenter)
        self.layout.addWidget(self.directoryLabel)
        self.create_info_layout()
        self.create_control_layout()
        self.displayLayout = QHBoxLayout()
        self.layout.addLayout(self.displayLayout)
        self.setLayout(self.layout)
        self.params = {}
        self.loading = False
        
    def create_info_layout(self):
        self.infoLayout = QHBoxLayout()
        leftLayout = QVBoxLayout()
        rightLayout = QVBoxLayout()
        self.infoLayout.addLayout(leftLayout)
        self.infoLayout.addLayout(rightLayout)
        self.pulsarLabel = QLabel()
        self.pulsarLabel.setAlignment(Qt.AlignHCenter)
        leftLayout.addWidget(self.pulsarLabel)
        self.fileCountLabel = QLabel()
        self.fileCountLabel.setAlignment(Qt.AlignHCenter)
        rightLayout.addWidget(self.fileCountLabel)
        self.layout.addLayout(self.infoLayout)

    def set_labels(self):
        self.directoryLabel.setText(f"Directory: {self.folder}")
        if self.pulsars.shape[0]==0:
            self.pulsarLabel.setText('Pulsars : --')
        elif self.pulsars.shape[0]==1:
            self.pulsarLabel.setText(f"Pulsar: {self.pulsars[0]}")
        else:
            self.pulsarLabel.setText(f"Pulsars: {self.pulsars}")
        self.fileCountLabel.setText(f'Files: {self.fileList.shape[0]}')

    def create_control_layout(self):
        self.flagsLayout = QHBoxLayout()
        self.flagsLayout.setAlignment(Qt.AlignHCenter)
        self.goodCheck = QCheckBox('Good')
        self.interestingCheck = QCheckBox('Interesting')
        self.svdCheck = QCheckBox('SVD')
        self.workCheck = QCheckBox('Needs Work')
        self.flagsLayout.addWidget(self.goodCheck)
        self.flagsLayout.addWidget(self.interestingCheck)
        self.flagsLayout.addWidget(self.svdCheck)
        self.flagsLayout.addWidget(self.workCheck)
        self.layout.addLayout(self.flagsLayout)
        self.buttonsLayout = QHBoxLayout()
        self.buttonsLayout.setAlignment(Qt.AlignHCenter)
        self.nextButton = QPushButton("Next")
        self.prevButton = QPushButton("Previous")
        self.saveButton = QPushButton('Save')
        self.buttonsLayout.addWidget(self.prevButton)
        self.buttonsLayout.addWidget(self.saveButton)
        self.buttonsLayout.addWidget(self.nextButton)
        self.layout.addLayout(self.buttonsLayout)

        self.nextButton.clicked.connect(self.next_push)
        self.prevButton.clicked.connect(self.prev_push)
        self.saveButton.clicked.connect(self.save_push)

        self.goodCheck.stateChanged.connect(self.toggle_good)
        self.interestingCheck.stateChanged.connect(self.toggle_interesting)
        self.workCheck.stateChanged.connect(self.toggle_work)
        self.svdCheck.stateChanged.connect(self.toggle_svd)

        self.goodCheck.setEnabled(False)
        self.svdCheck.setEnabled(False)
        self.workCheck.setEnabled(False)
        self.interestingCheck.setEnabled(False)
        self.nextButton.setEnabled(False)
        self.prevButton.setEnabled(False)
        self.saveButton.setEnabled(False)

    def loadpush(self):
        oldFolder = self.folder
        newFolder = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        if not newFolder == '':
            self.set_folder(newFolder)
        

    def set_folder(self,folder):
        self.folder = folder
        self.fileList = np.sort(np.array([f for f in os.listdir(self.folder) if f[-5:]=='.fits']))
        self.pulsars = np.unique(np.array([f.split('_')[0] for f in self.fileList]))
        self.set_labels()
        self.parent.clear_plots()
        self.currentID = 0
        self.params = {}
        if self.fileList.shape[0]>0:
            self.parent.tabs.setTabEnabled(1,True)
            if os.path.isfile(os.path.join(self.folder,'pirate_params.pkl')):
                with open(os.path.join(self.folder,'pirate_params.pkl'), 'rb') as f:
                    self.params = pkl.load(f)
            self.set_file()
            self.goodCheck.setEnabled(True)
            self.svdCheck.setEnabled(True)
            self.workCheck.setEnabled(True)
            self.interestingCheck.setEnabled(True)
            if self.fileList.shape[0]>1:
                self.nextButton.setEnabled(True)
                self.prevButton.setEnabled(True)
            self.saveButton.setEnabled(True)
        else:
            self.emptyFolder.emit()
            self.goodCheck.setEnabled(False)
            self.svdCheck.setEnabled(False)
            self.workCheck.setEnabled(False)
            self.interestingCheck.setEnabled(False)
            self.nextButton.setEnabled(False)
            self.prevButton.setEnabled(False)
            self.saveButton.setEnabled(False)
        for file in self.fileList:
            if file not in self.params.keys():
                self.params.update({file : {'good' : False, 'svd': False, 'work' : False, 'interesting' : False}})

    def set_file(self):
        self.loading = True
        self.goodCheck.setChecked(False)
        self.svdCheck.setChecked(False)
        self.interestingCheck.setChecked(False)
        self.workCheck.setChecked(False)
        self.parent.clear_plots()
        if self.fileList.shape[0]>0:
            if self.fileList[self.currentID] in self.params.keys():
                if 'good' in self.params[self.fileList[self.currentID]].keys():
                    self.goodCheck.setChecked(self.params[self.fileList[self.currentID]]['good'])
                if 'work' in self.params[self.fileList[self.currentID]].keys():
                    self.workCheck.setChecked(self.params[self.fileList[self.currentID]]['work'])
                if 'interesting' in self.params[self.fileList[self.currentID]].keys():
                    self.interestingCheck.setChecked(self.params[self.fileList[self.currentID]]['interesting'])
                if 'svd' in self.params[self.fileList[self.currentID]].keys():
                    self.svdCheck.setChecked(self.params[self.fileList[self.currentID]]['svd'])
            else:
                self.params.update({self.fileList[self.currentID] : {'good' : False, 'work' : False, 'svd' : False, 'interesting' : False}})
            self.parent.dyn.load(os.path.join(self.folder,self.fileList[self.currentID]),
                                self.svdCheck.isChecked(),
                                self.goodCheck.isChecked(),
                                self.workCheck.isChecked(),
                                self.interestingCheck.isChecked()
                                )
        self.parent.update_plots()
        self.loading = False

    def next_push(self):
        if self.fileList.shape[0]>0:
            self.currentID = np.mod(self.currentID+1,self.fileList.shape[0])
            self.set_file()

    def prev_push(self):
        if self.fileList.shape[0]>0:
            self.currentID = np.mod(self.currentID-1,self.fileList.shape[0])
            self.set_file()

    def save_push(self):
        with open(os.path.join(self.folder,'pirate_params.pkl'), 'wb') as f:
            pkl.dump(self.params, f, pkl.HIGHEST_PROTOCOL)

    def toggle_good(self,state):
        if not self.loading:
            if self.fileList.shape[0]>0:
                self.params[self.fileList[self.currentID]].update({'good' : self.goodCheck.isChecked()})

    def toggle_work(self,state):
        if not self.loading:
            if self.fileList.shape[0]>0:
                self.params[self.fileList[self.currentID]].update({'work' : self.workCheck.isChecked()})

    def toggle_interesting(self,state):
        if not self.loading:
            if self.fileList.shape[0]>0:
                self.params[self.fileList[self.currentID]].update({'interesting' : self.interestingCheck.isChecked()})

    def toggle_svd(self,state):
        if not self.loading:
            if self.fileList.shape[0]>0:
                self.params[self.fileList[self.currentID]].update({'svd' : self.svdCheck.isChecked()})
                self.parent.dyn.set_svd(self.svdCheck.isChecked())
                self.parent.clear_plots()
                self.parent.update_plots()

    def on_dyn_load(self):
        pass

class DynWidget(QWidget):
    newLoaded = pyqtSignal()
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.dyn = Dynspec()
        self.filename = ''
        self.singleThetaTheta = False
        
        self.good = False
        self.work = False
        self.interesting = False
        self.svd = False
        self.set = False

    def load(self,filename,svd,good,work,interesting):
        if self.filename !=filename.split('/')[-1]:
            self.set = True
            self.dspec,self.time,self.freq = dspec_from_fits(filename)
            self.singleThetaTheta = False
            self.filename = filename.split('/')[-1]
            self.work = work
            self.interesting = interesting
            self.good = good
            self.set_svd(svd)
            
        
    def set_svd(self,value):
        self.dyn = Dynspec()
        if value:
            dspec = do_svd(self.dspec)
            self.dyn.load_dyn_obj(BasicDyn(dyn = dspec,times = self.time.value,freqs=self.freq.value,dt=(self.time[1]-self.time[0]).value,df=(self.freq[1]-self.freq[0]).value),process=False)
        else:
            self.dyn.load_dyn_obj(BasicDyn(dyn = self.dspec,times = self.time.value,freqs=self.freq.value,dt=(self.time[1]-self.time[0]).value,df=(self.freq[1]-self.freq[0]).value),process=False)
        self.svd = value
        self.newLoaded.emit()


    def basic_plot(self):
        dsExtent = thth.ext_find(self.time,self.freq)
        fd = thth.fft_axis(self.time,u.mHz)
        tau=thth.fft_axis(self.freq,u.us)
        ssExtent = thth.ext_find(fd,tau)
        cs = np.fft.fftshift(np.fft.fft2(np.nan_to_num(self.dyn.dyn-np.nanmean(self.dyn.dyn))))
        ss = np.abs(cs)**2
        return(self.dyn.dyn,dsExtent,ss,ssExtent)

class ThThTab(QWidget):
    parametersSet = pyqtSignal()
    saveParams = pyqtSignal()
    def __init__(self,parent,general) -> None:
        super().__init__()
        self.fastUpdate = False
        self.layout = QHBoxLayout(parent)
        self.parent = parent
        self.general = general

        self.leftLayout = QVBoxLayout()
        self.leftLayout.setAlignment(Qt.AlignHCenter)
        self._make_filter_box()
        self.fileSelect = QComboBox()
        self.fileSelect.currentTextChanged.connect(self.load)
        self.leftLayout.addWidget(self.fileSelect)
        self._make_options_box()
        self.displayLayout = QHBoxLayout()
        self.leftLayout.addLayout(self.displayLayout)
        self._make_info_box()
        self.layout.addLayout(self.leftLayout)
        saveButton = QPushButton('Save Parameters')
        saveButton.clicked.connect(self.save_params)
        self.leftLayout.addWidget(saveButton)

        
        self.rightLayout = QVBoxLayout()
        self.rightLayout.setAlignment(Qt.AlignHCenter)
        self._create_diagnostic_layout()
        self.layout.addLayout(self.rightLayout)
        self.setLayout(self.layout)

    def _create_diagnostic_layout(self):
        self.diagnosticFigure = MplCanvas()
        fig = self.diagnosticFigure.fig
        grid = matplotlib.pyplot.GridSpec(4, 2)
        self.dataDspecPlot = fig.add_subplot(grid[0, 0])
        self.modelDspecPlot = fig.add_subplot(grid[0,1])
        self.dataSspecPlot = fig.add_subplot(grid[1, 0])
        self.modelSspecPlot = fig.add_subplot(grid[1,1])
        self.dataTHTHPlot = fig.add_subplot(grid[2,0])
        self.modelTHTHPlot = fig.add_subplot(grid[2,1])
        self.evoPlot = fig.add_subplot(grid[3,:])

        self.singleButton = QPushButton('Run Single')
        self.singleButton.setEnabled(False)
        self.singleButton.clicked.connect(self._run_single)
        self.rightLayout.addWidget(self.singleButton)
        self.rightLayout.addWidget(self.diagnosticFigure)
        self.savePlotButton = QPushButton('Save Figure')
        self.savePlotButton.clicked.connect(self._save_figure)
        self.savePlotButton.setEnabled(False)
        self.rightLayout.addWidget(self.savePlotButton)

    def _save_figure(self):
        self.diagnosticFigure.save(self.general.folder,self.parent.dyn.filename.split('.')[0]+'_single.png')

    def clear_plots(self):
        self.dataDspecPlot.cla()
        self.modelDspecPlot.cla()
        self.dataSspecPlot.cla()
        self.modelSspecPlot.cla()
        self.dataTHTHPlot.cla()
        self.modelTHTHPlot.cla()
        self.evoPlot.cla()
    
    def _run_single(self):
        self.clear_plots()
        self.singleButton.setEnabled(False)
        self.savePlotButton.setEnabled(False)
        if self.parent.dyn.singleThetaTheta:
            etas,eigs,popt = self.parent.dyn.dyn.thetatheta_single(plot=False,arrays=True)
            self._make_plots(etas,eigs,popt)
        self.singleButton.setEnabled(True)
        self.savePlotButton.setEnabled(True)

    def _make_plots(self,etas,measure,fit_res):
        ncf = self.ncfControl.value()
        cwf = self.parent.dyn.freq.shape[0]//ncf
        cwf -= np.mod(cwf,2)
        nct = self.nctControl.value()
        cwt = self.parent.dyn.time.shape[0]//nct
        cwt -= np.mod(cwt,2)

        dspec = self.parent.dyn.dyn.dyn[:cwf,:cwt]
        time = self.parent.dyn.time[:cwt]
        freq = self.parent.dyn.freq[:cwf]
        fd = thth.fft_axis(time,u.mHz, self.parent.dyn.dyn.npad)
        tau = thth.fft_axis(freq,u.us, self.parent.dyn.dyn.npad)

        dspec2 = np.copy(dspec)
        mn = np.nanmean(dspec2)
        dspec2 -= mn
        dspec_pad = np.pad(np.nan_to_num(dspec2), ((
            0, self.parent.dyn.dyn.npad*cwf), (0, self.parent.dyn.dyn.npad*cwt)), mode='constant',
            constant_values=0)
        CS = np.fft.fftshift(np.fft.fft2(dspec_pad))
        CS[np.abs(tau)<self.parent.dyn.dyn.thth_tau_mask]=0
        
        edges = self.parent.dyn.dyn.edges

        etas0 = etas[np.isfinite(measure)]
        eigs0 = measure[np.isfinite(measure)]

        etas_fit = etas0[np.abs(etas0 - etas0[eigs0 == eigs0.max()])
                        < self.parent.dyn.dyn.fw * etas0[eigs0 == eigs0.max()]]
        eigs_fit = eigs0[np.abs(etas0 - etas0[eigs0 == eigs0.max()])
                        < self.parent.dyn.dyn.fw * etas0[eigs0 == eigs0.max()]]

        if fit_res is None:
            eta_fit = etas.mean()
            eta_sig = etas.mean()/100.
        else:
            eta_fit = fit_res[1]*u.us/u.mHz**2
            eta_sig = np.sqrt((eigs_fit - thth.chi_par(etas_fit.value,
                                *fit_res)).std() / np.abs(fit_res[0]))*u.us/u.mHz**2

        # Verify units
        tau = thth.unit_checks(tau, "tau", u.us)
        fd = thth.unit_checks(fd, "fd", u.mHz)
        edges = thth.unit_checks(edges, "edges", u.mHz)
        eta_fit = thth.unit_checks(eta_fit, "eta_fit", u.s**3)
        eta_sig = thth.unit_checks(eta_sig, "eta_sig", u.s**3)
        etas = thth.unit_checks(etas, "etas", u.s**3)
        etas_fit = thth.unit_checks(etas_fit, "etas_fit", u.s**3)

        tau_lim = tau.max()

        # Determine fd limits
        fd_lim = min(2 * edges.max(), fd.max()).value

        # Determine TH-TH and model
        if np.isnan(eta_fit):
            eta = etas.mean()
            thth_red, thth2_red, recov, model, edges_red, w, V = thth.modeler(
                CS, tau, fd, etas.mean(), edges
            )
        else:
            eta = eta_fit
            thth_red, thth2_red, recov, model, edges_red, w, V = thth.modeler(
                CS, tau, fd, eta_fit, edges
            )

        # Create model Wavefield and Conjugate Wavefield
        ththE_red = thth_red * 0
        ththE_red[ththE_red.shape[0] // 2, :] = np.conjugate(V) * np.sqrt(w)
        # Map back to time/frequency space
        recov_E = thth.rev_map(ththE_red, tau, fd, eta, edges_red, hermetian=False)
        model_E = np.fft.ifft2(np.fft.ifftshift(recov_E))[
            : dspec.shape[0], : dspec.shape[1]
        ]
        model_E *= dspec.shape[0] * dspec.shape[1] / 4
        model_E[dspec > 0] = np.sqrt(dspec[dspec > 0]) * np.exp(
            1j * np.angle(model_E[dspec > 0])
        )
        model_E = np.pad(
            model_E,
            (
                (0, CS.shape[0] - model_E.shape[0]),
                (0, CS.shape[1] - model_E.shape[1]),
            ),
            mode="constant",
            constant_values=0,
        )
        recov_E = np.abs(np.fft.fftshift(np.fft.fft2(model_E))) ** 2
        model_E = model_E[: dspec.shape[0], : dspec.shape[1]]
        N_E = recov_E[: recov_E.shape[0] // 4, :].mean()

        model = model[: dspec.shape[0], : dspec.shape[1]]
        model -= model.mean()
        model *= np.nanstd(dspec) / np.std(model)
        model += np.nanmean(dspec)

        # Data Dynamic Spectrum
        self.dataDspecPlot.imshow(
            dspec,
            aspect="auto",
            extent=thth.ext_find(time.to(u.min), freq),
            origin="lower",
            vmin=np.nanmean(dspec) - 5 * np.nanstd(dspec),
            vmax=np.nanmean(dspec) + 5 * np.nanstd(dspec),
            cmap='magma',
        )
        self.dataDspecPlot.set_xlabel("Time (min)")
        self.dataDspecPlot.set_ylabel("Freq (MHz)")
        self.dataDspecPlot.set_title("Data Dynamic Spectrum")

        # Model Dynamic Spectrum
        self.modelDspecPlot.imshow(
            model[: dspec.shape[0], : dspec.shape[1]],
            aspect="auto",
            extent=thth.ext_find(time.to(u.min), freq),
            origin="lower",
            vmin=np.nanmean(dspec) - 5 * np.nanstd(dspec),
            vmax=np.nanmean(dspec) + 5 * np.nanstd(dspec),
            cmap='magma',
        )
        self.modelDspecPlot.set_xlabel("Time (min)")
        self.modelDspecPlot.set_yticks([])
        self.modelDspecPlot.set_title("Model Dynamic Spectrum")

        # Data Secondary Spectrum
        SS = np.abs(CS) ** 2
        self.dataSspecPlot.imshow(
            SS,
            norm=LogNorm(
                vmin=np.median(SS[SS>0]), vmax=SS.max()
            ),
            origin="lower",
            aspect="auto",
            extent=thth.ext_find(fd, tau),
            cmap='magma'
        )
        self.dataSspecPlot.set_xlim((-fd_lim, fd_lim))
        self.dataSspecPlot.set_ylim((0, tau_lim.value))
        self.dataSspecPlot.set_xlabel(r"$f_D$ (mHz)")
        self.dataSspecPlot.set_ylabel(r"$\tau$ (us)")
        self.dataSspecPlot.set_title("Data Secondary Spectrum")
        self.dataSspecPlot.plot(fd, eta * (fd**2), "r", alpha=0.7,lw=1)
        

        # Model Secondary Spectrum
        self.modelSspecPlot.imshow(
            np.abs(recov) ** 2,
            norm=LogNorm(
                vmin=np.median(SS[SS>0]), vmax=SS.max()
            ),
            origin="lower",
            aspect="auto",
            extent=thth.ext_find(fd, tau),
            cmap='magma',
        )
        self.modelSspecPlot.set_xlim((-fd_lim, fd_lim))
        self.modelSspecPlot.set_ylim((0, tau_lim.value))
        self.modelSspecPlot.set_title("Model Secondary Spectrum")
        self.modelSspecPlot.set_xlabel(r"$f_D$ (mHz)")
        self.modelSspecPlot.set_yticks([])

        Sthth = np.abs(thth_red) ** 2
        # Data TH-TH
        self.dataTHTHPlot.imshow(
            Sthth,
            norm=LogNorm(
                vmin=np.median(Sthth[Sthth>0]),
                vmax=Sthth.max(),
            ),
            origin="lower",
            # aspect="auto",
            extent=[
                edges_red[0].value,
                edges_red[-1].value,
                edges_red[0].value,
                edges_red[-1].value,
            ],
            cmap='magma',
        )
        self.dataTHTHPlot.set_xlabel(r"$\theta_1$")
        self.dataTHTHPlot.set_ylabel(r"$\theta_2$")
        self.dataTHTHPlot.set_title(r"Data $\theta-\theta$")

        # Model TH-TH
        self.modelTHTHPlot.imshow(
            np.abs(thth2_red) ** 2,
            norm=LogNorm(
                vmin=np.median(Sthth[Sthth>0]),
                vmax=Sthth.max(),
            ),
            origin="lower",
            # aspect="auto",
            extent=[
                edges_red[0].value,
                edges_red[-1].value,
                edges_red[0].value,
                edges_red[-1].value,
            ],
            cmap='magma'
        )
        self.modelTHTHPlot.set_xlabel(r"$\theta_1$")
        self.modelTHTHPlot.set_yticks([])
        self.modelTHTHPlot.set_title(r"Model $\theta-\theta$")

        self.evoPlot.plot(etas, measure,lw=1)
        if not fit_res is None:
            fit_string, err_string = thth.errString(eta_fit, eta_sig)
            self.evoPlot.plot(
                etas_fit,
                thth.chi_par(etas_fit.value, *fit_res),
                label=r"$\eta$ = %s $\pm$ %s $s^3$" % (fit_string, err_string),lw=1
            )
            self.evoPlot.legend()
        self.evoPlot.set_title("Eigenvalue Search")
        self.evoPlot.set_ylabel(r"Largest Eigenvalue")
        self.evoPlot.set_xlabel(r"$\eta$ ($s^3$)")

        self.diagnosticFigure.fig.tight_layout()
        self.diagnosticFigure.draw()

    def _make_filter_box(self):
        filterBox = QHBoxLayout()
        self.goodCheck = QCheckBox('Good')
        self.interestingCheck = QCheckBox('Interesting')
        self.svdCheck = QCheckBox('SVD')
        self.workCheck = QCheckBox('Needs Work')
        filterBox.addWidget(self.goodCheck)
        filterBox.addWidget(self.interestingCheck)
        filterBox.addWidget(self.svdCheck)
        filterBox.addWidget(self.workCheck)
        self.leftLayout.addLayout(filterBox)
        self.goodCheck.stateChanged.connect(self.update_items)
        self.svdCheck.stateChanged.connect(self.update_items)
        self.interestingCheck.stateChanged.connect(self.update_items)
        self.workCheck.stateChanged.connect(self.update_items)

    def filter_from_dyn(self):
        self.goodCheck.setChecked(self.parent.dyn.good)
        self.svdCheck.setChecked(self.parent.dyn.svd)
        self.workCheck.setChecked(self.parent.dyn.work)
        self.interestingCheck.setChecked(self.parent.dyn.interesting)

    def update_items(self):
        if not self.fastUpdate:
            self.fastUpdate = True
            prevSelection = self.fileSelect.currentText()
            self.fileSelect.clear()
            allFiles = np.copy(self.general.fileList)
            if self.goodCheck.isChecked():
                mask = [self.general.params[f]['good'] for f in allFiles]
                allFiles = allFiles[mask]
            if self.interestingCheck.isChecked():
                mask = [self.general.params[f]['interesting'] for f in allFiles]
                allFiles = allFiles[mask]
            if self.workCheck.isChecked():
                mask = [self.general.params[f]['work'] for f in allFiles]
                allFiles = allFiles[mask]
            if self.svdCheck.isChecked():
                mask = [self.general.params[f]['svd'] for f in allFiles]
                allFiles = allFiles[mask]  
            self.fileSelect.addItems(allFiles)
            oldIndex = self.fileSelect.findText(prevSelection)
            if oldIndex >= 0:
                self.fileSelect.setCurrentIndex(oldIndex)
            elif len(self.fileSelect.currentText())>0:
                self.fastUpdate = False
                self.load(self.fileSelect.currentText())
            self.fastUpdate = False

    def _make_options_box(self):
        optionsBox = QVBoxLayout()
        splitLayout = QHBoxLayout()
        dspecOptions = QHBoxLayout()
        freqOptions = QVBoxLayout()
        timeOptions = QVBoxLayout()
        sspecOption = QVBoxLayout()

        self.ncfControl = QSpinBox()
        self.ncfControl.setSuffix(' frequency chunks')
        self.ncfControl.setMinimum(1)
        self.ncfControl.setAlignment(Qt.AlignHCenter)
        self.nctControl = QSpinBox()
        self.nctControl.setMinimum(1)
        self.nctControl.setSuffix(' time chunks')
        self.nctControl.setAlignment(Qt.AlignHCenter)

        self.ncfControl.valueChanged.connect(self._calc_cwf)
        self.nctControl.valueChanged.connect(self._calc_cwt)

        self.cwfLabel = QLabel()
        self.cwfLabel.setText('___ Chanels per Chunk')
        self.cwtLabel = QLabel()
        self.cwtLabel.setText('___ Bins per Chunk')

        limitBox = QHBoxLayout()
        limitBox.setAlignment(Qt.AlignHCenter)
        limitLabel = QLabel()
        limitLabel.setText('Edges Limit :')
        self.limitControl = QDoubleSpinBox()
        self.limitControl.setMinimum(0)
        self.limitControl.setSuffix(' mHz')
        limitBox.addWidget(limitLabel)
        limitBox.addWidget(self.limitControl)

        maskBox = QHBoxLayout()
        maskBox.setAlignment(Qt.AlignHCenter)
        maskLabel = QLabel()
        maskLabel.setTextFormat(Qt.RichText)
        maskLabel.setText('Mask \u03C4 below')
        self.maskControl = QDoubleSpinBox()
        self.maskControl.setMinimum(0)
        self.maskControl.setSuffix(' \u03BCs')
        maskBox.addWidget(maskLabel)
        maskBox.addWidget(self.maskControl)

        etaBox = QHBoxLayout()
        etaBox.setAlignment(Qt.AlignHCenter)
        minBox = QHBoxLayout()
        minBox.setAlignment(Qt.AlignHCenter)
        maxBox = QHBoxLayout()
        maxBox.setAlignment(Qt.AlignHCenter)

        minLabel = QLabel()
        minLabel.setTextFormat(Qt.RichText)
        minLabel.setText('Min \u03B7:')
        self.minControl = QDoubleSpinBox()
        self.minControl.setMinimum(0)
        self.minControl.setSuffix(' s\u00b3')
        self.minControl.setDecimals(4)
        minBox.addWidget(minLabel)
        minBox.addWidget(self.minControl)

        maxLabel = QLabel()
        maxLabel.setTextFormat(Qt.RichText)
        maxLabel.setText('Max \u03B7:')
        self.maxControl = QDoubleSpinBox()
        self.maxControl.setMinimum(0)
        self.maxControl.setSuffix(' s\u00b3')
        self.maxControl.setDecimals(4)
        maxBox.addWidget(maxLabel)
        maxBox.addWidget(self.maxControl)

        etaBox.addLayout(minBox)
        etaBox.addLayout(maxBox)

        self.maskControl.valueChanged.connect(self.update_overlay)
        self.limitControl.valueChanged.connect(self.update_overlay)
        self.minControl.valueChanged.connect(self.update_overlay)
        self.maxControl.valueChanged.connect(self.update_overlay)

        timeOptions.setAlignment(Qt.AlignHCenter)
        freqOptions.setAlignment(Qt.AlignHCenter)
        timeOptions.addWidget(self.nctControl)
        freqOptions.addWidget(self.ncfControl)
        timeOptions.addWidget(self.cwtLabel)
        freqOptions.addWidget(self.cwfLabel)
        dspecOptions.addLayout(freqOptions)
        dspecOptions.addLayout(timeOptions)
        splitLayout.addLayout(dspecOptions)

        sspecOption.addLayout(limitBox)
        sspecOption.addLayout(maskBox)
        sspecOption.addLayout(etaBox)
        splitLayout.addLayout(sspecOption)
        optionsBox.addLayout(splitLayout)

        fwHolder = QHBoxLayout()
        self.fwSlider = QSlider(Qt.Orientation.Horizontal)
        self.fwSlider.setMinimum(1)
        self.fwSlider.setMaximum(99)
        self.fwSlider.setValue(10)
        self.fwSlider.valueChanged.connect(self.update_thetatheta_params)
        self.fwSlider.valueChanged.connect(self.update_fw_display)
        fwLabel = QLabel()
        # fwLabel.setAlignment(Qt.AlignHStart)
        fwLabel.setText('Fractional Fitting Width')
        self.fwValue = QLabel()
        # self.fwValue.setAlignment(Qt.AlignHRight)
        self.fwValue.setText('10 %')
        fwHolder.addWidget(fwLabel)
        fwHolder.addWidget(self.fwSlider)
        fwHolder.addWidget(self.fwValue)
        optionsBox.addLayout(fwHolder)
        
        
        self.leftLayout.addLayout(optionsBox)

    def update_fw_display(self):
        self.fwValue.setText(f'{str(self.fwSlider.value()).zfill(2)} %')

    def _make_info_box(self):
        infoBox = QVBoxLayout()
        infoBox.setAlignment(Qt.AlignHCenter)
        self.netaLabel = QLabel()
        self.netaLabel.setText(f'Number of \u03B7 to search : ---')
        self.nedgeLabel = QLabel()
        self.nedgeLabel.setText(f'\u03B8-\u03B8 size : ---')
        self.nedgeLabel.setAlignment(Qt.AlignHCenter)
        self.netaLabel.setAlignment(Qt.AlignHCenter)

        infoBox.addWidget(self.netaLabel)
        infoBox.addWidget(self.nedgeLabel)
        self.leftLayout.addLayout(infoBox)

    def limits_from_dyn(self):
        self.ncfControl.setMaximum(self.parent.dyn.freq.shape[0]//2)
        self.nctControl.setMaximum(self.parent.dyn.time.shape[0]//2)
        fd = thth.fft_axis(self.parent.dyn.time,u.mHz)
        tau = thth.fft_axis(self.parent.dyn.time,u.mHz)
        self.limitControl.setMaximum(fd.max().value/2)

    def on_open(self):
        pass

    def _calc_cwf(self,value):
        target = self.parent.dyn.freq.shape[0]//value
        target -= np.mod(target,2)
        self.cwfLabel.setText(f'{target} Chanels per Chunk')
        self.update_overlay()

    def _calc_cwt(self,value):
        target = self.parent.dyn.time.shape[0]//value
        target -= np.mod(target,2)
        self.cwtLabel.setText(f'{target} Bins per Chunk')
        self.update_overlay()

    def update_overlay(self):
        if not self.fastUpdate:
            self.parent.remove_overlay()
            ncf = self.ncfControl.value()
            cwf = self.parent.dyn.freq.shape[0]//ncf
            cwf -= np.mod(cwf,2)
            nct = self.nctControl.value()
            cwt = self.parent.dyn.time.shape[0]//nct
            cwt -= np.mod(cwt,2)
            edges_lim = self.limitControl.value()
            tau_mask = self.maskControl.value()
            eta_min = self.minControl.value()
            eta_max = self.maxControl.value()
            self.parent.add_overlay(ncf,cwf,nct,cwt,edges_lim,tau_mask,eta_min,eta_max)
            self.update_thetatheta_params()         

    def update_thetatheta_params(self):
        if not self.fastUpdate:
            ncf = self.ncfControl.value()
            cwf = self.parent.dyn.freq.shape[0]//ncf
            cwf -= np.mod(cwf,2)
            nct = self.nctControl.value()
            cwt = self.parent.dyn.time.shape[0]//nct
            cwt -= np.mod(cwt,2)
            edges_lim = self.limitControl.value()
            tau_mask = self.maskControl.value()
            eta_min = self.minControl.value()
            eta_max = self.maxControl.value()
            fw = self.fwSlider.value()/100.
            if eta_max > eta_min and eta_min > 0 and edges_lim >0:
                    self.parent.dyn.dyn.prep_thetatheta(cwf=cwf,cwt = cwt,
                                                        edges_lim = edges_lim*u.mHz,
                                                        tau_mask = tau_mask*u.us,
                                                        eta_min = eta_min*u.s**3,eta_max = eta_max*u.s**3,
                                                        fw = fw,
                                                        verbose=False)
                    self.parent.dyn.singleThetaTheta = True
                    self.singleButton.setEnabled(True)
                    self.netaLabel.setText(f'Number of \u03B7 to search : {self.parent.dyn.dyn.neta}')
                    self.nedgeLabel.setText(f'\u03B8-\u03B8 size : {self.parent.dyn.dyn.edges.shape[0]-1}')
                    self.parametersSet.emit()

    def load(self,text):
        if not self.fastUpdate:
            self.parent.dyn.load(os.path.join(self.general.folder,text),
            self.general.params[text]['svd'],
            self.general.params[text]['good'],
            self.general.params[text]['work'],
            self.general.params[text]['interesting'],
            )

    def on_dyn_load(self):
        if self.parent.dyn.set:
            self.fastUpdate = True
            self.singleButton.setEnabled(False)
            self.savePlotButton.setEnabled(False)
            self.clear_plots()
            self.limits_from_dyn()
            
            if 'theta-theta' in self.general.params[self.parent.dyn.filename].keys():
                ncf = self.general.params[self.parent.dyn.filename]['theta-theta']['ncf']
                nct = self.general.params[self.parent.dyn.filename]['theta-theta']['nct']
                fw = int(self.general.params[self.parent.dyn.filename]['theta-theta']['fw']*100)
                edges_lim = self.general.params[self.parent.dyn.filename]['theta-theta']['edges_lim']
                tau_mask = self.general.params[self.parent.dyn.filename]['theta-theta']['tau_mask']
                eta_min = self.general.params[self.parent.dyn.filename]['theta-theta']['eta_min']
                eta_max = self.general.params[self.parent.dyn.filename]['theta-theta']['eta_max']
                self.ncfControl.setValue(ncf)
                self.nctControl.setValue(nct)
                self.minControl.setValue(eta_min)
                self.maxControl.setValue(eta_max)
                self.limitControl.setValue(edges_lim)
                self.maskControl.setValue(tau_mask)
                self.fwSlider.setValue(fw)
                self._calc_cwf(self.ncfControl.value())
                self.fastUpdate = False
                self._calc_cwt(self.nctControl.value())
                self.update_fw_display()
                self.singleButton.setEnabled(True)
                self.parametersSet.emit()
            else:
                self.netaLabel.setText(f'Number of \u03B7 to search : ---')
                self.nedgeLabel.setText(f'\u03B8-\u03B8 size : ---')
                self.ncfControl.setValue(1)
                self.nctControl.setValue(1)
                self.minControl.setValue(0)
                self.maxControl.setValue(0)
                self.limitControl.setValue(0)
                self.maskControl.setValue(0)
                self.fwSlider.setValue(10)
                self.update_fw_display()
                self._calc_cwf(1)
                self._calc_cwt(1)
                self.diagnosticFigure.fig.tight_layout()
                self.diagnosticFigure.draw()
            self.fastUpdate = False

    def save_params(self):
        ncf = self.ncfControl.value()
        cwf = self.parent.dyn.freq.shape[0]//ncf
        cwf -= np.mod(cwf,2)
        nct = self.nctControl.value()
        cwt = self.parent.dyn.time.shape[0]//nct
        cwt -= np.mod(cwt,2)
        edges_lim = self.limitControl.value()
        tau_mask = self.maskControl.value()
        eta_min = self.minControl.value()
        eta_max = self.maxControl.value()
        fw = self.fwSlider.value()/100.
        self.general.params[self.parent.dyn.filename].update({'theta-theta' : { 'ncf' : ncf,
                                                                                'nct' : nct,
                                                                                'edges_lim' : edges_lim,
                                                                                'fw' : fw,
                                                                                'eta_min' : eta_min,
                                                                                'eta_max' : eta_max,
                                                                                'tau_mask' : tau_mask
                                                                                }})
        self.saveParams.emit()

class PhaseRetrievalTab(QWidget):
    def __init__(self,parent) -> None:
        self.parent=parent
        super().__init__()
        self.layout = QVBoxLayout(self)
        self.runButton = QPushButton('Run')
        self.saveButton = QPushButton('Save Figure')
        self.mainPlots = MplCanvas()
        self.wavefieldHolder = MplCanvas()
        self.layout.addWidget(self.runButton)
        self.layout.addWidget(self.mainPlots)
        self.layout.addWidget(self.wavefieldHolder)
        self.saveButton.setEnabled(False)
        

        self.tauSlider = QSlider(Qt.Orientation.Horizontal)
        self.tauSlider.setMinimum(1)
        self.tauSlider.setMaximum(100)
        self.tauSlider.setValue(100)
        tauLabel = QLabel()
        tauLabel.setAlignment(Qt.AlignHCenter)
        tauLabel.setText('\u03C4 Extent')
        self.layout.addWidget(tauLabel)
        self.layout.addWidget(self.tauSlider)

        self.fdSlider = QSlider(Qt.Orientation.Horizontal)
        self.fdSlider.setMinimum(1)
        self.fdSlider.setMaximum(100)
        self.fdSlider.setValue(100)
        fdLabel = QLabel()
        fdLabel.setAlignment(Qt.AlignHCenter)
        fdLabel.setText('fD Extent')
        self.layout.addWidget(fdLabel)
        self.layout.addWidget(self.fdSlider)

        self.fdSlider.sliderMoved.connect(self.fd_scale)
        self.tauSlider.sliderMoved.connect(self.tau_scale)

        self.layout.addWidget(self.saveButton)
        grid = matplotlib.pyplot.GridSpec(2,2)

        self.evoPlot = self.mainPlots.fig.add_subplot(grid[0,:])
        self.dataDspec = self.mainPlots.fig.add_subplot(grid[1,0])
        self.modelDspec = self.mainPlots.fig.add_subplot(grid[1,1])
        self.wavefield = self.wavefieldHolder.fig.add_subplot(111)

        self.runButton.clicked.connect(self.run)

    def fd_scale(self):
        dyn = self.parent.dyn
        fd = thth.fft_axis(dyn.time[:dyn.dyn.wavefield.shape[1]],u.mHz)
        extent = thth.ext_find(fd,fd)
        mx= extent[1]
        fraction = self.fdSlider.value()/100.
        self.wavefield.set_xlim((-fraction*mx,fraction*mx))
        self.wavefieldHolder.draw()

    def tau_scale(self):
        dyn = self.parent.dyn
        tau = thth.fft_axis(dyn.freq[:dyn.dyn.wavefield.shape[0]],u.us)
        extent = thth.ext_find(tau,tau)
        mx= extent[1]
        fraction = self.tauSlider.value()/100.
        self.wavefield.set_ylim((0,fraction*mx))
        self.wavefieldHolder.draw()

    def run(self):
        self.clear_plots()
        self.parent.dyn.dyn.fit_thetatheta()
        self.make_evo_plot()
        self.parent.dyn.dyn.calc_wavefield()
        self.make_dspec_plots()
        self.parent.dyn.dyn.calc_wavefield(gs=True,niter=10)
        self.mainPlots.fig.tight_layout()
        self.mainPlots.draw()
        self.make_wavefield_plot()
        self.wavefieldHolder.fig.tight_layout()
        self.wavefieldHolder.draw()
        
        self.saveButton.setEnabled(True)

    def make_evo_plot(self):
        time_avg = True
        dyn = self.parent.dyn.dyn
        if time_avg and dyn.nct_fit>1:
            eta_avg = np.nanmean(dyn.eta_evo,1)
            eta_count = np.nansum(dyn.eta_evo,1)/eta_avg
            avg_err = np.nanstd(dyn.eta_evo,1)/np.sqrt(eta_count-1)
            tofit =  np.isfinite(eta_avg)*np.isfinite(avg_err)
            A = (np.sum(eta_avg[tofit] / (dyn.f0s * 
                                          avg_err)[tofit] ** 2) /
                np.sum(1 / (dyn.f0s**2 * avg_err)[tofit] ** 2)).to(u.s**3 * u.MHz**2)
            A_err = np.sqrt(1 / np.sum(2 / ((dyn.f0s**2) * avg_err)[tofit] ** 2)).to(u.s**3 * u.MHz**2)
        else:
            tofit =  np.isfinite(dyn.eta_evo)*np.isfinite(dyn.eta_evo_err)
            A = (np.sum(dyn.eta_evo[tofit] / (dyn.f0s[:, np.newaxis] *
                                               dyn.eta_evo_err)[tofit] ** 2) /
                np.sum(1 / ((dyn.f0s[:, np.newaxis]**2) *
                            dyn.eta_evo_err)[tofit] ** 2)).to(u.s**3 * u.MHz**2)
            A_err = np.sqrt(
                1 / np.sum(2 /
                           ((dyn.f0s[:, np.newaxis]**2) *
                            dyn.eta_evo_err)[tofit] ** 2)).to(u.s**3 * u.MHz**2)
        dyn.ththeta = A/dyn.fref**2
        dyn.ththetaerr = A_err/dyn.fref**2

        fr = dyn.fref
        fit_string, err_string = \
            thth.errString(dyn.ththeta*(fr/np.floor(fr))**2,
                            dyn.ththetaerr*(fr/np.floor(fr))**2)
        if time_avg and dyn.nct_fit>1:
            self.evoPlot.errorbar(dyn.f0s.value, eta_avg, yerr=avg_err, fmt='.')
        else:
            self.evoPlot.errorbar(np.ravel(dyn.f0s.value[:,np.newaxis] *
                                    np.ones(dyn.eta_evo.shape)),
                            np.ravel(dyn.eta_evo.value),
                            yerr=np.ravel(dyn.eta_evo_err.value), fmt='.')
        self.evoPlot.plot(dyn.f0s,A/dyn.f0s**2,label = r'$\eta_{%s}$ = %s $\pm$ %s $s^3$' %
        (np.floor(dyn.fref),fit_string, err_string),lw=1)
        self.evoPlot.set_xlabel(r'$\rm{Freq}~\left(\rm{MHz}\right)$')
        self.evoPlot.set_ylabel(r'$\eta~\left(\rm{s}^3\right)$')
        self.evoPlot.legend()

    def make_dspec_plots(self):
        dyn = self.parent.dyn.dyn
        wf = dyn.wavefield
        tcut = dyn.times[:wf.shape[1]]*u.s
        fcut = dyn.freqs[:wf.shape[0]]*u.MHz
        modelDS = np.abs(wf)**2
        modelDS /= modelDS.mean()
        dataDS = np.copy(dyn.dyn[:wf.shape[0],:wf.shape[1]])
        dataDS/=np.nanmean(dataDS)

        self.dataDspec.imshow(dataDS,
                                origin='lower',aspect='auto',
                                extent=thth.ext_find(tcut,fcut),
                                vmin = np.percentile(modelDS,1),vmax = np.percentile(modelDS,99),
                                cmap='magma')
        self.modelDspec.imshow(modelDS,
                                origin='lower',aspect='auto',
                                extent=thth.ext_find(tcut,fcut),
                                vmin = np.percentile(modelDS,1),vmax = np.percentile(modelDS,99),
                                cmap='magma')
        self.dataDspec.set_ylabel(r'$\nu~\left(\rm{MHz}\right)$')
        self.dataDspec.set_xlabel(r'$t~\left(\rm{s}\right)$')
        self.modelDspec.set_xlabel(r'$t~\left(\rm{s}\right)$')
        self.dataDspec.set_yticks([])

    def make_wavefield_plot(self):
        dyn = self.parent.dyn.dyn
        wf = dyn.wavefield
        tcut = dyn.times[:wf.shape[1]]*u.s
        fcut = dyn.freqs[:wf.shape[0]]*u.MHz
        tau=thth.fft_axis(fcut,u.us)
        fd=thth.fft_axis(tcut,u.mHz)
        swf = np.abs(np.fft.fftshift(np.fft.fft2(wf)))**2
        extent = thth.ext_find(fd,tau)
        self.wavefield.imshow(swf,
                                origin='lower',aspect='auto',
                                extent=extent,
                                norm=LogNorm(vmin=np.percentile(swf,75)),
                                cmap='magma'
                                )
        self.wavefield.set_ylabel(r'$\tau~\left(\mu\rm{s}\right)$')
        self.wavefield.set_xlabel(r'$f_D~\left(\rm{mHz}\right)$')
        self.wavefield.set_ylim((0,extent[-1]))
        self.wavefield.plot(fd,(dyn.ththeta*fd**2),lw=1)

    def clear_plots(self):
        self.dataDspec.cla()
        self.modelDspec.cla()
        self.evoPlot.cla()
        self.wavefield.cla()
        self.saveButton.setEnabled(False)

    def on_dyn_load(self):
        self.clear_plots
        self.saveButton.setEnabled(False)

class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.left = 0
        self.top = 0
        self.width = 300
        self.height = 200
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.setWindowTitle(r"P.I.R.A.T.E.")
        self.table_widget = MyTableWidget(self)
        self.setCentralWidget(self.table_widget)
        self.showMaximized()