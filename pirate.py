## Pulsar Intesnity R Archive T Explorer

## QT Imports
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from widgets import MainWindow
import os

import warnings
warnings.filterwarnings("ignore", module="pyqt5")

app = QApplication([])
w = MainWindow()
app.exec_()