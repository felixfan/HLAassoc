#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyQt4 import QtGui, QtCore
import pandas as pd
import sys
import os

class HLAassocWin(QtGui.QWidget):
	def __init__(self):
		QtGui.QMainWindow.__init__(self)
		self.width = 700
		self.height =500
		self.setFixedSize(self.width, self.height)
		self.move(200, 200)
		self.setWindowTitle('HLAassoc')
		#######################################################
		self.gfile = ''
		self.covfile = ''
		self.outfile = ''
		#######################################################
		self.gLab = QtGui.QLabel('Genotype file')
		self.gEdit = QtGui.QLineEdit()
		self.gButton = QtGui.QPushButton("Open")
		self.gvButton = QtGui.QPushButton("View")

		self.vLab = QtGui.QLabel('Covariates file')
		self.vEdit = QtGui.QLineEdit()
		self.vButton = QtGui.QPushButton("Open")
		self.vvButton = QtGui.QPushButton("View")

		self.oLab = QtGui.QLabel('Output file')
		self.oEdit = QtGui.QLineEdit()
		self.oButton = QtGui.QPushButton("Open")

		self.traitLab = QtGui.QLabel('Trait')
		self.traitCombo = QtGui.QComboBox()
		self.traitCombo.addItem("disease trait/case-control study")
		self.traitCombo.addItem("quantitative trait")

		self.testLab = QtGui.QLabel('Test')
		self.testCombo = QtGui.QComboBox()
		self.testCombo.addItem("Pearson chi-squared test",)
		self.testCombo.addItem("Fisher's exact test")
		self.testCombo.addItem("Logistic regression")
		self.testCombo.addItem("Raw Test")
		self.testCombo.addItem("Score Test")

		self.digitLab = QtGui.QLabel('Digits')
		self.digitCombo = QtGui.QComboBox()
		self.digitCombo.addItems(["4","2","6"])

		self.modelLab = QtGui.QLabel('Model')
		self.modelCombo = QtGui.QComboBox()
		self.modelCombo.addItems(["allelic", "dominant", "recessive"])

		self.adjLab = QtGui.QLabel('adjustment')
		self.adjCombo = QtGui.QComboBox()
		self.adjCombo.addItems(["FDR", "Bonferroni", "Holm"])

		self.frqLab = QtGui.QLabel('allele frequency')
		self.frqEdit = QtGui.QLineEdit()

		self.permLab = QtGui.QLabel('Number of permutation')
		self.permEdit = QtGui.QLineEdit()

		self.covLab = QtGui.QLabel('Covariates name')
		self.covEdit = QtGui.QLineEdit()

		self.runButton = QtGui.QPushButton("Run")
		#######################################################
		self.grid = QtGui.QGridLayout()
		self.grid.setSpacing(10)
		self.grid.addWidget(self.gLab, 0, 0) # first row
		self.grid.addWidget(self.gEdit, 0, 1)
		self.grid.addWidget(self.gButton, 0, 2)
		self.grid.addWidget(self.gvButton, 0, 3)
		self.grid.addWidget(self.vLab, 1, 0) # second row
		self.grid.addWidget(self.vEdit, 1, 1)
		self.grid.addWidget(self.vButton, 1, 2)
		self.grid.addWidget(self.vvButton, 1, 3)
		self.grid.addWidget(self.oLab, 2, 0) # third row
		self.grid.addWidget(self.oEdit, 2, 1)
		self.grid.addWidget(self.oButton, 2, 2)
		self.grid.addWidget(self.traitLab, 3, 0) # fourth row
		self.grid.addWidget(self.traitCombo, 3, 1, 1, 3)
		self.grid.addWidget(self.testLab, 4, 0) # 5
		self.grid.addWidget(self.testCombo, 4, 1, 1, 3)
		self.grid.addWidget(self.digitLab, 5, 0) # 6
		self.grid.addWidget(self.digitCombo, 5, 1, 1, 3)
		self.grid.addWidget(self.modelLab, 6, 0) # 7
		self.grid.addWidget(self.modelCombo, 6, 1, 1, 3)
		self.grid.addWidget(self.adjLab, 7, 0) # 8
		self.grid.addWidget(self.adjCombo, 7, 1, 1, 3)
		self.grid.addWidget(self.frqLab, 8, 0) # 9
		self.grid.addWidget(self.frqEdit, 8, 1, 1, 3)
		self.grid.addWidget(self.permLab, 9, 0) # 10
		self.grid.addWidget(self.permEdit, 9, 1, 1, 3)
		self.grid.addWidget(self.covLab, 10, 0) # 11
		self.grid.addWidget(self.covEdit, 10, 1, 1, 3)
		self.grid.addWidget(self.runButton, 11, 1) # 12
		########################################################
		### init
		self.vLab.setEnabled(False)
		self.vEdit.setEnabled(False)
		self.vButton.setEnabled(False)
		self.vvButton.setEnabled(False)
		self.gvButton.setEnabled(False)
		self.frqEdit.setText('0.05')
		self.covLab.setEnabled(False)
		self.covEdit.setEnabled(False)
		### Events and signals
		self.connect(self.traitCombo, QtCore.SIGNAL('activated(QString)'), self.traitCombo_chosen)
		self.connect(self.testCombo, QtCore.SIGNAL('activated(QString)'), self.testCombo_chosen)
		self.gButton.clicked.connect(self.gButtonClicked)
		self.vButton.clicked.connect(self.vButtonClicked)
		self.oButton.clicked.connect(self.oButtonClicked)
		self.runButton.clicked.connect(self.runButtonClicked)
		self.gvButton.clicked.connect(self.gvButtonClicked)
		self.vvButton.clicked.connect(self.vvButtonClicked)
		#######################################################
		self.setLayout(self.grid)
		#######################################################
	def traitCombo_chosen(self, text):
		if text == 'quantitative trait':
			self.vLab.setEnabled(True)
			self.vEdit.setEnabled(True)
			self.vButton.setEnabled(True)
			self.covLab.setEnabled(True)
			self.covEdit.setEnabled(True)
			self.testCombo.clear()
			self.testCombo.addItem("Linear regression")
			self.modelLab.setEnabled(False)
			self.modelCombo.setEnabled(False)
			self.covLab.setEnabled(True)
			self.covEdit.setEnabled(True)
		else:
			self.vLab.setEnabled(False)
			self.vEdit.setEnabled(False)
			self.vButton.setEnabled(False)
			self.vvButton.setEnabled(False)
			self.covLab.setEnabled(False)
			self.covEdit.setEnabled(False)
			self.testCombo.clear()
			self.testCombo.addItem("Pearson chi-squared test")
			self.testCombo.addItem("Fisher's exact test")
			self.testCombo.addItem("Logistic regression")
			self.testCombo.addItem("Raw Test")
			self.testCombo.addItem("Score Test")
			if str(self.testCombo.currentText()) == 'Logistic regression':
				self.vLab.setEnabled(True)
				self.vEdit.setEnabled(True)
				self.vButton.setEnabled(True)
				self.vvButton.setEnabled(True)
				self.covLab.setEnabled(True)
				self.covEdit.setEnabled(True)
			if str(self.testCombo.currentText()) == "Pearson chi-squared test" or str(self.testCombo.currentText()) == "Fisher's exact test":
				self.modelLab.setEnabled(True)
				self.modelCombo.setEnabled(True)
	def testCombo_chosen(self, text):
		if text == 'Logistic regression' or text == 'Linear regression':
			self.modelLab.setEnabled(False)
			self.modelCombo.setEnabled(False)
			self.vLab.setEnabled(True)
			self.vEdit.setEnabled(True)
			self.vButton.setEnabled(True)
			self.covLab.setEnabled(True)
			self.covEdit.setEnabled(True)
		else:
			self.modelLab.setEnabled(True)
			self.modelCombo.setEnabled(True)
			self.vLab.setEnabled(False)
			self.vEdit.setEnabled(False)
			self.vButton.setEnabled(False)
			self.vvButton.setEnabled(False)
			self.covLab.setEnabled(False)
			self.covEdit.setEnabled(False)	
		if text == 'Raw Test' or text == 'Score Test':
			self.modelLab.setEnabled(False)
			self.modelCombo.setEnabled(False)
	def gButtonClicked(self):
		self.gfile = QtGui.QFileDialog.getOpenFileName(self,'Open File', '.')
		self.gEdit.setText(self.gfile)
		if self.gfile:
			self.gvButton.setEnabled(True)
	def vButtonClicked(self):
		self.covfile = QtGui.QFileDialog.getOpenFileName(self,'Open File', '.')
		self.vEdit.setText(self.covfile)
		if self.covfile:
			self.vvButton.setEnabled(True)
	def oButtonClicked(self):
		self.outfile = QtGui.QFileDialog.getSaveFileName(self,'Open File', '.')
		self.oEdit.setText(self.outfile)
	def gvButtonClicked(self):
		self.gdialog = QtGui.QDialog()
		self.gdialog.resize(self.width, self.height)
		self.gdialog.setWindowTitle("view of genotype data")
		df  = pd.read_csv(str(self.gfile), delim_whitespace= True, index_col = None, header = None)
		gdatatable = QtGui.QTableWidget(parent=self.gdialog)
		gdatatable.setColumnCount(len(df.columns))
		gdatatable.setRowCount(len(df.index))
		labels = list(df.iloc[0])
		labels = labels[2:]
		lheader = ['IID', 'PHT']
		lindex = 0
		for label in labels:
			temp = label.split('*')
			tlabel = 'HLA-' + temp[0]
			if lindex == 0:
				tlabel += '-1'
				lindex = 1
			else:
				tlabel += '-2'
				lindex = 0
			lheader.append(tlabel)
		gdatatable.setHorizontalHeaderLabels(lheader)
		for i in range(len(df.index)):
			for j in range(len(df.columns)):
				gdatatable.setItem(i,j,QtGui.QTableWidgetItem(str(df.iget_value(i, j))))
		gdatatable.resizeColumnsToContents()
		ghbox = QtGui.QHBoxLayout()
		ghbox.addWidget(gdatatable)
		self.gdialog.setLayout(ghbox)
		self.gdialog.show()
	def vvButtonClicked(self):
		self.vdialog = QtGui.QDialog()
		self.vdialog.resize(self.width, self.height)
		self.vdialog.setWindowTitle("view of covariates data")
		df  = pd.read_csv(str(self.covfile), delim_whitespace= True, index_col = None, header = 0)
		vdatatable = QtGui.QTableWidget(parent=self.vdialog)
		vdatatable.setColumnCount(len(df.columns))
		vdatatable.setRowCount(len(df.index))
		vdatatable.setHorizontalHeaderLabels(list(df.columns.values))
		for i in range(len(df.index)):
			for j in range(len(df.columns)):
				vdatatable.setItem(i,j,QtGui.QTableWidgetItem(str(df.iget_value(i, j))))
		vdatatable.resizeColumnsToContents()

		vhbox = QtGui.QHBoxLayout()
		vhbox.addWidget(vdatatable)
		self.vdialog.setLayout(vhbox)
		self.vdialog.show()
	def runButtonClicked(self):
		comm = 'python HLAassoc.py --file '
		comm += str(self.gfile)
		digi = str(self.digitCombo.currentText())
		comm += ' --digits '
		comm += digi
		testStr = str(self.testCombo.currentText())
		testM = ''
		if testStr == 'Pearson chi-squared test':
			testM = 'chisq'
		elif testStr == "Fisher's exact test":
			testM = 'fisher'
		elif testStr == 'Logistic regression':
			testM = 'logistic'
		elif testStr == 'Linear regression':
			testM = 'linear'
		elif testStr == 'Raw Test':
			testM = 'raw'
		elif testStr == 'Score Test':
			testM = 'score'
		comm += ' --test '
		comm += testM
		if testM == 'chisq' or testM == 'fisher':
			model = str(self.modelCombo.currentText())
			comm += ' --model '
			comm += model
		freq = str(self.frqEdit.text())
		comm += ' --freq '
		comm += freq
		adjust = str(self.adjCombo.currentText())
		if testM != 'score' and testM != 'raw':
			comm += ' --adjust '
			comm += adjust
		if testM == 'logistic' or testM == 'linear':
			if self.covfile:
				comm += ' --covar '
				comm += str(self.covfile)
				if str(self.covEdit.text()):
					comm += ' --covarname '
					comm += str(self.covEdit.text())
		if str(self.permEdit.text()) and int(str(self.permEdit.text())) > 0:
			comm += ' --perm '
			comm += str(self.permEdit.text())
		comm += ' --out '
		comm += str(self.outfile)
		os.system(comm)
		# show results
		self.odialog = QtGui.QDialog()
		self.odialog.resize(self.width, self.height)
		self.odialog.setWindowTitle("view of output")
		df  = pd.read_csv(str(self.outfile), delim_whitespace= True, index_col = None, header = 0)
		odatatable = QtGui.QTableWidget(parent=self.odialog)
		odatatable.setColumnCount(len(df.columns))
		odatatable.setRowCount(len(df.index))
		odatatable.setHorizontalHeaderLabels(list(df.columns.values))
		for i in range(len(df.index)):
			for j in range(len(df.columns)):
				odatatable.setItem(i,j,QtGui.QTableWidgetItem(str(df.iget_value(i, j))))
		odatatable.resizeColumnsToContents()
		ohbox = QtGui.QHBoxLayout()
		ohbox.addWidget(odatatable)
		self.odialog.setLayout(ohbox)
		self.odialog.show()
#####################################################

if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	gui = HLAassocWin()
	gui.show()
	sys.exit(app.exec_())

