import string
import sys
import os
import math


# Some constants
binInterval = 20  			# The size of bin (20 millisecond)

forceStartMessage = 0;			# Whether the output should start from a specific message (0).
					# If this is set to 1, the analysis will set the time of the message as the start time of the plot.
					# If you have already created the output report by using the interest periods in Data Viewer, set
					# this to 0.

syncMessage = "SYNCTIME"		# If forceStartMessage = 1, need to specify the start message; no data output will be available 
					# if the syncMessage is not found in the file, 
						
includeSamplesInFixationOnly = 0;	# Whether we should consider those samples in a fixation only (1); this option should be set to 0 if samples 
					# during saccade should also be considered in the report, .

excludeEntireBinIfInSaccade = 0;	# Whether the entire bin should be discarded if some samples in the interval are in saccade.
					# For a fixation-only analysis, this option should be 1.

maxNumberofIAs = 4;			# Maximum number of interest areas used in each trial (ID can be 1, 2, 3, 4).

forceOutsideIATo0 = 1			# If the eye sample is not within any interest areas; force it to being ID = 0
					# Otherwise, the sample is simply ignored.

collapseAcrossTrials = 0		# Collapse the frequency of fixations to each of the interest area across trials;
					# IMPORTANT: THIS OPTION IS ONLY AVAILABLE IF THE INTEREST AREA ID IN EACH OF THE TRIALS
					# CORRESPONDS TO THE SAME THING (e.g., ID = 1 refers to the target across all of the trial,
					# ID = 2 refers to the semantic distractor, etc).  
					# If the interest area ID refers to different items across trials, it does not make any 
					# sense to perform this collapsing analysis.  You will need to redo the ID of the interest 
					# areas in Data Viewer before performing this analysis.  
					# Alternatively,use the following "recodeInterestAreaBasedOnLabels" option to recode the ID of interest areas. 

maxBinSize = 2000;			# Maximum number of bin size in the collapse analysis.  This corresponds to
					# maxBinSize * binInterval millseconds of data in each trial (2000* 20 ms of recording after the sync message)

groupingVariableLabel = "target"	# Label of the column you planned to do trial grouping; the data in the same condition will be collapsed.
					# Note that this label is case-sensitive
					# The trial grouping will be ignored if the variable is not among the list.


recodeInterestAreaBasedOnLabels =  0	# Check whether the interest area ID used in the bin output is organized by the "LEFT_INTEREST_AREA_ID" 
					# and/or "RIGHT_INTEREST_AREA_ID".
					# If you have a consistent interest area mapping (e.g., target is always mapped to ID 1 of the interest areas,
					# and distract_1 is mapped to 2, etc), you don't need to change the default settings (0).
					# If you code the interest areas based on the labels (although the ID can vary across trials), turn 
					# this to 1 and additionally recode the lines following line #656 "if recodeInterestAreaBasedOnLabels:"


#########  Check data content;  #########
def isInteger(astring):
	try: int(astring)
	except ValueError: return 0
	else: return 1


def isLong(astring):
	try: long(astring)
	except ValueError: return 0
	else: return 1
	
	
def isFloat(astring):
	try: float(astring)
	except ValueError: return 0
	else: return 1
	
	
def isMissingValue(astring):
	if astring == '.':
		return 1
	else:
		return 0;
		


#########  Structure of the bin  #########
class BinData:
	def __init__(self, startTime=0):
		self.sampleCount = 0;		# Total number of samples in the bin
		self.iaLeftEye = list();	# list used to contain fixation frequency to each of the interest area
		self.iaRightEye = list();
		self.timeBinEnd = 0;		# Time of the start sample in the bin
		self.timeBinStart = startTime;	# Time of the end sample in the bin
		
		self.hasLeftSaccadeOrBlinkOrNoData = 0;		# Whether there is a saccade, blink, or no data within the current bin;
		self.hasRightSaccadeOrBlinkOrNoData = 0;
		
		for i in range(maxNumberofIAs+1):
			self.iaLeftEye.append(0);
			self.iaRightEye.append(0);	


#########  Structure of the data file  #########
class DataFile:
	def __init__(self, fileString):
		self.fileName = fileString	# Input data file
		self.inputFile = None
		self.outputFile = None;
		
		self.hasIALabel = 0;		# Check whether the input file contains the labels for the interest areas
		self.iaLabels = list();		# Variable used to contain the interest area labels
		self.hasRecordingSessionLabel = 0;

		self.variableLabelPrinted = 0;	# Check whether the label of the columns has been printed to the output file	
		self.variableLabels = list();	# List used to contain the labels of the data columns 
		self.myDict = dict();		# Dictionary for mapping the data column label and value.
		
		self.trialStartTime = 0;	# Start time of the trial
		self.previousTrial = None;	# Label of the previous trial
		self.previousSession = None 	# Label of the previous session
		self.currentBin = 0;		# Current bin of the samples in a trial; 
		self.previousBin= 0;
		
		self.eyeUsed = -1;		# Check whether the data file contains the left or right eye information.
		self.currentLine =0;		
		
		self.currentCondition ="null"	# the value of the current trial condition; this is recorded for trial grouping later
		
		self.binAcrossTrialsLeftEye = [] # trial data for each of the bin
		self.binAcrossTrialsRightEye = []
			
		self.performTrialGrouping = 0;	# Whether data collapsing should be based on the different groups of trials;
		self.trialCondition = list();	# List used to store the condition information for each trial	
		
		self.extraUserDefinedTrialConditionLabels = list();
		self.extraUserDefinedTrialConditionValues = list();
		self.trialConditionExtraUserDefinedVariables = list();
		# Todo; append your own variables;
		self.extraUserDefinedTrialConditionLabels.append('List');
		self.extraUserDefinedTrialConditionLabels.append('verb');
		self.extraUserDefinedTrialConditionLabels.append('cond');
			
			
	### Open the input file and create an output file according to each trial (no collapsing of the data across trials)
	def openFiles(self):
		self.inputFile = open(self.fileName, 'r')
		if self.inputFile == None:
			print "cannot open file ", self.fileName
			return 1;
	
		baseFilename = os.path.splitext(os.path.basename(self.fileName));
		self.outputFile = open(baseFilename[0]+'_trial.xls', 'wt')
		if self.outputFile == None:
			print "cannot open file ", self.fileName
			return 1; 	
		
		return 0;
		
		
	### Check the content of the input file
	def checkInputDataFile(self, line):
		print "1.\tBin interval ", binInterval
		print "2.\tmaxBinSize ", maxBinSize, "(Maximum analysis duration for each trial ", maxBinSize*binInterval, ")"

		if line.find('RECORDING_SESSION_LABEL') == -1:
			print "Input data doesn't contain the 'RECORDING_SESSION_LABEL' column!"
			return 1;

		if line.find('TRIAL_LABEL') == -1:
			print "Input data doesn't contain the 'TRIAL_LABEL' column!"
			return 1;

		if line.find('TIMESTAMP') == -1:
			print "Input data doesn't contain the 'TIMESTAMP' column!"
			return 1;

		if forceStartMessage and line.find('SAMPLE_MESSAGE') == -1:
			print "Input data doesn't contain the 'SAMPLE_MESSAGE' column!"
			return 1;

		# Check the eye and interest area label information
		self.hasIALabel = 0;
		if line.find('LEFT_INTEREST_AREA_ID') > -1 and line.find('RIGHT_INTEREST_AREA_ID') > -1:
			self.eyeUsed = 2
			if line.find('LEFT_INTEREST_AREA_LABEL') > -1 and line.find('RIGHT_INTEREST_AREA_LABEL') > -1:
				self.hasIALabel = 1;
		elif line.find('RIGHT_INTEREST_AREA_ID') > -1:
			self.eyeUsed = 1
			if line.find('RIGHT_INTEREST_AREA_LABEL') > -1:
				self.hasIALabel = 1;
		elif line.find('LEFT_INTEREST_AREA_ID'):
			self.eyeUsed = 0		
			if line.find('LEFT_INTEREST_AREA_LABEL') > -1:
				self.hasIALabel = 1;			
		else:
			print "Input data doesn't contain the 'LEFT_INTEREST_AREA_ID' or 'RIGHT_INTEREST_AREA_ID'column!"
			return 1;						

		# Check whether the report should consider those samples within a fixation only
		# If so, the input file needs to contain those columns indicating whether the sample is in a saccade or blink
		if includeSamplesInFixationOnly:
			print "Ignore the samples during blink or saccade";

			if self.eyeUsed == 0 or self.eyeUsed == 2:
				if line.find('LEFT_IN_BLINK') == -1 or line.find('LEFT_IN_SACCADE') == -1:
					print "Input data doesn't contain the 'LEFT_IN_BLINK' or 'LEFT_IN_SACCADE' column!"
					return 1;

			if self.eyeUsed == 1 or self.eyeUsed == 2:
				if line.find('RIGHT_IN_BLINK') == -1 or line.find('RIGHT_IN_SACCADE') == -1:
					print "Input data doesn't contain the 'RIGHT_IN_BLINK' or 'RIGHT_IN_SACCADE' column!"
					return 1;
					
			
		# If need to collapse across trials, we need to allocate buffer
		if collapseAcrossTrials:
			for i in range(maxBinSize +1):
				self.binAcrossTrialsLeftEye.append(list());
				self.binAcrossTrialsRightEye.append(list());
		
		# If need to group the trials based on one variable, need to check that
		if line.find(groupingVariableLabel) == -1:
			print "3.\tCannot find the grouping variable '", groupingVariableLabel, "'\nNo trial grouping is performed!"
			self.performTrialGrouping = 0;
		else:
			print "3.\tTrial grouping is based on the variable '", groupingVariableLabel, "'"
			self.performTrialGrouping = 1;

		# Check whether variables are actually in the list.
		if len(self.extraUserDefinedTrialConditionLabels) >0:
			length = len(self.extraUserDefinedTrialConditionLabels)
			for i in range(length):
				variable = self.extraUserDefinedTrialConditionLabels[length-1 - i]
				if line.find(variable) == -1:
					print "******* ", variable, " not found; dropped for user-defined variable list"
					self.extraUserDefinedTrialConditionLabels.remove(variable)
		print "4.\tExtra User Defined Trial Condition Labels in the trial-based report "
		print "\t", self.extraUserDefinedTrialConditionLabels 
		
		print "5.\tincludeSamplesInFixationOnly = ", includeSamplesInFixationOnly
		print "\texcludeEntireBinIfInSaccade = ", excludeEntireBinIfInSaccade
		if excludeEntireBinIfInSaccade:
			print "\tEntire bin discarded if any samples in the interval is in \n\tsaccade or blink"
		else:	
			print "\tSome samples still kept even if any samples in the interval \n\tis in saccade or blink (default settings)"
		
		
		print "6.\trecodeInterestAreaBasedOnLabels ", recodeInterestAreaBasedOnLabels
		if recodeInterestAreaBasedOnLabels:
			print "\tRecode the Index of the interest areas based on the interest area labels";
		else:
			print "\tIndex of the interest areas are based on the \n\tLEFT_INTEREST_AREA_ID/RIGHT_INTEREST_AREA_ID (default settings)";
		
		
		print "7.\tMaximum number of interest areas used in each trial ", maxNumberofIAs
		
		print "8.\tforceOutsideIATo0 = ", forceOutsideIATo0
		if forceOutsideIATo0:
			print "\tIf the eye sample is not within any interest areas; \n\tforce it to being ID = 0 (default settings)"
		else:
			print "\tIf the eye sample is not within any interest areas, \n\tthat sample is ignored and dropped from the analysis"
		
		
		print "9.\tforceStartMessage = ", forceStartMessage
		if forceStartMessage:
			print "\tThe first bin starts on the onset of message ", syncMessage
		else:
			print "\tThe first bin of each trial starts on the first sample in the trial"
			print "\t(this assumes that you have created the interest period and created"
			print "\tthe sample report for this period; default settings)"
		
		print "-------------------------------------------------------------------\n"
		return 0;						
	
	
	### Outputs the bin analysis to the files. 
	def binReportToOutput(self, currentSampleBin, sessionLabel, trialLabel, condition):
		# Check whether it has exceeds the bin size.
		if self.previousBin >= maxBinSize:
			print "\tSkip data output because the maximum bin size has been reached"
			return 0;
	
		# Check whether the column labels have been printed out for the trials 		
		if not self.variableLabelPrinted:
			self.variableLabelPrinted = 1;
			text = "RECORDING_SESSION_LABEL\tTRIAL_LABEL\tTRIAL_CONDITION\tTRIAL_START_TIME\t";
			text +="CURRENT_BIN\tBIN_START_TIME\tBIN_END_TIME\tSAMPLE_COUNT\t"
			
			# Print out the left eye data
			if self.eyeUsed == 0 or self.eyeUsed == 2:
				for i in range(maxNumberofIAs +1):
					text +="LEFT_%d_C\t"%(i)
				for i in range(maxNumberofIAs +1):
					text +="LEFT_%d_P\t"%(i)
			
			# Print out the right eye data
			if self.eyeUsed == 1 or self.eyeUsed == 2:
				for i in range(maxNumberofIAs +1):
					text +="RIGHT_%d_C\t"%(i)
				for i in range(maxNumberofIAs +1):
					text +="RIGHT_%d_P\t"%(i)
	
			# Print out the labels for the columns
			for i in range(maxNumberofIAs +1):
				text +="IA_LABEL_%d\t"%(i)
				
			for label in self.extraUserDefinedTrialConditionLabels:
				text += label + "\t"
			self.outputFile.write(text)


		text  ="\n%s\t%s\t%s\t%d\t"%(str(sessionLabel), str(trialLabel), str(self.currentCondition), self.trialStartTime)
		text +="%d\t%d\t%d\t"%(self.previousBin, currentSampleBin.timeBinStart, currentSampleBin.timeBinEnd);
		text +="%d\t"%(currentSampleBin.sampleCount);
	
		count = float(currentSampleBin.sampleCount)
	
		# Print out the left eye sample frequency to each interest area
		# and then calculate the proportion to each interest area 
			
	
		if self.eyeUsed == 0 or self.eyeUsed == 2:
			binDataFromTheCurrentTrial = [];
			if currentSampleBin.hasLeftSaccadeOrBlinkOrNoData == 1:
				for i in range(maxNumberofIAs +1):
					currentSampleBin.iaLeftEye[i] = 0;

			for i in range(maxNumberofIAs +1):
				text +="%d\t"%(currentSampleBin.iaLeftEye[i]);

			if count > 0:
				for i in range(maxNumberofIAs +1):
					text +="%.2f\t"%(currentSampleBin.iaLeftEye[i]/count);
					binDataFromTheCurrentTrial.append(currentSampleBin.iaLeftEye[i])
			else:
				for i in range(maxNumberofIAs +1):
					text +="%.2f\t"%(0.0);
					binDataFromTheCurrentTrial.append(0);
					
			if collapseAcrossTrials and self.previousBin < maxBinSize:
				binDataFromTheCurrentTrial.append(currentSampleBin.sampleCount*currentSampleBin.hasLeftSaccadeOrBlinkOrNoData)
				binDataFromTheCurrentTrial.append(condition)
				self.binAcrossTrialsLeftEye[self.previousBin].append(binDataFromTheCurrentTrial)
				#print "=", binDataFromTheCurrentTrial	
		
		if self.eyeUsed == 1 or self.eyeUsed == 2:
			binDataFromTheCurrentTrial = [];
			if currentSampleBin.hasRightSaccadeOrBlinkOrNoData == 1:
				for i in range(maxNumberofIAs +1):
					currentSampleBin.iaRightEye[i] = 0;

			for i in range(maxNumberofIAs +1):
				text +="%d\t"%(currentSampleBin.iaRightEye[i]);

			if count > 0:
				for i in range(maxNumberofIAs +1):
					text +="%.2f\t"%(currentSampleBin.iaRightEye[i]/count);
					binDataFromTheCurrentTrial.append(currentSampleBin.iaRightEye[i])
			else:
				for i in range(maxNumberofIAs +1):
					text +="%.2f\t"%(0.0);
					binDataFromTheCurrentTrial.append(0);
					
			if collapseAcrossTrials and self.previousBin < maxBinSize:
				binDataFromTheCurrentTrial.append(currentSampleBin.sampleCount*currentSampleBin.hasRightSaccadeOrBlinkOrNoData)
				binDataFromTheCurrentTrial.append(condition)
				self.binAcrossTrialsRightEye[self.previousBin].append(binDataFromTheCurrentTrial)
				#print "=", binDataFromTheCurrentTrial			
							
		# Print out the label for each interest areas
		for i in range(maxNumberofIAs +1):
			text +="%s\t"%(str(self.iaLabels[i]));
	
		for value in self.extraUserDefinedTrialConditionValues:
			text += str(value) + "\t"
	
		self.outputFile.write(text)
		return 0



	### Outputs the bin analysis to the files. 
	def collapseTrialOutput(self):
		# Check how many conditions we have;
		self.trialGroupLabel = list();
		for i in range(len(self.trialCondition)):
			if not (self.trialCondition[i] in self.trialGroupLabel):
				self.trialGroupLabel.append(self.trialCondition[i])
	
		print "\n\nCollapsing data across trials"
		#Going through each of the groups;
		for group in self.trialGroupLabel:
			if self.performTrialGrouping:
				print groupingVariableLabel, "->", group
				
			baseFilename = os.path.splitext(os.path.basename(self.fileName));
			fileName = baseFilename[0]+'_collapsed_' +  str(group) + '.xls'
			self.collapseOutputFile = open(fileName, 'wt')
			if self.collapseOutputFile == None:
				print "cannot open file ", fileName
				return 1; 	

			#Print out the column headings
			text = "RECORDING_SESSION_LABEL\tCURRENT_BIN\tGROUP\tTRIAL_COUNT\t";

			# Print out the left eye data
			if self.eyeUsed == 0 or self.eyeUsed == 2:
				for i in range(maxNumberofIAs +1):
					text +="LEFT_%d_C\t"%(i)

				text +="LEFT_SACCADE_COUNT\tLEFT_VALID_COUNT\t"

				for i in range(maxNumberofIAs +1):
					text +="LEFT_%d_P\t"%(i)

			# Print out the right eye data
			if self.eyeUsed == 1 or self.eyeUsed == 2:
				for i in range(maxNumberofIAs +1):
					text +="RIGHT_%d_C\t"%(i)

				text +="RIGHT_SACCADE_COUNT\tRIGHT_VALID_COUNT\t"

				for i in range(maxNumberofIAs +1):
					text +="RIGHT_%d_P\t"%(i)

			# Print out the labels for the columns
			for i in range(maxNumberofIAs +1):
				text +="IA_LABEL_%d\t"%(i)

			# This output is actually not applicable to the collapsed trial reports
			#for label in self.extraUserDefinedTrialConditionLabels:
			#	text += label + "\t"

			self.collapseOutputFile.write(text)
			
			
			# Going through each of the data bin;
			for individualBin in range(len(self.binAcrossTrialsRightEye)):				
				# Counts the number of trials in the group;
				leftGroupTrialCount = 0;
				trialCount = len(self.binAcrossTrialsLeftEye[individualBin]) 
				if trialCount> 0:
					for trial in range(trialCount):
						if  self.binAcrossTrialsLeftEye[individualBin][trial][maxNumberofIAs+2] == group:
							leftGroupTrialCount += 1;
				rightGroupTrialCount = 0;
				trialCount = len(self.binAcrossTrialsRightEye[individualBin]) 
				if trialCount> 0:
					for trial in range(trialCount):
						if  self.binAcrossTrialsRightEye[individualBin][trial][maxNumberofIAs+2] == group:
							rightGroupTrialCount += 1;
			
				#Skip the output if no data is available
				groupTrialCount = max(leftGroupTrialCount, rightGroupTrialCount)
				if groupTrialCount == 0:
					break;
				
				text  ="\n%s\t"%(str(self.myDict['RECORDING_SESSION_LABEL']));
				text +="%d\t"%(individualBin);
				text +="%s\t"%(str(group));
				text +="%d\t"%(groupTrialCount)

				# If the file contains the left-eye data
				if self.eyeUsed == 0 or self.eyeUsed == 2:
					iaDistribution = list();
					total = 0 
					for ia in range(maxNumberofIAs+2):	# ID = 0, 1, 2, 3, 4, and saccade/blink flag
						iaDistribution.append(0);

						trialCount = len(self.binAcrossTrialsLeftEye[individualBin]) 
						if trialCount> 0:
							for trial in range(trialCount):
								# Check whether the label of the trial matches the label of the group
								if  self.binAcrossTrialsLeftEye[individualBin][trial][maxNumberofIAs+2] == group:
									iaDistribution[ia] += self.binAcrossTrialsLeftEye[individualBin][trial][ia]
							
							if ia < maxNumberofIAs+1:
								total += iaDistribution[ia]
						else:
							iaDistribution[ia] = 0;
						text +="%.2f\t"%(iaDistribution[ia])

					text +="%d\t"%(total)

					for ia in range(maxNumberofIAs+1):
						if total > 0:
							iaDistribution[ia] = (iaDistribution[ia]*100.0)/float(total)
						else:
							iaDistribution[ia] = 0;
						text +="%.2f\t"%(iaDistribution[ia])			
				
				# If the file contains the right-eye data
				if self.eyeUsed == 1 or self.eyeUsed == 2:
					iaDistribution = list();
					total = 0 
					for ia in range(maxNumberofIAs+2):
						iaDistribution.append(0);

						trialCount = len(self.binAcrossTrialsRightEye[individualBin]) 
						if trialCount> 0:
							for trial in range(trialCount):
								# Check whether the label of the trial matches the label of the group
								if  self.binAcrossTrialsRightEye[individualBin][trial][maxNumberofIAs+2] == group:
									iaDistribution[ia] += self.binAcrossTrialsRightEye[individualBin][trial][ia]

							if ia < maxNumberofIAs+1:
								total += iaDistribution[ia]
						else:
							iaDistribution[ia] = 0;
						text +="%.2f\t"%(iaDistribution[ia])

					text +="%d\t"%(total)

					for ia in range(maxNumberofIAs+1):
						if total > 0:
							iaDistribution[ia] = (iaDistribution[ia]*100.0)/float(total)
						else:
							iaDistribution[ia] = 0;
						text +="%.2f\t"%(iaDistribution[ia])			


				# Print out the label for each interest areas
				for i in range(maxNumberofIAs +1):
					text +="%s\t"%(str(self.iaLabels[i]));

				#This output is actually not applicable to the collapsed trial reports.
				#index = self.trialCondition.index(group)
				#for value in self.trialConditionExtraUserDefinedVariables[index]:
				#	text += str(value) + "\t"

				self.collapseOutputFile.write(text)		

			# Close the file at the end;
			self.collapseOutputFile.close();
		# End of writing files;
		return 0;



	### Method used to perform data analysis 
	def extractData(self):
		# Open the input data file
		if self.openFiles():
			return 1;
	
		# Initialize the sample bin
		currentSampleBin = BinData();
		
		if int(binInterval) <= 0:
			print "Error in the bin interval value ", binInterval
			return 1;			
				
				
		######### Read the content of the TEXT file #########
		while True:
			line=self.inputFile.readline()
			if len(line) <= 0:
				break;
			else:
				if "\n" in line:
					line = line.strip('\n')
				stringCache = line.split("\t");
				self.currentLine += 1;
				
				# If the first line, we need to check the column headings.
				if self.currentLine == 1:					
					if self.checkInputDataFile(line) == 1:
						return 1
					
					for label in stringCache: # Remove quotes in the labels
						if "\n" in label:
							label = label.strip('\n')
						label = label.strip()
						self.variableLabels.append(label.strip('"'))
					
				else:
					# Read in the content of for each data line;
					# Create a dictionary for the data line
			
					# Populate the value of for each column variable
					for i in range(len(self.variableLabels)):
						
						# Code used to convert "3,12" to "3.12" for some regional settings
						# (Control Panel -> Date, Time, Lanugage, and Regional Options  -> Regional and Language Options -> Customize -> Decimal Symbol)
						temp = stringCache[i]
						if "," in temp:
							stringCache[i] = temp.replace(",", ".");
					
						if isInteger( stringCache[i]):	
							self.myDict[self.variableLabels[i]] = int(stringCache[i])
						elif isLong(stringCache[i]):
							self.myDict[self.variableLabels[i]] = long(stringCache[i])
						elif isFloat(stringCache[i]):
							self.myDict[self.variableLabels[i]] = float(stringCache[i])
						else:	
							self.myDict[self.variableLabels[i]] = stringCache[i]
					
					
					# Is this a new trial?
					if self.myDict['TRIAL_LABEL'] != self.previousTrial:
						print "Processing ", self.myDict['RECORDING_SESSION_LABEL'], "  ", self.myDict['TRIAL_LABEL']
						#Write out the last bin for the previous trial;
						if self.previousTrial and currentSampleBin.sampleCount > 0:
							self.binReportToOutput(currentSampleBin, self.previousSession, self.previousTrial, self.currentCondition);
						
						self.iaLabels = list();
						for i in range(maxNumberofIAs+1):
							self.iaLabels.append(".");
						
						#Reset the trial data;
						currentSampleBin = BinData(long(self.myDict['TIMESTAMP']));
						if forceStartMessage == 1:
							self.trialStartTime = 0;
						else:
							self.trialStartTime = long(self.myDict['TIMESTAMP'])
						
						self.previousTrial = self.myDict['TRIAL_LABEL'] 
						self.previousSession = self.myDict['RECORDING_SESSION_LABEL'];
					
						# Need to check the trial condition value if we need to do the trial grouping.	
						if self.performTrialGrouping:
							self.currentCondition = self.myDict[groupingVariableLabel]
							print "\tgrouping variable :", groupingVariableLabel, " -> ", self.myDict[groupingVariableLabel]
						else:
							self.currentCondition = "null"
						self.trialCondition.append(self.currentCondition);

						self.extraUserDefinedTrialConditionValues = list();
						for label in self.extraUserDefinedTrialConditionLabels:
							self.extraUserDefinedTrialConditionValues.append(self.myDict[label])
						
						if len(self.extraUserDefinedTrialConditionValues):
							print "\t",self.extraUserDefinedTrialConditionLabels
							print "\t",self.extraUserDefinedTrialConditionValues
						
						self.trialConditionExtraUserDefinedVariables.append(self.extraUserDefinedTrialConditionValues)
							

					#Check for the trial start message for each trial;
					#The message is simply ignored if you do not have to force checking the start message.
					if forceStartMessage and self.myDict['SAMPLE_MESSAGE'] != '.' and self.myDict['SAMPLE_MESSAGE'].find(syncMessage) != -1:
						self.trialStartTime = long (self.myDict['TIMESTAMP'])
						currentSampleBin = BinData(long(self.myDict['TIMESTAMP']));


					# Check whether the samples are within the interest period;
					# (ie., after the sync message if forceStartMessage is on)
					if self.trialStartTime > 0:
						self.currentBin =long(self.myDict['TIMESTAMP'] - self.trialStartTime)/int(binInterval)
						#print self.myDict['TRIAL_LABEL'],  "  ->  ", self.myDict['TIMESTAMP'], " -> ", self.currentBin
						
						#Output the data if a new bin starts
						if self.currentBin > self.previousBin: 
							self.binReportToOutput(currentSampleBin, self.previousSession, self.myDict['TRIAL_LABEL'], self.currentCondition);
							currentSampleBin = BinData(long(self.myDict['TIMESTAMP']));
							
						currentSampleBin.timeBinEnd  = long(self.myDict['TIMESTAMP'])
						currentSampleBin.sampleCount += 1;			
						self.previousBin = self.currentBin
						
						# The following code segment is only useful if you need to recode the index of the
						# interest areas based on the label of the interest areas "LEFT_INTEREST_AREA_LABEL"
						# and/or 'RIGHT_INTEREST_AREA_ID'
						if recodeInterestAreaBasedOnLabels:
							if self.myDict['LEFT_INTEREST_AREA_ID'] != '.':
								temp = self.myDict['LEFT_INTEREST_AREA_LABEL']
								temp = temp.strip()

								if temp == 'target':
									self.myDict['LEFT_INTEREST_AREA_ID'] = 1
								elif temp == '2nd_targ':
									self.myDict['LEFT_INTEREST_AREA_ID'] = 2
								elif temp == 'comp_1':
									self.myDict['LEFT_INTEREST_AREA_ID'] = 3
								elif temp == 'comp_2':
									self.myDict['LEFT_INTEREST_AREA_ID'] = 4
								elif temp == 'fix':
									self.myDict['LEFT_INTEREST_AREA_ID'] = 5
								else:
									self.myDict['LEFT_INTEREST_AREA_ID'] = 0;

							if self.myDict['RIGHT_INTEREST_AREA_ID'] != '.':
								temp = self.myDict['RIGHT_INTEREST_AREA_LABEL']
								temp = temp.strip()

								if temp == 'target':
									self.myDict['RIGHT_INTEREST_AREA_ID'] = 1
								elif temp == '2nd_targ':
									self.myDict['RIGHT_INTEREST_AREA_ID'] = 2
								elif temp == 'comp_1':
									self.myDict['RIGHT_INTEREST_AREA_ID'] = 3
								elif temp == 'comp_2':
									self.myDict['RIGHT_INTEREST_AREA_ID'] = 4
								elif temp == 'fix':
									self.myDict['RIGHT_INTEREST_AREA_ID'] = 5
								else:
									self.myDict['RIGHT_INTEREST_AREA_ID'] = 0;
								
							
						# Now tally samples within each interest area.
						if self.eyeUsed == 0 or self.eyeUsed == 2:
							# If the eye sample is not in any of the interest areas; do we need to convert
							# it to 0? This is useful if for example you just define the foreground interest area
							# You can then treat the rest samples as fixations to the background 
							if forceOutsideIATo0 and self.myDict['LEFT_INTEREST_AREA_ID'] == '.':
								self.myDict['LEFT_INTEREST_AREA_ID'] = 0;
						
							if self.myDict['LEFT_INTEREST_AREA_ID'] != '.':
								currentID = int(self.myDict['LEFT_INTEREST_AREA_ID'])
								if  currentID > maxNumberofIAs:
									print "error in the LEFT interest area ID value (%d), which must be within %d" % (currentID, maxNumberofIAs)
									return 1;
								else:
									if excludeEntireBinIfInSaccade and (self.myDict['LEFT_IN_BLINK'] or self.myDict["LEFT_IN_SACCADE"]):
										currentSampleBin.hasLeftSaccadeOrBlinkOrNoData = 1;
									
									# We determine whether we only includes those samples within a fixation only;
									# If so, we need to remove those samples during blink or during saccade.
									if includeSamplesInFixationOnly and (self.myDict['LEFT_IN_BLINK'] or self.myDict["LEFT_IN_SACCADE"]):
										pass;
									else:
										currentSampleBin.iaLeftEye[currentID] += 1;
										if self.hasIALabel == 1:
											self.iaLabels[currentID] = self.myDict['LEFT_INTEREST_AREA_LABEL']

						if self.eyeUsed == 1 or self.eyeUsed == 2:
							# If the eye sample is not in any of the interest areas; do we need to convert
							# it to 0? This is useful if for example you just define the foreground interest area
							# You can then treat the rest samples as fixations to the background 
							if forceOutsideIATo0 and self.myDict['RIGHT_INTEREST_AREA_ID'] == '.':
								self.myDict['RIGHT_INTEREST_AREA_ID'] = 0;

							if self.myDict['RIGHT_INTEREST_AREA_ID'] != '.':
								currentID = int(self.myDict['RIGHT_INTEREST_AREA_ID'])

								if  currentID> maxNumberofIAs:
									print "error in the RIGHT interest area ID value (%d), which must be within %d" %(currentID, maxNumberofIAs)
									return 1;
								else:
									if excludeEntireBinIfInSaccade and (self.myDict['RIGHT_IN_BLINK'] or self.myDict["RIGHT_IN_SACCADE"]):
										currentSampleBin.hasRightSaccadeOrBlinkOrNoData = 1;
									
									if includeSamplesInFixationOnly and (self.myDict['RIGHT_IN_BLINK'] or self.myDict["RIGHT_IN_SACCADE"]):
										pass;
									else:
										currentSampleBin.iaRightEye[currentID] += 1;
										if self.hasIALabel == 1:
											self.iaLabels[currentID] = self.myDict['RIGHT_INTEREST_AREA_LABEL']


				#print self.currentLine, self.myDict['TRIAL_LABEL'], self.trialStartTime, "bin ", self.currentBin


		# end of .txt file.
		if self.previousTrial and currentSampleBin.sampleCount > 0:
			self.binReportToOutput(currentSampleBin, self.previousSession, self.previousTrial, self.currentCondition);
		
		# If each of the interest area ID corresponds to the same thing,
		# You may collapse the data across trials;
		if collapseAcrossTrials:
			self.collapseTrialOutput();
		
		self.outputFile.write("\n");
		self.outputFile.close();
		self.inputFile.close()
		return 0;
		

print "\n"
print "** VisualWorld_YYMMDD.py                                            **"
print "** This analysis code aims to bin proportion of samples to interest **"
print "** areas within a fixed interval (e.g., in visual-world paradigm).  **"
print "** To run this code, type the following in the command Prompt       **"
print "** visualworld.py [input file name]                                 **"
print "** For example, 'visualworld.py sample.txt'                         **"
print "** The [input file name] is an output file from EyeLink Data Viewer **"
print "** Sample Report output ('Analysis -> Reports -> Sample Report')    **"
print "** It's recommended that you first create an interest period in     **"
print "** Data Viewer.  Following this, create the sample report           **"
print "** you may consider including all of the variables into the report. **"
print "** Some constants that can be modified at the beginning of the code **"
print "**                                                                  **"
print "** DISCLAIMER:                                                      **" 
print "** This example code was developed by an SR Research Applications   **"
print "** Engineer for technical support purposes. The code is supported   **"
print "** by SR Research Ltd., but may not have been completely verified   **"
print "** and tested.  Email support@sr-research.com to report bugs.       **"
print "**                                                                  **"
print "** HISTORY:                                                         **"
print "** mm/dd/yy                                                         **"
print "** 02/27/08 js  Initial version                                     **"
print "** 10/03/08 js  Supports trial grouping by one variable             **" 
print "** 02/05/09 js  Bug fixes                                           **" 
print "** 09/20/11 js  Allow to output custom variables                    **" 
print "** 09/20/11 js  Allow to recode the interest areas based on labels  **" 


######### Main #########
fileName = None;

if len(sys.argv) > 1:
	fileName = sys.argv[1];		
	print "\n-------------------------------------------------------------------"
	print "Input Data File: ", fileName, "\n\nConfigurations:"
else:
	print "\n-------------------------------------------------------------------"
	fileName = raw_input("Please enter the input file: ");
	print "Received input is : ", fileName

if fileName and os.path.exists(fileName) == True:
	data = DataFile(fileName)
	data.extractData();
else:
	print "Invalid Input File !!! ", fileName
	
