#!/usr/bin/env python
# This is the base of the regression tester
import imp, sys, os
import argparse
import subprocess

# Arguments
ap = argparse.ArgumentParser()
ap.add_argument('-v', '--verbose', help="When passing this flag, all printing will not be redirected to /dev/null.  This will create a lot of spew.", action="store_true", required=False, default=False)
args = ap.parse_args()


# Test Variables
tests = []
GSRevision = "60b9664"

# Open the file defining which tests to run
baseDir = "testsuites"
inpFile = "regression_tests.inp"
with open("%s/%s"%(baseDir,inpFile),"r") as f:
	for l in f:
		# ignore comments
		if l[0] == "#":
			continue
		ls = l.rstrip().split("\t")
		tests.append((ls[0], ls[1]))

# Print frontmatter
str1 = " | Lattice Microbes Regression Tester |"
rev = -1
try:
    p = subprocess.Popen(["bash","-c","git log --abbrev-commit --pretty=oneline | head -1 | awk '{print $1}'"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    rev = out.decode("utf-8").rstrip().split("\n")[0]
except:
	pass
str2 = " Git Revision:           %s"%(rev)
str3 = " Gold Standard Revision: %s"%(GSRevision)

print("")
print("\\"*len(str1))
print(str1)
print("/"*len(str1))
print(str2)
print(str3)
print("")
print("")



# Run tests
totalTests  = 0
passedTests = 0

count = 1
for testDir, testName in tests:
	startStr = "%d. Running %s"%(count, testDir)
	print("-"*len(startStr))
	print(startStr)
	print("-"*len(startStr))

	# Turn off output
	if not args.verbose:
		sys.stdout.flush()
		sys.stderr.flush()
		newStdOut = os.dup(1) # Copy StdOut
		newStdErr = os.dup(2) # Copy StdErr
		devNull = os.open('/dev/null', os.O_WRONLY)
		os.dup2(devNull, 1)
		os.dup2(devNull, 2)
	
	# Import test file
	testMod = imp.load_source("%s"%(testName), "%s/%s/%s"%(baseDir,testDir,testName))

	# Run the tests
	passed, total, timings, resDict = testMod.LMRegressionTester()

	# Reenable output
	if not args.verbose:
		os.close(devNull)
		sys.stdout = os.fdopen(newStdOut, 'w')
		sys.stderr = os.fdopen(newStdErr, 'w')

	# Update statistics
	totalTests += total
	passedTests += passed

	# Write results
	count2 = 0
	for k,v in resDict.items():
		print("{0:<16} - {1} - {2}s".format(k,v[0], timings[count2]))
		for err in v[1]:
			print(" "*20,err)
		count2 += 1
	print("")

	# Update test number
	count += 1



# Print overall Summary
print("")
print("=============")
print("|| Summary ||")
print("=============")
print("Total test:\t",totalTests)
print("Passed test:\t",passedTests)
print(round(100.0*float(passedTests)/float(totalTests),1),"% correct")
print("")

