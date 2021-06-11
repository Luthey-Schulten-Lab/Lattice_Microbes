import time

# A test for passing test
#   0 -> bool, [error strings]
def test1():
	return True, []

# A test for failed test from return
#   0 -> bool, [error strings]
def test2():
	return False, ["Test failed","By construction"]

# A test for failed test from Exception
#   0 -> bool, [error strings]
def test3():
	raise Exception("Excepted by construction")



# This function is imported
#  0 -> int : number passed, 
#       int : total tests, 
#       []  : test times, 
#       {str : testname -> (str: Passed/Failed, [] : a list of error strings/exceptions)} : Test results 
def LMRegressionTester():
	# Return variables
	passedTests = 0
	testResults = {}
	timings = []

	# Tests array
	tests = []

	# Add all the tests in this file
	tests.append(("testPassed",test1))
	tests.append(("testFailed",test2))
	tests.append(("testException",test3))

	# Run tests
	for testName, test in tests:
		startT = time.time()
		try:
			res, reason = test()
			if res:
				testResults[testName] = ("Passed",[])
				passedTests += 1
			else:
				testResults[testName] = ("Failed",reason)
		except Exception as err:
			testResults[testName] = ("Failed",[err.message])
		endT = time.time()
		timings.append(endT-startT)

	return passedTests, len(tests), timings, testResults
	

