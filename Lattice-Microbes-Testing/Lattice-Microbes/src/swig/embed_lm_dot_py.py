import sys

with open(sys.argv[1]) as inf, open(sys.argv[2],"w") as outf:
    print('static const char lm_swig_wrapper[] = ', file=outf)
    for line in inf:
        l = line.rstrip().replace('\\', '\\\\').replace('"', '\\"')
        print('"{}\\n"'.format(l), file=outf)
    print(';', file=outf)

