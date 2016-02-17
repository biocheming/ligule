import sys
l = sys.argv[1].split(',')
i = int(sys.argv[2])
try:
    print l[i]
except IndexError:
    print "Indexing error"
finally:
    print l[-1]


print "doing some other stuff"
