import multired as mr
import sys


if len(sys.argv) < 2:
    print "Usage: %s <layer_list>" % sys.argv[0]
    sys.exit(1)

print "Loading layers...",
m = mr.multiplex_red(sys.argv[1])
print "[DONE]"

print "Computing layer entropies...",
m.compute_layer_entropies()
print "[DONE]"

print "Computing JSD matrix...",
m.compute_JSD_matrix()
print "[DONE]"

print "Performing reduction...",
m.reduce()
print "[DONE]"

print "Getting partitons...",
part = m.compute_partitions()
print "[DONE]"

print "Partitions:...",
m.dump_partitions()
print "[DONE]"



