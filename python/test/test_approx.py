import multired as mr
import sys


if len(sys.argv) < 2:
    print "Usage: %s <layer_list>" % sys.argv[0]
    sys.exit(1)

print "Loading layers...",
m = mr.multiplex_red(sys.argv[1])
print "[DONE]"

print "Computing layer entropies (approx)...",
m.compute_layer_entropies_approx()
print "[DONE]"

print "Computing JSD matrix (approx)...",
m.compute_JSD_matrix_approx()
print "[DONE]"

print "Performing reduction (approx)...",
m.reduce_approx()
print "[DONE]"


print "Getting partitons...",
part = m.compute_partitions_approx()
print "[DONE]"

print "Partitions:..."
m.dump_partitions_approx()
print "[DONE]"


