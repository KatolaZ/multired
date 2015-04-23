import multired as mr
import sys


if len(sys.argv) < 2:
    print "Usage: %s <layer_list>" % sys.argv[0]
    sys.exit(1)

print "Loading layers...",
m = mr.multiplex_red(sys.argv[1], verbose=True)
print "[DONE]"

print "Getting partitons...",
part = m.compute_partitions()
print "[DONE]"

print "Partitions:..."
m.dump_partitions()

m.draw_dendrogram()


