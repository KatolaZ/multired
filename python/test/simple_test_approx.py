import multired as mr
import sys


if len(sys.argv) < 2:
    print "Usage: %s <layer_list>" % sys.argv[0]
    sys.exit(1)

m = mr.multiplex_red(sys.argv[1], verbose=True, fit_degree=20)
part = m.compute_partitions_approx()

print "Partitions:..."
m.dump_partitions_approx()

m.draw_dendrogram_approx()

