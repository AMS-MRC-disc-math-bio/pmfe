# Subdirectories to recurse in build order
SUBDIR = rnascoring pmfe bin

.PHONY: $(SUBDIR) recurse

$(MAKECMDGOALS) recurse: $(SUBDIR)

rnascoring: # no dependencies

pmfe: rnascoring

bin: pmfe rnascoring

$(SUBDIR):
	$(MAKE) -C $@ $(MAKECMDGOALS)
