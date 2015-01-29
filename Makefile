# Subdirectories to recurse in build order
SUBDIR = bin pmfe rnascoring

.PHONY: $(SUBDIR) recurse

$(MAKECMDGOALS) recurse: $(SUBDIR)

rnascoring: # no dependencies

pmfe: rnascoring

bin: pmfe rnascoring

$(SUBDIR):
	$(MAKE) -C $@ $(MAKECMDGOALS)
