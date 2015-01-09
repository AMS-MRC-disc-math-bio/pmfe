# Subdirectories to recurse in build order
SUBDIR = rnascoring gtmfe bin

.PHONY: $(SUBDIR) recurse

$(MAKECMDGOALS) recurse: $(SUBDIR)

$(SUBDIR):
	$(MAKE) -C $@ $(MAKECMDGOALS)
