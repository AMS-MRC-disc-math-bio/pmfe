# Subdirectories to recurse in build order
SUBDIR = bin pmfe rnascoring

# Install prefix
export PREFIX = /usr/local
export PREFIXSHARE = $(PREFIX)/share/pmfe/
export PREFIXBIN = $(PREFIX)/bin/

.PHONY: $(SUBDIR) recurse
.PHONY: install

$(MAKECMDGOALS) recurse: $(SUBDIR)

rnascoring: # no dependencies

pmfe: rnascoring

bin: pmfe rnascoring

install: bin
	-mkdir -p $(PREFIXSHARE)
	-cp -rv Turner99 $(PREFIXSHARE)

uninstall:
	-rm -vrf $(PREFIXSHARE)

clean:
	-rm -vf parametrizer # Delete any old binaries lying around

$(SUBDIR):
	$(MAKE) -C $@ $(MAKECMDGOALS)
