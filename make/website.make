################################################################
# website (without the rest of the docs)

################################################################
#####  SECURITY -- check these values for lilypond.org #########
################################################################
ifeq ($(WEBSITE_ONLY_BUILD),1)
  ### for lilypond.org
  TOP_SRC_DIR=$(HOME)/src/lilypond
  TRUSTED_DIR=$(HOME)/lilypond/trusted-scripts
  top-src-dir=$(TOP_SRC_DIR)
  depth=.
  trusted-dir=$(TRUSTED_DIR)
  script-dir=$(trusted-dir)
  texi2html-init-file=$(trusted-dir)/lilypond-texi2html.init
  top-htaccess=$(trusted-dir)/lilypond.org.htaccess
  dir-htaccess=$(trusted-dir)/website-dir.htaccess
  TEXI2HTML_PROGRAM=$(HOME)/usr/bin/texi2html
  EXAMPLES=$(HOME)/media/ly-examples/
  PICTURES=$(HOME)/media/pictures
else
  ### for normal git
  script-dir=$(top-src-dir)/scripts/build/
  texi2html-init-file=$(top-src-dir)/Documentation/lilypond-texi2html.init
  top-htaccess=$(top-src-dir)/Documentation/web/server/lilypond.org.htaccess
  dir-htaccess=$(top-src-dir)/Documentation/web/server/website-dir.htaccess
  include $(config_make)
  # I assume this is run from top-build-dir
  EXAMPLES=Documentation/web/ly-examples/out-www/
  PICTURES=Documentation/pictures/out-www/
endif


################################################################
OUT=out-website

### only update this when the language compiles correctly!
# LANGUAGES = (site, de, es, fr, hu, it, ja, nl)
WEB_LANGS = es fr nl

TEXI2HTML=ONLY_WEB=1 TOP_SRC_DIR=$(top-src-dir) DEPTH=$(depth) PERL_UNICODE=SD $(TEXI2HTML_PROGRAM)

EXTRACT_TEXI_FILENAMES=python $(script-dir)/extract_texi_filenames.py
CREATE_VERSION=python $(script-dir)/create-version-itexi.py
CREATE_WEBLINKS=python $(script-dir)/create-weblinks-itexi.py
MASS_LINK=python $(script-dir)/mass-link.py
WEB_POST=python $(script-dir)/website_post.py

SERVER_FILES=$(top-src-dir)/Documentation/web/server/

# don't include web
MANUALS=$(wildcard $(top-src-dir)/Documentation/*.tely)
MANUALS+=$(top-src-dir)/Documentation/contributor.texi

website-test:
	echo $(TEXI2HTML)

website-version:
	mkdir -p $(OUT)
	$(CREATE_VERSION) $(top-src-dir) > $(OUT)/version.itexi
	$(CREATE_WEBLINKS) $(top-src-dir) > $(OUT)/weblinks.itexi

website-xrefs: website-version
	for l in '' $(WEB_LANGS); do \
		$(EXTRACT_TEXI_FILENAMES) \
			-I $(top-src-dir)/Documentation/ \
			-I $(top-src-dir)/Documentation/"$$l" \
			-I $(OUT) -o $(OUT) --split=node \
			$(top-src-dir)/Documentation/"$$l"/web.texi ;\
		for m in $(MANUALS); do \
			n=`echo "$$m" | sed 's/Documentation/Documentation\/'$$l'/'` ; \
			b=`basename "$$n" .texi`; \
			d=`basename "$$b" .tely`; \
			if [ -e "$$n" ] ; then \
				$(EXTRACT_TEXI_FILENAMES) \
				-I $(top-src-dir)/Documentation/ \
				-I $(top-src-dir)/Documentation/"$$l" \
				-I $(top-src-dir)/Documentation/"$$l"/"$$d"/ \
				-I $(OUT) -o $(OUT) "$$n" ; \
			fi ; \
		done; \
	done;



website-texinfo: website-version website-xrefs
	for l in '' $(WEB_LANGS); do \
	        if test -n "$$l"; then \
			langopt=--lang="$$l"; \
			langsuf=.$$l; \
		fi; \
		$(TEXI2HTML) --prefix=index \
			--split=section \
			--I=$(top-src-dir)/Documentation/"$$l" \
			--I=$(top-src-dir)/Documentation/ \
			--I=$(OUT) \
			$$langopt \
			--init-file=$(texi2html-init-file) \
			-D web_version \
			--output=$(OUT)/"$$l" \
			$(top-src-dir)/Documentation/"$$l"/web.texi ; \
		find $(OUT)/$$l/ -name '*.html' | xargs grep -L 'UNTRANSLATED NODE: IGNORE ME' | sed 's!$(OUT)/'$$l'/!!g' | xargs $(MASS_LINK) --prepend-suffix="$$langsuf" hard $(OUT)/$$l/ $(OUT)/website/ ; \
	done


website-css:
	cp $(top-src-dir)/Documentation/css/*.css $(OUT)/website/

website-pictures:
	mkdir -p $(OUT)/website/pictures/
	cp $(PICTURES)/* $(OUT)/website/pictures/
	ln -sf website/pictures $(OUT)/pictures

website-examples:
	mkdir -p $(OUT)/website/ly-examples
	cp $(EXAMPLES)/* $(OUT)/website/ly-examples

web-post:
	$(WEB_POST) $(OUT)/website/

website: website-texinfo website-css website-pictures website-examples web-post
	cp $(SERVER_FILES)/favicon.ico $(OUT)/website/
	cp $(SERVER_FILES)/robots.txt $(OUT)/website/
	cp $(top-htaccess) $(OUT)/.htaccess
	cp $(dir-htaccess) $(OUT)/website/.htaccess


