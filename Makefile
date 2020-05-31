PREFIX=/media/jeksterslab/scripts/r
PKG=$(PREFIX)/jeksterslabRsem

.PHONY: all clean rm

all : rm
	Rscript -e 'jeksterslabRpkg::pkg_build("$(PKG)", git = TRUE, github = TRUE, commit_msg = "Automated build")'

clean : rm
	-git add --all
	-git commit -m "[skip ci] Automated clean."
	-git push

rm :
	-rm -rf ${PKG}/docs/*
	-rm -rf ${PKG}/man/*
	-rm -rf ${PKG}/tests/testthat/*.html
	-rm -rf ${PKG}/tests/testthat/*.md
	-rm -rf ${PKG}/vignettes/*.html
	-rm -rf ${PKG}/vignettes/*.md
	-rm -rf ${PKG}/vignettes/notes/*.html
	-rm -rf ${PKG}/vignettes/notes/*_files
	-rm -rf ${PKG}/vignettes/notes/*_cache
	-rm -rf ${PKG}/vignettes/notes/*.md
	-rm -rf ${PKG}/vignettes/tests/*
