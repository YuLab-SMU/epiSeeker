html:
	quarto render --to pdf

pdf:
	quarto render --to pdf

clean:
	Rscript -e 'bookdown::clean_book()';\
	rm -rf _bookdown_files epiSeeker_cache epiSeeker_files

serve:
	quarto preview

push_doc:
	cd ../epiSeeker_gh_pages;\
	cp -r gh-pages/* ../epiSeeker_gh_pages