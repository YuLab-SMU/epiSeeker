html:
	quarto render --to html

pdf:
	quarto render --to pdf

serve:
	quarto preview

push_doc:
	cp -r gh-pages/* ../epiSeeker_gh_pages