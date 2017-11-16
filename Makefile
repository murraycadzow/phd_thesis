FILES = index.Rmd 01-Introduction.Rmd 02-Methods.Rmd 03-SelectionPipeline.Rmd 04-Clustering.Rmd 05-Selection_and_association.Rmd 06-Conclusion_discussions.Rmd A1-Unimputed_coreExome_anlaysis.Rmd A5-Supplementary_Tables.Rmd references.bib

render_book : compile_book pdf_book

preview_pdf: pdf

compile_book : $(FILES)
	Rscript make_book.R	

pdf_book : compile_book
	mv _book/Thesis.pdf _book/Thesis_bookdown.pdf
	sed 's/\\includegraphics{/\\includegraphics{_bookdown_files\//g'< _book/Thesis.tex  > Thesis.tex
	#mv _bookdown_files/* .
	cp -r images	_bookdown_files
	pdflatex Thesis
	bibtex Thesis
	makeglossaries Thesis	
	pdflatex -bib Thesis
	pdflatex Thesis
	mv Thesis.pdf _book/Thesis.pdf
	#cd ../
	#rm -r Thesis_files
	xdg-open _book/Thesis.pdf



preview : $(FILES)	
	Rscript preview_thesis.R
	
pdf : preview
	mv _book/Thesis.pdf _book/Thesis_bookdown.pdf
	sed 's/\\bibliography/\\fontsize{10}{12}\n\\linespread\{1}\\selectfont\n\\\bibliography/' < _book/Thesis.tex > Thesis.tex
	ln -s _bookdown_files/Thesis_files .
	pdflatex Thesis
	bibtex Thesis
	makeglossaries Thesis	
	pdflatex -bib Thesis
	pdflatex Thesis
	rm Thesis_files
	mv Thesis.pdf _book/Thesis.pdf
	xdg-open _book/Thesis.pdf

.PHONY : clean
clean : 
	@rm -f  Thesis.Rmd Thesis.md Thesis.rds Thesis.tex Thesis.log Thesis.lof Thesis.bbl Thesis.toc Thesis.out Thesis.lot Thesis.ist Thesis.glsdefs Thesis.gls Thesis.glo Thesis.glg Thesis.blg Thesis.aux Thesis.alg Thesis.acr Thesis.acn tmp-*.tex
	@rm -r _book/ _bookdown_files/
 
