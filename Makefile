FILES = index.Rmd 01-Introduction.Rmd 02-Methods.Rmd 03-SelectionPipeline.Rmd 04-Clustering.Rmd 05-Selection_and_association.Rmd 06-Conclusion_discussions.Rmd references.bib

preview_pdf: pdf

pdf_book : $(FILES)
	

preview : $(FILES)	
	Rscript preview_thesis.R
	
pdf : preview
	cp _book/Thesis.tex .
	ln -s _bookdown_files/Thesis_files .
	pdflatex Thesis
	bibtex Thesis
	makeglossaries Thesis	
	pdflatex -bib Thesis
	pdflatex Thesis
	rm Thesis_files
	mv Thesis.pdf _book/

.PHONY : clean
clean : 
	@rm -f  Thesis.Rmd Thesis.tex Thesis.log Thesis.lof Thesis.bbl Thesis.toc Thesis.out Thesis.lot Thesis.ist Thesis.glsdefs Thesis.gls Thesis.glo Thesis.glg Thesis.blg Thesis.aux Thesis.alg Thesis.acr Thesis.acn 
	rm -r _book/ _bookdown_files/
	rm tmp-*.tex
