preview : index.Rmd 01-Introduction.Rmd 02-Methods.Rmd 03-SelectionPipeline.Rmd 04-Clustering.Rmd 05-Selection_and_association.Rmd 06-Conclusion_discussions.Rmd references.bib
	Rscript preview_thesis.R
	cp _book/Thesis.tex .
	makeglossaries Thesis	
	Rscript preview_thesis.R

.PHONY : clean
clean: 
	rm Thesis.Rmd Thesis.tex _book/*
