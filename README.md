# dendro_hybrid_std

## Testing Bayesian climate modeling with hierarchical standardization

Daniel J Hocking and Laura G. Smith   

## Abstract


## Objectives


## Project Summary

**Project summary can be found at: [http://djhocking.github.io/dendro_bayes/](http://djhocking.github.io/dendro_bayes/)**

## Updating the Project Webpage

run `git pull origin master` to ensure the local project is synched with the GitHub version.

1. Edit the `index.Rmd` file on the master branch
2. knit the index file in RStudio. The following YAML code at the top of the Rmd file calls to the html template to incorporate css and javascript during the knitting:

```
output: 
  html_document: 
    keep_md: yes
    template: sheds-template.html
```

3. Edit any table details such as column and row names in the resulting `index.html` file.
4. git add, commit, and push the master branch files.
5. On the command line run `git checkout gh-pages`
5. Now in the gh-pages branch that generates the webpage, run `git checkout master --index.html`. This brings the file over from the master branch.
7. Add and commit the changes
8. `git push origin gh-pages`
9. `git checkout master` to get back on the master branch.

