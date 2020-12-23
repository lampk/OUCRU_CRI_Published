# 2020-12-22 --------------------------------------------------------------

library(workflowr)
wflow_build()
wflow_status()
wflow_publish()

# set up Git
##> # initiate the upstream tracking of the project on the GitHub repo
##> git remote add origin https://github.com/lampk/OUCRU_CRI_Published.git

##> # pull all files from the GitHub repo (typically just readme, license, gitignore)
##> git pull origin master

##> # set up GitHub repo to track changes on local machine
##> git push -u origin master

wflow_git_push()
wflow_publish(c("analysis/interim.Rmd", "analysis/about.Rmd", "analysis/index.Rmd", "analysis/labbook.Rmd"), "update")
