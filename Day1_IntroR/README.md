# Intro to R and Data Handling

Follow the instructions below to prepare for the workshop.

## Install R and RStudio 

### macOS

1. Install R
	1.	Under Download and Install R, click Download R for macOS. 
	2.	Choose the latest R-x.x.x.pkg file (top of the page). 
	3.	Open the .pkg file once it finishes downloading. 
	4.	Follow the installer prompts → Continue → Install. 
	5.	Enter your Mac password if asked. 

2. Install RStudio
	1.	Go to the RStudio download [page](https://posit.co/download/rstudio-desktop/)
	2.	Under All Installers, download the macOS version (.dmg).
	3.	Open the downloaded .dmg file.
	4.	Drag the RStudio icon into the Applications folder.
	5.	Open RStudio from Applications.

RStudio will automatically detect your R installation.

### Windows

1. Install R
	1.	Go to the CRAN R project [website](https://cran.r-project.org)
	2.	Click Download R for Windows.
	3.	Click base.
	4.	Download the R-x.x.x-win.exe installer.
	5.	Run the installer:
	•	Choose default options unless you need something specific.
	•	Click Next until Install.

2. Install RStudio
	1.	Go to the RStudio download [page](https://posit.co/download/rstudio-desktop/)
	2.	Under All Installers, download the Windows version (.exe).
	3.	Run the installer and click through the default installation steps.
	4.	Open RStudio from the Start Menu.

RStudio will automatically detect your R installation.

You can orient yourself to RStudio using this [guide](RStudio.pdf).

## Clone the repo

You may have been provided a Google Drive link to download the files from this repository. You should be sure you know where those files are, as we will use them in the instructions below. 

Another option is to clone this repository. Open Rstudio and in the bottom left corner you'll see the Files pane. Navigate to a directory where you would like to keep these files. In the bottom left pane, select the Terminal pane and execute the following command.

```
git clone https://github.com/BioDataAnalysisLab/BioData-Training-School-2025.git
```

Now you can navigate to Day1_IntroR within the new directory that was just created.

## Setup RStudio environment

1. Open the [Setup.R](Setup.R) script in RStudio. On mac, you can do this by File > Open File. 
2. In the top left pane, you will see a button that says Source in the top right corner. Press this and it will run the Setup Script.
3. In the bottom left console, read the output. If you did not pass all of the checks, something may be wrong.

## Let's code

The [IntroR.Rmd](Intro.Rmd) file has all of the code that we will go through together. However, I encourage you to open the [IntroR_working.Rmd](IntroR_working.Rmd) file, which has blanks for you to type the code yourself. THis will make for more effective learning. However, it may be helpful for you to have both open, so that you can copy and paste some code if you get lost. 

## Challenge yourself

At the end of the R Markdown file there is a Challenge. You can view the answers [here](ChallengeAnswers.Rmd)
