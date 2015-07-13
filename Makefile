#
# A very basic Makefile to compile the plugin
#
# For compatibility with old(er) ImageJ installations,
# ensure that 'java6' points to a 1.6-version java compiler
#
#

# Options for the java compiler
JFLAGS = -g -Xlint:unchecked -extdirs ./external -d ./ 
JC = javac6
JAR = jar

# remove command to clean up
RM = rm -rvf

.PHONY: clean doc

# Build the program
ESI_Analysis:	ESI_Analysis.class

ESI_Analysis.class: $(wildcard *.java)
	$(JC) $(JFLAGS) $(wildcard *.java)

# create jar file
jar	: ESI_Analysis
	$(JAR) -mcvf Manifest.txt ESI_$(shell date +%Y%m%d-%H%M).jar plugins.config \
	de/bio_photonics/*/*.class  

# create jar file (w. source)
jarsrc	: ESI_Analysis
	$(JAR) -mcvf Manifest.txt ESIsrc_$(shell date +%Y%m%d-%H%M).jar plugins.config \
	de/bio_photonics/*/*.class \
	*.java


# create javadoc
doc:
	javadoc -d doc/ -classpath ./ -extdirs ./external -subpackages de.bio_photonics.esi *.java 

# clean
clean :
	$(RM) de/bio_photonics/esi/*.class *.jar doc/*
