molgenis-genotype-reader
========================

**Downloading latested jar**
http://www.molgenis.org/jenkins/job/systemsgenetics/nl.systemsgenetics$molgenis-genotype-reader/lastBuild/

Make sure to download the stand alone jar: molgenis-genotype-reader-******-stand-alone.jar	

**Examples**
The org.molgenis.genotype.examples package contains some basic examples. This should get you started


**How to build**
<ul>
	<li>Install Maven (http://maven.apache.org/) if you haven't done so already</li>
	<li>Follow instructions in the general systemsgenetics repro on how to instal local jars</li>
	<li>Build molgenis-genotype-reader (<code>mvn clean install)</code></li>
</ul>
<br/>

**How to open in Eclipse**
<ul>
	<li>
		Install the m2e Maven plugin; This is standard installed in the more recent versions of Eclipse, if
		you got an older version you probably need to install it.
	</li>
	<li>
		Import the project in Eclipse with <code>File/Import -> Existing Maven Projects</code>
		and select the molgenis-genotype-reader root folder
	</li>
</ul>

