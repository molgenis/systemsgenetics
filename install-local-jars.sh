mvn install:install-file -Dfile=./lib/jsc.jar -DgroupId=jsc -DartifactId=jsc -Dversion=1.0 -Dpackaging=jar
mvn install:install-file -Dfile=./lib/molgenis-core-0.0.1.jar -DgroupId=nl.systemsgenetics -DartifactId=molgenis-core -Dversion=0.0.1-SNAPSHOT -Dpackaging=jar
mvn install:install-file -Dfile=./lib/sam-1.83.jar -DgroupId=net.sf -DartifactId=samtools -Dversion=1.83 -Dpackaging=jar
