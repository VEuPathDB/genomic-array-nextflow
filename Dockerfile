FROM eclipse-temurin:17-jdk-jammy

WORKDIR /app
COPY src ./src
RUN javac src/java/org/apidb/ggtools/array/ChIP_Chip_Peak_Finder.java -d .
RUN jar cf ChIPChipPeakFinder.jar org/apidb/ggtools/array/ChIP_Chip_Peak_Finder*.class
#CMD ["java", "YourClass"]
