package com.sg.fastq;

import java.io.IOException;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.DirectoryStream;

import com.sg.fastq.GZipUnzip;
import com.sg.fastq.CountPrimers;

/*
 * ProcessFastq
 *
 * This is the main class for the com.sg.fastq package.
 * The arguments passed to the main function are 
 * 
 * fastqDir : the directory holding all fastq files to be processed
 * primerDir : the directory holding all of the primer-list.csv and primer-route.csv
 * clusterDir : the directory we are saving all child fastq files produced by primer-based clustering
 *
 */
class ProcessFastq {

	private static final String CLUSTER_DIR = "clusters";
	private static final String PRIMER_COUNT_DIR = "primercounts";
	private static final String TARGET_COUNT_DIR = "targetcounts";

	public static void main (String[] args) throws IOException {

		String fastqDir = args[0];		// This dir is full of fastq files
		String primerDir = args[1];		// This dir contains csv primer lists and routes

		// Get the path separator
		String separator = System.getProperty("file.separator");
		Path fastqPath = Paths.get(fastqDir);
		Path primerPath = Paths.get(primerDir);
		Path fastqParent = fastqPath.getParent();		

		Path clusterPath = Paths.get(fastqParent + separator + CLUSTER_DIR);
		Path primerCountPath = Paths.get(fastqParent + separator + PRIMER_COUNT_DIR);
		Path targetCountPath = Paths.get(fastqParent + separator + TARGET_COUNT_DIR);

		// Create an empty directory for the fastq clusters, primer and target counts
		if (!Files.exists(clusterPath)) { Files.createDirectory(clusterPath); }
		if (!Files.exists(primerCountPath)) { Files.createDirectory(primerCountPath); }
		if (!Files.exists(targetCountPath)) { Files.createDirectory(targetCountPath); }

		// Unzip all the fastq.gz files
		try (
			DirectoryStream<Path> stream = Files.newDirectoryStream(fastqPath,"*.gz")
		) {
			for (Path p : stream) {
				GZipUnzip.gunzipFastq(p.toString());
			}
		}
		
		// Count primers in all .fastq files
		try (
			DirectoryStream<Path> stream = Files.newDirectoryStream(fastqPath, "*.fastq")
		) {
			Path firstOfPair = null;
			for (Path p : stream) {
				if (firstOfPair == null) {
					firstOfPair = p;
				} else if (isAPair(firstOfPair.toString(), p.toString())) {
					CountPrimers.countPrimers(primerPath, firstOfPair, p, primerCountPath, targetCountPath);
					DivideFastq.divideFastq(primerPath, primerPath, firstOfPair, p, clusterPath);
					firstOfPair = null;
				}
			}
		}

		// Delete sample fastq files
		try (
			DirectoryStream<Path> stream = Files.newDirectoryStream(fastqPath,"*.fastq")
		) {
			for (Path p : stream) {
				Files.delete(p);
			}
		}
	}
	private static boolean isAPair(String read1, String read2) {
		int index1;
		int index2;
		if ((index1 = read1.indexOf("_R")) >= 0 && (index2 = read2.indexOf("_R")) >= 0) {
			if (read1.substring(0, index1).equals(read2.substring(0, index2))) {
				return true;
			}	
		}
		return false;
	}
}
