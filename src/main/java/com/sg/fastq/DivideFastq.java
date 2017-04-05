package com.sg.fastq;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;

import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.DirectoryStream;

import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;

class DivideFastq {	
	
	private final static Charset ENCODING = StandardCharsets.UTF_8;
	private final static int PRIMER_REGION = 30;
	private final static String NONE = "none";

	public static void divideFastq (Path primerListPath, Path primerRoutesPath, Path path1, Path path2, Path clusterPath) throws IOException {

		// Read routes into a map
		// The map key is a target, e.g., A-1, and the list contains all destinations and allowable pairings
		Map<String, List<String>> routes = new HashMap<String, List<String>>();
		try (
			DirectoryStream<Path> stream = Files.newDirectoryStream(primerRoutesPath, "*routes.csv")
		) {
			for (Path p : stream) {
				try (
					BufferedReader br = Files.newBufferedReader(p, ENCODING);
				) {
					String[] line = null;
					while ((line = br.readLine()) != null) {
						String[] items = line.split(",");
						String target = items[0];
						routes.put(target, new ArrayList<String>());
						for (int i = 1; i < items.length; i++) {
							routes.get(target).add(items[i]);
						}
					}
				}
			}
		}

		// Create a path for each cluster file
		// Each cluster has a R1 and R2 file
		String separator = System.getProperty("file.separator");
		String clusterPrefix = clusterPath.toString();
		String id1 = path1.getFileName().toString().substring(0, path1.getFileName().toString().indexOf(".fastq"));
		String id2 = path2.getFileName().toString().substring(0, path2.getFileName().toString().indexOf(".fastq"));
		Map<String, Path[]> pathMap = new HashMap<String, Path[]>();
		for (String target : routes) {
			Path p1 = Paths.get(clusterPrefix + separator + id1 + target + ".fastq");
			Path p2 = Paths.get(clusterPrefix + separator + id2 + target + ".fastq");
			pathMap.put(target, new Path[]{p1,p2});
		}

		// Read the primer files data into arrays
		List<String> names = new ArrayList<String>();
		List<String> sequences = new ArrayList<String>();
		try (
			DirectoryStream<Path> stream = Files.newDirectoryStream(primerPath, "*list.csv")
		) {
			for (Path p : stream) {
				try (
					BufferedReader br = Files.newBufferedReader(p, ENCODING);
				) {
					String line = null;
					while ((line = br.readLine()) != null) {
						String[] items = line.split(",");
						names.add(items[0]);
						sequences.add(items[1]);
					}
				}
			}
		}		
			
		// Read files, looking at the corresponding R1 and R2 sequences together 
		System.out.println();
		System.out.println("------Start clustering: " + path1.getFileName().toString());
		System.out.println("------Start clustering: " + path2.getFileName().toString());
		try (
			BufferedReader br1 = Files.newBufferedReader(path1, ENCODING);
			BufferedReader br2 = Files.newBufferedReader(path2, ENCODING);
		) {
			String line1 = null;
			String line2 = null;
			while ((line1 = br1.readLine()) != null && (line2 = br2.readLine()) != null) {
				
				String header1 = line1;
				String header2 = line2;			
	
				String seq1 = br1.readLine();
				String seq2 = br2.readLine();
				
				// Check for exact matches to primer sequences
				String p1 = checkForPrimer(names, sequences, counts, seq1);
				String p2 = checkForPrimer(names, sequences, counts, seq2);

				String target1 = !p1.equals(NONE) ? p1.split("-")[0] + "-" + p1.split("-")[1] : NONE;
				String target2 = !p2.equals(NONE) ? p2.split("-")[0] + "-" + p2.split("-")[1] : NONE;

				// Assign plus signs
				String plus1 = br1.readLine();
				String plus2 = br2.readLine();

				// Assign quality scores
				String quality1 = br1.readLine();
				String quality2 = br2.readLine();

				count++;
			}
		}
		System.out.println("----------Clustering Done: " + path1.getFileName().toString());
		System.out.println("----------Clustering Done: " + path2.getFileName().toString());
		System.out.println("--------Sequences Checked: " + count);
		System.out.println();

		// Write the count file.txt
		String primerPrefix = primerCountPath.toString();
		String targetPrefix = targetCountPath.toString();
		String byPrimer = primerPrefix + separator + id + "-primers.csv";
		String byTarget = targetPrefix + separator + id + "-targets.csv";
		Path byPrimerPath = Paths.get(byPrimer);
		Path byTargetPath = Paths.get(byTarget);
		try (
			BufferedWriter bwPrimer = Files.newBufferedWriter(byPrimerPath, ENCODING);
			BufferedWriter bwTarget = Files.newBufferedWriter(byTargetPath, ENCODING);
		) {
			String newLine = System.getProperty("line.separator");

			// Put the target counts into a hashmap, keyed by target, "e.g., A-2"
			Map<String, Integer> targetCounts = new HashMap<String, Integer>();

			for (int i = 1; i < counts.length; i++) {
				String name = names.get(i);
				String line = name + "," + sequences.get(i) + "," + counts[i];
				line += newLine;
				// Write to bwPrimer
				bwPrimer.write(line, 0, line.length());

				// Populate targetCounts map
				// Target is the String preceding the second "-"
				String[] items = name.split("-");
				String target = items[0] + "-" + items[1];
				if (!targetCounts.containsKey(target)) {
					targetCounts.put(target, 0);
				}
				// Add individual primer count to total
				targetCounts.put(target, targetCounts.get(target)+counts[i]);
			}
			// Write to bwTarget
			// Iterate the entries in sorted order
			Map<String, Integer> sortedTargetCounts = new TreeMap<String, Integer>(targetCounts);
			for (Map.Entry<String, Integer> entry : targetCounts.entrySet()) {
				String line = entry.getKey() + "," + entry.getValue() + newLine;
				bwTarget.write(line, 0, line.length());
			}
		}
	}
	
	// Search for exact matches and increment the count array
	private static String checkForPrimer(List<String> names, List<String> sequences, int[] counts, String read) {
		for (int i = 0; i < sequences.size(); i++) {
			if (read.substring(0, PRIMER_REGION).contains(sequences.get(i))
					|| read.substring(0, PRIMER_REGION).contains(revComp(sequences.get(i)))) {
				counts[i]++;
				return names.get(i);
			}
		}
		return NONE;
	}
	
	// Move this to another class
	private static String revComp(String original) {
		String in = original.toUpperCase();
		String out = "";
		for (int i = original.length()-1; i >= 0; i--) {
			char base = original.charAt(i);
			switch (base) {
				case 'A': out += 'T';
									break; 
				case 'C': out += 'G';
									break; 
				case 'G': out += 'C';
									break; 
				case 'T': out += 'A';
									break; 
				default: out += base;
									break;
			}
		}
		return out;
	}
}
