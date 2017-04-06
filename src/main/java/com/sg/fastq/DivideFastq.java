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
import java.nio.file.StandardOpenOption;

import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;

class DivideFastq {	
	
	private final static Charset ENCODING = StandardCharsets.UTF_8;
	private final static int PRIMER_REGION = 30;
	private final static String NONE = "NONE";

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
					String line = null;
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
		String path1Name = path1.getFileName().toString();
		String path2Name = path2.getFileName().toString();
		String samplePrefix = path1Name.substring(0, path1Name.indexOf("_R"));
		Path samplePath = Paths.get(clusterPrefix + separator + samplePrefix);
		if (!Files.exists(samplePath)) { Files.createDirectory(samplePath); }
		String id1 = path1Name.substring(0, path1Name.indexOf(".fastq"));
		String id2 = path2Name.substring(0, path2Name.indexOf(".fastq"));
		Map<String, Path[]> pathMap = new HashMap<String, Path[]>();
		for (String target : routes.keySet()) {
			Path p1 = Paths.get(
								clusterPrefix 
								+ separator 
								+ samplePrefix 
								+ separator 
								+ id1.substring(0, id1.indexOf("_S"))
								+ "-" 
								+ target 
								+ id1.substring(id1.indexOf("_S"))
								+ ".fastq"
			);
			Path p2 = Paths.get(
								clusterPrefix 
								+ separator 
								+ samplePrefix 
								+ separator 
								+ id2.substring(0, id2.indexOf("_S"))
								+ "-" 
								+ target 
								+ id2.substring(id2.indexOf("_S"))
								+ ".fastq"
			);
			pathMap.put(target, new Path[]{p1,p2});
		}

		// Read the primer files data into arrays
		List<String> names = new ArrayList<String>();
		List<String> sequences = new ArrayList<String>();
		try (
			DirectoryStream<Path> stream = Files.newDirectoryStream(primerListPath, "*list.csv")
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
		int count = 0;
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
				
				// Assign headers
				String header1 = line1;
				String header2 = line2;			
	
				// Assign sequences	
				String seq1 = br1.readLine();
				String seq2 = br2.readLine();
				
				// Assign plus signs
				String plus1 = br1.readLine();
				String plus2 = br2.readLine();

				// Assign quality scores
				String quality1 = br1.readLine();
				String quality2 = br2.readLine();

				// Check for exact matches to primer sequences
				String p1 = checkForPrimer(names, sequences, seq1);
				String p2 = checkForPrimer(names, sequences, seq2);

				// Read target in primer name
				String target1 = !p1.equals(NONE) ? p1.split("-")[0] + "-" + p1.split("-")[1] : NONE;
				String target2 = !p2.equals(NONE) ? p2.split("-")[0] + "-" + p2.split("-")[1] : NONE;

				String[] quad1 = { header1, seq1, plus1, quality1 };
				String[] quad2 = { header2, seq2, plus2, quality2 };				
		
				// If targets are both NONE, then no more 
				if (target1.equals(NONE) && target2.equals(NONE)) {
					continue;
				}				

				// Calculate the destination(s)
				List<String> destinations = getDestinations(routes, target1, target2);

				// Write read pairs to each destination
				for (String dest : destinations) {
					StandardOpenOption oo = StandardOpenOption.APPEND;
					if (!Files.exists(pathMap.get(dest)[0])) {
						oo = StandardOpenOption.CREATE;
					}
					try (
						BufferedWriter bw1 = Files.newBufferedWriter(pathMap.get(dest)[0], ENCODING, oo);				
						BufferedWriter bw2 = Files.newBufferedWriter(pathMap.get(dest)[1], ENCODING, oo);				
					) {
						for (String s : quad1) {
							bw1.write(s);
							bw1.newLine();
						}
						for (String s : quad2) {
							bw2.write(s);
							bw2.newLine();
						}
					}
				}
				count++;
			}
		}
		System.out.println("-------Clustering Done: " + path1.getFileName().toString());
		System.out.println("-------Clustering Done: " + path2.getFileName().toString());
		System.out.println("-----Sequences Checked: " + count);
		System.out.println();
	
		// GZip the clusters
		try (
			DirectoryStream<Path> stream = Files.newDirectoryStream(samplePath, "*.fastq")
		) {
			for (Path p : stream) {
				GZipUnzip.gzipFastq(p.toString());
				Files.delete(p);
			}
		}		
	}

	// Use the routes map to generate a list of destinations for each read pair
	// Nested for loops parse each comma-separated destination spec, e.g.,
	// A-1/NONE/A-1/B-1 -- from this forward slash-delimited string
	// the destination is target A-1, and the allowed buddy reads
	// can have targets of NONE, A-1, or B-1 in order for the pair to be 
	// written to the A-1 cluster file
	private static List<String> getDestinations(Map<String, List<String>> routes, String target1, String target2) {
		List<String> out = new ArrayList<String>();
		String target = target1.equals(NONE) ? target2 : target1;
		String buddy = target.equals(target1) ? target2 : target1;
		if (routes.containsKey(target)) {
			List<String> destinations = routes.get(target);
			for (String line : destinations) { 
				String[] items = line.split(",");
				for (String item : items) {
					String[] fields = item.split("/");
					String d = fields[0];
					for (int i = 1; i < fields.length; i++) {
						if (buddy.equals(fields[i])) {
							out.add(d);
							break;
						}
					}
				}
			}
		}
		return out;
	}

	// Search for exact matches and increment the count array
	private static String checkForPrimer(List<String> names, List<String> sequences, String read) {
		for (int i = 0; i < sequences.size(); i++) {
			if (read.substring(0, PRIMER_REGION).contains(sequences.get(i))
					|| read.substring(0, PRIMER_REGION).contains(revComp(sequences.get(i)))) {
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
