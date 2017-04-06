package com.sg.fastq;

import java.nio.file.Files;
import java.nio.file.Paths;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

class GZipUnzip {

	public static void gunzipFastq(String filename) throws IOException {

		// Extra extension check - remove?
		if (!filename.endsWith(".gz")) { return; }

		byte[] buffer = new byte[1024];

		try (
			FileInputStream fis = new FileInputStream(filename);
			GZIPInputStream gis = new GZIPInputStream(fis);
			FileOutputStream fos = new FileOutputStream(filename.replaceAll(".gz",""));			
		) {
			int len;
			while ((len = gis.read(buffer)) > 0) {
				fos.write(buffer, 0, len);
			}	
		}
	}

	public static void gzipFastq(String filename) throws IOException {

		// Extra extension check - remove?
		if (!filename.endsWith(".fastq")) { return; }

		byte[] buffer = new byte[1024];

		try (
			FileInputStream fis = new FileInputStream(filename);
			FileOutputStream fos = new FileOutputStream(filename.replaceAll(".fastq",".fastq.gz"));			
			GZIPOutputStream gos = new GZIPOutputStream(fos);
		) {
			int len;
			while ((len = fis.read(buffer)) > 0) {
				gos.write(buffer, 0, len);
			}
		}
	}
}
