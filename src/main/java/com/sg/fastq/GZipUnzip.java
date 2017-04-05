package com.sg.fastq;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import java.util.zip.GZIPInputStream;

class GZipUnzip {

	public static void gunzipFastq(String filename) throws IOException {

		// Extra extension check - remove?
		if (!filename.endsWith(".gz")) { return; }

		System.out.println();
		System.out.println("*****************************************************");
		System.out.println("Unzipping: " + filename);

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
			gis.close();
			fos.close();
			System.out.println("*****************************************************");
			System.out.println();
	
		}
	}
}
