package com.sg.fastq;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import java.util.zip.GZIPInputStream;

public class GZipUnzip {

	public static void main(String[] args) {

		String filename = args[0];

		System.out.println("Input file name: " + filename);

		try {

			System.out.println();
			System.out.println("*****************************************************");
			System.out.println();

			System.out.println("File attributes...");

			File file = new File(filename);

			System.out.println("exists: " + file.exists());
			System.out.println("getAbsoluteFile: " + file.getAbsoluteFile());
			System.out.println("getAbsolutePath: " + file.getAbsolutePath());
			System.out.println("getCanonicalFile: " + file.getCanonicalFile());
			System.out.println("getCanonicalPath: " + file.getCanonicalPath());
			System.out.println("getName: " + file.getName());
			System.out.println("getParent: " + file.getParent());
			System.out.println("getParentFile: " + file.getParentFile());
			System.out.println("getPath: " + file.getPath());
			System.out.println("hashCode: " + file.hashCode());
			System.out.println("isAbsolute: " + file.isAbsolute());
			System.out.println("isDirectory: " + file.isDirectory());
			System.out.println("isFile: " + file.isFile());
			System.out.println("isHidden: " + file.isHidden());
			System.out.println("lastModified: " + file.lastModified());
			System.out.println("length: " + file.length());

			System.out.println();
			System.out.println("*****************************************************");
			System.out.println();
			
			System.out.println("FileInputStream attributes");

			FileInputStream fis = new FileInputStream(filename);

			System.out.println("Class name: " + fis.getClass().getName());
			System.out.println("Available before: " + fis.available());
			
			long start = System.currentTimeMillis();
			int count = 0;
			while (fis.read() != -1) {
				count++;
			}
			
			long timeToRead = System.currentTimeMillis() - start;
			timeToRead /= 1000;
			System.out.println("Seconds to read: " + timeToRead); 
			System.out.println("Byte count: " + count);
			System.out.println("Available after: " + fis.available());
			System.out.println("File descriptor (getFD): " + fis.getFD());

			System.out.println();
			System.out.println("*****************************************************");
			System.out.println();

		} catch (IOException ex) {

			ex.printStackTrace();

		}
	}

}
