/*
 * Stores and manipulates relevant exon information
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

public class Exon
{
	int chr;
	int posLeft;
	int posRight;
	String geneName;
	double FPKM;
	String referenceSequence;
	double SNPs;
	int numberOfOverlappingReads;
	ArrayList<MappedReads> reads;
	HashMap<Integer, Integer> SNPPositions;
	State state;

	public static final int READLENGTH = 50;

	/*
	 *  parse a line of data (modified SAM format) and input appropriate fields
	 */
	Exon(String dataLine)
	{
		String[] fields = dataLine.split("\t");
		try
		{
			this.chr  = Integer.parseInt(fields[0].substring(3));
		}catch(NumberFormatException e)
		{
			this.chr = 0; // for now, ignore X, Y, M chromosomes
		}
		this.posLeft = Integer.parseInt(fields[5]);
		this.posRight = Integer.parseInt(fields[6]);
		this.geneName = fields[7];
		this.FPKM = 0;
		this.SNPs = 0;
		this.reads = new ArrayList<MappedReads>();
		this.numberOfOverlappingReads = 0;
		this.SNPPositions = new HashMap<Integer, Integer>();
		this.state = new State(); 
	}

	// Pretty printing
	public String toString()
	{
		return("chr" + this.chr + "\t" + this.posLeft + "\t" + this.posRight + "\t" + this.geneName + "\t" + 
				this.SNPs + "\t" + this.FPKM + "\t" + this.state.stateName);

		//return(this.FPKM + "\t" + this.SNPs + "\t" + this.state.stateName);

		/*	return("chr" + this.chr + " (" + this.posLeft + ", " + this.posRight + ") " + 
				this.geneName + "\t" + this.FPKM + "\t" + this.SNPs + "\t" + 
				this.numberOfOverlappingReads + "\t" + this.SNPPositions.keySet().size() + "\t" + this.referenceSequence);
		 */
	}

	/*
	 *  Parses the file (br)
	 *  Reads file until the next chromosome
	 *  Loads into an array of Exon objects
	 */
	public static ArrayList<Exon> readExon(String exonFileName, int chromosome)
	{
		ArrayList<Exon> exons = new ArrayList<Exon>();
		String line = null; 
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(exonFileName));
			while( (line = br.readLine()) != null)
			{
				Exon exon = new Exon(line);
				exons.add(exon);
			}
		} catch (IOException e)
		{
			System.err.println("Error: Unable to process exon file");
			e.printStackTrace();
			System.exit(0);
		}

		return exons;
	}

	/*
	 *  Given an array of exons and expression data, calculate the FPKM levels per exon
	 *  Expressions & Exons should be sorted in order of beginning position
	 */
	static void getFPKM(ArrayList<Expression> Expressions, ArrayList<Exon> Exons)
	{
		int count = 0;
		int exprIndex = 0;

		// Helps to account for transcripts which are present in multiple exons
		int numberOfMovedExprPositions = 0;
		int numberOfMovedExonPositions = 0;

		for (int exonIndex = 0; exonIndex < Exons.size(); exonIndex++)
		{
			// until we are indexing the same chromosome, increment indices
			// get to the same chromosome in exons & expressions
			while(exprIndex < Expressions.size() 
					&& Exons.get(exonIndex).chr != Expressions.get(exprIndex).chr)
			{
				if(Exons.get(exonIndex).chr > Expressions.get(exprIndex).chr)
				{
					exprIndex++;
				}
				if(Exons.get(exonIndex).chr < Expressions.get(exprIndex).chr)
				{
					exonIndex++;	
					numberOfMovedExonPositions++;
				}
			}

			// go to the first transcript that in the exon
			while( exprIndex < Expressions.size() 
					&& Exons.get(exonIndex).chr == Expressions.get(exprIndex).chr 
					&& Exons.get(exonIndex).posLeft > Expressions.get(exprIndex).transcriptStart)
			{
				exprIndex++;
				numberOfMovedExprPositions++;
			}


			// for all the transcripts within the exon
			// add up the expression levels as long is its within the boundaries of the exon
			while ( exprIndex < Expressions.size() 
					&& Exons.get(exonIndex).chr  == Expressions.get(exprIndex).chr 
					&& Exons.get(exonIndex).posLeft <= Expressions.get(exprIndex).transcriptStart
					&& Exons.get(exonIndex).posRight >= Expressions.get(exprIndex).transcriptEnd )
			{
				numberOfMovedExprPositions++;
				count += Expressions.get(exprIndex).rpkm;
				exprIndex++;
			}

			Exons.get(exonIndex).FPKM = count;

			exonIndex = exonIndex - numberOfMovedExonPositions;
			exprIndex = exprIndex - numberOfMovedExprPositions;
			numberOfMovedExprPositions=0;
			numberOfMovedExonPositions=0;
			count = 0;

		}
	}

	/*
	 * Reads reference file and computes the reference sequence for the exon
	 */
	public void getReferenceSequence(RandomAccessFile inputReference)
	{
		String refRead = "";
		try
		{
			inputReference.seek(this.posLeft - 1);
			for(int i=0; i<(this.posRight - this.posLeft + 1); i++)
			{
				byte b = inputReference.readByte();
				char c = (char)b;
				if (c == '\n')
					i--;
				else if (c == '>') // ignore the header lines
				{
					i--;
					inputReference.readLine();
				}
				else
					refRead += ((char)b);
			}
		} catch (IOException e)
		{
			System.err.println("Error: Unable to process reference genome file: " + e.getMessage());
			e.printStackTrace();
			System.exit(0);
		}

		this.referenceSequence = refRead;
	}

	/*
	 * Gets the number of differences (ie SNPs) between the two sequences
	 */
	public static int getNumberOfDifferences(String s1, String s2)
	{
		if (s1.length() > READLENGTH)
			System.out.println("Read length > expected read length of " + READLENGTH); 
		// check if the lengths are equal
		if (s1.length() != s2.length())
			return -1;

		int count = 0;
		for(int pos = 0; pos < s1.length(); pos++)
		{
			if (s1.charAt(pos) != s2.charAt(pos))
				count++;
		}
		return count;
	}

	/*
	 * Gets the number of SNPs per exon
	 */
	public static void getSNPs(ArrayList<Exon> exons, ArrayList<MappedReads> mappedReads)
	{
		for (int exonIndex = 0; exonIndex < exons.size(); exonIndex++)
		{
			Exon exon = exons.get(exonIndex);
			if (mappedReads.isEmpty())
			{
				System.err.println("No mapped reads present");
				return;
			}

			int readsIndex = 0;
			MappedReads read = null;
			while(!mappedReads.isEmpty())
			{
				read = mappedReads.get(readsIndex);
				if (exon.posLeft > (read.startPosition + read.read.length()))
				{
					//System.out.println("Read removed: " + read.startPosition);
					mappedReads.remove(readsIndex);
				}
				else
					break;
			}

			while ( exon.overlaps(read) ) // mappedRead is within exon
			{
				calculateSNPPositions(exon, read);
				exon.numberOfOverlappingReads ++;
				readsIndex ++;
				if (readsIndex >= mappedReads.size())
					break;
				read = mappedReads.get(readsIndex);
				//System.out.println( "Read position: " + read.startPosition );
			}
			exon.SNPs = (double)exon.SNPPositions.keySet().size()/exon.length();
		}
	}

	int length()
	{
		return (this.posRight-this.posLeft+1);
	}

	private static void calculateSNPPositions(Exon exon, MappedReads read)
	{
		String referenceSeq = "";
		String readSeq = "";
		int exonStart = Math.max(exon.posLeft, read.startPosition) - exon.posLeft;
		int exonEnd = Math.min(exon.posRight, read.startPosition+read.read.length()) - exon.posLeft;
		int readStart = Math.max(exon.posLeft, read.startPosition) - read.startPosition;
		int readEnd = Math.min(exon.posRight, read.startPosition+read.read.length()) - read.startPosition;

		referenceSeq = exon.referenceSequence.toLowerCase().substring(exonStart, exonEnd);
		readSeq = read.read.toLowerCase().substring(readStart, readEnd);

		for(int pos = 0; pos < referenceSeq.length(); pos++)
		{
			if (referenceSeq.charAt(pos) != readSeq.charAt(pos))
				exon.SNPPositions.put(pos+Math.max(exon.posLeft, read.startPosition), 1);
		}

	}

	// Does exon (this) overlap with mapped reads?
	private boolean overlaps(MappedReads mappedRead)
	{
		if (this.chr == mappedRead.chr)
		{
			int overlapDistance = Math.min(this.posRight, mappedRead.startPosition+mappedRead.read.length())
			- Math.max(this.posLeft, mappedRead.startPosition);
			if (overlapDistance > 0)
				return true;
			else
				return false;

		}
		else
			return false;
	}

	// Sort exons in order of their genomic position
	static void sortExons(ArrayList<Exon> data)
	{
		final Comparator<Exon> order =
			new Comparator<Exon>() 
			{
			public int compare(Exon e1, Exon e2) {
				// chr diff & position diff
				int positionDiff = e1.chr - e2.chr;
				if (positionDiff == 0)
					positionDiff = e1.posLeft - e2.posLeft;
				if (positionDiff == 0)
					positionDiff = e1.posRight - e2.posRight;
				return positionDiff;
			}
			};
			Collections.sort(data, order);
	}

	public boolean containsPosition(int position)
	{
		if (this.posLeft <= position && position <= this.posRight)
			return true;
		else
			return false;
	}

}
