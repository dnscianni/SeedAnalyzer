/**
 * Research provided by Science Educational Enhancement Services, and the Hearst Foundation.  
 * Work is advised by Dr. Amar Raheja, of California State Polytechnic University, Pomona.  
 * Overall supervision is done by Dr. Steve Alas, the director of SEES at Cal Poly.  
 * This research project is written by David Scianni, a student at Cal Poly, Pomona.
 */
package edu.csupomona.cs.sees.research;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;

import org.opencv.*;
import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfByte;
import org.opencv.highgui.Highgui;
import org.opencv.imgproc.Imgproc;

import javax.imageio.ImageIO;

/**
 * Seed Analyzer will take in an image of seeds, and modify it in order to
 * analyze and deduce the number of seeds, and the length of each seed.
 * 
 * @author David Scianni
 * 
 */
public class SeedAnalyzer {

	// final private static int FILTER_SIZE = 9;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.loadLibrary(Core.NATIVE_LIBRARY_NAME);

		BufferedImage img = null;

		try {
			img = ImageIO.read(new File("seed1.JPG"));
		} catch (IOException e) {
			e.printStackTrace();
		}
		File outputfile = new File("lines	.JPG");
		try {
			//ImageIO.write(convertImg(fillIn(ZhangSuenThinning(noiseFilter(noiseFilter(noiseFilter(noiseFilter(noiseFilter(cutImage(threshold2(isolateChannel(img, 2),isolateChannel(img, 1), 48)),5),5),7),7),7))), img.getType()), "JPG", outputfile);
			printInfo(calcSeedLength(fillIn(ZhangSuenThinning(noiseFilter(noiseFilter(noiseFilter(noiseFilter(noiseFilter(cutImage(threshold2(isolateChannel(img, 2),isolateChannel(img, 1), 48)),5),5),7),7),7)))), calcPixelToMM(ZhangSuenThinning(noiseFilter(threshold(isolateChannel(img, 3), 20),5))));
			ImageIO.write(convertImg(ZhangSuenThinning(noiseFilter(threshold(isolateChannel(img, 3), 20),5)), img.getType()), "JPG", outputfile);
			//ImageIO.write(convertImg(fillIn(ZhangSuenThinning(noiseFilter(noiseFilter(noiseFilter(noiseFilter(noiseFilter(cutImage(threshold2(isolateChannel(img, 2),isolateChannel(img, 1), 48)),5),5),7),7),7))),img.getType()), "JPG", outputfile);
			/*
			 * ImageIO.write( cutImage(threshold2(isolateChannel(img, 2),
			 * isolateChannel(img, 1), 45)), "JPG", outputfile);
			 */
			//ImageIO.write(convertImg((noiseFilter(noiseFilter(noiseFilter(noiseFilter(cutImage(threshold2(isolateChannel(img, 2),isolateChannel(img, 1), 47)),5), 5), 7), 7)), img.getType()), "JPG", outputfile);
			//ImageIO.write(convertImg((noiseFilter(noiseFilter(noiseFilter(cutImage(threshold2(isolateChannel(img, 2),isolateChannel(img, 1), 50)),5), 5), 5)), img.getType()), "JPG", outputfile);
			/*
			 * noiseFilter(cutImage(threshold(isolateChannel(img, 2), 60))),
			 * "JPG", outputfile);
			 */
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static BufferedImage convertImg(int[][] img, int type) {
		BufferedImage resultImg = new BufferedImage(img.length,
				img[0].length, type);
		
		for(int i = 0; i < resultImg.getWidth(); i++) {
			for(int j = 0; j < resultImg.getHeight(); j++) {
				if(img[i][j] == 0) {
					resultImg.setRGB(i, j, Color.BLACK.getRGB());
				} else {
					resultImg.setRGB(i, j, Color.WHITE.getRGB());
				}
			}
		}
		return resultImg;
	}

	/**
	 * isolateChannel will take in a buffered image, and will go through each
	 * pixel of the image, calculating the color values of each pixel for its
	 * red, green, and blue channels. Then based on the option given, it will
	 * change the color value of that pixel to only the red, green, or blue
	 * channel. The resulting image is returned as the same image, but with only
	 * one color channel, instead of the three original channels.
	 * 
	 * @param img
	 *            The image to be manipulated.
	 * @param rgb
	 *            1 if isolating red, 2 if isolating green, 3 if isolating blue;
	 *            defaults to blue.
	 * @return a new image with isolated color channel.
	 */
	public static BufferedImage isolateChannel(BufferedImage img, int rgb) {
		BufferedImage resultImg = new BufferedImage(img.getWidth(),
				img.getHeight(), img.getType());
		int red, green, blue, pixel = 0;

		for (int i = 0; i < img.getWidth(); i++) {
			for (int j = 0; j < img.getHeight(); j++) {
				int color = img.getRGB(i, j);
				pixel = 0;

				red = (color & 0x00ff0000) >> 16;
				green = (color & 0x0000ff00) >> 8;
				blue = (color & 0x000000ff) >> 0;

				if (rgb == 1) {
					pixel = pixel | red << 16;
				} else if (rgb == 2) {
					pixel = pixel | green << 8;
				} else {
					pixel = pixel | blue;
				}

				resultImg.setRGB(i, j, pixel);
			}
		}
		return resultImg;
	}

	/**
	 * threshold will take in an image, going through each pixel of the image
	 * and comparing its rgb value to a number given as a parameter. If the
	 * value is greater than or equal to the number given, then the pixel's
	 * color value will be changed to black. Otherwise, it will be changed to
	 * white. This will turn the image into a binary image, focusing on part of
	 * the image but turning the important sections white, while turning the
	 * unimportant areas black.
	 * 
	 * @param img
	 *            The image to be manipulated.
	 * @param threshNum
	 *            The number that is to be compared with the rgb value
	 * @return a new binary image that consists of black and white pixels.
	 */
	public static int[][] threshold(BufferedImage img, int threshNum) {
		int[][] resultImg = new int[img.getWidth()][img.getHeight()];
		int red, green, blue = 0;

		for (int i = 0; i < img.getWidth(); i++) {
			for (int j = 0; j < img.getHeight(); j++) {
				int color = img.getRGB(i, j);

				red = (color & 0x00ff0000) >> 16;
				green = (color & 0x0000ff00) >> 8;
				blue = (color & 0x000000ff) >> 0;

				if ((red + green + blue) / 3 < threshNum) {
					resultImg[i][j] = 1;
				} else {
					resultImg[i][j] = 0;
				}
			}
		}

		return resultImg;
	}

	public static int[][] threshold2(BufferedImage img,
			BufferedImage img2, int threshNum) {
		int[][] resultImg = new int[img.getWidth()][img.getHeight()];
		int red, green, blue, red2, green2, blue2 = 0;

		for (int i = 0; i < img.getWidth(); i++) {
			for (int j = 0; j < img.getHeight(); j++) {
				int color = img.getRGB(i, j);
				int color2 = img2.getRGB(i, j);

				red = (color & 0x00ff0000) >> 16;
				green = (color & 0x0000ff00) >> 8;
				blue = (color & 0x000000ff) >> 0;

				red2 = (color2 & 0x00ff0000) >> 16;
				green2 = (color2 & 0x0000ff00) >> 8;
				blue2 = (color2 & 0x000000ff) >> 0;

				if ((red + green + blue) / 3 >= threshNum
						&& (red2 + green2 + blue2) / 3 >= threshNum) {
					resultImg[i][j] = 1;
				} else {
					resultImg[i][j] = 0;
				}
			}
		}

		return resultImg;
	}

	public static int[][] cutImage(int[][] img) {
		int[][] resultImg = new int[img.length - 150][img[0].length-500];
		for (int i = 0, a = 99; i < resultImg.length; i++, a++) {
			for (int j = 0, b = 99; j < resultImg[0].length; j++, b++) {
				resultImg[i][j] = img[a][b];
			}
		}
		return resultImg;
	}

	public static int median(int[] a) {
		mergeSort(a);
		if (a.length % 2 == 1)
			return a[a.length / 2];
		else
			return ((a[a.length / 2] + a[a.length / 2 - 1]) / 2);
	}

	public static int[] getArray(int[][] img, int x, int y, int z) {
		int[] a;
		int xmin, xmax, ymin, ymax;

		xmin = x - z / 2;
		xmax = x + z / 2;
		ymin = y - z / 2;
		ymax = y + z / 2;

		if (xmin < 0)
			xmin = 0;
		if (xmax > (img.length - 1))
			xmax = img.length - 1;
		if (ymin < 0)
			ymin = 0;
		if (ymax > (img[0].length - 1))
			ymax = img[0].length - 1;
		a = new int[(xmax - xmin + 1) * (ymax - ymin + 1)];
		for (int i = xmin, k = 0; i <= xmax; i++) {
			for (int j = ymin; j <= ymax; j++, k++) {
				a[k] = img[i][j];
			}
		}
		return a;
	}

	public static int[][] noiseFilter(int[][] img, int filterSize) {
		int[][] resultImg = new int[img.length][img[0].length];

		int[] a;
		int R, G, B;

		for (int k = 0; k < img[0].length; k++) {
			for (int j = 0; j < img.length; j++) {
				a = getArray(img, j, k, filterSize);
				int[] red, green, blue;
				red = new int[a.length];
				green = new int[a.length];
				blue = new int[a.length];
				for (int i = 0; i < a.length; i++) {
					red[i] = (a[i] & 0x00ff0000) >> 16;
					green[i] = (a[i] & 0x0000ff00) >> 8;
					blue[i] = (a[i] & 0x000000ff) >> 0;
				}
				R = median(red);
				G = median(green);
				B = median(blue);
				resultImg[j][k] = (R << 16 | G << 8 | B);
			}
		}
		return resultImg;
	}

	public static void mergeSort(int[] a) {
		if (a.length > 1) {
			int[] left = leftHalf(a);
			int[] right = rightHalf(a);

			mergeSort(left);
			mergeSort(right);

			merge(a, left, right);
		}
	}

	public static int[] leftHalf(int[] a) {
		int size1 = a.length / 2;
		int[] left = new int[size1];
		for (int i = 0; i < size1; i++) {
			left[i] = a[i];
		}
		return left;
	}

	public static int[] rightHalf(int[] a) {
		int size1 = a.length / 2;
		int size2 = a.length - size1;
		int[] right = new int[size2];
		for (int i = 0; i < size2; i++) {
			right[i] = a[i + size1];
		}
		return right;
	}

	public static void merge(int[] result, int[] left, int[] right) {
		int i1 = 0;
		int i2 = 0;

		for (int i = 0; i < result.length; i++) {
			if (i2 >= right.length
					|| (i1 < left.length && left[i1] <= right[i2])) {
				result[i] = left[i1];
				i1++;
			} else {
				result[i] = right[i2];
				i2++;
			}
		}
	}
	
	public static int[][] ZhangSuenThinning(int[][] img) {
        int a, b;
 
        LinkedList<Point> pointsToChange = new LinkedList<Point>();
        boolean hasChange;
 
        do {
 
            hasChange = false;
            for (int y = 1; y + 1 < img.length; y++) {
                for (int x = 1; x + 1 < img[y].length; x++) {
                    a = getA(img, y, x);
                    b = getB(img, y, x);
                    if ( img[y][x]==1 && 2 <= b && b <= 6 && a == 1
                            && (img[y - 1][x] * img[y][x + 1] * img[y + 1][x] == 0)
                            && (img[y][x + 1] * img[y + 1][x] * img[y][x - 1] == 0)) {
                        pointsToChange.add(new Point(x, y));
                        //img[y][x] = 0;
                        hasChange = true;
                    }
                }
            }
 
            for (Point point : pointsToChange) {
                img[point.getY()][point.getX()] = 0;
            }
 
            pointsToChange.clear();
 
            for (int y = 1; y + 1 < img.length; y++) {
                for (int x = 1; x + 1 < img[y].length; x++) {
                    a = getA(img, y, x);
                    b = getB(img, y, x);
                    if ( img[y][x]==1 && 2 <= b && b <= 6 && a == 1
                            && (img[y - 1][x] * img[y][x + 1] * img[y][x - 1] == 0)
                            && (img[y - 1][x] * img[y + 1][x] * img[y][x - 1] == 0)) {
                        pointsToChange.add(new Point(x, y));
 
                        hasChange = true;
                    }
                }
            }
 
            for (Point point : pointsToChange) {
                img[point.getY()][point.getX()] = 0;
            }
 
            pointsToChange.clear();
 
        } while (hasChange);
 
        return img;
    }
 
    private static int getA(int[][] img, int y, int x) {
 
        int count = 0;
        //p2 p3
        if (img[y - 1][x] == 0 && img[y - 1][x + 1] == 1) {
            count++;
        }
        //p3 p4
        if (img[y - 1][x + 1] == 0 && img[y][x + 1] == 1) {
            count++;
        }
        //p4 p5
        if (img[y][x + 1] == 0 && img[y + 1][x + 1] == 1) {
            count++;
        }
        //p5 p6
        if (img[y + 1][x + 1] == 0 && img[y + 1][x] == 1) {
            count++;
        }
        //p6 p7
        if (img[y + 1][x] == 0 && img[y + 1][x - 1] == 1) {
            count++;
        }
        //p7 p8
        if (img[y + 1][x - 1] == 0 && img[y][x - 1] == 1) {
            count++;
        }
        //p8 p9
        if (img[y][x - 1] == 0 && img[y - 1][x - 1] == 1) {
            count++;
        }
        //p9 p2
        if (img[y - 1][x - 1] == 0 && img[y - 1][x] == 1) {
            count++;
        }
 
        return count;
    }
 
    private static int getB(int[][] img, int y, int x) {
 
        return img[y - 1][x] + img[y - 1][x + 1] + img[y][x + 1]
                + img[y + 1][x + 1] + img[y + 1][x] + img[y + 1][x - 1]
                + img[y][x - 1] + img[y - 1][x - 1];
    }
	
    public static int[][] HilditchsThinning(int[][] img) {
        int a, b;
        boolean hasChange;
        do {
            hasChange = false;
            for (int y = 1; y + 1 < img.length; y++) {
                for (int x = 1; x + 1 < img[y].length; x++) {
                    a = getA2(img, y, x);
                    b = getB2(img, y, x);
                    if (img[y][x]==1 && 2 <= b && b <= 6 && a == 1
                            && ((img[y - 1][x] * img[y][x + 1] * img[y][x - 1] == 0) || (getA2(img, y - 1, x) != 1))
                            && ((img[y - 1][x] * img[y][x + 1] * img[y + 1][x] == 0) || (getA2(img, y, x + 1) != 1)))
                    {
                        img[y][x] = 0;
                        hasChange = true;
                    }
                }
            }
        } while (hasChange);
 
        return img;
    }
 
    private static int getA2(int [][] img, int y, int x){
 
        int count =0;
        //p2 p3
        if(y-1 >= 0 && x+1 < img[y].length && img[y-1][x]==0 && img[y-1][x+1]==1 ){
        count++;
        }
        //p3 p4
        if(y-1>= 0 && x+1 < img[y].length && img[y-1][x+1]==0 && img[y][x+1]==1){
        count++;
        }
        //p4 p5
        if(y+1<img.length && x+1 < img[y].length && img[y][x+1]== 0 && img[y+1][x+1]==1){
        count++;
        }
        //p5 p6
        if(y+1<img.length && x+1 < img[y].length && img[y+1][x+1]==0 && img[y+1][x]==1){
        count++;
        }
        //p6 p7
        if(y+1<img.length && x-1 >= 0 && img[y+1][x]==0 && img[y+1][x-1]==1){
        count++;
        }
        //p7 p8
        if(y+1<img.length && x-1 >= 0 && img[y+1][x-1]== 0 && img[y][x-1]==1){
        count++;
        }
        //p8 p9
        if(y-1>=0 && x-1 >= 0 && img[y][x-1]==0 && img[y-1][x-1]==1){
        count++;
        }
        //p9 p2
        if(y-1>=0 && x-1>= 0 &&img[y-1][x-1]==0 && img[y-1][x]==1){
        count++;
        }
 
    return count;
    }
 
    private static int getB2(int[][] img, int y, int x) {
 
        return img[y - 1][x] + img[y - 1][x + 1] + img[y][x + 1]
                + img[y + 1][x + 1] + img[y + 1][x] + img[y + 1][x - 1]
                + img[y][x - 1] + img[y - 1][x - 1];
    }
    
//    public static int[][] close(int[][] im) {
//		return erode( dilate( im));
//	}
    
    public static int[][] dilate(int[][] im) {
    	
		return percentileFilter(im, 100.0);
	}
    
    public static int[][] erode(int[][] im) {
		
		return percentileFilter(im, 0.0);
	}
    
    public static int[][] percentileFilter(int[][] im, double perc) {
		return percentileFilter(im, perc, true);
	}
    
    public static int[][] percentileFilter(int[][] im, double perc, boolean symmetric) {

		int[][] out = im.clone();

		int width = out.length;
		int height = out[0].length;

		// Send image to structuring element for speed and set symmetry.
		for(int i=0; i<width; i++) {
			for (int j=0; j<height; j++) {
				out[i][j] = (int)Math.round(
					percentile( 
					out[i], perc ) );
//				System.out.print(x);System.out.print(" ");System.out.print(y);System.out.println();
			}
		}
		return out;
	}
    
    private static double percentile(int[] values, double perc) {
		Arrays.sort(values);

		// n is the rank of the percentile value in the sorted array.
		double n = (perc/100.0)*(values.length - 1) + 1.0;

		if (Math.round(n)==1) {
			return values[0];
		} else if (Math.round(n)==values.length) {
			return values[values.length-1];
		} else {
			int k = (int)n;
			double d = n-k;
			return values[k-1] + d*(values[k] - values[k-1]);
		}
	}
    
    public static int[][] close(int[][] img, int type) {
    	Mat initialImg = matify(convertImg(img,type));
    	Mat kernel = new Mat();
    	Mat resultImg = new Mat();
    	Imgproc.morphologyEx(initialImg, resultImg, Imgproc.MORPH_CLOSE, kernel);
    	
    	int[][] result = BuffToArray(getImage(resultImg));
    	return result;
    }
    
    public static int[][] dilate(int[][] img, int type) {
    	Mat initialImg = matify(convertImg(img,type));
    	Mat kernel = new Mat();
    	Mat resultImg = new Mat();
    	Imgproc.morphologyEx(initialImg, resultImg, Imgproc.MORPH_DILATE, kernel);
    	
    	int[][] result = BuffToArray(getImage(resultImg));
    	return result;
    }
    
    public static int[][] BuffToArray(BufferedImage image) {
    	int[][] resultImg = new int[image.getWidth()][image.getHeight()];
    	for(int i = 0; i < resultImg.length; i++) {
			for(int j = 0; j < resultImg[0].length; j++) {
				if(image.getRGB(i, j) == Color.BLACK.getRGB()) {
					resultImg[i][j] = 0;
				} else {
					resultImg[i][j] = 1;
				}
			}
		}
    	return resultImg;
	}

	public static Mat matify(BufferedImage im) {
        // Convert INT to BYTE
        //im = new BufferedImage(im.getWidth(), im.getHeight(),BufferedImage.TYPE_3BYTE_BGR);
        // Convert bufferedimage to byte array
        byte[] pixels = ((DataBufferByte) im.getRaster().getDataBuffer())
                .getData();

        // Create a Matrix the same size of image
        Mat image = new Mat(im.getHeight(), im.getWidth(), CvType.CV_8UC3);
        // Fill Matrix with image values
        image.put(0, 0, pixels);

        return image;

    }
    public static BufferedImage getImage(Mat img){
    	
    	MatOfByte mob = new MatOfByte();
		//convert the matrix into a matrix of bytes appropriate for
		//this file extension
		Highgui.imencode(".JPG", img ,mob); 
		//convert the "matrix of bytes" into a byte array
		 byte[] byteArray = mob.toArray();
		 BufferedImage bufImage = null;
		 try {
		        InputStream in = new ByteArrayInputStream(byteArray);
		        bufImage = ImageIO.read(in);
		    } catch (Exception e) {
		        e.printStackTrace();
		    }
		 return bufImage;
	}
    
    public static int[][] fillIn(int[][] img) {
		for(int i = img[0].length-1; i > 400; i--) {
			for(int j = 0; j < img.length; j++){
				if(img[j][i] == 1) {
					if(img[j][i+1] != 1 && img[j-1][i+1] != 1 && img[j+1][i+1] != 1) {
					fillInHelper(img, i, j);
					}
				}
			}
		}
		return img;
	}
	
	public static void fillInHelper(int[][] img, int i, int j) {
		while(img[j][i] == 1) {
			if(img[j][i-1] == 1) {
				i = i-1;
			} else if(img[j-1][i-1] == 1) {
				j = j-1;
				i = i-1;
			} else if(img[j+1][i-1] == 1) {
				j = j+1;
				i = i-1;
			} else {
				break;
			}
		}
		int bottomJ = j;
		int bottomI = i;
		int topJ = -1;
		int topI = -1;
		for(int k = i-1; k >= 0; k--) {
			if(img[j][k] == 1) {
				topJ = j;
				topI = k;
				break;
			} else if(img[j-1][k] == 1) {
				topJ = j-1;
				topI = k;
				break;
			} else if(img[j+1][k] == 1) {
				topJ = j+1;
				topI = k;
				break;
			} else if(img[j-2][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+2][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-3][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+3][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-4][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+4][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-5][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+5][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-6][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+6][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-7][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+7][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-8][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+8][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-9][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+9][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-10][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+10][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-11][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+11][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-12][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+12][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-13][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+13][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-14][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+14][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			} else if(img[j-15][k] == 1) {
				topJ = j-2;
				topI = k;
				break;
			} else if(img[j+15][k] == 1) {
				topJ = j+2;
				topI = k;
				break;
			}
		}
		if(topJ >= 0) {
			int height = Math.abs(bottomI - topI);
			int test;
			if (bottomJ-topJ != 0) {
				test = height / Math.abs(bottomJ - topJ);
			} else {
				test = height;
			}
			for(int k = bottomI, l = bottomJ, m = 1; k > topI; k--, m++) {
				if(m == test) {
					if(topJ > bottomJ) {
						l++;
					} else if (topJ < bottomJ) {
						l--;
					}
					m = 1;
				}
				img[l][k] = 1;
			}
		}
	}
	
	public static float calcPixelToMM(int[][] img) {
		int result = 1;
		int i = img[0].length - 1;
		int j = img.length/2;
		while(true) {
			if(img[j][i] == 1) {
				i--;
				break;
			} else {
				i--;
			}
		}
		
		while(img[j][i] == 0) {
			result++;
			i--;
		}
		return result;
	}
	
	public static ArrayList<Seed> calcSeedLength(int[][] img) {
		ArrayList<Seed> list = new ArrayList<Seed>();
		
		for(int i = img[0].length-1; i > 450; i--) {
			for(int j = 0; j < img.length; j++){
				if(img[j][i] == 1) {
					if(img[j][i+1] != 1 && img[j-1][i+1] != 1 && img[j+1][i+1] != 1 && img[j-1][i] != 1) {
						list.add(new Seed(j, calcSeedHelper(img, i , j)));
					}
				}
			}
		}
		Collections.sort(list);
		return list;
	}

	private static float calcSeedHelper(int[][] img, int i, int j) {
		
		float result = 1;
		
		while(img[j][i] == 1) {
			if(img[j][i-1] == 1) {
				i = i-1;
				result++;
			} else if(img[j-1][i-1] == 1) {
				j = j-1;
				i = i-1;
				result++;
			} else if(img[j+1][i-1] == 1) {
				j = j+1;
				i = i-1;
				result++;
			} else if(img[j-2][i-1] == 1) {
				j = j-1;
				i = i-1;
				result++;
			} else if(img[j+2][i-1] == 1) {
				j = j+1;
				i = i-1;
				result++;
			} else if(img[j-3][i-1] == 1) {
				j = j-1;
				i = i-1;
				result++;
			} else if(img[j+3][i-1] == 1) {
				j = j+1;
				i = i-1;
				result++;
			} else if(img[j-4][i-1] == 1) {
				j = j-1;
				i = i-1;
				result++;
			} else if(img[j+4][i-1] == 1) {
				j = j+1;
				i = i-1;
				result++;
			} else if(img[j-1][i-1] == 1) {
				j = j-1;
				i = i-1;
				result++;
			} else if(img[j+5][i-5] == 1) {
				j = j+1;
				i = i-1;
				result++;
			} else {
				break;
			}
		}
		return result;
	}
	
	public static void printInfo(ArrayList<Seed> list, float mm) {
		int i = 1;
		float size;
		for(Seed s: list) {
			size = (s.getLength() / mm) * 10;
			System.out.println("Seed " + i +": ");
			System.out.println("Location: " + s.getStartLocation());
			System.out.println("\tPixel length: " + s.getLength() + " px");
			System.out.println("\tLength of seed: " + size + " mm");
			i++;
		}
	}

}
