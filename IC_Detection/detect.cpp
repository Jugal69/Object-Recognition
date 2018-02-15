//
// detect.cpp : Detect integrated circuits in printed circuit board (PCB) images.
//
// Based on skeleton code by D. Crandall, Spring 2018
//
// Authors -->  Shreyas Fadnavis(ID:shfadn) & Khusaal Giri(ID:kngiri)
//
//

#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc.
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose.

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);

    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedBox class may be helpful!
//  Feel free to modify.
//
class DetectedBox {
public:
  int row, col, width, height;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<DetectedBox> &ics)
{
  ofstream ofs(filename.c_str());

  for(vector<DetectedBox>::const_iterator s=ics.begin(); s != ics.end(); ++s)
    ofs << s->row << " " << s->col << " " << s->width << " " << s->height << " " << s->confidence << endl;
}

// Function that outputs a visualization of detected boxes
void  write_detection_image(const string &filename, const vector<DetectedBox> &ics, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];

  for(int p=0; p<3; p++)
    {
      output_planes[p] = input;
      for(vector<DetectedBox>::const_iterator s=ics.begin(); s != ics.end(); ++s)
	overlay_rectangle(output_planes[p], s->row, s->col, s->row+s->height-1, s->col+s->width-1, p==2?255:0, 2);
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  int padding_bit_row= row_filter.cols()-1,temp,padding_bit_col= col_filter.rows()-1;
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane input_temp(input.rows()+ padding_bit_col,input.cols()+  padding_bit_row);
  SDoublePlane output_temp1(input.rows()+  padding_bit_col,input.cols());

  //Padding the image boundaries with k bits from all the sides
  for (int i=0;i<input.rows();i++){
	for(int j=0;j<input.cols();j++)
		input_temp[i+(padding_bit_col%2?padding_bit_col/2+1:padding_bit_col/2)][j+(padding_bit_row%2?padding_bit_row/2+1:padding_bit_row/2)]=input[i][j];
	}

  //Flipping the row filter
  for(int i=0;i<row_filter.cols()/2;i++)
  {
	temp=row_filter[0][i];
	row_filter[0][i]=row_filter[0][row_filter.cols()-i-1];
	row_filter[0][row_filter.cols()-i-1]=temp;
  }

  //Convolution with a row filter
  for(int i=0;i<input_temp.rows();i++){
  	for(int j=0;j<output_temp1.cols();j++){
  		for(int k=0;k<row_filter.cols();k++)
  			output_temp1[i][j]+=row_filter[0][k]*input_temp[i][k+j];
	  }
	}

  //Flipping the column filter
  for(int i=0;i<col_filter.rows()/2;i++)
  {
	temp=col_filter[i][0];
	col_filter[i][0]=col_filter[col_filter.rows()-i-1][0];
	col_filter[col_filter.rows()-i-1][0]=temp;
  }

  //Convolution with a column vector
  for(int i=0;i<output_temp1.cols();i++){
  		for(int j=0;j<output.rows();j++){
  			for(int k=0;k<col_filter.rows();k++)
  				output[j][i]+=col_filter[k][0]*output_temp1[k+j][i];
	  }
	}
	return output;
}

// Convolve an image with a  convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
  int padding_bit_row= filter.cols()-1,padding_bit_col= filter.rows()-1;
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane input_temp(input.rows()+ padding_bit_col,input.cols()+  padding_bit_row);
  SDoublePlane flipped_filter(filter.rows(),filter.cols());

  //Flipping the filter horizontally & vertically
  for (int i=filter.rows()-1;i>=0;i--){
	for(int j=filter.cols()-1;j>=0;j--)
		flipped_filter[filter.rows()-i-1][filter.cols()-j-1]=filter[i][j];
	}

  //Padding the image boundaries with k bits from all the sides
  for (int i=0;i<input.rows();i++){
	for(int j=0;j<input.cols();j++)
		input_temp[i+(padding_bit_col%2?padding_bit_col/2+1:padding_bit_col/2)][j+(padding_bit_row%2?padding_bit_row/2+1:padding_bit_row/2)]=input[i][j];
	}

  //Convolution step
  int k_limit=(padding_bit_col%2?padding_bit_col/2:padding_bit_col/2-1),l_limit=(padding_bit_row%2?padding_bit_row/2:padding_bit_row/2-1);
  int i_pad=(padding_bit_col%2?padding_bit_col/2+1:padding_bit_col/2),j_pad=(padding_bit_row%2?padding_bit_row/2+1:padding_bit_row/2);
  for(int i=i_pad;i<input_temp.rows()-padding_bit_col/2;i++){
	  for(int j=j_pad;j<input_temp.cols()-padding_bit_row/2;j++){
		  for (int k = -(padding_bit_col/2); k <=k_limit ; k++) {
			  for (int l = -(padding_bit_row/2); l <= l_limit; l++)
				  output[i - i_pad][j - j_pad] += input_temp[k + i][l + j] * flipped_filter[k + padding_bit_col/2][l + padding_bit_row/2];
		  }
	  }
  }
  return output;
}

// Apply a sobel operator to an image, returns the result
//
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  SDoublePlane output(input.rows(), input.cols());

  // Has been implemented in another function named as "convolve_seperable"

  return output;
}

// Apply an edge detector to an image, returns the binary edge map
//
SDoublePlane find_edges(const SDoublePlane &input, double thresh)
{
  const double pi=3.142;
  int count=0,diag_len,k;
  SDoublePlane output(input.rows(), input.cols());

  //Thresholding the input image
  for(int i=0;i<input.rows();i++){
  	for(int j=0;j<input.cols();j++){
  		if(input[i][j]<thresh)
  			input[i][j]=0;
  		else{
  			input[i][j]=255;
  			count++;
		  }
	  }
  }

  int n_r,bin_theta=1,n_theta=180/bin_theta;
  double u_c=floor(input.rows()/2),v_c=floor(input.cols()/2),del_theta=pi/n_theta,del_r;
  diag_len=ceil(sqrt(pow(u_c,2)+pow(v_c,2)));
  n_r=int(diag_len);
  del_r=2*diag_len/n_r;
  
  //Hough accumulator array of theta vs rho
  SDoublePlane accumulator(int(n_r)+1,n_theta),accumulator_x(int(n_r)+1,n_theta),accumulator_y(int(n_r)+1,n_theta);
  for(int i=0;i<input.rows();i++){
  	for(int j=0;j<input.cols();j++){
  		if(input[i][j]){
  			int x=i-u_c,y=j-v_c;
  			for(int k=0;k<n_theta;k++){
				  double theta=del_theta*k;
				  double r=x*cos(theta)+y*sin(theta);
				  int i_r=n_r/2+round(r/del_r);
        if(i_r>=0&&i_r<n_r){
          accumulator[i_r][k]++;
          accumulator_x[i_r][k]=i;
          accumulator_y[i_r][k]=j;
          if(i_r){                       //To avoid inaccurate point co-ordinates & rounding errors
            accumulator[i_r-1][k]++;
            accumulator_x[i_r-1][k]=i;
            accumulator_y[i_r-1][k]=j;
          }
          if(i_r!=accumulator.rows()){
            accumulator[i_r+1][k]++;
            accumulator_x[i_r+1][k]=i;
            accumulator_y[i_r+1][k]=j;
          }
        }
		  }
		}
	}
}
  
  //Thresholding of the accumulator array in the Hough space 
  for (int i=0;i<accumulator.rows();i++){
		for(int j=0;j<accumulator.cols();j++){
			if(accumulator[i][j]<input.rows())
				accumulator[i][j]=0;
			}
	}

  //Non-Maximal suppression 
  for(int i=1;i<accumulator.rows()-1;i++){
    for(int j=1;j<accumulator.cols()-1;j++){
      if(accumulator[i][j]){
        for(int k=-1;k<2;k++){
          for(int l=-1;l<2;l++)
            accumulator[i+k][j+l]=(accumulator[i+k][j+l]<accumulator[i][j]?0:accumulator[i+k][j+l]);
        }
      }
    }
  }

  //Peak finding in the Hough accumulator array
  SDoublePlane peaks(100000000,4);
  SDoublePlane ortho_check(600000000,3);
  int index=0,index_ortho=0;
  float t_theta=3,t_alpha=3,t_rho=3,t_l=0.4,d_theta;
  for(int i=0;i<accumulator.rows();i++){
  	for(int j=0;j<accumulator.cols();j++){
  		if(accumulator[i][j]){
  			for(int k=0;k<accumulator.rows();k++){
			  	for(int l=0;l<accumulator.cols();l++){
			  		if(accumulator[k][l]&&i!=k&&j!=l){
			  			double r_i=(2*i-n_r),r_k=(2*k-n_r);
						  if(del_theta*(180/pi)*abs(j-l)<t_theta && abs(r_i+r_k)<t_rho && (accumulator[i][j]-accumulator[k][l])<t_l*0.5*(accumulator[i][j]+accumulator[k][l])){
							  peaks[index][1]=del_theta*(180/pi)*0.5*(j+l);
							  peaks[index][2]=del_theta*(180/pi)*(j-l);
							  peaks[index][3]=r_i-r_k;
							  peaks[index++][0]=0.5*abs(r_i-r_k);
						  }
					  }
				  }
			  }
		  }
	  }
  }
  
  //Orthogonality check for all the peaks of the Hough Accumulator array
  for(int i=0;i<index;i++){
  	for(int j=i+1;j<index;j++){
  		if((del_theta*(180/pi)*abs(peaks[i][1]-peaks[j][1])-90)<t_alpha){
  			ortho_check[index_ortho][0]=i;
  			ortho_check[index_ortho][2]=sqrt(pow(peaks[i][2],2)+pow(peaks[j][2],2)+pow(del_theta*(180/pi)*abs(peaks[i][1]-peaks[j][1])-90,2)+4*pow(peaks[i][3],2)+4*pow(peaks[j][3],2));
  			ortho_check[index_ortho++][1]=j;
		  }
	  }
  }
 
  //Removal of Duplicated Rectangles
  SDoublePlane final_ortho(100000,3);
  int final_index=0,i=1,temp1=ortho_check[0][0],temp3=ortho_check[0][1];
  double temp2=ortho_check[0][2];
  while(i<index_ortho){
  	if(ortho_check[i][0]==ortho_check[i-1][0]){
  		if(ortho_check[i][2]<temp2){
			temp3=ortho_check[i][1];
  			temp2=ortho_check[i][2];
  		}
		  i++;
	  }
	else{
		final_ortho[final_index][0]=temp1;
		final_ortho[final_index][1]=temp3;
		final_ortho[final_index++][2]=temp2;
		temp1=ortho_check[i][0];
		temp3=ortho_check[i][1];
		temp2=ortho_check[i++][2];
	}
  }
  
  SDoublePlane final_ortho_col(50000,3);
  int final_index_col=0;
  i=1,temp1=final_ortho[0][0],temp3=final_ortho[0][1];
  temp2=final_ortho[0][2];
  while(i<final_index){
  	if(final_ortho[i][1]==final_ortho[i-1][1]){
  		if(final_ortho[i][2]<temp2){
			temp1=final_ortho[i][0];
  			temp2=final_ortho[i][2];
  		}
		  i++;
	  }
	else{
		final_ortho_col[final_index_col][0]=temp1;
		final_ortho_col[final_index_col][1]=temp3;
		final_ortho_col[final_index_col++][2]=temp2;
		temp1=final_ortho[i][0];
		temp3=final_ortho[i][1];
		temp2=final_ortho[i++][2];
	}
  }
  
  SDoublePlane rectangles(final_index_col,5);
  for(int i=0;i<final_index_col;i++){
    rectangles[i][0]=accumulator_x[int(peaks[final_ortho_col[i][0]][0])][int(peaks[final_ortho_col[i][1]][1])]*(1-cos(peaks[final_ortho_col[i][0]][1]*180/pi));
    rectangles[i][1]=accumulator_y[int(peaks[final_ortho_col[i][0]][0])][int(peaks[final_ortho_col[i][1]][1])]*(1-sin(peaks[final_ortho_col[i][1]][1]*180/pi));
    rectangles[i][2]=2*peaks[final_ortho_col[i][0]][0];
    rectangles[i][3]=2*peaks[final_ortho_col[i][1]][0];
    rectangles[i][4]=final_ortho_col[i][2]/10000;
  }
  return rectangles;
}

//
// This main file just outputs a few test images. You'll want to change it to do
//  something more interesting!
//
int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
  SDoublePlane output_image(input_image.rows(),input_image.cols());
  // test step 2 by applying mean filters to the input image
  int fil_row=3,fil_col=3;
  SDoublePlane mean_filter(fil_row,fil_col);
  /*for(int i=0; i<mean_filter.rows(); i++)
    for(int j=0; j<mean_filter.cols(); j++)
      mean_filter[i][j] = 1/9.0;*/

    //Gaussian Filter
	mean_filter[0][0]=1/16.0;
  mean_filter[0][1]=1/8.0;
	mean_filter[0][2]=1/16.0;
	mean_filter[1][0]=1/8.0;
	mean_filter[1][1]=1/4.0;
	mean_filter[1][2]=1/8.0;
	mean_filter[2][0]=1/16.0;
	mean_filter[2][1]=1/8.0;
	mean_filter[2][2]=1/16.0;

  int xrow=3,xcol=3,ycol=3,yrow=3;
  SDoublePlane sobel_xrow(1,xrow);
  SDoublePlane sobel_xcol(xcol,1);
  SDoublePlane sobel_yrow(1,yrow);
  SDoublePlane sobel_ycol(ycol,1);

  sobel_xrow[0][0]=1;
  sobel_xrow[0][1]=0;
  sobel_xrow[0][2]=-1;

  sobel_xcol[0][0]=1;
  sobel_xcol[1][0]=2;
  sobel_xcol[2][0]=1;

  sobel_ycol[0][0]=1;
  sobel_ycol[1][0]=0;
  sobel_ycol[2][0]=-1;

  sobel_yrow[0][0]=1;
  sobel_yrow[0][1]=2;
  sobel_yrow[0][2]=1;

  SDoublePlane output_image_x=convolve_separable(convolve_general(input_image,mean_filter),sobel_xrow,sobel_xcol);
  SDoublePlane output_image_y=convolve_separable(convolve_general(input_image,mean_filter),sobel_yrow,sobel_ycol);
  
  //Computing the magnitude of the gradient obtained from the sobel filter
  for(int i=0;i<output_image.rows();i++){
  	for(int j=0;j<output_image.cols();j++)
  		output_image[i][j]=sqrt(pow(output_image_x[i][j],2)+pow(output_image_y[i][j],2));
  }
  float avg=0;
  for (int i = 0; i < output_image.rows(); i++) {
	  for (int j = 0; j < output_image.cols(); j++)
		  avg += output_image[i][j];
  }
  avg /= (output_image.rows()*output_image.cols());
  
	SDoublePlane rectangles;
	rectangles=find_edges(output_image,avg*3);
	
  vector<DetectedBox> ics;
  write_detection_image("edges.png", ics,output_image);
  for(int i=0; i<rectangles.rows(); i++)
    {
      if((rectangles[i][0]<10||rectangles[i][0]>output_image.rows()-10)||(rectangles[i][1]<10||rectangles[i][1]>output_image.cols()-10))
        continue;
      DetectedBox s;
      s.row=int(rectangles[i][0]);
      s.col=int(rectangles[i][1]);
	    s.width = int(rectangles[i][2]);
      s.height = int(rectangles[i][3]);
      s.confidence = rectangles[i][4];
      ics.push_back(s);
    }

  write_detection_txt("detected.txt", ics);
  write_detection_image("detected.png", ics, output_image);
}
