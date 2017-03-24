#include "../interface/CommonTools.h"

bool CommonTools::hardwareSimulation=true;

/* Function which simulate the HardWare representation of the values : manage UNSIGNED and SIGNED (2's complement) overflows and accuracy according to the available dynamic of the binary word */
double CommonTools::binning(double fNumber, int nMSBpowOfTwo, int nBits, HW_SIGN_TYPE signType)
{
  if (!hardwareSimulation)
    //If the Hardware binning simulation is not asked, return directly the original number
    return fNumber;

  if (signType == UNSIGNED && fNumber < 0)
    {
      //Bad interpretation, a negative number is stored in an UNSIGNED format (sign lost)
      fNumber = -fNumber;
    }
  
  int nLSBpowOfTwo;
	
  //Process the power of two of the LSB for the binary representation
  if (signType == UNSIGNED)
    {
      //If UNSIGNED
      nLSBpowOfTwo = nMSBpowOfTwo - (nBits-1);
    }	
  else
    {
      //If SIGNED, 1 bit is used for the sign
      nLSBpowOfTwo = nMSBpowOfTwo - (nBits-2);
    }

  /* Accuracy Simulation */

  //Divide the number by the power of two of the LSB => the integer part of the new number is the value we are looking for
  fNumber = fNumber / pow(2, nLSBpowOfTwo);
	
  //Remove the fractionnal part by rounding down (for both positive and negative values), this simulate the HW truncature
  fNumber = floor(fNumber);
	
  //Multiply the number by the power of two of the LSB to get the correct float value
  fNumber = fNumber * pow(2, nLSBpowOfTwo);


  double fBinnedNumber = fNumber;

  /* Overflow Simulation */

  if (signType == UNSIGNED)
    {
      //If the number is in UNSIGNED representation
      fNumber = fmod(fNumber, pow(2, nMSBpowOfTwo+1));
    }
  else
    {
      //If the number is in SIGNED representation (2's complement)
      
      double fTempResult = fNumber - pow(2, nMSBpowOfTwo+1); //substract the possible range to the number

      if (fTempResult >= 0)
        {
          //If there is an overflow, it's a positive one
          fNumber = fmod(fTempResult, pow(2, nMSBpowOfTwo+2)) - pow(2, nMSBpowOfTwo+1);
        }
      else
        {
          //If there is an overflow, it's a negative one (2's complement format has an asymetric range for positive and negative values)
          fNumber = fmod(fTempResult + pow(2, nLSBpowOfTwo), pow(2, nMSBpowOfTwo+2)) - pow(2, nLSBpowOfTwo) + pow(2, nMSBpowOfTwo+1);
        }
    }

  //If the new number is different from the previous one, an HW overflow occured
  if (fNumber != fBinnedNumber)
    {
      cout<<"WARNING HW overflow for the value : "<<fBinnedNumber<<" resulting value : "<<fNumber<<" (diff= "<<fBinnedNumber-fNumber<<")"<<endl;
    }
	
  return fNumber;
}

/*Bitwise emulation of the firmware CORDIC module (process the translation between cartesian and polar coordinates)*/
void CommonTools::binCordic(double X, double Y, double &result_R, double &result_PHI)
{

  if (!hardwareSimulation){
    //Floating point computing
    result_R = sqrt(X*X + Y*Y);
    result_PHI = atan(Y/X);
    return;
  }

  //Number of iterations for the CORDIC algorithm
  int nIter = 17;
  
  //Factor used to process R (depends on the number of iterations)
  double fScaleFactor = 0.60725289583206176757812500;
  
  //Micro-rotation values for the iterations
  double tfCorrec[17] = {
    0.78540039062500000000,
    0.46365356445312500000,
    0.24497985839843750000,
    0.12435913085937500000,
    0.06242370605468750000,
    0.03123474121093750000,
    0.01562500000000000000,
    0.00781250000000000000,
    0.00390625000000000000,
    0.00195312500000000000,
    0.00097656250000000000,
    0.00048828125000000000,
    0.00024414062500000000,
    0.00012207031250000000,
    0.00006103515625000000,
    0.00003051757812500000,
    0.00001525878906250000};

  //Initialization
  double nextX, nextY;
  double Xcordic = binning(X, 6, 18, SIGNED);
  double Ycordic = binning(Y, 6, 18, SIGNED);
  double angle = 0.0;

  //Iterations
  for (int i =0; i<nIter; i++)
    {
      if (Ycordic >= 0.0)
	{
	  nextX = Xcordic + binning(Ycordic * pow(2,-i), 7, 26, SIGNED);
	  nextY = Ycordic - binning(Xcordic * pow(2,-i), 7, 26, SIGNED);
	  angle += tfCorrec[i];
	}
      else
	{
	  nextX = Xcordic - binning(Ycordic * pow(2,-i), 7, 26, SIGNED);
	  nextY = Ycordic + binning(Xcordic * pow(2,-i), 7, 26, SIGNED);
	  angle -= tfCorrec[i];
	}
      Xcordic = binning(nextX, 7, 26, SIGNED);
      Ycordic = binning(nextY, 7, 26, SIGNED);
    }
	
  //Apply the scale factor
  result_R = Xcordic * fScaleFactor;
	
  //Binning of the results
  result_R = binning(result_R, 6, 18, SIGNED);
  result_PHI = binning(angle, 0, 18, SIGNED);
}
