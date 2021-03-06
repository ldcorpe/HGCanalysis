#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TAxis.h"
#include "TStyle.h"

#include <iostream>

using namespace std;

//
std::pair<float,float> getLocalCoords(int cell, float cellSize, float h, float bl, float tl)
{
  //linear parameterization of the trapezoid
  float a=2*h/(tl-bl);
  float b=-h*(tl+bl)/(tl-bl);
  
  //find the y-row iteratively
  int maxKy=floor(2*h/cellSize);
  int ky(0),testCell(0);
  for(int iky=0; iky<maxKy; iky++)
    {
      int deltay( floor( (iky*cellSize-h-b)/(a*cellSize) ) );
      if(testCell+deltay > cell) break;
      testCell+=deltay;
      ky++;
    }
  
  //find the x-column
  int kx=cell-testCell;

  //all done here (return centered at cell)
  return std::pair<float,float>((kx+0.5)*cellSize,(ky+0.5)*cellSize-h);
}


//
int assignCell(float x, float y, float cellSize, float h, float bl, float tl)
{
  //linear parameterization of the trapezoid
  float a=2*h/(tl-bl);
  float b=-h*(tl+bl)/(tl-bl);
  
  //this is the cell # in the row and column
  int kx=floor( fabs(x)/cellSize );
  int ky=floor((y+h)/cellSize);
  
  //find the cell sequentially in the trapezoid
  //notice the arithmetic sum can't be used as \sum floor(x) != floor( \sum x )
  int icell(0);
  for(int iky=0; iky<ky; iky++)
    icell += floor( (iky*cellSize-h-b)/(a*cellSize) );
  icell += kx;

  //all done here
  return icell;
}


void divideTrapezoid() 
{
  float cellSize(9.5);
  //   float bl=58.3041 ; float tl=283.979 ; float h=639.934 ; float y=-18.1448 ; float x=-568.768;
  //   float bl=58.4537 ; float tl=284.45  ; float h=640.844 ; float y=-14.4645 ; float x=-570.561;
  //   float bl=58.6438 ; float tl=285.048 ; float h=642.001 ; float y=-9.78826 ; float x=-572.838;
  //   float bl=58.8339 ; float tl=285.646 ; float h=643.157 ; float y=-5.11207 ; float x=-575.115;
  //   float bl=59.024  ; float tl=286.244 ; float h=644.314 ; float y=-0.435881; float x=-577.392;
  //   float bl=59.2141 ; float tl=286.842 ; float h=645.47  ; float y=4.24031  ; float x=-579.669;
  //   float bl=59.4042 ; float tl=287.44  ; float h=646.627 ; float y=8.9165   ; float x=-581.946;
  //   float bl=59.5943 ; float tl=288.038 ; float h=647.784 ; float y=13.5927  ; float x=-584.223;
  //float bl=57.556  ; float tl=281.626 ; float h=635.383 ;      float x=-30.902  ; float y=-561.339;
  //  float bl=56.2377; float tl=265.942 ; float h=594.645;      float x=-156.75; float y=513.767;
  float bl=60.4919; float tl=285.782 ; float h=638.841;      float x=432.25; float y=606.631;


// Layer 26 Sector 2
// h=638.841 b=60.4919 t=285.782 cell=9.5
//   (x0,y0,z0)=(981.908,-2.15798e-12,3423.55)
// rho=981.908phi=4.59163e-41
//   |xx,xy,xz|=|1.22461e-161-6.12303e-16|
//   |yx,yy,yz|=|-1-3.67382e-16-3.67382e-16|
//   |zx,zy,zz|=|1.22461e-16-2.44921e-161|
//  Can't accumulate @ (432.25 606.631)
// 432.25,606.631,432.25,606.631
// 2160
  //23.75,572.409
  
  //x=cellSize*(0+1./2);
  //  y=-h+cellSize*(0+1./2);

  //linear parameterization of the trapezoid
  float a=2*h/(tl-bl);
  float b=-h*(tl+bl)/(tl-bl);

  //this is the cell # in the row and column
  int kx=TMath::Floor( fabs(x)/cellSize );
  int ky=TMath::Floor((y+h)/cellSize);

  //find the cell sequentially in the trapezoid
  //notice the arithmetic sum can't be used as \sum floor(x) != floor( \sum x )
  int icell=assignCell(x,y,cellSize,h,bl,tl);
  std::cout <<"(" << kx << "," << ky << ") -> " << icell << std::endl;
  
  //  std::pair<float,float> localxy=getLocalCoords(icell,cellSize,h,bl,tl);
  std::pair<float,float> localxy=getLocalCoords(2160,cellSize,h,bl,tl);
  std::cout << localxy.first << "," << localxy.second << std::endl;
  std::cout << x << "," << y << " <- orig" << std::endl;

  //show all of this

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *c=new TCanvas("c","c",600,600);
  c->cd();

  TGraph *gr=new TGraph;
  gr->SetPoint(0,-bl,-h);
  gr->SetPoint(1,bl,-h);
  gr->SetPoint(2,tl,h);
  gr->SetPoint(3,-tl,h);
  gr->SetPoint(4,-bl,-h);
  gr->SetFillColor(kGray);
  gr->Draw("afl");
  gr->GetXaxis()->SetTitle("x [mm]");
  gr->GetYaxis()->SetTitle("y [mm]");

  TGraphErrors *grerr=new TGraphErrors;
  grerr->SetMarkerColor(kRed);
  grerr->SetFillStyle(0);
  grerr->SetLineColor(kRed);
  grerr->SetLineWidth(4);
  grerr->SetMarkerStyle(20);
  grerr->SetPoint(0,x,y);
  grerr->SetPointError(0,cellSize,cellSize);
  grerr->SetPoint(1,-x,y);
  grerr->SetPointError(1,cellSize,cellSize);
  grerr->Draw("e2p");

    
  Int_t ncells(0);
  Int_t Ny = TMath::Floor(2*h/cellSize);
  TLatex *ptcell=new TLatex;
  ptcell->SetTextFont(42);
  ptcell->SetTextAlign(12);
  ptcell->SetTextSize(0.015);
  for(Int_t i=0; i<Ny; i++)
    {
      float iy=i*cellSize-h;

      float xMax=(iy-b)/a;
      Int_t Nx = TMath::Floor(xMax/cellSize);

      for(Int_t j=0; j<Nx; j++)
	{
	  float ix=j*cellSize;

	  if(fabs(x)>ix && 
	     fabs(x)<ix+cellSize &&
	     y>iy && 
	     y<iy+cellSize) 
	    cout << ncells << endl;

	  TGraph *pad=new TGraph;
	  pad->SetPoint(0,ix,iy);
	  pad->SetPoint(1,ix+cellSize,iy);
	  pad->SetPoint(2,ix+cellSize,iy+cellSize);
	  pad->SetPoint(3,ix,iy+cellSize);
	  pad->SetPoint(4,ix,iy);
	  pad->Draw("l");

	  if(ncells>icell-10 && ncells<icell+10){
	    char buf[200];
	    sprintf(buf,"%d",ncells);
	    ptcell->DrawLatex(ix+cellSize/4,iy+cellSize/4,buf);	  
	  }
	  ncells++;
	}
    }


  TPaveText *pt=new TPaveText(0.12,0.9,0.5,0.95,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  char buf[200];
  sprintf(buf,"%3.1f x %3.1f mm^{2} cells",cellSize,cellSize);
  pt->AddText(buf);
  sprintf(buf,"# cells / 10 deg sub-sector = %d",ncells);
  pt->AddText(buf);
  pt->SetTextAlign(12);
  pt->Draw();
}
