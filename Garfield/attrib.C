void   Attrib    ( Double_t siz=20, Double_t ts=20, Double_t yts=20,
                   Double_t xto=1.5, Double_t yto=1.0 );
//void   AttribH   ( TH1D *hist, Int_t col, Double_t alpha=0.6 );
void   AttribH   ( TH1F *hist, Int_t col, Double_t alpha=0.6, Double_t marker=0.6);
void   AttribL   ( TLegend *leg, Double_t mgn=0.15, Double_t size=24, Int_t font=23 );
//void   AttribG   ( TGraph *g, Int_t col, Int_t msty, Int_t lsty );
void   AttribG   ( TGraph *g, Int_t col, Int_t msty, Int_t lsty, Int_t lwdy,Double_t marker=0.6);
void   AttribF   ( TF1 *func, Int_t col, Int_t lwid, Int_t lsty );
void   MarkerDraw( Double_t x, Double_t y,
                   Int_t sty, Int_t col, Double_t si, Bool_t NDC = kFALSE );
void   TextDraw  ( Double_t x, Double_t y, TString t,
                   Double_t s, Int_t col = kBlack, Bool_t fNDC = kFALSE );
void   TextDraw2 ( Double_t x, Double_t y, TString t,
                   Double_t s, Int_t col = kBlack, Bool_t fNDC = kFALSE,
                   Double_t fid = 23 );
TLine* LineDraw  ( Double_t x1, Double_t x2, Double_t y1, Double_t y2,
                   Int_t col, Int_t sty, Int_t wid );
void   SetTag    ( void );

// example  Attrib( 20, 24, 24, 0, 0 );
void Attrib( Double_t siz, Double_t xts, Double_t yts,
             Double_t xto, Double_t yto )
{
  //gStyle->SetTitleFontSize( 1000000.0*siz );
  gStyle->SetTitleSize  ( xts, "X" );
  gStyle->SetTitleSize  ( yts, "Y" );
  gStyle->SetLabelSize  ( siz, "X" );
  gStyle->SetLabelSize  ( siz, "Y" );
  gStyle->SetTitleOffset( xto, "X" );
  gStyle->SetTitleOffset( yto, "Y" );
  gStyle->SetNdivisions( 505, "X" );
  gStyle->SetNdivisions( 505, "Y" );
  gStyle->SetOptStat( 0 );
  gStyle->SetPalette( 55 );
  gStyle->SetNumberContours( 100 );
  gStyle->SetTextFont( 23 );
  gStyle->SetLabelFont( 23, "XY" );
  gStyle->SetTitleFont( 23, "XY" );
}

void SetTag( void )
{
  TText *note = new TText( 0.21,0.74,"Should be corrected");
  note->SetNDC();
  note->SetTextSize(80);
  note->SetTextColorAlpha( kRed, 0.5 );
  note->SetTextAngle(-45);
  note->Draw();
}

void MarkerDraw( Double_t x, Double_t y,
                 Int_t sty, Int_t col, Double_t si, Bool_t NDC )
{
  TMarker *mk = new TMarker( x, y, sty );
  mk->SetMarkerColor( col );
  mk->SetMarkerSize ( si  );
  if( NDC ) mk->SetNDC();
  mk->Draw();
}

void TextDraw( Double_t x, Double_t y, TString t,
               Double_t s, Int_t col, Bool_t fNDC )
{
  TLatex *text = new TLatex( x, y, t );
  text->SetTextSize( s );
  text->SetTextColor( col );
  if( fNDC ) text->SetNDC();
  text->Draw();
}

void TextDraw2( Double_t x, Double_t y, TString t,
                Double_t s, Int_t col, Bool_t fNDC, Double_t fid )
{
  TLatex *text = new TLatex( x, y, t );
  text->SetTextFont( fid );
  text->SetTextSize( s );
  text->SetTextColor( col );
  if( fNDC ) text->SetNDC();
  text->Draw();
}

//void AttribH( TH1D *hist, Int_t col, Double_t alpha )
void AttribH( TH1F *hist, Int_t col, Double_t alpha , Double_t marker)
{
  hist->SetMarkerStyle( 24 );
  hist->SetMarkerSize ( marker );
  hist->SetMarkerColor( col );
  hist->SetLineColor( col );
  hist->SetFillStyle( 1000 );
  hist->SetFillColorAlpha( col, alpha);
}

void AttribL( TLegend *leg, Double_t mgn, Double_t size, Int_t font )
{
  leg->SetFillStyle ( 0 );
  leg->SetBorderSize( 0 );
  leg->SetTextSize  ( size );
  leg->SetTextFont  ( font );
  leg->SetMargin    ( mgn );
}

void AttribG( TGraph *g, Int_t col, Int_t msty, Int_t lsty , Int_t lwdy, Double_t marker)
{
  g->SetMarkerStyle( msty );
  g->SetMarkerColor( col );
  g->SetLineStyle( lsty );
  g->SetLineColor( col );
  g->SetLineWidth( lwdy );
  g->SetMarkerSize ( marker );
}

void AttribF( TF1 *func, Int_t col, Int_t lwid, Int_t lsty )
{
  func->SetLineColor( col  );
  func->SetLineWidth( lwid );
  func->SetLineStyle( lsty );
}

TLine* LineDraw( Double_t x1, Double_t x2, Double_t y1, Double_t y2,
                 Int_t col, Int_t sty, Int_t wid )
{

  Double_t x1c = (x1 <= x2)? x1 : x2;
  Double_t x2c = (x1 <= x2)? x2 : x1;
  Double_t y1c = (x1 <= x2)? y1 : y2;
  Double_t y2c = (x1 <= x2)? y2 : y1;

  TLine *line = new TLine( x1c, y1c, x2c, y2c );
  line->SetLineColor ( col );
  line->SetLineStyle( sty );
  line->SetLineWidth( wid );
  line->Draw();

  return line;
}
  