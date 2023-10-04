TGraph **peakg;
TGraph **peakg2;
TGraph **peakg3;
const int nbExponentialPar=4;
Double_t resp(TGraph *gr, Double_t *x, Double_t *par) {return gr->Eval(x[0]+par[0])*par[1];}

Double_t respf(Double_t *x, Double_t *par) { return resp(peakg[(int)par[2]],x,par) ;}
//Double_t respf2(Double_t *x, Double_t *par) { return resp(peakg2[(int)par[2]],x,par) ;}
//Double_t respf3(Double_t *x, Double_t *par) { return resp(peakg3[(int)par[2]],x,par) ;}

Double_t expf1(Double_t *x, Double_t *par) { return TMath::Exp(par[0]+par[1]*x[0]) + TMath::Exp( par[2]+par[3]*x[0] );}

Double_t ex_respf1(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 1 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}
Double_t ex_respf2(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 2 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
  return rValue;
}
Double_t ex_respf3(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 3 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
  return rValue;
}

Double_t ex_respf4(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 4 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}
Double_t ex_respf5(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 5 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}
Double_t ex_respf6(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 6 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}
Double_t ex_respf7(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 7 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}
Double_t ex_respf8(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 8 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}

Double_t ex_respf9(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 9 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}
Double_t ex_respf10(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 10 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}
Double_t ex_respf11(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 11 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}
Double_t ex_respf12(Double_t *x, Double_t *par) {
  double rValue = expf1(x,par);
  for(Int_t j = 0 ; j < 12 ;j++)
    rValue += resp(peakg[j],x,par+nbExponentialPar+2*j);
    return rValue;
}
