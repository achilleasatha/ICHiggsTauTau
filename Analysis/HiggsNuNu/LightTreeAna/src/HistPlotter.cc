#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/HistPlotter.h"
#include <iostream>
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <map>
//#include "CommonTools/Utils/interface/TFileDirectory.h"
#include <boost/algorithm/string.hpp>
#include "TDirectory.h"
#include "TFile.h"

namespace ic{

  LTPlotElement::LTPlotElement(){
    unit_="GeV";
    in_stack_=false;
  };

  LTPlotElement::~LTPlotElement(){ ;};

  void LTPlotElement::ApplyStyle(){
    hist_ptr_->SetName(legname_.c_str());
    hist_ptr_->Scale(scale_);
    if(draw_marker_){
      hist_ptr_->SetMarkerColor(marker_color_);
      hist_ptr_->SetMarkerStyle(marker_style_);
      hist_ptr_->SetMarkerSize(marker_size_);
      drawopts_+="P";
      legopts_+="p";
    }
    if(draw_fill_){
      hist_ptr_->SetFillColor(fill_color_);
      hist_ptr_->SetFillStyle(fill_style_);
      hist_ptr_->SetLineColor(line_color_);
      drawopts_+="hist";
      legopts_+="f";
    }
    if(draw_line_){
      hist_ptr_->SetLineColor(line_color_);
      hist_ptr_->SetLineStyle(line_style_);
      hist_ptr_->SetLineWidth(line_width_);
      legopts_+="l";
    }
    if(draw_stat_error_y_==true){
      drawopts_+="E1";
      legopts_+="e";
    }
  };

  void HistPlotter::SetMCStackStyle(ic::LTPlotElement* ele) {
    ele->set_fill_color(ele->color());
    ele->set_fill_style(1001);
    ele->set_draw_fill(true);
    ele->set_draw_marker(false);
    ele->set_draw_line(false);
    ele->set_line_width(2);
    ele->set_draw_stat_error_y(false);
    ele->set_line_color(1);
    return;
  }
  void HistPlotter::SetSignalStyle(ic::LTPlotElement* ele) {
    ele->set_fill_style(1001);
    ele->set_draw_fill(true);
    ele->set_draw_marker(false);
    ele->set_draw_line(true);
    ele->set_draw_stat_error_y(false);
    ele->set_draw_fill_in_legend(false);
    ele->set_line_style(11);
    ele->set_fill_color(0);
    ele->set_line_color(ele->color());
    ele->set_line_width(2);
    return;
  }
  void HistPlotter::SetDataStyle(ic::LTPlotElement* ele) {
    ele->set_marker_color(1);
    ele->set_line_color(1);
    ele->set_fill_color(1);
    ele->set_fill_style(0);
    ele->set_draw_fill(false);
    ele->set_line_width(2);
    ele->set_draw_marker(true);
    ele->set_draw_line(true);
    ele->set_marker_style(20);
    ele->set_draw_stat_error_y(true);
    ele->set_marker_size(1.1);
    return;
  }

  HistPlotter::HistPlotter(std::string name) : LTModule(name){
  };

  HistPlotter::~HistPlotter(){ ;};

  int HistPlotter::Init(TFile* fs){
    fs_=fs;
    std::cout<<"Initialisation info for "<<module_name_<<":"<<std::endl;

    return 0;
  };

  int HistPlotter::Run(LTFiles* filemanager){
    std::cout<<module_name_<<":"<<std::endl;
    TFile* file=fs_;

    //GET DIRECTORY TO WRITE TO
    TDirectory* writedir;
    if(dirname_==""){
      writedir=file->mkdir("controlplots");
    }
    else if(!fs_->GetDirectory(dirname_.c_str())){
      writedir=file->mkdir(dirname_.c_str());
    }
    else{
      writedir=file->GetDirectory(dirname_.c_str());
    }
    writedir->cd();

    //SETUP STYLE
    for(unsigned iElement=0;iElement<elements_.size();iElement++){
    }


    //LOOP OVER ALL THE VARIABLES TO PLOT
    for(unsigned iShape=0;iShape<shapes_.size();iShape++){
      std::cout<<"  Drawing plot for "<<shapes_[iShape]<<std::endl;
      THStack *stack=new THStack("stack","stacked plots");
      bool stackempty=true;
      stack->SetTitle(shapes_[iShape].c_str());
      //stack->GetXaxis()->SetTitle(shapes_[iShape].c_str());

      //EXTRACT ALL THE HISTOS AND PUT THEM IN STACKED OR UNSTACKED GROUPS
      std::cout<<"    Getting histograms.."<<std::endl;
      for(unsigned iElement=0;iElement<elements_.size();iElement++){
	if(elements_[iElement].sample()==""){
	  std::cout<<"ERROR: Element with empty name exiting with status 1"<<std::endl;
	  return 1;
	}
	//std::cout<<"      "<<elements_[iElement].sample()<<std::endl;
	if(!fs_->GetDirectory(elements_[iElement].sample().c_str())){
	  std::cout<<"ERROR: No directory with name: "<<elements_[iElement].sample()<<std::endl;
	  std::cout<<"Exiting with status 1"<<std::endl;
	  return 1;
	}
	writedir->cd();
	TH1F* histo =dynamic_cast<TH1F*>(file->Get((elements_[iElement].sample()+"/"+shapes_[iShape]).c_str()));
	writedir->cd();
	elements_[iElement].set_hist_ptr(histo);
	
	//SETUP STYLE
	if(elements_[iElement].is_data()){
	  SetDataStyle(&(elements_[iElement]));
	}
	else if(elements_[iElement].in_stack()){
	  SetMCStackStyle(&(elements_[iElement]));
	  stackempty=false;
	}
	else{
	  SetSignalStyle(&(elements_[iElement]));
	}
	elements_[iElement].ApplyStyle();

	//ADD STACKED HISTOS TO STACK
	if((!elements_[iElement].is_data())&&elements_[iElement].in_stack()){
	  stack->Add(elements_[iElement].hist_ptr());
	}
      }

      //SETUP THE CANVAS
      TCanvas *c1=new TCanvas(shapes_[iShape].c_str(),shapes_[iShape].c_str());
      c1->cd();
      bool first=true;
      double ymax=0;
      if(!stackempty) ymax=stack->GetMaximum();
      for(unsigned iElement=0;iElement<elements_.size();iElement++){
	if(elements_[iElement].hist_ptr()->GetMaximum()>ymax) ymax=elements_[iElement].hist_ptr()->GetMaximum();
      }
      
      //DRAW THE STACK AND ALL THE UNSTACKED HISTOS INCLUDING THE DATA
      if(!stackempty){ 
	std::cout<<"    Drawing Stack.."<<std::endl;
	if(first){
	  stack->SetMaximum(ymax);
	  stack->Draw("hist");
	  c1->Update();
	  first=false;
	}
	else stack->Draw("histsame");
      }
      std::cout<<"    Drawing Unstacked.."<<std::endl;
      for(unsigned iElement=0;iElement<elements_.size();iElement++){
	if(!(elements_[iElement].in_stack())){
	  if(first){
	    elements_[iElement].hist_ptr()->Draw(elements_[iElement].drawopts().c_str());
	    elements_[iElement].hist_ptr()->GetYaxis()->SetRangeUser(0.,ymax+10);
	    elements_[iElement].hist_ptr()->Draw(elements_[iElement].drawopts().c_str());
	    c1->Update();
	    first=false;
	  }
	  else elements_[iElement].hist_ptr()->Draw(("same"+elements_[iElement].drawopts()).c_str());
	}
      }

      //SETUP AND DRAW THE LEGEND
      TLegend* leg =new TLegend(0.6,0.6,0.9,0.9);
      for(unsigned iElement=0;iElement<elements_.size();iElement++){
	leg->AddEntry(elements_[iElement].hist_ptr(),elements_[iElement].hist_ptr()->GetName(),elements_[iElement].legopts().c_str());
      }
      leg->Draw("same");
      c1->Update();
      
      //WRITE TO FILE
      writedir->cd();
      c1->Write();
      c1->Close();
    }
    writedir->Close();

    return 0;
  };

}
