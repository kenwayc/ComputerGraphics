// generated by Fast Light User Interface Designer (fluid) version 1.0107

#ifndef SceneviewInterface_h
#define SceneviewInterface_h
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include "SceneviewUI.h"
#include <FL/Fl_Button.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Output.H>

class SceneviewInterface {
public:
  Fl_Double_Window* make_window();
  Fl_Double_Window *m_sceneviewWindow;
private:
  void cb_Load_i(Fl_Button*, void*);
  static void cb_Load(Fl_Button*, void*);
public:
  Fl_Check_Button *m_bLighting;
private:
  void cb_m_bLighting_i(Fl_Check_Button*, void*);
  static void cb_m_bLighting(Fl_Check_Button*, void*);
public:
  Fl_Check_Button *m_bInteractive;
private:
  void cb_m_bInteractive_i(Fl_Check_Button*, void*);
  static void cb_m_bInteractive(Fl_Check_Button*, void*);
public:
  Fl_Output *m_txtStatus;
  Fl_Check_Button *m_bIBar;
private:
  void cb_m_bIBar_i(Fl_Check_Button*, void*);
  static void cb_m_bIBar(Fl_Check_Button*, void*);
public:
  Fl_Check_Button *m_bIBarHelp;
private:
  void cb_m_bIBarHelp_i(Fl_Check_Button*, void*);
  static void cb_m_bIBarHelp(Fl_Check_Button*, void*);
public:
  SceneviewInterface();
private:
  SceneviewUI sceneviewUI;
public:
  UIInterface * getUI();
private:
  std::string strSceneName;
};
#endif
