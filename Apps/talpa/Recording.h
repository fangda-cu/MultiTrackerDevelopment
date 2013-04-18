//
//  Recording.h
//  BASim
//
//  Created by Fang Da on 4/18/13.
//
//

#ifndef BASim_Recording_h
#define BASim_Recording_h

#include <string>
#include <eltopo.h>
#include <fstream>
#include <surftrack.h>
#include <sstream>
#include <iostream>

class Recording
{
public:
    Recording() : m_recording_name("rec"), m_current_frame(0), m_current_step(0), m_recording(false) { }
    
    void setRecordingName(const std::string & name) { m_recording_name = name; }
    const std::string & recordingName() const { return m_recording_name; }
    
    void setCurrentFrame(int frame) { m_current_frame = frame; m_current_step = 0; m_of.close(); m_if.close(); }
    int currentFrame() const { return m_current_frame; }
    void setCurrentStep(int step) { m_current_step = step % m_step_pos.size(); }
    int currentStep() const { return m_current_step; }
    
    void recordSurfTrack(const ElTopo::SurfTrack & st);
    void loadRecording(ElTopo::SurfTrack & st, int next = 0);
    
    void turnOnRecording() { m_recording = true; m_playback = false; }
    void turnOffRecording() { m_recording = false; }
    bool isRecording() const { return m_recording; }
    
    void turnOnPlayback() { m_playback = true; m_recording = false; }
    void turnOffPlayback() { m_playback = false; }
    bool isPlaybackOn() const { return m_playback; }
    
public:
    static void writeSurfTrack(std::ostream & os, const ElTopo::SurfTrack & st);
    static void readSurfTrack(std::istream & is, ElTopo::SurfTrack & st);
    
public:
    std::ostream & log() { return m_log; }
    
protected:
    std::string m_recording_name;
    int m_current_frame;  // frame number
    int m_current_step;   // step number within the frame
    
    bool m_recording;
    bool m_playback;
    
    std::vector<std::streampos> m_step_pos;
    
    std::ofstream m_of;
    std::ifstream m_if;
    
    std::stringstream m_log;
};

extern Recording g_recording;



#endif
