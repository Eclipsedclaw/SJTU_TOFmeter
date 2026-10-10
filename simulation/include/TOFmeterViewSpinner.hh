//
// Slow automatic rotation of the current viewer in the Qt GUI: /tof/vis/spin <deg/s> turns the
// view about its up vector until the first mouse click, scroll or key press in a viewer
// (0 stops it). Only built when Geant4 has Qt (TOF_USE_QT, see CMakeLists.txt); the command
// exists in interactive GUI sessions only.
//

#ifndef TOFmeterViewSpinner_h
#define TOFmeterViewSpinner_h 1

#ifdef TOF_USE_QT

#include "globals.hh"

#include <QElapsedTimer>
#include <QObject>

class G4GenericMessenger;
class QTimer;

class TOFmeterViewSpinner : public QObject
{
  public:
    TOFmeterViewSpinner();
    ~TOFmeterViewSpinner() override;

    void Spin(G4double degreesPerSecond);   // 0 stops

  protected:
    // watches the whole application for the first interaction with a viewer
    bool eventFilter(QObject* watched, QEvent* event) override;

  private:
    void Step();

    QTimer* fTimer;
    QElapsedTimer fClock;
    G4double fSpeed;   // rad/s
    G4GenericMessenger* fMessenger;
};

#endif
#endif
