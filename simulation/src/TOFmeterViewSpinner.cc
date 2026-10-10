#include "TOFmeterViewSpinner.hh"

#ifdef TOF_USE_QT

#include "G4GenericMessenger.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSceneHandler.hh"
#include "G4VViewer.hh"
#include "G4ViewParameters.hh"
#include "G4VisManager.hh"

#include <QCoreApplication>
#include <QEvent>
#include <QTimer>

#include <algorithm>

TOFmeterViewSpinner::TOFmeterViewSpinner()
  : fTimer(new QTimer(this)), fSpeed(0.), fMessenger(nullptr)
{
  fTimer->setInterval(40);   // 25 frames per second
  QObject::connect(fTimer, &QTimer::timeout, [this]() { Step(); });

  fMessenger = new G4GenericMessenger(this, "/tof/vis/", "Qt viewer helpers");
  fMessenger->DeclareMethod("spin", &TOFmeterViewSpinner::Spin,
                            "Turn the current viewer about its up vector at this speed (deg/s) until "
                            "the first mouse click, scroll or key press in a viewer; 0 stops")
    .SetParameterName("degPerSecond", true)
    .SetDefaultValue("10");
}

TOFmeterViewSpinner::~TOFmeterViewSpinner()
{
  if (QCoreApplication::instance()) QCoreApplication::instance()->removeEventFilter(this);
  delete fMessenger;
}

void TOFmeterViewSpinner::Spin(G4double degreesPerSecond)
{
  QCoreApplication* app = QCoreApplication::instance();
  fSpeed = degreesPerSecond * deg;
  if (fSpeed == 0. || !app) {
    fTimer->stop();
    if (app) app->removeEventFilter(this);
    return;
  }
  app->installEventFilter(this);
  fClock.start();
  fTimer->start();
}

bool TOFmeterViewSpinner::eventFilter(QObject* watched, QEvent* event)
{
  switch (event->type()) {
    case QEvent::MouseButtonPress:
    case QEvent::Wheel:
    case QEvent::KeyPress:
    case QEvent::TouchBegin:
    case QEvent::NativeGesture:
      // OpenGL areas of the viewers: QGLWidget (Qt5) or QOpenGLWidget
      if (watched->inherits("QGLWidget") || watched->inherits("QOpenGLWidget")) Spin(0.);
      break;
    default:
      break;
  }
  return false;   // the viewer still gets the event
}

void TOFmeterViewSpinner::Step()
{
  const G4double dt = std::min(fClock.restart() * 1e-3, 0.2);   // s since the last frame
  if (G4StateManager::GetStateManager()->GetCurrentState() != G4State_Idle) return;
  auto visManager = dynamic_cast<G4VisManager*>(G4VVisManager::GetConcreteInstance());
  G4VViewer* viewer = visManager ? visManager->GetCurrentViewer() : nullptr;
  if (!viewer || !viewer->GetSceneHandler() || !viewer->GetSceneHandler()->GetScene()) return;

  // as /vis/viewer/set/viewpointVector followed by /vis/viewer/refresh, without the echo
  G4ViewParameters vp = viewer->GetViewParameters();
  G4Vector3D direction = vp.GetViewpointDirection();
  direction.rotate(fSpeed * dt, vp.GetUpVector());
  vp.SetViewAndLights(direction);
  viewer->SetViewParameters(vp);
  viewer->SetView();
  viewer->ClearView();
  viewer->DrawView();
}

#endif
