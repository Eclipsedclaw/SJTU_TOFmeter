#include "TOFmeterAppIcon.hh"

#ifdef TOF_HAS_LOGO

#include "TOFmeterLogo.hh"

#include <QApplication>
#include <QByteArray>
#include <QIcon>
#include <QImage>
#include <QPainter>
#include <QPalette>
#include <QPixmap>
#include <QSvgRenderer>

#include <algorithm>

void SetTOFmeterAppIcon()
{
  if (!qobject_cast<QApplication*>(QCoreApplication::instance())) return;

  // light marks on a dark desktop
  const bool dark = QApplication::palette().color(QPalette::Window).lightness() < 128;
  QSvgRenderer renderer(QByteArray(dark ? TOFmeterLogo::kMarkReversed : TOFmeterLogo::kMark));
  if (!renderer.isValid()) return;

  // Pixmaps from 16 to 1024 px, the mark centred with a margin as for macOS app icons
  const QSizeF mark = renderer.defaultSize();
  QIcon icon;
  for (int size : {16, 32, 64, 128, 256, 512, 1024}) {
    QImage image(size, size, QImage::Format_ARGB32_Premultiplied);
    image.fill(Qt::transparent);
    QPainter painter(&image);
    painter.setRenderHint(QPainter::Antialiasing);
    const double scale = 0.84 * size / std::max(mark.width(), mark.height());
    const QSizeF drawn = mark * scale;
    renderer.render(&painter, QRectF(0.5 * (size - drawn.width()), 0.5 * (size - drawn.height()),
                                     drawn.width(), drawn.height()));
    painter.end();
    icon.addPixmap(QPixmap::fromImage(image));
  }
  QApplication::setWindowIcon(icon);
}

#else

void SetTOFmeterAppIcon() {}

#endif
