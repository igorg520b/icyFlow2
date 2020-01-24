#ifndef MESH_H
#define MESH_H

#include <QObject>
#include <vector>
#include "node.h"
#include "element.h"
#include "cz.h"
#include "face.h"
#include "surfacefragment.h"
#include "translation.h"

namespace icy { class Mesh; }

class icy::Mesh : public QObject
{
    Q_OBJECT

public:
    Mesh();
    std::vector<Node*> nodes;
    std::vector<Element*> elems;
    std::vector<CZ*> czs;
    std::vector<Face*> faces;
    std::vector<SurfaceFragment*> surfaceFragments;
    std::vector<Element*> surfaceElements; // elements that can potentially come in contact

      //  public TranslationCollection translationCollection = new TranslationCollection();


};

#endif // MESH_H
