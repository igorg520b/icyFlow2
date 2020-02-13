#include "czinsertiontool.h"


icy::ExtendedFace::ExtendedFace(std::tuple<int,int,int> ftuple) {
    vrts[0] = std::get<0>(ftuple);
    vrts[1] = std::get<1>(ftuple);
    vrts[2] = std::get<2>(ftuple);
    this->ftuple = ftuple;
}


const int icy::ExtendedElement::myConvention[4][3] = {
{3,1,2},
{0,3,2},
{0,1,3},
{1,0,2}};

std::tuple<int,int,int> icy::ExtendedElement::make_sorted_tuple(int a, int b, int c) {
    int arr[3]={a,b,c};
    std::sort(arr, arr+3);
    return std::make_tuple(arr[0],arr[1],arr[2]);
}

void icy::ExtendedElement::GenerateFaces()
{
    for(int i=0;i<4;i++) {
        faces_as_tuples[i] = make_sorted_tuple(
                    vrts[myConvention[i][0]],
                vrts[myConvention[i][1]],
                vrts[myConvention[i][2]]);
    }
}

void icy::ExtendedElement::SubstituteNode(int oldNode, int newNode) {
    bool found = false;
    for(int i=0;i<4;i++)
        if(vrts[i]==oldNode) {
            found = true;
            vrts[i] = newNode;
            break;
        }
    if(!found) throw std::runtime_error("substituted node not found");
}

int icy::ExtendedElement::WhichFace(ExtendedFace &fc) {
    if(fc.ftuple == make_sorted_tuple(vrts[1],vrts[2],vrts[3])) return 0;
    else if(fc.ftuple == make_sorted_tuple(vrts[0],vrts[2],vrts[3])) return 1;
    else if(fc.ftuple == make_sorted_tuple(vrts[0],vrts[1],vrts[3])) return 2;
    else if(fc.ftuple == make_sorted_tuple(vrts[0],vrts[1],vrts[2])) return 3;
    else throw std::runtime_error("whichface error");
}

bool icy::ExtendedElement::ContainsFace(ExtendedFace &fc)
{
    for(int i=0;i<3;i++) {
        int v = fc.vrts[i];
        if(v != vrts[0] && v != vrts[1] && v != vrts[2] && v != vrts[3]) return false;
    }
    return true;
}

void icy::ExtendedElement::AddFaceIfContains(ExtendedFace &fc)
{
    if(ContainsFace(fc)) faces.push_back(fc.id);
}



icy::ExtendedFace icy::ExtendedElement::CreateFaceFromIdx(int idx)
{
    ExtendedFace r;

    r.vrts[0] = vrts[myConvention[idx][0]];
    r.vrts[1] = vrts[myConvention[idx][1]];
    r.vrts[2] = vrts[myConvention[idx][2]];

    return r;
}



//============================================================




void icy::CZInsertionTool::Extend(Mesh &mg)
{
    nodes.resize(mg.nodes.size());

    for(size_t i= 0;i<mg.nodes.size();i++) {
        Node &nd = mg.nodes[i];
        ExtendedNode &nd2 = nodes[i];
        nd2.id = i;
        if(nd.id != (int)i) throw std::runtime_error("node index error");
        nd2.x0 = nd.x0;
        nd2.y0 = nd.y0;
        nd2.z0 = nd.z0;
    }


    elems.resize(mg.elems.size());
    for(size_t i=0;i<mg.elems.size();i++) {
        Element &elem = mg.elems[i];
        ExtendedElement &elem2 = elems[i];
        elem2.tag = elem.tag;
        elem2.id = i;
        for(int j=0;j<4;j++) elem2.vrts[j] = elem.vrts[j]->id;
        elem2.GenerateFaces();
        for(int j=0;j<4;j++)
        {
            std::tuple<int,int,int> ftuple = elem2.faces_as_tuples[j];
            if(face_map.find(ftuple) == face_map.end()) {
                // add ExtendedFace
                face_map[ftuple] = ExtendedFace(ftuple);
            }
            face_map[ftuple].connectedElements.push_back(elem2.id);
            if(face_map[ftuple].connectedElements.size() > 2)
                throw std::runtime_error("too many elements");
        }
    }

    int count = 0;
    for(auto &val : face_map) {
        ExtendedFace &fc = val.second;
        if(fc.connectedElements.size() == 1 ||
                (fc.connectedElements.size() ==2 &&
                 elems[fc.connectedElements[0]].tag != elems[fc.connectedElements[1]].tag)) {
            fc.id = count++;
            faces.push_back(fc);
        }
        else if(fc.connectedElements.size() == 0 || fc.connectedElements.size() > 2)
            throw std::runtime_error("face with incorrect number of adj tetra");

    }

    czs.clear();

}


void icy::CZInsertionTool::ConvertBack(Mesh &mg)
{

}



void icy::CZInsertionTool::InsertCohesiveElements(Mesh &mg)
{
    Extend(mg);
//    IdentifyParentsOfTriangles();

    // store surface faces and inner faces
    std::vector<int> surface;
    std::vector<int> innerTris;
    for(auto &fc : faces){
        if(fc.connectedElements.size() == 1) surface.push_back(fc.id);
        else if(fc.connectedElements.size() == 2) innerTris.push_back(fc.id);
        else throw std::runtime_error("invalid element count");
    }
    std::cout << "surface fcs " << surface.size() << std::endl;
    std::cout << "innter fcs " << innerTris.size() << std::endl;

    for(auto &fcidx : surface) {
        ExtendedFace &fc = faces[fcidx];
        for(int i=0;i<3;i++) nodes[fc.vrts[i]].isSurface = true;
    }

    std::vector<std::pair<int, int>> surfaceFaces; // <elemIdx, faceIdx>

    for(auto &fc : faces) {
        for(auto &elemIdx : fc.connectedElements) {
            ExtendedElement &elem = elems[elemIdx];
            int which = elem.WhichFace(fc);
            if(fc.connectedElements.size()==1) surfaceFaces.push_back(std::make_pair(elemIdx,which));
            elem.faces.push_back(which);
        }
    }

    // create czs
    for(auto &fcidx : innerTris) {
        ExtendedFace &fc = faces[fcidx];
        ExtendedCZ cz;
        cz.e0 = fc.connectedElements[0];
        cz.e1 = fc.connectedElements[1];
        ExtendedElement &elem0 = elems[cz.e0];
        cz.fidx0 = elem0.WhichFace(fc);
        ExtendedElement &elem1 = elems[cz.e1];
        cz.fidx1 = elem1.WhichFace(fc);

        int first_nd_idx = elem0.vrts[icy::ExtendedElement::myConvention[cz.fidx0][0]];
        cz.relativeOrientationIndex = -1;
        for(int i=0;i<3;i++) {
            int second_nd_idx = elem1.vrts[icy::ExtendedElement::myConvention[cz.fidx1][i]];
            if(first_nd_idx == second_nd_idx) {cz.relativeOrientationIndex=i; break;}
        }
        if(cz.relativeOrientationIndex == -1) throw std::runtime_error("relative orientation idx");
    }


    // prepare linked list (all nodes)


    // list the nodes that are connected to CZs
    for(auto &fidx : innerTris) {
        ExtendedFace &fc = faces[fidx];
        for(int i=0;i<3;i++) {
            ExtendedNode &nd = nodes[fc.vrts[i]];
            nd.belongs_to_cz = true;
        }
    }

    std::vector<int> nczs; // nodes connected to cohesive zones
    for(auto &nd : nodes)
        if(nd.belongs_to_cz) nczs.push_back(nd.id);


    // split the nodes, which belong to cohesive elements
//    foreach (ExtendedNode nd in nczs) nd.SplitNode(ll);


    // linked list becomes the new list of nodes; resequence


}
