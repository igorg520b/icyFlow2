#include "czinsertiontool.h"


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
    const int face0 = (8|2|4);
    const int face1 = (1|8|4);
    const int face2 = (1|2|8);
    const int face3 = (1|2|4);

    auto findIdx = [](int ndId, int (&vrts)[4]) {
        for(int i=0;i<4;i++) if(vrts[i]==ndId) return i;
        throw std::runtime_error("idx not found"); };

    int fid = 0;
    for(int i=0;i<3;i++) fid |= findIdx(fc.vrts[i], vrts);

    if(fid==face0) return 0;
    else if(fid==face1) return 1;
    else if(fid==face2) return 2;
    else if(fid==face3) return 3;
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

const int icy::ExtendedElement::myConvention[4][3] = {
{3,1,2},
{0,3,2},
{0,1,3},
{1,0,2}};

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
    faces.resize(mg.faces.size());
    elems.resize(mg.elems.size());

    for(size_t i= 0;i<mg.nodes.size();i++) {
        Node &nd = mg.nodes[i];
        ExtendedNode &nd2 = nodes[i];
        nd2.id = i;
        if(nd.id != (int)i) throw std::runtime_error("node index error");
        nd2.x0 = nd.x0;
        nd2.y0 = nd.y0;
        nd2.z0 = nd.z0;
    }

    for(size_t i=0;i<mg.faces.size();i++) {
        Face &fc = mg.faces[i];
        ExtendedFace &fc2 = faces[i];
        fc2.id = i;
        for(int j=0;j<3;j++) fc2.vrts[j] = fc.vrts[j]->id;
    }

    for(size_t i=0;i<mg.elems.size();i++) {
        Element &elem = mg.elems[i];
        ExtendedElement &elem2 = elems[i];
        elem2.tag = elem.tag;
        elem2.id = i;
        for(int j=0;j<4;j++) elem2.vrts[j] = elem.vrts[j]->id;
    }

    czs.clear();
}


void icy::CZInsertionTool::ConvertBack(Mesh &mg)
{

}


void icy::CZInsertionTool::IdentifyParentsOfTriangles()
{
//    for(auto &nd : nodes) nd.elementsOfNode.clear();
//    for(auto &fc : faces) fc.elementsOfFace.clear();

    for(auto &elem : elems)
        for(int i=0;i<4;i++){
            int ndidx = elem.vrts[i];
            ExtendedNode &nd = nodes[ndidx];
            nd.elementsOfNode.push_back(elem.id);
            nd.grains.insert(elem.tag);
        }

    for(auto &fc : faces) {
        ExtendedNode &nd = nodes[fc.vrts[0]];
        for(auto &elemidx : nd.elementsOfNode) {
            ExtendedElement &elem = elems[elemidx];
            if(elem.ContainsFace(fc)) {
                fc.elementsOfFace.push_back(elemidx);
                if(fc.elementsOfFace.size()==2) break;
            }
        }
    }
}


void icy::CZInsertionTool::InsertCohesiveElements(Mesh &mg)
{
    Extend(mg);
    IdentifyParentsOfTriangles();

    // store surface faces and inner faces
    std::vector<int> surface;
    std::vector<int> innerTris;
    for(auto &fc : faces){
        if(fc.elementsOfFace.size() == 1) surface.push_back(fc.id);
        else if(fc.elementsOfFace.size() == 2) innerTris.push_back(fc.id);
        else throw std::runtime_error("invalid element count");
    }
    std::cout << "surface fcs " << surface.size() << std::endl;
    std::cout << "innter fcs " << innerTris.size() << std::endl;
}
