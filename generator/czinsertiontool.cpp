#include "czinsertiontool.h"

// ====================================== node
void icy::ExtendedNode::SplitNode(
        std::list<int> &allNodesAsLL,
        std::vector<ExtendedNode> &allNodes,
        std::vector<ExtendedElement> &elems,
        std::vector<ExtendedNode> &additionalNodes)
{
    std::list<int>::iterator itLL = std::find(allNodesAsLL.begin(), allNodesAsLL.end(), this->id);
    if(itLL == allNodesAsLL.end()) throw std::runtime_error("node not found in LL");

    std::vector<int> grains_as_vector(grains.begin(), grains.end());

//    std::map<int,bool> localNewNodes; // for double-check
//    localNewNodes[this->id] = false;

    // add (nGranules-1) new nodes and update connected elements
    int nGrains = (int)grains.size();
//    std::cout << "split nd " << id << "; grains " << nGrains << "; vec size " << grains_as_vector.size()
//    std::cout << "; elems of node " << elementsOfNode.size() << std::endl;
    for(int i=0;i<nGrains-1;i++) {
        int currentGrain = grains_as_vector[i];
//        std::cout << "current grain " << currentGrain << std::endl;
        ExtendedNode newNode = *this;
        newNode.id = (int)allNodes.size()+additionalNodes.size();
        additionalNodes.push_back(newNode);
        allNodesAsLL.insert(itLL, newNode.id);
//        localNewNodes[newNode.id]=false;

        for(auto &elemIdx : elementsOfNode) {
            ExtendedElement &elem = elems[elemIdx];
            if(elem.tag == currentGrain)
                elem.SubstituteNode(this->id, newNode.id);
        }
    }
/*
    // verify that all localNewNodes are connected
    for(auto &elemIdx : elementsOfNode) {
        for(int i=0;i<4;i++) {
            int ndIdx = elems[elemIdx].vrts[i];
            if(localNewNodes.find(ndIdx)!=localNewNodes.end())
                localNewNodes[ndIdx] = true;
        }
    }

    bool disconnect_found = false;
    for(auto &pair : localNewNodes)
        if(pair.second == false)
        {
            disconnect_found = true;
            std::cout << "disconnected node " << pair.first << " is present when splitting " << this->id <<std::endl;
        }

    if(disconnect_found) {
        std::cout << "  nGrains " << nGrains << "; localNewNodes: " << localNewNodes.size();
        std::cout << "; elems of nds: " << elementsOfNode.size() << std::endl;
    }
    */
}


// ================================== face
icy::ExtendedFace::ExtendedFace(std::tuple<int,int,int> ftuple) {
    vrts[0] = std::get<0>(ftuple);
    vrts[1] = std::get<1>(ftuple);
    vrts[2] = std::get<2>(ftuple);
    this->ftuple = ftuple;
}

//==================================== element
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
//    std::cout << "elem " << id << "; subst " << oldNode << " for " << newNode << std::endl;
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

//=========================================== CZ

void icy::ExtendedCZ::ReinitializeVerticeArrays(std::vector<ExtendedElement> &elems)
{
    for(int i=0;i<3;i++) {
        ExtendedElement &elem0 = elems[e0];
        ExtendedElement &elem1 = elems[e1];

        vrts[i] = elem0.vrts[icy::ExtendedElement::myConvention[fidx0][i]];
        vrts[i+3] = elem1.vrts[icy::ExtendedElement::myConvention[fidx1][(i + relativeOrientationIndex) % 3]];
    }
    int temp = vrts[5];
    vrts[5] = vrts[4];
    vrts[4] = temp;
}



//============================================================




void icy::CZInsertionTool::Extend(Mesh &mg)
{
    // count disconnected nds in mg
    std::vector<bool> double_check;
    double_check.resize(mg.nodes.size());
    for(size_t i=0;i<double_check.size();i++) double_check[i] = false;

    for(auto &elem : mg.elems)
        for(int i=0;i<4;i++)
            double_check[elem.vrts[i]->id] = true;

    int count = 0;
    for(size_t i=0;i<double_check.size();i++) {
        if(!double_check[i]) count++;
    }
    std::cout << "count disconnected nds in mg: "<< count << std::endl;

    czs.clear();
    faces.clear();
    elems.clear();
    face_map.clear();

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
        for(int j=0;j<4;j++) {
            int nodeIdx = elem.vrts[j]->id;
            ExtendedNode &nd = nodes[nodeIdx];
            nd.grains.insert(elem.tag); // populate grains of the node
            nd.elementsOfNode.push_back(elem2.id);
            elem2.vrts[j] = nodeIdx;
        }
        elem2.GenerateFaces();

        // create a "map" of faces; each face connecting to 1-2 elements
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

    // only pick faces that are either (1) exterior or (2) separate grains
    count = 0;
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
}



void icy::CZInsertionTool::InsertCohesiveElements(Mesh &mg)
{
    Extend(mg);

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

//    for(auto &fcidx : surface) {
//        ExtendedFace &fc = faces[fcidx];
//        for(int i=0;i<3;i++) nodes[fc.vrts[i]].isSurface = true;
//    }

//    std::vector<std::pair<int, int>> surfaceFaces; // <elemIdx, faceIdx>

    for(auto &fc : faces) {
        for(auto &elemIdx : fc.connectedElements) {
            ExtendedElement &elem = elems[elemIdx];
            int which = elem.WhichFace(fc);
            if(fc.connectedElements.size()==1) //surfaceFaces.push_back(std::make_pair(elemIdx,which));
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
        czs.push_back(cz);
    }

    // list the nodes that are connected to CZs
    for(auto &fidx : innerTris) {
        ExtendedFace &fc = faces[fidx];
        for(int i=0;i<3;i++) {
            ExtendedNode &nd = nodes[fc.vrts[i]];
            nd.belongs_to_cz = true;
        }
    }

    // prepare linked list and a list of nodes that need to split
    std::vector<int> nczs; // nodes connected to cohesive zones
    allNodesAsLL.clear();
    for(auto &nd : nodes) {
        allNodesAsLL.push_back(nd.id);
        if(nd.belongs_to_cz) nczs.push_back(nd.id);
    }

    // split the nodes that belong to cohesive elements
    std::cout << "nodes before " << nodes.size() << std::endl;

    std::vector<ExtendedNode> additionalNodes;
    for(auto &ndIdx : nczs) {
        ExtendedNode &nd = nodes[ndIdx];
        nd.SplitNode(allNodesAsLL, nodes, elems, additionalNodes);
    }

    std::cout << "additional nodes size " << additionalNodes.size() << std::endl;
    // additionalNodes -> nodes
    nodes.insert(nodes.end(), additionalNodes.begin(),additionalNodes.end());

    std::cout << "nodes after split " << nodes.size() << std::endl;
    std::cout << "idxs in allNodesAsLL " << allNodesAsLL.size() << std::endl;
    std::cout << "number of czs " << czs.size() << std::endl;

    // linked list becomes the new list of nodes
    // create a new "vector" of nodes from allNodesAsLL and update elem.vrts
    std::map<int, int> resequence;
    std::vector<ExtendedNode> nodes2;
    int count = 0;
    for(auto &ndIdx : allNodesAsLL) {
        ExtendedNode &nd = nodes[ndIdx];
        nd.newId = count++;
        nodes2.push_back(nd);
    }

    // replace node indices in elems
    for(auto &elem : elems)
        for(int i=0;i<4;i++)
            elem.vrts[i] = nodes[elem.vrts[i]].newId;

    // replace the node array
    nodes = nodes2;
    for(auto &nd : nodes) nd.id = nd.newId;


    // infer cz.vrts[] from fidx
    for(auto &cz : czs) cz.ReinitializeVerticeArrays(elems);

    // FACES
    faces.clear();

    // surface faces are re-created from element.faces
    for(auto &elem : elems)
        for(int &fcIdx : elem.faces)
            faces.push_back(elem.CreateFaceFromIdx(fcIdx));
    std::cout << "exposed faces " << faces.size() << std::endl;

    // cz faces are re-created from czs
    count = 0;
    for(auto &cz : czs) {
        ExtendedElement &elem0 = elems[cz.e0];
        ExtendedElement &elem1 = elems[cz.e1];

        ExtendedFace fc0 = elem0.CreateFaceFromIdx(cz.fidx0);
        ExtendedFace fc1 = elem1.CreateFaceFromIdx(cz.fidx1);
        cz.fid0 = fc0.id = (int)faces.size();
        faces.push_back(fc0);
        cz.fid1 = fc1.id = (int)faces.size();
        faces.push_back(fc1);
        fc0.exposed = false;
        fc1.exposed = false;
    }
    std::cout << "total faces " << faces.size() << std::endl;

    // CONVERT BACK
    mg.Clear();

    // convert nodes
    mg.nodes.resize(nodes.size());
    for(int i=0;i<(int)nodes.size();i++) {
        Node &nd = mg.nodes[i];
        ExtendedNode &nde = nodes[i];
        nd.x0 = nd.cx = nd.tx = nde.x0;
        nd.y0 = nd.cy = nd.ty = nde.y0;
        nd.z0 = nd.cz = nd.tz = nde.z0;
        nd.id = nde.id;
    }

    // elements
    mg.elems.resize(elems.size());
    for(int i=0;i<(int)elems.size();i++) {
        Element &elem = mg.elems[i];
        ExtendedElement &elemx = elems[i];
        elem.id = i;
        elem.tag = elemx.tag;
        for(int j=0;j<4;j++) elem.vrts[j] = &mg.nodes[elemx.vrts[j]];
    }

    // faces
    mg.faces.resize(faces.size());
    for(int i=0;i<(int)faces.size();i++) {
        Face &fc = mg.faces[i];
        ExtendedFace &fcx = faces[i];
        fc.id = i;
        fc.exposed = fcx.exposed;

        for(int j=0;j<3;j++) fc.vrts[j] = &mg.nodes[fcx.vrts[j]];
    }

    // czs
    mg.czs.resize(czs.size());
    for(int i=0;i<(int)czs.size();i++) {
        ExtendedCZ &czx = czs[i];
        CZ &cz = mg.czs[i];
        for(int j=0;j<6;j++) cz.vrts[j] = &mg.nodes[czx.vrts[j]];
        cz.faces[0] = &mg.faces[czx.fid0];
        cz.faces[1] = &mg.faces[czx.fid1];
    }

    // verify that each node is connected to either element or cz

    std::vector<bool> double_check;
    double_check.resize(mg.nodes.size());
    for(size_t i=0;i<double_check.size();i++) double_check[i] = false;

    for(auto &elem : mg.elems)
        for(int i=0;i<4;i++)
            double_check[elem.vrts[i]->id] = true;

    count = 0;
    for(size_t i=0;i<double_check.size();i++) {
        if(!double_check[i]) count++;
    }
    std::cout << "count disconnected nds: "<< count << std::endl;

}


