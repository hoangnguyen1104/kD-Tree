#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */

kDTree::kDTree(int k)
{
    this->k = k;
    root = nullptr;
}

kDTree::~kDTree()
{
    deleteTree(root);
}

kDTreeNode *kDTree::copyTree(const kDTreeNode *node)
{
    if (node == nullptr)
    {
        return nullptr;
    }

    kDTreeNode *newNode = new kDTreeNode(node->data);
    newNode->left = copyTree(node->left);
    newNode->right = copyTree(node->right);

    return newNode;
}

void kDTree::deleteTree(kDTreeNode *node)
{
    if (node == nullptr)
    {
        return;
    }

    deleteTree(node->left);
    deleteTree(node->right);
    delete node;
}

const kDTree &kDTree::operator=(const kDTree &other)
{
    if (this != &other)
    {
        k = other.k;
        deleteTree(root);
        root = copyTree(other.root);
    }

    return *this;
}

kDTree::kDTree(const kDTree &other)
{
    k = other.k;
    root = copyTree(other.root);
}

void kDTree::inorderHelper(const kDTreeNode *node) const
{
    if (node == nullptr)
    {
        return;
    }
    inorderHelper(node->left);
    std::cout << *node << " ";
    inorderHelper(node->right);
}

void kDTree::inorderTraversal() const
{
    inorderHelper(root);
}

void kDTree::preorderHelper(const kDTreeNode *node) const
{
    if (node == nullptr)
    {
        return;
    }

    std::cout << *node << " ";

    preorderHelper(node->left);
    preorderHelper(node->right);
}

void kDTree::preorderTraversal() const
{
    preorderHelper(root);
}

void kDTree::postorderHelper(const kDTreeNode *node) const
{
    if (node == nullptr)
    {
        return;
    }

    postorderHelper(node->left);
    postorderHelper(node->right);
    std::cout << *node << " ";
}

void kDTree::postorderTraversal() const
{
    postorderHelper(root);
}

int kDTree::countNodes(const kDTreeNode *node) const
{
    if (node == nullptr)
    {
        return 0;
    }

    int leftCount = countNodes(node->left);
    int rightCount = countNodes(node->right);

    return leftCount + rightCount + 1;
}
int kDTree::countLeafNodes(const kDTreeNode *node) const
{
    if (node == nullptr)
    {
        return 0;
    }

    if (node->left == nullptr && node->right == nullptr)
    {
        return 1;
    }

    int leftCount = countLeafNodes(node->left);
    int rightCount = countLeafNodes(node->right);

    return leftCount + rightCount;
}
int kDTree::calculateHeight(const kDTreeNode *node) const
{
    if (node == nullptr)
    {
        return 0;
    }

    int leftHeight = calculateHeight(node->left);
    int rightHeight = calculateHeight(node->right);

    return std::max(leftHeight, rightHeight) + 1;
}

int kDTree::height() const
{
    return calculateHeight(root);
}

int kDTree::nodeCount() const
{
    return countNodes(root);
}

int kDTree::leafCount() const
{
    return countLeafNodes(root);
}

kDTreeNode *kDTree::insertHelper(kDTreeNode *node, const std::vector<int> &point, int depth)
{
    if (node == nullptr)
    {
        return new kDTreeNode(point);
    }

    int currentDimension = depth % k;

    if (point[currentDimension] < node->data[currentDimension])
    {
        node->left = insertHelper(node->left, point, depth + 1);
    }
    else
    {
        node->right = insertHelper(node->right, point, depth + 1);
    }

    return node;
}

void kDTree::insert(const std::vector<int> &point)
{
    root = insertHelper(root, point, 0);
}

kDTreeNode *kDTree::removeHelper(kDTreeNode *node, const std::vector<int> &point, int depth)
{
    if (node == nullptr)
    {
        return nullptr;
    }
    int currentDimension = depth % k;

    if (point == node->data)
    {
        if (node->right != nullptr)
        {
            kDTreeNode *minNode = findMinNode(node->right, currentDimension);
            node->data = minNode->data;
            node->right = removeHelper(node->right, minNode->data, depth + 1);
        }
        else if (node->left != nullptr)
        {
            kDTreeNode *minNode = findMinNode(node->left, currentDimension);
            node->data = minNode->data;
            node->right = removeHelper(node->left, minNode->data, depth + 1);
            node->left = nullptr;
        }
        else
        {
            delete node;
            return nullptr;
        }
    }
    else if (point[currentDimension] < node->data[currentDimension])
    {
        node->left = removeHelper(node->left, point, depth + 1);
    }
    else
    {
        node->right = removeHelper(node->right, point, depth + 1);
    }

    return node;
}

kDTreeNode *kDTree::findMinNode(kDTreeNode *node, int dimension) const
{
    if (node == nullptr)
    {
        return nullptr;
    }

    kDTreeNode *minNode = node;
    kDTreeNode *leftMin = findMinNode(node->left, dimension);
    kDTreeNode *rightMin = findMinNode(node->right, dimension);

    if (leftMin != nullptr && leftMin->data[dimension] < minNode->data[dimension])
    {
        minNode = leftMin;
    }

    if (rightMin != nullptr && rightMin->data[dimension] < minNode->data[dimension])
    {
        minNode = rightMin;
    }

    return minNode;
}

void kDTree::remove(const std::vector<int> &point)
{
    root = removeHelper(root, point, 0);
}

bool kDTree::searchHelper(kDTreeNode *node, const std::vector<int> &point, int depth)
{
    if (node == nullptr)
    {
        return false;
    }

    int currentDimension = depth % k;

    if (point == node->data)
    {
        return true;
    }

    if (point[currentDimension] < node->data[currentDimension])
    {
        return searchHelper(node->left, point, depth + 1);
    }
    else
    {
        return searchHelper(node->right, point, depth + 1);
    }
}

bool kDTree::search(const std::vector<int> &point)
{
    return searchHelper(root, point, 0);
}

void kDTree::bubbleSort(std::vector<std::vector<int>> &pointList, int dimension)
{
    int n = pointList.size();
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - i - 1; j++)
        {
            if (pointList[j][dimension] > pointList[j + 1][dimension])
            {
                std::swap(pointList[j], pointList[j + 1]);
            }
        }
    }
}

kDTreeNode *kDTree::buildTreeHelper(std::vector<std::vector<int>> pointList, int depth)
{
    if (pointList.size() == 0)
    {
        return nullptr;
    }
    if (pointList.size() == 1)
    {
        kDTreeNode *node = new kDTreeNode(pointList[0]);
        return node;
    }

    int currentDimension = depth % k;
    int n = pointList.size();
    // for (int i=0; i<n; i++){
    //     kDTreeNode *node2 = new kDTreeNode(pointList[i]);
    //     cout << node2 << " ";
    // }
    // cout<<"\n";
    vector<vector<int>> point_left, point_right;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - i - 1; j++)
        {
            if (pointList[j][currentDimension] > pointList[j + 1][currentDimension])
            {
                std::swap(pointList[j], pointList[j + 1]);
            }
        }
    }

    // for (int i=0; i<n; i++){
    //     kDTreeNode *node2 = new kDTreeNode(pointList[i]);
    //     cout<<*node2<<" ";
    // }
    // cout<<"\n";

    int mid = (n + 1) / 2 - 1;
    // cout<<mid<<"\n";
    for (int i = 0; i < mid; i++)
        point_left.push_back(pointList[i]);
    for (int i = mid + 1; i < n; i++)
        point_right.push_back(pointList[i]);
    kDTreeNode *node = new kDTreeNode(pointList[mid]);

    node->left = buildTreeHelper(point_left, depth + 1);
    node->right = buildTreeHelper(point_right, depth + 1);

    return node;
}

void kDTree::buildTree(const std::vector<std::vector<int>> &pointList)
{
    std::vector<std::vector<int>> tmp = pointList;
    root = buildTreeHelper(tmp, 0);
}

double kDTree::calculateDistance(const std::vector<int> &point1, const std::vector<int> &point2)
{
    double sum = 0.0;
    for (size_t i = 0; i < point1.size(); ++i)
    {
        double diff = static_cast<double>(point1[i] - point2[i]);
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

void kDTree::nearestNeighbourHelper(kDTreeNode *node, const std::vector<int> &target, kDTreeNode *&best, double &bestDistance, int depth)
{
    if (node == nullptr)
        return;

    double currentDistance = calculateDistance(node->data, target);
    // cout<<*node<<"\n";

    if (currentDistance < bestDistance)
    {
        best = node;        
        bestDistance = currentDistance;
        // cout<<bestDistance<<"oo\n";
    }

    int currentDimension = depth % k;
    int targetValue = target[currentDimension];
    int nodeValue = node->data[currentDimension];
    // cout<<targetValue<<" "<<nodeValue<<"\n";

    if (targetValue < nodeValue)
    {
        nearestNeighbourHelper(node->left, target, best, bestDistance, depth + 1);
        if (nodeValue - targetValue <= bestDistance)
            nearestNeighbourHelper(node->right, target, best, bestDistance, depth + 1);
    }
    else
    {
        nearestNeighbourHelper(node->right, target, best, bestDistance, depth + 1);
        if (targetValue - nodeValue <= bestDistance)
            nearestNeighbourHelper(node->left, target, best, bestDistance, depth + 1);
    }
}

void kDTree::nearestNeighbour(const std::vector<int> &target, kDTreeNode *&best)
{
    double bestDistance = 999999;
    nearestNeighbourHelper(root, target, best, bestDistance, 0);
}

void kDTree::kNearestNeighbourHelper(kDTreeNode *node, const std::vector<int> &target, int k,
                                     std::vector<std::pair<double, kDTreeNode *>> &nearestList,
                                     double &maxDistance, int depth)
{
    // for (int i=0; i<nearestList.size(); i++) cout<<nearestList[i].first<<" "<<*nearestList[i].second<<"\n";    
    
    // cout<<nearestList.size()<<"\n";
    if (node == nullptr)
        return;

    double currentDistance = calculateDistance(node->data, target);        
    if (node->left == nullptr && node->right == nullptr){
        //cout<<*node<<" "<<"===\n";
        if (currentDistance <= maxDistance || nearestList.size() < k){
            nearestList.push_back({currentDistance, node});        
            int insertPos = nearestList.size() - 1;
            while (insertPos > 0 && nearestList[insertPos - 1].first > currentDistance)
            {
                nearestList[insertPos] = nearestList[insertPos - 1];
                insertPos--;
            }
            nearestList[insertPos] = {currentDistance, node};
            if (nearestList.size() > k) nearestList.pop_back();
            maxDistance = nearestList.back().first;
        }
        return;
    }

    int currentDimension = depth % k;
    int targetValue = target[currentDimension];
    int nodeValue = node->data[currentDimension];

    if (targetValue < nodeValue)
    {
        kNearestNeighbourHelper(node->left, target, k, nearestList, maxDistance, depth + 1);
        // cout<<maxDistance<<" ";
        if (nearestList.size() > 0) maxDistance = nearestList.back().first;
        // cout<<maxDistance<<"\n";
        if (currentDistance <= maxDistance || nearestList.size() < k){
            nearestList.push_back({currentDistance, node});        
            int insertPos = nearestList.size() - 1;
            while (insertPos > 0 && nearestList[insertPos - 1].first > currentDistance)
            {
                nearestList[insertPos] = nearestList[insertPos - 1];
                insertPos--;
            }
            nearestList[insertPos] = {currentDistance, node};
            if (nearestList.size() > k) nearestList.pop_back();
            maxDistance = nearestList.back().first;
        }
        if (nodeValue - targetValue <= maxDistance)
            kNearestNeighbourHelper(node->right, target, k, nearestList, maxDistance, depth + 1);
    }
    else
    {
        kNearestNeighbourHelper(node->right, target, k, nearestList, maxDistance, depth + 1);
        // cout<<maxDistance<<" ";
        if (nearestList.size() > 0) maxDistance = nearestList.back().first;
        // cout<<maxDistance<<"\n";
        if (currentDistance <= maxDistance || nearestList.size() < k){
            nearestList.push_back({currentDistance, node});        
            int insertPos = nearestList.size() - 1;
            while (insertPos > 0 && nearestList[insertPos - 1].first > currentDistance)
            {
                nearestList[insertPos] = nearestList[insertPos - 1];
                insertPos--;
            }
            nearestList[insertPos] = {currentDistance, node};
            if (nearestList.size() > k) nearestList.pop_back();
            maxDistance = nearestList.back().first;
        }
        if (targetValue - nodeValue <= maxDistance)
            kNearestNeighbourHelper(node->left, target, k, nearestList, maxDistance, depth + 1);
    }
}

void kDTree::kNearestNeighbour(const std::vector<int> &target, int k, std::vector<kDTreeNode *> &bestList)
{
    std::vector<std::pair<double, kDTreeNode *>> nearestList;

    double maxDistance = 999999;

    kNearestNeighbourHelper(root, target, k, nearestList, maxDistance, 0);
    // for (int i=0; i<nearestList.size(); i++) cout<<nearestList[i].first<<" "<<*nearestList[i].second<<"\n";    

    bestList.reserve(k);

    for (size_t i = 0; i < k && i < nearestList.size(); i++)
    {
        kDTreeNode *tmp = new kDTreeNode(nearestList[i].second->data);                
        bestList.push_back(tmp);
    }    
}

kNN::kNN(int k)
{
    this->k = k;
}

void kNN::fit(Dataset &X_train, Dataset &y_train)
{
    this->X_train = X_train;
    this->y_train = y_train;
    list<list<int>> data_y = y_train.data;
    std::vector<int> labels;
    std::vector<std::vector<int>> pointList;
    for (const std::list<int> &innerList : data_y)
    {
        for (int element : innerList)
        {
            labels.push_back(element);
        }
    }
    list<list<int>> data_x = X_train.data;
    int row = 0;
    for (const std::list<int> &innerList : data_x)
    {
        std::vector<int> features;
        features.push_back(labels[row]);
        for (int element : innerList)
        {
            features.push_back(element);
        }
        pointList.push_back(features);
        row += 1;
    }
    tree.buildTree(pointList);
}

Dataset kNN::predict(Dataset &X_test)
{
    Dataset y_pred;
    list<list<int>> data_x_test = X_test.data;
    // int nRows = data_x_test.size();
    // int nCols = 1;
    // y_pred.getShape(nRows, nCols);
    list<list<int>> data;
    for (const std::list<int> &innerList : data_x_test)
    {
        std::vector<int> features;
        vector<kDTreeNode *> bestList;
        for (int element : innerList)
        {
            features.push_back(element);
        }
        tree.kNearestNeighbour(features, k, bestList);
        std::vector<int> labels;
        for (kDTreeNode *node : bestList)
        {
            labels.push_back(node->data[0]);
        }
        // for (int i=0; i<labels.size(); i++) std::cout<<labels[i]<<" ";
        // std::cout<<std::endl;
        int max_d = 0, d = -1, dem;
        for (int i = 0; i < labels.size(); i++)
        {
            dem = 0;
            for (int j = 0; j < labels.size(); j++)
                if (labels[j] == labels[i])
                    dem += 1;
            if (max_d < dem)
            {
                max_d = dem;
                d = labels[i];
            }
        }
        // std::cout<<d;
        // std::cout<<std::endl;
        list<int> cols;
        cols.push_back(d);
        data.push_back(cols);
    }
    y_pred.data = data;
    y_pred.columnName.push_back("label");
    return y_pred;
}

double kNN::score(const Dataset &y_test, const Dataset &y_pred)
{
    vector<int> labels_test, labels_pred;
    list<list<int>> data = y_test.data;
    for (const std::list<int> &innerList : data)
    {
        for (int element : innerList)
        {
            labels_test.push_back(element);
        }
    }
    data = y_pred.data;
    for (const std::list<int> &innerList : data)
    {
        for (int element : innerList)
        {
            labels_pred.push_back(element);
        }
    }
    int numCorrect = 0;
    int nRows = labels_pred.size();
    for (int i = 0; i < nRows; i++){
        if (labels_pred[i] == labels_test[i])
            numCorrect += 1;
    }
        
    return static_cast<double>(numCorrect) / nRows;
}