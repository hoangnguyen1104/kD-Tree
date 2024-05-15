#include "main.hpp"
#include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */
struct kDTreeNode
{
    vector<int> data;
    kDTreeNode *left;
    kDTreeNode *right;
    int label;
    kDTreeNode(vector<int> data, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
    }

    friend ostream &operator<<(ostream &os, const kDTreeNode &node)
    {
        os << "(";
        for (int i = 0; i < node.data.size(); i++)
        {
            os << node.data[i];
            if (i != node.data.size() - 1)
            {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }
};

class kDTree
{
private:
    int k;
    kDTreeNode *root;
    kDTreeNode *copyTree(const kDTreeNode *node);
    void deleteTree(kDTreeNode *node);
    void inorderHelper(const kDTreeNode *node) const;
    void preorderHelper(const kDTreeNode *node) const;
    void postorderHelper(const kDTreeNode *node) const;
    int countNodes(const kDTreeNode *node) const;
    int countLeafNodes(const kDTreeNode *node) const;
    int calculateHeight(const kDTreeNode *node) const;
    kDTreeNode *insertHelper(kDTreeNode *node, const std::vector<int> &point, int depth);
    kDTreeNode *removeHelper(kDTreeNode *node, const std::vector<int> &point, int depth);
    kDTreeNode *findMinNode(kDTreeNode *node, int dimension) const;
    bool searchHelper(kDTreeNode *node, const std::vector<int> point, int depth);
    kDTreeNode *buildTreeHelper(std::vector<std::vector<int>> pointList, int depth);
    void bubbleSort(std::vector<std::vector<int>> &pointList, int dimension);
    void nearestNeighbourHelper(kDTreeNode *node, const std::vector<int> &target, kDTreeNode *&best, double &bestDistance, int depth);
    double calculateDistance(const std::vector<int> &point1, const std::vector<int> &point2);
    void kNearestNeighbourHelper(kDTreeNode *node, const std::vector<int> &target, int k,
                                 std::vector<std::pair<double, kDTreeNode *>> &nearestList,
                                 double &maxDistance, int depth);

public:
    kDTree(int k = 2);
    ~kDTree();

    const kDTree &operator=(const kDTree &other);
    kDTree(const kDTree &other);
    void setK(int vk){
        k = vk;
    }

    void inorderTraversal() const;
    void preorderTraversal() const;
    void postorderTraversal() const;
    int height() const;
    int nodeCount() const;
    int leafCount() const;

    void insert(const vector<int> &point);
    void remove(const vector<int> &point);
    bool search(const vector<int> &point);
    void buildTree(const vector<vector<int>> &pointList);
    void nearestNeighbour(const vector<int> &target, kDTreeNode *&best);
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList);
};

class kNN
{
private:
    int k;
    Dataset X_train;
    Dataset y_train;
    kDTree tree;
    kDTree treeLabels;

public:
    kNN(int k = 5);
    void fit(Dataset &X_train, Dataset &y_train);
    Dataset predict(Dataset &X_test);
    double score(const Dataset &y_test, const Dataset &y_pred);
};

// Please add more or modify as needed
