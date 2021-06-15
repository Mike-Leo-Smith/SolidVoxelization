#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

template<typename T> struct Vec3
{
    T x;
    T y;
    T z;
    Vec3(T _x = 0, T _y = 0, T _z = 0): x(_x), y(_y), z(_z) {}

    T dot(Vec3<T> right) { return x * right.x + y * right.y + z * right.z; }
    Vec3<T> operator-(Vec3<T> right) { return Vec3<T>(x - right.x, y - right.y, z - right.z); }
    Vec3<T> operator+(Vec3<T> right) { return Vec3<T>(x + right.x, y + right.y, z + right.z); }
    Vec3<T> operator*(T right) { return Vec3<T>(x * right, y * right, z * right); }

    friend std::ostream& operator<<(std::ostream& os, Vec3<T> vec)
    {
        os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
        return os;
    }
};

struct Segment
{
    int x;
    int y;
    int z_min;
    int z_max;

    Segment(int _x = 0, int _y = 0, int _z_min = 0, int _z_max = 0) : x(_x), y(_y), z_min(_z_min), z_max(_z_max) {}
};

struct EmbreeRay {
    Vec3<float> origin;
    Vec3<float> direction;
    EmbreeRay(Vec3<float> _origin = Vec3<float>(), Vec3<float> _direction = Vec3<float>()): origin(_origin), direction(_direction) {}
};

struct OctNode
{
    OctNode* child_xpypzp; //xp is x >= midpoint.x, xn is x < midpoint.x, others the same
    OctNode* child_xpypzn;
    OctNode* child_xpynzp;
    OctNode* child_xpynzn;
    OctNode* child_xnypzp;
    OctNode* child_xnypzn;
    OctNode* child_xnynzp;
    OctNode* child_xnynzn;

    int x_min, x_max; //[x_min, x_max), others the same
    int y_min, y_max;
    int z_min, z_max;

    Vec3<int> midpoint;

    enum Status { EMPTY, HALF_FULL, FULL } status;

    OctNode() {}

    friend std::ostream& operator<<(std::ostream& os, OctNode& node)
    {
        os << "(" << node.x_min << " to " << node.x_max << ", " << node.y_min << " to " << node.y_max << ", " << node.z_min << " to " << node.z_max << ", ";
        switch (node.status)
        {
        case OctNode::EMPTY:
            std::cout << "EMPTY)";
            break;
        case OctNode::HALF_FULL:
            std::cout << "HALF)";
            break;
        case OctNode::FULL:
            std::cout << "FULL)";
            break;
        default:
            std::cout << "SHIT!)";
            break;
        }
        return os;
    }
};

OctNode* build_octnode(std::vector<Vec3<int>>& points, int x_min, int x_max, int y_min, int y_max, int z_min, int z_max)
{
    OctNode* node = new OctNode;

    node->x_min = x_min;
    node->x_max = x_max;
    node->y_min = y_min;
    node->y_max = y_max;
    node->z_min = z_min;
    node->z_max = z_max;

    //std::cout << "build point ";
    //for (auto p : points)
    //    std::cout << p << " ";
    //std::cout << " whose bounding box is " << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << "\n";

    if (points.empty() || x_min >= x_max || y_min >= y_max || z_min >= z_max)
    {
        node->status = OctNode::EMPTY;
        return node;
    }
    if (points.size() == (x_max - x_min) * (y_max - y_min) * (z_max - z_min))
    {
        node->status = OctNode::FULL;
        return node;
    }

    node->status = OctNode::HALF_FULL;

    Vec3<int> midpoint((x_min + x_max) / 2, (y_min + y_max) / 2, (z_min + z_max) / 2);
    node->midpoint = midpoint;

    std::vector<Vec3<int>> points_xpypzp;
    std::vector<Vec3<int>> points_xpypzn;
    std::vector<Vec3<int>> points_xpynzp;
    std::vector<Vec3<int>> points_xpynzn;
    std::vector<Vec3<int>> points_xnypzp;
    std::vector<Vec3<int>> points_xnypzn;
    std::vector<Vec3<int>> points_xnynzp;
    std::vector<Vec3<int>> points_xnynzn;

    for (auto& point : points)
    {
        if (point.x >= midpoint.x && point.y >= midpoint.y && point.z >= midpoint.z)
            points_xpypzp.push_back(point);
        else if (point.x >= midpoint.x && point.y >= midpoint.y && point.z < midpoint.z)
            points_xpypzn.push_back(point);
        else if (point.x >= midpoint.x && point.y < midpoint.y && point.z >= midpoint.z)
            points_xpynzp.push_back(point);
        else if (point.x >= midpoint.x && point.y < midpoint.y && point.z < midpoint.z)
            points_xpynzn.push_back(point);
        else if (point.x < midpoint.x && point.y >= midpoint.y && point.z >= midpoint.z)
            points_xnypzp.push_back(point);
        else if (point.x < midpoint.x && point.y >= midpoint.y && point.z < midpoint.z)
            points_xnypzn.push_back(point);
        else if (point.x < midpoint.x && point.y < midpoint.y && point.z >= midpoint.z)
            points_xnynzp.push_back(point);
        else if (point.x < midpoint.x && point.y < midpoint.y && point.z < midpoint.z)
            points_xnynzn.push_back(point);
    }

    node->child_xpypzp = build_octnode(points_xpypzp, midpoint.x, x_max, midpoint.y, y_max, midpoint.z, z_max);
    node->child_xpypzn = build_octnode(points_xpypzn, midpoint.x, x_max, midpoint.y, y_max, z_min, midpoint.z);
    node->child_xpynzp = build_octnode(points_xpynzp, midpoint.x, x_max, y_min, midpoint.y, midpoint.z, z_max);
    node->child_xpynzn = build_octnode(points_xpynzn, midpoint.x, x_max, y_min, midpoint.y, z_min, midpoint.z);
    node->child_xnypzp = build_octnode(points_xnypzp, x_min, midpoint.x, midpoint.y, y_max, midpoint.z, z_max);
    node->child_xnypzn = build_octnode(points_xnypzn, x_min, midpoint.x, midpoint.y, y_max, z_min, midpoint.z);
    node->child_xnynzp = build_octnode(points_xnynzp, x_min, midpoint.x, y_min, midpoint.y, midpoint.z, z_max);
    node->child_xnynzn = build_octnode(points_xnynzn, x_min, midpoint.x, y_min, midpoint.y, z_min, midpoint.z);

    return node;
}

std::vector<Vec3<int>> segment_to_points(std::vector<Segment>& segments)
{
    std::vector<Vec3<int>> points;

    for (auto& segment : segments)
    {
        for (int z = segment.z_min; z <= segment.z_max; z++)
            points.push_back(Vec3<int>(segment.x, segment.y, z));
    }

    return points;
}

OctNode* build_octtree(std::vector<Segment>& segments, int res)
{
    auto points = segment_to_points(segments);

//    int x_min = INT32_MAX, y_min = INT32_MAX, z_min = INT32_MAX;
//    int x_max = INT32_MIN, y_max = INT32_MIN, z_max = INT32_MIN;
//
//    for (auto& point : points)
//    {
//        if (point.x < x_min)
//            x_min = point.x;
//        if (point.x > x_max)
//            x_max = point.x;
//        if (point.y < y_min)
//            y_min = point.y;
//        if (point.y > y_max)
//            y_max = point.y;
//        if (point.z < z_min)
//            z_min = point.z;
//        if (point.z > z_max)
//            z_max = point.z;
//    }

    return build_octnode(points, 0, res, 0, res, 0, res);
}

void show_octtree(OctNode* node)
{
    if (node->status == OctNode::EMPTY)
    {
        std::cout << *node << " is empty\n";
        return;
    }
    if (node->status == OctNode::FULL)
    {
        std::cout << *node << " is full\n";
        return;
    }

    std::cout << *node << " has child " << *(node->child_xpypzp) << "\n";
    std::cout << *node << " has child " << *(node->child_xpypzn) << "\n";
    std::cout << *node << " has child " << *(node->child_xpynzp) << "\n";
    std::cout << *node << " has child " << *(node->child_xpynzn) << "\n";
    std::cout << *node << " has child " << *(node->child_xnypzp) << "\n";
    std::cout << *node << " has child " << *(node->child_xnypzn) << "\n";
    std::cout << *node << " has child " << *(node->child_xnynzp) << "\n";
    std::cout << *node << " has child " << *(node->child_xnynzn) << "\n\n";

    show_octtree(node->child_xpypzp);
    show_octtree(node->child_xpypzn);
    show_octtree(node->child_xpynzp);
    show_octtree(node->child_xpynzn);
    show_octtree(node->child_xnypzp);
    show_octtree(node->child_xnypzn);
    show_octtree(node->child_xnynzp);
    show_octtree(node->child_xnynzn);
}

struct obj_v
{
    int x;
    int y;
    int z;
    obj_v(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {}
};

struct obj_f
{
    int v1;
    int v2;
    int v3;
    int v4;
    obj_f(int _v1, int _v2, int _v3, int _v4) : v1(_v1), v2(_v2), v3(_v3), v4(_v4) {}
};

void octtree_to_obj(OctNode* node, std::vector<obj_v>& vertices, std::vector<obj_f>& faces)
{
    if (node->status == OctNode::EMPTY)
        return;

    if (node->status == OctNode::FULL)
    {
        int base = vertices.size();

        vertices.push_back(obj_v(node->x_min, node->y_min, node->z_min)); //base + 1
        vertices.push_back(obj_v(node->x_min, node->y_min, node->z_max)); //base + 2
        vertices.push_back(obj_v(node->x_min, node->y_max, node->z_min)); //base + 3
        vertices.push_back(obj_v(node->x_min, node->y_max, node->z_max)); //base + 4
        vertices.push_back(obj_v(node->x_max, node->y_min, node->z_min)); //base + 5
        vertices.push_back(obj_v(node->x_max, node->y_min, node->z_max)); //base + 6
        vertices.push_back(obj_v(node->x_max, node->y_max, node->z_min)); //base + 7
        vertices.push_back(obj_v(node->x_max, node->y_max, node->z_max)); //base + 8

        faces.push_back(obj_f(base + 1, base + 2, base + 3, base + 4));
        faces.push_back(obj_f(base + 5, base + 8, base + 7, base + 8));
        faces.push_back(obj_f(base + 1, base + 3, base + 5, base + 7));
        faces.push_back(obj_f(base + 2, base + 4, base + 6, base + 8));
        faces.push_back(obj_f(base + 1, base + 2, base + 5, base + 6));
        faces.push_back(obj_f(base + 3, base + 4, base + 7, base + 8));

        return;
    }

    //HALF_FULL
    octtree_to_obj(node->child_xpypzp, vertices, faces);
    octtree_to_obj(node->child_xpypzn, vertices, faces);
    octtree_to_obj(node->child_xpynzp, vertices, faces);
    octtree_to_obj(node->child_xpynzn, vertices, faces);
    octtree_to_obj(node->child_xnypzp, vertices, faces);
    octtree_to_obj(node->child_xnypzn, vertices, faces);
    octtree_to_obj(node->child_xnynzp, vertices, faces);
    octtree_to_obj(node->child_xnynzn, vertices, faces);
}

void obj(OctNode* node, std::string filename)
{
    std::ofstream out(filename);

    std::vector<obj_v> vertices;
    std::vector<obj_f> faces;

    octtree_to_obj(node, vertices, faces);

    for (auto& vertice : vertices)
        out << "v " << vertice.x << " " << vertice.y << " " << vertice.z << "\n";

    for (auto& face : faces)
        out << "f " << face.v1 << " " << face.v2 << " " << face.v3 << " " << face.v4 << "\n";

    out.close();
}

float intersection_f(EmbreeRay ray, EmbreeRay normal)
{
    Vec3<float> t = ray.origin - normal.origin;
    Vec3<float>& n = normal.direction;
    Vec3<float>& dir = ray.direction;
    return -t.dot(n) / dir.dot(n);
}

float intersection_f_x(EmbreeRay ray, int x, int y_min, int y_max, int z_min, int z_max)
{
    EmbreeRay normal(Vec3<float>(x, y_min, z_min), Vec3<float>(1, 0, 0));
    float f = intersection_f(ray, normal);
    
    Vec3<float> point = ray.origin + ray.direction * f;
    auto point_y = point.y;
    auto point_z = point.z;

    if (y_min <= point_y && point_y < y_max && z_min <= point_z && point_z < z_max)
        return f;
    else
        return -1;
}

float intersection_f_y(EmbreeRay ray, int x_min, int x_max, int y, int z_min, int z_max)
{
    EmbreeRay normal(Vec3<float>(x_min, y, z_min), Vec3<float>(0, 1, 0));
    float f = intersection_f(ray, normal);

    Vec3<float> point = ray.origin + ray.direction * f;
    auto point_x = point.x;
    auto point_z = point.z;

    if (x_min <= point_x && point_x < x_max && z_min <= point_z && point_z < z_max)
        return f;
    else
        return -1;
}

float intersection_f_z(EmbreeRay ray, int x_min, int x_max, int y_min, int y_max, int z)
{
    EmbreeRay normal(Vec3<float>(x_min, y_min, z), Vec3<float>(0, 0, 1));
    float f = intersection_f(ray, normal);

    Vec3<float> point = ray.origin + ray.direction * f;
    auto point_x = point.x;
    auto point_y = point.y;

    if (x_min <= point_x && point_x < x_max && y_min <= point_y && point_y < y_max)
        return f;
    else
        return -1;
}

Vec3<float> intersection(EmbreeRay ray, OctNode* node)
{
    if (node->status == OctNode::EMPTY)
        return Vec3<float>(-1, -1, -1);

    std::vector<float> fs;

    fs.push_back(intersection_f_x(ray, node->x_min, node->y_min, node->y_max, node->z_min, node->z_max));
    fs.push_back(intersection_f_x(ray, node->midpoint.x, node->y_min, node->y_max, node->z_min, node->z_max));
    fs.push_back(intersection_f_x(ray, node->x_max, node->y_min, node->y_max, node->z_min, node->z_max));
    fs.push_back(intersection_f_y(ray, node->x_min, node->x_max, node->y_min, node->z_min, node->z_max));
    fs.push_back(intersection_f_y(ray, node->x_min, node->x_max, node->midpoint.y, node->z_min, node->z_max));
    fs.push_back(intersection_f_y(ray, node->x_min, node->x_max, node->y_max, node->z_min, node->z_max));
    fs.push_back(intersection_f_z(ray, node->x_min, node->x_max, node->y_min, node->y_max, node->z_min));
    fs.push_back(intersection_f_z(ray, node->x_min, node->x_max, node->y_min, node->y_max, node->midpoint.z));
    fs.push_back(intersection_f_z(ray, node->x_min, node->x_max, node->y_min, node->y_max, node->z_max));
    
    for (auto it = fs.begin(); it != fs.end();)
    {
        if (*it < 0)
            it = fs.erase(it);
        else
            it++;
    }

    if (fs.empty())
        return Vec3<float>(-1, -1, -1);

    std::sort(fs.begin(), fs.end());

    std::cout << ray.origin << " " << ray.direction << " intersects " << *node << "at:\n";
    for (auto f : fs)
    {
        auto point = ray.origin + ray.direction * f;
        std::cout << "\t" << point << "\n";
    }

    if (node->status == OctNode::FULL)
    {
        return ray.origin + ray.direction * fs[0];
    }

    Vec3<float> closest_point;

    for (int i = 0; i < fs.size() - 1; i++)
    {
        Vec3<float> mid = ray.origin + ray.direction * ((fs[i] + fs[i + 1]) / 2);
        auto mid_x = mid.x;
        auto mid_y = mid.y;
        auto mid_z = mid.z;

        if (mid_x >= node->midpoint.x && mid_y >= node->midpoint.y && mid_z >= node->midpoint.z && (closest_point = intersection(ray, node->child_xpypzp)).x >= 0)
            return closest_point;
        if (mid_x >= node->midpoint.x && mid_y >= node->midpoint.y && mid_z < node->midpoint.z && (closest_point = intersection(ray, node->child_xpypzn)).x >= 0)
            return closest_point;
        if (mid_x >= node->midpoint.x && mid_y < node->midpoint.y && mid_z >= node->midpoint.z && (closest_point = intersection(ray, node->child_xpynzp)).x >= 0)
            return closest_point;
        if (mid_x >= node->midpoint.x && mid_y < node->midpoint.y && mid_z < node->midpoint.z && (closest_point = intersection(ray, node->child_xpynzn)).x >= 0)
            return closest_point;
        if (mid_x < node->midpoint.x && mid_y >= node->midpoint.y && mid_z >= node->midpoint.z && (closest_point = intersection(ray, node->child_xnypzp)).x >= 0)
            return closest_point;
        if (mid_x < node->midpoint.x && mid_y >= node->midpoint.y && mid_z < node->midpoint.z && (closest_point = intersection(ray, node->child_xnypzn)).x >= 0)
            return closest_point;
        if (mid_x < node->midpoint.x && mid_y < node->midpoint.y && mid_z >= node->midpoint.z && (closest_point = intersection(ray, node->child_xnynzp)).x >= 0)
            return closest_point;
        if (mid_x < node->midpoint.x && mid_y < node->midpoint.y && mid_z < node->midpoint.z && (closest_point = intersection(ray, node->child_xnynzn)).x >= 0)
            return closest_point;
    }

    return Vec3<float>(-1, -1, -1);
}


int main()
{
    std::ifstream in("test-128.txt");

    int res;
    in >> res;
    
    int segment_num;
    in >> segment_num;

    std::cout << "resolution: " << res << ", segment_num: " << segment_num << "\n";

    std::vector<Segment> segments;
    for (int i = 0; i < segment_num; i++)
    {
        int x, y, z_min, z_max;
        in >> x >> y >> z_min >> z_max;
        segments.push_back(Segment(x, y, z_min, z_max));
    }

    OctNode* root = build_octtree(segments, res);

    //show_octtree(root);

    obj(root, "out.obj");

//    while (1)
//    {
//        int origin_x, origin_y, origin_z;
//        int direction_x, direction_y, direction_z;
//        std::cin >> origin_x >> origin_y >> origin_z >> direction_x >> direction_y >> direction_z;
//        Ray ray(Vec3<float>(origin_x, origin_y, origin_z), Vec3<float>(direction_x, direction_y, direction_z));
//        auto point = intersection(ray, root);
//        std::cout << point.x << " " << point.y << " " << point.z << "\n";
//    }
}
