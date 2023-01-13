#include "trojanmap.h"

//-----------------------------------------------------
// TODO: Student should implement the following:
//-----------------------------------------------------
/**
 * GetLat: Get the latitude of a Node given its id. If id does not exist, return
 * -1.
 *
 * @param  {std::string} id : location id
 * @return {double}         : latitude
 */
double TrojanMap::GetLat(const std::string &id) { 
  auto res = data.find(id);
  if (res != data.end()) {
    return res->second.lat;
  } else {
    return -1;
  }
}

/**
 * GetLon: Get the longitude of a Node given its id. If id does not exist,
 * return -1.
 *
 * @param  {std::string} id : location id
 * @return {double}         : longitude
 */
double TrojanMap::GetLon(const std::string &id) {
  auto res = data.find(id);
  if (res != data.end()) {
    return res->second.lon;
  } else {
    return -1;
  }
}

/**
 * GetName: Get the name of a Node given its id. If id does not exist, return
 * "NULL".
 *
 * @param  {std::string} id : location id
 * @return {std::string}    : name
 */
std::string TrojanMap::GetName(const std::string &id) {
  auto res = data.find(id);
  if (res != data.end()) {
    return res->second.name;
  } else {
    return "NULL";
  }
}

/**
 * GetNeighborIDs: Get the neighbor ids of a Node. If id does not exist, return
 * an empty vector.
 *
 * @param  {std::string} id            : location id
 * @return {std::vector<std::string>}  : neighbor ids
 */
std::vector<std::string> TrojanMap::GetNeighborIDs(const std::string &id) {
  auto res = data.find(id);
  if (res != data.end()) {
    return res->second.neighbors;
  } else {
    return {};
  }
}

/**
 * GetID: Given a location name, return the id.
 * If the node does not exist, return an empty string.
 *
 * @param  {std::string} name          : location name
 * @return {int}  : id
 */
std::string TrojanMap::GetID(const std::string &name) {
  std::string res = "";
  for (auto it = data.begin(); it != data.end(); ++it) {
    if (name == it->second.name) {
      res = it->second.id;
      break;
    }
  }
  return res;
}

/**
 * GetPosition: Given a location name, return the position. If id does not
 * exist, return (-1, -1).
 *
 * @param  {std::string} name          : location name
 * @return {std::pair<double,double>}  : (lat, lon)
 */
std::pair<double, double> TrojanMap::GetPosition(std::string name) {
  std::pair<double, double> results(-1, -1);
  std::string id = GetID(name);
  if (!id.empty()) {
    results.first = GetLat(id);
    results.second = GetLon(id);
  }
  return results;
  
}

/**
 * CalculateEditDistance: Calculate edit distance between two location names
 *
 */
int TrojanMap::CalculateEditDistance(std::string a, std::string b) {
  int m = a.size(), n = b.size();
  std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
  for (int i = 0; i <= m; ++i) {
    dp[i][0] = i;
  }     

  for (int j = 0; j <= n; ++j) {
    dp[0][j] = j;
  }

  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
      if (std::tolower(a.at(i - 1)) != std::tolower(b.at(j - 1))) {
        dp[i][j] = 1 + std::min(std::min(dp[i - 1][j], dp[i][j - 1]), dp[i - 1][j - 1]);
      } else {
        dp[i][j] = 1 + std::min(std::min(dp[i - 1][j], dp[i][j - 1]), dp[i - 1][j - 1] - 1);
      }
    }
  }

  return dp[m][n];
}

/**
 * FindClosestName: Given a location name, return the name with smallest edit
 * distance.
 *
 * @param  {std::string} name          : location name
 * @return {std::string} tmp           : similar name
 */
std::string TrojanMap::FindClosestName(std::string name) {
  std::string res = "";
  int distance = INT_MAX;
  for (auto it = data.begin(); it != data.end(); ++it) {
    if (it->second.name.empty()) {
      continue;
    }

    int distanceTemp = CalculateEditDistance(name, it->second.name);  
    if (distance > distanceTemp) {
      distance = distanceTemp;
      res = it->second.name;
      
    }
  }

  return res;
}

/**
 * Autocomplete: Given a parital name return all the possible locations with
 * partial name as the prefix. The function should be case-insensitive.
 *
 * @param  {std::string} name          : partial name
 * @return {std::vector<std::string>}  : a vector of full names
 */
std::vector<std::string> TrojanMap::Autocomplete(std::string name) {
  if (name.empty()) {
    return {};
  }
  std::vector<std::string> results;
  std::transform(name.begin(), name.end(), name.begin(), ::tolower);
  for (auto it = data.begin(); it != data.end(); ++it) {
    if (it->second.name.empty()) {
      continue;
    }
    std::string locName = it->second.name.substr(0, name.size());
    std::transform(locName.begin(), locName.end(), locName.begin(), ::tolower);
    if (name == locName) {
      results.push_back(it->second.name);
    }
  }
  return results;
}

/**
 * GetAllCategories: Return all the possible unique location categories, i.e.
 * there should be no duplicates in the output.
 *
 * @return {std::vector<std::string>}  : all unique location categories
 */
std::vector<std::string> TrojanMap::GetAllCategories() {
  std::set<std::string> s;
  for (auto it = data.begin(); it != data.end(); ++it) {
    if (it->second.attributes.empty()) {
      continue;
    } else {
      s.insert(*it->second.attributes.begin());
    }
  }

  std::vector<std::string> res;
  for (auto it = s.begin(); it != s.end(); ++it) {
    res.push_back(*it);
  }
  
  return res;
}

/**
 * GetAllLocationsFromCategory: Return all the locations of the input category (i.e.
 * 'attributes' in data.csv). If there is no location of that category, return
 * (-1, -1). The function should be case-insensitive.
 *
 * @param  {std::string} category          : category name (attribute)
 * @return {std::vector<std::string>}     : ids
 */
std::vector<std::string> TrojanMap::GetAllLocationsFromCategory(
    std::string category) {
  std::transform(category.begin(), category.end(), category.begin(), ::tolower);
  std::vector<std::string> res;
  for (auto it = data.begin(); it != data.end(); ++it) {
    if (it->second.attributes.empty()) {
      continue;
    } else {
      std::string str = *it->second.attributes.begin();
      std::transform(str.begin(), str.end(), str.begin(), ::tolower);
      if(str == category) {
        res.push_back(it->first);
      }
    }
  }

  if (res.empty()) {
    return {"-1, -1"};
  }

  return res;
}

/**
 * GetLocationRegex: Given the regular expression of a location's name, your
 * program should first check whether the regular expression is valid, and if so
 * it returns all locations that match that regular expression.
 *
 * @param  {std::regex} location name      : the regular expression of location
 * names
 * @return {std::vector<std::string>}     : ids
 */
std::vector<std::string> TrojanMap::GetLocationRegex(std::regex location) {
  std::vector<std::string> res;
  for (auto it = data.begin(); it != data.end(); ++it) {
    if (std::regex_match(it->second.name, location)) {
      res.push_back(it->first);
    }
  }
  
  return res;
}

/**
 * CalculateDistance: Get the distance between 2 nodes.
 *
 * @param  {std::string} a  : a_id
 * @param  {std::string} b  : b_id
 * @return {double}  : distance in mile
 */
double TrojanMap::CalculateDistance(const std::string &a_id,
                                    const std::string &b_id) {
  // Do not change this function
  Node a = data[a_id];
  Node b = data[b_id];
  double dlon = (b.lon - a.lon) * M_PI / 180.0;
  double dlat = (b.lat - a.lat) * M_PI / 180.0;
  double p = pow(sin(dlat / 2), 2.0) + cos(a.lat * M_PI / 180.0) *
                                           cos(b.lat * M_PI / 180.0) *
                                           pow(sin(dlon / 2), 2.0);
  double c = 2 * asin(std::min(1.0, sqrt(p)));
  return c * 3961;
}

/**
 * CalculatePathLength: Calculates the total path length for the locations
 * inside the vector.
 *
 * @param  {std::vector<std::string>} path : path
 * @return {double}                        : path length
 */
double TrojanMap::CalculatePathLength(const std::vector<std::string> &path) {
  // Do not change this function
  double sum = 0;
  for (int i = 0; i < int(path.size()) - 1; i++) {
    sum += CalculateDistance(path[i], path[i + 1]);
  }
  return sum;
}

/**
 * CalculateShortestPath_Dijkstra: Given 2 locations, return the shortest path
 * which is a list of id. Hint: Use priority queue.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_Dijkstra(
    std::string location1_name, std::string location2_name) {
  std::string location1_id = GetID(location1_name), location2_id = GetID(location2_name);
  std::priority_queue<std::vector<std::pair<double, std::string>>, std::vector<std::pair<double, std::string>>, std::greater<std::pair<double, std::string>>> pq;
  std::unordered_map<std::string, bool> visited;
  std::unordered_map<std::string, std::string> prev;
  std::unordered_map<std::string, double> distance;
  std::vector<std::string> path;
  for (auto d : data) {
    distance[d.first] = INT_MAX;
  }

  distance[location1_id] = 0;
  pq.push({0, location1_id});
  while (!pq.empty() && !visited[location2_id]) {
    auto u = pq.top().second;
    pq.pop();
    if (visited[u]) {
      continue;
    }

    visited[u] = true;
    for (auto v : data[u].neighbors) {
      double weight = CalculateDistance(u, v);
      if (distance[v] > distance[u] + weight) {
        distance[v] = distance[u] + weight;
        pq.push({distance[v], v});
        prev[v] = u;
      }
    }
  }

  path.push_back(location2_id);
  while (location2_id != location1_id) {
    location2_id = prev[location2_id];
    path.push_back(location2_id);
  }
  
  std::reverse(path.begin(), path.end());
  return path;
}

/**
 * CalculateShortestPath_Bellman_Ford: Given 2 locations, return the shortest
 * path which is a list of id. Hint: Do the early termination when there is no
 * change on distance.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_Bellman_Ford(
    std::string location1_name, std::string location2_name) {
  std::string location1_id = GetID(location1_name), location2_id = GetID(location2_name);
  std::unordered_map<std::string, std::string> prev;
  std::unordered_map<std::string, double> distance;
  std::vector<std::string> path;
  for (auto loc : data) {
    distance[loc.first] = INT_MAX;
  }
  distance[location1_id] = 0;
  for (int i = 0; i < distance.size() - 1; ++i) {
    bool flag = true;
    for (auto u : data) {
      for (auto adj : u.second.neighbors) {
        if (distance[adj] == INT_MAX && distance[u.first] == INT_MAX) {
          continue;
        } else {
          double weight = CalculateDistance(u.first, adj);
          if (distance[adj] > distance[u.first] + weight) {
            distance[adj] = distance[u.first] + weight;
            prev[adj] = u.first;
            flag = false;
          }
        }
      }
    }
    if (flag) break;
  }

  path.push_back(location2_id);
  while (location2_id != location1_id) {
    location2_id = prev[location2_id];
    path.push_back(location2_id);
  }
  
  std::reverse(path.begin(), path.end());
  return path;
}

/**
 * Traveling salesman problem: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path
 */
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_Brute_force(
                                    std::vector<std::string> location_ids) {
  std::pair<double, std::vector<std::vector<std::string>>> records;
  return records;
}

std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_Backtracking(
                                    std::vector<std::string> location_ids) {
  std::pair<double, std::vector<std::vector<std::string>>> records;
  return records;
}

// Hint: https://en.wikipedia.org/wiki/2-opt
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_2opt(
      std::vector<std::string> location_ids){
  std::pair<double, std::vector<std::vector<std::string>>> records;
  return records;
}

/**
 * Given CSV filename, it read and parse locations data from CSV file,
 * and return locations vector for topological sort problem.
 *
 * @param  {std::string} locations_filename     : locations_filename
 * @return {std::vector<std::string>}           : locations
 */
std::vector<std::string> TrojanMap::ReadLocationsFromCSVFile(
    std::string locations_filename) {
  std::vector<std::string> location_names_from_csv;
  std::fstream fin;
  fin.open(locations_filename, std::ios::in);
  std::string line, word;
  getline(fin, line);
  while (getline(fin, word)) {
    location_names_from_csv.push_back(word);
  }
  fin.close();
  return location_names_from_csv;
}

/**
 * Given CSV filenames, it read and parse dependencise data from CSV file,
 * and return dependencies vector for topological sort problem.
 *
 * @param  {std::string} dependencies_filename     : dependencies_filename
 * @return {std::vector<std::vector<std::string>>} : dependencies
 */
std::vector<std::vector<std::string>> TrojanMap::ReadDependenciesFromCSVFile(
    std::string dependencies_filename) {
  std::vector<std::vector<std::string>> dependencies_from_csv;
  std::fstream fin;
  fin.open(dependencies_filename, std::ios::in);
  std::string line, word;
  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);
    std::vector<std::string> dependency;
    while (getline(s, word, ',')) {
      dependency.push_back(word);
    }
    dependencies_from_csv.push_back(dependency);
  }
  fin.close();
  return dependencies_from_csv;
}

/**
 * DeliveringTrojan: Given a vector of location names, it should return a
 * sorting of nodes that satisfies the given dependencies. If there is no way to
 * do it, return a empty vector.
 *
 * @param  {std::vector<std::string>} locations                     : locations
 * @param  {std::vector<std::vector<std::string>>} dependencies     :
 * prerequisites
 * @return {std::vector<std::string>} results                       : results
 */
std::vector<std::string> TrojanMap::DeliveringTrojan(
    std::vector<std::string> &locations,
    std::vector<std::vector<std::string>> &dependencies) {
  std::vector<std::string> result;
  std::unordered_map<std::string, std::unordered_set<std::string>> g;
  std::unordered_map<std::string, int> indegree;
  std::queue<std::string> q;
  for (auto loc : locations) {
    g[loc] = {};
    indegree[loc] = 0;
  }

  for (auto dep : dependencies) {
    g[dep[0]].insert(dep[1]);
  }

  for (auto node : g) {
    for (auto adj : node.second) {
      indegree[adj]++;
    }
  }

  for (auto in : indegree) {
    if (in.second == 0) {
      q.push(in.first);
    }
  }

  while (!q.empty()) {
    auto u = q.front();
    q.pop();
    result.push_back(u);
    for (auto adj : g[u]) {
      if (--indegree[adj] == 0) {
        q.push(adj);
      }
    }
  }

  return result;
}

/**
 * inSquare: Give a id retunr whether it is in square or not.
 *
 * @param  {std::string} id            : location id
 * @param  {std::vector<double>} square: four vertexes of the square area
 * @return {bool}                      : in square or not
 */
bool TrojanMap::inSquare(std::string id, std::vector<double> &square) {
  auto n = data[id];
  if (n.lon > square[0] && n.lon < square[1] && n.lat < square[2] && n.lat > square[3]) {
    return true;
  } else {
    return false;
  }
  
}

/**
 * GetSubgraph: Give four vertexes of the square area, return a list of location
 * ids in the squares
 *
 * @param  {std::vector<double>} square         : four vertexes of the square
 * area
 * @return {std::vector<std::string>} subgraph  : list of location ids in the
 * square
 */
std::vector<std::string> TrojanMap::GetSubgraph(std::vector<double> &square) {
  // include all the nodes in subgraph
  std::vector<std::string> subgraph;
  for (auto it = data.begin(); it != data.end(); ++it) {
    if (inSquare(it->first, square)) {
      subgraph.push_back(it->first);
    }
  }
  
  return subgraph;
}

/**
 * Cycle Detection: Given four points of the square-shape subgraph, return true
 * if there is a cycle path inside the square, false otherwise.
 *
 * @param {std::vector<std::string>} subgraph: list of location ids in the
 * square
 * @param {std::vector<double>} square: four vertexes of the square area
 * @return {bool}: whether there is a cycle or not
 */
bool TrojanMap::CycleDetection(std::vector<std::string> &subgraph, std::vector<double> &square) {
  std::unordered_map<std::string, int> status;
  std::unordered_map<std::string, std::string> parent;
  for (auto n : subgraph) {
    status[n] = 0;
  }
  for (auto start_node : subgraph) {
    if (ReachesACycleHelper(start_node, status, parent, square)) {
      return true;
    }

    parent.clear();
    status.clear();
  }

  return false;
}

bool TrojanMap::ReachesACycleHelper(std::string start_node, std::unordered_map<std::string, int> &status, std::unordered_map<std::string, std::string> &parent, std::vector<double> &square) {
  status[start_node] = 1;
  for (auto adj : data[start_node].neighbors) {
    if (!inSquare(adj, square)) {
      continue;
    }

    if (status[adj] == 1 && parent[start_node] != adj) {
      return true;
    } else {
      if (status[adj] == 0) {
        parent[adj] = start_node;
        ReachesACycleHelper(adj, status, parent, square);
      }
    }
  }
  status[start_node] = 2;
  return false;
}

/**
 * FindNearby: Given a class name C, a location name L and a number r,
 * find all locations in class C on the map near L with the range of r and
 * return a vector of string ids
 *
 * @param {std::string} className: the name of the class
 * @param {std::string} locationName: the name of the location
 * @param {int} r: search radius
 * @param {int} k: search numbers
 * @return {std::vector<std::string>}: location name that meets the requirements
 */
std::vector<std::string> TrojanMap::FindNearby(std::string attributesName, std::string name, double r, int k) {
  std::vector<std::string> res;
  return res;
}

/**
 * Shortest Path to Visit All Nodes: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of
 * total distance and the all the progress to get final path
 */
std::vector<std::string> TrojanMap::TrojanPath(
      std::vector<std::string> &location_names) {
    std::vector<std::string> res;
    return res;
}

/**
 * CreateGraphFromCSVFile: Read the map data from the csv file
 *
 */
void TrojanMap::CreateGraphFromCSVFile() {
  // Do not change this function
  std::fstream fin;
  fin.open("src/lib/data.csv", std::ios::in);
  std::string line, word;

  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);

    Node n;
    int count = 0;
    while (getline(s, word, ',')) {
      word.erase(std::remove(word.begin(), word.end(), '\''), word.end());
      word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '{'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '}'), word.end());
      if (count == 0)
        n.id = word;
      else if (count == 1)
        n.lat = stod(word);
      else if (count == 2)
        n.lon = stod(word);
      else if (count == 3)
        n.name = word;
      else {
        word.erase(std::remove(word.begin(), word.end(), ' '), word.end());
        if (isalpha(word[0])) n.attributes.insert(word);
        if (isdigit(word[0])) n.neighbors.push_back(word);
      }
      count++;
    }
    data[n.id] = n;
  }
  fin.close();
}
