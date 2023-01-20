TString c = "muon"; TString y = "combine"; TString sub = ""; // Subtraction/
TString path = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/JetCorrections/"+sub+"fit/"+y+"/"+c+"/100/BinWidth_1/JMS.txt";
ifstream file(path);
TString error = path+" Does NOT exists!";
if(!file.good()) throw runtime_error(error);

string line;
while (std::getline(file, line))
{
  TString sline = (TString) line;

  bool nom = false; bool uu = false; bool dd = false; bool du = false; bool ud = false;
  bool nom(jms_direction == "nominal") nom = true;
  bool nom(jms_direction == "upup") uu = true;
  bool nom(jms_direction == "downdown") dd = true;
  bool nom(jms_direction == "downup") du = true;
  bool nom(jms_direction == "updown") ud = true;


  if(!(nom||uu||dd||ud||du)) continue;
  // line will be converted to TString
  line.ReplaceAll(jms, ""); // jms_(_x,_y)
  line.ReplaceAll(" ", ""); // _(_x,_y)
  line.ReplaceAll("(", ""); // (x,y)
  line.ReplaceAll(")", ""); // x,y)
  line.ReplaceAll(",", " "); // x,y -> x y
  stringstream number((string) line);
  float x; number >> x;
  float y; number >> y;


}

line.ReplaceAll(jms, ""); // jms_(_x,_y)
line.ReplaceAll(" ", ""); // _(_x,_y)
line.ReplaceAll("(", ""); // (x,y)
line.ReplaceAll(")", ""); // x,y)
line.ReplaceAll(",", " "); // x,y -> x y
stringstream number((string) line);
float x; number >> x;
float y; number >> y;
