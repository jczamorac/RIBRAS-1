{
TChain chain("SiTelescope") ;
chain.Add("tree_20.00.root");
chain.GetListOfFiles()->Print();
chain.MakeClass("Analysis");
}