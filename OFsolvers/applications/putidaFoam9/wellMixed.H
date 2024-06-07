{
    std::vector<double> solnvec(microbemodel::nvars);
    double advance_time=0.5; //hours
    double dt_react=advance_time/50.0; //hours
    double final_time=24.0; //hours
    
    solnvec[microbemodel::X]=0.5;
    solnvec[microbemodel::O2]=0.214;
    solnvec[microbemodel::G]=500.0;
    solnvec[microbemodel::M]=0.0;

    std::ofstream outfile("wellmixed_soln.dat");

    double current_time=0.0;
    outfile<<current_time<<"\t";
    for(int i=0;i<microbemodel::nvars;i++)
    {
        outfile<<solnvec[i]<<"\t";
    }
    outfile<<"\n";
    while(current_time < final_time)
    {
        current_time += advance_time;
        microbemodel::advance(solnvec,microbemodel::nvars,current_time,advance_time,dt_react);
        outfile<<current_time<<"\t";
        for(int i=0;i<microbemodel::nvars;i++)
        {
            outfile<<solnvec[i]<<"\t";
        }
        outfile<<"\n";
    }
    outfile.close();
}
