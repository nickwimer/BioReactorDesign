{
    std::vector<double> solnvec(microbemodel::nvars);
   /*
    double advance_time=5.0; //hours
    double dt_react=advance_time/50.0; //hours
    double final_time=30.0; //hours
    
    solnvec[microbemodel::X]=0.5;
    solnvec[microbemodel::O2]=0.214;
    solnvec[microbemodel::G]=73.0;
    solnvec[microbemodel::M]=0.0;
 */
    
    double advance_time=2.0; //hours
    double dt_react=0.01; //hours
    double final_time=34.0; //hours

    solnvec[microbemodel::X]=0.329754101;
    solnvec[microbemodel::O2]=0.05;
    solnvec[microbemodel::G]=547.2541796;
    solnvec[microbemodel::M]=0.173782;
    

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
        microbemodel::advance(solnvec,microbemodel::nvars,current_time,advance_time,dt_react, false);
        outfile<<current_time<<"\t";
        for(int i=0;i<microbemodel::nvars;i++)
        {
            outfile<<solnvec[i]<<"\t";
        }
        outfile<<"\n";
    }
    outfile.close();


    std::ofstream outfile_ml("wellmixed_soln_ml.dat");

    solnvec[microbemodel::X]=0.329754101;
    solnvec[microbemodel::O2]=0.05;
    solnvec[microbemodel::G]=547.2541796;
    solnvec[microbemodel::M]=0.173782;

    current_time=0.0;
    outfile_ml<<current_time<<"\t";
    for(int i=0;i<microbemodel::nvars;i++)
    {
        outfile_ml<<solnvec[i]<<"\t";
    }
    outfile_ml<<"\n";
    while(current_time < final_time)
    {
        current_time += advance_time;
        microbemodel::advance(solnvec,microbemodel::nvars,current_time,advance_time,dt_react, true);
        outfile_ml<<current_time<<"\t";
        for(int i=0;i<microbemodel::nvars;i++)
        {
            outfile_ml<<solnvec[i]<<"\t";
        }
        outfile_ml<<"\n";
    }
    outfile_ml.close();
}
