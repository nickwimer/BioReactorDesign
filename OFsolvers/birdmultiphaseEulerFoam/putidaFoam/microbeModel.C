#include"microbeModel.H"
#include<map>
#include<torch/torch.h>
#include<torch/script.h>

namespace microbemodel
{
    const int X=0;
    const int O2=1;
    const int G=2;
    const int M=3;
    // const int CO2=4;
    const int nvars=4;
    std::map<std::string, int> sp_keys = {{"putida.liquid", 0},
						  {"O2.liquid", 1},
						  {"C6H12O6.liquid", 2},
						  {"C6H4O4.liquid", 3}};
                                                  // {"CO2.liquid", 4}};

    // Parameters fit from Will 0d Fixed Stir Rate Data
    // Along with Initial parameter guess provided by Jamshidzadeh et al. (2025) 
    
    // const double y_xs = 0.04119077;
    // const double y_ms = 0.31743851;
    // const double y_os = 2.0319953;

    // const double x_max = 7.9;
    // const double q_max = 8.26529425;
    // const double o2_max = 0.188125;

    // const double K_o = 0.10375191; 
    // const double K_s = 0.40010819;

    // const double kLa = 540; // Not used here; kLa calculated in OpenFOAM and used with OTR modeling there

    // Parameters are now passed in from OpenFOAM as MM_params vectors
    std::map<std::string, int> MM_param_keys = {{"y_xs", 0},
                          {"y_ms", 1},
                          {"y_os", 2},
                          {"x_max", 3},
                          {"q_max", 4},
                          {"o2_max", 5},
                          {"K_o", 6},
                          {"K_s", 7},
                          {"kLa", 8}};

    double X_avg=0.0;
    double G_avg=0.0;
    double M_avg=0.0;
    double O2_avg=0.0;

    torch::jit::script::Module model;
    bool model_loaded = false;

    std::vector<double> load_scaler(const std::string& filename)
    {
        std::ifstream file(filename);
        std::vector<double> values;
        double value;
        while (file >> value)
        {
            values.push_back(value);
        }
        return values;
    }

    std::vector<std::vector<double>> load_points(const std::string& filename)
    {
        std::ifstream file(filename);
        std::vector<std::vector<double>> points;
        std::string line;

        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::vector<double> point;
            double value;
            while (iss >> value)
            {
                point.push_back(value);
            }
            points.push_back(point);
        }
        return points;
    }

    std::vector<std::vector<int>> load_simplices(const std::string& filename)
    {
        std::ifstream file(filename);
        std::vector<std::vector<int>> simplices;
        std::string line;

        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::vector<int> simplex;
            int index;
            while (iss >> index)
            {
                simplex.push_back(index);
            }
            simplices.push_back(simplex);
        }
        return simplices;
    }

    std::vector<double> x_mean = load_scaler("scaler_x_mean.txt");
    std::vector<double> x_scale = load_scaler("scaler_x_scale.txt");
    std::vector<double> y_mean = load_scaler("scaler_y_mean.txt");
    std::vector<double> y_scale = load_scaler("scaler_y_scale.txt");

    std::vector<std::vector<double>> tri_points = load_points("tri_points.txt");
    std::vector<std::vector<int>> tri_simplices = load_simplices("tri_simplices.txt");

    std::vector<double> scale_input(const std::vector<double>& input,
                                const std::vector<double>& mean,
                                const std::vector<double>& scale)
    {
        std::vector<double> scaled;
        for (size_t i = 0; i < input.size(); ++i)
        {
            scaled.push_back((input[i] - mean[i]) / scale[i]);
        }
        return scaled;
    }

    std::vector<double> unscale_output(const std::vector<double>& output,
                                   const std::vector<double>& mean,
                                   const std::vector<double>& scale)
    {
        std::vector<double> unscaled;
        for (size_t i = 0; i < output.size(); ++i)
        {
            unscaled.push_back(output[i] * scale[i] + mean[i]);
        }
        return unscaled;
    }

    bool is_point_in_triangle(const std::vector<double>& p,
                          const std::vector<double>& a,
                          const std::vector<double>& b,
                          const std::vector<double>& c)
    {
        double detT = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        double alpha = ((b[0] - p[0]) * (c[1] - p[1]) - (b[1] - p[1]) * (c[0] - p[0])) / detT;
        double beta = ((c[0] - p[0]) * (a[1] - p[1]) - (c[1] - p[1]) * (a[0] - p[0])) / detT;
        double gamma = 1.0 - alpha - beta;

        return (alpha >= 0 && beta >= 0 && gamma >= 0);
    }

    bool is_point_in_convex_hull(const std::vector<double>& point,
                                const std::vector<std::vector<double>>& points,
                                const std::vector<std::vector<int>>& simplices)
    {
        for (const auto& simplex : simplices)
        {
            if (is_point_in_triangle(point,
                                    points[simplex[0]],
                                    points[simplex[1]],
                                    points[simplex[2]]))
            {
                return true;
            }
        }
        return false;
    }

    std::vector<double> eval_torch_model(const std::vector<double>& inputs)
    {
        if (!model_loaded)
        {
            try
            {
                std::cout << "Attempting to load model..." << std::endl;
                model = torch::jit::load("scripted_model_biomass_h_muconate.pt");
                std::cout << "Model Loaded..." << std::endl;
                model.eval();
                model_loaded = true;
            }
            catch (const c10::Error& e)
            {
                std::cerr << "Error loading the PyTorch model: " << e.what() << std::endl;
                return {};
            }
        }

        if (is_point_in_convex_hull(inputs, tri_points, tri_simplices))
        {
            std::vector<double> scaled_input = scale_input(inputs, x_mean, x_scale);
            torch::Tensor input_tensor = torch::tensor(scaled_input).unsqueeze(0).to(torch::kFloat32);
            torch::Tensor output_tensor = model.forward({input_tensor}).toTensor();

            std::vector<double> raw_output = {
                output_tensor[0][0].item<double>(),
                output_tensor[0][1].item<double>(),
                output_tensor[0][2].item<double>()
            };

            std::vector<double> outputs = unscale_output(raw_output, y_mean, y_scale);

        return outputs;
        } else {
            std::cout << "Point is outside the convex hull." << std::endl;
            std::cout << "Inputs: " << inputs << std::endl;
            return std::vector<double>(3, 0.0);

        }
    }

    void get_sp_id(std::string name, std::vector<int>& id_map, int foam_id)
    {
      auto it = sp_keys.find(name);
      if(it != sp_keys.end())
	{
	  id_map[it->second] = foam_id;
	}
    }


    void get_sp_id(std::string name, std::vector<int>& id_map, int foam_id)
    {
      auto it = sp_keys.find(name);
      if(it != sp_keys.end())
	{
	  id_map[it->second] = foam_id;
	}
    }

  
    void get_rhs(std::vector<double>& rhs, std::vector<double> solnvec, double t, int nvars, std::vector<double>& MM_params)
    {
 
        double K_s = MM_params[MM_param_keys["K_s"]];
        double K_o = MM_params[MM_param_keys["K_o"]];
        double q_max = MM_params[MM_param_keys["q_max"]];
        double x_max = MM_params[MM_param_keys["x_max"]];
        double y_xs = MM_params[MM_param_keys["y_xs"]];
        double y_ms = MM_params[MM_param_keys["y_ms"]];
        double y_os = MM_params[MM_param_keys["y_os"]];

        // calculate q_s
        double F_s = solnvec[G]/(solnvec[G] + K_s);
        double F_e = solnvec[O2]/(solnvec[O2] + K_o);
        double q_s = q_max*F_s*F_e;

        // calculate final rates
        rhs[X] = y_xs*q_s*solnvec[X]*(1 - solnvec[X]/x_max);

        //set RHS of o2 to 0 as it is solved in CFD
        rhs[O2] = 0.0;
        rhs[G] = -q_s*solnvec[X];
        rhs[M] = y_ms*q_s*solnvec[X];
	// rhs[CO2] = 0.0;

    }

    void get_rhs_ml(std::vector<double>& rhs, std::vector<double> solnvec, double t, int nvars, std::vector<double>& MM_params)
    {
        double K_s = MM_params[MM_param_keys["K_s"]];
        double K_o = MM_params[MM_param_keys["K_o"]];
        double q_max = MM_params[MM_param_keys["q_max"]];
        double x_max = MM_params[MM_param_keys["x_max"]];
        double y_xs = MM_params[MM_param_keys["y_xs"]];
        double y_ms = MM_params[MM_param_keys["y_ms"]];
        double y_os = MM_params[MM_param_keys["y_os"]];

	    // calculate q_s
        double F_s = solnvec[G]/(solnvec[G] + K_s);
        double F_e = solnvec[O2]/(solnvec[O2] + K_o);
        double q_s = q_max*F_s*F_e; // mol/(kg*hr)

        double rglu = -q_s*solnvec[X]; // mol / (m^3 * hr)
        double our = y_os * q_s * solnvec[X];

        // Prepare input for the model
        std::vector<double> inputs = {rglu / solnvec[X], -our / solnvec[X]}; // mol/(kg*hr)
        // Evaluate model
        std::vector<double> outputs = eval_torch_model(inputs);
        // Extract the outputs
        double mu_bio = outputs[0];
        double r_muc = outputs[2];
        double rbio_ml = mu_bio * solnvec[X];
        double rmuc_ml = r_muc * solnvec[X];

        // Calculate final rates
        rhs[X] = rbio_ml;
        rhs[O2] = 0.0;
        rhs[G] = rglu;
        rhs[M] = rmuc_ml;
        rhs[CO2] = 0.0;

    }

    void advance(std::vector<double>& solnvec, int nvars, double t_now, double t_adv, double dt, std::vector<double>& MM_params, bool do_ml)
    {
        double current_time=t_now;
        double final_time=t_now+t_adv;

        std::vector<double> rhs(nvars);
        std::vector<double> solnvec_n(nvars);

        while(current_time < final_time)
        {
            current_time += dt;

            //at current time level n
            solnvec_n=solnvec;

            //Doing RK23

            //stage 1
            if (do_ml) {
                get_rhs_ml(rhs, solnvec, current_time, nvars, MM_params);
            } else {
                get_rhs(rhs, solnvec, current_time, nvars, MM_params);
            }
            for(int i=0;i<nvars;i++)
            {
                solnvec[i] = solnvec_n[i] + 0.5*rhs[i]*dt;
            }

            //stage 2
            if (do_ml) {
                get_rhs_ml(rhs, solnvec, current_time, nvars, MM_params);
            } else {
                get_rhs(rhs, solnvec, current_time, nvars, MM_params);
            }
            for(int i=0;i<nvars;i++)
            {
                solnvec[i] = solnvec_n[i] + rhs[i]*dt;
            }
        }
    }

    double get_our(std::vector<double> solnvec,int nvars, std::vector<double>& MM_params)
    {
        double K_s = MM_params[MM_param_keys["K_s"]];
        double K_o = MM_params[MM_param_keys["K_o"]];
        double q_max = MM_params[MM_param_keys["q_max"]];
        double y_os = MM_params[MM_param_keys["y_os"]];

        // calculate q_s
        double F_s = solnvec[G]/(solnvec[G] + K_s);
        double F_e = solnvec[O2]/(solnvec[O2] + K_o);
        double q_s = q_max*F_s*F_e;

        double our = y_os*q_s*solnvec[X];

        return our;
    }
}
