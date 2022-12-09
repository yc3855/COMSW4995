#ifndef WAVESOLVER_H_
#define WAVESOLVER_H_
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <fstream>
#include <iostream>
#include "CustomLibraries/np.hpp"
#include "CustomLibraries/np_to_matplot.hpp"

/*!
 *  \addtogroup wavePlotter
 *  @{
 */

//! Custom plotter class
namespace wavePlotter
{
    //! Plotter class
    //! \brief This class is used to plot the wave field
    //! TODO: make it multithreaded
    class Plotter
    {
    public:
        //! Constructor
        Plotter(const boost::multi_array<double, 3> &u, const matplot::vector_2d &Xp, const matplot::vector_2d &Zp, int num_levels, int nt)
        {
            this->u.resize(boost::extents[u.shape()[0]][u.shape()[1]][u.shape()[2]]);
            this->u = u;
            this->Xp = Xp;
            this->Zp = Zp;
            this->num_levels = num_levels;
            this->nt = nt;
            double min_u = np::min(u);
            double max_u = np::max(u);
            std::cout << "min_u = " << min_u << " max_u = " << max_u << "\n";
            this->levels = matplot::linspace(min_u, max_u, num_levels);
        }
        //! Renders a frame of the wave field to a image on disk
        void renderFrame(int index)
        {
            matplot::vector_2d Up = np::convert_to_matplot(this->u[index]);
            matplot::contourf(this->Xp, this->Zp, Up, this->levels);
            matplot::save(save_directory + "/contourf_" + format_num(index) + ".png");
        }

        //! Renders all frames of the wave field to form an animation to be saved on disk
        void renderAllFrames(int begin_frame_index, int end_frame_index)
        {
            for (int i = begin_frame_index; i < end_frame_index; i++)
            {
                renderFrame(i);
            }
        }

        //! Renders a complete video animation of the wave field
        void animate(std::string output_file_name, int begin_frame_index, int end_frame_index, int frame_rate)
        {
            renderAllFrames(begin_frame_index, end_frame_index);
            std::string ffmpeg_render_command = "ffmpeg -framerate " + std::to_string(frame_rate) + " -pattern_type glob -i '" + save_directory + "/*.png' -c:v libx264 -pix_fmt yuv420p " + output_file_name;
            std::system(ffmpeg_render_command.c_str());
        }
        //! Export a frame of the wave field to a .csv format for external use
        void exportFrame(int index)
        {
            matplot::vector_2d Up = np::convert_to_matplot(this->u[index]);
            std::ofstream outfile;
            outfile.open(save_directory + "/frame_" + format_num(index) + ".csv");
            for (std::size_t i = 0; i < Up.size(); i++)
            {
                for (std::size_t j = 0; j < Up[i].size(); j++)
                {
                    outfile << Up[i][j];
                    if (j != Up[i].size() - 1)
                        outfile << ",";
                }
                outfile << "\n";
            }
            outfile.close();
        }

        //! Export all frames of the wave field to a .csv format for external use
        void exportAllFrames(int begin_frame_index, int end_frame_index)
        {
            for (int i = begin_frame_index; i < end_frame_index; i++)
            {
                exportFrame(i);
            }
        }

        //! Set the save directory for the rendered frames
        void setSaveDirectory(std::string save_directory)
        {
            this->save_directory = save_directory;
        }

    private:
        std::string format_num(int num, int length = 8)
        {
            std::string str_num = std::to_string(num);

            int str_length = str_num.length();
            for (int i = 0; i < length - str_length; i++)
                str_num = "0" + str_num;
            return str_num;
        }

        boost::multi_array<double, 3> u;
        matplot::vector_2d Xp;
        matplot::vector_2d Zp;
        int num_levels;
        int nt;
        std::vector<double> levels;
        std::string save_directory = "output";
    };

}
#endif