#ifndef __PEAK__
#define __PEAK__

//#include <Rcpp.h>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/max_cardinality_matching.hpp>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<tuple>
#include<cstdlib>
#include<cmath>


// [[Rcpp::depends(BH)]]

//using namespace Rcpp;
typedef std::basic_string<char> string;

class peak{

public:
	
		typedef std::tuple<int, double, double> data_Type;
		
		//Default constructor		
		peak()=default;
		
		//Costruttore con valori e lunghezza
		peak(const std::vector<int> & _abscissa, const std::vector<double> & _value, const std::vector<double> & _deriv){
			M_data.resize(_abscissa.size());
			for(unsigned int i=0; i<M_data.size(); i++)
			{
				M_data[i]=std::make_tuple(_abscissa[i], _value[i], _deriv[i]);
			}
	
		};
	
		
		peak(const std::vector<int> &_abscissa, const double &y, const double &y_1){
			M_data.resize(_abscissa.size());
			for(unsigned int i=0; i<M_data.size(); i++)
			{
				M_data[i]=std::make_tuple(_abscissa[i], y, y_1);
			}
		}
				
		//Copy Constructor
		peak(const peak& _p): M_data(_p.M_data)
			{ };
		
		//Copy Assignment
		peak & operator= (const peak & _p){
		
			if(this!=& _p){  //beware of self-assignment
				
				M_data.resize(_p.getData().size());
			
				for(unsigned int i=0; i<M_data.size(); i++)
				{
					M_data[i] = _p.M_data[i];
				}
								
			}
			
			return *this;
		};	
		
		//Distructor
		~peak() = default;
		//Peak features
		
		//Length of the peak
		unsigned  int size() const{
			return M_data.size();
		}
		
		std::vector<int> getRange() const{
			
			std::vector<int> range(2);
			range[0] = std::get<0>(M_data[0]);
			range[1] = std::get<0>(M_data[M_data.size()-1]);
			return range;
		
		}
		
		const std::vector<data_Type> & getData( ) const{return M_data;}
		
		std::vector<data_Type> & getData( ) {return M_data;}
		
		std::vector<int> getAbscissa() const {
			
			std::vector<int> x(M_data.size());
			
			for(unsigned int i=0; i<M_data.size(); i++ )
			{
				x[i] = std::get<0>(M_data[i]);
			}
			return x;
		}
		
		std::vector<double> getY() const {
			
			std::vector<double> x(M_data.size());
			
			for(unsigned int i=0; i<M_data.size(); i++ )
			{
				x[i] = std::get<1>(M_data[i]);
			}
			return x;
		}
		
		std::vector<double> getDerivative() const {
			
			std::vector<double> x(M_data.size());
			
			for(unsigned int i=0; i<M_data.size(); i++ )
			{
				x[i] = std::get<2>(M_data[i]);
			}
			return x;
		}
		
		std::vector<double> area(int order_Lp, char normalize) const{
			
			std::vector<double> area(2,0);
            std::vector<double> area_def(2,0);
			
            int D;
            if (normalize=='T')
            {
                D = std::get<0>(M_data[M_data.size()-1]) - std::get<0>(M_data[0]) + 1;
            }
            if (normalize=='F')
            {
                D = 1;
            }
            
			if (order_Lp==2) 
            {
                for (unsigned int i=0; i<M_data.size(); i++) // oppure fino a M_data.size()-1
			    {
			        area[0] += 2*(std::get<1>(M_data[i]))*(std::get<1>(M_data[i]));
			        area[1] += 2*(std::get<2>(M_data[i]))*(std::get<2>(M_data[i]));
			    }
			
                area[0] = area[0] - 
                    (std::get<1>(M_data[0]))*(std::get<1>(M_data[0])) - 
                        (std::get<1>(M_data[M_data.size()-1]))*(std::get<1>(M_data[M_data.size()-1]));
               area[1] = area[1] - 
                    (std::get<2>(M_data[0]))*(std::get<2>(M_data[0])) - 
                         (std::get<2>(M_data[M_data.size()-1]))*(std::get<2>(M_data[M_data.size()-1]));  
                area_def[0] = sqrt(area[0]/2)/D;
                area_def[1] = sqrt(area[1]/2)/D;
			}
       
            if (order_Lp==1)
            {
                for (unsigned int i=0; i<M_data.size(); i++)
			    {
			        area[0] += std::fabs(std::get<1>(M_data[i]));
			        area[1] += std::fabs(std::get<2>(M_data[i]));
			    }
                
                area_def[0] = area[0]/D;
                area_def[1] = area[1]/D;
            }
            
            if (order_Lp==0)
            {
                for (unsigned int i=0; i<M_data.size(); i++)
			    {
			        if (std::fabs(std::get<1>(M_data[i])) > area[0]) 
                        {
                            area[0] = std::fabs(std::get<1>(M_data[i]));
                        }
    			    if (std::fabs(std::get<2>(M_data[i])) > area[1]) 
                        {
                            area[1] = std::fabs(std::get<2>(M_data[i]));
                        } 
			    }
                
                area_def[0] = area[0]/D;
                area_def[1] = area[1]/D;
            }
            return area_def;
			
		}
		
		//Overload <<
		
		friend std::ostream& operator<<(std::ostream& OUT, const peak& _p){
			
			OUT.precision(2);
			
			if (_p.M_data.size()==0)
			{
				OUT<<"EMPTY";
			}
			else{
				for(int i=0; i<5; i++ ) //_p.M_data.size()
				{
					OUT<<std::get<0>(_p.M_data[i])<<' ';
				}
			
				OUT<<'\n';
			
				for(int i=0; i<5; i++ ) //_p.M_data.size()
				{
					OUT<<std::get<1>(_p.M_data[i])<<' ';
				}
			
				OUT<<'\n';
			
				for(int i=0; i<5; i++ )//_p.M_data.size()
				{
					OUT<<std::get<2>(_p.M_data[i])<<' ';
				}		
				
			}
			
			return OUT;
		}
		
		// Overload -
		// for peaks WITH THE SAME ASCISSA!! 
		friend peak operator - (const peak & _p1, const peak & _p2 ){
			
			std::vector<int> x(_p1.M_data.size(),0);
			std::vector<double> y(_p1.M_data.size(),0);
			std::vector<double> y_1(_p1.M_data.size(),0);
			
			for (unsigned int i=0; i<_p1.M_data.size(); i++)
			{
				x[i]= std::get<0>(_p1.M_data[i]);
				y[i]= std::get<1>(_p1.M_data[i]) - std::get<1>(_p2.M_data[i]);
				y_1[i]= std::get<2>(_p1.M_data[i]) - std::get<2>(_p2.M_data[i]);
			}
			
			return peak(x, y, y_1);
		}
		
		
private:
		
	std::vector<data_Type> M_data;
		
		
		
		
};
#endif
