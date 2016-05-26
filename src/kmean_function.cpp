#include"peak.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<tuple>
#include<map>
#include<cmath>
#include<iterator>
#include<cstdlib>
#include<R.h>
#include<Rdefines.h>
#include<Rmath.h>
#include<Rcpp.h>


typedef std::pair<int, double> shift_dist;

peak shift(const peak &_start, const int t);
void intersection(const peak &v1,const peak &v2, peak &_out);
void distance_L2_intersected(const peak &v1, const peak &v2, std::vector<double> &d, int order_Lp, char normalize);
void unione(const peak &v1,const peak &v2, peak &_out);
void distance_L2_union(const peak &v1, const peak &v2, const double &all_dato, 
                       const double &all_deriv, std::vector<double> &d, int order_Lp, char normalize);
std::pair<int, double> optim_shift (const peak &_first, const peak &_second, int LB, int UB, const char &tipo_dominio, const  double &all,
                                    const  double &all_der, const double & aplha, const double & peso, int order_Lp, char normalize);
void one_iteration(const std::vector<peak> &dati, const std::vector<peak> &templates, 
                   std::vector<int> &coef_shift, std::vector<double> & simil, std::vector<int> & template_assigned,
                   const char &tipo_dom, const double all, const double all_der, int & LB, int & UB,
                   const double & alpha, const double & peso, int order_Lp, char normalize);
void distingui_shift_dist(const std::vector<shift_dist> &coppie, std::vector<int> & coef_shift, std::vector<double> & dist);
void subset_picchi(const std::vector<peak> & dati, std::vector<peak> & extract, std::vector<int> &chosen);
std::map<int, std::vector<int>> distinguish_cluster(const std::vector<int> & clusters, std::vector<int> & num_elem_cluster);
void find_template(const std::vector<peak> & dati, const std::vector<int> & cluster, std::vector<peak> & templates, 
                   std::vector<double> & similarity_total,
                   const char & tipo_dominio, const double &all_dato, const double &all_deriv, const double & alpha,
                   const double & peso, int order_Lp, char normalize);
void shift_all(std::vector<peak> & dati, const std::vector<int> & coef_shift);
void kma_discrete(std::vector<peak> & dati, const int & n_clust, std::vector<int> & seeds, const double & t_max, 
                  const double & toll, const int & iter_max, const char & domain, const double & allung_dato, const double & allung_der,
                  std::vector<int> & labels, std::vector<double> & distances, std::vector<int> & shift, const double & alpha, 
                  const double & peso, int order_Lp, char normalize, const char verbose);  
void normalize_data(std::vector<int> & shift, std::vector<int>  & cluster, const int & n_clust);    
std::vector<double> global_distance(std::vector<peak> & dati,const double & alpha, const double & peso, const char & domain, 
                                    const double & allung_dato, const double & allung_der, int order_Lp, char normalize);


extern "C"{
  
  SEXP kmean_function( SEXP x, SEXP spline, SEXP spline_der, 
                           SEXP lenght, SEXP seeds, SEXP align, 
                           SEXP k, SEXP weight_input, SEXP alpha_input,
                           SEXP p_input, SEXP t_max_input, SEXP verbose)
  {
    
    // check in R that the dimentions of x, spline and spline_der are coherent
    // All the matrixes must have the same number of rows (num_data)
    // and columns (num_points)
    
    char norm= 'F'; //the distance will not be normalized with the length of the domain
    
    int num_data = INTEGER(GET_DIM(x))[0];
    int num_points = INTEGER(GET_DIM(x)) [1];
    
    Rcpp::NumericMatrix x_cpp(x);
    Rcpp::NumericMatrix spline_cpp(spline);
    Rcpp::NumericMatrix spline_der_cpp(spline_der);
    
    Rcpp::NumericVector length_cpp(lenght);
    
    std::vector<peak> dati(num_data);
    
    for (unsigned int i =0 ; i<num_data; i++){
      
      Rcpp::NumericVector tmp_x(num_points); 
      Rcpp::NumericVector tmp_spline(num_points);
      Rcpp::NumericVector tmp_spline_der(num_points);
      
      for (unsigned int t =0 ; t < num_points; t++)
      {
        tmp_x(t) = x_cpp(i,t);
        tmp_spline(t) = spline_cpp(i,t);
        tmp_spline_der(t) = spline_der_cpp(i,t);
      }
      
      std::vector<int> tmp_x_vect = Rcpp::as<std::vector<int>>(tmp_x);
      std::vector<double> tmp_spline_vect = Rcpp::as<std::vector<double>>(tmp_spline);
      std::vector<double> tmp_spline_der_vect = Rcpp::as<std::vector<double>>(tmp_spline_der);
      
      tmp_x_vect.resize(length_cpp(i));
      tmp_spline_vect.resize(length_cpp(i));
      tmp_spline_der_vect.resize(length_cpp(i));
      
      dati[i] = peak(tmp_x_vect, tmp_spline_vect, tmp_spline_der_vect);
    }
    
    char verb= Rcpp::as<char>(verbose);
    
    int numero_clusters = Rcpp::as<int>(k);
    
    if (verb == 'T')
    {
      Rcpp::Rcout<<"Definition of k="<<numero_clusters<<" clusters"<<'\n'; 
    }
    
    Rcpp::NumericVector seeds_cpp(seeds);
    std::vector<int> template_scelti = Rcpp::as<std::vector<int>>(seeds_cpp);
    
    
    std::vector<peak> templates(numero_clusters);
    subset_picchi(dati, templates, template_scelti);
    
    std::vector<int> coef_tra(dati.size());
    std::vector<double> dist_tra(dati.size());
    std::vector<int> clusters_new(dati.size());
    
    double all=0;
    double all_der=0;
    
    char alignment= Rcpp::as<char>(align);
  
    
    double t_max;
    if (alignment == 'T')
    {
      t_max= Rcpp::as<double>(t_max_input);
    }
    if (alignment == 'F')
    {
      t_max=0;
    }
    
    double toll=0.01;
    int iter_max=10;
    
    double peso = Rcpp::as<double>(weight_input); 
    
    char tipo_dom='U';
    
    double alpha=Rcpp::as<double>(alpha_input);
    
    if (verb == 'T')
    {
      Rcpp::Rcout<<"With weight "<<peso<<" you have fixed the value of alpha for the distance : ((1-alpha) * d0L2 + alpha * weight * d1L2)  to:"<<alpha<<'\n';
    }
    
    
    int order_p= Rcpp::as<int>(p_input);
    
    
    if (verb == 'T')
    {
      Rcpp::Rcout<<"The value of p for the order of the metric is "<<order_p<<'\n';
    } 
   
    
    std::vector<int> labels(num_data);
    std::vector<double> distances(num_data);
    std::vector<int> shift(num_data);
    kma_discrete(dati, numero_clusters, template_scelti, t_max, toll,  iter_max, 
                 tipo_dom, all, all_der, labels, distances, shift, alpha, peso, order_p, norm, verb);
    
    
    // Return the output in a list
    
    Rcpp::List res;
    
    res["labels"] = labels;
    
    
    res["shift"] = shift;
    
    res["distances"] = distances;
    
    
    return res;
  }
  
  
  
  
}




extern "C"{
    
    SEXP distance_matrix( SEXP x, SEXP spline, SEXP spline_der,
                          SEXP len, 
                          SEXP p_input)
    {
        int num_data = INTEGER(GET_DIM(x))[0];
        int num_points = INTEGER(GET_DIM(x)) [1];
        
        int order_p= Rcpp::as<int>(p_input);
        
        double all=0; //value used to extend the data
        double all_der=0; //value used to extend the derivative
        
        Rcpp::NumericMatrix x_cpp(x);
        Rcpp::NumericMatrix spline_cpp(spline);
        Rcpp::NumericMatrix spline_der_cpp(spline_der);
        
        Rcpp::NumericVector length_cpp(len);
        
        std::vector<peak> dati(num_data);
        
        for (unsigned int i =0 ; i<num_data; i++){
            
            Rcpp::NumericVector tmp_x(num_points); 
            Rcpp::NumericVector tmp_spline(num_points);
            Rcpp::NumericVector tmp_spline_der(num_points);
            
            for (unsigned int t =0 ; t < num_points; t++)
            {
                tmp_x(t) = x_cpp(i,t);
                tmp_spline(t) = spline_cpp(i,t);
                tmp_spline_der(t) = spline_der_cpp(i,t);
            }
            
            std::vector<int> tmp_x_vect = Rcpp::as<std::vector<int>>(tmp_x);
            std::vector<double> tmp_spline_vect = Rcpp::as<std::vector<double>>(tmp_spline);
            std::vector<double> tmp_spline_der_vect = Rcpp::as<std::vector<double>>(tmp_spline_der);
            
            int elem = length_cpp(i);
            
            tmp_x_vect.resize(elem);
            tmp_spline_vect.resize(elem);
            tmp_spline_der_vect.resize(elem);
            
            dati[i] = peak(tmp_x_vect, tmp_spline_vect, tmp_spline_der_vect);
            
        }
        
        char tipo_dominio = 'U';
        
        char normalize ='F';
        
        std::vector < std::vector<double> > D0(dati.size(), std::vector<double>(dati.size())); 
        std::vector < std::vector<double> > D1(dati.size(), std::vector<double>(dati.size())); 
        
        for (unsigned int k=0; k<dati.size(); k++)
        {
            
            for (unsigned int p=k ; p<dati.size(); p++)
            {
                
                std::vector<double> dist(2);
                if (tipo_dominio == 'I') 
                {
                    distance_L2_intersected(dati[k], dati[p], dist, order_p, normalize);
                }
                if (tipo_dominio == 'U')
                {
                    distance_L2_union(dati[k], dati[p], all, all_der, dist, order_p, normalize);
                }
                
                D0[p][k]= dist[0];
                D0[k][p]=D0[p][k];
                
                D1[p][k]= dist[1];
                D1[k][p]=D1[p][k];
                
            }
        }
        
        Rcpp::List res;
        
        res["dist_D0"] = D0;
        res["dist_D1"] = D1; 
        
        return res;
    }
    
}


// computes the intersection between peaks. The result is the update of _out.
// ABSCISSA : the common abcissa beween v1 and v2 peaks
// DATA and DERIVATIVE : of the FIRST peak (v1)

void intersection(const peak &v1,const peak &v2, peak &_out)
{
    
    int primo_v1, primo_v2, ultimo_v1, ultimo_v2;
    primo_v1= v1.getRange()[0];
    ultimo_v1= v1.getRange()[1];
    primo_v2= v2.getRange()[0];
    ultimo_v2= v2.getRange()[1];
    if (ultimo_v1 >= primo_v2 and ultimo_v2 >= primo_v1)  
    {
        set_intersection( std::begin(v1.getData()),std::end(v1.getData()),
                          std::begin(v2.getData()),std::end(v2.getData()),
                          std::back_inserter(_out.getData()), 
                          [](const peak::data_Type & first, const peak::data_Type & second) {                                          
                              return std::get<0>(first) < std::get<0>(second);});
    }
    
    
}

// computes the distance (in d) between ordinate and derivative of v1 and v2 in the overlapping part

void distance_L2_intersected(const peak &v1, const peak &v2, std::vector<double> &d, int order_Lp, char normalize)
{
    peak intersect_1;
    intersection(v1, v2, intersect_1);
    if (intersect_1.getData().size()!=0)
    {
        peak intersect_2;
        intersection(v2, intersect_1, intersect_2);
        peak differenza = intersect_1 - intersect_2;
        d=differenza.area(order_Lp, normalize);   
    }
    else{
        d[0] = 1000000;
        d[1] = 1000000;
    }
}

void unione(const peak &v1,const peak &v2, peak &_out)
{
    
    int primo_v1, primo_v2, ultimo_v1, ultimo_v2;
    primo_v1= v1.getRange()[0];
    ultimo_v1= v1.getRange()[1];
    primo_v2= v2.getRange()[0];
    ultimo_v2= v2.getRange()[1];
    if (ultimo_v1 >= primo_v2 and ultimo_v2 >= primo_v1)
    {
        std::set_union(std::begin(v1.getData()),std::end(v1.getData()),
                       std::begin(v2.getData()),std::end(v2.getData()),
                       std::back_inserter(_out.getData()), 
                       [](const peak::data_Type & first, const peak::data_Type & second) {                                          
                           return std::get<0>(first) < std::get<0>(second);});
    }
    
    
    
}


void distance_L2_union(const peak &v1, const peak &v2, const double &all_dato, const double &all_deriv, std::vector<double> &d, int order_Lp, char normalize)
{
    peak union_1;
    unione(v1, v2, union_1);
    if (union_1.getData().size()!=0)
    {
        peak union_2(union_1.getAbscissa(), all_dato, all_deriv);
        peak v1_long, v2_long;
        unione(v1, union_2, v1_long);
        unione(v2, union_2, v2_long);
        
        peak differenza= v1_long - v2_long;
        
        d=differenza.area(order_Lp, normalize); 
    }
    else{
        d[0] = 1000000;
        d[1] = 1000000;
    }
}







// shift of a peak of a INTEGER value t

peak shift(const peak &_start, const int _t){
  
  peak copia(_start);
  
  for (auto & i : copia.getData())
  {
    std::get<0>(i) += _t;
  }
  
  return copia;
}

// optimization function: 
// detects the best shift between LB and UB and returns a pair with as first element the coefficient of the shift
// and as second the distance

std::pair<int, double> optim_shift (const peak &_first, const peak &_second, int LB, int UB, const char &tipo_dominio, const double & all, const double & all_der, const double & alpha, const double & peso, int order_Lp, char normalize)
{
  std::vector<double> distance(UB-LB+1);
  int position=0;
  for (int t = LB; t<=UB; t++)
  {
    peak traslato=shift(_second, t);
    std::vector<double> dist(2);
    if (tipo_dominio == 'I') 
    {
      distance_L2_intersected(_first, traslato, dist, order_Lp, normalize);
    }
    if (tipo_dominio == 'U')
    {
      distance_L2_union(_first, traslato, all, all_der, dist, order_Lp, normalize); //all and all_der are the values used to extend the data and the derivative
    }
    distance[position]= (1-alpha) * dist[0]+ alpha * peso * dist[1];
   
    position++; //in the end position contains the lenght: UB-LB+1
  }
  int pto_min=LB;
  double min=distance[0];
  for (int i=1; i<position; i++)
  {
    if (distance[i] < min){
      min=distance[i];
      pto_min=LB+i;
    }
    if ((distance[i] ==  min) and (std::fabs(LB+i) < std::fabs(pto_min))) {
      min=distance[i];
      pto_min=LB+i;
    }
    
  }
  std::pair<int, double> OUT (pto_min, min);
  return(OUT);
}


void one_iteration(const std::vector<peak> &dati, const std::vector<peak> &templates, 
                   std::vector<int> &coef_shift, std::vector<double> & simil, std::vector<int> & template_assigned, //template_assigned will contain the number of clusters
                   const char &tipo_dom, const double all, const double all_der, int & LB, int & UB, const double & alpha, const double & peso, int order_Lp, char normalize)
{
  std::vector < std::vector<shift_dist> > info(dati.size(), std::vector<shift_dist>(templates.size())); 
  for (unsigned int k=0; k<templates.size(); k++)
  {
    
    for (unsigned int p=0 ; p<dati.size(); p++)
    {
       info[p][k]= optim_shift( templates[k], dati[p], LB, UB , tipo_dom, all, all_der, alpha, peso, order_Lp, normalize);  
      //the second element is shifted. It is NOT the template.

    }
  }
  
  std::vector<shift_dist> assegnamento(dati.size());
  for (unsigned int i=0; i<dati.size(); i++)
  {
    int template_min_dist=0;
    double min_dist=std::get<1>(info[i][0]);
    for (unsigned int j=1; j<templates.size(); j++)
    {
      if (std::get<1>(info[i][j])<min_dist)
      {
        min_dist=std::get<1>(info[i][j]);
        template_min_dist=j;
      }
    }
    template_assigned[i]=template_min_dist;
    assegnamento[i]=info[i][template_min_dist];
  }
  distingui_shift_dist(assegnamento, coef_shift, simil);
  
}


void distingui_shift_dist(const std::vector<shift_dist> &coppie, std::vector<int> & coef_shift, std::vector<double> & dist)
{
  for (unsigned int i=0; i<coef_shift.size();i++)
  {
    coef_shift[i]=std::get<0>(coppie[i]);
    dist[i]=std::get<1>(coppie[i]);
  }
}

// indexes listed in chosen. extract will contain dati[chosen]
void subset_picchi(const std::vector<peak> & dati, std::vector<peak> & extract, std::vector<int> &chosen)
{
  for ( unsigned int k=0; k<chosen.size(); k++)
  {
    extract[k] = dati[chosen[k]];
  }
}

std::map<int, std::vector<int>> distinguish_cluster(const std::vector<int> & clusters, std::vector<int> & num_elem_cluster)
{
  std::map<int, std::vector<int>> lista_cluster;
  for (unsigned int k=0; k<num_elem_cluster.size(); k++)
  {
    num_elem_cluster[k]=0;
  }
  for (unsigned int i=0; i<clusters.size(); i++)
  {
    lista_cluster[clusters[i]].push_back(i);
    num_elem_cluster[clusters[i]]++;
  }
  return lista_cluster;
}


// from dati (vector of all the 'peak' elements) and a vector of subdivisions of dati, find the best tempaltes
void find_template(const std::vector<peak> & dati, const std::vector<int> & cluster, std::vector<peak> & templates, //similarity_total will contain the overall distance of the km clusters.
                   std::vector<double> & similarity_total,  const char & tipo_dominio,const double &all_dato,
                   const double &all_deriv, const double & alpha, const double & peso, int order_Lp, char normalize)
{
  std::vector<int> number_elem(templates.size());
  std::map<int, std::vector<int>> elementi_cluster;
  elementi_cluster=distinguish_cluster(cluster, number_elem);
  
  std::vector<int> templates_new(templates.size());
  
  int cluster_vuoti=0;
  
  for (unsigned int k=0; k<templates.size(); k++)
  {
    std::vector<double> dist_totali(number_elem[k],0);
    
    //     std::cout<<"cluster "<<k<<"con elementi"<<number_elem[k]<<'\n';
    if (number_elem[k]!=0)
    {
      cluster_vuoti ++;
      for (int t=0; t<number_elem[k]; t++)
      {
        for (int q=t+1; q<number_elem[k]; q++)
        {
          if ( tipo_dominio == 'I')
          {
            std::vector<double> dist(2);
            distance_L2_intersected(dati[elementi_cluster[k][t]], dati[elementi_cluster[k][q]], dist,
                                    order_Lp, normalize );
            dist_totali[t] += (1-alpha) * dist[0] + alpha * peso * dist[1];
            dist_totali[q] += (1-alpha) * dist[0] + alpha * peso * dist[1];
          }
          if ( tipo_dominio == 'U')
          {
            std::vector<double> dist(2);
            distance_L2_union(dati[elementi_cluster[k][t]], dati[elementi_cluster[k][q]], all_dato, all_deriv, 
                              dist, order_Lp, normalize);
            dist_totali[t] += (1-alpha) * dist[0] + alpha * peso * dist[1];
            dist_totali[q] += (1-alpha) * dist[0] + alpha * peso * dist[1];
          }
        }
      }
      
      
      
      templates_new[k]=elementi_cluster[k][0];
      similarity_total[k]=1000000;
      
      for (int t=0; t<number_elem[k]; t++)
      {
        if (dist_totali[t]<=similarity_total[k])
        {
          templates_new[k] = elementi_cluster[k][t];
          similarity_total[k] = dist_totali[t];
        }
      }
    }
  }
  subset_picchi(dati, templates, templates_new);
}

void shift_all(std::vector<peak> & dati, const std::vector<int> & coef_shift)
{
  for (unsigned int i=0 ; i<dati.size(); i++)
  {
    peak copia;
    copia = shift(dati[i], coef_shift[i]);
    dati[i] = copia;
  }
}

// it performs k-medoids algorithm 
// warping.method = shift
// similarity method d0.L2 (but easily generaliz with d0ed1.L2)
void kma_discrete(std::vector<peak> & dati, const int & n_clust, std::vector<int> & seeds, const double & t_max, 
                  const double & toll, const int & iter_max, const char & domain, const double & allung_dato, const double & allung_der, //parameters of the model 
                  std::vector<int> & labels, std::vector<double> & distances, std::vector<int> & shift, const double &  alpha, 
                  const double & peso, int order_Lp, char normalize, const char verbose) 
{
  // original data
  std::vector<peak> dati_orig(dati.size());
  for (unsigned int i=0; i<dati.size(); i++)
  {
    dati_orig[i]=dati[i];
  }
  // defintion of the maximum shitf
  std::vector<int> range_1 = dati[0].getRange();
  double coeff_shift = t_max * (range_1[1]-range_1[0]);
  // cast truncate the double, and I want the floor so it is ok, if I want the right conversion shift = shift + 0.5; 
  // int trasl = (int) shift;  or the ceiling +1
  int trasl = (int) coeff_shift;
  std::vector<int> range_new(2);
  double shift_coeff_new;
  int trasl_new;
  for (unsigned int i=1; i<dati.size(); i++)
  {
    range_new = dati[i].getRange();
    shift_coeff_new = t_max * (range_new[1]-range_new[0]);
    trasl_new = (int) shift_coeff_new;
    if (trasl_new < trasl)
    {
      trasl = trasl_new;
    }
  }
  int LB = -trasl;
  int UB = trasl;
  
  if (verbose == 'T') 
  {
    Rcpp::Rcout<<"The value of t_max "<<t_max<<" corresponds to a maximum shift of "<< UB<<'\n';
  }
  
  if ( (unsigned int)n_clust != seeds.size())
  {
    Rcpp::Rcerr<<"Number of clusters "<<n_clust<<" differnt from the lenght of vector seeds: PROVIDE A CORRECT NUMBER OF SEEDS FOR INITIALIZATION"<<std::endl;
  }
  

  std::vector<peak> templates(n_clust);
  subset_picchi(dati, templates, seeds);
  
  // initialization
  int iter = 1;
  std::vector<double> distances_old(dati.size(), 0);
  std::vector<int> labels_old(dati.size(), 10);
  
  int number_distances_low=0; // will contain the number of data with a SIGNIFICANT(i.e <toll) decrease of the distance with the iteration
  int number_clusters_different=dati.size(); //will contain the number of data belonging to differnt clusters
  
  std::vector<double> distances_new(dati.size());
  std::vector<int> labels_new(dati.size());
  std::vector<int> shift_new(dati.size());
  
  std::vector<int> shift_global(dati.size(),0);
  
  int cluster_vuoti=0;
  
  while( iter < iter_max and number_distances_low < dati.size()  and cluster_vuoti==0){  //and number_clusters_different > 0
    
    if (verbose=='T')
    {
      Rcpp::Rcout<< "Iteration "<<iter<<'\n'<<"------------------------------------"<<'\n';
    }
    
    distances_new.clear();
    distances_new.resize(dati.size());
    labels_new.clear();
    labels_new.resize(dati.size());
    shift_new.clear();
    shift_new.resize(dati.size());
    
    number_distances_low=0;
    number_clusters_different=0;
    one_iteration(dati, templates, shift_new, distances_new, labels_new, domain, allung_dato, 
                  allung_der, LB, UB, alpha, peso, order_Lp, normalize);
    
    shift_all(dati, shift_new);
    for (unsigned int i=0; i< dati.size(); i++)
    {
      shift_global[i] += shift_new[i];
      if ( (distances_new[i] - distances_old[i]) < toll)
      {
        number_distances_low++;
      }
      if ( (labels_new[i] - labels_old[i]) != 0)
      {
        number_clusters_different ++;
      }
      distances_old[i]=distances_new[i];
      distances_new[i]=100000;
      labels_old[i]=labels_new[i];
      labels_new[i]=10;
    }
    
    if (verbose == 'T')
    {
      if ( (unsigned int)number_clusters_different == dati.size() & iter != 1)
      {
        Rcpp::Rcout<<"In iteration "<<iter<<" all the data change cluster: WARNING: possible loop"<<'\n';
      }
      if (number_distances_low== dati.size())
      {
        Rcpp::Rcout<<"number of data which decrease the similarity is dati.size()"<<'\n';
      }
      if ( (unsigned int)number_clusters_different == 0)
      {
        Rcpp::Rcout<<"no variation in the clusters"<<'\n';
      }
    }
    
    std::vector<double> similarity_total(n_clust);
    find_template(dati, labels_old, templates, similarity_total, domain, allung_dato, allung_der, alpha, peso, order_Lp, normalize);
   
    std::vector<int> number_elem_clusters(templates.size());
    for (unsigned int i=0; i<dati.size(); i++)
    {
      number_elem_clusters[labels_old[i]]++;
    }
    for (unsigned int k=0; k<templates.size(); k++)
    {
      if (verbose== 'T')
      {
        Rcpp::Rcout<<"Number elementes cluster "<<k<<" : "<<number_elem_clusters[k]<<'\n'; 
      }
      if (number_elem_clusters[k]==0)
      {
        cluster_vuoti++;
      }  
    }
    if (verbose =='T')
    {
      if (cluster_vuoti!=0) 
      {
        Rcpp::Rcout<<"empty cluster detected => algorithm stops: try to decrease number of clusters or change seeds"<<'\n';
      }
    }
   
    iter++;
  }
  for (unsigned int i=0; i<dati.size(); i++)
  {
    labels[i]=labels_old[i];
    shift[i]=shift_global[i];
  }
  
  normalize_data(shift, labels, n_clust);
  
  for (unsigned int i=0; i<dati.size(); i++)
  {
    std::vector<double> dist(2);
    if (domain=='I')
    {
      distance_L2_intersected(dati[i], templates[labels[i]], dist , order_Lp, normalize);
    }
    if (domain == 'U')
    {
      distance_L2_union(dati[i], templates[labels[i]], allung_dato, allung_der, dist, order_Lp, normalize);
    }
    distances[i]= (1-alpha) * dist[0] + alpha * peso * dist[1];
    
  }
  shift_all(dati_orig,shift); //now I have the correct tanslation of the data and I store it in dati.
  for (unsigned int i=0; i<dati.size(); i++)
  {
    dati[i]=dati_orig[i];
  }
}

void normalize_data(std::vector<int> & shift, std::vector<int>  & cluster, const int & n_clust) // make shift = shift - mean(shift, k) with mean(shift, k) the mean of the shift of cluster k
{
  std::vector<int> sum_shift_clusters(n_clust, 0);
  std::vector<int> number_elements_clusters(n_clust, 0);
  for(unsigned int i=0; i<shift.size(); i++)
  {
    sum_shift_clusters[cluster[i]] += shift[i];
    number_elements_clusters[cluster[i]]++;
  }
  std::vector<int> average_shift(n_clust,0);
  for (unsigned int i=0; i<n_clust; i++)
  {
    if (number_elements_clusters[i]!=0)
    {
      double average;
      //std::cout<<"numero elementi cluster "<<i<<": "<<number_elements_clusters[i]<<'\n';
      average = sum_shift_clusters[i] / number_elements_clusters[i] + 0.5; //+0.5 to have the correct round and not only cut off with the cast
      average_shift[i]= (int) average;
    }  
  }
  for (unsigned int i = 0; i<shift.size(); i++)
  {
    shift[i] = shift[i]-average_shift[cluster[i]];
  }
}

std::vector<double> global_distance(std::vector<peak> & dati,const double & alpha, const double & peso, const char & domain, const double & allung_dato, const double & allung_der, int order_Lp, char normalize)
{
  std::vector<double> dist_tot(dati.size(),0);
  for (unsigned int i =0 ; i<dati.size(); i++)
  {
    for (unsigned int j = i+1; j<dati.size(); j++)
    {
      std::vector<double> dist(2);
      if (domain=='I')
      {
        distance_L2_intersected(dati[i], dati[j], dist, order_Lp, normalize );
      }
      if (domain == 'U')
      {
        distance_L2_union(dati[i], dati[j], allung_dato, allung_der, dist, order_Lp, normalize);
      }
      dist_tot[i] += (1-alpha) * dist[0] + alpha * peso * dist[1];
      dist_tot[j] += (1-alpha) * dist[0] + alpha * peso * dist[1];
    }
  }
  return dist_tot;
}
