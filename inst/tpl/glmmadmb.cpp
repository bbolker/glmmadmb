  #include <admodel.h>
  #include <df1b2fun.h>
  ofstream ofs11("b1");
  ofstream ofs12("b2");
  ofstream ofs13("s1");
  ofstream ofs14("s2");
  //void add_slave_suffix(const adstring & s){}
  dvariable betaln(const prevariable& a,const prevariable& b )
  {
    return gammln(a)+gammln(b)-gammln(a+b);
  }
  dvariable ln_beta_density(double y,const prevariable & mu,
    const prevariable& phi)
  {
    dvariable omega=mu*phi;
    dvariable tau=phi-mu*phi;
    dvariable lb=betaln(omega,tau);
    dvariable d=(omega-1)*log(y)+(tau-1)*log(1.0-y)-lb;
    return d;
  }
  df1b2variable betaln(const df1b2variable& a,const df1b2variable& b )
  {
    return gammln(a)+gammln(b)-gammln(a+b);
  }
  df1b2variable ln_beta_density(double y,const df1b2variable & mu,
    const df1b2variable& phi)
  {
    df1b2variable omega=mu*phi;
    df1b2variable tau=phi-mu*phi;
    df1b2variable lb=betaln(omega,tau);
    df1b2variable d=(omega-1)*log(y)+(tau-1)*log(1.0-y)-lb;
    return d;
  }
  dvariable random_bound(const prevariable& u,double a)
  {
    if (fabs(value(u))<=a)
      return u;
    else if (value(u)>a)
    {
      dvariable y=u-a;
      return a+y/(1+y);
    }
    else if (value(u)<-a)
    {
      dvariable y=-a-u;
      return -a-y/(1+y);
    }
  }
  df1b2variable random_bound(const df1b2variable& u,double a)
  {
    if (fabs(value(u))<=a)
      return u;
    else if (value(u)>a)
    {
      df1b2variable y=u-a;
      return a+y/(1+y);
    }
    else if (value(u)<-a)
    {
      df1b2variable y=-a-u;
      return -a-y/(1+y);
    }
  }
  dvar_vector random_bound(const dvar_vector& v,double a)
  {
    int mmin=v.indexmin();
    int mmax=v.indexmax();
    dvar_vector tmp(mmin,mmax);
    for (int i=mmin;i<=mmax;i++)
    {
      tmp(i)=random_bound(v(i),a);
    }
    return tmp;
  }
  df1b2vector random_bound(const df1b2vector& v,double a)
  {
    int mmin=v.indexmin();
    int mmax=v.indexmax();
    df1b2vector tmp(mmin,mmax);
    for (int i=mmin;i<=mmax;i++)
    {
      tmp(i)=random_bound(v(i),a);
    }
    return tmp;
  }
#include <admodel.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <glmmadmb.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  p_y.allocate("p_y");
  y.allocate(1,n,1,p_y,"y");
  p.allocate("p");
  X.allocate(1,n,1,p,"X");
  M.allocate("M");
  q.allocate(1,M,"q");
  m.allocate(1,M,"m");
  ncolZ.allocate("ncolZ");
  Z.allocate(1,n,1,ncolZ,"Z");
  I.allocate(1,n,1,ncolZ,"I");
  cor_flag.allocate(1,M,"cor_flag");
  cor_block_start.allocate(1,M,"cor_block_start");
  cor_block_stop.allocate(1,M,"cor_block_stop");
  numb_cor_params.allocate("numb_cor_params");
  like_type_flag.allocate("like_type_flag");
  link_type_flag.allocate("link_type_flag");
  rlinkflag.allocate("rlinkflag");
  no_rand_flag.allocate("no_rand_flag");
  zi_flag.allocate("zi_flag");
  zi_kluge.allocate("zi_kluge");
  poiss_prob_bound.allocate("poiss_prob_bound");
  nbinom1_flag.allocate("nbinom1_flag");
  intermediate_maxfn.allocate("intermediate_maxfn");
  has_offset.allocate("has_offset");
  offset.allocate(1,n,"offset");
  rr.allocate(1,n,1,6);
  phi.allocate(1,p,1,p);
  int i,j;
  phi.initialize();
  ymax=log(15.0*max(y)+1);
  for (i=1;i<=p;i++)
  {
    phi(i,i)=1.0;
  }
  dmatrix trr=trans(rr);
  trr(6).fill_seqadd(1,1);
  rr=trans(trr);
  dmatrix TX(1,p,1,n);
  TX=trans(X);
  for (i=1;i<=p;i++)
  {
    double tmp=norm(TX(i));
    TX(i)/=tmp;
    phi(i)/=tmp;
    for (j=i+1;j<=p;j++)
    {
      double a=TX(j)*TX(i);
      TX(j)-=a*TX(i);
      phi(j)-=a*phi(i);
    }
  }
  X=trans(TX);
  sum_mq = 0;
  for (i=1;i<=M;i++)
    sum_mq += m(i)*q(i);
  ofstream ofs("phi.rep");
  for (i=1; i<=p; i++)
  {
    for (j=1; j<=p; j++)
    {
      ofs << phi(i,j) << " ";
    }
    ofs << endl;
  }
  ofs << endl;
}

void model_parameters::initializationfunction(void)
{
  tmpL.set_initial_value(1.0);
  tmpL1.set_initial_value(0.0);
  log_alpha.set_initial_value(1);
  pz.set_initial_value(.001);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  // BMB: FIXME: do we need this?  formerly disallowed for binomial (was like_type_flag 2);
  //     would be problematic for binary data but otherwise OK.  Should test in R code
  //  if(zi_flag && (like_type_flag>=2))
  // {
  //  cerr << "Zero inflation not allowed for this response type" << endl;
  //  ad_exit(1);
  // }
  // Determines "phases", i.e. when the various parameters  becomes active in the optimization process
  int pctr = 2;		// "Current phase" in the stagewise procedure
  // FIXME: move trunc_poisson, logistic earlier in like_type hierarchy (or add a scale parameter flag vector)
  int alpha_phase = like_type_flag>1 && like_type_flag!=6 ? pctr++ : -1;        // Phase 2 if active
  int zi_phase = zi_flag ? pctr++ : -1;                      			// After alpha
  int rand_phase = no_rand_flag==0 ? pctr++ : -1;    				// SD of RE's
  int cor_phase = (rand_phase>0) && (sum(cor_flag)>0) ? pctr++ : -1 ; 		// Correlations of RE's
  // Count the number of variance/correlation parameters to be estimated
  ivector ncolS(1,M);
  double log_alpha_lowerbound = nbinom1_flag==1 ? 0.001 : -5.0 ;
  ncolS = m;                       	// Uncorrelated random effects
  for (int i=1;i<=M;i++)                // Modifies the correlated ones
    if(cor_flag(i)>0)
      ncolS(i) = m(i)*(m(i)+1)/2;
  int nS = sum(ncolS);             	//  Total number
  pz.allocate(.000001,0.999,zi_phase,"pz");
  beta.allocate(1,p,1,"beta");
  real_beta.allocate(1,p,"real_beta");
  tmpL.allocate(1,ncolZ,-10,10.5,rand_phase,"tmpL");
  tmpL1.allocate(1,numb_cor_params,-10,10.5,cor_phase,"tmpL1");
  log_alpha.allocate(log_alpha_lowerbound,6.,alpha_phase,"log_alpha");
  alpha.allocate("alpha");
  S.allocate(1,nS,"S");
  u.allocate(1,sum_mq,rand_phase,"u");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  g.allocate("g");  /* ADOBJECTIVEFUNCTION */
  for (int i=I.indexmin() ; i<= I.indexmax(); i++)
  {
    if (min(I(i)) < u.indexmin() || max(I(i)) > u.indexmax())
    {
      cerr << "Bounds error in I for i = " << i << endl
         << " I(" << i << ") = " << I(i) << endl
         <<  "u.indexmin() = " << u.indexmin() << endl
          <<  "u.indexmax() = " << u.indexmax() << endl;
      ad_exit(1);
    }
  }
}
void model_parameters::userfunction(void)
{
  g =0.0;
  g=0.0;
  int i;
  if(!no_rand_flag)
    for (i=1;i<=sum_mq;i++)
      n01_prior(u(i));			// u's are N(0,1) distributed
  if (rlinkflag && !last_phase()) 
  {
    betapen(beta);
  }
  for(i=1;i<=n;i++)
    log_lik(i,tmpL,tmpL1,u(I(i)),beta,log_alpha,pz);
  if (sd_phase())
  {
    alpha = exp(log_alpha);
    real_beta=beta*phi;
    int i,j,i_m;
    int i1=1, i2=1, ii=1;
    for (i_m=1;i_m<=M;i_m++)
    {
      dvar_matrix L(1,m(i_m),1,m(i_m));
      L.initialize();
      dvar_matrix tmpS(1,m(i_m),1,m(i_m));
      tmpS.initialize();
      if(cor_flag(i_m)>0) // full cor structure
      {
        int ii=1;
        L(1,1)=1;
        for (i=1;i<=m(i_m);i++)
        {
          L(i,i)=1; // set diagonal
          for (int j=1;j<i;j++)
            L(i,j)=tmpL1(i2++); // fill in off-diag
          L(i)(1,i)/=norm(L(i)(1,i)); // ???
        }
        for (i=1;i<=m(i_m);i++)
          L(i)*=exp(tmpL(i1++)); //scale row
      }
      else  /// diagonal cor structure
      {
        for (i=1;i<=m(i_m);i++)
          L(i,i)=exp(tmpL(i1++));
      }
      tmpS=L*trans(L);  // square
      for (i=1;i<=m(i_m);i++) // fill in var-cov vector
      {
          if(cor_flag(i_m)>0)
            for(j=1;j<i;j++)
              S(ii++) = tmpS(i,j);
          S(ii++) = tmpS(i,i);
      }
    }
  }
			
}

void SEPFUN1  model_parameters::kludgepen(const prevariable&  v)
{
  begin_df1b2_funnel();
 g +=.5*square(v);
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::betapen(const dvar_vector&  v)
{
  begin_df1b2_funnel();
  g+=0.5*norm2(v);
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::n01_prior(const prevariable&  u)
{
  begin_df1b2_funnel();
 g -= -0.5*log(2.0*M_PI) - 0.5*square(u);
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::log_lik(int _i,const dvar_vector& tmpL,const dvar_vector& tmpL1,const dvar_vector& _ui, const dvar_vector& beta,const prevariable& log_alpha, const prevariable& pz)
{
  begin_df1b2_funnel();
  
  ADUNCONST(dvar_vector,ui)
  int i,j, i_m, Ni;
  double e1=1e-8; // formerly 1.e-20; current agrees with nbmm.tpl
  double e2=1e-8; // formerly 1.e-20; current agrees with nbmm.tpl
  double e3=1e-6; // Poisson prob=0 hack
  double e4=1e-10; // Poisson prob=0 hack for last phase
  dvariable alpha = e2+exp(log_alpha);
  // Construct random effects vector with proper var-covar structure from u
  dvar_vector b(1,ncolZ);
  int i1=1, i2=1;
  for (i_m=1;i_m<=M;i_m++)
  {
    dvar_matrix L(1,m(i_m),1,m(i_m));
    L.initialize();
    if(cor_flag(i_m)>0)
    {
      L(1,1)=1;
      for (i=1;i<=m(i_m);i++)
      {
        L(i,i)=1;
        for (int j=1;j<i;j++)
          L(i,j)=tmpL1(i2++);
        L(i)(1,i)/=norm(L(i)(1,i));
      }
      for (i=1;i<=m(i_m);i++)
        L(i)*=exp(tmpL(i1++));
    }
    else
    {
      for (i=1;i<=m(i_m);i++)
        L(i,i)=exp(tmpL(i1++));
    }
    int upper = sum(m(1,i_m));
    int lower = upper-m(i_m)+1;
    dvar_vector tmp1(1,m(i_m));
    //tmp1 = ui(lower,upper).shift(1);
 
    // FIXME: re-introduce rlinkflag here?
    if (initial_params::current_phase < initial_params::max_number_phases-1)
    {
      tmp1 = random_bound(ui(lower,upper).shift(1),5);
    }
    else if 
      (initial_params::current_phase == initial_params::max_number_phases-1)
    {
      tmp1 = random_bound(ui(lower,upper).shift(1),20);
    }
    else
    {
    
      tmp1 = ui(lower,upper).shift(1);
    }
    tmp1 = L*tmp1;
    b(lower,upper) = tmp1.shift(lower);
  }
  // fudge factors for inverse link
  double eps=1.e-2;
  double eps1=1.e-2;
  double eps2=1.e-4; // works on cloglog test in link.R; FAILS when eps2=1.e-6
  switch (current_phase())
  {
  case 1:
    eps=.01;
    eps1=.05;
    break;
  default:
    eps=.01;
    eps1=.01;
  }
  dvariable eta = X(_i)*beta + Z(_i)*b;
  if(has_offset)
    eta += offset(_i);
  dvariable lambda;
  dvariable fpen=0.0;
  switch(link_type_flag) 
  {
     case 0:    // [robust] log
       if (rlinkflag) {
          lambda = e1+mfexp(eta);
       } else {
          lambda = exp(eta);
       }
       break;
     case 1:    // [robust] logistic
       if (rlinkflag) {
        if (value(eta)<10)
          {
            dvariable eeta=mfexp(eta);
            lambda=eeta/(1.0+eeta);
          }
          else
         {
            dvariable eeta=mfexp(-eta);
            lambda=1.0/(1.0+eeta);
         }
       } else {
          lambda = 1.0/(1.0+exp(-eta));
       }
       if (initial_params::current_phase < initial_params::max_number_phases-1)
       {
         lambda=0.999*lambda+0.0005;
       }
       else if 
         (initial_params::current_phase == initial_params::max_number_phases-1)
       {
         lambda=0.999999*lambda+0.0000005;
       }
       break;
     case 2:   // probit (cum norm)
       lambda = cumd_norm(eta);
       break;
     case 3:   // inverse link
       if (rlinkflag) {
         lambda = 1.0/(eps+posfun(eta-eps,eps1,fpen));
          g+=fpen;
        } else {
         lambda = 1.0/eta;
        }
       break;
     case 4: // cloglog
	 {
       // FIXME: document/clarify epsilon values  Add rlinkflag?
       dvariable eeta;
       
       if (rlinkflag) {
          eeta = -mfexp(eta);
       } else {
          eeta = -exp(eta);
       }
       const double onesixth=1.0/6.0;
       const double one24=1.0/24.0;
       const double one120=1.0/120.0;
       // safe (1-exp()); tip from http://www.johndcook.com/cpp_expm1.html
       if (fabs(value(eeta))<1e-5) {
	   lambda = -eeta*
             (1.0+eeta*(0.5+eeta*(onesixth+eeta*(one24+eeta*one120))));
       } else {
	   lambda = 1-mfexp(eeta);
       }
       // restrict to (0,1)
      /*
       if (rlinkflag) {
          lambda = posfun(lambda,eps2,fpen);
          lambda = (1.0-posfun(1.0-lambda,eps2,fpen));
       }
      */
       if (rlinkflag && !last_phase()) 
       {
         lambda=0.999999*lambda+0.0000005;
       }
       }
       break;
     case 5: // identity
        lambda = eta;
	break;
     default:
       cerr << "Illegal value for link_type_flag" << endl;
       ad_exit(1);
  }
  dvariable  tau = nbinom1_flag ? alpha : 1.0 + e1 + lambda/alpha ;
  dvariable tmpl; 				// Log likelihood
  int cph=current_phase();
  // FIXME: does having lots of choices (for like_type_flag and link_flag) slow things down?
  // Is there any advantage to doing this stuff in a vectorized way?
  // Is there some other better approach to the per-point switch() ?
  switch(like_type_flag)
  {
    case 0:   // Poisson
	if (poiss_prob_bound==0) { 
	    tmpl =  log_density_poisson(y(_i,1),lambda);
	} else {
          if (cph<5)  
	    tmpl = log(e3+exp(log_density_poisson(y(_i,1),lambda)));  // DF hack May 2013
          else
	    tmpl = log(e4+exp(log_density_poisson(y(_i,1),lambda)));  // DF hack May 2013
	}
      break;
    case 1:   // Binomial: y(_i,1)=#successes, y(_i,2)=#failures, 
      if (p_y==1) {
         if (y(_i,1)==1) {
           tmpl = log(e1+lambda); //BMB: bvprobit.tpl uses e1=1e-10 here
         } else {
           tmpl = log(e1+(1.0-lambda));
         }
      }	 else {
         Ni = sum(y(_i));			// Number of trials
         tmpl = log_comb(Ni,y(_i,1)) + y(_i,1)*log(lambda) + (Ni-y(_i,1))*log(1.0-lambda);
      }
      break;
    case 2:   // neg binomial
      if (cph<2)  // would like to use alpha_phase rather than 2 but it's in local_calcs
      	   tmpl = -square(log(1.0+y(_i,1))-log(1.0+lambda));
      if (cph<4)  
	  tmpl = -(1.0+y(_i,1))*square(log(1.0+y(_i,1))-log(1.0+lambda));
      else
	  tmpl = log_negbinomial_density(y(_i,1),lambda,tau);
      break;
    case 3: // Gamma 
        tmpl = log_gamma_density(y(_i,1),alpha,alpha/lambda);
	break;
    case 4: // Beta 
      // FIXME: "log_beta_density" seems more consistent but changing name
      //       causes problems -- already exists somewhere?
      tmpl = ln_beta_density(y(_i,1),lambda,alpha);
      break;
    case 5: // Gaussian
      tmpl = -0.5*(log(2.0*M_PI))-log(alpha)-0.5*square((y(_i,1)-lambda)/alpha);	
      break; 
    case 6:   // truncated Poisson
      // FIXME: check somewhere (here, or preferably in R code) for trunc poisson + not ZI + 0 in response
      if (value(lambda) > 1.0e-10) {
          tmpl = log_density_poisson(y(_i,1),lambda)-log(1.0-exp(-lambda));
      } else {
          tmpl = log_density_poisson(y(_i,1),lambda)-log(lambda);
      }
      break;
    case 7:  // truncated NB
    // NB(0) = p^alpha = (alpha/(alpha+lambda))^alpha = (1+lambda/alpha)^(-alpha)
    //   -> exp(-lambda) as alpha -> infty
      if (cph<2)  // ignore zero-inflation for first phase
        tmpl = -square(log(1.0+y(_i,1))-log(1.0+lambda));
      else
        tmpl = log_negbinomial_density(y(_i,1),lambda,tau)-log(1.0-pow(1.0+lambda/alpha,-alpha));
      break;
    case 8: // logistic 
	tmpl = -log(alpha) + (y(_i,1)-lambda)/alpha - 2*log(1+exp((y(_i,1)-lambda)/alpha));
      break;
    case 9: // beta-binomial
	Ni = sum(y(_i));
	tmpl = log_comb(Ni,y(_i,1)) + // log(C(Ni,y(_i,1)))
	    gammln(y(_i,2)+alpha*(1-lambda))+gammln(y(_i,1)+alpha*lambda)-gammln(Ni+alpha) + // lbeta(...)
	    -(gammln(alpha*(1-lambda))+gammln(alpha*lambda)-gammln(alpha)); // lbeta(...)
	break;
  default:
      cerr << "Illegal value for like_type_flag" << endl;
      ad_exit(1);
  }
	
  // Zero inflation part 
  // zi_kluge: apply ZI whether or not zi_flag is true
  if(zi_flag || zi_kluge) {
    // if (zi_model_flag) {
    // 
    // }
    if(y(_i,1)==0)
      g -= log(e2+pz+(1.0-pz)*mfexp(tmpl));
    else
      g -= log(e2+(1.0-pz)*mfexp(tmpl));
   } else g -= tmpl;
  end_df1b2_funnel();
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << beta*phi << endl;
}

void model_parameters::preliminary_calculations(void) 
{

  admaster_slave_variable_interface(*this);
  cout << setprecision(4);
}
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
  ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize=40000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  gradient_structure::set_MAX_NVAR_OFFSET(2000000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    df1b2variable::noallocate=1;
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
    mp.preliminary_calculations();
    initial_df1b2params::separable_flag=1;
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

void df1b2_parameters::user_function(void)
{
  g =0.0;
  g=0.0;
  int i;
  if(!no_rand_flag)
    for (i=1;i<=sum_mq;i++)
      n01_prior(u(i));			// u's are N(0,1) distributed
  if (rlinkflag && !last_phase()) 
  {
    betapen(beta);
  }
  for(i=1;i<=n;i++)
    log_lik(i,tmpL,tmpL1,u(I(i)),beta,log_alpha,pz);
  if (sd_phase())
  {
    alpha = exp(log_alpha);
    real_beta=beta*phi;
    int i,j,i_m;
    int i1=1, i2=1, ii=1;
    for (i_m=1;i_m<=M;i_m++)
    {
      df1b2matrix L(1,m(i_m),1,m(i_m));
      L.initialize();
      df1b2matrix tmpS(1,m(i_m),1,m(i_m));
      tmpS.initialize();
      if(cor_flag(i_m)>0) // full cor structure
      {
        int ii=1;
        L(1,1)=1;
        for (i=1;i<=m(i_m);i++)
        {
          L(i,i)=1; // set diagonal
          for (int j=1;j<i;j++)
            L(i,j)=tmpL1(i2++); // fill in off-diag
          L(i)(1,i)/=norm(L(i)(1,i)); // ???
        }
        for (i=1;i<=m(i_m);i++)
          L(i)*=exp(tmpL(i1++)); //scale row
      }
      else  /// diagonal cor structure
      {
        for (i=1;i<=m(i_m);i++)
          L(i,i)=exp(tmpL(i1++));
      }
      tmpS=L*trans(L);  // square
      for (i=1;i<=m(i_m);i++) // fill in var-cov vector
      {
          if(cor_flag(i_m)>0)
            for(j=1;j<i;j++)
              S(ii++) = tmpS(i,j);
          S(ii++) = tmpS(i,i);
      }
    }
  }
			
}

void   df1b2_pre_parameters::kludgepen(const funnel_init_df1b2variable&  v)
{
  begin_df1b2_funnel();
 g +=.5*square(v);
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::betapen(const funnel_init_df1b2vector&  v)
{
  begin_df1b2_funnel();
  g+=0.5*norm2(v);
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::n01_prior(const funnel_init_df1b2variable&  u)
{
  begin_df1b2_funnel();
 g -= -0.5*log(2.0*M_PI) - 0.5*square(u);
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::log_lik(int _i,const funnel_init_df1b2vector& tmpL,const funnel_init_df1b2vector& tmpL1,const funnel_init_df1b2vector& _ui, const funnel_init_df1b2vector& beta,const funnel_init_df1b2variable& log_alpha, const funnel_init_df1b2variable& pz)
{
  begin_df1b2_funnel();
  
  ADUNCONST(df1b2vector,ui)
  int i,j, i_m, Ni;
  double e1=1e-8; // formerly 1.e-20; current agrees with nbmm.tpl
  double e2=1e-8; // formerly 1.e-20; current agrees with nbmm.tpl
  double e3=1e-6; // Poisson prob=0 hack
  double e4=1e-10; // Poisson prob=0 hack for last phase
  df1b2variable alpha = e2+exp(log_alpha);
  // Construct random effects vector with proper var-covar structure from u
  df1b2vector b(1,ncolZ);
  int i1=1, i2=1;
  for (i_m=1;i_m<=M;i_m++)
  {
    df1b2matrix L(1,m(i_m),1,m(i_m));
    L.initialize();
    if(cor_flag(i_m)>0)
    {
      L(1,1)=1;
      for (i=1;i<=m(i_m);i++)
      {
        L(i,i)=1;
        for (int j=1;j<i;j++)
          L(i,j)=tmpL1(i2++);
        L(i)(1,i)/=norm(L(i)(1,i));
      }
      for (i=1;i<=m(i_m);i++)
        L(i)*=exp(tmpL(i1++));
    }
    else
    {
      for (i=1;i<=m(i_m);i++)
        L(i,i)=exp(tmpL(i1++));
    }
    int upper = sum(m(1,i_m));
    int lower = upper-m(i_m)+1;
    df1b2vector tmp1(1,m(i_m));
    //tmp1 = ui(lower,upper).shift(1);
 
    // FIXME: re-introduce rlinkflag here?
    if (initial_params::current_phase < initial_params::max_number_phases-1)
    {
      tmp1 = random_bound(ui(lower,upper).shift(1),5);
    }
    else if 
      (initial_params::current_phase == initial_params::max_number_phases-1)
    {
      tmp1 = random_bound(ui(lower,upper).shift(1),20);
    }
    else
    {
    
      tmp1 = ui(lower,upper).shift(1);
    }
    tmp1 = L*tmp1;
    b(lower,upper) = tmp1.shift(lower);
  }
  // fudge factors for inverse link
  double eps=1.e-2;
  double eps1=1.e-2;
  double eps2=1.e-4; // works on cloglog test in link.R; FAILS when eps2=1.e-6
  switch (current_phase())
  {
  case 1:
    eps=.01;
    eps1=.05;
    break;
  default:
    eps=.01;
    eps1=.01;
  }
  df1b2variable eta = X(_i)*beta + Z(_i)*b;
  if(has_offset)
    eta += offset(_i);
  df1b2variable lambda;
  df1b2variable fpen=0.0;
  switch(link_type_flag) 
  {
     case 0:    // [robust] log
       if (rlinkflag) {
          lambda = e1+mfexp(eta);
       } else {
          lambda = exp(eta);
       }
       break;
     case 1:    // [robust] logistic
       if (rlinkflag) {
        if (value(eta)<10)
          {
            df1b2variable eeta=mfexp(eta);
            lambda=eeta/(1.0+eeta);
          }
          else
         {
            df1b2variable eeta=mfexp(-eta);
            lambda=1.0/(1.0+eeta);
         }
       } else {
          lambda = 1.0/(1.0+exp(-eta));
       }
       if (initial_params::current_phase < initial_params::max_number_phases-1)
       {
         lambda=0.999*lambda+0.0005;
       }
       else if 
         (initial_params::current_phase == initial_params::max_number_phases-1)
       {
         lambda=0.999999*lambda+0.0000005;
       }
       break;
     case 2:   // probit (cum norm)
       lambda = cumd_norm(eta);
       break;
     case 3:   // inverse link
       if (rlinkflag) {
         lambda = 1.0/(eps+posfun(eta-eps,eps1,fpen));
          g+=fpen;
        } else {
         lambda = 1.0/eta;
        }
       break;
     case 4: // cloglog
	 {
       // FIXME: document/clarify epsilon values  Add rlinkflag?
       df1b2variable eeta;
       
       if (rlinkflag) {
          eeta = -mfexp(eta);
       } else {
          eeta = -exp(eta);
       }
       const double onesixth=1.0/6.0;
       const double one24=1.0/24.0;
       const double one120=1.0/120.0;
       // safe (1-exp()); tip from http://www.johndcook.com/cpp_expm1.html
       if (fabs(value(eeta))<1e-5) {
	   lambda = -eeta*
             (1.0+eeta*(0.5+eeta*(onesixth+eeta*(one24+eeta*one120))));
       } else {
	   lambda = 1-mfexp(eeta);
       }
       // restrict to (0,1)
      /*
       if (rlinkflag) {
          lambda = posfun(lambda,eps2,fpen);
          lambda = (1.0-posfun(1.0-lambda,eps2,fpen));
       }
      */
       if (rlinkflag && !last_phase()) 
       {
         lambda=0.999999*lambda+0.0000005;
       }
       }
       break;
     case 5: // identity
        lambda = eta;
	break;
     default:
       cerr << "Illegal value for link_type_flag" << endl;
       ad_exit(1);
  }
  df1b2variable  tau = nbinom1_flag ? alpha : 1.0 + e1 + lambda/alpha ;
  df1b2variable tmpl; 				// Log likelihood
  int cph=current_phase();
  // FIXME: does having lots of choices (for like_type_flag and link_flag) slow things down?
  // Is there any advantage to doing this stuff in a vectorized way?
  // Is there some other better approach to the per-point switch() ?
  switch(like_type_flag)
  {
    case 0:   // Poisson
	if (poiss_prob_bound==0) { 
	    tmpl =  log_density_poisson(y(_i,1),lambda);
	} else {
          if (cph<5)  
	    tmpl = log(e3+exp(log_density_poisson(y(_i,1),lambda)));  // DF hack May 2013
          else
	    tmpl = log(e4+exp(log_density_poisson(y(_i,1),lambda)));  // DF hack May 2013
	}
      break;
    case 1:   // Binomial: y(_i,1)=#successes, y(_i,2)=#failures, 
      if (p_y==1) {
         if (y(_i,1)==1) {
           tmpl = log(e1+lambda); //BMB: bvprobit.tpl uses e1=1e-10 here
         } else {
           tmpl = log(e1+(1.0-lambda));
         }
      }	 else {
         Ni = sum(y(_i));			// Number of trials
         tmpl = log_comb(Ni,y(_i,1)) + y(_i,1)*log(lambda) + (Ni-y(_i,1))*log(1.0-lambda);
      }
      break;
    case 2:   // neg binomial
      if (cph<2)  // would like to use alpha_phase rather than 2 but it's in local_calcs
      	   tmpl = -square(log(1.0+y(_i,1))-log(1.0+lambda));
      if (cph<4)  
	  tmpl = -(1.0+y(_i,1))*square(log(1.0+y(_i,1))-log(1.0+lambda));
      else
	  tmpl = log_negbinomial_density(y(_i,1),lambda,tau);
      break;
    case 3: // Gamma 
        tmpl = log_gamma_density(y(_i,1),alpha,alpha/lambda);
	break;
    case 4: // Beta 
      // FIXME: "log_beta_density" seems more consistent but changing name
      //       causes problems -- already exists somewhere?
      tmpl = ln_beta_density(y(_i,1),lambda,alpha);
      break;
    case 5: // Gaussian
      tmpl = -0.5*(log(2.0*M_PI))-log(alpha)-0.5*square((y(_i,1)-lambda)/alpha);	
      break; 
    case 6:   // truncated Poisson
      // FIXME: check somewhere (here, or preferably in R code) for trunc poisson + not ZI + 0 in response
      if (value(lambda) > 1.0e-10) {
          tmpl = log_density_poisson(y(_i,1),lambda)-log(1.0-exp(-lambda));
      } else {
          tmpl = log_density_poisson(y(_i,1),lambda)-log(lambda);
      }
      break;
    case 7:  // truncated NB
    // NB(0) = p^alpha = (alpha/(alpha+lambda))^alpha = (1+lambda/alpha)^(-alpha)
    //   -> exp(-lambda) as alpha -> infty
      if (cph<2)  // ignore zero-inflation for first phase
        tmpl = -square(log(1.0+y(_i,1))-log(1.0+lambda));
      else
        tmpl = log_negbinomial_density(y(_i,1),lambda,tau)-log(1.0-pow(1.0+lambda/alpha,-alpha));
      break;
    case 8: // logistic 
	tmpl = -log(alpha) + (y(_i,1)-lambda)/alpha - 2*log(1+exp((y(_i,1)-lambda)/alpha));
      break;
    case 9: // beta-binomial
	Ni = sum(y(_i));
	tmpl = log_comb(Ni,y(_i,1)) + // log(C(Ni,y(_i,1)))
	    gammln(y(_i,2)+alpha*(1-lambda))+gammln(y(_i,1)+alpha*lambda)-gammln(Ni+alpha) + // lbeta(...)
	    -(gammln(alpha*(1-lambda))+gammln(alpha*lambda)-gammln(alpha)); // lbeta(...)
	break;
  default:
      cerr << "Illegal value for like_type_flag" << endl;
      ad_exit(1);
  }
	
  // Zero inflation part 
  // zi_kluge: apply ZI whether or not zi_flag is true
  if(zi_flag || zi_kluge) {
    // if (zi_model_flag) {
    // 
    // }
    if(y(_i,1)==0)
      g -= log(e2+pz+(1.0-pz)*mfexp(tmpl));
    else
      g -= log(e2+(1.0-pz)*mfexp(tmpl));
   } else g -= tmpl;
  end_df1b2_funnel();
}
   
void df1b2_pre_parameters::setup_quadprior_calcs(void) 
{ 
  df1b2_gradlist::set_no_derivatives(); 
  quadratic_prior::in_qp_calculations=1; 
}  
  
void df1b2_pre_parameters::begin_df1b2_funnel(void) 
{ 
  (*re_objective_function_value::pobjfun)=0; 
  other_separable_stuff_begin(); 
  f1b2gradlist->reset();  
  if (!quadratic_prior::in_qp_calculations) 
  { 
    df1b2_gradlist::set_yes_derivatives();  
  } 
  funnel_init_var::allocate_all();  
}  
 
void df1b2_pre_parameters::end_df1b2_funnel(void) 
{  
  lapprox->do_separable_stuff(); 
  other_separable_stuff_end(); 
} 
  
void model_parameters::begin_df1b2_funnel(void) 
{ 
  if (lapprox)  
  {  
    {  
      begin_funnel_stuff();  
    }  
  }  
}  
 
void model_parameters::end_df1b2_funnel(void) 
{  
  if (lapprox)  
  {  
    end_df1b2_funnel_stuff();  
  }  
} 

void df1b2_parameters::allocate(void) 
{
  // BMB: FIXME: do we need this?  formerly disallowed for binomial (was like_type_flag 2);
  //     would be problematic for binary data but otherwise OK.  Should test in R code
  //  if(zi_flag && (like_type_flag>=2))
  // {
  //  cerr << "Zero inflation not allowed for this response type" << endl;
  //  ad_exit(1);
  // }
  // Determines "phases", i.e. when the various parameters  becomes active in the optimization process
  int pctr = 2;		// "Current phase" in the stagewise procedure
  // FIXME: move trunc_poisson, logistic earlier in like_type hierarchy (or add a scale parameter flag vector)
  int alpha_phase = like_type_flag>1 && like_type_flag!=6 ? pctr++ : -1;        // Phase 2 if active
  int zi_phase = zi_flag ? pctr++ : -1;                      			// After alpha
  int rand_phase = no_rand_flag==0 ? pctr++ : -1;    				// SD of RE's
  int cor_phase = (rand_phase>0) && (sum(cor_flag)>0) ? pctr++ : -1 ; 		// Correlations of RE's
  // Count the number of variance/correlation parameters to be estimated
  ivector ncolS(1,M);
  double log_alpha_lowerbound = nbinom1_flag==1 ? 0.001 : -5.0 ;
  ncolS = m;                       	// Uncorrelated random effects
  for (int i=1;i<=M;i++)                // Modifies the correlated ones
    if(cor_flag(i)>0)
      ncolS(i) = m(i)*(m(i)+1)/2;
  int nS = sum(ncolS);             	//  Total number
  pz.allocate(.000001,0.999,zi_phase,"pz");
  beta.allocate(1,p,1,"beta");
  real_beta.allocate(1,p,"real_beta");
  tmpL.allocate(1,ncolZ,-10,10.5,rand_phase,"tmpL");
  tmpL1.allocate(1,numb_cor_params,-10,10.5,cor_phase,"tmpL1");
  log_alpha.allocate(log_alpha_lowerbound,6.,alpha_phase,"log_alpha");
  alpha.allocate("alpha");
  S.allocate(1,nS,"S");
  u.allocate(1,sum_mq,rand_phase,"u");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  g.allocate("g");  /* ADOBJECTIVEFUNCTION */
  for (int i=I.indexmin() ; i<= I.indexmax(); i++)
  {
    if (min(I(i)) < u.indexmin() || max(I(i)) > u.indexmax())
    {
      cerr << "Bounds error in I for i = " << i << endl
         << " I(" << i << ") = " << I(i) << endl
         <<  "u.indexmin() = " << u.indexmin() << endl
          <<  "u.indexmax() = " << u.indexmax() << endl;
      ad_exit(1);
    }
  }
}
