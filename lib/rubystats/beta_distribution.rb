require 'rubystats/probability_distribution'

module Rubystats
  class BetaDistribution < Rubystats::ProbabilityDistribution 
    include Rubystats::SpecialMath

    attr_reader :p, :q

    #dgr_p = degrees of freedom p
    #dgr_q = degrees of freedom q
    def initialize(dgr_p, dgr_q)
      if dgr_p <= 0 || dgr_q <= 0
        raise ArgumentError.new("Paramters must be greater than zero.")
      end
      @p = dgr_p.to_f
      @q = dgr_q.to_f
    end

    def mean
      @p.to_f / (@p.to_f + @q.to_f)
    end

    def standard_deviation
      Math.sqrt(@p * @q / ((@p + @q)**2 * (@p + @q + 1)))
    end

    def pdf(x) 
      if x.class == Array 
        pdf_vals = []
        for i in (0 ... x.size) 
          check_range(x[i])
          if x[i] == 0.0 || x[i] == 1.0 
            pdf_vals[i] = 0.0
          else 
            pdf_vals[i] = Math.exp( - log_beta(@p, @q) + (@p - 1.0) * Math.log(x[i]) + (@q - 1.0) * Math.log(1.0 - x[i]))
          end
        end
        return pdf_vals
      else 
        check_range(x)
        if  (x == 0.0) || (x == 1.0)  
          return 0.0
        else 
          return Math.exp( -1.0 * log_beta(@p, @q) + 
          (@p - 1.0) * Math.log(x) + 
          (@q - 1.0) * Math.log(1.0 - x))
        end
      end
    end

    def cdf(x) 
      if x.class == Array 
        cdf_vals = Array.new
        for i in 0 ... x.size 
          check_range(x[i])
          cdf_vals[i] = incomplete_beta(x[i], @p, @q)
        end
        return cdf_vals
      else 
        check_range(x)
        cdf_val = incomplete_beta(x, @p, @q)
        return cdf_val
      end
    end

    def icdf(prob) 
      if prob.class == Array 
        inv_vals = Array.new
        for i in 0 ... prob.size
          check_range(prob[i])
          if prob[i] == 0.0 
            inv_vals[i] = 0.0
          end
          if prob[i] == 1.0
            inv_vals[i] = 1.0
          end
          inv_vals[i] = find_root(prob[i], 0.5, 0.0, 1.0)
        end
        return inv_vals
      else 
        check_range(prob)
        return 0.0 if prob == 0.0
        return 1.0 if prob == 1.0
        return find_root(prob, 0.5, 0.0, 1.0)
      end
    end

    # From http://www.xycoon.com/beta_randomnumbers.htm
    # Let Gam(1,a) ~ ln CartProduct(1,a) U(0,1)
    # Let Gam(1,b) ~ ln CartProduct(1,b) U(0,1)
    # Then Beta(a,b) ~ Gam(1,a)/(Gam(1,a)+Gam(1,b))

    def get_rng # p and q must be integers
      p_int, q_int = p.to_i, q.to_i
      raise ArgumentError, "P's and Q's must be integers for this random number generator." if p_int!=p||q_int!=q
      raise ArgumentError, "P's and Q's must be greater than 1 for this random number generator." if p_int < 1 || q_int < 1
      gamma_p = gamma_approx(p_int)
      gamma_q = gamma_approx(q_int)
      return gamma_p/(gamma_p+gamma_q)
    end

    private

    def gamma_approx(x)
      Math.log(Array(1..x).inject(1){|cart_prod, v| cart_prod*rand})
    end

  end
end
