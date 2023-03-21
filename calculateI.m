function I = calculateI(tf, hf, bf)

    I = (hf^3 * tf)/12 + 2*((tf^3 * bf)/12) + 2*(bf*tf*(hf/2 -tf/2)^2);

end