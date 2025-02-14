OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4053722) q[0];
sx q[0];
rz(-0.020981941) q[0];
sx q[0];
rz(-1.1809281) q[0];
rz(-0.48038545) q[1];
sx q[1];
rz(-1.9863702) q[1];
sx q[1];
rz(2.6748599) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48343998) q[0];
sx q[0];
rz(-2.3077093) q[0];
sx q[0];
rz(-1.9574653) q[0];
x q[1];
rz(-2.4619815) q[2];
sx q[2];
rz(-0.70290297) q[2];
sx q[2];
rz(1.694569) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3148823) q[1];
sx q[1];
rz(-0.2383543) q[1];
sx q[1];
rz(1.5520067) q[1];
rz(-pi) q[2];
rz(-2.1331514) q[3];
sx q[3];
rz(-0.97259023) q[3];
sx q[3];
rz(1.4149208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0878318) q[2];
sx q[2];
rz(-1.6873282) q[2];
sx q[2];
rz(1.9161179) q[2];
rz(2.7671704) q[3];
sx q[3];
rz(-1.1332847) q[3];
sx q[3];
rz(-2.5882914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9541009) q[0];
sx q[0];
rz(-2.2611389) q[0];
sx q[0];
rz(0.53400293) q[0];
rz(-2.7502637) q[1];
sx q[1];
rz(-2.6983039) q[1];
sx q[1];
rz(2.2543529) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0011883) q[0];
sx q[0];
rz(-2.6713704) q[0];
sx q[0];
rz(-2.4109439) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0897691) q[2];
sx q[2];
rz(-0.88406056) q[2];
sx q[2];
rz(-2.9406282) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56480592) q[1];
sx q[1];
rz(-0.24165711) q[1];
sx q[1];
rz(0.56743001) q[1];
rz(-pi) q[2];
rz(1.0530472) q[3];
sx q[3];
rz(-1.8705006) q[3];
sx q[3];
rz(-0.24417711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.54681626) q[2];
sx q[2];
rz(-1.8211326) q[2];
sx q[2];
rz(2.7776264) q[2];
rz(-1.4173896) q[3];
sx q[3];
rz(-0.28984362) q[3];
sx q[3];
rz(-2.3072306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.904838) q[0];
sx q[0];
rz(-0.098243864) q[0];
sx q[0];
rz(2.8114317) q[0];
rz(2.5579021) q[1];
sx q[1];
rz(-0.81283641) q[1];
sx q[1];
rz(-2.6469753) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0000293) q[0];
sx q[0];
rz(-2.7589679) q[0];
sx q[0];
rz(-2.3429246) q[0];
x q[1];
rz(1.0801267) q[2];
sx q[2];
rz(-1.3441836) q[2];
sx q[2];
rz(-1.8267711) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51926326) q[1];
sx q[1];
rz(-2.0262782) q[1];
sx q[1];
rz(-1.0002656) q[1];
x q[2];
rz(1.4026138) q[3];
sx q[3];
rz(-0.59051149) q[3];
sx q[3];
rz(-1.8174949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3459449) q[2];
sx q[2];
rz(-2.4730885) q[2];
sx q[2];
rz(3.0188959) q[2];
rz(-3.0986541) q[3];
sx q[3];
rz(-2.0624845) q[3];
sx q[3];
rz(0.19449657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73404679) q[0];
sx q[0];
rz(-2.6341944) q[0];
sx q[0];
rz(0.88974446) q[0];
rz(0.45097688) q[1];
sx q[1];
rz(-0.86800066) q[1];
sx q[1];
rz(0.45423347) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12338598) q[0];
sx q[0];
rz(-2.7834547) q[0];
sx q[0];
rz(-0.30787719) q[0];
rz(1.5726015) q[2];
sx q[2];
rz(-1.1434525) q[2];
sx q[2];
rz(-1.3722562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.52773038) q[1];
sx q[1];
rz(-2.3017028) q[1];
sx q[1];
rz(-1.1682603) q[1];
rz(-2.9353981) q[3];
sx q[3];
rz(-2.7795994) q[3];
sx q[3];
rz(-0.21807204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1526327) q[2];
sx q[2];
rz(-2.5071414) q[2];
sx q[2];
rz(-0.75759849) q[2];
rz(2.9518413) q[3];
sx q[3];
rz(-1.8228143) q[3];
sx q[3];
rz(0.54722133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309288) q[0];
sx q[0];
rz(-2.3543816) q[0];
sx q[0];
rz(-1.2934562) q[0];
rz(-2.5594607) q[1];
sx q[1];
rz(-1.4385834) q[1];
sx q[1];
rz(-0.17939803) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6027733) q[0];
sx q[0];
rz(-1.7742549) q[0];
sx q[0];
rz(-0.18360965) q[0];
x q[1];
rz(-2.8371528) q[2];
sx q[2];
rz(-2.4477262) q[2];
sx q[2];
rz(2.2562698) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.071956858) q[1];
sx q[1];
rz(-1.3838125) q[1];
sx q[1];
rz(0.93902875) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87406335) q[3];
sx q[3];
rz(-1.1723601) q[3];
sx q[3];
rz(1.5521727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3680215) q[2];
sx q[2];
rz(-2.7715235) q[2];
sx q[2];
rz(0.5698815) q[2];
rz(2.701345) q[3];
sx q[3];
rz(-1.8972242) q[3];
sx q[3];
rz(-0.35278916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2920947) q[0];
sx q[0];
rz(-0.26420132) q[0];
sx q[0];
rz(-2.0000892) q[0];
rz(-1.796272) q[1];
sx q[1];
rz(-2.1514386) q[1];
sx q[1];
rz(0.17123953) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6692176) q[0];
sx q[0];
rz(-2.6504189) q[0];
sx q[0];
rz(-0.21531658) q[0];
rz(-1.0467347) q[2];
sx q[2];
rz(-1.9643133) q[2];
sx q[2];
rz(-1.524802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1223733) q[1];
sx q[1];
rz(-1.7369441) q[1];
sx q[1];
rz(-2.5801587) q[1];
rz(-1.5909252) q[3];
sx q[3];
rz(-0.55680767) q[3];
sx q[3];
rz(2.0272154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.77219999) q[2];
sx q[2];
rz(-0.93595305) q[2];
sx q[2];
rz(-0.13747036) q[2];
rz(-1.8933659) q[3];
sx q[3];
rz(-2.5745001) q[3];
sx q[3];
rz(1.9198157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0565979) q[0];
sx q[0];
rz(-0.63768142) q[0];
sx q[0];
rz(1.6531264) q[0];
rz(0.78530637) q[1];
sx q[1];
rz(-1.4005125) q[1];
sx q[1];
rz(-1.3135501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5851597) q[0];
sx q[0];
rz(-1.6939511) q[0];
sx q[0];
rz(0.086104546) q[0];
x q[1];
rz(-1.4573757) q[2];
sx q[2];
rz(-1.7608425) q[2];
sx q[2];
rz(2.6148877) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.056712778) q[1];
sx q[1];
rz(-3.0482349) q[1];
sx q[1];
rz(-2.8356524) q[1];
rz(0.87487674) q[3];
sx q[3];
rz(-1.531344) q[3];
sx q[3];
rz(0.46515282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58275783) q[2];
sx q[2];
rz(-1.02966) q[2];
sx q[2];
rz(-2.6666857) q[2];
rz(-0.20259914) q[3];
sx q[3];
rz(-0.83297268) q[3];
sx q[3];
rz(-1.9995662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1982034) q[0];
sx q[0];
rz(-3.0084963) q[0];
sx q[0];
rz(2.8357847) q[0];
rz(-2.7022779) q[1];
sx q[1];
rz(-0.82249928) q[1];
sx q[1];
rz(1.8261212) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02365774) q[0];
sx q[0];
rz(-2.2542037) q[0];
sx q[0];
rz(1.5562431) q[0];
rz(1.1313296) q[2];
sx q[2];
rz(-0.60053289) q[2];
sx q[2];
rz(0.10490049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2478359) q[1];
sx q[1];
rz(-0.79036056) q[1];
sx q[1];
rz(-1.7414581) q[1];
rz(-pi) q[2];
rz(-2.148904) q[3];
sx q[3];
rz(-0.6738014) q[3];
sx q[3];
rz(-1.0579619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7266015) q[2];
sx q[2];
rz(-1.5433658) q[2];
sx q[2];
rz(-2.7169363) q[2];
rz(-0.031938227) q[3];
sx q[3];
rz(-1.5141124) q[3];
sx q[3];
rz(-2.4922075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53751078) q[0];
sx q[0];
rz(-0.6518971) q[0];
sx q[0];
rz(-2.5867468) q[0];
rz(-2.8765053) q[1];
sx q[1];
rz(-1.4684497) q[1];
sx q[1];
rz(-1.2010942) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5145668) q[0];
sx q[0];
rz(-1.1794588) q[0];
sx q[0];
rz(-2.5364801) q[0];
rz(-pi) q[1];
rz(2.8657929) q[2];
sx q[2];
rz(-1.7364304) q[2];
sx q[2];
rz(1.5716108) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7857901) q[1];
sx q[1];
rz(-1.5340163) q[1];
sx q[1];
rz(2.793502) q[1];
rz(-0.026215629) q[3];
sx q[3];
rz(-2.391572) q[3];
sx q[3];
rz(1.9103844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25192866) q[2];
sx q[2];
rz(-0.98968518) q[2];
sx q[2];
rz(0.84452334) q[2];
rz(2.0097513) q[3];
sx q[3];
rz(-2.2364538) q[3];
sx q[3];
rz(1.7822781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97245401) q[0];
sx q[0];
rz(-0.3903946) q[0];
sx q[0];
rz(-1.1727232) q[0];
rz(-0.68924618) q[1];
sx q[1];
rz(-2.0447562) q[1];
sx q[1];
rz(-0.26395878) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.496752) q[0];
sx q[0];
rz(-0.94366108) q[0];
sx q[0];
rz(1.0154614) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9969292) q[2];
sx q[2];
rz(-0.80748122) q[2];
sx q[2];
rz(-1.5308183) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86706272) q[1];
sx q[1];
rz(-1.8856876) q[1];
sx q[1];
rz(2.9462183) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12216025) q[3];
sx q[3];
rz(-2.2111243) q[3];
sx q[3];
rz(-1.8807008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3802203) q[2];
sx q[2];
rz(-1.994588) q[2];
sx q[2];
rz(2.0210361) q[2];
rz(-1.5348966) q[3];
sx q[3];
rz(-0.94529072) q[3];
sx q[3];
rz(1.631261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0705538) q[0];
sx q[0];
rz(-0.65503913) q[0];
sx q[0];
rz(-0.10792637) q[0];
rz(-2.4212266) q[1];
sx q[1];
rz(-1.9279059) q[1];
sx q[1];
rz(2.7005213) q[1];
rz(0.60463641) q[2];
sx q[2];
rz(-2.2199134) q[2];
sx q[2];
rz(1.8668957) q[2];
rz(1.4955487) q[3];
sx q[3];
rz(-0.9552707) q[3];
sx q[3];
rz(0.09494119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
