OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95028967) q[0];
sx q[0];
rz(6.0130881) q[0];
sx q[0];
rz(10.313378) q[0];
rz(-1.3130045) q[1];
sx q[1];
rz(-1.5994025) q[1];
sx q[1];
rz(1.7607652) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61533538) q[0];
sx q[0];
rz(-1.7563631) q[0];
sx q[0];
rz(0.19677563) q[0];
x q[1];
rz(-2.1031453) q[2];
sx q[2];
rz(-0.84538904) q[2];
sx q[2];
rz(0.40276819) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9819698) q[1];
sx q[1];
rz(-1.9062942) q[1];
sx q[1];
rz(1.5606176) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9350151) q[3];
sx q[3];
rz(-1.7813588) q[3];
sx q[3];
rz(-0.59277804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8218653) q[2];
sx q[2];
rz(-1.3250019) q[2];
sx q[2];
rz(2.3925609) q[2];
rz(0.15394112) q[3];
sx q[3];
rz(-1.0356244) q[3];
sx q[3];
rz(3.0416987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6946436) q[0];
sx q[0];
rz(-1.4780937) q[0];
sx q[0];
rz(-1.9248167) q[0];
rz(1.0617537) q[1];
sx q[1];
rz(-0.80445015) q[1];
sx q[1];
rz(2.7144576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9574643) q[0];
sx q[0];
rz(-1.2250568) q[0];
sx q[0];
rz(-2.5347802) q[0];
rz(-pi) q[1];
rz(-2.8982694) q[2];
sx q[2];
rz(-0.93987018) q[2];
sx q[2];
rz(2.6532946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4915575) q[1];
sx q[1];
rz(-1.8582398) q[1];
sx q[1];
rz(0.14658714) q[1];
x q[2];
rz(1.8187567) q[3];
sx q[3];
rz(-0.21020791) q[3];
sx q[3];
rz(1.5402286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0624258) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(1.742935) q[2];
rz(0.22377293) q[3];
sx q[3];
rz(-2.2327435) q[3];
sx q[3];
rz(1.5980501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021521213) q[0];
sx q[0];
rz(-1.3503617) q[0];
sx q[0];
rz(-0.045510005) q[0];
rz(-1.1687763) q[1];
sx q[1];
rz(-1.2380506) q[1];
sx q[1];
rz(2.5659836) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92514491) q[0];
sx q[0];
rz(-0.63918221) q[0];
sx q[0];
rz(-3.1228288) q[0];
rz(-pi) q[1];
rz(-3.124442) q[2];
sx q[2];
rz(-0.91020012) q[2];
sx q[2];
rz(-1.9018381) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82390416) q[1];
sx q[1];
rz(-1.0270734) q[1];
sx q[1];
rz(-0.25441092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91797249) q[3];
sx q[3];
rz(-1.550972) q[3];
sx q[3];
rz(-1.4696079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0098972926) q[2];
sx q[2];
rz(-2.3991149) q[2];
sx q[2];
rz(-0.95727813) q[2];
rz(2.5668528) q[3];
sx q[3];
rz(-1.5301907) q[3];
sx q[3];
rz(-1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15605536) q[0];
sx q[0];
rz(-2.1888581) q[0];
sx q[0];
rz(-1.8804469) q[0];
rz(3.0592697) q[1];
sx q[1];
rz(-1.0876834) q[1];
sx q[1];
rz(-1.3538768) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60786965) q[0];
sx q[0];
rz(-1.3114616) q[0];
sx q[0];
rz(1.813867) q[0];
x q[1];
rz(2.8026695) q[2];
sx q[2];
rz(-0.51033516) q[2];
sx q[2];
rz(-1.5477382) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8396047) q[1];
sx q[1];
rz(-1.0763554) q[1];
sx q[1];
rz(1.5093263) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1534377) q[3];
sx q[3];
rz(-2.372962) q[3];
sx q[3];
rz(0.60291327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7202619) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(1.2565695) q[2];
rz(-0.71150696) q[3];
sx q[3];
rz(-1.810775) q[3];
sx q[3];
rz(0.28212696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8307777) q[0];
sx q[0];
rz(-2.2956235) q[0];
sx q[0];
rz(2.7984483) q[0];
rz(0.061773069) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(1.4168581) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8648711) q[0];
sx q[0];
rz(-1.7407932) q[0];
sx q[0];
rz(1.2890588) q[0];
rz(-pi) q[1];
rz(-0.70059641) q[2];
sx q[2];
rz(-1.3178409) q[2];
sx q[2];
rz(2.7299936) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2572391) q[1];
sx q[1];
rz(-1.4352918) q[1];
sx q[1];
rz(-2.0550904) q[1];
rz(-2.9067578) q[3];
sx q[3];
rz(-1.9855301) q[3];
sx q[3];
rz(1.0339586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5477649) q[2];
sx q[2];
rz(-1.6739028) q[2];
sx q[2];
rz(-1.0859547) q[2];
rz(-0.91935277) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(2.5277188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74349657) q[0];
sx q[0];
rz(-1.9323823) q[0];
sx q[0];
rz(-0.79291517) q[0];
rz(2.4010557) q[1];
sx q[1];
rz(-0.99895993) q[1];
sx q[1];
rz(-0.83121306) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7949267) q[0];
sx q[0];
rz(-1.3627317) q[0];
sx q[0];
rz(1.1870334) q[0];
rz(-pi) q[1];
rz(0.0074499091) q[2];
sx q[2];
rz(-1.2937355) q[2];
sx q[2];
rz(1.3598833) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38922849) q[1];
sx q[1];
rz(-2.3550451) q[1];
sx q[1];
rz(-1.1334555) q[1];
x q[2];
rz(1.6115723) q[3];
sx q[3];
rz(-2.6705461) q[3];
sx q[3];
rz(2.149793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0703766) q[2];
sx q[2];
rz(-1.2717609) q[2];
sx q[2];
rz(2.0224723) q[2];
rz(1.8708771) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(1.3814111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0041466) q[0];
sx q[0];
rz(-0.16682145) q[0];
sx q[0];
rz(-1.5995837) q[0];
rz(2.1265325) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(-2.9687845) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61758608) q[0];
sx q[0];
rz(-0.79389555) q[0];
sx q[0];
rz(2.2909597) q[0];
rz(-pi) q[1];
rz(2.4577499) q[2];
sx q[2];
rz(-1.2798066) q[2];
sx q[2];
rz(3.0802205) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0513251) q[1];
sx q[1];
rz(-2.6410612) q[1];
sx q[1];
rz(-1.0330908) q[1];
x q[2];
rz(2.1234346) q[3];
sx q[3];
rz(-2.2595836) q[3];
sx q[3];
rz(-1.3017201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.479594) q[2];
sx q[2];
rz(-0.84344232) q[2];
sx q[2];
rz(2.1638828) q[2];
rz(-0.57502037) q[3];
sx q[3];
rz(-1.9366555) q[3];
sx q[3];
rz(-1.2303111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32790312) q[0];
sx q[0];
rz(-0.29569018) q[0];
sx q[0];
rz(3.1258702) q[0];
rz(1.188259) q[1];
sx q[1];
rz(-2.5091722) q[1];
sx q[1];
rz(1.1994919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2000113) q[0];
sx q[0];
rz(-2.0111548) q[0];
sx q[0];
rz(2.8811127) q[0];
rz(-pi) q[1];
rz(3.0891782) q[2];
sx q[2];
rz(-2.0611066) q[2];
sx q[2];
rz(0.744396) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6576865) q[1];
sx q[1];
rz(-0.71345854) q[1];
sx q[1];
rz(-0.76303996) q[1];
x q[2];
rz(-2.8335081) q[3];
sx q[3];
rz(-1.020806) q[3];
sx q[3];
rz(-2.1717482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99379313) q[2];
sx q[2];
rz(-1.9027998) q[2];
sx q[2];
rz(1.7564868) q[2];
rz(-2.5177453) q[3];
sx q[3];
rz(-1.8892989) q[3];
sx q[3];
rz(-1.0998211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95595908) q[0];
sx q[0];
rz(-0.82056844) q[0];
sx q[0];
rz(-0.69325915) q[0];
rz(2.8765826) q[1];
sx q[1];
rz(-0.82844228) q[1];
sx q[1];
rz(1.5230491) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6054305) q[0];
sx q[0];
rz(-2.3443065) q[0];
sx q[0];
rz(-1.5553586) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5184143) q[2];
sx q[2];
rz(-2.0668536) q[2];
sx q[2];
rz(-0.11962275) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0761145) q[1];
sx q[1];
rz(-1.5468742) q[1];
sx q[1];
rz(0.33965276) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4632439) q[3];
sx q[3];
rz(-0.58379025) q[3];
sx q[3];
rz(1.3046169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.30274621) q[2];
sx q[2];
rz(-1.6281717) q[2];
sx q[2];
rz(-1.6839074) q[2];
rz(-0.95528209) q[3];
sx q[3];
rz(-2.6421319) q[3];
sx q[3];
rz(-2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7569698) q[0];
sx q[0];
rz(-2.5769825) q[0];
sx q[0];
rz(2.6589174) q[0];
rz(-2.2118498) q[1];
sx q[1];
rz(-1.0044121) q[1];
sx q[1];
rz(0.39628705) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.020655) q[0];
sx q[0];
rz(-1.9299986) q[0];
sx q[0];
rz(2.6976556) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70102923) q[2];
sx q[2];
rz(-2.9091479) q[2];
sx q[2];
rz(-1.6384517) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0925818) q[1];
sx q[1];
rz(-1.3900458) q[1];
sx q[1];
rz(-0.84047079) q[1];
x q[2];
rz(-2.73783) q[3];
sx q[3];
rz(-1.9039394) q[3];
sx q[3];
rz(3.1333095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79260176) q[2];
sx q[2];
rz(-1.7020117) q[2];
sx q[2];
rz(-0.3271884) q[2];
rz(-2.3606825) q[3];
sx q[3];
rz(-1.1612929) q[3];
sx q[3];
rz(-1.4769295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1031716) q[0];
sx q[0];
rz(-0.99089834) q[0];
sx q[0];
rz(0.25767576) q[0];
rz(-2.7720263) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(0.69330458) q[2];
sx q[2];
rz(-1.605576) q[2];
sx q[2];
rz(-1.3645542) q[2];
rz(1.0450324) q[3];
sx q[3];
rz(-0.67787328) q[3];
sx q[3];
rz(-1.9280435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
