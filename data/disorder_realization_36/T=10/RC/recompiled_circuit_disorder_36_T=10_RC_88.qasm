OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(-2.6262992) q[0];
rz(4.4858785) q[1];
sx q[1];
rz(2.9872515) q[1];
sx q[1];
rz(6.8607688) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.261895) q[0];
sx q[0];
rz(-1.4059773) q[0];
sx q[0];
rz(0.66189712) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0787813) q[2];
sx q[2];
rz(-1.0928109) q[2];
sx q[2];
rz(-0.21131549) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9211728) q[1];
sx q[1];
rz(-2.6487659) q[1];
sx q[1];
rz(2.1652031) q[1];
rz(-0.9382117) q[3];
sx q[3];
rz(-1.0217561) q[3];
sx q[3];
rz(3.0905746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43705964) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(-0.68721592) q[2];
rz(-2.1263188) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17094831) q[0];
sx q[0];
rz(-1.0785372) q[0];
sx q[0];
rz(1.8815536) q[0];
rz(-1.0062224) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(2.2959183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2226919) q[0];
sx q[0];
rz(-1.0534304) q[0];
sx q[0];
rz(2.0738686) q[0];
rz(0.13375608) q[2];
sx q[2];
rz(-1.4417366) q[2];
sx q[2];
rz(1.918902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3413275) q[1];
sx q[1];
rz(-2.6086573) q[1];
sx q[1];
rz(2.76782) q[1];
x q[2];
rz(-0.91652292) q[3];
sx q[3];
rz(-2.407981) q[3];
sx q[3];
rz(-0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(2.6611924) q[2];
rz(-1.7885615) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(1.9539179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1903494) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-2.0551596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.438293) q[0];
sx q[0];
rz(-2.6801531) q[0];
sx q[0];
rz(-1.567054) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2181182) q[2];
sx q[2];
rz(-1.4787276) q[2];
sx q[2];
rz(1.1084686) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.8311027) q[1];
sx q[1];
rz(-1.576014) q[1];
sx q[1];
rz(0.28880854) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17351405) q[3];
sx q[3];
rz(-1.1713235) q[3];
sx q[3];
rz(-0.80843335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(1.3712937) q[2];
rz(-2.7584372) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(-0.048359811) q[0];
rz(0.16391779) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(1.6960467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018054124) q[0];
sx q[0];
rz(-2.2211383) q[0];
sx q[0];
rz(-1.3440078) q[0];
x q[1];
rz(0.21638685) q[2];
sx q[2];
rz(-1.4828223) q[2];
sx q[2];
rz(-0.83425922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2913937) q[1];
sx q[1];
rz(-1.7899917) q[1];
sx q[1];
rz(1.6367957) q[1];
x q[2];
rz(1.0501782) q[3];
sx q[3];
rz(-1.8225267) q[3];
sx q[3];
rz(2.7057735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(-1.6131489) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4399453) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.2874999) q[0];
rz(-0.31907407) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(-2.2873926) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720649) q[0];
sx q[0];
rz(-0.85692353) q[0];
sx q[0];
rz(-3.1307427) q[0];
x q[1];
rz(-0.87128432) q[2];
sx q[2];
rz(-1.185002) q[2];
sx q[2];
rz(0.7427578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0336049) q[1];
sx q[1];
rz(-1.0228979) q[1];
sx q[1];
rz(-2.4461436) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0075931) q[3];
sx q[3];
rz(-1.0922722) q[3];
sx q[3];
rz(1.0790881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.95191082) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(-2.5642776) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(2.234941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1381056) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(-0.072120897) q[0];
rz(1.1068608) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(3.0153826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.135658) q[0];
sx q[0];
rz(-1.5392443) q[0];
sx q[0];
rz(-2.9289782) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77504471) q[2];
sx q[2];
rz(-1.2942874) q[2];
sx q[2];
rz(-1.0724049) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0446484) q[1];
sx q[1];
rz(-1.7587887) q[1];
sx q[1];
rz(2.2915927) q[1];
rz(-pi) q[2];
rz(1.0522271) q[3];
sx q[3];
rz(-0.37399451) q[3];
sx q[3];
rz(-0.75043375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90298992) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(0.827968) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(-2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5489952) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(0.098408498) q[0];
rz(1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-0.55955204) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7647117) q[0];
sx q[0];
rz(-1.6983713) q[0];
sx q[0];
rz(-0.60140951) q[0];
rz(-pi) q[1];
rz(-0.26142188) q[2];
sx q[2];
rz(-1.6209507) q[2];
sx q[2];
rz(0.23471552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2110062) q[1];
sx q[1];
rz(-1.6582656) q[1];
sx q[1];
rz(-1.7006111) q[1];
rz(-1.7667174) q[3];
sx q[3];
rz(-2.19176) q[3];
sx q[3];
rz(1.2697112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.45903912) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(0.34379488) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664117) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-2.4108316) q[0];
rz(-2.9991951) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(0.87160814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14411892) q[0];
sx q[0];
rz(-1.5796356) q[0];
sx q[0];
rz(-3.0661422) q[0];
rz(0.34198728) q[2];
sx q[2];
rz(-2.7139398) q[2];
sx q[2];
rz(-0.74117408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.826088) q[1];
sx q[1];
rz(-1.2102038) q[1];
sx q[1];
rz(-2.4396067) q[1];
rz(-2.1731247) q[3];
sx q[3];
rz(-1.6864711) q[3];
sx q[3];
rz(2.065421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(-2.0020206) q[2];
rz(-1.6428044) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(-2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20294872) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(1.2217481) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(1.6171914) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4964543) q[0];
sx q[0];
rz(-1.5997412) q[0];
sx q[0];
rz(-1.653198) q[0];
x q[1];
rz(0.35877123) q[2];
sx q[2];
rz(-2.6371187) q[2];
sx q[2];
rz(-0.82746738) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9164239) q[1];
sx q[1];
rz(-1.4125707) q[1];
sx q[1];
rz(1.9041063) q[1];
rz(-pi) q[2];
rz(-0.63906007) q[3];
sx q[3];
rz(-1.3165054) q[3];
sx q[3];
rz(1.1921079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.021585492) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(2.7900556) q[2];
rz(1.0567788) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4979424) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(-1.836401) q[0];
rz(0.38048831) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(0.25340733) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5288552) q[0];
sx q[0];
rz(-0.40898541) q[0];
sx q[0];
rz(-1.8324052) q[0];
rz(2.9741653) q[2];
sx q[2];
rz(-1.2360459) q[2];
sx q[2];
rz(1.6850922) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0903783) q[1];
sx q[1];
rz(-1.2284632) q[1];
sx q[1];
rz(1.315209) q[1];
x q[2];
rz(2.5305383) q[3];
sx q[3];
rz(-1.9285893) q[3];
sx q[3];
rz(-0.65294453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(2.6386476) q[2];
rz(2.2425966) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(-1.9780654) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1702561) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(-2.4304216) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(-0.74264991) q[2];
sx q[2];
rz(-0.37336083) q[2];
sx q[2];
rz(0.20867418) q[2];
rz(2.3861804) q[3];
sx q[3];
rz(-2.073954) q[3];
sx q[3];
rz(3.0975773) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
