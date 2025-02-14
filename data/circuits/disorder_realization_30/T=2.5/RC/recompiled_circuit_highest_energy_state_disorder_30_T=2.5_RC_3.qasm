OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27424115) q[0];
sx q[0];
rz(-2.7354136) q[0];
sx q[0];
rz(0.22845536) q[0];
rz(-2.7574381) q[1];
sx q[1];
rz(-1.6698807) q[1];
sx q[1];
rz(2.6928112) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3224027) q[0];
sx q[0];
rz(-2.7459641) q[0];
sx q[0];
rz(2.7725793) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2972791) q[2];
sx q[2];
rz(-0.43028545) q[2];
sx q[2];
rz(1.6543433) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8463414) q[1];
sx q[1];
rz(-1.59473) q[1];
sx q[1];
rz(1.5071391) q[1];
x q[2];
rz(-0.21073802) q[3];
sx q[3];
rz(-3.11298) q[3];
sx q[3];
rz(-0.60625854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7479129) q[2];
sx q[2];
rz(-2.9069275) q[2];
sx q[2];
rz(2.5133384) q[2];
rz(2.0755532) q[3];
sx q[3];
rz(-0.99848905) q[3];
sx q[3];
rz(2.0121241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5265441) q[0];
sx q[0];
rz(-0.55416179) q[0];
sx q[0];
rz(0.91914415) q[0];
rz(-2.9650086) q[1];
sx q[1];
rz(-1.44839) q[1];
sx q[1];
rz(-2.6939556) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95631188) q[0];
sx q[0];
rz(-3.1321648) q[0];
sx q[0];
rz(2.5032237) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39645529) q[2];
sx q[2];
rz(-1.8282991) q[2];
sx q[2];
rz(-0.56397179) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3474975) q[1];
sx q[1];
rz(-1.4117876) q[1];
sx q[1];
rz(2.4560132) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3354946) q[3];
sx q[3];
rz(-0.67538031) q[3];
sx q[3];
rz(-0.71682978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8944019) q[2];
sx q[2];
rz(-1.2195419) q[2];
sx q[2];
rz(-0.88376898) q[2];
rz(2.4023458) q[3];
sx q[3];
rz(-1.2284307) q[3];
sx q[3];
rz(-1.8993186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3833375) q[0];
sx q[0];
rz(-2.1012335) q[0];
sx q[0];
rz(-3.1357159) q[0];
rz(-2.8351496) q[1];
sx q[1];
rz(-1.076661) q[1];
sx q[1];
rz(1.8773016) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41691859) q[0];
sx q[0];
rz(-1.8643799) q[0];
sx q[0];
rz(0.39263607) q[0];
x q[1];
rz(3.1074171) q[2];
sx q[2];
rz(-1.4219936) q[2];
sx q[2];
rz(3.0224295) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6578232) q[1];
sx q[1];
rz(-0.59952901) q[1];
sx q[1];
rz(-0.35783135) q[1];
rz(-pi) q[2];
rz(2.131553) q[3];
sx q[3];
rz(-1.2298785) q[3];
sx q[3];
rz(-0.775767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82187247) q[2];
sx q[2];
rz(-1.7629273) q[2];
sx q[2];
rz(2.482448) q[2];
rz(1.871073) q[3];
sx q[3];
rz(-0.29265413) q[3];
sx q[3];
rz(-1.6856153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4996516) q[0];
sx q[0];
rz(-0.15436509) q[0];
sx q[0];
rz(2.2041712) q[0];
rz(1.591466) q[1];
sx q[1];
rz(-2.4731686) q[1];
sx q[1];
rz(2.392427) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1099656) q[0];
sx q[0];
rz(-0.82200371) q[0];
sx q[0];
rz(1.9796994) q[0];
rz(-2.5116483) q[2];
sx q[2];
rz(-1.0039181) q[2];
sx q[2];
rz(0.34467372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3377825) q[1];
sx q[1];
rz(-0.69548464) q[1];
sx q[1];
rz(1.2616322) q[1];
rz(-pi) q[2];
x q[2];
rz(2.375256) q[3];
sx q[3];
rz(-1.4141937) q[3];
sx q[3];
rz(-2.5559154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5472827) q[2];
sx q[2];
rz(-0.82188934) q[2];
sx q[2];
rz(1.8853356) q[2];
rz(2.3516583) q[3];
sx q[3];
rz(-0.94090763) q[3];
sx q[3];
rz(1.1676211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0867406) q[0];
sx q[0];
rz(-0.60844839) q[0];
sx q[0];
rz(2.9229981) q[0];
rz(-0.78544468) q[1];
sx q[1];
rz(-1.1050858) q[1];
sx q[1];
rz(-1.2931664) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6712286) q[0];
sx q[0];
rz(-0.73743907) q[0];
sx q[0];
rz(2.295275) q[0];
rz(-1.3082383) q[2];
sx q[2];
rz(-1.6505604) q[2];
sx q[2];
rz(-2.3924116) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.67261117) q[1];
sx q[1];
rz(-2.1124351) q[1];
sx q[1];
rz(-0.85503135) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.038083548) q[3];
sx q[3];
rz(-1.3553737) q[3];
sx q[3];
rz(2.0652536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.834632) q[2];
sx q[2];
rz(-2.0768276) q[2];
sx q[2];
rz(-1.628423) q[2];
rz(-2.4026134) q[3];
sx q[3];
rz(-1.2836974) q[3];
sx q[3];
rz(1.6966604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6936679) q[0];
sx q[0];
rz(-1.4883) q[0];
sx q[0];
rz(-0.36201763) q[0];
rz(-1.1297049) q[1];
sx q[1];
rz(-1.2672707) q[1];
sx q[1];
rz(2.3806908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7689038) q[0];
sx q[0];
rz(-2.3279705) q[0];
sx q[0];
rz(1.3074101) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40973969) q[2];
sx q[2];
rz(-2.271816) q[2];
sx q[2];
rz(-1.0052731) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.43266549) q[1];
sx q[1];
rz(-2.558518) q[1];
sx q[1];
rz(-1.1770171) q[1];
rz(-2.6057556) q[3];
sx q[3];
rz(-0.69147516) q[3];
sx q[3];
rz(2.7590318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1165983) q[2];
sx q[2];
rz(-2.8572539) q[2];
sx q[2];
rz(-1.9653758) q[2];
rz(1.0233915) q[3];
sx q[3];
rz(-2.0659645) q[3];
sx q[3];
rz(0.37337506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1171595) q[0];
sx q[0];
rz(-0.42851448) q[0];
sx q[0];
rz(-1.8350711) q[0];
rz(3.040763) q[1];
sx q[1];
rz(-1.3918624) q[1];
sx q[1];
rz(-0.94591013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8656509) q[0];
sx q[0];
rz(-1.8748594) q[0];
sx q[0];
rz(1.7070387) q[0];
rz(0.51351764) q[2];
sx q[2];
rz(-1.1217204) q[2];
sx q[2];
rz(-2.8359536) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.57179831) q[1];
sx q[1];
rz(-1.6299953) q[1];
sx q[1];
rz(-3.0059391) q[1];
x q[2];
rz(-1.1710555) q[3];
sx q[3];
rz(-0.75467602) q[3];
sx q[3];
rz(3.1295329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1176318) q[2];
sx q[2];
rz(-1.6098216) q[2];
sx q[2];
rz(-0.35759932) q[2];
rz(-0.050497342) q[3];
sx q[3];
rz(-2.7183618) q[3];
sx q[3];
rz(2.7549506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24001515) q[0];
sx q[0];
rz(-1.3192663) q[0];
sx q[0];
rz(0.59558076) q[0];
rz(-1.9161842) q[1];
sx q[1];
rz(-1.775368) q[1];
sx q[1];
rz(0.43232408) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5744988) q[0];
sx q[0];
rz(-1.7364952) q[0];
sx q[0];
rz(-0.49543799) q[0];
rz(-pi) q[1];
rz(0.064266845) q[2];
sx q[2];
rz(-1.7241577) q[2];
sx q[2];
rz(-0.72976412) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9185168) q[1];
sx q[1];
rz(-1.0035512) q[1];
sx q[1];
rz(-1.5722797) q[1];
x q[2];
rz(2.6981399) q[3];
sx q[3];
rz(-0.93741527) q[3];
sx q[3];
rz(-1.2284281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7633742) q[2];
sx q[2];
rz(-1.32064) q[2];
sx q[2];
rz(-1.0279559) q[2];
rz(-1.6592735) q[3];
sx q[3];
rz(-1.4488528) q[3];
sx q[3];
rz(-0.75215522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1886275) q[0];
sx q[0];
rz(-1.231622) q[0];
sx q[0];
rz(2.6787483) q[0];
rz(-1.2340744) q[1];
sx q[1];
rz(-1.2214829) q[1];
sx q[1];
rz(3*pi/8) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5355204) q[0];
sx q[0];
rz(-1.7458445) q[0];
sx q[0];
rz(-0.025786215) q[0];
rz(-pi) q[1];
rz(2.0527564) q[2];
sx q[2];
rz(-0.83363998) q[2];
sx q[2];
rz(-0.23689517) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14069362) q[1];
sx q[1];
rz(-1.447666) q[1];
sx q[1];
rz(-0.85832125) q[1];
rz(2.7998802) q[3];
sx q[3];
rz(-1.6843421) q[3];
sx q[3];
rz(0.36817238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6356293) q[2];
sx q[2];
rz(-1.3869945) q[2];
sx q[2];
rz(-2.3734234) q[2];
rz(0.65028894) q[3];
sx q[3];
rz(-0.30935973) q[3];
sx q[3];
rz(0.55997854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77785093) q[0];
sx q[0];
rz(-0.44054458) q[0];
sx q[0];
rz(-1.1458696) q[0];
rz(-2.5762985) q[1];
sx q[1];
rz(-2.4259613) q[1];
sx q[1];
rz(0.40564767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34590301) q[0];
sx q[0];
rz(-2.3064605) q[0];
sx q[0];
rz(2.664734) q[0];
rz(-pi) q[1];
rz(-0.27743199) q[2];
sx q[2];
rz(-1.2068958) q[2];
sx q[2];
rz(2.8958447) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5989117) q[1];
sx q[1];
rz(-0.57334954) q[1];
sx q[1];
rz(-0.47853761) q[1];
x q[2];
rz(-0.43934699) q[3];
sx q[3];
rz(-1.6682819) q[3];
sx q[3];
rz(-2.0078703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0946577) q[2];
sx q[2];
rz(-1.2311225) q[2];
sx q[2];
rz(-2.6222353) q[2];
rz(2.6804067) q[3];
sx q[3];
rz(-2.0824194) q[3];
sx q[3];
rz(0.69380277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2431348) q[0];
sx q[0];
rz(-1.5903789) q[0];
sx q[0];
rz(-0.90573885) q[0];
rz(1.657919) q[1];
sx q[1];
rz(-1.9762194) q[1];
sx q[1];
rz(-1.2414052) q[1];
rz(2.1853191) q[2];
sx q[2];
rz(-1.6884138) q[2];
sx q[2];
rz(-1.5547167) q[2];
rz(-2.6775581) q[3];
sx q[3];
rz(-1.1190718) q[3];
sx q[3];
rz(-3.0116871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
