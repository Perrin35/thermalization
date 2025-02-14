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
rz(-2.8673515) q[0];
sx q[0];
rz(-0.40617901) q[0];
sx q[0];
rz(2.9131373) q[0];
rz(0.38415456) q[1];
sx q[1];
rz(1.6698807) q[1];
sx q[1];
rz(8.9759965) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2160546) q[0];
sx q[0];
rz(-1.9384697) q[0];
sx q[0];
rz(1.4212763) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8443135) q[2];
sx q[2];
rz(-2.7113072) q[2];
sx q[2];
rz(1.6543433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0579867) q[1];
sx q[1];
rz(-0.068002105) q[1];
sx q[1];
rz(1.9307095) q[1];
x q[2];
rz(3.1136127) q[3];
sx q[3];
rz(-1.5767808) q[3];
sx q[3];
rz(0.75388349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3936798) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61504859) q[0];
sx q[0];
rz(-2.5874309) q[0];
sx q[0];
rz(0.91914415) q[0];
rz(2.9650086) q[1];
sx q[1];
rz(-1.6932026) q[1];
sx q[1];
rz(0.44763705) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31792163) q[0];
sx q[0];
rz(-1.5783675) q[0];
sx q[0];
rz(1.5764144) q[0];
x q[1];
rz(-2.7451374) q[2];
sx q[2];
rz(-1.8282991) q[2];
sx q[2];
rz(-0.56397179) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3474975) q[1];
sx q[1];
rz(-1.4117876) q[1];
sx q[1];
rz(-0.68557941) q[1];
rz(-2.3354946) q[3];
sx q[3];
rz(-2.4662123) q[3];
sx q[3];
rz(0.71682978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24719079) q[2];
sx q[2];
rz(-1.2195419) q[2];
sx q[2];
rz(2.2578237) q[2];
rz(0.73924685) q[3];
sx q[3];
rz(-1.913162) q[3];
sx q[3];
rz(1.242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7582551) q[0];
sx q[0];
rz(-2.1012335) q[0];
sx q[0];
rz(3.1357159) q[0];
rz(0.30644304) q[1];
sx q[1];
rz(-1.076661) q[1];
sx q[1];
rz(-1.2642911) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41691859) q[0];
sx q[0];
rz(-1.8643799) q[0];
sx q[0];
rz(2.7489566) q[0];
rz(-pi) q[1];
rz(-0.034175506) q[2];
sx q[2];
rz(-1.719599) q[2];
sx q[2];
rz(0.11916313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2325461) q[1];
sx q[1];
rz(-1.0139483) q[1];
sx q[1];
rz(1.805748) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39671582) q[3];
sx q[3];
rz(-1.0457888) q[3];
sx q[3];
rz(0.58806149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82187247) q[2];
sx q[2];
rz(-1.7629273) q[2];
sx q[2];
rz(2.482448) q[2];
rz(-1.2705196) q[3];
sx q[3];
rz(-0.29265413) q[3];
sx q[3];
rz(-1.6856153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6419411) q[0];
sx q[0];
rz(-2.9872276) q[0];
sx q[0];
rz(0.9374215) q[0];
rz(1.5501267) q[1];
sx q[1];
rz(-0.66842404) q[1];
sx q[1];
rz(2.392427) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5431108) q[0];
sx q[0];
rz(-2.3078663) q[0];
sx q[0];
rz(2.7373256) q[0];
x q[1];
rz(-2.238041) q[2];
sx q[2];
rz(-2.090881) q[2];
sx q[2];
rz(-1.5423216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8038102) q[1];
sx q[1];
rz(-2.446108) q[1];
sx q[1];
rz(1.8799604) q[1];
rz(2.9177298) q[3];
sx q[3];
rz(-2.3626257) q[3];
sx q[3];
rz(1.1457486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5472827) q[2];
sx q[2];
rz(-2.3197033) q[2];
sx q[2];
rz(1.8853356) q[2];
rz(0.78993434) q[3];
sx q[3];
rz(-2.200685) q[3];
sx q[3];
rz(1.1676211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0867406) q[0];
sx q[0];
rz(-2.5331443) q[0];
sx q[0];
rz(-0.21859455) q[0];
rz(-2.356148) q[1];
sx q[1];
rz(-1.1050858) q[1];
sx q[1];
rz(1.2931664) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4703641) q[0];
sx q[0];
rz(-0.73743907) q[0];
sx q[0];
rz(-2.295275) q[0];
rz(-1.3082383) q[2];
sx q[2];
rz(-1.6505604) q[2];
sx q[2];
rz(0.74918109) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4689815) q[1];
sx q[1];
rz(-1.0291576) q[1];
sx q[1];
rz(2.2865613) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3985211) q[3];
sx q[3];
rz(-2.9228811) q[3];
sx q[3];
rz(2.2416473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30696067) q[2];
sx q[2];
rz(-1.0647651) q[2];
sx q[2];
rz(-1.5131697) q[2];
rz(0.73897922) q[3];
sx q[3];
rz(-1.8578953) q[3];
sx q[3];
rz(-1.6966604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6936679) q[0];
sx q[0];
rz(-1.4883) q[0];
sx q[0];
rz(-2.779575) q[0];
rz(-1.1297049) q[1];
sx q[1];
rz(-1.2672707) q[1];
sx q[1];
rz(-0.76090181) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99859649) q[0];
sx q[0];
rz(-0.79299295) q[0];
sx q[0];
rz(-2.8727813) q[0];
rz(-pi) q[1];
rz(0.40973969) q[2];
sx q[2];
rz(-0.86977661) q[2];
sx q[2];
rz(2.1363195) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89448358) q[1];
sx q[1];
rz(-1.037408) q[1];
sx q[1];
rz(-0.24786149) q[1];
rz(0.53583709) q[3];
sx q[3];
rz(-2.4501175) q[3];
sx q[3];
rz(-2.7590318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1165983) q[2];
sx q[2];
rz(-2.8572539) q[2];
sx q[2];
rz(1.9653758) q[2];
rz(2.1182012) q[3];
sx q[3];
rz(-1.0756282) q[3];
sx q[3];
rz(-2.7682176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0244331) q[0];
sx q[0];
rz(-2.7130782) q[0];
sx q[0];
rz(-1.8350711) q[0];
rz(0.10082968) q[1];
sx q[1];
rz(-1.7497302) q[1];
sx q[1];
rz(-0.94591013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8656509) q[0];
sx q[0];
rz(-1.2667333) q[0];
sx q[0];
rz(1.434554) q[0];
x q[1];
rz(-2.3657799) q[2];
sx q[2];
rz(-2.472942) q[2];
sx q[2];
rz(-2.5324627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57179831) q[1];
sx q[1];
rz(-1.6299953) q[1];
sx q[1];
rz(-3.0059391) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9705371) q[3];
sx q[3];
rz(-0.75467602) q[3];
sx q[3];
rz(0.012059742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0239608) q[2];
sx q[2];
rz(-1.5317711) q[2];
sx q[2];
rz(-2.7839933) q[2];
rz(-3.0910953) q[3];
sx q[3];
rz(-2.7183618) q[3];
sx q[3];
rz(-2.7549506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9015775) q[0];
sx q[0];
rz(-1.3192663) q[0];
sx q[0];
rz(-0.59558076) q[0];
rz(1.2254084) q[1];
sx q[1];
rz(-1.3662246) q[1];
sx q[1];
rz(-0.43232408) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5670938) q[0];
sx q[0];
rz(-1.7364952) q[0];
sx q[0];
rz(-0.49543799) q[0];
rz(-1.4171227) q[2];
sx q[2];
rz(-1.5072848) q[2];
sx q[2];
rz(-2.2907298) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9185168) q[1];
sx q[1];
rz(-2.1380414) q[1];
sx q[1];
rz(1.569313) q[1];
x q[2];
rz(-0.8882167) q[3];
sx q[3];
rz(-1.9239263) q[3];
sx q[3];
rz(-0.61643657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7633742) q[2];
sx q[2];
rz(-1.8209527) q[2];
sx q[2];
rz(2.1136368) q[2];
rz(1.6592735) q[3];
sx q[3];
rz(-1.6927398) q[3];
sx q[3];
rz(-0.75215522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9529652) q[0];
sx q[0];
rz(-1.9099706) q[0];
sx q[0];
rz(-2.6787483) q[0];
rz(-1.2340744) q[1];
sx q[1];
rz(-1.2214829) q[1];
sx q[1];
rz(-5*pi/8) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45904392) q[0];
sx q[0];
rz(-0.17691806) q[0];
sx q[0];
rz(1.7155619) q[0];
rz(-pi) q[1];
rz(-2.0527564) q[2];
sx q[2];
rz(-2.3079527) q[2];
sx q[2];
rz(2.9046975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.000899) q[1];
sx q[1];
rz(-1.447666) q[1];
sx q[1];
rz(-0.85832125) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7998802) q[3];
sx q[3];
rz(-1.4572506) q[3];
sx q[3];
rz(-2.7734203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.50596333) q[2];
sx q[2];
rz(-1.3869945) q[2];
sx q[2];
rz(-0.76816922) q[2];
rz(0.65028894) q[3];
sx q[3];
rz(-0.30935973) q[3];
sx q[3];
rz(-2.5816141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3637417) q[0];
sx q[0];
rz(-0.44054458) q[0];
sx q[0];
rz(1.1458696) q[0];
rz(2.5762985) q[1];
sx q[1];
rz(-2.4259613) q[1];
sx q[1];
rz(-0.40564767) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89116468) q[0];
sx q[0];
rz(-1.9180204) q[0];
sx q[0];
rz(2.3655212) q[0];
rz(-pi) q[1];
rz(-0.27743199) q[2];
sx q[2];
rz(-1.9346969) q[2];
sx q[2];
rz(-2.8958447) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.010505761) q[1];
sx q[1];
rz(-1.0684135) q[1];
sx q[1];
rz(1.8598063) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9156014) q[3];
sx q[3];
rz(-0.44934326) q[3];
sx q[3];
rz(-2.9087272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.046935) q[2];
sx q[2];
rz(-1.2311225) q[2];
sx q[2];
rz(0.51935736) q[2];
rz(-2.6804067) q[3];
sx q[3];
rz(-1.0591732) q[3];
sx q[3];
rz(-2.4477899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8984579) q[0];
sx q[0];
rz(-1.5512137) q[0];
sx q[0];
rz(2.2358538) q[0];
rz(1.657919) q[1];
sx q[1];
rz(-1.9762194) q[1];
sx q[1];
rz(-1.2414052) q[1];
rz(-2.1853191) q[2];
sx q[2];
rz(-1.4531789) q[2];
sx q[2];
rz(1.5868759) q[2];
rz(-1.0736856) q[3];
sx q[3];
rz(-1.9852255) q[3];
sx q[3];
rz(-1.6559813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
