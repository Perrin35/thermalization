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
rz(0.089224815) q[0];
sx q[0];
rz(-0.74639809) q[0];
sx q[0];
rz(0.69190061) q[0];
rz(0.46299419) q[1];
sx q[1];
rz(-1.6166592) q[1];
sx q[1];
rz(2.1114299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2041516) q[0];
sx q[0];
rz(-0.61249706) q[0];
sx q[0];
rz(0.4093616) q[0];
x q[1];
rz(0.18333034) q[2];
sx q[2];
rz(-1.9631443) q[2];
sx q[2];
rz(3.0073056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9249484) q[1];
sx q[1];
rz(-2.7949667) q[1];
sx q[1];
rz(-1.2216755) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6923156) q[3];
sx q[3];
rz(-2.2560824) q[3];
sx q[3];
rz(-2.6916162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87972155) q[2];
sx q[2];
rz(-0.77360952) q[2];
sx q[2];
rz(-0.36641463) q[2];
rz(2.053818) q[3];
sx q[3];
rz(-0.78947869) q[3];
sx q[3];
rz(0.68939775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.825603) q[0];
sx q[0];
rz(-3.0593384) q[0];
sx q[0];
rz(-0.91973037) q[0];
rz(-2.3986744) q[1];
sx q[1];
rz(-1.8835604) q[1];
sx q[1];
rz(-1.119335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9237083) q[0];
sx q[0];
rz(-0.15734921) q[0];
sx q[0];
rz(1.1697392) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3886003) q[2];
sx q[2];
rz(-2.2815653) q[2];
sx q[2];
rz(-1.7882535) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85038919) q[1];
sx q[1];
rz(-0.77050401) q[1];
sx q[1];
rz(0.34363665) q[1];
rz(0.51471424) q[3];
sx q[3];
rz(-1.052894) q[3];
sx q[3];
rz(1.6165773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6433581) q[2];
sx q[2];
rz(-1.9469807) q[2];
sx q[2];
rz(0.052113459) q[2];
rz(2.0335967) q[3];
sx q[3];
rz(-2.2823157) q[3];
sx q[3];
rz(0.064854709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68010083) q[0];
sx q[0];
rz(-1.2315741) q[0];
sx q[0];
rz(-1.8517866) q[0];
rz(2.3884804) q[1];
sx q[1];
rz(-1.6878004) q[1];
sx q[1];
rz(0.49711102) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4418934) q[0];
sx q[0];
rz(-1.9632247) q[0];
sx q[0];
rz(2.5470069) q[0];
rz(-pi) q[1];
rz(2.4646453) q[2];
sx q[2];
rz(-2.1008889) q[2];
sx q[2];
rz(2.9924336) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0935082) q[1];
sx q[1];
rz(-1.9447548) q[1];
sx q[1];
rz(1.611534) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0158832) q[3];
sx q[3];
rz(-2.4851228) q[3];
sx q[3];
rz(-1.8905115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46811732) q[2];
sx q[2];
rz(-2.868729) q[2];
sx q[2];
rz(-0.72431481) q[2];
rz(2.916548) q[3];
sx q[3];
rz(-1.8658172) q[3];
sx q[3];
rz(-2.3465274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1555136) q[0];
sx q[0];
rz(-2.2561095) q[0];
sx q[0];
rz(2.8610863) q[0];
rz(1.6579423) q[1];
sx q[1];
rz(-1.2018485) q[1];
sx q[1];
rz(-2.6703506) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6408117) q[0];
sx q[0];
rz(-2.7735595) q[0];
sx q[0];
rz(1.5076007) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6837101) q[2];
sx q[2];
rz(-1.5051923) q[2];
sx q[2];
rz(-1.124231) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.9556826) q[1];
sx q[1];
rz(-1.5430527) q[1];
sx q[1];
rz(1.300059) q[1];
rz(-2.332451) q[3];
sx q[3];
rz(-2.005072) q[3];
sx q[3];
rz(-0.44444042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.19807854) q[2];
sx q[2];
rz(-2.5835865) q[2];
sx q[2];
rz(-1.0151218) q[2];
rz(-2.4070168) q[3];
sx q[3];
rz(-1.7769122) q[3];
sx q[3];
rz(-0.53810292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839142) q[0];
sx q[0];
rz(-1.9240802) q[0];
sx q[0];
rz(-0.96018803) q[0];
rz(-3.0958815) q[1];
sx q[1];
rz(-2.263133) q[1];
sx q[1];
rz(0.033230573) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35141459) q[0];
sx q[0];
rz(-2.040967) q[0];
sx q[0];
rz(2.1432502) q[0];
rz(-pi) q[1];
x q[1];
rz(0.015093283) q[2];
sx q[2];
rz(-1.303267) q[2];
sx q[2];
rz(-2.044673) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.14151621) q[1];
sx q[1];
rz(-1.7011397) q[1];
sx q[1];
rz(-1.2356067) q[1];
rz(-2.700272) q[3];
sx q[3];
rz(-1.1989824) q[3];
sx q[3];
rz(-1.9723434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1632605) q[2];
sx q[2];
rz(-0.99908081) q[2];
sx q[2];
rz(-2.8283289) q[2];
rz(-2.7836109) q[3];
sx q[3];
rz(-1.0397747) q[3];
sx q[3];
rz(0.09859214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98704308) q[0];
sx q[0];
rz(-1.6313666) q[0];
sx q[0];
rz(-0.48598591) q[0];
rz(1.9588574) q[1];
sx q[1];
rz(-2.4458838) q[1];
sx q[1];
rz(-2.8114496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2898588) q[0];
sx q[0];
rz(-2.8427474) q[0];
sx q[0];
rz(-1.6720222) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5146944) q[2];
sx q[2];
rz(-0.90993249) q[2];
sx q[2];
rz(3.0783184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15974074) q[1];
sx q[1];
rz(-1.0969437) q[1];
sx q[1];
rz(0.57094806) q[1];
x q[2];
rz(-2.7948912) q[3];
sx q[3];
rz(-1.7809521) q[3];
sx q[3];
rz(-2.8488248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6351472) q[2];
sx q[2];
rz(-0.61656419) q[2];
sx q[2];
rz(0.35663024) q[2];
rz(-2.5470274) q[3];
sx q[3];
rz(-1.6584572) q[3];
sx q[3];
rz(2.4026292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67721382) q[0];
sx q[0];
rz(-0.76736275) q[0];
sx q[0];
rz(0.58694029) q[0];
rz(-2.7691973) q[1];
sx q[1];
rz(-1.3601235) q[1];
sx q[1];
rz(-0.24758235) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7149087) q[0];
sx q[0];
rz(-1.539402) q[0];
sx q[0];
rz(-2.035729) q[0];
rz(0.32682287) q[2];
sx q[2];
rz(-2.4832209) q[2];
sx q[2];
rz(3.1382552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.048760895) q[1];
sx q[1];
rz(-0.93834651) q[1];
sx q[1];
rz(-2.7238242) q[1];
rz(-pi) q[2];
rz(2.7917737) q[3];
sx q[3];
rz(-1.8050151) q[3];
sx q[3];
rz(1.8893591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26329142) q[2];
sx q[2];
rz(-1.9204488) q[2];
sx q[2];
rz(2.9325874) q[2];
rz(-0.59669295) q[3];
sx q[3];
rz(-2.7946819) q[3];
sx q[3];
rz(1.6056812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40808943) q[0];
sx q[0];
rz(-2.1812289) q[0];
sx q[0];
rz(2.7485513) q[0];
rz(-0.6706925) q[1];
sx q[1];
rz(-1.1436983) q[1];
sx q[1];
rz(1.447698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.873998) q[0];
sx q[0];
rz(-1.0769258) q[0];
sx q[0];
rz(0.14125342) q[0];
rz(-2.8301881) q[2];
sx q[2];
rz(-2.1620353) q[2];
sx q[2];
rz(-0.55252749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8116074) q[1];
sx q[1];
rz(-2.0714272) q[1];
sx q[1];
rz(-0.58793427) q[1];
rz(-1.0974698) q[3];
sx q[3];
rz(-0.83126634) q[3];
sx q[3];
rz(-0.39287469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.712901) q[2];
sx q[2];
rz(-0.8608326) q[2];
sx q[2];
rz(3.0248896) q[2];
rz(1.0137001) q[3];
sx q[3];
rz(-1.2319535) q[3];
sx q[3];
rz(1.8359756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882465) q[0];
sx q[0];
rz(-1.4660864) q[0];
sx q[0];
rz(-1.4353132) q[0];
rz(0.995579) q[1];
sx q[1];
rz(-1.8228056) q[1];
sx q[1];
rz(-2.6011655) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2453702) q[0];
sx q[0];
rz(-0.55573744) q[0];
sx q[0];
rz(-0.51387365) q[0];
x q[1];
rz(-0.58806555) q[2];
sx q[2];
rz(-2.1508475) q[2];
sx q[2];
rz(-1.5813476) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41758834) q[1];
sx q[1];
rz(-0.88694983) q[1];
sx q[1];
rz(-1.0761912) q[1];
x q[2];
rz(2.5164817) q[3];
sx q[3];
rz(-1.2338936) q[3];
sx q[3];
rz(2.5663216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.871375) q[2];
sx q[2];
rz(-1.3048708) q[2];
sx q[2];
rz(2.3649141) q[2];
rz(3.1089879) q[3];
sx q[3];
rz(-1.6083345) q[3];
sx q[3];
rz(1.5672654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8107574) q[0];
sx q[0];
rz(-2.7726655) q[0];
sx q[0];
rz(2.708013) q[0];
rz(0.66917229) q[1];
sx q[1];
rz(-1.1829665) q[1];
sx q[1];
rz(0.42632595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45352061) q[0];
sx q[0];
rz(-1.768508) q[0];
sx q[0];
rz(-1.7889678) q[0];
rz(-pi) q[1];
rz(-2.3947515) q[2];
sx q[2];
rz(-1.3051118) q[2];
sx q[2];
rz(2.8248252) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5475216) q[1];
sx q[1];
rz(-1.2202377) q[1];
sx q[1];
rz(2.622753) q[1];
x q[2];
rz(2.4678585) q[3];
sx q[3];
rz(-2.4387283) q[3];
sx q[3];
rz(-1.4386615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7508037) q[2];
sx q[2];
rz(-2.8533253) q[2];
sx q[2];
rz(-2.0354347) q[2];
rz(-0.11354167) q[3];
sx q[3];
rz(-0.93324408) q[3];
sx q[3];
rz(-0.043702628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15585598) q[0];
sx q[0];
rz(-2.2423797) q[0];
sx q[0];
rz(2.5806497) q[0];
rz(0.38324311) q[1];
sx q[1];
rz(-1.1226729) q[1];
sx q[1];
rz(2.9164006) q[1];
rz(0.55632527) q[2];
sx q[2];
rz(-1.3090324) q[2];
sx q[2];
rz(-0.39568452) q[2];
rz(2.6970277) q[3];
sx q[3];
rz(-2.2722383) q[3];
sx q[3];
rz(2.1376192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
