OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7810516) q[0];
sx q[0];
rz(3.7563503) q[0];
sx q[0];
rz(10.496245) q[0];
rz(2.7331424) q[1];
sx q[1];
rz(-0.93745679) q[1];
sx q[1];
rz(-1.6931005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14872257) q[0];
sx q[0];
rz(-2.5186484) q[0];
sx q[0];
rz(2.2882266) q[0];
rz(-1.9117457) q[2];
sx q[2];
rz(-2.1181137) q[2];
sx q[2];
rz(2.7173619) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.422157) q[1];
sx q[1];
rz(-1.4220752) q[1];
sx q[1];
rz(1.4805111) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86032773) q[3];
sx q[3];
rz(-1.9091879) q[3];
sx q[3];
rz(2.8380054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94413269) q[2];
sx q[2];
rz(-2.2977273) q[2];
sx q[2];
rz(-2.3485363) q[2];
rz(2.872725) q[3];
sx q[3];
rz(-1.3258508) q[3];
sx q[3];
rz(-0.9602921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31084138) q[0];
sx q[0];
rz(-0.38250592) q[0];
sx q[0];
rz(0.88019669) q[0];
rz(2.510732) q[1];
sx q[1];
rz(-0.81268251) q[1];
sx q[1];
rz(1.0989443) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99115935) q[0];
sx q[0];
rz(-2.9637103) q[0];
sx q[0];
rz(-3.1129897) q[0];
rz(-2.9756594) q[2];
sx q[2];
rz(-2.182462) q[2];
sx q[2];
rz(2.2076904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0365043) q[1];
sx q[1];
rz(-1.1394115) q[1];
sx q[1];
rz(-2.9626767) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1707525) q[3];
sx q[3];
rz(-2.4106541) q[3];
sx q[3];
rz(-0.39492861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38689303) q[2];
sx q[2];
rz(-1.6563481) q[2];
sx q[2];
rz(-2.5999542) q[2];
rz(-0.70332876) q[3];
sx q[3];
rz(-0.43916217) q[3];
sx q[3];
rz(0.37731236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90848732) q[0];
sx q[0];
rz(-2.2623514) q[0];
sx q[0];
rz(0.096906699) q[0];
rz(0.58755177) q[1];
sx q[1];
rz(-0.89043003) q[1];
sx q[1];
rz(2.5403835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9158949) q[0];
sx q[0];
rz(-0.89381274) q[0];
sx q[0];
rz(-1.3381132) q[0];
rz(-0.16232441) q[2];
sx q[2];
rz(-1.6337724) q[2];
sx q[2];
rz(1.2895101) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.435531) q[1];
sx q[1];
rz(-1.5630504) q[1];
sx q[1];
rz(-2.6813497) q[1];
x q[2];
rz(2.2279068) q[3];
sx q[3];
rz(-1.2877712) q[3];
sx q[3];
rz(-2.117827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2931557) q[2];
sx q[2];
rz(-1.6907254) q[2];
sx q[2];
rz(-0.71451521) q[2];
rz(-1.4826639) q[3];
sx q[3];
rz(-2.8667993) q[3];
sx q[3];
rz(-0.37884918) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191384) q[0];
sx q[0];
rz(-1.9062573) q[0];
sx q[0];
rz(-1.6374913) q[0];
rz(-2.0488886) q[1];
sx q[1];
rz(-1.8605109) q[1];
sx q[1];
rz(2.5561996) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71676895) q[0];
sx q[0];
rz(-2.1930074) q[0];
sx q[0];
rz(1.2347925) q[0];
rz(3.131065) q[2];
sx q[2];
rz(-1.4665571) q[2];
sx q[2];
rz(-0.57360211) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55951872) q[1];
sx q[1];
rz(-2.2041956) q[1];
sx q[1];
rz(2.5410209) q[1];
rz(-2.500072) q[3];
sx q[3];
rz(-1.975394) q[3];
sx q[3];
rz(-1.4610964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66712159) q[2];
sx q[2];
rz(-0.67255628) q[2];
sx q[2];
rz(2.7453864) q[2];
rz(-0.19043663) q[3];
sx q[3];
rz(-0.93947828) q[3];
sx q[3];
rz(2.6207391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025295479) q[0];
sx q[0];
rz(-2.681356) q[0];
sx q[0];
rz(0.39475557) q[0];
rz(-1.0074298) q[1];
sx q[1];
rz(-1.9837244) q[1];
sx q[1];
rz(1.6068858) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7459864) q[0];
sx q[0];
rz(-1.853038) q[0];
sx q[0];
rz(-0.76947468) q[0];
rz(-pi) q[1];
rz(1.9362279) q[2];
sx q[2];
rz(-0.88813587) q[2];
sx q[2];
rz(-2.3548369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4940961) q[1];
sx q[1];
rz(-2.7733884) q[1];
sx q[1];
rz(-1.0009264) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0232361) q[3];
sx q[3];
rz(-1.1471738) q[3];
sx q[3];
rz(2.771559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4452303) q[2];
sx q[2];
rz(-2.4424489) q[2];
sx q[2];
rz(-0.17808476) q[2];
rz(2.1961424) q[3];
sx q[3];
rz(-2.8036696) q[3];
sx q[3];
rz(0.43615714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6984542) q[0];
sx q[0];
rz(-2.4789424) q[0];
sx q[0];
rz(0.15289256) q[0];
rz(1.0180417) q[1];
sx q[1];
rz(-0.94296229) q[1];
sx q[1];
rz(2.5315888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40326443) q[0];
sx q[0];
rz(-0.77304196) q[0];
sx q[0];
rz(-2.9570262) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1334396) q[2];
sx q[2];
rz(-0.9969396) q[2];
sx q[2];
rz(-1.2795841) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.568012) q[1];
sx q[1];
rz(-1.8046772) q[1];
sx q[1];
rz(-2.3877445) q[1];
x q[2];
rz(-2.445137) q[3];
sx q[3];
rz(-3.0142733) q[3];
sx q[3];
rz(2.6501973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4623922) q[2];
sx q[2];
rz(-1.2827164) q[2];
sx q[2];
rz(-2.882615) q[2];
rz(2.9955043) q[3];
sx q[3];
rz(-0.77272213) q[3];
sx q[3];
rz(-2.5025388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884969) q[0];
sx q[0];
rz(-0.86576068) q[0];
sx q[0];
rz(-0.75188941) q[0];
rz(-0.73196661) q[1];
sx q[1];
rz(-1.3804133) q[1];
sx q[1];
rz(0.52925777) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89443356) q[0];
sx q[0];
rz(-0.34810796) q[0];
sx q[0];
rz(-2.0624119) q[0];
rz(-2.0257607) q[2];
sx q[2];
rz(-1.6935529) q[2];
sx q[2];
rz(-1.3467186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.098245278) q[1];
sx q[1];
rz(-1.4895413) q[1];
sx q[1];
rz(3.1286865) q[1];
rz(-pi) q[2];
rz(-1.115836) q[3];
sx q[3];
rz(-1.2833724) q[3];
sx q[3];
rz(-1.5446203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8818714) q[2];
sx q[2];
rz(-1.0292116) q[2];
sx q[2];
rz(-0.79505801) q[2];
rz(1.2231539) q[3];
sx q[3];
rz(-1.7984248) q[3];
sx q[3];
rz(-2.3651626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.935598) q[0];
sx q[0];
rz(-1.5482276) q[0];
sx q[0];
rz(-3.026631) q[0];
rz(0.089275442) q[1];
sx q[1];
rz(-2.142579) q[1];
sx q[1];
rz(-2.9866536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1316589) q[0];
sx q[0];
rz(-1.7007052) q[0];
sx q[0];
rz(-3.0402501) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2346588) q[2];
sx q[2];
rz(-2.0217381) q[2];
sx q[2];
rz(-0.622657) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.36679249) q[1];
sx q[1];
rz(-1.8272562) q[1];
sx q[1];
rz(2.7406997) q[1];
rz(-pi) q[2];
rz(2.2959034) q[3];
sx q[3];
rz(-1.6250623) q[3];
sx q[3];
rz(-2.3783663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0838919) q[2];
sx q[2];
rz(-1.8996779) q[2];
sx q[2];
rz(-2.4532301) q[2];
rz(-0.56811959) q[3];
sx q[3];
rz(-1.0024242) q[3];
sx q[3];
rz(0.68305558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883009) q[0];
sx q[0];
rz(-0.81128565) q[0];
sx q[0];
rz(-1.2055093) q[0];
rz(-1.0909117) q[1];
sx q[1];
rz(-2.5542407) q[1];
sx q[1];
rz(-1.6039414) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086370416) q[0];
sx q[0];
rz(-2.0049262) q[0];
sx q[0];
rz(2.7403797) q[0];
rz(-pi) q[1];
rz(1.0028935) q[2];
sx q[2];
rz(-1.5712103) q[2];
sx q[2];
rz(1.7985207) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1343775) q[1];
sx q[1];
rz(-0.79979169) q[1];
sx q[1];
rz(-0.15534955) q[1];
x q[2];
rz(1.0991715) q[3];
sx q[3];
rz(-1.0256919) q[3];
sx q[3];
rz(2.4555912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9566112) q[2];
sx q[2];
rz(-2.0272171) q[2];
sx q[2];
rz(1.9967009) q[2];
rz(1.6635118) q[3];
sx q[3];
rz(-0.76668113) q[3];
sx q[3];
rz(-2.0293106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0803273) q[0];
sx q[0];
rz(-0.63617951) q[0];
sx q[0];
rz(-1.7058477) q[0];
rz(1.1371293) q[1];
sx q[1];
rz(-0.69157332) q[1];
sx q[1];
rz(2.2045076) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6244753) q[0];
sx q[0];
rz(-1.269425) q[0];
sx q[0];
rz(-1.4791489) q[0];
x q[1];
rz(-0.50855277) q[2];
sx q[2];
rz(-1.1932696) q[2];
sx q[2];
rz(3.0005531) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4623066) q[1];
sx q[1];
rz(-0.7711773) q[1];
sx q[1];
rz(-0.036840082) q[1];
rz(-pi) q[2];
rz(-0.4611909) q[3];
sx q[3];
rz(-0.66890016) q[3];
sx q[3];
rz(1.2707641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.81637853) q[2];
sx q[2];
rz(-2.4983695) q[2];
sx q[2];
rz(2.8312259) q[2];
rz(-0.54272932) q[3];
sx q[3];
rz(-2.2118745) q[3];
sx q[3];
rz(-2.824596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19001374) q[0];
sx q[0];
rz(-1.0931451) q[0];
sx q[0];
rz(-2.6059294) q[0];
rz(-0.71470064) q[1];
sx q[1];
rz(-1.8938046) q[1];
sx q[1];
rz(-1.6751777) q[1];
rz(2.428029) q[2];
sx q[2];
rz(-0.55222558) q[2];
sx q[2];
rz(-1.4520558) q[2];
rz(1.6313995) q[3];
sx q[3];
rz(-2.2987859) q[3];
sx q[3];
rz(1.4684341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
