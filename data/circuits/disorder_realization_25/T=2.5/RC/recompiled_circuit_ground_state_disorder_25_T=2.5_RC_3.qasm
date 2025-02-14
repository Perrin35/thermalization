OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9659757) q[0];
sx q[0];
rz(-1.5953925) q[0];
sx q[0];
rz(1.6311128) q[0];
rz(-2.053396) q[1];
sx q[1];
rz(-1.8367986) q[1];
sx q[1];
rz(0.78707492) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45739503) q[0];
sx q[0];
rz(-1.4589063) q[0];
sx q[0];
rz(2.4613613) q[0];
rz(-pi) q[1];
rz(-0.12949582) q[2];
sx q[2];
rz(-0.39948764) q[2];
sx q[2];
rz(-2.0129888) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8382227) q[1];
sx q[1];
rz(-1.7821687) q[1];
sx q[1];
rz(1.6081393) q[1];
rz(0.98104279) q[3];
sx q[3];
rz(-1.0372122) q[3];
sx q[3];
rz(-0.55280815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2934908) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(2.9254986) q[2];
rz(-2.6381093) q[3];
sx q[3];
rz(-1.2516021) q[3];
sx q[3];
rz(-3.0751198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1450495) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(0.37522069) q[0];
rz(-0.61620617) q[1];
sx q[1];
rz(-2.1245978) q[1];
sx q[1];
rz(-2.9527051) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8804042) q[0];
sx q[0];
rz(-1.1587901) q[0];
sx q[0];
rz(2.3212577) q[0];
rz(-pi) q[1];
rz(-2.1470931) q[2];
sx q[2];
rz(-0.81953555) q[2];
sx q[2];
rz(-0.34162765) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9465883) q[1];
sx q[1];
rz(-0.71041162) q[1];
sx q[1];
rz(-1.2571774) q[1];
rz(-pi) q[2];
rz(1.1718719) q[3];
sx q[3];
rz(-1.0780943) q[3];
sx q[3];
rz(-2.1239831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45794332) q[2];
sx q[2];
rz(-1.8234437) q[2];
sx q[2];
rz(1.9981492) q[2];
rz(2.5055366) q[3];
sx q[3];
rz(-3.0885185) q[3];
sx q[3];
rz(1.5422356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77883333) q[0];
sx q[0];
rz(-2.9017359) q[0];
sx q[0];
rz(-1.1676316) q[0];
rz(-3.0150343) q[1];
sx q[1];
rz(-1.626222) q[1];
sx q[1];
rz(2.5680241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.180998) q[0];
sx q[0];
rz(-0.74686909) q[0];
sx q[0];
rz(-1.7043714) q[0];
rz(-pi) q[1];
rz(0.36678183) q[2];
sx q[2];
rz(-1.8213118) q[2];
sx q[2];
rz(2.9652852) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.998497) q[1];
sx q[1];
rz(-1.5228038) q[1];
sx q[1];
rz(2.001862) q[1];
rz(-pi) q[2];
rz(-2.4819991) q[3];
sx q[3];
rz(-1.9232188) q[3];
sx q[3];
rz(3.0291758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6557189) q[2];
sx q[2];
rz(-1.8334917) q[2];
sx q[2];
rz(-1.6386848) q[2];
rz(1.8466628) q[3];
sx q[3];
rz(-2.86627) q[3];
sx q[3];
rz(0.82726971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5591705) q[0];
sx q[0];
rz(-2.9952413) q[0];
sx q[0];
rz(2.225112) q[0];
rz(0.16042635) q[1];
sx q[1];
rz(-1.285099) q[1];
sx q[1];
rz(2.4442576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.682919) q[0];
sx q[0];
rz(-1.3240755) q[0];
sx q[0];
rz(1.9015067) q[0];
rz(0.38390145) q[2];
sx q[2];
rz(-2.4426201) q[2];
sx q[2];
rz(-2.0052103) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8749344) q[1];
sx q[1];
rz(-1.8968079) q[1];
sx q[1];
rz(-2.3951016) q[1];
x q[2];
rz(1.55682) q[3];
sx q[3];
rz(-1.4406573) q[3];
sx q[3];
rz(1.997245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15630284) q[2];
sx q[2];
rz(-2.6142945) q[2];
sx q[2];
rz(-1.1227597) q[2];
rz(1.5924234) q[3];
sx q[3];
rz(-1.4607818) q[3];
sx q[3];
rz(2.7482225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3093981) q[0];
sx q[0];
rz(-0.10646146) q[0];
sx q[0];
rz(2.0124117) q[0];
rz(-2.26217) q[1];
sx q[1];
rz(-1.8303454) q[1];
sx q[1];
rz(-2.5130491) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.382269) q[0];
sx q[0];
rz(-0.17062561) q[0];
sx q[0];
rz(2.3294446) q[0];
rz(-pi) q[1];
rz(-1.0209924) q[2];
sx q[2];
rz(-0.98787427) q[2];
sx q[2];
rz(1.1958808) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6595093) q[1];
sx q[1];
rz(-1.0190367) q[1];
sx q[1];
rz(-0.047288069) q[1];
rz(-pi) q[2];
rz(0.83129779) q[3];
sx q[3];
rz(-1.4834838) q[3];
sx q[3];
rz(1.0333453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92155543) q[2];
sx q[2];
rz(-2.3735789) q[2];
sx q[2];
rz(1.3735695) q[2];
rz(0.31275648) q[3];
sx q[3];
rz(-2.0650568) q[3];
sx q[3];
rz(-3.055174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3304928) q[0];
sx q[0];
rz(-2.436315) q[0];
sx q[0];
rz(0.16264859) q[0];
rz(-1.9074408) q[1];
sx q[1];
rz(-1.1470497) q[1];
sx q[1];
rz(2.2208234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2049157) q[0];
sx q[0];
rz(-1.7584118) q[0];
sx q[0];
rz(0.102392) q[0];
rz(-pi) q[1];
rz(-0.92451325) q[2];
sx q[2];
rz(-1.1072031) q[2];
sx q[2];
rz(-1.6789503) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.28345767) q[1];
sx q[1];
rz(-2.3871536) q[1];
sx q[1];
rz(-0.094875022) q[1];
rz(-pi) q[2];
rz(0.1699888) q[3];
sx q[3];
rz(-2.3521479) q[3];
sx q[3];
rz(-1.1657749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2938701) q[2];
sx q[2];
rz(-2.4203478) q[2];
sx q[2];
rz(3.0976307) q[2];
rz(-1.4219159) q[3];
sx q[3];
rz(-2.7094262) q[3];
sx q[3];
rz(-3.1066084) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111572) q[0];
sx q[0];
rz(-0.80428094) q[0];
sx q[0];
rz(-0.093611896) q[0];
rz(1.0253819) q[1];
sx q[1];
rz(-1.5954115) q[1];
sx q[1];
rz(-0.78790087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87374114) q[0];
sx q[0];
rz(-2.2722167) q[0];
sx q[0];
rz(2.3021477) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4008994) q[2];
sx q[2];
rz(-2.133103) q[2];
sx q[2];
rz(0.0058836939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9331197) q[1];
sx q[1];
rz(-1.6115687) q[1];
sx q[1];
rz(-0.59894125) q[1];
x q[2];
rz(-1.4282088) q[3];
sx q[3];
rz(-1.4514006) q[3];
sx q[3];
rz(0.52055955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0660144) q[2];
sx q[2];
rz(-0.54328537) q[2];
sx q[2];
rz(-0.96599609) q[2];
rz(-2.994359) q[3];
sx q[3];
rz(-1.4598264) q[3];
sx q[3];
rz(-2.537263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052893355) q[0];
sx q[0];
rz(-1.3857144) q[0];
sx q[0];
rz(1.2526441) q[0];
rz(-3.0839651) q[1];
sx q[1];
rz(-1.7689972) q[1];
sx q[1];
rz(2.7431814) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3301901) q[0];
sx q[0];
rz(-2.5017362) q[0];
sx q[0];
rz(2.8485203) q[0];
rz(-pi) q[1];
rz(-0.49968991) q[2];
sx q[2];
rz(-1.7202639) q[2];
sx q[2];
rz(2.2563344) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1507453) q[1];
sx q[1];
rz(-1.7707033) q[1];
sx q[1];
rz(0.0041109494) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6774954) q[3];
sx q[3];
rz(-2.1114285) q[3];
sx q[3];
rz(-1.65834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.56732279) q[2];
sx q[2];
rz(-1.8696573) q[2];
sx q[2];
rz(2.2197913) q[2];
rz(-0.70720339) q[3];
sx q[3];
rz(-2.0965529) q[3];
sx q[3];
rz(2.8309477) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1929753) q[0];
sx q[0];
rz(-3.0204168) q[0];
sx q[0];
rz(0.90712547) q[0];
rz(2.6616197) q[1];
sx q[1];
rz(-2.0860806) q[1];
sx q[1];
rz(2.5352246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1760178) q[0];
sx q[0];
rz(-1.8794354) q[0];
sx q[0];
rz(-1.465331) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53177436) q[2];
sx q[2];
rz(-2.0645879) q[2];
sx q[2];
rz(2.6419248) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.231016) q[1];
sx q[1];
rz(-1.2779932) q[1];
sx q[1];
rz(2.3412555) q[1];
rz(-pi) q[2];
rz(-2.4828047) q[3];
sx q[3];
rz(-1.4843936) q[3];
sx q[3];
rz(-2.1571546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85764641) q[2];
sx q[2];
rz(-1.5479156) q[2];
sx q[2];
rz(0.66961163) q[2];
rz(-2.5745463) q[3];
sx q[3];
rz(-1.0457057) q[3];
sx q[3];
rz(-1.4001575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46212101) q[0];
sx q[0];
rz(-1.7993878) q[0];
sx q[0];
rz(-2.2456428) q[0];
rz(-0.040291928) q[1];
sx q[1];
rz(-0.94172421) q[1];
sx q[1];
rz(1.5641854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9624845) q[0];
sx q[0];
rz(-0.024027457) q[0];
sx q[0];
rz(-1.3732617) q[0];
rz(2.4559458) q[2];
sx q[2];
rz(-1.9138991) q[2];
sx q[2];
rz(-2.1941663) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9312694) q[1];
sx q[1];
rz(-1.0627373) q[1];
sx q[1];
rz(-2.6711051) q[1];
rz(-pi) q[2];
rz(2.0730998) q[3];
sx q[3];
rz(-0.66876679) q[3];
sx q[3];
rz(-0.34688244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8335874) q[2];
sx q[2];
rz(-2.8218125) q[2];
sx q[2];
rz(1.8316815) q[2];
rz(-1.9561907) q[3];
sx q[3];
rz(-2.2535321) q[3];
sx q[3];
rz(-1.6185224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0014521) q[0];
sx q[0];
rz(-1.8403213) q[0];
sx q[0];
rz(2.1885827) q[0];
rz(-3.0630655) q[1];
sx q[1];
rz(-1.086906) q[1];
sx q[1];
rz(-0.79798098) q[1];
rz(-2.2453515) q[2];
sx q[2];
rz(-1.7755388) q[2];
sx q[2];
rz(-0.15646738) q[2];
rz(-1.4137951) q[3];
sx q[3];
rz(-1.3678322) q[3];
sx q[3];
rz(2.598138) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
