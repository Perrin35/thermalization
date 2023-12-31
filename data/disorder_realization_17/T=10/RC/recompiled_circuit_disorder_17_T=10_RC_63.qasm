OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9085812) q[0];
sx q[0];
rz(4.3281875) q[0];
sx q[0];
rz(8.2789658) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(2.8491128) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16380331) q[0];
sx q[0];
rz(-1.6262494) q[0];
sx q[0];
rz(-1.1763014) q[0];
rz(-3.0897065) q[2];
sx q[2];
rz(-0.78009006) q[2];
sx q[2];
rz(2.8442596) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9521192) q[1];
sx q[1];
rz(-1.6166302) q[1];
sx q[1];
rz(0.97656753) q[1];
rz(1.0969767) q[3];
sx q[3];
rz(-1.572527) q[3];
sx q[3];
rz(2.2497821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73074377) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(-0.4494108) q[2];
rz(0.64569008) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(-1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21355024) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(-2.399562) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(1.7785633) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3516386) q[0];
sx q[0];
rz(-0.75766701) q[0];
sx q[0];
rz(-0.58654465) q[0];
rz(-pi) q[1];
rz(2.6724733) q[2];
sx q[2];
rz(-1.763478) q[2];
sx q[2];
rz(2.0872781) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9957673) q[1];
sx q[1];
rz(-1.9103423) q[1];
sx q[1];
rz(1.550436) q[1];
rz(2.1915073) q[3];
sx q[3];
rz(-2.8182497) q[3];
sx q[3];
rz(-2.8107373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.30119511) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(-1.2724686) q[2];
rz(-0.33723351) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(-3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067588016) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(0.2628251) q[0];
rz(-2.7858531) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(-1.9062818) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0316186) q[0];
sx q[0];
rz(-1.5556766) q[0];
sx q[0];
rz(-1.5605687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.026651816) q[2];
sx q[2];
rz(-1.3304552) q[2];
sx q[2];
rz(-1.9453366) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1097818) q[1];
sx q[1];
rz(-2.8019252) q[1];
sx q[1];
rz(-2.8207645) q[1];
rz(1.9800817) q[3];
sx q[3];
rz(-2.403879) q[3];
sx q[3];
rz(-1.1324594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1245023) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(-2.9612605) q[2];
rz(2.5056433) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76698774) q[0];
sx q[0];
rz(-2.7174482) q[0];
sx q[0];
rz(-0.71907991) q[0];
rz(2.6285697) q[1];
sx q[1];
rz(-1.9388371) q[1];
sx q[1];
rz(1.0346574) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67544115) q[0];
sx q[0];
rz(-3.0794641) q[0];
sx q[0];
rz(-0.32390578) q[0];
rz(1.4356394) q[2];
sx q[2];
rz(-1.617859) q[2];
sx q[2];
rz(2.1152209) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0532916) q[1];
sx q[1];
rz(-0.75690714) q[1];
sx q[1];
rz(-2.7234368) q[1];
x q[2];
rz(-0.33200522) q[3];
sx q[3];
rz(-0.49626353) q[3];
sx q[3];
rz(2.9018324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2085312) q[2];
sx q[2];
rz(-0.26991093) q[2];
sx q[2];
rz(0.2362403) q[2];
rz(1.158372) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(-2.2896144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029595705) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(-2.239256) q[0];
rz(-2.2205655) q[1];
sx q[1];
rz(-0.61704707) q[1];
sx q[1];
rz(-0.62228084) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3233933) q[0];
sx q[0];
rz(-2.8022031) q[0];
sx q[0];
rz(-2.1470977) q[0];
rz(-pi) q[1];
rz(1.4786517) q[2];
sx q[2];
rz(-2.2170057) q[2];
sx q[2];
rz(-0.19942936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2942549) q[1];
sx q[1];
rz(-1.4811185) q[1];
sx q[1];
rz(1.7404106) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47982521) q[3];
sx q[3];
rz(-1.2180426) q[3];
sx q[3];
rz(1.6640116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92419147) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(1.4422669) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(-0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6570046) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(-0.044629991) q[0];
rz(-2.1394829) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1761985) q[0];
sx q[0];
rz(-1.584443) q[0];
sx q[0];
rz(1.3854881) q[0];
rz(-pi) q[1];
rz(1.7897286) q[2];
sx q[2];
rz(-0.5133252) q[2];
sx q[2];
rz(-2.9190612) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9958482) q[1];
sx q[1];
rz(-1.8954344) q[1];
sx q[1];
rz(-1.0041298) q[1];
x q[2];
rz(-1.4061635) q[3];
sx q[3];
rz(-0.87222404) q[3];
sx q[3];
rz(-0.64100953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78952152) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(-0.96674031) q[2];
rz(3.1294075) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(3.0325082) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0063342) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(2.1475041) q[0];
rz(-0.055123568) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(-2.4023043) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2544884) q[0];
sx q[0];
rz(-1.4758037) q[0];
sx q[0];
rz(-1.7083005) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57619608) q[2];
sx q[2];
rz(-1.9857166) q[2];
sx q[2];
rz(-0.31838271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48078295) q[1];
sx q[1];
rz(-1.0479095) q[1];
sx q[1];
rz(1.7899662) q[1];
rz(-pi) q[2];
rz(-1.8592632) q[3];
sx q[3];
rz(-0.97086421) q[3];
sx q[3];
rz(2.695042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58549515) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(-2.5778256) q[2];
rz(-0.24027763) q[3];
sx q[3];
rz(-0.53301817) q[3];
sx q[3];
rz(1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4450842) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(-2.5575496) q[0];
rz(0.42075992) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(0.27841321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70008343) q[0];
sx q[0];
rz(-2.1372876) q[0];
sx q[0];
rz(-0.14871116) q[0];
rz(-2.7483447) q[2];
sx q[2];
rz(-1.7448145) q[2];
sx q[2];
rz(2.4207123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.003493017) q[1];
sx q[1];
rz(-1.3521191) q[1];
sx q[1];
rz(1.2117282) q[1];
x q[2];
rz(-0.24448963) q[3];
sx q[3];
rz(-2.1049307) q[3];
sx q[3];
rz(0.29230803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2723096) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(-0.38796866) q[2];
rz(1.7913473) q[3];
sx q[3];
rz(-2.6781121) q[3];
sx q[3];
rz(-0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85912722) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(0.7318837) q[0];
rz(-2.9751119) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(0.073721185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4972992) q[0];
sx q[0];
rz(-1.4696346) q[0];
sx q[0];
rz(2.0072719) q[0];
rz(-1.45103) q[2];
sx q[2];
rz(-1.6511859) q[2];
sx q[2];
rz(-2.5188841) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8718865) q[1];
sx q[1];
rz(-1.5333813) q[1];
sx q[1];
rz(-2.5579631) q[1];
rz(2.5684486) q[3];
sx q[3];
rz(-1.4775652) q[3];
sx q[3];
rz(1.4179109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59447294) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(-2.7129042) q[2];
rz(-1.8244913) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(-1.930442) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73228943) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(2.0987341) q[0];
rz(1.6304784) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(-1.7932549) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3955584) q[0];
sx q[0];
rz(-1.4575882) q[0];
sx q[0];
rz(3.1157225) q[0];
x q[1];
rz(0.080532455) q[2];
sx q[2];
rz(-2.3057612) q[2];
sx q[2];
rz(2.1554975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7235782) q[1];
sx q[1];
rz(-1.2472767) q[1];
sx q[1];
rz(-0.07467204) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7616368) q[3];
sx q[3];
rz(-1.8075426) q[3];
sx q[3];
rz(2.5773347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(0.094816118) q[2];
rz(1.173165) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(-2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.512758) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(1.5402773) q[1];
sx q[1];
rz(-1.5581144) q[1];
sx q[1];
rz(1.5099572) q[1];
rz(1.9376471) q[2];
sx q[2];
rz(-0.22782142) q[2];
sx q[2];
rz(0.08624764) q[2];
rz(-2.6125828) q[3];
sx q[3];
rz(-0.88376868) q[3];
sx q[3];
rz(0.39467011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
