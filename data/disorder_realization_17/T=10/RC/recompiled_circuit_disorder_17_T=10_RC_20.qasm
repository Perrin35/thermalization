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
rz(-1.9549978) q[0];
sx q[0];
rz(-1.1458122) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(-0.29247984) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2746409) q[0];
sx q[0];
rz(-2.7434218) q[0];
sx q[0];
rz(1.7142332) q[0];
x q[1];
rz(-0.77941676) q[2];
sx q[2];
rz(-1.5343108) q[2];
sx q[2];
rz(-1.2365637) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18947345) q[1];
sx q[1];
rz(-1.5249624) q[1];
sx q[1];
rz(-0.97656753) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0446159) q[3];
sx q[3];
rz(-1.572527) q[3];
sx q[3];
rz(-0.8918106) q[3];
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
rz(2.6921819) q[2];
rz(0.64569008) q[3];
sx q[3];
rz(-0.68104762) q[3];
sx q[3];
rz(1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21355024) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(2.399562) q[0];
rz(-1.4713326) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(1.7785633) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3516386) q[0];
sx q[0];
rz(-2.3839256) q[0];
sx q[0];
rz(2.555048) q[0];
rz(-pi) q[1];
rz(-0.46911932) q[2];
sx q[2];
rz(-1.763478) q[2];
sx q[2];
rz(2.0872781) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.431753) q[1];
sx q[1];
rz(-1.589994) q[1];
sx q[1];
rz(0.33961105) q[1];
x q[2];
rz(2.9491049) q[3];
sx q[3];
rz(-1.8322332) q[3];
sx q[3];
rz(-2.1646433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30119511) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(-1.8691241) q[2];
rz(2.8043591) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(-0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0740046) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(-2.8787676) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(1.9062818) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4609769) q[0];
sx q[0];
rz(-1.5810228) q[0];
sx q[0];
rz(-0.015120487) q[0];
rz(-pi) q[1];
rz(-3.1149408) q[2];
sx q[2];
rz(-1.8111374) q[2];
sx q[2];
rz(-1.196256) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0318109) q[1];
sx q[1];
rz(-2.8019252) q[1];
sx q[1];
rz(2.8207645) q[1];
x q[2];
rz(2.7945307) q[3];
sx q[3];
rz(-2.2357781) q[3];
sx q[3];
rz(2.5393328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0170903) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(-2.9612605) q[2];
rz(-0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3746049) q[0];
sx q[0];
rz(-2.7174482) q[0];
sx q[0];
rz(2.4225127) q[0];
rz(-2.6285697) q[1];
sx q[1];
rz(-1.9388371) q[1];
sx q[1];
rz(2.1069353) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9229139) q[0];
sx q[0];
rz(-1.551034) q[0];
sx q[0];
rz(-0.058905525) q[0];
x q[1];
rz(1.2345418) q[2];
sx q[2];
rz(-0.14306919) q[2];
sx q[2];
rz(-0.87749315) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0883011) q[1];
sx q[1];
rz(-2.3846855) q[1];
sx q[1];
rz(-2.7234368) q[1];
x q[2];
rz(0.47311584) q[3];
sx q[3];
rz(-1.72662) q[3];
sx q[3];
rz(-1.5161878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2085312) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(-2.9053524) q[2];
rz(-1.9832206) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(-2.2896144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029595705) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(0.90233666) q[0];
rz(-2.2205655) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-2.5193118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392004) q[0];
sx q[0];
rz(-1.7532187) q[0];
sx q[0];
rz(-1.2829885) q[0];
x q[1];
rz(2.4933382) q[2];
sx q[2];
rz(-1.4972685) q[2];
sx q[2];
rz(1.3157805) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8497148) q[1];
sx q[1];
rz(-1.4018702) q[1];
sx q[1];
rz(3.0506163) q[1];
x q[2];
rz(0.673224) q[3];
sx q[3];
rz(-2.554318) q[3];
sx q[3];
rz(-0.67929635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.92419147) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(-1.4422669) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6570046) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(3.0969627) q[0];
rz(1.0021098) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(0.034428509) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39715595) q[0];
sx q[0];
rz(-1.3855055) q[0];
sx q[0];
rz(-3.1277083) q[0];
rz(0.12182932) q[2];
sx q[2];
rz(-2.0707154) q[2];
sx q[2];
rz(3.1140529) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.039636314) q[1];
sx q[1];
rz(-2.4974681) q[1];
sx q[1];
rz(1.0108175) q[1];
rz(-0.19272007) q[3];
sx q[3];
rz(-0.7145213) q[3];
sx q[3];
rz(-0.89380985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78952152) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(-2.1748523) q[2];
rz(3.1294075) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(3.0325082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0063342) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(2.1475041) q[0];
rz(0.055123568) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(0.73928839) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3294322) q[0];
sx q[0];
rz(-1.433916) q[0];
sx q[0];
rz(0.095892266) q[0];
rz(2.0544858) q[2];
sx q[2];
rz(-1.0488044) q[2];
sx q[2];
rz(-0.99624485) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6608097) q[1];
sx q[1];
rz(-1.0479095) q[1];
sx q[1];
rz(1.7899662) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61974157) q[3];
sx q[3];
rz(-1.8078139) q[3];
sx q[3];
rz(-1.2902416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.58549515) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(0.56376702) q[2];
rz(-2.901315) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(-1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4450842) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(-0.58404303) q[0];
rz(-2.7208327) q[1];
sx q[1];
rz(-1.0412443) q[1];
sx q[1];
rz(-0.27841321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4415092) q[0];
sx q[0];
rz(-2.1372876) q[0];
sx q[0];
rz(-0.14871116) q[0];
rz(-pi) q[1];
rz(-2.7483447) q[2];
sx q[2];
rz(-1.7448145) q[2];
sx q[2];
rz(2.4207123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1380996) q[1];
sx q[1];
rz(-1.3521191) q[1];
sx q[1];
rz(-1.2117282) q[1];
rz(-2.897103) q[3];
sx q[3];
rz(-2.1049307) q[3];
sx q[3];
rz(-0.29230803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.86928308) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(-2.753624) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-2.6781121) q[3];
sx q[3];
rz(0.2750245) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(-0.7318837) q[0];
rz(-2.9751119) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(-3.0678715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4972992) q[0];
sx q[0];
rz(-1.6719581) q[0];
sx q[0];
rz(2.0072719) q[0];
x q[1];
rz(-0.97754064) q[2];
sx q[2];
rz(-2.9974555) q[2];
sx q[2];
rz(1.6050715) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8970866) q[1];
sx q[1];
rz(-2.5569041) q[1];
sx q[1];
rz(-3.0737682) q[1];
rz(2.9708462) q[3];
sx q[3];
rz(-0.57983825) q[3];
sx q[3];
rz(-0.29614007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5471197) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(-2.7129042) q[2];
rz(-1.8244913) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.73228943) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(2.0987341) q[0];
rz(1.5111142) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(-1.7932549) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139075) q[0];
sx q[0];
rz(-1.5450918) q[0];
sx q[0];
rz(-1.4575507) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6595608) q[2];
sx q[2];
rz(-0.73854337) q[2];
sx q[2];
rz(0.86631394) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41801449) q[1];
sx q[1];
rz(-1.2472767) q[1];
sx q[1];
rz(0.07467204) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3799558) q[3];
sx q[3];
rz(-1.8075426) q[3];
sx q[3];
rz(2.5773347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3412791) q[2];
sx q[2];
rz(-2.6276402) q[2];
sx q[2];
rz(-3.0467765) q[2];
rz(1.173165) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(0.15869424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.512758) q[0];
sx q[0];
rz(-2.6678968) q[0];
sx q[0];
rz(-2.0976023) q[0];
rz(-1.5402773) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(1.3576635) q[2];
sx q[2];
rz(-1.6518946) q[2];
sx q[2];
rz(1.2988731) q[2];
rz(0.81090609) q[3];
sx q[3];
rz(-1.1699642) q[3];
sx q[3];
rz(-1.5311833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
