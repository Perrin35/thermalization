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
rz(-0.29247984) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7115294) q[0];
sx q[0];
rz(-1.1769413) q[0];
sx q[0];
rz(-3.0815365) q[0];
rz(-pi) q[1];
rz(0.77941676) q[2];
sx q[2];
rz(-1.5343108) q[2];
sx q[2];
rz(1.2365637) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3136183) q[1];
sx q[1];
rz(-2.5458114) q[1];
sx q[1];
rz(1.4890563) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0019449751) q[3];
sx q[3];
rz(-2.0446152) q[3];
sx q[3];
rz(0.67809826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4108489) q[2];
sx q[2];
rz(-1.0935254) q[2];
sx q[2];
rz(-0.4494108) q[2];
rz(0.64569008) q[3];
sx q[3];
rz(-0.68104762) q[3];
sx q[3];
rz(-1.2908363) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21355024) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(0.74203062) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(1.7785633) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3516386) q[0];
sx q[0];
rz(-2.3839256) q[0];
sx q[0];
rz(-2.555048) q[0];
rz(-pi) q[1];
rz(-1.3554553) q[2];
sx q[2];
rz(-1.1110348) q[2];
sx q[2];
rz(2.721867) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14582536) q[1];
sx q[1];
rz(-1.2312504) q[1];
sx q[1];
rz(-1.5911566) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3046673) q[3];
sx q[3];
rz(-1.3849272) q[3];
sx q[3];
rz(-0.64418018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30119511) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(1.2724686) q[2];
rz(2.8043591) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067588016) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(-0.2628251) q[0];
rz(0.35573959) q[1];
sx q[1];
rz(-1.6825312) q[1];
sx q[1];
rz(1.9062818) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6567898) q[0];
sx q[0];
rz(-3.1233388) q[0];
sx q[0];
rz(0.59469964) q[0];
rz(-1.3303731) q[2];
sx q[2];
rz(-1.5449107) q[2];
sx q[2];
rz(-2.7607069) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.23535145) q[1];
sx q[1];
rz(-1.6760577) q[1];
sx q[1];
rz(-0.32354849) q[1];
rz(-pi) q[2];
rz(2.2658306) q[3];
sx q[3];
rz(-1.8417629) q[3];
sx q[3];
rz(0.7489487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1245023) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(-0.18033218) q[2];
rz(-0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(-0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.2027556) q[1];
sx q[1];
rz(-1.0346574) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67544115) q[0];
sx q[0];
rz(-0.062128566) q[0];
sx q[0];
rz(-0.32390578) q[0];
rz(3.0940975) q[2];
sx q[2];
rz(-1.7058027) q[2];
sx q[2];
rz(2.6035655) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0532916) q[1];
sx q[1];
rz(-0.75690714) q[1];
sx q[1];
rz(-2.7234368) q[1];
rz(-pi) q[2];
rz(1.7454809) q[3];
sx q[3];
rz(-1.1038728) q[3];
sx q[3];
rz(0.13388453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93306142) q[2];
sx q[2];
rz(-0.26991093) q[2];
sx q[2];
rz(-2.9053524) q[2];
rz(1.158372) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(2.2896144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029595705) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(2.239256) q[0];
rz(2.2205655) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-0.62228084) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.92684) q[0];
sx q[0];
rz(-1.8536957) q[0];
sx q[0];
rz(-0.19006417) q[0];
x q[1];
rz(3.0201966) q[2];
sx q[2];
rz(-0.65181323) q[2];
sx q[2];
rz(2.7898942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29187782) q[1];
sx q[1];
rz(-1.7397225) q[1];
sx q[1];
rz(3.0506163) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6617674) q[3];
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
rz(-2.2174012) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(-0.021281555) q[2];
rz(-1.6993258) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(2.2309979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6570046) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(0.044629991) q[0];
rz(-1.0021098) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-0.034428509) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39715595) q[0];
sx q[0];
rz(-1.7560871) q[0];
sx q[0];
rz(-0.013884355) q[0];
rz(3.0197633) q[2];
sx q[2];
rz(-1.0708772) q[2];
sx q[2];
rz(-0.027539754) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.039636314) q[1];
sx q[1];
rz(-0.64412457) q[1];
sx q[1];
rz(-2.1307751) q[1];
rz(-pi) q[2];
rz(1.4061635) q[3];
sx q[3];
rz(-0.87222404) q[3];
sx q[3];
rz(-2.5005831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78952152) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(2.1748523) q[2];
rz(-3.1294075) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(3.0325082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0063342) q[0];
sx q[0];
rz(-0.2671347) q[0];
sx q[0];
rz(2.1475041) q[0];
rz(-0.055123568) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(-0.73928839) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3294322) q[0];
sx q[0];
rz(-1.433916) q[0];
sx q[0];
rz(3.0457004) q[0];
x q[1];
rz(-2.5653966) q[2];
sx q[2];
rz(-1.9857166) q[2];
sx q[2];
rz(-2.8232099) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97923219) q[1];
sx q[1];
rz(-1.3812961) q[1];
sx q[1];
rz(2.6081671) q[1];
x q[2];
rz(1.8592632) q[3];
sx q[3];
rz(-2.1707284) q[3];
sx q[3];
rz(-0.44655061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5560975) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(-0.24027763) q[3];
sx q[3];
rz(-0.53301817) q[3];
sx q[3];
rz(1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4450842) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(0.58404303) q[0];
rz(-2.7208327) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(0.27841321) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4415092) q[0];
sx q[0];
rz(-2.1372876) q[0];
sx q[0];
rz(-0.14871116) q[0];
rz(-2.7483447) q[2];
sx q[2];
rz(-1.3967782) q[2];
sx q[2];
rz(-2.4207123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6555357) q[1];
sx q[1];
rz(-1.2206519) q[1];
sx q[1];
rz(-0.23305594) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0233381) q[3];
sx q[3];
rz(-1.3609144) q[3];
sx q[3];
rz(1.4048214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2723096) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(-0.38796866) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(-0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(-2.9751119) q[1];
sx q[1];
rz(-2.6961168) q[1];
sx q[1];
rz(-0.073721185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6442935) q[0];
sx q[0];
rz(-1.4696346) q[0];
sx q[0];
rz(2.0072719) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97754064) q[2];
sx q[2];
rz(-0.14413713) q[2];
sx q[2];
rz(1.6050715) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8718865) q[1];
sx q[1];
rz(-1.5333813) q[1];
sx q[1];
rz(-2.5579631) q[1];
rz(-pi) q[2];
rz(0.1707465) q[3];
sx q[3];
rz(-0.57983825) q[3];
sx q[3];
rz(0.29614007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59447294) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(0.42868844) q[2];
rz(-1.3171014) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(-1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4093032) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(1.0428585) q[0];
rz(-1.6304784) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(-1.3483378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3955584) q[0];
sx q[0];
rz(-1.4575882) q[0];
sx q[0];
rz(3.1157225) q[0];
rz(-0.83421631) q[2];
sx q[2];
rz(-1.6305106) q[2];
sx q[2];
rz(-0.63876736) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.18689449) q[1];
sx q[1];
rz(-2.8098626) q[1];
sx q[1];
rz(-1.3518672) q[1];
x q[2];
rz(0.66626196) q[3];
sx q[3];
rz(-2.8386335) q[3];
sx q[3];
rz(-0.12476866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3412791) q[2];
sx q[2];
rz(-2.6276402) q[2];
sx q[2];
rz(-0.094816118) q[2];
rz(1.9684277) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(-0.15869424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6288347) q[0];
sx q[0];
rz(-2.6678968) q[0];
sx q[0];
rz(-2.0976023) q[0];
rz(1.5402773) q[1];
sx q[1];
rz(-1.5581144) q[1];
sx q[1];
rz(1.5099572) q[1];
rz(-1.2039456) q[2];
sx q[2];
rz(-0.22782142) q[2];
sx q[2];
rz(0.08624764) q[2];
rz(2.6125828) q[3];
sx q[3];
rz(-2.257824) q[3];
sx q[3];
rz(-2.7469225) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
