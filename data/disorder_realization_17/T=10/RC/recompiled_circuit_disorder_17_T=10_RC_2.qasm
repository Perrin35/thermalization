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
rz(0.97283483) q[1];
sx q[1];
rz(-1.4714779) q[1];
sx q[1];
rz(0.29247984) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4300633) q[0];
sx q[0];
rz(-1.9646514) q[0];
sx q[0];
rz(-3.0815365) q[0];
x q[1];
rz(2.3621759) q[2];
sx q[2];
rz(-1.6072818) q[2];
sx q[2];
rz(-1.905029) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3136183) q[1];
sx q[1];
rz(-2.5458114) q[1];
sx q[1];
rz(-1.6525364) q[1];
x q[2];
rz(-1.5745893) q[3];
sx q[3];
rz(-0.47382254) q[3];
sx q[3];
rz(2.4592318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73074377) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(-0.4494108) q[2];
rz(2.4959026) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9280424) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(-0.74203062) q[0];
rz(1.6702601) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(-1.3630294) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3516386) q[0];
sx q[0];
rz(-2.3839256) q[0];
sx q[0];
rz(-0.58654465) q[0];
rz(-pi) q[1];
rz(-1.7861373) q[2];
sx q[2];
rz(-1.1110348) q[2];
sx q[2];
rz(-2.721867) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.084761707) q[1];
sx q[1];
rz(-2.8014604) q[1];
sx q[1];
rz(-3.0840193) q[1];
x q[2];
rz(-1.3046673) q[3];
sx q[3];
rz(-1.7566655) q[3];
sx q[3];
rz(2.4974125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30119511) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(-1.2724686) q[2];
rz(2.8043591) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(-3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0740046) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(-2.8787676) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.6825312) q[1];
sx q[1];
rz(-1.9062818) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4609769) q[0];
sx q[0];
rz(-1.5810228) q[0];
sx q[0];
rz(-3.1264722) q[0];
rz(-pi) q[1];
rz(-1.4624865) q[2];
sx q[2];
rz(-2.8998067) q[2];
sx q[2];
rz(-2.0568648) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.77093) q[1];
sx q[1];
rz(-1.2491033) q[1];
sx q[1];
rz(1.4598203) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34706195) q[3];
sx q[3];
rz(-2.2357781) q[3];
sx q[3];
rz(2.5393328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0170903) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(-0.18033218) q[2];
rz(-0.63594937) q[3];
sx q[3];
rz(-1.9790117) q[3];
sx q[3];
rz(0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76698774) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(2.4225127) q[0];
rz(0.51302296) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(-2.1069353) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35095222) q[0];
sx q[0];
rz(-1.5119023) q[0];
sx q[0];
rz(1.5509997) q[0];
rz(3.0940975) q[2];
sx q[2];
rz(-1.4357899) q[2];
sx q[2];
rz(-2.6035655) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3116341) q[1];
sx q[1];
rz(-1.8533851) q[1];
sx q[1];
rz(-0.71210536) q[1];
rz(-2.8095874) q[3];
sx q[3];
rz(-0.49626353) q[3];
sx q[3];
rz(-2.9018324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93306142) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(-0.2362403) q[2];
rz(-1.158372) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(-0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1119969) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(-2.239256) q[0];
rz(0.92102712) q[1];
sx q[1];
rz(-0.61704707) q[1];
sx q[1];
rz(-0.62228084) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.92684) q[0];
sx q[0];
rz(-1.2878969) q[0];
sx q[0];
rz(2.9515285) q[0];
x q[1];
rz(2.4933382) q[2];
sx q[2];
rz(-1.4972685) q[2];
sx q[2];
rz(1.3157805) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29187782) q[1];
sx q[1];
rz(-1.4018702) q[1];
sx q[1];
rz(0.090976322) q[1];
rz(-2.6617674) q[3];
sx q[3];
rz(-1.92355) q[3];
sx q[3];
rz(1.477581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2174012) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(-3.1203111) q[2];
rz(-1.4422669) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(-2.2309979) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845881) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(-3.0969627) q[0];
rz(2.1394829) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8196626) q[0];
sx q[0];
rz(-0.18580431) q[0];
sx q[0];
rz(1.4968605) q[0];
x q[1];
rz(-1.0677412) q[2];
sx q[2];
rz(-1.6776553) q[2];
sx q[2];
rz(1.5397132) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1019563) q[1];
sx q[1];
rz(-0.64412457) q[1];
sx q[1];
rz(-2.1307751) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4061635) q[3];
sx q[3];
rz(-2.2693686) q[3];
sx q[3];
rz(-0.64100953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78952152) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(2.1748523) q[2];
rz(3.1294075) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
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
rz(-1.8693285) q[1];
sx q[1];
rz(-0.73928839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4261739) q[0];
sx q[0];
rz(-2.9746375) q[0];
sx q[0];
rz(0.963361) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0871068) q[2];
sx q[2];
rz(-1.0488044) q[2];
sx q[2];
rz(2.1453478) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48078295) q[1];
sx q[1];
rz(-2.0936831) q[1];
sx q[1];
rz(-1.7899662) q[1];
x q[2];
rz(2.7474653) q[3];
sx q[3];
rz(-2.4836802) q[3];
sx q[3];
rz(-3.1042299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5560975) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(2.5778256) q[2];
rz(2.901315) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(-1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69650841) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(0.58404303) q[0];
rz(2.7208327) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(-0.27841321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4415092) q[0];
sx q[0];
rz(-1.004305) q[0];
sx q[0];
rz(2.9928815) q[0];
rz(-pi) q[1];
rz(-2.7114696) q[2];
sx q[2];
rz(-0.42818907) q[2];
sx q[2];
rz(2.6870514) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0502187) q[1];
sx q[1];
rz(-0.41793567) q[1];
sx q[1];
rz(1.0068847) q[1];
rz(-pi) q[2];
rz(-1.9592459) q[3];
sx q[3];
rz(-0.58245814) q[3];
sx q[3];
rz(2.9782481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.86928308) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(2.753624) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85912722) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(-0.7318837) q[0];
rz(2.9751119) q[1];
sx q[1];
rz(-2.6961168) q[1];
sx q[1];
rz(0.073721185) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8548944) q[0];
sx q[0];
rz(-0.44730967) q[0];
sx q[0];
rz(-1.8064503) q[0];
x q[1];
rz(3.0606255) q[2];
sx q[2];
rz(-1.690174) q[2];
sx q[2];
rz(2.2031684) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.24450609) q[1];
sx q[1];
rz(-0.58468854) q[1];
sx q[1];
rz(-0.067824407) q[1];
x q[2];
rz(-1.4599667) q[3];
sx q[3];
rz(-1.0004527) q[3];
sx q[3];
rz(-0.092872083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59447294) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(-2.7129042) q[2];
rz(1.3171014) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(1.2111506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4093032) q[0];
sx q[0];
rz(-1.4855054) q[0];
sx q[0];
rz(1.0428585) q[0];
rz(1.6304784) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(1.7932549) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74603426) q[0];
sx q[0];
rz(-1.4575882) q[0];
sx q[0];
rz(3.1157225) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4820319) q[2];
sx q[2];
rz(-2.4030493) q[2];
sx q[2];
rz(-0.86631394) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7235782) q[1];
sx q[1];
rz(-1.8943159) q[1];
sx q[1];
rz(3.0669206) q[1];
rz(-pi) q[2];
rz(-0.24095778) q[3];
sx q[3];
rz(-1.3853419) q[3];
sx q[3];
rz(-1.0518187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(0.094816118) q[2];
rz(-1.173165) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(1.7839292) q[2];
sx q[2];
rz(-1.4896981) q[2];
sx q[2];
rz(-1.8427195) q[2];
rz(-2.1223162) q[3];
sx q[3];
rz(-0.84001361) q[3];
sx q[3];
rz(2.7915814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
