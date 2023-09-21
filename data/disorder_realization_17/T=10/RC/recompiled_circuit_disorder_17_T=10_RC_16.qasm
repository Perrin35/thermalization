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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2746409) q[0];
sx q[0];
rz(-0.39817087) q[0];
sx q[0];
rz(1.7142332) q[0];
rz(-2.3621759) q[2];
sx q[2];
rz(-1.6072818) q[2];
sx q[2];
rz(1.905029) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7293207) q[1];
sx q[1];
rz(-2.1643157) q[1];
sx q[1];
rz(0.055298474) q[1];
x q[2];
rz(2.0446159) q[3];
sx q[3];
rz(-1.572527) q[3];
sx q[3];
rz(-2.2497821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73074377) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(0.4494108) q[2];
rz(2.4959026) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(1.8507563) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21355024) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(-0.74203062) q[0];
rz(-1.4713326) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(-1.3630294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5308967) q[0];
sx q[0];
rz(-0.96141059) q[0];
sx q[0];
rz(-2.0531528) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6724733) q[2];
sx q[2];
rz(-1.763478) q[2];
sx q[2];
rz(1.0543146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.084761707) q[1];
sx q[1];
rz(-0.34013224) q[1];
sx q[1];
rz(3.0840193) q[1];
rz(-pi) q[2];
rz(0.95008534) q[3];
sx q[3];
rz(-2.8182497) q[3];
sx q[3];
rz(2.8107373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8403975) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(-1.8691241) q[2];
rz(0.33723351) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(-0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067588016) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(0.2628251) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(1.9062818) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4848029) q[0];
sx q[0];
rz(-3.1233388) q[0];
sx q[0];
rz(2.546893) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8112196) q[2];
sx q[2];
rz(-1.5966819) q[2];
sx q[2];
rz(0.38088574) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23535145) q[1];
sx q[1];
rz(-1.4655349) q[1];
sx q[1];
rz(-0.32354849) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9800817) q[3];
sx q[3];
rz(-0.73771362) q[3];
sx q[3];
rz(-2.0091332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0170903) q[2];
sx q[2];
rz(-1.7586781) q[2];
sx q[2];
rz(-2.9612605) q[2];
rz(-2.5056433) q[3];
sx q[3];
rz(-1.9790117) q[3];
sx q[3];
rz(2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
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
rz(-pi) q[2];
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
rz(1.0346574) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9229139) q[0];
sx q[0];
rz(-1.551034) q[0];
sx q[0];
rz(0.058905525) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7059533) q[2];
sx q[2];
rz(-1.617859) q[2];
sx q[2];
rz(-2.1152209) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0883011) q[1];
sx q[1];
rz(-2.3846855) q[1];
sx q[1];
rz(-0.41815586) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3961117) q[3];
sx q[3];
rz(-1.1038728) q[3];
sx q[3];
rz(-0.13388453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(2.5193118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392004) q[0];
sx q[0];
rz(-1.3883739) q[0];
sx q[0];
rz(-1.8586041) q[0];
x q[1];
rz(0.64825443) q[2];
sx q[2];
rz(-1.6443242) q[2];
sx q[2];
rz(1.3157805) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9363074) q[1];
sx q[1];
rz(-0.19166066) q[1];
sx q[1];
rz(1.081341) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.673224) q[3];
sx q[3];
rz(-2.554318) q[3];
sx q[3];
rz(-2.4622963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2174012) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(-3.1203111) q[2];
rz(-1.6993258) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(2.2309979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6570046) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(-0.044629991) q[0];
rz(1.0021098) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-3.1071641) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653942) q[0];
sx q[0];
rz(-1.584443) q[0];
sx q[0];
rz(-1.3854881) q[0];
rz(3.0197633) q[2];
sx q[2];
rz(-1.0708772) q[2];
sx q[2];
rz(3.1140529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.039636314) q[1];
sx q[1];
rz(-0.64412457) q[1];
sx q[1];
rz(2.1307751) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19272007) q[3];
sx q[3];
rz(-0.7145213) q[3];
sx q[3];
rz(-0.89380985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3520711) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(0.96674031) q[2];
rz(3.1294075) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352585) q[0];
sx q[0];
rz(-0.2671347) q[0];
sx q[0];
rz(-2.1475041) q[0];
rz(-0.055123568) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(-0.73928839) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3294322) q[0];
sx q[0];
rz(-1.433916) q[0];
sx q[0];
rz(0.095892266) q[0];
rz(-pi) q[1];
rz(0.67989345) q[2];
sx q[2];
rz(-2.4455564) q[2];
sx q[2];
rz(-1.807715) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.90034396) q[1];
sx q[1];
rz(-2.5785891) q[1];
sx q[1];
rz(-0.36069936) q[1];
rz(0.61974157) q[3];
sx q[3];
rz(-1.3337787) q[3];
sx q[3];
rz(-1.851351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58549515) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(2.901315) q[3];
sx q[3];
rz(-0.53301817) q[3];
sx q[3];
rz(1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.4450842) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(2.5575496) q[0];
rz(0.42075992) q[1];
sx q[1];
rz(-1.0412443) q[1];
sx q[1];
rz(-0.27841321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42785545) q[0];
sx q[0];
rz(-0.58361485) q[0];
sx q[0];
rz(-1.7996656) q[0];
rz(-pi) q[1];
rz(2.7483447) q[2];
sx q[2];
rz(-1.3967782) q[2];
sx q[2];
rz(2.4207123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0502187) q[1];
sx q[1];
rz(-2.723657) q[1];
sx q[1];
rz(1.0068847) q[1];
rz(0.24448963) q[3];
sx q[3];
rz(-2.1049307) q[3];
sx q[3];
rz(2.8492846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85912722) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(-0.16648079) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(3.0678715) q[1];
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
rz(-1.1343207) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.164052) q[2];
sx q[2];
rz(-2.9974555) q[2];
sx q[2];
rz(-1.6050715) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24450609) q[1];
sx q[1];
rz(-0.58468854) q[1];
sx q[1];
rz(3.0737682) q[1];
x q[2];
rz(-2.9708462) q[3];
sx q[3];
rz(-2.5617544) q[3];
sx q[3];
rz(-0.29614007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59447294) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(0.42868844) q[2];
rz(-1.8244913) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4093032) q[0];
sx q[0];
rz(-1.4855054) q[0];
sx q[0];
rz(-1.0428585) q[0];
rz(1.5111142) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(-1.7932549) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139075) q[0];
sx q[0];
rz(-1.5450918) q[0];
sx q[0];
rz(1.684042) q[0];
rz(-2.3073763) q[2];
sx q[2];
rz(-1.6305106) q[2];
sx q[2];
rz(-2.5028253) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1765602) q[1];
sx q[1];
rz(-1.6415879) q[1];
sx q[1];
rz(-1.2464347) q[1];
rz(-pi) q[2];
rz(-2.4753307) q[3];
sx q[3];
rz(-2.8386335) q[3];
sx q[3];
rz(-0.12476866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80031359) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(-3.0467765) q[2];
rz(1.9684277) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-0.512758) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(-1.6013153) q[1];
sx q[1];
rz(-1.5581144) q[1];
sx q[1];
rz(1.5099572) q[1];
rz(-1.9376471) q[2];
sx q[2];
rz(-2.9137712) q[2];
sx q[2];
rz(-3.055345) q[2];
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
