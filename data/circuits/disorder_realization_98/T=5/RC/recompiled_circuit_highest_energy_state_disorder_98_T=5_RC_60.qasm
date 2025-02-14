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
rz(-0.67796081) q[0];
sx q[0];
rz(-0.27422658) q[0];
sx q[0];
rz(1.8507313) q[0];
rz(-0.21386799) q[1];
sx q[1];
rz(-2.8254421) q[1];
sx q[1];
rz(1.4996127) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363269) q[0];
sx q[0];
rz(-1.8417497) q[0];
sx q[0];
rz(-1.3265557) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8617282) q[2];
sx q[2];
rz(-2.3290344) q[2];
sx q[2];
rz(1.9951374) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8655244) q[1];
sx q[1];
rz(-1.8929947) q[1];
sx q[1];
rz(2.404271) q[1];
rz(-0.47874449) q[3];
sx q[3];
rz(-2.1447721) q[3];
sx q[3];
rz(0.51306242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9700254) q[2];
sx q[2];
rz(-2.1000523) q[2];
sx q[2];
rz(0.901326) q[2];
rz(-2.6775635) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(-3.071781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-2.0483911) q[0];
sx q[0];
rz(-3.1049325) q[0];
sx q[0];
rz(1.2777591) q[0];
rz(-0.094206421) q[1];
sx q[1];
rz(-2.6779046) q[1];
sx q[1];
rz(1.5346079) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8746895) q[0];
sx q[0];
rz(-0.97123607) q[0];
sx q[0];
rz(0.039681704) q[0];
rz(-pi) q[1];
rz(-1.7530551) q[2];
sx q[2];
rz(-1.6267952) q[2];
sx q[2];
rz(0.1268498) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.31904768) q[1];
sx q[1];
rz(-1.3101868) q[1];
sx q[1];
rz(-2.6316363) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3393163) q[3];
sx q[3];
rz(-0.55219383) q[3];
sx q[3];
rz(1.4331762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.719912) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(-2.9812532) q[2];
rz(0.73873377) q[3];
sx q[3];
rz(-0.84638798) q[3];
sx q[3];
rz(1.7140478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.5147603) q[0];
sx q[0];
rz(-1.7034096) q[0];
sx q[0];
rz(-2.6318188) q[0];
rz(-1.431142) q[1];
sx q[1];
rz(-1.9606083) q[1];
sx q[1];
rz(-2.3950155) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1603161) q[0];
sx q[0];
rz(-2.123194) q[0];
sx q[0];
rz(-0.023560087) q[0];
rz(0.72991972) q[2];
sx q[2];
rz(-1.0906719) q[2];
sx q[2];
rz(1.4381222) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4164734) q[1];
sx q[1];
rz(-1.9320357) q[1];
sx q[1];
rz(-2.3010985) q[1];
x q[2];
rz(-1.1034562) q[3];
sx q[3];
rz(-0.73690542) q[3];
sx q[3];
rz(-1.74685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15472445) q[2];
sx q[2];
rz(-0.26686033) q[2];
sx q[2];
rz(1.9179087) q[2];
rz(-1.0736505) q[3];
sx q[3];
rz(-1.4370388) q[3];
sx q[3];
rz(0.8980208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680442) q[0];
sx q[0];
rz(-2.2548964) q[0];
sx q[0];
rz(0.37937382) q[0];
rz(0.1768449) q[1];
sx q[1];
rz(-1.481448) q[1];
sx q[1];
rz(0.79536974) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.11957) q[0];
sx q[0];
rz(-1.9384465) q[0];
sx q[0];
rz(-3.0421542) q[0];
rz(-pi) q[1];
rz(0.81669871) q[2];
sx q[2];
rz(-1.3743322) q[2];
sx q[2];
rz(2.8981371) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4702372) q[1];
sx q[1];
rz(-2.0908326) q[1];
sx q[1];
rz(-2.8440471) q[1];
rz(-pi) q[2];
rz(-1.627911) q[3];
sx q[3];
rz(-1.6962255) q[3];
sx q[3];
rz(-0.40298395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7103601) q[2];
sx q[2];
rz(-1.7962339) q[2];
sx q[2];
rz(-1.7809407) q[2];
rz(-2.1121173) q[3];
sx q[3];
rz(-1.6607213) q[3];
sx q[3];
rz(-1.3220538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65115702) q[0];
sx q[0];
rz(-1.54162) q[0];
sx q[0];
rz(1.5486451) q[0];
rz(-0.40052888) q[1];
sx q[1];
rz(-1.488204) q[1];
sx q[1];
rz(1.4097479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23435322) q[0];
sx q[0];
rz(-2.2924404) q[0];
sx q[0];
rz(-1.3035151) q[0];
rz(-pi) q[1];
rz(-2.7422041) q[2];
sx q[2];
rz(-2.7542973) q[2];
sx q[2];
rz(-3.0687817) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61297902) q[1];
sx q[1];
rz(-2.3720212) q[1];
sx q[1];
rz(-2.2320094) q[1];
rz(-pi) q[2];
rz(0.94878133) q[3];
sx q[3];
rz(-1.5153441) q[3];
sx q[3];
rz(-2.2353719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3132402) q[2];
sx q[2];
rz(-1.7007549) q[2];
sx q[2];
rz(-0.43133119) q[2];
rz(0.055015419) q[3];
sx q[3];
rz(-0.31156817) q[3];
sx q[3];
rz(-2.5855248) q[3];
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
rz(2.3738275) q[0];
sx q[0];
rz(-2.8887833) q[0];
sx q[0];
rz(2.1164236) q[0];
rz(-0.45285666) q[1];
sx q[1];
rz(-0.57944524) q[1];
sx q[1];
rz(-2.2377009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65176455) q[0];
sx q[0];
rz(-1.8790885) q[0];
sx q[0];
rz(-2.2405008) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3329817) q[2];
sx q[2];
rz(-0.8166733) q[2];
sx q[2];
rz(-2.6220235) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.55906193) q[1];
sx q[1];
rz(-0.85147038) q[1];
sx q[1];
rz(1.2079617) q[1];
x q[2];
rz(-0.76670209) q[3];
sx q[3];
rz(-0.23582102) q[3];
sx q[3];
rz(0.23877777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8338833) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(-0.21635381) q[2];
rz(2.2458535) q[3];
sx q[3];
rz(-2.4560865) q[3];
sx q[3];
rz(1.640813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0093339) q[0];
sx q[0];
rz(-2.0348771) q[0];
sx q[0];
rz(-2.4427781) q[0];
rz(1.9578594) q[1];
sx q[1];
rz(-1.7203169) q[1];
sx q[1];
rz(-2.2898477) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9601933) q[0];
sx q[0];
rz(-1.6087247) q[0];
sx q[0];
rz(1.0344124) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57735195) q[2];
sx q[2];
rz(-1.5492348) q[2];
sx q[2];
rz(-0.93761629) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.56543437) q[1];
sx q[1];
rz(-1.9322188) q[1];
sx q[1];
rz(0.71908497) q[1];
x q[2];
rz(2.5301039) q[3];
sx q[3];
rz(-0.73966714) q[3];
sx q[3];
rz(2.8304493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8375887) q[2];
sx q[2];
rz(-1.8303266) q[2];
sx q[2];
rz(-2.6336929) q[2];
rz(1.0287644) q[3];
sx q[3];
rz(-0.84380904) q[3];
sx q[3];
rz(2.540551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69865882) q[0];
sx q[0];
rz(-1.826094) q[0];
sx q[0];
rz(0.0096631924) q[0];
rz(-1.5254321) q[1];
sx q[1];
rz(-1.3776255) q[1];
sx q[1];
rz(-2.756871) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5949109) q[0];
sx q[0];
rz(-1.718504) q[0];
sx q[0];
rz(-1.3315931) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1037154) q[2];
sx q[2];
rz(-0.42755238) q[2];
sx q[2];
rz(0.84168692) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5009969) q[1];
sx q[1];
rz(-1.2944376) q[1];
sx q[1];
rz(-0.21548693) q[1];
rz(-pi) q[2];
x q[2];
rz(2.557005) q[3];
sx q[3];
rz(-0.86170022) q[3];
sx q[3];
rz(-0.60230909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5821417) q[2];
sx q[2];
rz(-0.95953512) q[2];
sx q[2];
rz(-0.01586308) q[2];
rz(-1.9150241) q[3];
sx q[3];
rz(-2.3358986) q[3];
sx q[3];
rz(-2.5198643) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3847619) q[0];
sx q[0];
rz(-2.4767196) q[0];
sx q[0];
rz(1.0913947) q[0];
rz(2.3233844) q[1];
sx q[1];
rz(-1.3312157) q[1];
sx q[1];
rz(0.6689201) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92907897) q[0];
sx q[0];
rz(-0.90796472) q[0];
sx q[0];
rz(0.79848358) q[0];
rz(1.0878272) q[2];
sx q[2];
rz(-0.7730128) q[2];
sx q[2];
rz(0.20394606) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.131625) q[1];
sx q[1];
rz(-2.5796081) q[1];
sx q[1];
rz(-1.704292) q[1];
x q[2];
rz(-2.8112765) q[3];
sx q[3];
rz(-1.013375) q[3];
sx q[3];
rz(-1.5239609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6341256) q[2];
sx q[2];
rz(-0.52906817) q[2];
sx q[2];
rz(2.1779306) q[2];
rz(1.4108747) q[3];
sx q[3];
rz(-1.4286634) q[3];
sx q[3];
rz(1.3655519) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1249579) q[0];
sx q[0];
rz(-0.5683012) q[0];
sx q[0];
rz(-1.6868663) q[0];
rz(-1.1085054) q[1];
sx q[1];
rz(-1.6051555) q[1];
sx q[1];
rz(-2.813521) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011400819) q[0];
sx q[0];
rz(-1.6284124) q[0];
sx q[0];
rz(-1.685623) q[0];
rz(-0.180822) q[2];
sx q[2];
rz(-1.4563113) q[2];
sx q[2];
rz(0.78071981) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1197443) q[1];
sx q[1];
rz(-1.9520063) q[1];
sx q[1];
rz(-0.66701835) q[1];
rz(-pi) q[2];
rz(1.397473) q[3];
sx q[3];
rz(-1.4651235) q[3];
sx q[3];
rz(2.8901576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46128094) q[2];
sx q[2];
rz(-1.842247) q[2];
sx q[2];
rz(2.7247562) q[2];
rz(0.69878116) q[3];
sx q[3];
rz(-2.4650033) q[3];
sx q[3];
rz(-2.7509287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9417435) q[0];
sx q[0];
rz(-1.9244292) q[0];
sx q[0];
rz(1.695965) q[0];
rz(-2.4327714) q[1];
sx q[1];
rz(-0.91239057) q[1];
sx q[1];
rz(-0.053587996) q[1];
rz(1.3016635) q[2];
sx q[2];
rz(-2.4872645) q[2];
sx q[2];
rz(-0.80254868) q[2];
rz(0.82984701) q[3];
sx q[3];
rz(-2.0086858) q[3];
sx q[3];
rz(0.84926844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
