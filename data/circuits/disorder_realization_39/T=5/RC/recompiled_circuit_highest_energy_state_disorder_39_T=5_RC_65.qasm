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
rz(-3.0012335) q[0];
sx q[0];
rz(-1.4599414) q[0];
sx q[0];
rz(-0.85293823) q[0];
rz(1.9536904) q[1];
sx q[1];
rz(3.3878769) q[1];
sx q[1];
rz(9.090957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1000017) q[0];
sx q[0];
rz(-2.0717588) q[0];
sx q[0];
rz(0.28613018) q[0];
rz(-pi) q[1];
rz(-0.63392459) q[2];
sx q[2];
rz(-1.156154) q[2];
sx q[2];
rz(-2.7716605) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8454518) q[1];
sx q[1];
rz(-0.23794623) q[1];
sx q[1];
rz(0.68745698) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4859441) q[3];
sx q[3];
rz(-1.6716262) q[3];
sx q[3];
rz(-2.3840897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21981257) q[2];
sx q[2];
rz(-2.2647936) q[2];
sx q[2];
rz(-0.63966695) q[2];
rz(-2.2031247) q[3];
sx q[3];
rz(-2.7191021) q[3];
sx q[3];
rz(-1.870702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4939782) q[0];
sx q[0];
rz(-3.0632126) q[0];
sx q[0];
rz(-1.6139503) q[0];
rz(2.899462) q[1];
sx q[1];
rz(-2.1277728) q[1];
sx q[1];
rz(-2.4193144) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15597181) q[0];
sx q[0];
rz(-1.9160144) q[0];
sx q[0];
rz(1.9916608) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79485915) q[2];
sx q[2];
rz(-1.7599918) q[2];
sx q[2];
rz(-2.4986588) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7178044) q[1];
sx q[1];
rz(-0.77149888) q[1];
sx q[1];
rz(-0.81591925) q[1];
x q[2];
rz(-2.9728209) q[3];
sx q[3];
rz(-0.95484551) q[3];
sx q[3];
rz(1.4142766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.051141288) q[2];
sx q[2];
rz(-2.800056) q[2];
sx q[2];
rz(3.0549468) q[2];
rz(-0.58049774) q[3];
sx q[3];
rz(-1.9167506) q[3];
sx q[3];
rz(-2.1897924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8101623) q[0];
sx q[0];
rz(-2.0482752) q[0];
sx q[0];
rz(-1.4587559) q[0];
rz(-3.0259865) q[1];
sx q[1];
rz(-1.9435725) q[1];
sx q[1];
rz(-0.16673949) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78555303) q[0];
sx q[0];
rz(-2.5305439) q[0];
sx q[0];
rz(-1.9370228) q[0];
rz(-1.9405792) q[2];
sx q[2];
rz(-1.175011) q[2];
sx q[2];
rz(-1.41768) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3003242) q[1];
sx q[1];
rz(-1.027123) q[1];
sx q[1];
rz(2.6946696) q[1];
rz(-pi) q[2];
rz(-2.3963967) q[3];
sx q[3];
rz(-1.4663854) q[3];
sx q[3];
rz(2.1946583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.705767) q[2];
sx q[2];
rz(-2.9220118) q[2];
sx q[2];
rz(1.1337918) q[2];
rz(-3.0886768) q[3];
sx q[3];
rz(-1.8688801) q[3];
sx q[3];
rz(2.4165966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13919203) q[0];
sx q[0];
rz(-0.63183689) q[0];
sx q[0];
rz(2.8644526) q[0];
rz(-2.3587522) q[1];
sx q[1];
rz(-1.506184) q[1];
sx q[1];
rz(-0.35071075) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6576516) q[0];
sx q[0];
rz(-2.4336877) q[0];
sx q[0];
rz(1.4483676) q[0];
rz(-pi) q[1];
rz(3.0530351) q[2];
sx q[2];
rz(-0.32677256) q[2];
sx q[2];
rz(1.8984924) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8102032) q[1];
sx q[1];
rz(-2.1955829) q[1];
sx q[1];
rz(2.5532755) q[1];
x q[2];
rz(0.043204149) q[3];
sx q[3];
rz(-1.4625878) q[3];
sx q[3];
rz(1.4747335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1002525) q[2];
sx q[2];
rz(-1.2790054) q[2];
sx q[2];
rz(1.9913541) q[2];
rz(2.9669115) q[3];
sx q[3];
rz(-0.29398578) q[3];
sx q[3];
rz(3.0512419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1133076) q[0];
sx q[0];
rz(-2.568013) q[0];
sx q[0];
rz(-1.9453402) q[0];
rz(-2.664227) q[1];
sx q[1];
rz(-0.82672516) q[1];
sx q[1];
rz(-0.54944077) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4065721) q[0];
sx q[0];
rz(-1.0112678) q[0];
sx q[0];
rz(-2.4491549) q[0];
rz(0.35713335) q[2];
sx q[2];
rz(-2.5741842) q[2];
sx q[2];
rz(1.0536989) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6345646) q[1];
sx q[1];
rz(-1.1643049) q[1];
sx q[1];
rz(-0.1312934) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1138238) q[3];
sx q[3];
rz(-1.9544056) q[3];
sx q[3];
rz(-2.3758604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4598733) q[2];
sx q[2];
rz(-2.7698066) q[2];
sx q[2];
rz(-2.8968774) q[2];
rz(-1.5185482) q[3];
sx q[3];
rz(-2.0270429) q[3];
sx q[3];
rz(-3.0486619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86843425) q[0];
sx q[0];
rz(-2.1491829) q[0];
sx q[0];
rz(2.7145845) q[0];
rz(1.2097516) q[1];
sx q[1];
rz(-1.6170343) q[1];
sx q[1];
rz(-0.0078113656) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2526379) q[0];
sx q[0];
rz(-1.7868306) q[0];
sx q[0];
rz(0.81058575) q[0];
x q[1];
rz(3.1321944) q[2];
sx q[2];
rz(-1.0236866) q[2];
sx q[2];
rz(-1.6449606) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.58578929) q[1];
sx q[1];
rz(-2.2206056) q[1];
sx q[1];
rz(1.5073677) q[1];
rz(0.2711556) q[3];
sx q[3];
rz(-0.89224766) q[3];
sx q[3];
rz(-3.1377047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4012287) q[2];
sx q[2];
rz(-0.87330356) q[2];
sx q[2];
rz(2.002423) q[2];
rz(-0.27979699) q[3];
sx q[3];
rz(-1.3236902) q[3];
sx q[3];
rz(-1.4626224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6497659) q[0];
sx q[0];
rz(-2.7597646) q[0];
sx q[0];
rz(-2.5873798) q[0];
rz(-3.0795433) q[1];
sx q[1];
rz(-2.4833312) q[1];
sx q[1];
rz(-2.2668692) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59959182) q[0];
sx q[0];
rz(-0.79757798) q[0];
sx q[0];
rz(2.9717658) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1360097) q[2];
sx q[2];
rz(-1.5632331) q[2];
sx q[2];
rz(2.7420704) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3530423) q[1];
sx q[1];
rz(-1.0354831) q[1];
sx q[1];
rz(-2.8393406) q[1];
rz(-1.0933541) q[3];
sx q[3];
rz(-0.72243172) q[3];
sx q[3];
rz(-2.5337766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41500652) q[2];
sx q[2];
rz(-1.0174624) q[2];
sx q[2];
rz(0.93878186) q[2];
rz(0.69432652) q[3];
sx q[3];
rz(-2.4735579) q[3];
sx q[3];
rz(-0.15650775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.071455) q[0];
sx q[0];
rz(-2.5818765) q[0];
sx q[0];
rz(2.8330084) q[0];
rz(-1.4831108) q[1];
sx q[1];
rz(-1.3456656) q[1];
sx q[1];
rz(2.9187091) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5102986) q[0];
sx q[0];
rz(-1.1985072) q[0];
sx q[0];
rz(2.3907803) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7168458) q[2];
sx q[2];
rz(-0.24492376) q[2];
sx q[2];
rz(-2.4670759) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.4207238) q[1];
sx q[1];
rz(-1.6096186) q[1];
sx q[1];
rz(-1.8160519) q[1];
x q[2];
rz(-2.127366) q[3];
sx q[3];
rz(-2.1750919) q[3];
sx q[3];
rz(-0.44318553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8181204) q[2];
sx q[2];
rz(-1.0264531) q[2];
sx q[2];
rz(0.66413122) q[2];
rz(-1.1431665) q[3];
sx q[3];
rz(-1.6616471) q[3];
sx q[3];
rz(2.0460879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42300647) q[0];
sx q[0];
rz(-1.9976595) q[0];
sx q[0];
rz(-0.018420694) q[0];
rz(-2.9711235) q[1];
sx q[1];
rz(-1.8959931) q[1];
sx q[1];
rz(-2.9439994) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.402408) q[0];
sx q[0];
rz(-1.5801593) q[0];
sx q[0];
rz(-1.4020756) q[0];
x q[1];
rz(-3.0891916) q[2];
sx q[2];
rz(-1.9366855) q[2];
sx q[2];
rz(3.0405557) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9925633) q[1];
sx q[1];
rz(-2.8125893) q[1];
sx q[1];
rz(2.002153) q[1];
rz(-pi) q[2];
rz(-0.13473265) q[3];
sx q[3];
rz(-2.7425457) q[3];
sx q[3];
rz(1.4302554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.0240872) q[2];
sx q[2];
rz(-2.0427637) q[2];
sx q[2];
rz(0.39829028) q[2];
rz(2.6309218) q[3];
sx q[3];
rz(-1.6528249) q[3];
sx q[3];
rz(-0.85810703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2354105) q[0];
sx q[0];
rz(-2.372083) q[0];
sx q[0];
rz(0.19943516) q[0];
rz(0.27345744) q[1];
sx q[1];
rz(-2.7077935) q[1];
sx q[1];
rz(-2.1308897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0293297) q[0];
sx q[0];
rz(-2.4284017) q[0];
sx q[0];
rz(0.40644706) q[0];
rz(-pi) q[1];
x q[1];
rz(1.550714) q[2];
sx q[2];
rz(-1.806136) q[2];
sx q[2];
rz(2.2724336) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1361183) q[1];
sx q[1];
rz(-1.6228383) q[1];
sx q[1];
rz(2.4050557) q[1];
x q[2];
rz(2.3673986) q[3];
sx q[3];
rz(-2.7280118) q[3];
sx q[3];
rz(-2.3516097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79591862) q[2];
sx q[2];
rz(-2.5359539) q[2];
sx q[2];
rz(-0.79552135) q[2];
rz(0.37122053) q[3];
sx q[3];
rz(-1.3859387) q[3];
sx q[3];
rz(-2.8871239) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026074792) q[0];
sx q[0];
rz(-1.4160897) q[0];
sx q[0];
rz(-1.8427451) q[0];
rz(-2.5904291) q[1];
sx q[1];
rz(-1.1977341) q[1];
sx q[1];
rz(-1.151998) q[1];
rz(-0.72533561) q[2];
sx q[2];
rz(-1.7400405) q[2];
sx q[2];
rz(-1.7766458) q[2];
rz(1.1887278) q[3];
sx q[3];
rz(-1.3466866) q[3];
sx q[3];
rz(-1.8849296) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
