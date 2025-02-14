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
rz(1.4242564) q[0];
sx q[0];
rz(1.9219226) q[0];
sx q[0];
rz(9.9998247) q[0];
rz(0.14248928) q[1];
sx q[1];
rz(-1.9228851) q[1];
sx q[1];
rz(-0.86427468) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40022093) q[0];
sx q[0];
rz(-1.7622041) q[0];
sx q[0];
rz(-1.7975397) q[0];
rz(-1.7792542) q[2];
sx q[2];
rz(-0.64178665) q[2];
sx q[2];
rz(-1.2520977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4839403) q[1];
sx q[1];
rz(-0.53002702) q[1];
sx q[1];
rz(0.0040667314) q[1];
rz(-pi) q[2];
rz(-0.56700403) q[3];
sx q[3];
rz(-1.1409014) q[3];
sx q[3];
rz(0.56741949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.076866604) q[2];
sx q[2];
rz(-1.4996424) q[2];
sx q[2];
rz(-2.1323252) q[2];
rz(0.81389728) q[3];
sx q[3];
rz(-0.32865694) q[3];
sx q[3];
rz(-2.8024659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0547884) q[0];
sx q[0];
rz(-1.4084933) q[0];
sx q[0];
rz(-0.24222294) q[0];
rz(1.9107266) q[1];
sx q[1];
rz(-2.5098398) q[1];
sx q[1];
rz(-2.3435074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7042328) q[0];
sx q[0];
rz(-0.95600545) q[0];
sx q[0];
rz(-1.4324709) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1801486) q[2];
sx q[2];
rz(-2.7728348) q[2];
sx q[2];
rz(0.40290305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68951571) q[1];
sx q[1];
rz(-0.59714666) q[1];
sx q[1];
rz(-2.4680016) q[1];
rz(0.35817082) q[3];
sx q[3];
rz(-1.4205975) q[3];
sx q[3];
rz(-1.1800486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.65608394) q[2];
sx q[2];
rz(-2.0723497) q[2];
sx q[2];
rz(2.3089224) q[2];
rz(2.5101856) q[3];
sx q[3];
rz(-2.4000945) q[3];
sx q[3];
rz(-3.1413063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0176508) q[0];
sx q[0];
rz(-2.232382) q[0];
sx q[0];
rz(2.9538474) q[0];
rz(0.84367696) q[1];
sx q[1];
rz(-1.8673106) q[1];
sx q[1];
rz(2.5086596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62719856) q[0];
sx q[0];
rz(-2.7717675) q[0];
sx q[0];
rz(1.373795) q[0];
x q[1];
rz(1.9203414) q[2];
sx q[2];
rz(-2.2270906) q[2];
sx q[2];
rz(-1.2906769) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.570185) q[1];
sx q[1];
rz(-1.270789) q[1];
sx q[1];
rz(2.6064998) q[1];
rz(1.9162779) q[3];
sx q[3];
rz(-2.7746183) q[3];
sx q[3];
rz(-0.93852321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26744276) q[2];
sx q[2];
rz(-2.1953857) q[2];
sx q[2];
rz(-1.2713185) q[2];
rz(1.6455796) q[3];
sx q[3];
rz(-1.4964024) q[3];
sx q[3];
rz(2.1700844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8895759) q[0];
sx q[0];
rz(-2.3887964) q[0];
sx q[0];
rz(2.4821607) q[0];
rz(-1.8239498) q[1];
sx q[1];
rz(-1.3392071) q[1];
sx q[1];
rz(0.79016322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58604974) q[0];
sx q[0];
rz(-1.9056742) q[0];
sx q[0];
rz(-1.4396776) q[0];
rz(-1.538058) q[2];
sx q[2];
rz(-1.5509836) q[2];
sx q[2];
rz(2.0231501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9875264) q[1];
sx q[1];
rz(-0.48168698) q[1];
sx q[1];
rz(-2.9996526) q[1];
rz(0.55410093) q[3];
sx q[3];
rz(-1.3375207) q[3];
sx q[3];
rz(1.3194039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0822175) q[2];
sx q[2];
rz(-1.5267812) q[2];
sx q[2];
rz(0.29602948) q[2];
rz(2.5791903) q[3];
sx q[3];
rz(-0.86807576) q[3];
sx q[3];
rz(-3.0276022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.247308) q[0];
sx q[0];
rz(-1.9534651) q[0];
sx q[0];
rz(-2.9344015) q[0];
rz(-2.1098792) q[1];
sx q[1];
rz(-1.9887911) q[1];
sx q[1];
rz(-1.656104) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15940672) q[0];
sx q[0];
rz(-2.3362118) q[0];
sx q[0];
rz(-2.5872487) q[0];
x q[1];
rz(-1.0096512) q[2];
sx q[2];
rz(-1.818294) q[2];
sx q[2];
rz(-1.1010608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7149026) q[1];
sx q[1];
rz(-1.443092) q[1];
sx q[1];
rz(-2.7903778) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4506102) q[3];
sx q[3];
rz(-0.73420364) q[3];
sx q[3];
rz(2.7199573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.038387211) q[2];
sx q[2];
rz(-0.81192553) q[2];
sx q[2];
rz(-1.1078328) q[2];
rz(-2.2191018) q[3];
sx q[3];
rz(-1.4100217) q[3];
sx q[3];
rz(-2.8973575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6237727) q[0];
sx q[0];
rz(-2.6277442) q[0];
sx q[0];
rz(2.9129831) q[0];
rz(2.8672583) q[1];
sx q[1];
rz(-0.81293303) q[1];
sx q[1];
rz(-2.6913604) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73295702) q[0];
sx q[0];
rz(-1.3191603) q[0];
sx q[0];
rz(0.77045124) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8000431) q[2];
sx q[2];
rz(-1.563213) q[2];
sx q[2];
rz(-2.6784865) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95388639) q[1];
sx q[1];
rz(-1.0774195) q[1];
sx q[1];
rz(1.8567371) q[1];
rz(-pi) q[2];
rz(2.8725876) q[3];
sx q[3];
rz(-1.169636) q[3];
sx q[3];
rz(-1.3022193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36394128) q[2];
sx q[2];
rz(-0.55299091) q[2];
sx q[2];
rz(-0.027776329) q[2];
rz(2.3751496) q[3];
sx q[3];
rz(-1.1602297) q[3];
sx q[3];
rz(0.69507039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9195093) q[0];
sx q[0];
rz(-1.5064025) q[0];
sx q[0];
rz(-2.5191504) q[0];
rz(-0.70603236) q[1];
sx q[1];
rz(-2.249554) q[1];
sx q[1];
rz(-2.5546254) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0737105) q[0];
sx q[0];
rz(-2.1176) q[0];
sx q[0];
rz(2.6220462) q[0];
rz(-pi) q[1];
rz(1.3038425) q[2];
sx q[2];
rz(-1.0994229) q[2];
sx q[2];
rz(1.4990569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1026969) q[1];
sx q[1];
rz(-1.4998815) q[1];
sx q[1];
rz(-2.6346579) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5587646) q[3];
sx q[3];
rz(-0.71206743) q[3];
sx q[3];
rz(0.67675096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.065757699) q[2];
sx q[2];
rz(-1.4328052) q[2];
sx q[2];
rz(-2.5471121) q[2];
rz(-2.8822656) q[3];
sx q[3];
rz(-1.8766873) q[3];
sx q[3];
rz(1.9243141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-3.1031621) q[0];
sx q[0];
rz(-2.2362464) q[0];
sx q[0];
rz(0.578798) q[0];
rz(1.831306) q[1];
sx q[1];
rz(-1.879004) q[1];
sx q[1];
rz(-2.328918) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2392501) q[0];
sx q[0];
rz(-2.6391811) q[0];
sx q[0];
rz(2.8397296) q[0];
rz(2.8825106) q[2];
sx q[2];
rz(-1.0244842) q[2];
sx q[2];
rz(-1.3457042) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95083773) q[1];
sx q[1];
rz(-0.79154888) q[1];
sx q[1];
rz(2.6218824) q[1];
x q[2];
rz(2.1401494) q[3];
sx q[3];
rz(-1.1318996) q[3];
sx q[3];
rz(-2.2341408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3506713) q[2];
sx q[2];
rz(-1.2988043) q[2];
sx q[2];
rz(-0.92362967) q[2];
rz(-2.1212497) q[3];
sx q[3];
rz(-0.89710051) q[3];
sx q[3];
rz(0.11788192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0807121) q[0];
sx q[0];
rz(-1.7055644) q[0];
sx q[0];
rz(1.4870148) q[0];
rz(0.05038536) q[1];
sx q[1];
rz(-2.3245508) q[1];
sx q[1];
rz(2.494716) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9084307) q[0];
sx q[0];
rz(-1.1457503) q[0];
sx q[0];
rz(2.2320497) q[0];
rz(-2.4755461) q[2];
sx q[2];
rz(-1.9589309) q[2];
sx q[2];
rz(-1.9789517) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.17740956) q[1];
sx q[1];
rz(-0.46176592) q[1];
sx q[1];
rz(2.6564471) q[1];
x q[2];
rz(2.991363) q[3];
sx q[3];
rz(-1.3597128) q[3];
sx q[3];
rz(-2.7631813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6369624) q[2];
sx q[2];
rz(-1.5865822) q[2];
sx q[2];
rz(-0.75616765) q[2];
rz(1.470083) q[3];
sx q[3];
rz(-0.78592891) q[3];
sx q[3];
rz(0.80529958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037755448) q[0];
sx q[0];
rz(-1.8917731) q[0];
sx q[0];
rz(2.0334429) q[0];
rz(0.26421079) q[1];
sx q[1];
rz(-1.6725531) q[1];
sx q[1];
rz(-0.60233751) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36985227) q[0];
sx q[0];
rz(-0.30769545) q[0];
sx q[0];
rz(1.0956531) q[0];
x q[1];
rz(2.9649023) q[2];
sx q[2];
rz(-1.7222705) q[2];
sx q[2];
rz(0.02515153) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8089167) q[1];
sx q[1];
rz(-2.4530468) q[1];
sx q[1];
rz(1.4946412) q[1];
x q[2];
rz(-0.65091316) q[3];
sx q[3];
rz(-1.3861361) q[3];
sx q[3];
rz(-2.1391275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5084874) q[2];
sx q[2];
rz(-0.54163951) q[2];
sx q[2];
rz(-2.2701021) q[2];
rz(-1.9320711) q[3];
sx q[3];
rz(-2.8612374) q[3];
sx q[3];
rz(-0.59396321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7021983) q[0];
sx q[0];
rz(-1.7560503) q[0];
sx q[0];
rz(-0.51881292) q[0];
rz(0.58204542) q[1];
sx q[1];
rz(-2.0379635) q[1];
sx q[1];
rz(1.4720974) q[1];
rz(-2.4051278) q[2];
sx q[2];
rz(-1.2869375) q[2];
sx q[2];
rz(1.7766042) q[2];
rz(0.68861674) q[3];
sx q[3];
rz(-1.4388765) q[3];
sx q[3];
rz(2.8394113) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
