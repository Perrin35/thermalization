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
rz(1.9873729) q[0];
sx q[0];
rz(-1.6183102) q[0];
sx q[0];
rz(-0.14259882) q[0];
rz(0.59347403) q[1];
sx q[1];
rz(-2.4886517) q[1];
sx q[1];
rz(-1.0452193) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7735239) q[0];
sx q[0];
rz(-2.0288101) q[0];
sx q[0];
rz(-0.59654838) q[0];
x q[1];
rz(-3.1047761) q[2];
sx q[2];
rz(-0.93339257) q[2];
sx q[2];
rz(2.3790145) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75946158) q[1];
sx q[1];
rz(-2.9540017) q[1];
sx q[1];
rz(-1.0008903) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7552283) q[3];
sx q[3];
rz(-2.5650555) q[3];
sx q[3];
rz(-0.54631305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1656701) q[2];
sx q[2];
rz(-2.133635) q[2];
sx q[2];
rz(-0.21318501) q[2];
rz(2.2255157) q[3];
sx q[3];
rz(-0.35917425) q[3];
sx q[3];
rz(0.38755125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85584545) q[0];
sx q[0];
rz(-0.88749945) q[0];
sx q[0];
rz(0.47072738) q[0];
rz(-2.0384906) q[1];
sx q[1];
rz(-2.6328937) q[1];
sx q[1];
rz(-0.57977605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.267201) q[0];
sx q[0];
rz(-1.6306595) q[0];
sx q[0];
rz(-0.52459985) q[0];
rz(1.2389555) q[2];
sx q[2];
rz(-2.08925) q[2];
sx q[2];
rz(1.4084852) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6439542) q[1];
sx q[1];
rz(-0.71566391) q[1];
sx q[1];
rz(0.96244754) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1768119) q[3];
sx q[3];
rz(-0.44108118) q[3];
sx q[3];
rz(1.6562605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9819928) q[2];
sx q[2];
rz(-2.2856568) q[2];
sx q[2];
rz(-3.0465872) q[2];
rz(0.57560086) q[3];
sx q[3];
rz(-2.3592981) q[3];
sx q[3];
rz(-1.6949722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59548241) q[0];
sx q[0];
rz(-0.70850104) q[0];
sx q[0];
rz(0.27221671) q[0];
rz(3.0369924) q[1];
sx q[1];
rz(-1.0403386) q[1];
sx q[1];
rz(-1.4629755) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7143216) q[0];
sx q[0];
rz(-0.84762379) q[0];
sx q[0];
rz(0.066825213) q[0];
rz(-pi) q[1];
x q[1];
rz(0.06426364) q[2];
sx q[2];
rz(-2.3552364) q[2];
sx q[2];
rz(-0.83640316) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.453664) q[1];
sx q[1];
rz(-1.2896104) q[1];
sx q[1];
rz(-2.4635876) q[1];
rz(-pi) q[2];
rz(-1.6312863) q[3];
sx q[3];
rz(-1.20245) q[3];
sx q[3];
rz(-0.18243901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0766803) q[2];
sx q[2];
rz(-1.400759) q[2];
sx q[2];
rz(-2.7317375) q[2];
rz(-3.0900433) q[3];
sx q[3];
rz(-0.99286538) q[3];
sx q[3];
rz(-0.73369098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9728397) q[0];
sx q[0];
rz(-2.3641455) q[0];
sx q[0];
rz(-2.4421316) q[0];
rz(-2.2109924) q[1];
sx q[1];
rz(-1.365463) q[1];
sx q[1];
rz(2.0420989) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4874026) q[0];
sx q[0];
rz(-0.64399566) q[0];
sx q[0];
rz(-2.0430768) q[0];
rz(-pi) q[1];
rz(2.1191988) q[2];
sx q[2];
rz(-1.7779997) q[2];
sx q[2];
rz(-2.0275379) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1395542) q[1];
sx q[1];
rz(-2.6894662) q[1];
sx q[1];
rz(0.76805784) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9774262) q[3];
sx q[3];
rz(-2.1726131) q[3];
sx q[3];
rz(0.072657771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76116556) q[2];
sx q[2];
rz(-2.0695504) q[2];
sx q[2];
rz(-0.94044828) q[2];
rz(-2.579651) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(2.565062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1171653) q[0];
sx q[0];
rz(-3.0553525) q[0];
sx q[0];
rz(-1.0166136) q[0];
rz(2.2158465) q[1];
sx q[1];
rz(-0.75030202) q[1];
sx q[1];
rz(-0.076676682) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8228418) q[0];
sx q[0];
rz(-2.7116576) q[0];
sx q[0];
rz(-1.1640401) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3252649) q[2];
sx q[2];
rz(-1.7707728) q[2];
sx q[2];
rz(2.1918346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8342585) q[1];
sx q[1];
rz(-1.7821719) q[1];
sx q[1];
rz(1.0049051) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2604976) q[3];
sx q[3];
rz(-0.7782794) q[3];
sx q[3];
rz(0.33184127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3735247) q[2];
sx q[2];
rz(-2.3536451) q[2];
sx q[2];
rz(-3.0184271) q[2];
rz(-2.4672616) q[3];
sx q[3];
rz(-0.47334039) q[3];
sx q[3];
rz(0.82790747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3097565) q[0];
sx q[0];
rz(-2.9530544) q[0];
sx q[0];
rz(-2.6605666) q[0];
rz(-2.9006529) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(-0.13490881) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3137993) q[0];
sx q[0];
rz(-2.302357) q[0];
sx q[0];
rz(0.93904943) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1730072) q[2];
sx q[2];
rz(-2.5998625) q[2];
sx q[2];
rz(2.2587905) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1926769) q[1];
sx q[1];
rz(-0.8444311) q[1];
sx q[1];
rz(2.7706258) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14860146) q[3];
sx q[3];
rz(-0.45431787) q[3];
sx q[3];
rz(-3.076072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9350819) q[2];
sx q[2];
rz(-0.93588459) q[2];
sx q[2];
rz(1.3235462) q[2];
rz(2.8694618) q[3];
sx q[3];
rz(-2.2836253) q[3];
sx q[3];
rz(0.14565295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018933522) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(2.4830699) q[0];
rz(1.5497442) q[1];
sx q[1];
rz(-2.1287983) q[1];
sx q[1];
rz(-0.50643593) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69909912) q[0];
sx q[0];
rz(-0.37537071) q[0];
sx q[0];
rz(-0.18113776) q[0];
rz(-0.030722458) q[2];
sx q[2];
rz(-2.0523768) q[2];
sx q[2];
rz(2.5831163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5750908) q[1];
sx q[1];
rz(-3.021214) q[1];
sx q[1];
rz(-0.68415227) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.562256) q[3];
sx q[3];
rz(-1.6054573) q[3];
sx q[3];
rz(-0.12753419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3353614) q[2];
sx q[2];
rz(-2.3956617) q[2];
sx q[2];
rz(-1.8719505) q[2];
rz(-1.3658203) q[3];
sx q[3];
rz(-3.1277872) q[3];
sx q[3];
rz(0.55240101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34614554) q[0];
sx q[0];
rz(-2.8988291) q[0];
sx q[0];
rz(-2.7224097) q[0];
rz(-3.0508793) q[1];
sx q[1];
rz(-2.2234629) q[1];
sx q[1];
rz(-1.027164) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5266383) q[0];
sx q[0];
rz(-1.1246343) q[0];
sx q[0];
rz(2.5102708) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33454169) q[2];
sx q[2];
rz(-3.1174264) q[2];
sx q[2];
rz(-0.4949567) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.724714) q[1];
sx q[1];
rz(-1.4644578) q[1];
sx q[1];
rz(1.7843397) q[1];
rz(-pi) q[2];
rz(-3.0869802) q[3];
sx q[3];
rz(-2.255618) q[3];
sx q[3];
rz(-0.38338654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1586228) q[2];
sx q[2];
rz(-1.0047487) q[2];
sx q[2];
rz(-0.99009222) q[2];
rz(2.8236735) q[3];
sx q[3];
rz(-2.3293994) q[3];
sx q[3];
rz(-3.1032491) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7541499) q[0];
sx q[0];
rz(-0.70873547) q[0];
sx q[0];
rz(2.5842174) q[0];
rz(-0.36224657) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(-0.16709669) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.635256) q[0];
sx q[0];
rz(-1.8765159) q[0];
sx q[0];
rz(1.8503667) q[0];
rz(-pi) q[1];
rz(-1.7542776) q[2];
sx q[2];
rz(-2.263732) q[2];
sx q[2];
rz(2.3132581) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9145979) q[1];
sx q[1];
rz(-2.216504) q[1];
sx q[1];
rz(1.1832994) q[1];
rz(-pi) q[2];
rz(1.4792419) q[3];
sx q[3];
rz(-0.6855489) q[3];
sx q[3];
rz(2.8753672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.96357137) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(0.42845217) q[2];
rz(2.4109449) q[3];
sx q[3];
rz(-1.080039) q[3];
sx q[3];
rz(0.60034269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44883248) q[0];
sx q[0];
rz(-1.4880143) q[0];
sx q[0];
rz(-2.0220508) q[0];
rz(2.4730832) q[1];
sx q[1];
rz(-0.57066494) q[1];
sx q[1];
rz(-2.8817435) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.665312) q[0];
sx q[0];
rz(-2.5710433) q[0];
sx q[0];
rz(-1.0744988) q[0];
rz(-pi) q[1];
rz(2.75801) q[2];
sx q[2];
rz(-2.0077939) q[2];
sx q[2];
rz(2.008977) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6653241) q[1];
sx q[1];
rz(-2.5905053) q[1];
sx q[1];
rz(2.4399806) q[1];
x q[2];
rz(2.9838461) q[3];
sx q[3];
rz(-2.254527) q[3];
sx q[3];
rz(-2.5585548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6707637) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(-1.0864351) q[2];
rz(-0.49232617) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(-0.72211784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54871854) q[0];
sx q[0];
rz(-1.508779) q[0];
sx q[0];
rz(-1.4887703) q[0];
rz(-1.8863574) q[1];
sx q[1];
rz(-1.8240758) q[1];
sx q[1];
rz(0.12771894) q[1];
rz(-1.7931425) q[2];
sx q[2];
rz(-2.7334474) q[2];
sx q[2];
rz(1.853142) q[2];
rz(0.90821785) q[3];
sx q[3];
rz(-1.5208018) q[3];
sx q[3];
rz(3.0326299) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
