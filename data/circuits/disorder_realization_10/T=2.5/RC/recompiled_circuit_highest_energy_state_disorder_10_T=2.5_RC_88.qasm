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
rz(-0.69775692) q[0];
sx q[0];
rz(-1.193576) q[0];
sx q[0];
rz(2.9370263) q[0];
rz(-0.76072955) q[1];
sx q[1];
rz(1.7525571) q[1];
sx q[1];
rz(11.166823) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8176048) q[0];
sx q[0];
rz(-1.7919375) q[0];
sx q[0];
rz(-0.093324109) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0206154) q[2];
sx q[2];
rz(-1.2024785) q[2];
sx q[2];
rz(-2.8415547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39069893) q[1];
sx q[1];
rz(-1.9913773) q[1];
sx q[1];
rz(1.4750359) q[1];
rz(-pi) q[2];
rz(-2.8329323) q[3];
sx q[3];
rz(-2.4081569) q[3];
sx q[3];
rz(1.4539469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43510398) q[2];
sx q[2];
rz(-0.032328345) q[2];
sx q[2];
rz(-2.6522563) q[2];
rz(-1.0232183) q[3];
sx q[3];
rz(-3.1232941) q[3];
sx q[3];
rz(2.0160915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-0.045687549) q[0];
sx q[0];
rz(-0.65653312) q[0];
sx q[0];
rz(-2.3387961) q[0];
rz(-0.071391694) q[1];
sx q[1];
rz(-0.26669058) q[1];
sx q[1];
rz(-0.057770483) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18081576) q[0];
sx q[0];
rz(-1.0313927) q[0];
sx q[0];
rz(2.8404854) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28654307) q[2];
sx q[2];
rz(-1.6725958) q[2];
sx q[2];
rz(1.2004146) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7747468) q[1];
sx q[1];
rz(-0.54301942) q[1];
sx q[1];
rz(-1.1096891) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3996734) q[3];
sx q[3];
rz(-2.4226396) q[3];
sx q[3];
rz(-2.3477682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.94566655) q[2];
sx q[2];
rz(-2.0464996) q[2];
sx q[2];
rz(1.2954953) q[2];
rz(0.96674353) q[3];
sx q[3];
rz(-0.77015489) q[3];
sx q[3];
rz(0.75743341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029723786) q[0];
sx q[0];
rz(-1.7787378) q[0];
sx q[0];
rz(-1.4335853) q[0];
rz(0.068610527) q[1];
sx q[1];
rz(-1.5741293) q[1];
sx q[1];
rz(-0.57919085) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60459701) q[0];
sx q[0];
rz(-1.662435) q[0];
sx q[0];
rz(0.058152278) q[0];
rz(-pi) q[1];
rz(0.4757232) q[2];
sx q[2];
rz(-1.5724725) q[2];
sx q[2];
rz(-2.7762108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9186185) q[1];
sx q[1];
rz(-2.9125305) q[1];
sx q[1];
rz(-2.983186) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3287401) q[3];
sx q[3];
rz(-0.62004706) q[3];
sx q[3];
rz(-0.64921415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27089831) q[2];
sx q[2];
rz(-1.2535932) q[2];
sx q[2];
rz(2.9581621) q[2];
rz(2.3373248) q[3];
sx q[3];
rz(-2.1427514) q[3];
sx q[3];
rz(-0.53406322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428228) q[0];
sx q[0];
rz(-0.11480055) q[0];
sx q[0];
rz(-0.59952366) q[0];
rz(2.7291258) q[1];
sx q[1];
rz(-3.1209374) q[1];
sx q[1];
rz(-1.0106769) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.855797) q[0];
sx q[0];
rz(-1.9205581) q[0];
sx q[0];
rz(1.4361022) q[0];
rz(-pi) q[1];
rz(0.050881906) q[2];
sx q[2];
rz(-2.2023099) q[2];
sx q[2];
rz(-2.8040407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0914826) q[1];
sx q[1];
rz(-0.98467876) q[1];
sx q[1];
rz(-0.32663235) q[1];
rz(0.75586163) q[3];
sx q[3];
rz(-1.6718277) q[3];
sx q[3];
rz(1.7441235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3707054) q[2];
sx q[2];
rz(-0.34660307) q[2];
sx q[2];
rz(2.8240805) q[2];
rz(-2.5192449) q[3];
sx q[3];
rz(-0.93572664) q[3];
sx q[3];
rz(-0.58421016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7790826) q[0];
sx q[0];
rz(-0.91479397) q[0];
sx q[0];
rz(2.1146178) q[0];
rz(-0.55038553) q[1];
sx q[1];
rz(-3.0774979) q[1];
sx q[1];
rz(1.9245573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7636722) q[0];
sx q[0];
rz(-0.91912133) q[0];
sx q[0];
rz(-2.4125189) q[0];
rz(-1.6486859) q[2];
sx q[2];
rz(-2.1078034) q[2];
sx q[2];
rz(-0.80918771) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0266307) q[1];
sx q[1];
rz(-1.3991881) q[1];
sx q[1];
rz(2.2725355) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86037614) q[3];
sx q[3];
rz(-1.851322) q[3];
sx q[3];
rz(-0.70395148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43651849) q[2];
sx q[2];
rz(-1.3631835) q[2];
sx q[2];
rz(-0.77318937) q[2];
rz(-3.0298722) q[3];
sx q[3];
rz(-1.3216647) q[3];
sx q[3];
rz(-2.2022061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47950995) q[0];
sx q[0];
rz(-0.22308068) q[0];
sx q[0];
rz(0.42698419) q[0];
rz(-0.93049479) q[1];
sx q[1];
rz(-3.1246779) q[1];
sx q[1];
rz(-0.46447909) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1980249) q[0];
sx q[0];
rz(-0.35145282) q[0];
sx q[0];
rz(0.60011421) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1709178) q[2];
sx q[2];
rz(-1.8541186) q[2];
sx q[2];
rz(-1.7194058) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0504004) q[1];
sx q[1];
rz(-0.37072966) q[1];
sx q[1];
rz(-3.0417473) q[1];
rz(-pi) q[2];
rz(1.5018015) q[3];
sx q[3];
rz(-2.7164) q[3];
sx q[3];
rz(2.0824279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6102607) q[2];
sx q[2];
rz(-1.320763) q[2];
sx q[2];
rz(2.8533234) q[2];
rz(2.098295) q[3];
sx q[3];
rz(-0.61560029) q[3];
sx q[3];
rz(-2.402795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4814602) q[0];
sx q[0];
rz(-1.2876502) q[0];
sx q[0];
rz(-0.75755358) q[0];
rz(-0.07269147) q[1];
sx q[1];
rz(-0.025765954) q[1];
sx q[1];
rz(-3.0933948) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91201997) q[0];
sx q[0];
rz(-2.1549112) q[0];
sx q[0];
rz(3.106227) q[0];
rz(-pi) q[1];
rz(-1.2343302) q[2];
sx q[2];
rz(-2.5162553) q[2];
sx q[2];
rz(-0.71695527) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0695677) q[1];
sx q[1];
rz(-1.9981019) q[1];
sx q[1];
rz(0.91520379) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4675619) q[3];
sx q[3];
rz(-2.6695731) q[3];
sx q[3];
rz(-1.9341521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4797719) q[2];
sx q[2];
rz(-1.4996108) q[2];
sx q[2];
rz(0.069395937) q[2];
rz(-1.5875459) q[3];
sx q[3];
rz(-0.78726751) q[3];
sx q[3];
rz(2.9230996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3875535) q[0];
sx q[0];
rz(-2.0856922) q[0];
sx q[0];
rz(1.4132502) q[0];
rz(2.3130401) q[1];
sx q[1];
rz(-3.0998402) q[1];
sx q[1];
rz(-2.5989596) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5722902) q[0];
sx q[0];
rz(-0.18768947) q[0];
sx q[0];
rz(-0.68861262) q[0];
x q[1];
rz(0.37658875) q[2];
sx q[2];
rz(-1.9053725) q[2];
sx q[2];
rz(-2.0221111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1664913) q[1];
sx q[1];
rz(-1.6246188) q[1];
sx q[1];
rz(3.0256998) q[1];
x q[2];
rz(-0.70809182) q[3];
sx q[3];
rz(-0.54881426) q[3];
sx q[3];
rz(-0.08591692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.98216206) q[2];
sx q[2];
rz(-2.7320778) q[2];
sx q[2];
rz(-2.8816667) q[2];
rz(-2.121117) q[3];
sx q[3];
rz(-2.8838938) q[3];
sx q[3];
rz(0.79403383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6771773) q[0];
sx q[0];
rz(-0.18615119) q[0];
sx q[0];
rz(-1.4790685) q[0];
rz(-1.5326477) q[1];
sx q[1];
rz(-1.0391935) q[1];
sx q[1];
rz(2.398568) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0940762) q[0];
sx q[0];
rz(-2.0819944) q[0];
sx q[0];
rz(3.009895) q[0];
x q[1];
rz(-2.1582339) q[2];
sx q[2];
rz(-0.62817103) q[2];
sx q[2];
rz(0.52631015) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2496532) q[1];
sx q[1];
rz(-1.6257451) q[1];
sx q[1];
rz(1.6366538) q[1];
rz(-pi) q[2];
rz(0.10492341) q[3];
sx q[3];
rz(-1.5399974) q[3];
sx q[3];
rz(1.9178176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6741901) q[2];
sx q[2];
rz(-0.82618606) q[2];
sx q[2];
rz(0.80545938) q[2];
rz(-1.6921267) q[3];
sx q[3];
rz(-1.2286681) q[3];
sx q[3];
rz(0.8031556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2529124) q[0];
sx q[0];
rz(-0.54179931) q[0];
sx q[0];
rz(2.3895277) q[0];
rz(-1.9883142) q[1];
sx q[1];
rz(-0.88429943) q[1];
sx q[1];
rz(0.28336743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8077791) q[0];
sx q[0];
rz(-2.5954843) q[0];
sx q[0];
rz(-2.9696652) q[0];
x q[1];
rz(2.9937075) q[2];
sx q[2];
rz(-0.82442946) q[2];
sx q[2];
rz(1.3651207) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9787406) q[1];
sx q[1];
rz(-1.8998441) q[1];
sx q[1];
rz(0.63905893) q[1];
rz(-pi) q[2];
x q[2];
rz(1.180319) q[3];
sx q[3];
rz(-2.1877398) q[3];
sx q[3];
rz(0.65333593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.89455426) q[2];
sx q[2];
rz(-3.0587695) q[2];
sx q[2];
rz(-1.7029597) q[2];
rz(2.8476207) q[3];
sx q[3];
rz(-3.1271264) q[3];
sx q[3];
rz(-1.0283874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.033584874) q[0];
sx q[0];
rz(-1.7243732) q[0];
sx q[0];
rz(1.617817) q[0];
rz(2.602018) q[1];
sx q[1];
rz(-2.3537666) q[1];
sx q[1];
rz(-2.981577) q[1];
rz(2.9931184) q[2];
sx q[2];
rz(-1.21429) q[2];
sx q[2];
rz(0.16333632) q[2];
rz(-2.5359375) q[3];
sx q[3];
rz(-0.45481053) q[3];
sx q[3];
rz(-0.74301471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
