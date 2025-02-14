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
rz(-2.0853618) q[0];
sx q[0];
rz(-0.27529278) q[0];
sx q[0];
rz(-3.0866301) q[0];
rz(-2.3995402) q[1];
sx q[1];
rz(-0.11062515) q[1];
sx q[1];
rz(-0.7551809) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7464712) q[0];
sx q[0];
rz(-0.81534213) q[0];
sx q[0];
rz(1.59207) q[0];
rz(-pi) q[1];
rz(1.0805425) q[2];
sx q[2];
rz(-0.95455652) q[2];
sx q[2];
rz(-2.3927671) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8557737) q[1];
sx q[1];
rz(-2.710418) q[1];
sx q[1];
rz(-1.5613501) q[1];
rz(-pi) q[2];
rz(-2.5489786) q[3];
sx q[3];
rz(-1.7120769) q[3];
sx q[3];
rz(-1.1326552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5893958) q[2];
sx q[2];
rz(-1.155921) q[2];
sx q[2];
rz(-1.4001621) q[2];
rz(0.33341148) q[3];
sx q[3];
rz(-2.4937778) q[3];
sx q[3];
rz(0.49745542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19631504) q[0];
sx q[0];
rz(-0.59756398) q[0];
sx q[0];
rz(-3.0481098) q[0];
rz(0.64158332) q[1];
sx q[1];
rz(-2.7020794) q[1];
sx q[1];
rz(-2.9494367) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7948611) q[0];
sx q[0];
rz(-2.4114176) q[0];
sx q[0];
rz(-0.037564354) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3619611) q[2];
sx q[2];
rz(-1.3614839) q[2];
sx q[2];
rz(-2.9787082) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.007602674) q[1];
sx q[1];
rz(-0.96293425) q[1];
sx q[1];
rz(-2.1278343) q[1];
rz(2.8175958) q[3];
sx q[3];
rz(-1.5662875) q[3];
sx q[3];
rz(1.5229932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41500348) q[2];
sx q[2];
rz(-0.64565349) q[2];
sx q[2];
rz(-2.1356706) q[2];
rz(3.1065324) q[3];
sx q[3];
rz(-0.56060767) q[3];
sx q[3];
rz(-2.3817271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9525035) q[0];
sx q[0];
rz(-1.8006421) q[0];
sx q[0];
rz(2.9553318) q[0];
rz(2.8798036) q[1];
sx q[1];
rz(-1.9354542) q[1];
sx q[1];
rz(0.50938767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4550393) q[0];
sx q[0];
rz(-1.5709366) q[0];
sx q[0];
rz(-1.554398) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7942717) q[2];
sx q[2];
rz(-2.2185549) q[2];
sx q[2];
rz(-1.8504686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68056831) q[1];
sx q[1];
rz(-2.1740434) q[1];
sx q[1];
rz(1.6955743) q[1];
rz(-pi) q[2];
rz(-2.4136606) q[3];
sx q[3];
rz(-1.0714487) q[3];
sx q[3];
rz(3.0161757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.95106) q[2];
sx q[2];
rz(-0.6742) q[2];
sx q[2];
rz(0.29779693) q[2];
rz(2.8262302) q[3];
sx q[3];
rz(-2.8331929) q[3];
sx q[3];
rz(-1.7409918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067658871) q[0];
sx q[0];
rz(-0.44026259) q[0];
sx q[0];
rz(-0.0010781188) q[0];
rz(0.06238097) q[1];
sx q[1];
rz(-1.9900813) q[1];
sx q[1];
rz(1.4952225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8366458) q[0];
sx q[0];
rz(-1.5016097) q[0];
sx q[0];
rz(2.0757282) q[0];
rz(-pi) q[1];
rz(2.0286331) q[2];
sx q[2];
rz(-1.0159696) q[2];
sx q[2];
rz(-1.2910281) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8431823) q[1];
sx q[1];
rz(-2.3376589) q[1];
sx q[1];
rz(1.1602379) q[1];
rz(0.92837628) q[3];
sx q[3];
rz(-1.4331665) q[3];
sx q[3];
rz(-2.9127667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8794787) q[2];
sx q[2];
rz(-0.72747362) q[2];
sx q[2];
rz(-1.0564085) q[2];
rz(-0.94111717) q[3];
sx q[3];
rz(-1.2818776) q[3];
sx q[3];
rz(-0.077805422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62517405) q[0];
sx q[0];
rz(-2.43483) q[0];
sx q[0];
rz(1.9551552) q[0];
rz(0.9390074) q[1];
sx q[1];
rz(-2.4133195) q[1];
sx q[1];
rz(0.77566385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5318308) q[0];
sx q[0];
rz(-1.8738835) q[0];
sx q[0];
rz(-2.9489452) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4294377) q[2];
sx q[2];
rz(-1.0593763) q[2];
sx q[2];
rz(-1.7101024) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6104148) q[1];
sx q[1];
rz(-1.9912331) q[1];
sx q[1];
rz(1.1320482) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9432582) q[3];
sx q[3];
rz(-2.9171067) q[3];
sx q[3];
rz(-0.62262229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0607854) q[2];
sx q[2];
rz(-0.60254097) q[2];
sx q[2];
rz(1.5134643) q[2];
rz(-1.3249409) q[3];
sx q[3];
rz(-2.5038268) q[3];
sx q[3];
rz(3.018124) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047664646) q[0];
sx q[0];
rz(-1.0075189) q[0];
sx q[0];
rz(-2.5061593) q[0];
rz(-2.0280929) q[1];
sx q[1];
rz(-0.43575382) q[1];
sx q[1];
rz(0.39438549) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8296426) q[0];
sx q[0];
rz(-1.4935635) q[0];
sx q[0];
rz(-3.0566099) q[0];
rz(-pi) q[1];
rz(-0.46994163) q[2];
sx q[2];
rz(-1.7423465) q[2];
sx q[2];
rz(-1.1217233) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.202521) q[1];
sx q[1];
rz(-0.58798169) q[1];
sx q[1];
rz(-0.67616868) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4448576) q[3];
sx q[3];
rz(-0.91154418) q[3];
sx q[3];
rz(-1.6492305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8013132) q[2];
sx q[2];
rz(-1.3542513) q[2];
sx q[2];
rz(0.24564965) q[2];
rz(2.4047325) q[3];
sx q[3];
rz(-0.95972484) q[3];
sx q[3];
rz(-2.5514742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0737632) q[0];
sx q[0];
rz(-1.2968061) q[0];
sx q[0];
rz(0.8514362) q[0];
rz(2.7379981) q[1];
sx q[1];
rz(-0.60786301) q[1];
sx q[1];
rz(3.0013989) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4969585) q[0];
sx q[0];
rz(-1.9352563) q[0];
sx q[0];
rz(-0.24043997) q[0];
rz(1.6824016) q[2];
sx q[2];
rz(-1.7652458) q[2];
sx q[2];
rz(0.60332509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.278001) q[1];
sx q[1];
rz(-0.74713445) q[1];
sx q[1];
rz(-1.1255582) q[1];
rz(-pi) q[2];
rz(1.1153489) q[3];
sx q[3];
rz(-2.1059016) q[3];
sx q[3];
rz(-0.29087152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95789528) q[2];
sx q[2];
rz(-0.87527466) q[2];
sx q[2];
rz(2.9389006) q[2];
rz(1.1787339) q[3];
sx q[3];
rz(-0.26114994) q[3];
sx q[3];
rz(2.0451666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.390585) q[0];
sx q[0];
rz(-3.106116) q[0];
sx q[0];
rz(-1.9417199) q[0];
rz(1.9203145) q[1];
sx q[1];
rz(-2.5883636) q[1];
sx q[1];
rz(-1.663307) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64740381) q[0];
sx q[0];
rz(-1.5664808) q[0];
sx q[0];
rz(-0.01135435) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9273527) q[2];
sx q[2];
rz(-2.5551559) q[2];
sx q[2];
rz(-2.9629493) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83542237) q[1];
sx q[1];
rz(-1.3886012) q[1];
sx q[1];
rz(-0.56791124) q[1];
rz(-pi) q[2];
rz(1.399089) q[3];
sx q[3];
rz(-1.2868361) q[3];
sx q[3];
rz(-0.79156706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.025909802) q[2];
sx q[2];
rz(-0.9240216) q[2];
sx q[2];
rz(-2.8795418) q[2];
rz(2.986749) q[3];
sx q[3];
rz(-0.20379977) q[3];
sx q[3];
rz(3.1018597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994951) q[0];
sx q[0];
rz(-0.69923788) q[0];
sx q[0];
rz(-0.56270885) q[0];
rz(-1.567125) q[1];
sx q[1];
rz(-1.2854311) q[1];
sx q[1];
rz(-0.80975103) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2520272) q[0];
sx q[0];
rz(-1.2148661) q[0];
sx q[0];
rz(-2.9019474) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2024683) q[2];
sx q[2];
rz(-1.0168827) q[2];
sx q[2];
rz(-0.46304747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0938125) q[1];
sx q[1];
rz(-1.8829104) q[1];
sx q[1];
rz(1.9563489) q[1];
x q[2];
rz(0.92110473) q[3];
sx q[3];
rz(-1.6739598) q[3];
sx q[3];
rz(-2.1589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8030545) q[2];
sx q[2];
rz(-0.43236098) q[2];
sx q[2];
rz(-1.7301249) q[2];
rz(2.9825397) q[3];
sx q[3];
rz(-1.8820857) q[3];
sx q[3];
rz(-0.66445309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14801046) q[0];
sx q[0];
rz(-0.2033041) q[0];
sx q[0];
rz(-2.5116442) q[0];
rz(-2.3155164) q[1];
sx q[1];
rz(-0.45336205) q[1];
sx q[1];
rz(2.1410543) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9634386) q[0];
sx q[0];
rz(-0.54433301) q[0];
sx q[0];
rz(3.0523249) q[0];
rz(2.6772772) q[2];
sx q[2];
rz(-1.3110187) q[2];
sx q[2];
rz(-0.73031727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7383768) q[1];
sx q[1];
rz(-0.94562247) q[1];
sx q[1];
rz(2.0887062) q[1];
x q[2];
rz(-2.431683) q[3];
sx q[3];
rz(-1.9485352) q[3];
sx q[3];
rz(-2.0722598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91653812) q[2];
sx q[2];
rz(-0.30320898) q[2];
sx q[2];
rz(-1.0603325) q[2];
rz(-0.66925085) q[3];
sx q[3];
rz(-2.4762912) q[3];
sx q[3];
rz(2.3745712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52436787) q[0];
sx q[0];
rz(-2.0414633) q[0];
sx q[0];
rz(2.166911) q[0];
rz(-0.928448) q[1];
sx q[1];
rz(-1.4310373) q[1];
sx q[1];
rz(1.9824082) q[1];
rz(0.84865271) q[2];
sx q[2];
rz(-2.3046222) q[2];
sx q[2];
rz(0.23474856) q[2];
rz(1.8228788) q[3];
sx q[3];
rz(-2.2889201) q[3];
sx q[3];
rz(-2.5826519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
