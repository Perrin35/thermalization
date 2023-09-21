OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73206168) q[0];
sx q[0];
rz(4.5067956) q[0];
sx q[0];
rz(11.542008) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(-1.9722809) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60351935) q[0];
sx q[0];
rz(-2.3581714) q[0];
sx q[0];
rz(2.6253683) q[0];
rz(-2.2157482) q[2];
sx q[2];
rz(-2.1241509) q[2];
sx q[2];
rz(0.64253053) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6572666) q[1];
sx q[1];
rz(-2.6063759) q[1];
sx q[1];
rz(-2.3372997) q[1];
rz(0.3124247) q[3];
sx q[3];
rz(-1.4989304) q[3];
sx q[3];
rz(1.7463328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(1.8189836) q[2];
rz(2.8406075) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(-1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(-0.47505501) q[0];
rz(1.7430199) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(2.1038726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71887165) q[0];
sx q[0];
rz(-1.0178716) q[0];
sx q[0];
rz(1.3119112) q[0];
rz(-pi) q[1];
rz(-0.77387626) q[2];
sx q[2];
rz(-2.4485588) q[2];
sx q[2];
rz(-0.80640031) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4251551) q[1];
sx q[1];
rz(-1.5967224) q[1];
sx q[1];
rz(-1.1643216) q[1];
rz(2.188835) q[3];
sx q[3];
rz(-2.2817094) q[3];
sx q[3];
rz(3.0666805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(-2.7627913) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(-2.7509007) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(0.57317615) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5432376) q[0];
sx q[0];
rz(-0.43524536) q[0];
sx q[0];
rz(1.4198562) q[0];
x q[1];
rz(2.3689752) q[2];
sx q[2];
rz(-1.3396016) q[2];
sx q[2];
rz(-0.85341838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0108311) q[1];
sx q[1];
rz(-1.2847932) q[1];
sx q[1];
rz(-2.7421013) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40804789) q[3];
sx q[3];
rz(-1.1432262) q[3];
sx q[3];
rz(1.6392631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(2.9476681) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(-0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28213421) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(2.0129054) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(2.7788924) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6837316) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(1.213221) q[0];
rz(-pi) q[1];
rz(-1.0417468) q[2];
sx q[2];
rz(-0.036465557) q[2];
sx q[2];
rz(-1.0066102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.73211654) q[1];
sx q[1];
rz(-0.64142694) q[1];
sx q[1];
rz(2.9684121) q[1];
rz(-1.4705212) q[3];
sx q[3];
rz(-2.5146211) q[3];
sx q[3];
rz(0.37563045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(3.1385699) q[2];
rz(-0.65888843) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9534849) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(3.127393) q[0];
rz(-0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(1.682122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1634954) q[0];
sx q[0];
rz(-2.8601544) q[0];
sx q[0];
rz(2.1241758) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2541788) q[2];
sx q[2];
rz(-1.1537342) q[2];
sx q[2];
rz(2.8487157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0121213) q[1];
sx q[1];
rz(-2.8399889) q[1];
sx q[1];
rz(2.4232037) q[1];
x q[2];
rz(-2.1305389) q[3];
sx q[3];
rz(-0.55941814) q[3];
sx q[3];
rz(1.884348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(-1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(-2.5573964) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(-0.13866436) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9753871) q[0];
sx q[0];
rz(-1.5741325) q[0];
sx q[0];
rz(3.1211981) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9406592) q[2];
sx q[2];
rz(-1.4954508) q[2];
sx q[2];
rz(-2.2546525) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0181959) q[1];
sx q[1];
rz(-2.1235848) q[1];
sx q[1];
rz(0.68560302) q[1];
rz(0.45913978) q[3];
sx q[3];
rz(-2.3448179) q[3];
sx q[3];
rz(1.2129984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77506322) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(1.8072051) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(-0.095120393) q[3];
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
rz(-0.60550624) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(2.6050674) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(0.62430635) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1034575) q[0];
sx q[0];
rz(-0.75224829) q[0];
sx q[0];
rz(-0.30970807) q[0];
rz(-0.46632669) q[2];
sx q[2];
rz(-2.1412686) q[2];
sx q[2];
rz(1.2517267) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1075322) q[1];
sx q[1];
rz(-2.1851087) q[1];
sx q[1];
rz(0.88453102) q[1];
x q[2];
rz(2.1950486) q[3];
sx q[3];
rz(-0.72580273) q[3];
sx q[3];
rz(-0.040369999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5252934) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(-3.110102) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87576762) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(-0.12891842) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3710204) q[0];
sx q[0];
rz(-0.41045529) q[0];
sx q[0];
rz(-1.8380941) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89739563) q[2];
sx q[2];
rz(-0.79332966) q[2];
sx q[2];
rz(-2.430254) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0205295) q[1];
sx q[1];
rz(-1.1464835) q[1];
sx q[1];
rz(-1.2709649) q[1];
rz(-pi) q[2];
rz(1.586732) q[3];
sx q[3];
rz(-2.3105379) q[3];
sx q[3];
rz(2.9908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63697469) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(-0.62409419) q[2];
rz(0.23877731) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794466) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.8918442) q[0];
rz(0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.3508266) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6791145) q[0];
sx q[0];
rz(-2.9318641) q[0];
sx q[0];
rz(2.1049343) q[0];
x q[1];
rz(2.5596041) q[2];
sx q[2];
rz(-0.66170035) q[2];
sx q[2];
rz(0.87994196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0235325) q[1];
sx q[1];
rz(-2.2242821) q[1];
sx q[1];
rz(0.94783028) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75307122) q[3];
sx q[3];
rz(-1.3967447) q[3];
sx q[3];
rz(-0.93833246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0525557) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(0.11432153) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91530144) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(2.7217216) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(0.54668033) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71781681) q[0];
sx q[0];
rz(-0.86200889) q[0];
sx q[0];
rz(2.8481759) q[0];
x q[1];
rz(-1.7231862) q[2];
sx q[2];
rz(-1.1001462) q[2];
sx q[2];
rz(-1.6413123) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9676799) q[1];
sx q[1];
rz(-0.2914857) q[1];
sx q[1];
rz(-3.013054) q[1];
x q[2];
rz(-0.22565266) q[3];
sx q[3];
rz(-1.435558) q[3];
sx q[3];
rz(-1.7547363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13835779) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(0.004301087) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99988408) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(-1.8757204) q[2];
sx q[2];
rz(-0.049449895) q[2];
sx q[2];
rz(0.54686875) q[2];
rz(0.32904939) q[3];
sx q[3];
rz(-1.1209189) q[3];
sx q[3];
rz(0.95595595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
