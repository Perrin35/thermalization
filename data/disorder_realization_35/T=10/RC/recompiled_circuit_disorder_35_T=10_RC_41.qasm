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
rz(-1.7763897) q[0];
sx q[0];
rz(2.1172297) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(-1.9722809) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60351935) q[0];
sx q[0];
rz(-0.78342122) q[0];
sx q[0];
rz(-2.6253683) q[0];
rz(-pi) q[1];
rz(-0.65806234) q[2];
sx q[2];
rz(-1.0339289) q[2];
sx q[2];
rz(1.8368349) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.642627) q[1];
sx q[1];
rz(-1.9470012) q[1];
sx q[1];
rz(0.39019231) q[1];
rz(0.23006769) q[3];
sx q[3];
rz(-2.821273) q[3];
sx q[3];
rz(2.7473118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8756276) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(-1.3226091) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(-1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(0.47505501) q[0];
rz(-1.3985727) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(-1.0377201) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25181928) q[0];
sx q[0];
rz(-2.5368241) q[0];
sx q[0];
rz(-2.7483727) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0449045) q[2];
sx q[2];
rz(-2.0453339) q[2];
sx q[2];
rz(1.431312) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2270826) q[1];
sx q[1];
rz(-0.40725476) q[1];
sx q[1];
rz(1.6362908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5492937) q[3];
sx q[3];
rz(-0.90511887) q[3];
sx q[3];
rz(-2.388282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(0.37880138) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(-1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5304853) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(0.57317615) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.976982) q[0];
sx q[0];
rz(-1.5073538) q[0];
sx q[0];
rz(-2.001686) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3689752) q[2];
sx q[2];
rz(-1.3396016) q[2];
sx q[2];
rz(-2.2881743) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0108311) q[1];
sx q[1];
rz(-1.2847932) q[1];
sx q[1];
rz(-0.39949135) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8543386) q[3];
sx q[3];
rz(-2.5594098) q[3];
sx q[3];
rz(2.4454988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0009784) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(0.097269416) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(-0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8594584) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(1.1286873) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(0.36270025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1283778) q[0];
sx q[0];
rz(-0.88071874) q[0];
sx q[0];
rz(2.8280558) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1231819) q[2];
sx q[2];
rz(-1.6022748) q[2];
sx q[2];
rz(1.6056431) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6994233) q[1];
sx q[1];
rz(-1.6740834) q[1];
sx q[1];
rz(0.63421722) q[1];
rz(-1.6710715) q[3];
sx q[3];
rz(-2.5146211) q[3];
sx q[3];
rz(2.7659622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68391934) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(-2.4827042) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(2.3390884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(3.127393) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(-1.682122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97809726) q[0];
sx q[0];
rz(-0.28143829) q[0];
sx q[0];
rz(2.1241758) q[0];
rz(-pi) q[1];
rz(0.51909165) q[2];
sx q[2];
rz(-0.95538288) q[2];
sx q[2];
rz(1.5965243) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38813218) q[1];
sx q[1];
rz(-1.7963444) q[1];
sx q[1];
rz(-1.3688341) q[1];
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
rz(-pi/2) q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(2.2996976) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1275948) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(-1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(3.0029283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9753871) q[0];
sx q[0];
rz(-1.5674601) q[0];
sx q[0];
rz(0.020394527) q[0];
x q[1];
rz(0.20093341) q[2];
sx q[2];
rz(-1.6461419) q[2];
sx q[2];
rz(-0.88694015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0181959) q[1];
sx q[1];
rz(-1.0180078) q[1];
sx q[1];
rz(-2.4559896) q[1];
rz(-pi) q[2];
rz(1.1451374) q[3];
sx q[3];
rz(-2.2666551) q[3];
sx q[3];
rz(1.828572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(1.8072051) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(-0.095120393) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(-0.62430635) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37492232) q[0];
sx q[0];
rz(-0.86219388) q[0];
sx q[0];
rz(1.8486345) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95958556) q[2];
sx q[2];
rz(-2.4215536) q[2];
sx q[2];
rz(2.6401273) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.034060409) q[1];
sx q[1];
rz(-0.95648396) q[1];
sx q[1];
rz(0.88453102) q[1];
rz(-pi) q[2];
rz(-0.94675605) q[3];
sx q[3];
rz(-1.1723926) q[3];
sx q[3];
rz(1.0362253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-2.7590511) q[2];
rz(0.031490695) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.265825) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(-1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(0.12891842) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55420586) q[0];
sx q[0];
rz(-1.4652068) q[0];
sx q[0];
rz(1.9681853) q[0];
rz(-2.5768186) q[2];
sx q[2];
rz(-2.161705) q[2];
sx q[2];
rz(1.5806944) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0205295) q[1];
sx q[1];
rz(-1.9951092) q[1];
sx q[1];
rz(-1.2709649) q[1];
rz(-pi) q[2];
rz(2.4017879) q[3];
sx q[3];
rz(-1.5590258) q[3];
sx q[3];
rz(1.7108325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(2.5174985) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(-1.0673267) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794466) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.8918442) q[0];
rz(-3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.3508266) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352167) q[0];
sx q[0];
rz(-1.750964) q[0];
sx q[0];
rz(-3.0336477) q[0];
x q[1];
rz(-2.5647854) q[2];
sx q[2];
rz(-1.2262605) q[2];
sx q[2];
rz(0.21201269) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13547922) q[1];
sx q[1];
rz(-2.0524426) q[1];
sx q[1];
rz(-0.7559795) q[1];
rz(-0.75307122) q[3];
sx q[3];
rz(-1.3967447) q[3];
sx q[3];
rz(0.93833246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0525557) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(-0.11432153) q[2];
rz(1.8814686) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(-1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(-0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(2.5949123) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71781681) q[0];
sx q[0];
rz(-2.2795838) q[0];
sx q[0];
rz(-2.8481759) q[0];
x q[1];
rz(1.7231862) q[2];
sx q[2];
rz(-1.1001462) q[2];
sx q[2];
rz(1.6413123) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9676799) q[1];
sx q[1];
rz(-0.2914857) q[1];
sx q[1];
rz(0.12853865) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22565266) q[3];
sx q[3];
rz(-1.435558) q[3];
sx q[3];
rz(1.7547363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(-0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1417086) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(1.5236241) q[2];
sx q[2];
rz(-1.5559559) q[2];
sx q[2];
rz(2.422239) q[2];
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
