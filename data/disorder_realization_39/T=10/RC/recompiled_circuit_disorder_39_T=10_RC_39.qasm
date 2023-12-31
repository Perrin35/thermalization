OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(-0.46407035) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(6.2160677) q[1];
sx q[1];
rz(10.709229) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0425134) q[0];
sx q[0];
rz(-0.49603841) q[0];
sx q[0];
rz(2.2201559) q[0];
rz(-pi) q[1];
rz(1.60381) q[2];
sx q[2];
rz(-1.0812949) q[2];
sx q[2];
rz(0.93544338) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1661108) q[1];
sx q[1];
rz(-1.9327285) q[1];
sx q[1];
rz(0.48717498) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5873317) q[3];
sx q[3];
rz(-1.4201122) q[3];
sx q[3];
rz(-2.2076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39711943) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(-0.31952566) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-0.47839034) q[3];
sx q[3];
rz(-0.67392504) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7330866) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(2.615036) q[0];
rz(2.5780442) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(-2.3449576) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8684727) q[0];
sx q[0];
rz(-0.94808775) q[0];
sx q[0];
rz(-2.7362583) q[0];
rz(0.57467069) q[2];
sx q[2];
rz(-1.7386912) q[2];
sx q[2];
rz(1.6811973) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7044201) q[1];
sx q[1];
rz(-2.17802) q[1];
sx q[1];
rz(-1.5864658) q[1];
x q[2];
rz(0.96062406) q[3];
sx q[3];
rz(-1.1517797) q[3];
sx q[3];
rz(-2.744439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18156302) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(-0.78655085) q[2];
rz(2.6484047) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(-2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(1.7012117) q[0];
rz(0.72021833) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(0.70297855) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9706791) q[0];
sx q[0];
rz(-1.850607) q[0];
sx q[0];
rz(3.0199416) q[0];
rz(2.643232) q[2];
sx q[2];
rz(-1.8866072) q[2];
sx q[2];
rz(-0.74893803) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30333334) q[1];
sx q[1];
rz(-1.5240182) q[1];
sx q[1];
rz(-0.10951885) q[1];
rz(-pi) q[2];
rz(-2.1941575) q[3];
sx q[3];
rz(-1.7121592) q[3];
sx q[3];
rz(2.3879104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8748223) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754958) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(1.874118) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(-0.56328303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.253809) q[0];
sx q[0];
rz(-0.94622181) q[0];
sx q[0];
rz(-0.33872351) q[0];
rz(0.40043719) q[2];
sx q[2];
rz(-2.4106328) q[2];
sx q[2];
rz(1.2618582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1174417) q[1];
sx q[1];
rz(-1.5772595) q[1];
sx q[1];
rz(-0.27177377) q[1];
x q[2];
rz(-2.4241583) q[3];
sx q[3];
rz(-1.3460025) q[3];
sx q[3];
rz(-0.7625398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(-1.9794827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6974555) q[0];
sx q[0];
rz(-1.1778957) q[0];
sx q[0];
rz(-2.5494954) q[0];
rz(-pi) q[1];
rz(1.048462) q[2];
sx q[2];
rz(-1.221721) q[2];
sx q[2];
rz(-2.366684) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1197549) q[1];
sx q[1];
rz(-0.86958414) q[1];
sx q[1];
rz(-0.55418684) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.01708548) q[3];
sx q[3];
rz(-2.4199329) q[3];
sx q[3];
rz(-2.7929896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6953485) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(-1.2472786) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(-2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(-0.75063467) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(-1.8355339) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36201492) q[0];
sx q[0];
rz(-1.0581731) q[0];
sx q[0];
rz(-0.4431475) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46361228) q[2];
sx q[2];
rz(-1.7525275) q[2];
sx q[2];
rz(-2.6169427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2893147) q[1];
sx q[1];
rz(-1.5631952) q[1];
sx q[1];
rz(-0.72035933) q[1];
rz(1.7197051) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(-2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3874454) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3951185) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(-2.0948998) q[0];
rz(-1.612161) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(2.7244862) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55506599) q[0];
sx q[0];
rz(-1.5952206) q[0];
sx q[0];
rz(-1.3506372) q[0];
x q[1];
rz(2.3796758) q[2];
sx q[2];
rz(-0.41315213) q[2];
sx q[2];
rz(-0.44943902) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4579791) q[1];
sx q[1];
rz(-2.2375467) q[1];
sx q[1];
rz(-0.97138202) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7583371) q[3];
sx q[3];
rz(-2.4517422) q[3];
sx q[3];
rz(0.81724973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1356915) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(1.3674412) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(1.7143543) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(-3.0292125) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0880786) q[0];
sx q[0];
rz(-1.7435762) q[0];
sx q[0];
rz(1.0900351) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8939395) q[2];
sx q[2];
rz(-0.34341771) q[2];
sx q[2];
rz(-2.0695956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5210133) q[1];
sx q[1];
rz(-2.3304906) q[1];
sx q[1];
rz(1.0902507) q[1];
rz(3.0197969) q[3];
sx q[3];
rz(-1.313901) q[3];
sx q[3];
rz(-1.6459873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(1.3486264) q[2];
rz(1.2049234) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-0.034974139) q[0];
rz(0.84683013) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(-2.2299178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23110403) q[0];
sx q[0];
rz(-1.3470874) q[0];
sx q[0];
rz(-1.847812) q[0];
x q[1];
rz(-2.9334925) q[2];
sx q[2];
rz(-2.5884183) q[2];
sx q[2];
rz(-2.2262239) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.94757838) q[1];
sx q[1];
rz(-1.7742426) q[1];
sx q[1];
rz(-2.6737763) q[1];
rz(-pi) q[2];
rz(-2.0613641) q[3];
sx q[3];
rz(-2.053223) q[3];
sx q[3];
rz(0.8182984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(2.8923477) q[2];
rz(-0.76672673) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578167) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(2.7695079) q[0];
rz(-2.5601939) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(-1.7262329) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85522643) q[0];
sx q[0];
rz(-0.33272538) q[0];
sx q[0];
rz(2.4871728) q[0];
x q[1];
rz(0.66843372) q[2];
sx q[2];
rz(-2.4032776) q[2];
sx q[2];
rz(0.43503209) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6485939) q[1];
sx q[1];
rz(-3.0393638) q[1];
sx q[1];
rz(-0.74536721) q[1];
rz(-0.38378999) q[3];
sx q[3];
rz(-1.0597611) q[3];
sx q[3];
rz(-2.8904861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8156478) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(-2.9369205) q[2];
rz(-1.4137911) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50080147) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(-1.5851371) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(-1.1164222) q[2];
sx q[2];
rz(-1.4618256) q[2];
sx q[2];
rz(1.7088919) q[2];
rz(-1.4428044) q[3];
sx q[3];
rz(-0.57590579) q[3];
sx q[3];
rz(-0.87344575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
