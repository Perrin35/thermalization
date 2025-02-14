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
rz(-1.2404233) q[0];
sx q[0];
rz(-1.197553) q[0];
sx q[0];
rz(2.9252606) q[0];
rz(0.95739111) q[1];
sx q[1];
rz(-2.4654145) q[1];
sx q[1];
rz(1.5114991) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0441598) q[0];
sx q[0];
rz(-2.5852647) q[0];
sx q[0];
rz(-3.1233643) q[0];
x q[1];
rz(2.6681603) q[2];
sx q[2];
rz(-1.8634999) q[2];
sx q[2];
rz(-2.6930489) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2547467) q[1];
sx q[1];
rz(-1.0751845) q[1];
sx q[1];
rz(1.9504471) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7483398) q[3];
sx q[3];
rz(-0.90688469) q[3];
sx q[3];
rz(0.25379405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.034885255) q[2];
sx q[2];
rz(-0.35197508) q[2];
sx q[2];
rz(-1.0484288) q[2];
rz(-0.18167051) q[3];
sx q[3];
rz(-0.964966) q[3];
sx q[3];
rz(-0.65339965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6492017) q[0];
sx q[0];
rz(-2.1972456) q[0];
sx q[0];
rz(-0.44678584) q[0];
rz(0.88042879) q[1];
sx q[1];
rz(-1.7767521) q[1];
sx q[1];
rz(2.3588038) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7715817) q[0];
sx q[0];
rz(-0.25213045) q[0];
sx q[0];
rz(1.8186411) q[0];
rz(-pi) q[1];
rz(0.61666416) q[2];
sx q[2];
rz(-1.9266124) q[2];
sx q[2];
rz(2.7105376) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35905169) q[1];
sx q[1];
rz(-0.56582574) q[1];
sx q[1];
rz(-0.58360175) q[1];
x q[2];
rz(1.0083593) q[3];
sx q[3];
rz(-0.76518067) q[3];
sx q[3];
rz(-2.6443114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.030152628) q[2];
sx q[2];
rz(-1.6609265) q[2];
sx q[2];
rz(2.1622369) q[2];
rz(2.7566946) q[3];
sx q[3];
rz(-1.9210457) q[3];
sx q[3];
rz(2.9773007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48591831) q[0];
sx q[0];
rz(-0.05412183) q[0];
sx q[0];
rz(2.3517877) q[0];
rz(0.18665953) q[1];
sx q[1];
rz(-1.4013545) q[1];
sx q[1];
rz(0.99367118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.073245) q[0];
sx q[0];
rz(-0.6093502) q[0];
sx q[0];
rz(0.44550271) q[0];
rz(-2.5285401) q[2];
sx q[2];
rz(-1.384735) q[2];
sx q[2];
rz(1.2836054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6216065) q[1];
sx q[1];
rz(-0.64711678) q[1];
sx q[1];
rz(-0.81955975) q[1];
x q[2];
rz(0.72088269) q[3];
sx q[3];
rz(-1.3647121) q[3];
sx q[3];
rz(-2.4997847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8197202) q[2];
sx q[2];
rz(-2.0237782) q[2];
sx q[2];
rz(-0.37128386) q[2];
rz(-0.34879455) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(2.546052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.45469859) q[0];
sx q[0];
rz(-2.1214387) q[0];
sx q[0];
rz(-1.4042847) q[0];
rz(0.67717254) q[1];
sx q[1];
rz(-1.9858457) q[1];
sx q[1];
rz(1.6489395) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66144511) q[0];
sx q[0];
rz(-2.8146525) q[0];
sx q[0];
rz(-0.72823712) q[0];
rz(1.1392966) q[2];
sx q[2];
rz(-1.8427263) q[2];
sx q[2];
rz(1.4828585) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3824635) q[1];
sx q[1];
rz(-1.8378403) q[1];
sx q[1];
rz(-1.8014531) q[1];
rz(-pi) q[2];
rz(-1.9203606) q[3];
sx q[3];
rz(-2.3858983) q[3];
sx q[3];
rz(-0.44103482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45381418) q[2];
sx q[2];
rz(-0.9684338) q[2];
sx q[2];
rz(2.7933534) q[2];
rz(-1.4549152) q[3];
sx q[3];
rz(-1.429052) q[3];
sx q[3];
rz(-1.3454364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9652047) q[0];
sx q[0];
rz(-0.56469733) q[0];
sx q[0];
rz(-0.89214605) q[0];
rz(-0.46420321) q[1];
sx q[1];
rz(-1.2449539) q[1];
sx q[1];
rz(-1.8468599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1071467) q[0];
sx q[0];
rz(-0.71104103) q[0];
sx q[0];
rz(2.3734295) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.928435) q[2];
sx q[2];
rz(-2.5937754) q[2];
sx q[2];
rz(2.167706) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6695822) q[1];
sx q[1];
rz(-1.246844) q[1];
sx q[1];
rz(-2.0355909) q[1];
rz(-pi) q[2];
rz(2.5312349) q[3];
sx q[3];
rz(-2.1564283) q[3];
sx q[3];
rz(-0.40025362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8190454) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(2.5757705) q[2];
rz(-3.0602509) q[3];
sx q[3];
rz(-2.1606725) q[3];
sx q[3];
rz(-0.62059039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6767947) q[0];
sx q[0];
rz(-0.066635266) q[0];
sx q[0];
rz(1.5860522) q[0];
rz(1.0606891) q[1];
sx q[1];
rz(-1.565275) q[1];
sx q[1];
rz(2.5097844) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87651265) q[0];
sx q[0];
rz(-1.8051462) q[0];
sx q[0];
rz(1.6523907) q[0];
rz(-pi) q[1];
rz(-0.49922322) q[2];
sx q[2];
rz(-2.2167335) q[2];
sx q[2];
rz(1.2690085) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7673847) q[1];
sx q[1];
rz(-0.65113089) q[1];
sx q[1];
rz(-0.99094772) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68354179) q[3];
sx q[3];
rz(-1.6850796) q[3];
sx q[3];
rz(-1.6268693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98171988) q[2];
sx q[2];
rz(-1.1184511) q[2];
sx q[2];
rz(-2.6427606) q[2];
rz(-1.3129129) q[3];
sx q[3];
rz(-2.4098318) q[3];
sx q[3];
rz(1.4341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.6071534) q[0];
sx q[0];
rz(-0.3648912) q[0];
sx q[0];
rz(0.69951192) q[0];
rz(2.6761159) q[1];
sx q[1];
rz(-2.270348) q[1];
sx q[1];
rz(2.2043998) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4995183) q[0];
sx q[0];
rz(-2.8694911) q[0];
sx q[0];
rz(2.4689552) q[0];
x q[1];
rz(-2.7363051) q[2];
sx q[2];
rz(-0.25988419) q[2];
sx q[2];
rz(2.0854307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22058567) q[1];
sx q[1];
rz(-1.3594419) q[1];
sx q[1];
rz(1.8615972) q[1];
x q[2];
rz(-1.8645511) q[3];
sx q[3];
rz(-1.2806727) q[3];
sx q[3];
rz(1.6725095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.84404868) q[2];
sx q[2];
rz(-1.301845) q[2];
sx q[2];
rz(-2.76827) q[2];
rz(-1.9518055) q[3];
sx q[3];
rz(-0.50896421) q[3];
sx q[3];
rz(0.51026195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74476403) q[0];
sx q[0];
rz(-2.0592392) q[0];
sx q[0];
rz(-2.3936791) q[0];
rz(2.3751936) q[1];
sx q[1];
rz(-0.26873573) q[1];
sx q[1];
rz(0.0029729923) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.620913) q[0];
sx q[0];
rz(-1.7629884) q[0];
sx q[0];
rz(0.8780794) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.019549088) q[2];
sx q[2];
rz(-1.4998933) q[2];
sx q[2];
rz(-0.30724684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.06719477) q[1];
sx q[1];
rz(-1.7247827) q[1];
sx q[1];
rz(3.0976553) q[1];
rz(-pi) q[2];
rz(1.7693172) q[3];
sx q[3];
rz(-1.7619507) q[3];
sx q[3];
rz(-1.6894345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11091867) q[2];
sx q[2];
rz(-1.3838394) q[2];
sx q[2];
rz(2.0745011) q[2];
rz(-3.057632) q[3];
sx q[3];
rz(-0.48560086) q[3];
sx q[3];
rz(-0.77897227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79754168) q[0];
sx q[0];
rz(-0.9386971) q[0];
sx q[0];
rz(-0.10636605) q[0];
rz(2.1616409) q[1];
sx q[1];
rz(-1.6500902) q[1];
sx q[1];
rz(-0.76593691) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87493103) q[0];
sx q[0];
rz(-0.70737544) q[0];
sx q[0];
rz(-0.69702638) q[0];
rz(-pi) q[1];
rz(1.4808876) q[2];
sx q[2];
rz(-1.3718296) q[2];
sx q[2];
rz(2.112191) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3801743) q[1];
sx q[1];
rz(-1.6055709) q[1];
sx q[1];
rz(-2.4654287) q[1];
x q[2];
rz(-2.8838653) q[3];
sx q[3];
rz(-3.1270087) q[3];
sx q[3];
rz(-2.7590883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47436675) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(-2.4412947) q[2];
rz(0.66655603) q[3];
sx q[3];
rz(-1.8066112) q[3];
sx q[3];
rz(-2.0535645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4661082) q[0];
sx q[0];
rz(-1.0270783) q[0];
sx q[0];
rz(1.745537) q[0];
rz(-1.667977) q[1];
sx q[1];
rz(-1.2584078) q[1];
sx q[1];
rz(2.4748763) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6121126) q[0];
sx q[0];
rz(-2.1938938) q[0];
sx q[0];
rz(0.77147958) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6595419) q[2];
sx q[2];
rz(-1.7983984) q[2];
sx q[2];
rz(1.2510117) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8620876) q[1];
sx q[1];
rz(-2.7976329) q[1];
sx q[1];
rz(-2.7441447) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3844149) q[3];
sx q[3];
rz(-0.72523967) q[3];
sx q[3];
rz(0.66981572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0746158) q[2];
sx q[2];
rz(-0.65698996) q[2];
sx q[2];
rz(-2.2508049) q[2];
rz(2.7095419) q[3];
sx q[3];
rz(-2.0712349) q[3];
sx q[3];
rz(-0.74705684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73410949) q[0];
sx q[0];
rz(-1.643184) q[0];
sx q[0];
rz(-1.2947422) q[0];
rz(-1.2474077) q[1];
sx q[1];
rz(-2.4220962) q[1];
sx q[1];
rz(-1.9069506) q[1];
rz(-2.9802889) q[2];
sx q[2];
rz(-1.9111173) q[2];
sx q[2];
rz(2.5572122) q[2];
rz(-2.6525146) q[3];
sx q[3];
rz(-1.5413956) q[3];
sx q[3];
rz(1.2914381) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
