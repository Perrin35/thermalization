OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(4.0868563) q[0];
sx q[0];
rz(9.950369) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68569505) q[0];
sx q[0];
rz(-1.0766317) q[0];
sx q[0];
rz(2.6463638) q[0];
rz(-pi) q[1];
rz(-0.71360795) q[2];
sx q[2];
rz(-2.8461694) q[2];
sx q[2];
rz(1.7507391) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.875784) q[1];
sx q[1];
rz(-1.335653) q[1];
sx q[1];
rz(-1.1633412) q[1];
x q[2];
rz(0.18890394) q[3];
sx q[3];
rz(-1.8763262) q[3];
sx q[3];
rz(0.73959914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(-0.096244372) q[2];
rz(-2.105666) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44089833) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(-0.76517117) q[0];
rz(-1.2922497) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(-0.66295019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0239319) q[0];
sx q[0];
rz(-1.5045325) q[0];
sx q[0];
rz(-0.06225417) q[0];
rz(2.2234369) q[2];
sx q[2];
rz(-1.895972) q[2];
sx q[2];
rz(-2.6174389) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.029164974) q[1];
sx q[1];
rz(-1.4883092) q[1];
sx q[1];
rz(2.6998181) q[1];
x q[2];
rz(-1.7327659) q[3];
sx q[3];
rz(-2.0521945) q[3];
sx q[3];
rz(-3.0961406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.092676) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(-2.4760903) q[3];
sx q[3];
rz(-2.9232959) q[3];
sx q[3];
rz(-1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(0.22234017) q[0];
rz(-2.1242583) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(-0.51868784) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4043536) q[0];
sx q[0];
rz(-0.81439942) q[0];
sx q[0];
rz(2.0703719) q[0];
x q[1];
rz(1.475004) q[2];
sx q[2];
rz(-1.3057858) q[2];
sx q[2];
rz(1.67213) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0573404) q[1];
sx q[1];
rz(-2.3358314) q[1];
sx q[1];
rz(1.8112103) q[1];
rz(-pi) q[2];
rz(-0.014702602) q[3];
sx q[3];
rz(-3.0869752) q[3];
sx q[3];
rz(2.2811449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8609994) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-3.0730491) q[2];
rz(2.5391501) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(-0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.65072) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(-2.9673476) q[0];
rz(2.6113367) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(-0.51309103) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1888694) q[0];
sx q[0];
rz(-0.9523305) q[0];
sx q[0];
rz(-1.4737211) q[0];
x q[1];
rz(2.1999173) q[2];
sx q[2];
rz(-2.2721014) q[2];
sx q[2];
rz(-0.47052449) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0780371) q[1];
sx q[1];
rz(-2.6647898) q[1];
sx q[1];
rz(-2.8366691) q[1];
x q[2];
rz(-0.030839132) q[3];
sx q[3];
rz(-1.3645002) q[3];
sx q[3];
rz(1.7145969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47485581) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(-0.22988698) q[2];
rz(-0.41904467) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(-2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6722365) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(0.67681926) q[0];
rz(-0.49304402) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(-2.5255323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6989243) q[0];
sx q[0];
rz(-1.6232423) q[0];
sx q[0];
rz(-0.035874493) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30454238) q[2];
sx q[2];
rz(-0.38123044) q[2];
sx q[2];
rz(2.9074557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9628145) q[1];
sx q[1];
rz(-1.5652834) q[1];
sx q[1];
rz(0.57761901) q[1];
x q[2];
rz(-2.8354007) q[3];
sx q[3];
rz(-2.2243143) q[3];
sx q[3];
rz(2.9432076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4074576) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(-2.8862254) q[2];
rz(-1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(-0.47250026) q[0];
rz(0.52608144) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(-2.2568259) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038267604) q[0];
sx q[0];
rz(-1.190435) q[0];
sx q[0];
rz(0.079770712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.066263513) q[2];
sx q[2];
rz(-2.931086) q[2];
sx q[2];
rz(-1.9312242) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6217475) q[1];
sx q[1];
rz(-1.2832844) q[1];
sx q[1];
rz(-1.9985755) q[1];
rz(-2.8956036) q[3];
sx q[3];
rz(-2.0912598) q[3];
sx q[3];
rz(0.99508475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3970967) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(-2.8302144) q[2];
rz(1.7729676) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(-0.5293203) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41806528) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(-0.73289245) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6717364) q[0];
sx q[0];
rz(-2.0977019) q[0];
sx q[0];
rz(1.1080452) q[0];
rz(-pi) q[1];
rz(-0.47809017) q[2];
sx q[2];
rz(-2.2149137) q[2];
sx q[2];
rz(0.23020506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0098795) q[1];
sx q[1];
rz(-0.94505802) q[1];
sx q[1];
rz(3.1335319) q[1];
rz(-pi) q[2];
rz(1.4634499) q[3];
sx q[3];
rz(-1.6119453) q[3];
sx q[3];
rz(0.51957182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(2.627009) q[2];
rz(-1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(-0.54491836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056203689) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(0.7094267) q[0];
rz(-1.6363232) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(2.8628796) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3942791) q[0];
sx q[0];
rz(-2.8907667) q[0];
sx q[0];
rz(-2.5670811) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8144572) q[2];
sx q[2];
rz(-2.116034) q[2];
sx q[2];
rz(-1.3512163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.41234327) q[1];
sx q[1];
rz(-1.7885498) q[1];
sx q[1];
rz(-1.2530243) q[1];
x q[2];
rz(-1.2392427) q[3];
sx q[3];
rz(-2.9058876) q[3];
sx q[3];
rz(1.2329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.30148208) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(-3.0855132) q[2];
rz(0.85514832) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(-0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1466325) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(-0.96889281) q[0];
rz(0.47337198) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(-1.1425346) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24507095) q[0];
sx q[0];
rz(-2.8327496) q[0];
sx q[0];
rz(1.8737428) q[0];
rz(-pi) q[1];
rz(-3.1387024) q[2];
sx q[2];
rz(-2.4382466) q[2];
sx q[2];
rz(3.0386915) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0513623) q[1];
sx q[1];
rz(-1.4453332) q[1];
sx q[1];
rz(-2.5999703) q[1];
rz(-0.90060602) q[3];
sx q[3];
rz(-1.4108676) q[3];
sx q[3];
rz(0.98774324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(-2.1208105) q[2];
rz(2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.97994119) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(-2.4401869) q[0];
rz(-0.91570634) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(1.9030301) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4558444) q[0];
sx q[0];
rz(-2.0861174) q[0];
sx q[0];
rz(2.8932163) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6221223) q[2];
sx q[2];
rz(-0.81846279) q[2];
sx q[2];
rz(-0.78879702) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4148256) q[1];
sx q[1];
rz(-1.6884202) q[1];
sx q[1];
rz(-2.2284501) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.978312) q[3];
sx q[3];
rz(-1.6654135) q[3];
sx q[3];
rz(-2.832151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4548268) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(-1.2010126) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(-1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6702406) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(0.025370601) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(1.17169) q[2];
sx q[2];
rz(-1.8532955) q[2];
sx q[2];
rz(2.4886139) q[2];
rz(2.3425441) q[3];
sx q[3];
rz(-1.4481164) q[3];
sx q[3];
rz(-1.8275402) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];