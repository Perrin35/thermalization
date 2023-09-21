OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3368971) q[0];
sx q[0];
rz(4.1788221) q[0];
sx q[0];
rz(9.0691789) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(-1.27966) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0227538) q[0];
sx q[0];
rz(-1.4316598) q[0];
sx q[0];
rz(1.3214146) q[0];
rz(-pi) q[1];
rz(0.42135294) q[2];
sx q[2];
rz(-1.3660649) q[2];
sx q[2];
rz(-2.8533964) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.76525926) q[1];
sx q[1];
rz(-1.1632803) q[1];
sx q[1];
rz(-2.5024662) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7198556) q[3];
sx q[3];
rz(-0.30656439) q[3];
sx q[3];
rz(2.7048064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0063643) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-2.485086) q[2];
rz(-2.4025829) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0083369) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(-3.0796675) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(-2.6541236) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9092642) q[0];
sx q[0];
rz(-3.0380913) q[0];
sx q[0];
rz(-2.2444025) q[0];
rz(-pi) q[1];
rz(-2.0080645) q[2];
sx q[2];
rz(-1.9568866) q[2];
sx q[2];
rz(-1.3646477) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4844955) q[1];
sx q[1];
rz(-0.24689455) q[1];
sx q[1];
rz(2.2730278) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8879714) q[3];
sx q[3];
rz(-1.4884236) q[3];
sx q[3];
rz(-1.2143283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(-0.3953735) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19668002) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(0.19038598) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(0.22274676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4653141) q[0];
sx q[0];
rz(-1.9440117) q[0];
sx q[0];
rz(0.42155427) q[0];
x q[1];
rz(2.9163755) q[2];
sx q[2];
rz(-1.2080964) q[2];
sx q[2];
rz(0.97937102) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.45914868) q[1];
sx q[1];
rz(-2.5249081) q[1];
sx q[1];
rz(-0.1174121) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1915216) q[3];
sx q[3];
rz(-2.7347703) q[3];
sx q[3];
rz(2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2591851) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(1.9474691) q[2];
rz(1.0549226) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(0.9202252) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.614914) q[0];
sx q[0];
rz(-0.81256142) q[0];
sx q[0];
rz(-0.99292361) q[0];
rz(-1.4843066) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(2.6038225) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3355545) q[1];
sx q[1];
rz(-1.1127377) q[1];
sx q[1];
rz(-0.6946509) q[1];
rz(2.782015) q[3];
sx q[3];
rz(-1.3782856) q[3];
sx q[3];
rz(1.0958835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(0.76830307) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(2.6884902) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26950726) q[0];
sx q[0];
rz(-3.0591024) q[0];
sx q[0];
rz(-2.0399658) q[0];
x q[1];
rz(1.4278973) q[2];
sx q[2];
rz(-1.5362527) q[2];
sx q[2];
rz(2.3847716) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0799775) q[1];
sx q[1];
rz(-1.4887759) q[1];
sx q[1];
rz(-1.7916937) q[1];
rz(-pi) q[2];
rz(-0.971332) q[3];
sx q[3];
rz(-1.1820275) q[3];
sx q[3];
rz(0.99861162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0129464) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(-1.5117234) q[2];
rz(0.72367469) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(-2.2912912) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1081651) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(-2.9034555) q[0];
rz(-2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.628081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0244004) q[0];
sx q[0];
rz(-1.2478561) q[0];
sx q[0];
rz(0.49844235) q[0];
x q[1];
rz(-2.020535) q[2];
sx q[2];
rz(-1.626484) q[2];
sx q[2];
rz(-2.4186717) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7442419) q[1];
sx q[1];
rz(-0.89351082) q[1];
sx q[1];
rz(2.1703297) q[1];
rz(-0.87168872) q[3];
sx q[3];
rz(-0.5629493) q[3];
sx q[3];
rz(-2.9339919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.05904077) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(1.9980105) q[2];
rz(0.13051662) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0156353) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(-1.827084) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(1.1869173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884739) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(0.32062809) q[0];
rz(2.6830964) q[2];
sx q[2];
rz(-3.0756604) q[2];
sx q[2];
rz(-2.3022848) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.742332) q[1];
sx q[1];
rz(-2.2982344) q[1];
sx q[1];
rz(-0.17952339) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2790473) q[3];
sx q[3];
rz(-2.1803133) q[3];
sx q[3];
rz(2.9677344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26178965) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(1.3558033) q[2];
rz(-1.6342182) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223406) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(0.85574714) q[0];
rz(3.1198655) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(-1.0303248) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5949769) q[0];
sx q[0];
rz(-0.82406509) q[0];
sx q[0];
rz(-2.0440408) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0439992) q[2];
sx q[2];
rz(-1.6128522) q[2];
sx q[2];
rz(-2.2823357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1460261) q[1];
sx q[1];
rz(-1.7221754) q[1];
sx q[1];
rz(1.8712908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1645032) q[3];
sx q[3];
rz(-1.1247375) q[3];
sx q[3];
rz(0.89531985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8490303) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940015) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(2.8253187) q[0];
rz(-1.0519741) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(-2.0064328) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59199698) q[0];
sx q[0];
rz(-1.9599008) q[0];
sx q[0];
rz(1.7745716) q[0];
rz(-pi) q[1];
rz(3.1124434) q[2];
sx q[2];
rz(-1.2505194) q[2];
sx q[2];
rz(-2.0707891) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3865349) q[1];
sx q[1];
rz(-0.48663501) q[1];
sx q[1];
rz(-1.4450032) q[1];
rz(-pi) q[2];
rz(1.5376066) q[3];
sx q[3];
rz(-2.5182708) q[3];
sx q[3];
rz(0.25191307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6241374) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(1.5448145) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.90821663) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(0.35183516) q[0];
rz(-0.31967638) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(2.9454254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0344714) q[0];
sx q[0];
rz(-0.84566085) q[0];
sx q[0];
rz(-0.79191533) q[0];
rz(-1.2162131) q[2];
sx q[2];
rz(-0.14419975) q[2];
sx q[2];
rz(-1.9229638) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6876467) q[1];
sx q[1];
rz(-1.5554264) q[1];
sx q[1];
rz(-2.3049298) q[1];
x q[2];
rz(0.84368002) q[3];
sx q[3];
rz(-1.0495249) q[3];
sx q[3];
rz(-2.664444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3907884) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(3.0604559) q[2];
rz(-1.1817415) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(-2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022973013) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(1.3148057) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(0.20028533) q[2];
sx q[2];
rz(-2.9030048) q[2];
sx q[2];
rz(1.5422274) q[2];
rz(1.3410939) q[3];
sx q[3];
rz(-1.615974) q[3];
sx q[3];
rz(1.6254197) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
