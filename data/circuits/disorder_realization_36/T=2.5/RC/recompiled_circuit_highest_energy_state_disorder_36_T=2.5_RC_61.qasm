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
rz(3.1004768) q[0];
sx q[0];
rz(-1.2007204) q[0];
sx q[0];
rz(-1.6706985) q[0];
rz(-0.38493758) q[1];
sx q[1];
rz(-0.5783143) q[1];
sx q[1];
rz(2.2338423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2157057) q[0];
sx q[0];
rz(-2.3952847) q[0];
sx q[0];
rz(0.67912905) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8706246) q[2];
sx q[2];
rz(-1.289164) q[2];
sx q[2];
rz(-2.9099438) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5564087) q[1];
sx q[1];
rz(-2.6165462) q[1];
sx q[1];
rz(-0.37918143) q[1];
rz(-pi) q[2];
rz(2.7192017) q[3];
sx q[3];
rz(-3.0085519) q[3];
sx q[3];
rz(0.62158004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3140728) q[2];
sx q[2];
rz(-1.4853518) q[2];
sx q[2];
rz(-0.32076389) q[2];
rz(-0.00094207923) q[3];
sx q[3];
rz(-0.92206803) q[3];
sx q[3];
rz(-0.48085406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28474057) q[0];
sx q[0];
rz(-0.72035235) q[0];
sx q[0];
rz(-2.5208933) q[0];
rz(1.2868098) q[1];
sx q[1];
rz(-1.4871253) q[1];
sx q[1];
rz(-2.1466045) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8061229) q[0];
sx q[0];
rz(-2.5988082) q[0];
sx q[0];
rz(-1.1837028) q[0];
rz(-2.8212655) q[2];
sx q[2];
rz(-1.1503845) q[2];
sx q[2];
rz(2.0589895) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7409119) q[1];
sx q[1];
rz(-1.3865292) q[1];
sx q[1];
rz(-1.1180941) q[1];
rz(-1.0424588) q[3];
sx q[3];
rz(-2.3064559) q[3];
sx q[3];
rz(-0.56921116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2113125) q[2];
sx q[2];
rz(-1.2418396) q[2];
sx q[2];
rz(2.9105183) q[2];
rz(1.5288345) q[3];
sx q[3];
rz(-1.4569747) q[3];
sx q[3];
rz(0.55135623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424161) q[0];
sx q[0];
rz(-1.1671952) q[0];
sx q[0];
rz(2.7447847) q[0];
rz(1.1948168) q[1];
sx q[1];
rz(-1.2620986) q[1];
sx q[1];
rz(1.4746812) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75735649) q[0];
sx q[0];
rz(-0.088619329) q[0];
sx q[0];
rz(-1.4411974) q[0];
rz(-1.0753158) q[2];
sx q[2];
rz(-1.0766509) q[2];
sx q[2];
rz(-1.0841498) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35016216) q[1];
sx q[1];
rz(-1.9099978) q[1];
sx q[1];
rz(-0.8703402) q[1];
rz(-pi) q[2];
rz(3.1297333) q[3];
sx q[3];
rz(-1.1766772) q[3];
sx q[3];
rz(0.32255641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84612238) q[2];
sx q[2];
rz(-2.2961871) q[2];
sx q[2];
rz(1.4836614) q[2];
rz(-2.1814003) q[3];
sx q[3];
rz(-0.63513297) q[3];
sx q[3];
rz(-0.6944164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8089777) q[0];
sx q[0];
rz(-1.806145) q[0];
sx q[0];
rz(0.72840011) q[0];
rz(1.2737466) q[1];
sx q[1];
rz(-1.961901) q[1];
sx q[1];
rz(1.5789998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9804763) q[0];
sx q[0];
rz(-0.022810629) q[0];
sx q[0];
rz(-1.7177607) q[0];
rz(-pi) q[1];
rz(1.4995826) q[2];
sx q[2];
rz(-0.54909583) q[2];
sx q[2];
rz(2.6656815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3167446) q[1];
sx q[1];
rz(-1.7051175) q[1];
sx q[1];
rz(-2.9243241) q[1];
rz(-pi) q[2];
rz(-1.9779491) q[3];
sx q[3];
rz(-0.71135697) q[3];
sx q[3];
rz(2.9748981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1100715) q[2];
sx q[2];
rz(-1.4036125) q[2];
sx q[2];
rz(2.9384379) q[2];
rz(-1.7848232) q[3];
sx q[3];
rz(-1.2099096) q[3];
sx q[3];
rz(1.8206966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(0.88306952) q[0];
sx q[0];
rz(-1.8061545) q[0];
sx q[0];
rz(-1.6254599) q[0];
rz(-2.7698703) q[1];
sx q[1];
rz(-2.5860791) q[1];
sx q[1];
rz(0.98977596) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3759339) q[0];
sx q[0];
rz(-1.759937) q[0];
sx q[0];
rz(2.2098757) q[0];
rz(-pi) q[1];
rz(0.7725352) q[2];
sx q[2];
rz(-1.4896817) q[2];
sx q[2];
rz(-0.81743956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6367148) q[1];
sx q[1];
rz(-1.9002007) q[1];
sx q[1];
rz(-0.18784951) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2618213) q[3];
sx q[3];
rz(-1.1011964) q[3];
sx q[3];
rz(2.0275786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61701361) q[2];
sx q[2];
rz(-2.9144574) q[2];
sx q[2];
rz(-3.0830834) q[2];
rz(1.2837563) q[3];
sx q[3];
rz(-0.82710281) q[3];
sx q[3];
rz(-1.1684928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44130317) q[0];
sx q[0];
rz(-1.0463725) q[0];
sx q[0];
rz(2.6562279) q[0];
rz(-2.6742588) q[1];
sx q[1];
rz(-0.58715564) q[1];
sx q[1];
rz(2.4553518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.195226) q[0];
sx q[0];
rz(-1.2513046) q[0];
sx q[0];
rz(0.736306) q[0];
rz(2.9250018) q[2];
sx q[2];
rz(-2.1859043) q[2];
sx q[2];
rz(-1.5647581) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.972922) q[1];
sx q[1];
rz(-1.0712364) q[1];
sx q[1];
rz(1.023175) q[1];
rz(-0.83250203) q[3];
sx q[3];
rz(-2.0568536) q[3];
sx q[3];
rz(1.488369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.301959) q[2];
sx q[2];
rz(-2.0531824) q[2];
sx q[2];
rz(-2.9097617) q[2];
rz(-0.92136446) q[3];
sx q[3];
rz(-1.5270343) q[3];
sx q[3];
rz(3.11619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63474083) q[0];
sx q[0];
rz(-1.0036108) q[0];
sx q[0];
rz(2.4941709) q[0];
rz(-1.0415152) q[1];
sx q[1];
rz(-1.2728649) q[1];
sx q[1];
rz(-2.0288859) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1739832) q[0];
sx q[0];
rz(-0.69794535) q[0];
sx q[0];
rz(0.83794727) q[0];
rz(-pi) q[1];
rz(-1.9693807) q[2];
sx q[2];
rz(-2.9867134) q[2];
sx q[2];
rz(1.9914371) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1061411) q[1];
sx q[1];
rz(-1.7454495) q[1];
sx q[1];
rz(-0.56366745) q[1];
x q[2];
rz(-1.1057439) q[3];
sx q[3];
rz(-1.5616284) q[3];
sx q[3];
rz(-2.0818215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.57704321) q[2];
sx q[2];
rz(-2.3304522) q[2];
sx q[2];
rz(0.88456336) q[2];
rz(-2.1029643) q[3];
sx q[3];
rz(-1.1218772) q[3];
sx q[3];
rz(1.5326726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7811979) q[0];
sx q[0];
rz(-2.9304657) q[0];
sx q[0];
rz(-1.1398844) q[0];
rz(-0.79849157) q[1];
sx q[1];
rz(-1.4356828) q[1];
sx q[1];
rz(3.0671157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89635623) q[0];
sx q[0];
rz(-2.470861) q[0];
sx q[0];
rz(1.1719431) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3788667) q[2];
sx q[2];
rz(-1.8755091) q[2];
sx q[2];
rz(-2.9305127) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0165538) q[1];
sx q[1];
rz(-1.2579664) q[1];
sx q[1];
rz(-2.4683663) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30280717) q[3];
sx q[3];
rz(-1.2330556) q[3];
sx q[3];
rz(2.1827121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8188339) q[2];
sx q[2];
rz(-0.61046159) q[2];
sx q[2];
rz(-2.9739213) q[2];
rz(2.3882315) q[3];
sx q[3];
rz(-1.2246917) q[3];
sx q[3];
rz(-1.0838449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7849671) q[0];
sx q[0];
rz(-1.6870808) q[0];
sx q[0];
rz(2.0941358) q[0];
rz(0.52182237) q[1];
sx q[1];
rz(-1.9597041) q[1];
sx q[1];
rz(0.80400115) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76011953) q[0];
sx q[0];
rz(-1.1520885) q[0];
sx q[0];
rz(2.9497854) q[0];
rz(1.8487156) q[2];
sx q[2];
rz(-2.5679595) q[2];
sx q[2];
rz(2.3282027) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8447574) q[1];
sx q[1];
rz(-1.3073107) q[1];
sx q[1];
rz(-2.6729463) q[1];
rz(-0.45553546) q[3];
sx q[3];
rz(-2.3856294) q[3];
sx q[3];
rz(2.8433133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1890586) q[2];
sx q[2];
rz(-1.6510094) q[2];
sx q[2];
rz(2.789759) q[2];
rz(-3.0686038) q[3];
sx q[3];
rz(-1.1806386) q[3];
sx q[3];
rz(-2.4019057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9214582) q[0];
sx q[0];
rz(-2.298482) q[0];
sx q[0];
rz(-1.5716918) q[0];
rz(-2.4193071) q[1];
sx q[1];
rz(-1.1415569) q[1];
sx q[1];
rz(-1.3977745) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.777228) q[0];
sx q[0];
rz(-1.9731063) q[0];
sx q[0];
rz(-1.8437632) q[0];
rz(1.7977925) q[2];
sx q[2];
rz(-2.2637061) q[2];
sx q[2];
rz(2.5347112) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.367041) q[1];
sx q[1];
rz(-1.4907717) q[1];
sx q[1];
rz(-2.5888279) q[1];
rz(-pi) q[2];
rz(-1.8595376) q[3];
sx q[3];
rz(-2.7061314) q[3];
sx q[3];
rz(-2.9277152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34652823) q[2];
sx q[2];
rz(-1.3115735) q[2];
sx q[2];
rz(1.2947003) q[2];
rz(-2.0857701) q[3];
sx q[3];
rz(-1.0603696) q[3];
sx q[3];
rz(1.0052217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9697363) q[0];
sx q[0];
rz(-2.2315401) q[0];
sx q[0];
rz(1.6730614) q[0];
rz(-0.83364529) q[1];
sx q[1];
rz(-1.2794762) q[1];
sx q[1];
rz(1.5247482) q[1];
rz(2.7387754) q[2];
sx q[2];
rz(-2.2765712) q[2];
sx q[2];
rz(1.2650875) q[2];
rz(0.41224319) q[3];
sx q[3];
rz(-1.1664248) q[3];
sx q[3];
rz(2.5965367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
