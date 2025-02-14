OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71896267) q[0];
sx q[0];
rz(-0.29932061) q[0];
sx q[0];
rz(0.49459767) q[0];
rz(1.142113) q[1];
sx q[1];
rz(-1.0057058) q[1];
sx q[1];
rz(-2.0118654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8324234) q[0];
sx q[0];
rz(-2.4465843) q[0];
sx q[0];
rz(0.2485991) q[0];
rz(0.45290516) q[2];
sx q[2];
rz(-1.1158841) q[2];
sx q[2];
rz(-0.95547966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59334457) q[1];
sx q[1];
rz(-1.3181837) q[1];
sx q[1];
rz(-2.3139364) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0315597) q[3];
sx q[3];
rz(-0.72537106) q[3];
sx q[3];
rz(1.6417208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8226681) q[2];
sx q[2];
rz(-1.319004) q[2];
sx q[2];
rz(-0.85282105) q[2];
rz(1.3301814) q[3];
sx q[3];
rz(-0.69555247) q[3];
sx q[3];
rz(-0.046796355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3332719) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(1.464123) q[0];
rz(1.4785712) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(1.2082072) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0532329) q[0];
sx q[0];
rz(-0.7506943) q[0];
sx q[0];
rz(0.73221598) q[0];
x q[1];
rz(1.9764429) q[2];
sx q[2];
rz(-2.6313836) q[2];
sx q[2];
rz(-1.134553) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3390926) q[1];
sx q[1];
rz(-2.8475757) q[1];
sx q[1];
rz(2.0014928) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3974819) q[3];
sx q[3];
rz(-2.0195761) q[3];
sx q[3];
rz(3.0822494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53604424) q[2];
sx q[2];
rz(-2.7216585) q[2];
sx q[2];
rz(1.4844683) q[2];
rz(-0.015497192) q[3];
sx q[3];
rz(-1.9280547) q[3];
sx q[3];
rz(-1.1789471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.991796) q[0];
sx q[0];
rz(-1.3301671) q[0];
sx q[0];
rz(2.8699744) q[0];
rz(0.89667165) q[1];
sx q[1];
rz(-0.50183693) q[1];
sx q[1];
rz(1.2976049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.04682) q[0];
sx q[0];
rz(-1.5808788) q[0];
sx q[0];
rz(-1.1077513) q[0];
rz(0.50432019) q[2];
sx q[2];
rz(-1.9078622) q[2];
sx q[2];
rz(-0.45318174) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0314762) q[1];
sx q[1];
rz(-1.3240485) q[1];
sx q[1];
rz(3.0241248) q[1];
x q[2];
rz(3.1362757) q[3];
sx q[3];
rz(-0.76197366) q[3];
sx q[3];
rz(3.0865106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4012332) q[2];
sx q[2];
rz(-0.77784246) q[2];
sx q[2];
rz(3.0878301) q[2];
rz(-1.8476123) q[3];
sx q[3];
rz(-2.3432178) q[3];
sx q[3];
rz(-2.9062041) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9409222) q[0];
sx q[0];
rz(-0.5680474) q[0];
sx q[0];
rz(-1.6773552) q[0];
rz(-1.1306521) q[1];
sx q[1];
rz(-2.407275) q[1];
sx q[1];
rz(3.0272223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371724) q[0];
sx q[0];
rz(-2.3995598) q[0];
sx q[0];
rz(-2.4993308) q[0];
rz(1.0724154) q[2];
sx q[2];
rz(-0.97626462) q[2];
sx q[2];
rz(-1.4630813) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3072171) q[1];
sx q[1];
rz(-0.67612069) q[1];
sx q[1];
rz(-1.5285049) q[1];
rz(-pi) q[2];
rz(-1.5642197) q[3];
sx q[3];
rz(-1.7688171) q[3];
sx q[3];
rz(-1.978156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6681246) q[2];
sx q[2];
rz(-1.6056085) q[2];
sx q[2];
rz(-0.075695666) q[2];
rz(0.43478742) q[3];
sx q[3];
rz(-1.1554759) q[3];
sx q[3];
rz(3.0619612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39744034) q[0];
sx q[0];
rz(-1.8295153) q[0];
sx q[0];
rz(-0.42006668) q[0];
rz(-2.9282667) q[1];
sx q[1];
rz(-1.7330287) q[1];
sx q[1];
rz(-1.2423645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39109502) q[0];
sx q[0];
rz(-0.89078442) q[0];
sx q[0];
rz(1.0155506) q[0];
x q[1];
rz(-2.5630997) q[2];
sx q[2];
rz(-2.6393386) q[2];
sx q[2];
rz(-1.945221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.99862) q[1];
sx q[1];
rz(-3.0398453) q[1];
sx q[1];
rz(1.4187901) q[1];
x q[2];
rz(2.9202634) q[3];
sx q[3];
rz(-1.5240655) q[3];
sx q[3];
rz(-2.5525301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3199978) q[2];
sx q[2];
rz(-2.2264806) q[2];
sx q[2];
rz(-0.37459174) q[2];
rz(-1.0125259) q[3];
sx q[3];
rz(-0.39207021) q[3];
sx q[3];
rz(-0.64468002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0300765) q[0];
sx q[0];
rz(-2.6426297) q[0];
sx q[0];
rz(0.43055713) q[0];
rz(2.997609) q[1];
sx q[1];
rz(-1.0824243) q[1];
sx q[1];
rz(-0.63124257) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66723204) q[0];
sx q[0];
rz(-0.96275126) q[0];
sx q[0];
rz(-1.505018) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8885884) q[2];
sx q[2];
rz(-0.43995198) q[2];
sx q[2];
rz(0.41340128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5938607) q[1];
sx q[1];
rz(-1.0775078) q[1];
sx q[1];
rz(2.4592722) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33511929) q[3];
sx q[3];
rz(-2.0010173) q[3];
sx q[3];
rz(-0.29725257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1708019) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(-1.7010752) q[2];
rz(-1.3614281) q[3];
sx q[3];
rz(-1.1037339) q[3];
sx q[3];
rz(2.6543999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68706566) q[0];
sx q[0];
rz(-2.8519958) q[0];
sx q[0];
rz(2.2173296) q[0];
rz(2.848792) q[1];
sx q[1];
rz(-1.2288789) q[1];
sx q[1];
rz(-2.3407095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0706884) q[0];
sx q[0];
rz(-1.8278121) q[0];
sx q[0];
rz(1.9563673) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34556324) q[2];
sx q[2];
rz(-0.93831944) q[2];
sx q[2];
rz(-1.8008302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85383254) q[1];
sx q[1];
rz(-1.4656855) q[1];
sx q[1];
rz(-1.5413324) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.597936) q[3];
sx q[3];
rz(-1.5487897) q[3];
sx q[3];
rz(2.1834971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.822829) q[2];
sx q[2];
rz(-1.9209361) q[2];
sx q[2];
rz(1.2804383) q[2];
rz(-2.473623) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(2.7282696) q[3];
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
rz(2.8845344) q[0];
sx q[0];
rz(-1.2385383) q[0];
sx q[0];
rz(1.9588233) q[0];
rz(-2.9044652) q[1];
sx q[1];
rz(-0.10219899) q[1];
sx q[1];
rz(3.0822486) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6366258) q[0];
sx q[0];
rz(-1.6904181) q[0];
sx q[0];
rz(1.2960394) q[0];
rz(1.0419215) q[2];
sx q[2];
rz(-0.9176853) q[2];
sx q[2];
rz(0.83335857) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4044734) q[1];
sx q[1];
rz(-1.307202) q[1];
sx q[1];
rz(3.0332159) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.24712) q[3];
sx q[3];
rz(-1.6442862) q[3];
sx q[3];
rz(-2.2761114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55988971) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(-1.3524559) q[2];
rz(-1.3672359) q[3];
sx q[3];
rz(-1.4227899) q[3];
sx q[3];
rz(0.8555612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2818114) q[0];
sx q[0];
rz(-2.3563522) q[0];
sx q[0];
rz(-2.6819041) q[0];
rz(3.1135318) q[1];
sx q[1];
rz(-1.1496081) q[1];
sx q[1];
rz(-1.185816) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81521124) q[0];
sx q[0];
rz(-1.5259169) q[0];
sx q[0];
rz(2.8759967) q[0];
rz(-pi) q[1];
x q[1];
rz(0.079396768) q[2];
sx q[2];
rz(-1.4697452) q[2];
sx q[2];
rz(3.0011645) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8227905) q[1];
sx q[1];
rz(-1.0981961) q[1];
sx q[1];
rz(-2.5986555) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.288663) q[3];
sx q[3];
rz(-0.97172874) q[3];
sx q[3];
rz(-1.7271068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49843732) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(-0.15677491) q[2];
rz(-2.5458941) q[3];
sx q[3];
rz(-0.34606338) q[3];
sx q[3];
rz(-3.116385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51782411) q[0];
sx q[0];
rz(-2.1110004) q[0];
sx q[0];
rz(0.18381707) q[0];
rz(0.078016438) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(-0.26430166) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28537946) q[0];
sx q[0];
rz(-1.0074248) q[0];
sx q[0];
rz(1.8895288) q[0];
x q[1];
rz(-0.68497505) q[2];
sx q[2];
rz(-1.994418) q[2];
sx q[2];
rz(0.26348235) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.42929998) q[1];
sx q[1];
rz(-1.4136864) q[1];
sx q[1];
rz(0.18204851) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6666344) q[3];
sx q[3];
rz(-2.0639827) q[3];
sx q[3];
rz(0.020190369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8514303) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(0.22872049) q[2];
rz(1.4043572) q[3];
sx q[3];
rz(-2.2018933) q[3];
sx q[3];
rz(3.0586045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1514773) q[0];
sx q[0];
rz(-2.293312) q[0];
sx q[0];
rz(1.5207416) q[0];
rz(2.4304541) q[1];
sx q[1];
rz(-1.3093206) q[1];
sx q[1];
rz(0.73285229) q[1];
rz(1.6522897) q[2];
sx q[2];
rz(-0.50028481) q[2];
sx q[2];
rz(-1.0505983) q[2];
rz(-1.0300954) q[3];
sx q[3];
rz(-1.0813011) q[3];
sx q[3];
rz(-2.4745221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
