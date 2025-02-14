OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.42263) q[0];
sx q[0];
rz(-2.842272) q[0];
sx q[0];
rz(2.646995) q[0];
rz(1.142113) q[1];
sx q[1];
rz(-1.0057058) q[1];
sx q[1];
rz(1.1297273) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30916921) q[0];
sx q[0];
rz(-0.69500837) q[0];
sx q[0];
rz(0.2485991) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45290516) q[2];
sx q[2];
rz(-2.0257086) q[2];
sx q[2];
rz(-2.186113) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2430302) q[1];
sx q[1];
rz(-2.3645325) q[1];
sx q[1];
rz(1.9352566) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2211773) q[3];
sx q[3];
rz(-1.2231881) q[3];
sx q[3];
rz(-2.7915366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8226681) q[2];
sx q[2];
rz(-1.319004) q[2];
sx q[2];
rz(-0.85282105) q[2];
rz(1.8114113) q[3];
sx q[3];
rz(-0.69555247) q[3];
sx q[3];
rz(0.046796355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3332719) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(1.464123) q[0];
rz(-1.4785712) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(1.9333855) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20044662) q[0];
sx q[0];
rz(-2.1028554) q[0];
sx q[0];
rz(-2.1284298) q[0];
x q[1];
rz(1.0958395) q[2];
sx q[2];
rz(-1.37687) q[2];
sx q[2];
rz(0.79481193) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18259163) q[1];
sx q[1];
rz(-1.6920857) q[1];
sx q[1];
rz(1.3022997) q[1];
rz(-1.3974819) q[3];
sx q[3];
rz(-2.0195761) q[3];
sx q[3];
rz(-3.0822494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53604424) q[2];
sx q[2];
rz(-2.7216585) q[2];
sx q[2];
rz(1.4844683) q[2];
rz(-3.1260955) q[3];
sx q[3];
rz(-1.213538) q[3];
sx q[3];
rz(-1.1789471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14979664) q[0];
sx q[0];
rz(-1.8114256) q[0];
sx q[0];
rz(-0.27161828) q[0];
rz(2.244921) q[1];
sx q[1];
rz(-0.50183693) q[1];
sx q[1];
rz(1.8439878) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51894278) q[0];
sx q[0];
rz(-2.033816) q[0];
sx q[0];
rz(0.011269022) q[0];
rz(-1.9515368) q[2];
sx q[2];
rz(-1.0972766) q[2];
sx q[2];
rz(-1.2981594) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5601968) q[1];
sx q[1];
rz(-0.27276688) q[1];
sx q[1];
rz(2.0062937) q[1];
rz(-3.1362757) q[3];
sx q[3];
rz(-2.379619) q[3];
sx q[3];
rz(3.0865106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4012332) q[2];
sx q[2];
rz(-0.77784246) q[2];
sx q[2];
rz(-0.05376251) q[2];
rz(-1.8476123) q[3];
sx q[3];
rz(-0.79837489) q[3];
sx q[3];
rz(2.9062041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2006705) q[0];
sx q[0];
rz(-0.5680474) q[0];
sx q[0];
rz(-1.4642375) q[0];
rz(1.1306521) q[1];
sx q[1];
rz(-2.407275) q[1];
sx q[1];
rz(0.11437036) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30442023) q[0];
sx q[0];
rz(-0.74203287) q[0];
sx q[0];
rz(-0.64226182) q[0];
rz(-pi) q[1];
rz(2.52621) q[2];
sx q[2];
rz(-2.3856731) q[2];
sx q[2];
rz(2.4494954) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3721766) q[1];
sx q[1];
rz(-1.5443364) q[1];
sx q[1];
rz(0.89511223) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.032764445) q[3];
sx q[3];
rz(-0.19812852) q[3];
sx q[3];
rz(-1.1968544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6681246) q[2];
sx q[2];
rz(-1.6056085) q[2];
sx q[2];
rz(-0.075695666) q[2];
rz(2.7068052) q[3];
sx q[3];
rz(-1.1554759) q[3];
sx q[3];
rz(-3.0619612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39744034) q[0];
sx q[0];
rz(-1.3120774) q[0];
sx q[0];
rz(0.42006668) q[0];
rz(0.21332598) q[1];
sx q[1];
rz(-1.408564) q[1];
sx q[1];
rz(1.2423645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5516226) q[0];
sx q[0];
rz(-1.1484572) q[0];
sx q[0];
rz(0.7606272) q[0];
rz(-pi) q[1];
rz(0.43102805) q[2];
sx q[2];
rz(-1.8371474) q[2];
sx q[2];
rz(0.14542994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27659163) q[1];
sx q[1];
rz(-1.5554155) q[1];
sx q[1];
rz(-1.4702142) q[1];
x q[2];
rz(-1.6186938) q[3];
sx q[3];
rz(-1.7918799) q[3];
sx q[3];
rz(2.1703699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8215948) q[2];
sx q[2];
rz(-2.2264806) q[2];
sx q[2];
rz(-0.37459174) q[2];
rz(-1.0125259) q[3];
sx q[3];
rz(-2.7495224) q[3];
sx q[3];
rz(0.64468002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11151611) q[0];
sx q[0];
rz(-0.498963) q[0];
sx q[0];
rz(-2.7110355) q[0];
rz(2.997609) q[1];
sx q[1];
rz(-2.0591683) q[1];
sx q[1];
rz(-2.5103501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66723204) q[0];
sx q[0];
rz(-2.1788414) q[0];
sx q[0];
rz(1.6365746) q[0];
rz(-pi) q[1];
rz(-1.8885884) q[2];
sx q[2];
rz(-0.43995198) q[2];
sx q[2];
rz(-0.41340128) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.54773195) q[1];
sx q[1];
rz(-2.0640848) q[1];
sx q[1];
rz(-0.6823205) q[1];
x q[2];
rz(0.9489551) q[3];
sx q[3];
rz(-0.53882155) q[3];
sx q[3];
rz(2.1486189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1708019) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(1.7010752) q[2];
rz(1.7801646) q[3];
sx q[3];
rz(-2.0378588) q[3];
sx q[3];
rz(-2.6543999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.454527) q[0];
sx q[0];
rz(-0.28959689) q[0];
sx q[0];
rz(2.2173296) q[0];
rz(-2.848792) q[1];
sx q[1];
rz(-1.9127138) q[1];
sx q[1];
rz(-2.3407095) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0706884) q[0];
sx q[0];
rz(-1.3137806) q[0];
sx q[0];
rz(1.1852253) q[0];
rz(-0.34556324) q[2];
sx q[2];
rz(-0.93831944) q[2];
sx q[2];
rz(1.3407624) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.71387163) q[1];
sx q[1];
rz(-1.6000976) q[1];
sx q[1];
rz(-0.10515611) q[1];
rz(-pi) q[2];
rz(-0.022014736) q[3];
sx q[3];
rz(-1.5979294) q[3];
sx q[3];
rz(-2.5282945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.822829) q[2];
sx q[2];
rz(-1.2206565) q[2];
sx q[2];
rz(-1.8611543) q[2];
rz(0.66796962) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(2.7282696) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25705826) q[0];
sx q[0];
rz(-1.9030544) q[0];
sx q[0];
rz(-1.9588233) q[0];
rz(0.23712748) q[1];
sx q[1];
rz(-0.10219899) q[1];
sx q[1];
rz(3.0822486) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33464319) q[0];
sx q[0];
rz(-0.29906267) q[0];
sx q[0];
rz(1.1537667) q[0];
x q[1];
rz(2.4164532) q[2];
sx q[2];
rz(-1.1584917) q[2];
sx q[2];
rz(0.39620846) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8000477) q[1];
sx q[1];
rz(-2.8570685) q[1];
sx q[1];
rz(-1.9519898) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0640934) q[3];
sx q[3];
rz(-1.8935673) q[3];
sx q[3];
rz(-0.68068824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.55988971) q[2];
sx q[2];
rz(-2.6010332) q[2];
sx q[2];
rz(-1.3524559) q[2];
rz(1.7743568) q[3];
sx q[3];
rz(-1.7188027) q[3];
sx q[3];
rz(-0.8555612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2818114) q[0];
sx q[0];
rz(-2.3563522) q[0];
sx q[0];
rz(0.45968858) q[0];
rz(-0.028060878) q[1];
sx q[1];
rz(-1.9919845) q[1];
sx q[1];
rz(1.185816) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3263814) q[0];
sx q[0];
rz(-1.6156757) q[0];
sx q[0];
rz(-0.26559592) q[0];
rz(3.0621959) q[2];
sx q[2];
rz(-1.4697452) q[2];
sx q[2];
rz(-3.0011645) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3188022) q[1];
sx q[1];
rz(-1.0981961) q[1];
sx q[1];
rz(2.5986555) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38717732) q[3];
sx q[3];
rz(-2.4868591) q[3];
sx q[3];
rz(0.93965215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.49843732) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(-2.9848177) q[2];
rz(0.59569851) q[3];
sx q[3];
rz(-2.7955293) q[3];
sx q[3];
rz(3.116385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(2.877291) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3027356) q[0];
sx q[0];
rz(-0.6386916) q[0];
sx q[0];
rz(-0.46052082) q[0];
rz(-pi) q[1];
rz(-2.0980706) q[2];
sx q[2];
rz(-2.1855178) q[2];
sx q[2];
rz(1.6312903) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7122927) q[1];
sx q[1];
rz(-1.7279062) q[1];
sx q[1];
rz(-0.18204851) q[1];
rz(-pi) q[2];
rz(-2.6464858) q[3];
sx q[3];
rz(-1.4864085) q[3];
sx q[3];
rz(-1.5051248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2901624) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(-2.9128722) q[2];
rz(1.4043572) q[3];
sx q[3];
rz(-0.93969932) q[3];
sx q[3];
rz(-3.0586045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9901154) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(-0.71113853) q[1];
sx q[1];
rz(-1.3093206) q[1];
sx q[1];
rz(0.73285229) q[1];
rz(3.0971211) q[2];
sx q[2];
rz(-1.0723249) q[2];
sx q[2];
rz(1.9981801) q[2];
rz(2.5855999) q[3];
sx q[3];
rz(-1.0992194) q[3];
sx q[3];
rz(1.9627375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
