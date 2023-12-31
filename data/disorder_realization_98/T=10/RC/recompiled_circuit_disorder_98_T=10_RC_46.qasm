OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(4.1376576) q[0];
sx q[0];
rz(7.1538038) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(-0.14970782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5538841) q[0];
sx q[0];
rz(-0.93548488) q[0];
sx q[0];
rz(2.5250838) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1120464) q[2];
sx q[2];
rz(-1.1272578) q[2];
sx q[2];
rz(-2.8461547) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.522361) q[1];
sx q[1];
rz(-1.751096) q[1];
sx q[1];
rz(0.038832263) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2655067) q[3];
sx q[3];
rz(-0.57170924) q[3];
sx q[3];
rz(1.3523462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8831138) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(-2.5090704) q[0];
rz(-2.6951492) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(2.4893563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95603847) q[0];
sx q[0];
rz(-1.599405) q[0];
sx q[0];
rz(-1.5584598) q[0];
rz(2.9035283) q[2];
sx q[2];
rz(-1.2859584) q[2];
sx q[2];
rz(2.970398) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10416874) q[1];
sx q[1];
rz(-1.2836604) q[1];
sx q[1];
rz(-1.8038521) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2964301) q[3];
sx q[3];
rz(-2.3361428) q[3];
sx q[3];
rz(2.1243387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5941045) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(-1.1616421) q[2];
rz(-1.15796) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(-0.85025775) q[0];
rz(-0.49750528) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(-1.7920378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3916546) q[0];
sx q[0];
rz(-0.91021252) q[0];
sx q[0];
rz(1.2664938) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0850052) q[2];
sx q[2];
rz(-3.0658256) q[2];
sx q[2];
rz(-2.5839992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.459356) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(-2.6911246) q[1];
rz(3.1317741) q[3];
sx q[3];
rz(-1.1447284) q[3];
sx q[3];
rz(-0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(0.30113014) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6501453) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(-1.4105463) q[0];
rz(2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(3.1052123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2934389) q[0];
sx q[0];
rz(-1.7833976) q[0];
sx q[0];
rz(2.2178177) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47469791) q[2];
sx q[2];
rz(-0.63650741) q[2];
sx q[2];
rz(-1.1324901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.42602793) q[1];
sx q[1];
rz(-1.4392816) q[1];
sx q[1];
rz(0.87042602) q[1];
x q[2];
rz(-0.39156885) q[3];
sx q[3];
rz(-0.25300004) q[3];
sx q[3];
rz(-2.5391425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(-1.1882163) q[2];
rz(0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.7017512) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(0.16648509) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(0.23434815) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7097965) q[0];
sx q[0];
rz(-2.095982) q[0];
sx q[0];
rz(1.1774506) q[0];
rz(-pi) q[1];
x q[1];
rz(2.064803) q[2];
sx q[2];
rz(-1.6550118) q[2];
sx q[2];
rz(-2.6780724) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25917945) q[1];
sx q[1];
rz(-0.83173527) q[1];
sx q[1];
rz(-1.3684567) q[1];
rz(-0.078951051) q[3];
sx q[3];
rz(-1.8936994) q[3];
sx q[3];
rz(-0.80432804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4328737) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(-2.9210572) q[2];
rz(0.43705127) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(-0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7261312) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(-0.81714001) q[0];
rz(0.56610402) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(1.1436499) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7276579) q[0];
sx q[0];
rz(-1.2000788) q[0];
sx q[0];
rz(1.6137705) q[0];
rz(1.3104865) q[2];
sx q[2];
rz(-0.91432768) q[2];
sx q[2];
rz(1.7973763) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5551344) q[1];
sx q[1];
rz(-2.0124334) q[1];
sx q[1];
rz(3.0416136) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51795824) q[3];
sx q[3];
rz(-1.3266139) q[3];
sx q[3];
rz(1.1740008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75366655) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(0.4894408) q[2];
rz(0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(-0.42603809) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5722826) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(1.2868767) q[0];
rz(0.6634179) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(1.9082665) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9702643) q[0];
sx q[0];
rz(-1.5597222) q[0];
sx q[0];
rz(-1.1867255) q[0];
rz(0.8823422) q[2];
sx q[2];
rz(-0.8493648) q[2];
sx q[2];
rz(-1.7989858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.308097) q[1];
sx q[1];
rz(-1.6819681) q[1];
sx q[1];
rz(1.584946) q[1];
rz(-pi) q[2];
rz(2.535378) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(-1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.037501637) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(1.6960779) q[0];
rz(-0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(-1.8833556) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1930775) q[0];
sx q[0];
rz(-0.45398871) q[0];
sx q[0];
rz(1.8817188) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3636742) q[2];
sx q[2];
rz(-0.81233835) q[2];
sx q[2];
rz(1.1493491) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82842365) q[1];
sx q[1];
rz(-1.0272044) q[1];
sx q[1];
rz(-1.6924752) q[1];
rz(-pi) q[2];
rz(-2.1636837) q[3];
sx q[3];
rz(-2.3330354) q[3];
sx q[3];
rz(-0.18329328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(-0.56274596) q[2];
rz(0.19966666) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(2.8905706) q[0];
rz(-2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(-0.02773157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8160307) q[0];
sx q[0];
rz(-1.0849909) q[0];
sx q[0];
rz(2.3150139) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8717732) q[2];
sx q[2];
rz(-1.36424) q[2];
sx q[2];
rz(2.2733462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7293538) q[1];
sx q[1];
rz(-2.5322399) q[1];
sx q[1];
rz(2.9669697) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1157007) q[3];
sx q[3];
rz(-1.6373487) q[3];
sx q[3];
rz(2.6641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.015908265) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(-0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(2.2843602) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(-2.8840816) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2902381) q[0];
sx q[0];
rz(-2.5332753) q[0];
sx q[0];
rz(0.71382199) q[0];
x q[1];
rz(2.4117878) q[2];
sx q[2];
rz(-0.68241718) q[2];
sx q[2];
rz(-0.44621106) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35044893) q[1];
sx q[1];
rz(-2.2669683) q[1];
sx q[1];
rz(-1.5894366) q[1];
rz(-2.0503644) q[3];
sx q[3];
rz(-1.3551095) q[3];
sx q[3];
rz(-1.8857764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8132849) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(-2.5349687) q[2];
rz(2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(2.9150302) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(0.67129927) q[2];
sx q[2];
rz(-0.57325267) q[2];
sx q[2];
rz(-0.43043955) q[2];
rz(-1.0874891) q[3];
sx q[3];
rz(-1.3369505) q[3];
sx q[3];
rz(-1.7575775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
