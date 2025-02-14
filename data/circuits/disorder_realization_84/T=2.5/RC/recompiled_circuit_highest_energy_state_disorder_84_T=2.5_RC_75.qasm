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
rz(-0.62491971) q[0];
sx q[0];
rz(4.9486296) q[0];
sx q[0];
rz(9.3715342) q[0];
rz(-2.6553395) q[1];
sx q[1];
rz(-3.0142205) q[1];
sx q[1];
rz(-1.706634) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20986461) q[0];
sx q[0];
rz(-1.7057944) q[0];
sx q[0];
rz(-2.854706) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8743319) q[2];
sx q[2];
rz(-1.4556985) q[2];
sx q[2];
rz(2.5989344) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9588152) q[1];
sx q[1];
rz(-2.7441886) q[1];
sx q[1];
rz(0.13840492) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1987183) q[3];
sx q[3];
rz(-0.74340313) q[3];
sx q[3];
rz(-1.2974031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3367553) q[2];
sx q[2];
rz(-1.3433604) q[2];
sx q[2];
rz(0.27377823) q[2];
rz(1.9233507) q[3];
sx q[3];
rz(-0.67353407) q[3];
sx q[3];
rz(0.28606733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186721) q[0];
sx q[0];
rz(-2.3781222) q[0];
sx q[0];
rz(0.77962312) q[0];
rz(-0.66462213) q[1];
sx q[1];
rz(-1.2088935) q[1];
sx q[1];
rz(1.2453311) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1561403) q[0];
sx q[0];
rz(-3.0110571) q[0];
sx q[0];
rz(-2.0640316) q[0];
rz(-pi) q[1];
rz(-2.2822218) q[2];
sx q[2];
rz(-1.6409163) q[2];
sx q[2];
rz(-1.9944348) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9960963) q[1];
sx q[1];
rz(-0.22095535) q[1];
sx q[1];
rz(1.2902106) q[1];
rz(1.0118724) q[3];
sx q[3];
rz(-1.3959612) q[3];
sx q[3];
rz(-1.2993985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7094884) q[2];
sx q[2];
rz(-0.69807845) q[2];
sx q[2];
rz(-0.049987642) q[2];
rz(-1.6677808) q[3];
sx q[3];
rz(-1.4155017) q[3];
sx q[3];
rz(0.93562359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98899984) q[0];
sx q[0];
rz(-2.279156) q[0];
sx q[0];
rz(-1.9631901) q[0];
rz(-1.1475457) q[1];
sx q[1];
rz(-0.18745628) q[1];
sx q[1];
rz(2.5386834) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99166223) q[0];
sx q[0];
rz(-1.4105182) q[0];
sx q[0];
rz(3.0383238) q[0];
x q[1];
rz(-2.2608537) q[2];
sx q[2];
rz(-1.1549486) q[2];
sx q[2];
rz(1.2643472) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6223095) q[1];
sx q[1];
rz(-1.9819248) q[1];
sx q[1];
rz(2.4101188) q[1];
rz(-pi) q[2];
rz(-2.1985377) q[3];
sx q[3];
rz(-2.3585417) q[3];
sx q[3];
rz(-0.32367009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19043645) q[2];
sx q[2];
rz(-2.3556605) q[2];
sx q[2];
rz(2.7583165) q[2];
rz(1.3112618) q[3];
sx q[3];
rz(-1.0878599) q[3];
sx q[3];
rz(3.0200628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.424778) q[0];
sx q[0];
rz(-2.1489693) q[0];
sx q[0];
rz(0.92856652) q[0];
rz(-0.83944744) q[1];
sx q[1];
rz(-1.8239559) q[1];
sx q[1];
rz(2.3323257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081079) q[0];
sx q[0];
rz(-1.9376457) q[0];
sx q[0];
rz(-2.1686694) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4014612) q[2];
sx q[2];
rz(-1.9953097) q[2];
sx q[2];
rz(-2.9237539) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41701554) q[1];
sx q[1];
rz(-0.25522403) q[1];
sx q[1];
rz(-1.8468813) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.094481) q[3];
sx q[3];
rz(-2.4274656) q[3];
sx q[3];
rz(-1.3830101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3889918) q[2];
sx q[2];
rz(-1.6099124) q[2];
sx q[2];
rz(-1.6947702) q[2];
rz(-1.1270479) q[3];
sx q[3];
rz(-0.73486745) q[3];
sx q[3];
rz(2.5126422) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22555722) q[0];
sx q[0];
rz(-1.1981523) q[0];
sx q[0];
rz(1.4300562) q[0];
rz(-0.23722181) q[1];
sx q[1];
rz(-2.1784541) q[1];
sx q[1];
rz(0.29475862) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4859516) q[0];
sx q[0];
rz(-1.427934) q[0];
sx q[0];
rz(-2.7014707) q[0];
rz(-0.90845256) q[2];
sx q[2];
rz(-1.2725432) q[2];
sx q[2];
rz(2.3191888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7511616) q[1];
sx q[1];
rz(-2.0195578) q[1];
sx q[1];
rz(0.46350355) q[1];
rz(-pi) q[2];
rz(0.57394694) q[3];
sx q[3];
rz(-1.6960295) q[3];
sx q[3];
rz(-0.67170152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1416867) q[2];
sx q[2];
rz(-0.20076951) q[2];
sx q[2];
rz(3.0086009) q[2];
rz(0.76987949) q[3];
sx q[3];
rz(-1.2917638) q[3];
sx q[3];
rz(2.8225115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8517476) q[0];
sx q[0];
rz(-2.802749) q[0];
sx q[0];
rz(-1.7293365) q[0];
rz(0.051636592) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(-0.19270611) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.76336) q[0];
sx q[0];
rz(-2.6274649) q[0];
sx q[0];
rz(-3.1030416) q[0];
x q[1];
rz(2.4102845) q[2];
sx q[2];
rz(-1.0885065) q[2];
sx q[2];
rz(2.363229) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6289857) q[1];
sx q[1];
rz(-0.64388212) q[1];
sx q[1];
rz(2.3888247) q[1];
rz(-pi) q[2];
rz(-2.5000743) q[3];
sx q[3];
rz(-1.2605485) q[3];
sx q[3];
rz(-1.1603242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.073033832) q[2];
sx q[2];
rz(-1.8517588) q[2];
sx q[2];
rz(1.3812836) q[2];
rz(-0.51501385) q[3];
sx q[3];
rz(-0.88740715) q[3];
sx q[3];
rz(1.380434) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16159049) q[0];
sx q[0];
rz(-1.1450293) q[0];
sx q[0];
rz(2.7951151) q[0];
rz(-2.1222291) q[1];
sx q[1];
rz(-0.66277021) q[1];
sx q[1];
rz(1.4541218) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1145059) q[0];
sx q[0];
rz(-2.1505304) q[0];
sx q[0];
rz(-2.113335) q[0];
rz(-pi) q[1];
rz(-2.2367661) q[2];
sx q[2];
rz(-2.3695393) q[2];
sx q[2];
rz(2.0029298) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7123966) q[1];
sx q[1];
rz(-2.0919261) q[1];
sx q[1];
rz(-2.8217535) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.725127) q[3];
sx q[3];
rz(-1.6701506) q[3];
sx q[3];
rz(-1.9973444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.60160294) q[2];
sx q[2];
rz(-1.3301962) q[2];
sx q[2];
rz(2.9023564) q[2];
rz(0.38419497) q[3];
sx q[3];
rz(-0.68238634) q[3];
sx q[3];
rz(2.189883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9397028) q[0];
sx q[0];
rz(-1.6738142) q[0];
sx q[0];
rz(1.3772759) q[0];
rz(1.9210303) q[1];
sx q[1];
rz(-2.2790597) q[1];
sx q[1];
rz(0.16622226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2234874) q[0];
sx q[0];
rz(-1.1150196) q[0];
sx q[0];
rz(2.8546643) q[0];
rz(-2.407309) q[2];
sx q[2];
rz(-1.3314651) q[2];
sx q[2];
rz(-0.83877968) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0609445) q[1];
sx q[1];
rz(-1.6014487) q[1];
sx q[1];
rz(0.21126195) q[1];
rz(-pi) q[2];
rz(-2.1020426) q[3];
sx q[3];
rz(-2.663718) q[3];
sx q[3];
rz(0.22051624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7938457) q[2];
sx q[2];
rz(-0.85024992) q[2];
sx q[2];
rz(-2.094685) q[2];
rz(1.2547803) q[3];
sx q[3];
rz(-2.051765) q[3];
sx q[3];
rz(1.8449529) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60085249) q[0];
sx q[0];
rz(-2.7678601) q[0];
sx q[0];
rz(-1.1131713) q[0];
rz(-2.3241849) q[1];
sx q[1];
rz(-1.4459041) q[1];
sx q[1];
rz(2.7730952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2042646) q[0];
sx q[0];
rz(-2.6952792) q[0];
sx q[0];
rz(2.242779) q[0];
x q[1];
rz(-0.13780528) q[2];
sx q[2];
rz(-2.3944582) q[2];
sx q[2];
rz(-1.8835619) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6121713) q[1];
sx q[1];
rz(-2.2670806) q[1];
sx q[1];
rz(1.2207635) q[1];
x q[2];
rz(-1.0115252) q[3];
sx q[3];
rz(-2.2314921) q[3];
sx q[3];
rz(-1.504577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.14785279) q[2];
sx q[2];
rz(-0.61823121) q[2];
sx q[2];
rz(-1.3573793) q[2];
rz(2.3944858) q[3];
sx q[3];
rz(-2.1785469) q[3];
sx q[3];
rz(-2.4660306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6152182) q[0];
sx q[0];
rz(-0.45658657) q[0];
sx q[0];
rz(2.0988462) q[0];
rz(-2.3710947) q[1];
sx q[1];
rz(-2.1310525) q[1];
sx q[1];
rz(2.3811293) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2337991) q[0];
sx q[0];
rz(-1.9007508) q[0];
sx q[0];
rz(2.8889662) q[0];
rz(3.1030948) q[2];
sx q[2];
rz(-1.7506934) q[2];
sx q[2];
rz(-0.28598374) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.14899602) q[1];
sx q[1];
rz(-1.1482052) q[1];
sx q[1];
rz(-1.484297) q[1];
x q[2];
rz(2.4861927) q[3];
sx q[3];
rz(-1.8701347) q[3];
sx q[3];
rz(1.6884402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9383135) q[2];
sx q[2];
rz(-0.61050296) q[2];
sx q[2];
rz(-2.7321613) q[2];
rz(-1.5583386) q[3];
sx q[3];
rz(-2.6810472) q[3];
sx q[3];
rz(-0.62622768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4938477) q[0];
sx q[0];
rz(-2.0997601) q[0];
sx q[0];
rz(0.44152015) q[0];
rz(-1.5589177) q[1];
sx q[1];
rz(-1.1987004) q[1];
sx q[1];
rz(-0.62335062) q[1];
rz(-1.6292439) q[2];
sx q[2];
rz(-1.4679906) q[2];
sx q[2];
rz(1.0124258) q[2];
rz(2.6473896) q[3];
sx q[3];
rz(-0.27402607) q[3];
sx q[3];
rz(-1.5370697) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
