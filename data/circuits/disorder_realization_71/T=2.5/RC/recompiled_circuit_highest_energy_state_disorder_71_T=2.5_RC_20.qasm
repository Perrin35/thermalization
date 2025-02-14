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
rz(-2.5833997) q[0];
sx q[0];
rz(-0.71891958) q[0];
sx q[0];
rz(-0.67607003) q[0];
rz(0.27322912) q[1];
sx q[1];
rz(4.0987045) q[1];
sx q[1];
rz(10.336032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9118496) q[0];
sx q[0];
rz(-1.5485244) q[0];
sx q[0];
rz(-0.020072083) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1315795) q[2];
sx q[2];
rz(-2.5160033) q[2];
sx q[2];
rz(2.6715476) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75685749) q[1];
sx q[1];
rz(-0.57418203) q[1];
sx q[1];
rz(0.059348182) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6966713) q[3];
sx q[3];
rz(-2.3231703) q[3];
sx q[3];
rz(-2.6361806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91997826) q[2];
sx q[2];
rz(-2.6875434) q[2];
sx q[2];
rz(-0.81217074) q[2];
rz(3.0657366) q[3];
sx q[3];
rz(-0.64801884) q[3];
sx q[3];
rz(-2.5465452) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6512063) q[0];
sx q[0];
rz(-2.6631329) q[0];
sx q[0];
rz(1.2745717) q[0];
rz(-1.6365341) q[1];
sx q[1];
rz(-2.7795364) q[1];
sx q[1];
rz(1.9777745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75117373) q[0];
sx q[0];
rz(-2.0295706) q[0];
sx q[0];
rz(1.0659046) q[0];
rz(2.4435448) q[2];
sx q[2];
rz(-0.62605941) q[2];
sx q[2];
rz(-1.2241838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3006176) q[1];
sx q[1];
rz(-1.3625486) q[1];
sx q[1];
rz(0.62842259) q[1];
x q[2];
rz(2.8234661) q[3];
sx q[3];
rz(-1.9788747) q[3];
sx q[3];
rz(-2.7216743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.290648) q[2];
sx q[2];
rz(-3.1182351) q[2];
sx q[2];
rz(1.5604431) q[2];
rz(-2.9110939) q[3];
sx q[3];
rz(-2.0932308) q[3];
sx q[3];
rz(0.44052625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5480963) q[0];
sx q[0];
rz(-0.31262147) q[0];
sx q[0];
rz(3.0189959) q[0];
rz(-0.7971881) q[1];
sx q[1];
rz(-1.0120579) q[1];
sx q[1];
rz(-3.0534993) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6428292) q[0];
sx q[0];
rz(-0.82422148) q[0];
sx q[0];
rz(2.7902664) q[0];
rz(-pi) q[1];
rz(1.7983059) q[2];
sx q[2];
rz(-1.6517793) q[2];
sx q[2];
rz(-1.1151659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56386197) q[1];
sx q[1];
rz(-2.9363765) q[1];
sx q[1];
rz(-0.83744253) q[1];
rz(1.9205536) q[3];
sx q[3];
rz(-1.0579374) q[3];
sx q[3];
rz(-2.8347391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.88283551) q[2];
sx q[2];
rz(-1.5828524) q[2];
sx q[2];
rz(1.8176414) q[2];
rz(2.6435408) q[3];
sx q[3];
rz(-0.96314722) q[3];
sx q[3];
rz(0.74670416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28904706) q[0];
sx q[0];
rz(-2.9883224) q[0];
sx q[0];
rz(-0.14719851) q[0];
rz(-0.64951253) q[1];
sx q[1];
rz(-1.5107061) q[1];
sx q[1];
rz(-1.6726327) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8513059) q[0];
sx q[0];
rz(-1.3477579) q[0];
sx q[0];
rz(-3.1342491) q[0];
rz(-pi) q[1];
rz(1.1873522) q[2];
sx q[2];
rz(-2.6818775) q[2];
sx q[2];
rz(-1.5292668) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8115639) q[1];
sx q[1];
rz(-2.4061476) q[1];
sx q[1];
rz(-1.0388908) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58261915) q[3];
sx q[3];
rz(-1.1938059) q[3];
sx q[3];
rz(2.2347799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9327675) q[2];
sx q[2];
rz(-0.494445) q[2];
sx q[2];
rz(2.8998568) q[2];
rz(2.406481) q[3];
sx q[3];
rz(-1.7416411) q[3];
sx q[3];
rz(-0.66489768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46297896) q[0];
sx q[0];
rz(-1.047387) q[0];
sx q[0];
rz(-3.0188634) q[0];
rz(2.2562476) q[1];
sx q[1];
rz(-1.0095936) q[1];
sx q[1];
rz(0.56040323) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3353676) q[0];
sx q[0];
rz(-2.0809552) q[0];
sx q[0];
rz(-0.2038184) q[0];
rz(1.3483062) q[2];
sx q[2];
rz(-1.280539) q[2];
sx q[2];
rz(-3.0477216) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0593157) q[1];
sx q[1];
rz(-0.53566414) q[1];
sx q[1];
rz(-2.0945063) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6953648) q[3];
sx q[3];
rz(-2.0110907) q[3];
sx q[3];
rz(-1.4780731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6178599) q[2];
sx q[2];
rz(-2.3475519) q[2];
sx q[2];
rz(-0.20475556) q[2];
rz(-2.2023885) q[3];
sx q[3];
rz(-1.3699646) q[3];
sx q[3];
rz(3.052616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8646249) q[0];
sx q[0];
rz(-0.11700103) q[0];
sx q[0];
rz(-0.65698874) q[0];
rz(0.14343801) q[1];
sx q[1];
rz(-1.5564432) q[1];
sx q[1];
rz(2.4030446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6929106) q[0];
sx q[0];
rz(-1.686272) q[0];
sx q[0];
rz(-0.67489745) q[0];
rz(1.9759278) q[2];
sx q[2];
rz(-2.6141593) q[2];
sx q[2];
rz(-2.1857174) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2966531) q[1];
sx q[1];
rz(-0.9228188) q[1];
sx q[1];
rz(-0.39877994) q[1];
rz(-pi) q[2];
rz(-1.012771) q[3];
sx q[3];
rz(-0.9415365) q[3];
sx q[3];
rz(0.93782872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2584381) q[2];
sx q[2];
rz(-0.65885764) q[2];
sx q[2];
rz(3.1385885) q[2];
rz(2.352412) q[3];
sx q[3];
rz(-0.3843669) q[3];
sx q[3];
rz(1.9743617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0810735) q[0];
sx q[0];
rz(-0.34927148) q[0];
sx q[0];
rz(2.9935167) q[0];
rz(-2.608346) q[1];
sx q[1];
rz(-2.8956469) q[1];
sx q[1];
rz(-1.6291133) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142198) q[0];
sx q[0];
rz(-0.54410579) q[0];
sx q[0];
rz(-2.6065488) q[0];
x q[1];
rz(0.81432366) q[2];
sx q[2];
rz(-2.0481128) q[2];
sx q[2];
rz(2.5392591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.80485) q[1];
sx q[1];
rz(-0.36174332) q[1];
sx q[1];
rz(-1.2233578) q[1];
rz(-pi) q[2];
rz(-1.809172) q[3];
sx q[3];
rz(-2.4052103) q[3];
sx q[3];
rz(0.045698085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44275722) q[2];
sx q[2];
rz(-1.7008702) q[2];
sx q[2];
rz(0.093078144) q[2];
rz(-0.062151521) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(1.2232346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020697866) q[0];
sx q[0];
rz(-2.5466205) q[0];
sx q[0];
rz(-2.4769532) q[0];
rz(-1.3497893) q[1];
sx q[1];
rz(-0.87769687) q[1];
sx q[1];
rz(2.452449) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0040379077) q[0];
sx q[0];
rz(-1.8144286) q[0];
sx q[0];
rz(1.451127) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87920636) q[2];
sx q[2];
rz(-1.1779281) q[2];
sx q[2];
rz(-0.61516064) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2011678) q[1];
sx q[1];
rz(-1.42527) q[1];
sx q[1];
rz(0.31647233) q[1];
x q[2];
rz(-1.6313305) q[3];
sx q[3];
rz(-1.7352967) q[3];
sx q[3];
rz(-1.9609083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.050921116) q[2];
sx q[2];
rz(-2.4854269) q[2];
sx q[2];
rz(1.7655168) q[2];
rz(-1.9164267) q[3];
sx q[3];
rz(-0.71198946) q[3];
sx q[3];
rz(-0.041697748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0390778) q[0];
sx q[0];
rz(-0.044487655) q[0];
sx q[0];
rz(-0.59120375) q[0];
rz(0.2419596) q[1];
sx q[1];
rz(-1.0958902) q[1];
sx q[1];
rz(2.718149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5553143) q[0];
sx q[0];
rz(-1.269377) q[0];
sx q[0];
rz(-0.31762808) q[0];
x q[1];
rz(2.9482163) q[2];
sx q[2];
rz(-2.6453553) q[2];
sx q[2];
rz(-1.9691182) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4462699) q[1];
sx q[1];
rz(-2.5188418) q[1];
sx q[1];
rz(-1.0458617) q[1];
x q[2];
rz(2.5245111) q[3];
sx q[3];
rz(-1.9071731) q[3];
sx q[3];
rz(1.468889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.55687904) q[2];
sx q[2];
rz(-0.60606474) q[2];
sx q[2];
rz(1.9752183) q[2];
rz(2.0139458) q[3];
sx q[3];
rz(-2.1731264) q[3];
sx q[3];
rz(2.6166272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.128085) q[0];
sx q[0];
rz(-0.43376827) q[0];
sx q[0];
rz(-2.4684913) q[0];
rz(2.290944) q[1];
sx q[1];
rz(-1.3530082) q[1];
sx q[1];
rz(2.8180715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7009525) q[0];
sx q[0];
rz(-1.8201168) q[0];
sx q[0];
rz(-1.1855769) q[0];
rz(-pi) q[1];
rz(-2.0704248) q[2];
sx q[2];
rz(-0.53348225) q[2];
sx q[2];
rz(-2.1738652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0464152) q[1];
sx q[1];
rz(-0.84734619) q[1];
sx q[1];
rz(1.1208388) q[1];
rz(-pi) q[2];
rz(0.011908493) q[3];
sx q[3];
rz(-2.5577196) q[3];
sx q[3];
rz(1.5425105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2181776) q[2];
sx q[2];
rz(-2.9604993) q[2];
sx q[2];
rz(-0.072048135) q[2];
rz(-1.0142903) q[3];
sx q[3];
rz(-2.225596) q[3];
sx q[3];
rz(-0.54916507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88076787) q[0];
sx q[0];
rz(-1.4243955) q[0];
sx q[0];
rz(-2.0212174) q[0];
rz(-0.33521677) q[1];
sx q[1];
rz(-2.4991279) q[1];
sx q[1];
rz(2.0229708) q[1];
rz(1.3041244) q[2];
sx q[2];
rz(-2.647807) q[2];
sx q[2];
rz(-0.37995445) q[2];
rz(-2.9670197) q[3];
sx q[3];
rz(-0.1452419) q[3];
sx q[3];
rz(-2.8066842) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
