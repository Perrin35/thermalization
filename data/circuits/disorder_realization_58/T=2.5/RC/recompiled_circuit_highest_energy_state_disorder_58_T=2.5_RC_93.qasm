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
rz(2.1287542) q[0];
sx q[0];
rz(-1.1829809) q[0];
sx q[0];
rz(-2.8943789) q[0];
rz(-1.5850868) q[1];
sx q[1];
rz(4.1559846) q[1];
sx q[1];
rz(13.520887) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5069141) q[0];
sx q[0];
rz(-1.5416814) q[0];
sx q[0];
rz(-1.1965189) q[0];
rz(-2.688301) q[2];
sx q[2];
rz(-1.6214091) q[2];
sx q[2];
rz(-2.5091189) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87434972) q[1];
sx q[1];
rz(-1.2552208) q[1];
sx q[1];
rz(-0.79488956) q[1];
rz(-pi) q[2];
rz(1.8778716) q[3];
sx q[3];
rz(-2.319984) q[3];
sx q[3];
rz(-1.6916569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9556094) q[2];
sx q[2];
rz(-2.535203) q[2];
sx q[2];
rz(-2.8601904) q[2];
rz(1.227281) q[3];
sx q[3];
rz(-1.809092) q[3];
sx q[3];
rz(1.9277771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94673741) q[0];
sx q[0];
rz(-2.7342716) q[0];
sx q[0];
rz(2.5545004) q[0];
rz(0.52014822) q[1];
sx q[1];
rz(-1.2112434) q[1];
sx q[1];
rz(-0.21336666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8321229) q[0];
sx q[0];
rz(-1.42064) q[0];
sx q[0];
rz(-3.0792159) q[0];
x q[1];
rz(-2.2198898) q[2];
sx q[2];
rz(-1.9919167) q[2];
sx q[2];
rz(-2.365961) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6345916) q[1];
sx q[1];
rz(-1.4774972) q[1];
sx q[1];
rz(-0.044528151) q[1];
rz(-1.1150241) q[3];
sx q[3];
rz(-0.57330647) q[3];
sx q[3];
rz(-2.440883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0352036) q[2];
sx q[2];
rz(-0.28855244) q[2];
sx q[2];
rz(2.8698548) q[2];
rz(-0.64618293) q[3];
sx q[3];
rz(-1.313442) q[3];
sx q[3];
rz(1.9338231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3216517) q[0];
sx q[0];
rz(-0.8701179) q[0];
sx q[0];
rz(0.2680378) q[0];
rz(-2.9053814) q[1];
sx q[1];
rz(-1.3583438) q[1];
sx q[1];
rz(-1.4150298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7572896) q[0];
sx q[0];
rz(-1.4688604) q[0];
sx q[0];
rz(1.524964) q[0];
x q[1];
rz(2.7162464) q[2];
sx q[2];
rz(-2.7470664) q[2];
sx q[2];
rz(2.2262642) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8658262) q[1];
sx q[1];
rz(-1.8477461) q[1];
sx q[1];
rz(-0.45487483) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2938817) q[3];
sx q[3];
rz(-0.91654449) q[3];
sx q[3];
rz(-0.89756075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1872824) q[2];
sx q[2];
rz(-0.3287181) q[2];
sx q[2];
rz(1.7598565) q[2];
rz(-2.7749744) q[3];
sx q[3];
rz(-1.6691875) q[3];
sx q[3];
rz(3.0439175) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5919375) q[0];
sx q[0];
rz(-0.20579919) q[0];
sx q[0];
rz(-0.016121443) q[0];
rz(-1.4664117) q[1];
sx q[1];
rz(-1.6203531) q[1];
sx q[1];
rz(-2.256567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67041677) q[0];
sx q[0];
rz(-2.6993641) q[0];
sx q[0];
rz(-2.7294374) q[0];
rz(-1.9375291) q[2];
sx q[2];
rz(-1.7842494) q[2];
sx q[2];
rz(0.17630028) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.72674561) q[1];
sx q[1];
rz(-1.9512842) q[1];
sx q[1];
rz(-2.6443015) q[1];
x q[2];
rz(-2.7376369) q[3];
sx q[3];
rz(-2.0519629) q[3];
sx q[3];
rz(1.4766525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5625988) q[2];
sx q[2];
rz(-0.4563953) q[2];
sx q[2];
rz(1.5835416) q[2];
rz(1.3715749) q[3];
sx q[3];
rz(-1.8669502) q[3];
sx q[3];
rz(-0.19545999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9850995) q[0];
sx q[0];
rz(-1.766196) q[0];
sx q[0];
rz(0.14955713) q[0];
rz(-1.043148) q[1];
sx q[1];
rz(-0.97277343) q[1];
sx q[1];
rz(-2.3057888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0319894) q[0];
sx q[0];
rz(-1.2471203) q[0];
sx q[0];
rz(-2.8780185) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2611375) q[2];
sx q[2];
rz(-1.5635075) q[2];
sx q[2];
rz(1.4237504) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0975102) q[1];
sx q[1];
rz(-2.9335576) q[1];
sx q[1];
rz(-1.5745622) q[1];
x q[2];
rz(-0.18330611) q[3];
sx q[3];
rz(-1.8170895) q[3];
sx q[3];
rz(-0.016427996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.122637) q[2];
sx q[2];
rz(-1.476172) q[2];
sx q[2];
rz(-0.41610757) q[2];
rz(0.082402669) q[3];
sx q[3];
rz(-0.96937886) q[3];
sx q[3];
rz(0.20127067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0788954) q[0];
sx q[0];
rz(-2.1128928) q[0];
sx q[0];
rz(2.0507574) q[0];
rz(-2.003032) q[1];
sx q[1];
rz(-1.9416315) q[1];
sx q[1];
rz(2.535215) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67865506) q[0];
sx q[0];
rz(-0.69742862) q[0];
sx q[0];
rz(0.3987958) q[0];
rz(-pi) q[1];
rz(-2.3396083) q[2];
sx q[2];
rz(-1.9030266) q[2];
sx q[2];
rz(1.3114245) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6003638) q[1];
sx q[1];
rz(-2.0906587) q[1];
sx q[1];
rz(-0.71345774) q[1];
x q[2];
rz(-0.13952556) q[3];
sx q[3];
rz(-2.4184113) q[3];
sx q[3];
rz(2.4833402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.828317) q[2];
sx q[2];
rz(-0.54855359) q[2];
sx q[2];
rz(-3.0307148) q[2];
rz(-1.5634792) q[3];
sx q[3];
rz(-2.9302247) q[3];
sx q[3];
rz(1.8160688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9368847) q[0];
sx q[0];
rz(-0.81975833) q[0];
sx q[0];
rz(-0.65716609) q[0];
rz(1.3232629) q[1];
sx q[1];
rz(-0.75138775) q[1];
sx q[1];
rz(-2.0164067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60356451) q[0];
sx q[0];
rz(-1.5865933) q[0];
sx q[0];
rz(1.5558916) q[0];
rz(-1.6680587) q[2];
sx q[2];
rz(-1.2583557) q[2];
sx q[2];
rz(-0.51431235) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4954056) q[1];
sx q[1];
rz(-2.0584724) q[1];
sx q[1];
rz(-2.0085006) q[1];
x q[2];
rz(2.4666305) q[3];
sx q[3];
rz(-1.8139903) q[3];
sx q[3];
rz(-1.4362818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40921673) q[2];
sx q[2];
rz(-1.9595307) q[2];
sx q[2];
rz(-2.1803975) q[2];
rz(0.0088648908) q[3];
sx q[3];
rz(-2.2564087) q[3];
sx q[3];
rz(1.9302906) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081414118) q[0];
sx q[0];
rz(-0.1135122) q[0];
sx q[0];
rz(-2.7139582) q[0];
rz(-1.112452) q[1];
sx q[1];
rz(-1.3149657) q[1];
sx q[1];
rz(-2.1449259) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8450981) q[0];
sx q[0];
rz(-2.0581492) q[0];
sx q[0];
rz(0.71655207) q[0];
rz(-1.2995424) q[2];
sx q[2];
rz(-0.82640663) q[2];
sx q[2];
rz(-2.8959993) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8045051) q[1];
sx q[1];
rz(-1.305991) q[1];
sx q[1];
rz(-0.9881641) q[1];
rz(0.57076591) q[3];
sx q[3];
rz(-1.3955355) q[3];
sx q[3];
rz(-1.4890763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6796278) q[2];
sx q[2];
rz(-2.6042016) q[2];
sx q[2];
rz(-2.4930387) q[2];
rz(-0.31351659) q[3];
sx q[3];
rz(-2.2532012) q[3];
sx q[3];
rz(3.0891109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.736883) q[0];
sx q[0];
rz(-2.9209904) q[0];
sx q[0];
rz(-2.1746461) q[0];
rz(-2.1838358) q[1];
sx q[1];
rz(-1.0916595) q[1];
sx q[1];
rz(-1.0831833) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5388018) q[0];
sx q[0];
rz(-1.1447009) q[0];
sx q[0];
rz(-1.2998796) q[0];
rz(-pi) q[1];
rz(0.73019256) q[2];
sx q[2];
rz(-2.0959404) q[2];
sx q[2];
rz(2.5333135) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5450648) q[1];
sx q[1];
rz(-1.178243) q[1];
sx q[1];
rz(-2.3190581) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5708709) q[3];
sx q[3];
rz(-1.0914472) q[3];
sx q[3];
rz(-2.0572995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0208685) q[2];
sx q[2];
rz(-0.82896295) q[2];
sx q[2];
rz(-1.6797569) q[2];
rz(-2.7802137) q[3];
sx q[3];
rz(-1.8166108) q[3];
sx q[3];
rz(-0.98226515) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7892889) q[0];
sx q[0];
rz(-2.4109349) q[0];
sx q[0];
rz(2.5590382) q[0];
rz(1.7763058) q[1];
sx q[1];
rz(-2.1047695) q[1];
sx q[1];
rz(0.22886151) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6270646) q[0];
sx q[0];
rz(-0.43958966) q[0];
sx q[0];
rz(2.702781) q[0];
rz(-pi) q[1];
rz(0.9791012) q[2];
sx q[2];
rz(-0.31394638) q[2];
sx q[2];
rz(0.57127956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6534253) q[1];
sx q[1];
rz(-1.5060802) q[1];
sx q[1];
rz(-1.72987) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39005847) q[3];
sx q[3];
rz(-1.164468) q[3];
sx q[3];
rz(1.4618946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.87054306) q[2];
sx q[2];
rz(-2.5898263) q[2];
sx q[2];
rz(1.40847) q[2];
rz(-0.72077858) q[3];
sx q[3];
rz(-1.9920789) q[3];
sx q[3];
rz(1.6225947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0558753) q[0];
sx q[0];
rz(-1.4376823) q[0];
sx q[0];
rz(1.8327376) q[0];
rz(1.5695288) q[1];
sx q[1];
rz(-2.5253898) q[1];
sx q[1];
rz(0.22402221) q[1];
rz(-1.189255) q[2];
sx q[2];
rz(-0.59902643) q[2];
sx q[2];
rz(0.35120102) q[2];
rz(-0.15275501) q[3];
sx q[3];
rz(-1.7293617) q[3];
sx q[3];
rz(0.21241906) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
