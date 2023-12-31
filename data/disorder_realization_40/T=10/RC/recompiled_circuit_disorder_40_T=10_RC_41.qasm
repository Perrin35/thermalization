OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5326795) q[0];
sx q[0];
rz(3.5182313) q[0];
sx q[0];
rz(9.3129899) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(0.15375528) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.114404) q[0];
sx q[0];
rz(-1.9986885) q[0];
sx q[0];
rz(1.2501636) q[0];
rz(-3.0646677) q[2];
sx q[2];
rz(-1.811532) q[2];
sx q[2];
rz(1.392729) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7003324) q[1];
sx q[1];
rz(-1.889519) q[1];
sx q[1];
rz(3.1096427) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93011232) q[3];
sx q[3];
rz(-2.6112587) q[3];
sx q[3];
rz(-1.9213284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6979606) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(1.704818) q[2];
rz(-2.4076961) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(0.51600391) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88965082) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(2.2170128) q[0];
rz(2.1444767) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(-1.8181713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01100563) q[0];
sx q[0];
rz(-2.4936094) q[0];
sx q[0];
rz(-2.381071) q[0];
rz(1.5263444) q[2];
sx q[2];
rz(-1.2215081) q[2];
sx q[2];
rz(-2.8291707) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4393649) q[1];
sx q[1];
rz(-2.339748) q[1];
sx q[1];
rz(-2.761809) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7534361) q[3];
sx q[3];
rz(-1.0309439) q[3];
sx q[3];
rz(0.29263465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9374342) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(-0.75817529) q[2];
rz(2.5126863) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(1.1531856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73308289) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(-2.4531903) q[0];
rz(0.06772659) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(2.6115131) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2358658) q[0];
sx q[0];
rz(-1.4639963) q[0];
sx q[0];
rz(1.279116) q[0];
x q[1];
rz(0.42874239) q[2];
sx q[2];
rz(-2.0566166) q[2];
sx q[2];
rz(-0.57810099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59280076) q[1];
sx q[1];
rz(-2.2081516) q[1];
sx q[1];
rz(-2.7194276) q[1];
rz(-pi) q[2];
rz(1.9534555) q[3];
sx q[3];
rz(-2.8842162) q[3];
sx q[3];
rz(-2.5352258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3337341) q[2];
sx q[2];
rz(-3.1299751) q[2];
sx q[2];
rz(-2.2401436) q[2];
rz(-0.83550134) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(1.3114595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19105844) q[0];
sx q[0];
rz(-1.5427417) q[0];
sx q[0];
rz(2.3572671) q[0];
rz(-0.061231881) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(3.004946) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2480859) q[0];
sx q[0];
rz(-2.6419905) q[0];
sx q[0];
rz(1.2953555) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6776671) q[2];
sx q[2];
rz(-1.2850962) q[2];
sx q[2];
rz(2.6505016) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9347958) q[1];
sx q[1];
rz(-3.0691642) q[1];
sx q[1];
rz(-1.9070684) q[1];
rz(-pi) q[2];
rz(2.1468048) q[3];
sx q[3];
rz(-1.4928865) q[3];
sx q[3];
rz(0.60914492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1293929) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(-2.5811035) q[2];
rz(-0.012332049) q[3];
sx q[3];
rz(-2.2380232) q[3];
sx q[3];
rz(-1.0906609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(2.7085768) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(-0.82114712) q[0];
rz(0.87617809) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(1.7339773) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6038937) q[0];
sx q[0];
rz(-1.1802117) q[0];
sx q[0];
rz(-1.9498528) q[0];
x q[1];
rz(1.1248963) q[2];
sx q[2];
rz(-2.5113528) q[2];
sx q[2];
rz(1.5843887) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.447532) q[1];
sx q[1];
rz(-0.58528712) q[1];
sx q[1];
rz(-1.8156169) q[1];
x q[2];
rz(1.9305265) q[3];
sx q[3];
rz(-1.779428) q[3];
sx q[3];
rz(-1.2071351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6901107) q[2];
sx q[2];
rz(-1.9268945) q[2];
sx q[2];
rz(-3.0991128) q[2];
rz(0.63043198) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(-0.4831627) q[0];
rz(-1.0522316) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(-2.5767456) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171377) q[0];
sx q[0];
rz(-2.0949445) q[0];
sx q[0];
rz(-0.15630396) q[0];
rz(3.0658423) q[2];
sx q[2];
rz(-1.1510013) q[2];
sx q[2];
rz(-0.23725739) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5833252) q[1];
sx q[1];
rz(-0.43154432) q[1];
sx q[1];
rz(-1.4285018) q[1];
rz(-pi) q[2];
rz(0.58638339) q[3];
sx q[3];
rz(-2.8660503) q[3];
sx q[3];
rz(-1.1875718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.57006449) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(-1.338039) q[2];
rz(1.3048874) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(0.6750955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17334443) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(0.54779732) q[0];
rz(-0.785218) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(2.8731667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9831532) q[0];
sx q[0];
rz(-1.8162677) q[0];
sx q[0];
rz(1.4863187) q[0];
rz(-0.78290501) q[2];
sx q[2];
rz(-2.1856538) q[2];
sx q[2];
rz(0.2851141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.35560289) q[1];
sx q[1];
rz(-0.59521788) q[1];
sx q[1];
rz(1.7463513) q[1];
rz(-1.0915756) q[3];
sx q[3];
rz(-2.5476533) q[3];
sx q[3];
rz(1.6483496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5238374) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(-2.3274373) q[2];
rz(2.7653149) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(-3.0686839) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5557264) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(-1.0193753) q[0];
rz(-0.85340071) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(0.44874915) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42496029) q[0];
sx q[0];
rz(-0.53629959) q[0];
sx q[0];
rz(2.0400356) q[0];
x q[1];
rz(1.2690527) q[2];
sx q[2];
rz(-1.8174968) q[2];
sx q[2];
rz(-1.6517284) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6905578) q[1];
sx q[1];
rz(-1.0447377) q[1];
sx q[1];
rz(1.2760389) q[1];
rz(-2.429871) q[3];
sx q[3];
rz(-2.2300643) q[3];
sx q[3];
rz(3.1275415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86924187) q[2];
sx q[2];
rz(-1.7607471) q[2];
sx q[2];
rz(-0.31420079) q[2];
rz(2.3172486) q[3];
sx q[3];
rz(-2.6894675) q[3];
sx q[3];
rz(2.3468988) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0712414) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(-1.8810133) q[0];
rz(2.4977327) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(-0.018928122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2429598) q[0];
sx q[0];
rz(-1.3071994) q[0];
sx q[0];
rz(1.5522869) q[0];
rz(-pi) q[1];
x q[1];
rz(1.47255) q[2];
sx q[2];
rz(-1.4155404) q[2];
sx q[2];
rz(1.9948024) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2889274) q[1];
sx q[1];
rz(-1.7551433) q[1];
sx q[1];
rz(1.6221371) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4631259) q[3];
sx q[3];
rz(-0.36558357) q[3];
sx q[3];
rz(1.4640704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5921322) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(2.4712759) q[2];
rz(0.51236764) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(-1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55384127) q[0];
sx q[0];
rz(-1.2441664) q[0];
sx q[0];
rz(-2.642139) q[0];
rz(-1.5746501) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(-1.0826899) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854541) q[0];
sx q[0];
rz(-0.98422613) q[0];
sx q[0];
rz(-1.6295208) q[0];
x q[1];
rz(0.033109025) q[2];
sx q[2];
rz(-0.81956714) q[2];
sx q[2];
rz(0.59226743) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.078552695) q[1];
sx q[1];
rz(-1.5895956) q[1];
sx q[1];
rz(0.54363721) q[1];
rz(-pi) q[2];
rz(-0.70268481) q[3];
sx q[3];
rz(-2.6583238) q[3];
sx q[3];
rz(1.0815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3988951) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(-2.6297074) q[2];
rz(-0.39294696) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.186541) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(-0.74116771) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(-1.0319866) q[2];
sx q[2];
rz(-2.451755) q[2];
sx q[2];
rz(-2.3103726) q[2];
rz(2.9560271) q[3];
sx q[3];
rz(-2.3864828) q[3];
sx q[3];
rz(-1.3226487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
