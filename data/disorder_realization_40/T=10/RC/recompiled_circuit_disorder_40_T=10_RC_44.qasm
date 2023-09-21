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
rz(-2.9878374) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0271887) q[0];
sx q[0];
rz(-1.1429042) q[0];
sx q[0];
rz(-1.2501636) q[0];
rz(-pi) q[1];
rz(1.2674238) q[2];
sx q[2];
rz(-0.25250013) q[2];
sx q[2];
rz(1.4361824) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7003324) q[1];
sx q[1];
rz(-1.2520737) q[1];
sx q[1];
rz(3.1096427) q[1];
rz(2.2114803) q[3];
sx q[3];
rz(-0.53033391) q[3];
sx q[3];
rz(1.9213284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6979606) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(-1.704818) q[2];
rz(2.4076961) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2519418) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(-2.2170128) q[0];
rz(-2.1444767) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(1.8181713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01100563) q[0];
sx q[0];
rz(-0.6479833) q[0];
sx q[0];
rz(2.381071) q[0];
rz(-pi) q[1];
rz(-1.6152482) q[2];
sx q[2];
rz(-1.2215081) q[2];
sx q[2];
rz(0.312422) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18113187) q[1];
sx q[1];
rz(-2.301553) q[1];
sx q[1];
rz(-1.204927) q[1];
x q[2];
rz(-1.0074535) q[3];
sx q[3];
rz(-0.65348071) q[3];
sx q[3];
rz(0.37936488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9374342) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(-2.3834174) q[2];
rz(-2.5126863) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(-1.1531856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.4085098) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(0.68840233) q[0];
rz(-0.06772659) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(-2.6115131) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63307525) q[0];
sx q[0];
rz(-1.2808262) q[0];
sx q[0];
rz(3.0301208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7128503) q[2];
sx q[2];
rz(-2.0566166) q[2];
sx q[2];
rz(0.57810099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2391795) q[1];
sx q[1];
rz(-1.2352408) q[1];
sx q[1];
rz(2.2526342) q[1];
rz(-pi) q[2];
rz(0.097966627) q[3];
sx q[3];
rz(-1.3324105) q[3];
sx q[3];
rz(1.0007678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3337341) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(-2.2401436) q[2];
rz(-0.83550134) q[3];
sx q[3];
rz(-1.6136026) q[3];
sx q[3];
rz(-1.3114595) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19105844) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(-2.3572671) q[0];
rz(-3.0803608) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(3.004946) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2480859) q[0];
sx q[0];
rz(-0.49960217) q[0];
sx q[0];
rz(1.2953555) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.253445) q[2];
sx q[2];
rz(-1.1270521) q[2];
sx q[2];
rz(1.2197942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.54388753) q[1];
sx q[1];
rz(-1.502431) q[1];
sx q[1];
rz(3.1176561) q[1];
x q[2];
rz(1.4284381) q[3];
sx q[3];
rz(-2.5609303) q[3];
sx q[3];
rz(-0.84238392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0121997) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(0.56048918) q[2];
rz(-0.012332049) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(-2.0509317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.2416302) q[1];
sx q[1];
rz(-1.7339773) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53769892) q[0];
sx q[0];
rz(-1.1802117) q[0];
sx q[0];
rz(-1.1917398) q[0];
x q[1];
rz(1.1248963) q[2];
sx q[2];
rz(-2.5113528) q[2];
sx q[2];
rz(1.5843887) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.40286139) q[1];
sx q[1];
rz(-2.1364473) q[1];
sx q[1];
rz(-2.9823142) q[1];
rz(-1.2110662) q[3];
sx q[3];
rz(-1.3621646) q[3];
sx q[3];
rz(1.2071351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.451482) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(3.0991128) q[2];
rz(-0.63043198) q[3];
sx q[3];
rz(-2.5094331) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(0.4831627) q[0];
rz(1.0522316) q[1];
sx q[1];
rz(-1.1455043) q[1];
sx q[1];
rz(0.56484708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72950596) q[0];
sx q[0];
rz(-0.54486638) q[0];
sx q[0];
rz(-1.3077523) q[0];
rz(-1.1499314) q[2];
sx q[2];
rz(-1.6399584) q[2];
sx q[2];
rz(1.7771306) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55826742) q[1];
sx q[1];
rz(-0.43154432) q[1];
sx q[1];
rz(1.4285018) q[1];
rz(0.58638339) q[3];
sx q[3];
rz(-2.8660503) q[3];
sx q[3];
rz(1.9540209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.57006449) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(-1.338039) q[2];
rz(-1.8367052) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(0.6750955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9682482) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(-2.5937953) q[0];
rz(-0.785218) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(-0.26842591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7086605) q[0];
sx q[0];
rz(-1.4888568) q[0];
sx q[0];
rz(-0.24631484) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3586876) q[2];
sx q[2];
rz(-0.9559388) q[2];
sx q[2];
rz(-2.8564786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5749579) q[1];
sx q[1];
rz(-0.98595) q[1];
sx q[1];
rz(-3.023874) q[1];
x q[2];
rz(0.3018474) q[3];
sx q[3];
rz(-1.0511304) q[3];
sx q[3];
rz(2.2083851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5238374) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(-0.37627775) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5557264) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(2.1222173) q[0];
rz(0.85340071) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(-0.44874915) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73483) q[0];
sx q[0];
rz(-1.8039628) q[0];
sx q[0];
rz(-1.083311) q[0];
rz(-pi) q[1];
rz(-0.86782311) q[2];
sx q[2];
rz(-2.7542369) q[2];
sx q[2];
rz(-0.58434904) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9948506) q[1];
sx q[1];
rz(-0.59616201) q[1];
sx q[1];
rz(-0.46390987) q[1];
x q[2];
rz(0.7738503) q[3];
sx q[3];
rz(-2.1132831) q[3];
sx q[3];
rz(-1.098793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86924187) q[2];
sx q[2];
rz(-1.7607471) q[2];
sx q[2];
rz(-0.31420079) q[2];
rz(-0.82434404) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(-2.3468988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070351275) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(1.2605793) q[0];
rz(0.64385995) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(3.1226645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2429598) q[0];
sx q[0];
rz(-1.8343933) q[0];
sx q[0];
rz(-1.5522869) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5817714) q[2];
sx q[2];
rz(-2.9580742) q[2];
sx q[2];
rz(0.57932094) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1259809) q[1];
sx q[1];
rz(-2.9503082) q[1];
sx q[1];
rz(0.26856883) q[1];
rz(-pi) q[2];
rz(0.67846672) q[3];
sx q[3];
rz(-2.7760091) q[3];
sx q[3];
rz(-1.4640704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5921322) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(-2.4712759) q[2];
rz(-0.51236764) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(-1.4204773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.55384127) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(0.49945369) q[0];
rz(-1.5669426) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(-1.0826899) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4944045) q[0];
sx q[0];
rz(-1.5218966) q[0];
sx q[0];
rz(-2.5542269) q[0];
rz(-pi) q[1];
rz(3.1084836) q[2];
sx q[2];
rz(-0.81956714) q[2];
sx q[2];
rz(2.5493252) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.618256) q[1];
sx q[1];
rz(-0.54392951) q[1];
sx q[1];
rz(-0.036332794) q[1];
rz(-0.70268481) q[3];
sx q[3];
rz(-2.6583238) q[3];
sx q[3];
rz(1.0815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7426976) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(-2.6297074) q[2];
rz(-0.39294696) q[3];
sx q[3];
rz(-1.4018551) q[3];
sx q[3];
rz(-1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95505161) q[0];
sx q[0];
rz(-1.825009) q[0];
sx q[0];
rz(0.64074989) q[0];
rz(0.74116771) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(2.7411186) q[2];
sx q[2];
rz(-2.1486712) q[2];
sx q[2];
rz(1.490544) q[2];
rz(-2.3951204) q[3];
sx q[3];
rz(-1.4440047) q[3];
sx q[3];
rz(0.11228893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
