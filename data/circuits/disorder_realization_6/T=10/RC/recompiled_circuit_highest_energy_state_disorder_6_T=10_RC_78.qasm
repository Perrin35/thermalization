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
rz(-2.1498635) q[0];
sx q[0];
rz(-2.5340134) q[0];
sx q[0];
rz(2.1460331) q[0];
rz(-2.6140656) q[1];
sx q[1];
rz(-1.0310643) q[1];
sx q[1];
rz(-0.48479015) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82370347) q[0];
sx q[0];
rz(-1.7021588) q[0];
sx q[0];
rz(-2.0844368) q[0];
rz(-pi) q[1];
rz(0.77930696) q[2];
sx q[2];
rz(-1.2504842) q[2];
sx q[2];
rz(-0.28398289) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0409702) q[1];
sx q[1];
rz(-0.61196152) q[1];
sx q[1];
rz(-1.4684567) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4977341) q[3];
sx q[3];
rz(-1.2803852) q[3];
sx q[3];
rz(-2.6892457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.63554865) q[2];
sx q[2];
rz(-2.1940993) q[2];
sx q[2];
rz(0.49723899) q[2];
rz(0.55073589) q[3];
sx q[3];
rz(-0.68998718) q[3];
sx q[3];
rz(-2.0042888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15744844) q[0];
sx q[0];
rz(-1.4009615) q[0];
sx q[0];
rz(0.76341158) q[0];
rz(-2.5570671) q[1];
sx q[1];
rz(-1.3817363) q[1];
sx q[1];
rz(-2.1302628) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4437067) q[0];
sx q[0];
rz(-2.9666565) q[0];
sx q[0];
rz(2.2921497) q[0];
x q[1];
rz(-0.23638921) q[2];
sx q[2];
rz(-1.4589785) q[2];
sx q[2];
rz(-0.37416247) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3608777) q[1];
sx q[1];
rz(-1.3643801) q[1];
sx q[1];
rz(1.2206627) q[1];
rz(-pi) q[2];
x q[2];
rz(2.743558) q[3];
sx q[3];
rz(-1.0938489) q[3];
sx q[3];
rz(1.8534891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4804907) q[2];
sx q[2];
rz(-1.2031518) q[2];
sx q[2];
rz(-1.0002381) q[2];
rz(-0.58146042) q[3];
sx q[3];
rz(-2.4088819) q[3];
sx q[3];
rz(1.940041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30705273) q[0];
sx q[0];
rz(-2.5426799) q[0];
sx q[0];
rz(0.55759984) q[0];
rz(2.8343976) q[1];
sx q[1];
rz(-0.59978825) q[1];
sx q[1];
rz(-2.6934326) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9014943) q[0];
sx q[0];
rz(-2.0848102) q[0];
sx q[0];
rz(-1.7230524) q[0];
rz(1.0501625) q[2];
sx q[2];
rz(-1.7562281) q[2];
sx q[2];
rz(0.38927573) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.243147) q[1];
sx q[1];
rz(-2.2013651) q[1];
sx q[1];
rz(3.10792) q[1];
rz(2.7072315) q[3];
sx q[3];
rz(-1.5944527) q[3];
sx q[3];
rz(-1.0756186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3192516) q[2];
sx q[2];
rz(-2.4288869) q[2];
sx q[2];
rz(-2.5437497) q[2];
rz(2.6939825) q[3];
sx q[3];
rz(-1.8748583) q[3];
sx q[3];
rz(-3.083057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41281259) q[0];
sx q[0];
rz(-0.030981177) q[0];
sx q[0];
rz(-2.431751) q[0];
rz(-0.9791044) q[1];
sx q[1];
rz(-2.282228) q[1];
sx q[1];
rz(-1.3213347) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0036572) q[0];
sx q[0];
rz(-2.4383579) q[0];
sx q[0];
rz(3.0556995) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3611907) q[2];
sx q[2];
rz(-1.4596887) q[2];
sx q[2];
rz(1.358913) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2817414) q[1];
sx q[1];
rz(-2.4415589) q[1];
sx q[1];
rz(1.9938049) q[1];
x q[2];
rz(2.1141601) q[3];
sx q[3];
rz(-1.8413944) q[3];
sx q[3];
rz(-0.71934127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1506302) q[2];
sx q[2];
rz(-1.7143098) q[2];
sx q[2];
rz(-0.70685753) q[2];
rz(-1.7957211) q[3];
sx q[3];
rz(-1.6790877) q[3];
sx q[3];
rz(-2.4367387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83955806) q[0];
sx q[0];
rz(-0.78080451) q[0];
sx q[0];
rz(1.5020405) q[0];
rz(1.5976277) q[1];
sx q[1];
rz(-2.2815956) q[1];
sx q[1];
rz(1.647321) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2236007) q[0];
sx q[0];
rz(-0.57044464) q[0];
sx q[0];
rz(-1.5738652) q[0];
rz(-1.5635033) q[2];
sx q[2];
rz(-2.1661758) q[2];
sx q[2];
rz(-0.0059758107) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9372796) q[1];
sx q[1];
rz(-1.8674074) q[1];
sx q[1];
rz(-1.2526172) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95559883) q[3];
sx q[3];
rz(-1.5892059) q[3];
sx q[3];
rz(1.9083202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8813701) q[2];
sx q[2];
rz(-1.8689195) q[2];
sx q[2];
rz(2.5831045) q[2];
rz(0.50885606) q[3];
sx q[3];
rz(-1.6130092) q[3];
sx q[3];
rz(1.2671027) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59638554) q[0];
sx q[0];
rz(-1.243243) q[0];
sx q[0];
rz(-0.16045706) q[0];
rz(0.67333737) q[1];
sx q[1];
rz(-1.9247749) q[1];
sx q[1];
rz(1.1268667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50235275) q[0];
sx q[0];
rz(-1.2099838) q[0];
sx q[0];
rz(2.8176184) q[0];
rz(-1.8931863) q[2];
sx q[2];
rz(-0.83288542) q[2];
sx q[2];
rz(-2.5423342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.018467598) q[1];
sx q[1];
rz(-1.7729802) q[1];
sx q[1];
rz(0.98337652) q[1];
rz(-pi) q[2];
x q[2];
rz(1.91332) q[3];
sx q[3];
rz(-2.5728101) q[3];
sx q[3];
rz(2.3186213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0141853) q[2];
sx q[2];
rz(-1.5994453) q[2];
sx q[2];
rz(-0.7737774) q[2];
rz(-0.989178) q[3];
sx q[3];
rz(-2.7343605) q[3];
sx q[3];
rz(-1.0729084) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2268552) q[0];
sx q[0];
rz(-1.6586774) q[0];
sx q[0];
rz(2.3684655) q[0];
rz(1.5559366) q[1];
sx q[1];
rz(-1.2890041) q[1];
sx q[1];
rz(1.9645436) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0814514) q[0];
sx q[0];
rz(-0.8340652) q[0];
sx q[0];
rz(-3.0990776) q[0];
x q[1];
rz(-0.99780166) q[2];
sx q[2];
rz(-1.3062451) q[2];
sx q[2];
rz(-0.71392347) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24304535) q[1];
sx q[1];
rz(-0.70705251) q[1];
sx q[1];
rz(1.9999595) q[1];
x q[2];
rz(0.66399948) q[3];
sx q[3];
rz(-0.2240847) q[3];
sx q[3];
rz(-0.9730207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3561463) q[2];
sx q[2];
rz(-2.7768504) q[2];
sx q[2];
rz(-2.6981603) q[2];
rz(-0.38608471) q[3];
sx q[3];
rz(-1.723879) q[3];
sx q[3];
rz(-2.940787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76483738) q[0];
sx q[0];
rz(-1.5926462) q[0];
sx q[0];
rz(-0.03373294) q[0];
rz(-2.4434166) q[1];
sx q[1];
rz(-0.59711421) q[1];
sx q[1];
rz(2.4765292) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9586104) q[0];
sx q[0];
rz(-1.1232875) q[0];
sx q[0];
rz(0.29194269) q[0];
rz(-pi) q[1];
rz(-0.84455281) q[2];
sx q[2];
rz(-1.9366656) q[2];
sx q[2];
rz(-1.4839811) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5678957) q[1];
sx q[1];
rz(-2.0665281) q[1];
sx q[1];
rz(0.8521518) q[1];
rz(-pi) q[2];
rz(0.31954664) q[3];
sx q[3];
rz(-1.5509774) q[3];
sx q[3];
rz(1.457959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4195121) q[2];
sx q[2];
rz(-2.2245912) q[2];
sx q[2];
rz(1.7783995) q[2];
rz(-0.52115399) q[3];
sx q[3];
rz(-1.8173953) q[3];
sx q[3];
rz(0.49191973) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9687965) q[0];
sx q[0];
rz(-1.7954614) q[0];
sx q[0];
rz(0.81934339) q[0];
rz(-2.4567545) q[1];
sx q[1];
rz(-1.2872773) q[1];
sx q[1];
rz(0.12231621) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93520892) q[0];
sx q[0];
rz(-2.1851086) q[0];
sx q[0];
rz(-2.6198316) q[0];
rz(-pi) q[1];
rz(2.7610131) q[2];
sx q[2];
rz(-2.1316445) q[2];
sx q[2];
rz(0.83319908) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6284991) q[1];
sx q[1];
rz(-2.4446396) q[1];
sx q[1];
rz(-0.51817399) q[1];
rz(-pi) q[2];
rz(0.64494452) q[3];
sx q[3];
rz(-2.44584) q[3];
sx q[3];
rz(2.6464045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43716064) q[2];
sx q[2];
rz(-1.8081343) q[2];
sx q[2];
rz(0.86755794) q[2];
rz(1.7433172) q[3];
sx q[3];
rz(-2.7531392) q[3];
sx q[3];
rz(-0.5790264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45359465) q[0];
sx q[0];
rz(-1.9168251) q[0];
sx q[0];
rz(2.0523742) q[0];
rz(0.047529686) q[1];
sx q[1];
rz(-2.0952974) q[1];
sx q[1];
rz(2.4542123) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7318667) q[0];
sx q[0];
rz(-2.0610524) q[0];
sx q[0];
rz(0.027564647) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1310523) q[2];
sx q[2];
rz(-0.49695914) q[2];
sx q[2];
rz(-0.69253773) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7043651) q[1];
sx q[1];
rz(-0.95672551) q[1];
sx q[1];
rz(1.006587) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93411719) q[3];
sx q[3];
rz(-1.6457594) q[3];
sx q[3];
rz(1.9270169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1462732) q[2];
sx q[2];
rz(-2.5226888) q[2];
sx q[2];
rz(1.6590365) q[2];
rz(0.65685529) q[3];
sx q[3];
rz(-1.992179) q[3];
sx q[3];
rz(2.759554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.252608) q[0];
sx q[0];
rz(-1.9883307) q[0];
sx q[0];
rz(-2.0345732) q[0];
rz(0.57869115) q[1];
sx q[1];
rz(-2.3042669) q[1];
sx q[1];
rz(-0.28490983) q[1];
rz(0.22004057) q[2];
sx q[2];
rz(-2.7773428) q[2];
sx q[2];
rz(-0.11014948) q[2];
rz(-2.9978776) q[3];
sx q[3];
rz(-2.6201709) q[3];
sx q[3];
rz(-1.3377671) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
