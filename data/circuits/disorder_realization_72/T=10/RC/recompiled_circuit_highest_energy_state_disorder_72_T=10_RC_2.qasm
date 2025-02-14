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
rz(-1.8347972) q[0];
sx q[0];
rz(-2.2936294) q[0];
sx q[0];
rz(0.037394878) q[0];
rz(-1.0132064) q[1];
sx q[1];
rz(-0.4816882) q[1];
sx q[1];
rz(1.4451292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28222154) q[0];
sx q[0];
rz(-1.9992775) q[0];
sx q[0];
rz(1.3708049) q[0];
rz(-pi) q[1];
rz(3.0710241) q[2];
sx q[2];
rz(-1.6733992) q[2];
sx q[2];
rz(3.0377394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7788275) q[1];
sx q[1];
rz(-2.3794075) q[1];
sx q[1];
rz(-2.306421) q[1];
rz(-pi) q[2];
rz(2.1358498) q[3];
sx q[3];
rz(-1.4742719) q[3];
sx q[3];
rz(-0.18894503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.41363132) q[2];
sx q[2];
rz(-3.0484338) q[2];
sx q[2];
rz(1.5067345) q[2];
rz(2.9480751) q[3];
sx q[3];
rz(-0.7842803) q[3];
sx q[3];
rz(0.94582742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0993318) q[0];
sx q[0];
rz(-1.7758545) q[0];
sx q[0];
rz(0.17901626) q[0];
rz(1.6049113) q[1];
sx q[1];
rz(-0.34114006) q[1];
sx q[1];
rz(1.5515597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66288469) q[0];
sx q[0];
rz(-2.835412) q[0];
sx q[0];
rz(-1.2423489) q[0];
rz(0.70945246) q[2];
sx q[2];
rz(-1.7691028) q[2];
sx q[2];
rz(2.545216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1074321) q[1];
sx q[1];
rz(-1.1723659) q[1];
sx q[1];
rz(-2.855019) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1158982) q[3];
sx q[3];
rz(-1.1602931) q[3];
sx q[3];
rz(0.59343597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67148525) q[2];
sx q[2];
rz(-1.476373) q[2];
sx q[2];
rz(3.1191077) q[2];
rz(0.89618987) q[3];
sx q[3];
rz(-0.78138566) q[3];
sx q[3];
rz(1.6270858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80980587) q[0];
sx q[0];
rz(-2.9349194) q[0];
sx q[0];
rz(1.608954) q[0];
rz(2.0017852) q[1];
sx q[1];
rz(-0.31875691) q[1];
sx q[1];
rz(-2.2679813) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5748567) q[0];
sx q[0];
rz(-2.2169161) q[0];
sx q[0];
rz(-2.2227051) q[0];
rz(-pi) q[1];
rz(0.64675764) q[2];
sx q[2];
rz(-1.1144179) q[2];
sx q[2];
rz(2.6025085) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1908299) q[1];
sx q[1];
rz(-2.0520981) q[1];
sx q[1];
rz(1.7442622) q[1];
x q[2];
rz(-0.26362147) q[3];
sx q[3];
rz(-2.859451) q[3];
sx q[3];
rz(-2.1745115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5736299) q[2];
sx q[2];
rz(-0.45845389) q[2];
sx q[2];
rz(-2.6652794) q[2];
rz(-1.025398) q[3];
sx q[3];
rz(-1.3566596) q[3];
sx q[3];
rz(-0.29388139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3697701) q[0];
sx q[0];
rz(-2.2875468) q[0];
sx q[0];
rz(0.59637946) q[0];
rz(0.75309938) q[1];
sx q[1];
rz(-1.785708) q[1];
sx q[1];
rz(0.3515884) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.387334) q[0];
sx q[0];
rz(-1.046858) q[0];
sx q[0];
rz(-1.7225527) q[0];
rz(0.99927665) q[2];
sx q[2];
rz(-0.7692996) q[2];
sx q[2];
rz(2.3623449) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5758613) q[1];
sx q[1];
rz(-1.8003776) q[1];
sx q[1];
rz(-1.1665535) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91686317) q[3];
sx q[3];
rz(-1.0427949) q[3];
sx q[3];
rz(2.7478086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4190462) q[2];
sx q[2];
rz(-1.498035) q[2];
sx q[2];
rz(-2.2124115) q[2];
rz(1.2292713) q[3];
sx q[3];
rz(-0.036673948) q[3];
sx q[3];
rz(1.2693955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1790328) q[0];
sx q[0];
rz(-2.5293009) q[0];
sx q[0];
rz(0.52781934) q[0];
rz(0.57012308) q[1];
sx q[1];
rz(-0.56990439) q[1];
sx q[1];
rz(-2.747587) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.67585) q[0];
sx q[0];
rz(-0.52211863) q[0];
sx q[0];
rz(0.39710775) q[0];
x q[1];
rz(-1.9377685) q[2];
sx q[2];
rz(-1.0093371) q[2];
sx q[2];
rz(0.16753627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7841508) q[1];
sx q[1];
rz(-1.3489745) q[1];
sx q[1];
rz(3.0206969) q[1];
rz(2.3218003) q[3];
sx q[3];
rz(-2.2283883) q[3];
sx q[3];
rz(2.4321604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3923308) q[2];
sx q[2];
rz(-1.8529842) q[2];
sx q[2];
rz(2.0070845) q[2];
rz(0.025731651) q[3];
sx q[3];
rz(-2.0870356) q[3];
sx q[3];
rz(0.64190763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5662017) q[0];
sx q[0];
rz(-1.8113149) q[0];
sx q[0];
rz(-2.6824685) q[0];
rz(0.8210012) q[1];
sx q[1];
rz(-0.88290015) q[1];
sx q[1];
rz(-2.6301036) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59508649) q[0];
sx q[0];
rz(-1.0985785) q[0];
sx q[0];
rz(-2.8123463) q[0];
x q[1];
rz(1.8051926) q[2];
sx q[2];
rz(-0.92051855) q[2];
sx q[2];
rz(2.5087207) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2177812) q[1];
sx q[1];
rz(-2.1379323) q[1];
sx q[1];
rz(1.7497211) q[1];
x q[2];
rz(1.4019764) q[3];
sx q[3];
rz(-2.4509015) q[3];
sx q[3];
rz(0.21571479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56200999) q[2];
sx q[2];
rz(-2.0333813) q[2];
sx q[2];
rz(-0.75135922) q[2];
rz(2.5747418) q[3];
sx q[3];
rz(-1.6040498) q[3];
sx q[3];
rz(-1.8142627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.421627) q[0];
sx q[0];
rz(-0.1816853) q[0];
sx q[0];
rz(-0.0075465329) q[0];
rz(-0.94002062) q[1];
sx q[1];
rz(-2.4766141) q[1];
sx q[1];
rz(0.13253458) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7081817) q[0];
sx q[0];
rz(-1.1601747) q[0];
sx q[0];
rz(-2.2020086) q[0];
x q[1];
rz(0.33905115) q[2];
sx q[2];
rz(-2.7586852) q[2];
sx q[2];
rz(2.6660181) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30250994) q[1];
sx q[1];
rz(-1.1801475) q[1];
sx q[1];
rz(1.7724228) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0293983) q[3];
sx q[3];
rz(-2.762714) q[3];
sx q[3];
rz(2.1357812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.14564766) q[2];
sx q[2];
rz(-1.4952679) q[2];
sx q[2];
rz(-3.0403467) q[2];
rz(-0.64949399) q[3];
sx q[3];
rz(-2.5443304) q[3];
sx q[3];
rz(-0.35972843) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7262909) q[0];
sx q[0];
rz(-1.8446209) q[0];
sx q[0];
rz(1.8743961) q[0];
rz(-1.2665117) q[1];
sx q[1];
rz(-0.15788618) q[1];
sx q[1];
rz(0.74323851) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948461) q[0];
sx q[0];
rz(-2.136793) q[0];
sx q[0];
rz(1.8106242) q[0];
rz(-pi) q[1];
rz(1.2879804) q[2];
sx q[2];
rz(-2.3415945) q[2];
sx q[2];
rz(1.9375506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91620244) q[1];
sx q[1];
rz(-1.7280726) q[1];
sx q[1];
rz(1.8099524) q[1];
rz(-pi) q[2];
x q[2];
rz(2.322213) q[3];
sx q[3];
rz(-2.1158822) q[3];
sx q[3];
rz(1.8784639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1242975) q[2];
sx q[2];
rz(-0.20169078) q[2];
sx q[2];
rz(-1.7151493) q[2];
rz(2.1258449) q[3];
sx q[3];
rz(-1.7044715) q[3];
sx q[3];
rz(-0.82971853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5947241) q[0];
sx q[0];
rz(-2.4877553) q[0];
sx q[0];
rz(0.015901707) q[0];
rz(-1.5025899) q[1];
sx q[1];
rz(-2.465261) q[1];
sx q[1];
rz(2.1471088) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54399207) q[0];
sx q[0];
rz(-2.3975537) q[0];
sx q[0];
rz(2.1617894) q[0];
rz(1.5076958) q[2];
sx q[2];
rz(-1.1931538) q[2];
sx q[2];
rz(-1.3444195) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2331824) q[1];
sx q[1];
rz(-0.841978) q[1];
sx q[1];
rz(0.89836095) q[1];
rz(-pi) q[2];
rz(-2.9904305) q[3];
sx q[3];
rz(-2.0796607) q[3];
sx q[3];
rz(0.14948949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9859621) q[2];
sx q[2];
rz(-1.4302284) q[2];
sx q[2];
rz(2.963781) q[2];
rz(-1.6832247) q[3];
sx q[3];
rz(-0.42176133) q[3];
sx q[3];
rz(1.2678857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47852248) q[0];
sx q[0];
rz(-1.2489742) q[0];
sx q[0];
rz(-1.857969) q[0];
rz(-0.78370699) q[1];
sx q[1];
rz(-0.77373928) q[1];
sx q[1];
rz(-1.3073889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9147025) q[0];
sx q[0];
rz(-0.44825867) q[0];
sx q[0];
rz(-2.9284555) q[0];
x q[1];
rz(-1.9587014) q[2];
sx q[2];
rz(-1.5989306) q[2];
sx q[2];
rz(-1.494734) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5518903) q[1];
sx q[1];
rz(-1.3968727) q[1];
sx q[1];
rz(-3.0822192) q[1];
rz(-pi) q[2];
rz(-0.41475716) q[3];
sx q[3];
rz(-2.0817753) q[3];
sx q[3];
rz(1.8243044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6002097) q[2];
sx q[2];
rz(-1.4098488) q[2];
sx q[2];
rz(0.48995885) q[2];
rz(-1.1146891) q[3];
sx q[3];
rz(-1.1763108) q[3];
sx q[3];
rz(2.8118242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212696) q[0];
sx q[0];
rz(-1.9135495) q[0];
sx q[0];
rz(1.3187153) q[0];
rz(-1.8439138) q[1];
sx q[1];
rz(-2.4083125) q[1];
sx q[1];
rz(-0.42793035) q[1];
rz(1.3135459) q[2];
sx q[2];
rz(-1.8218818) q[2];
sx q[2];
rz(2.8043361) q[2];
rz(-0.86434563) q[3];
sx q[3];
rz(-1.4180687) q[3];
sx q[3];
rz(-1.1579469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
