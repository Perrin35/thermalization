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
rz(1.2020741) q[0];
sx q[0];
rz(-2.095686) q[0];
sx q[0];
rz(2.3910971) q[0];
rz(-2.9208288) q[1];
sx q[1];
rz(-1.060744) q[1];
sx q[1];
rz(-1.9579252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1738244) q[0];
sx q[0];
rz(-0.90961752) q[0];
sx q[0];
rz(-2.1029766) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92249845) q[2];
sx q[2];
rz(-0.20805173) q[2];
sx q[2];
rz(1.7788943) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.466063) q[1];
sx q[1];
rz(-2.4616082) q[1];
sx q[1];
rz(3.1230981) q[1];
rz(-pi) q[2];
rz(-1.8791844) q[3];
sx q[3];
rz(-1.7636136) q[3];
sx q[3];
rz(-0.28379019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.815328) q[2];
sx q[2];
rz(-1.1417737) q[2];
sx q[2];
rz(0.092078837) q[2];
rz(-2.0103256) q[3];
sx q[3];
rz(-2.0720033) q[3];
sx q[3];
rz(0.85589516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.49694127) q[0];
sx q[0];
rz(-0.48321378) q[0];
sx q[0];
rz(-2.605873) q[0];
rz(0.016853111) q[1];
sx q[1];
rz(-0.87937513) q[1];
sx q[1];
rz(-0.0171612) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099875256) q[0];
sx q[0];
rz(-1.8922531) q[0];
sx q[0];
rz(-1.0519371) q[0];
x q[1];
rz(0.46673985) q[2];
sx q[2];
rz(-1.5108951) q[2];
sx q[2];
rz(-1.2546828) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8507599) q[1];
sx q[1];
rz(-2.8597745) q[1];
sx q[1];
rz(0.59808224) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5480367) q[3];
sx q[3];
rz(-0.16151229) q[3];
sx q[3];
rz(-2.7167729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0422334) q[2];
sx q[2];
rz(-1.6681654) q[2];
sx q[2];
rz(-2.8933706) q[2];
rz(-2.21375) q[3];
sx q[3];
rz(-2.35858) q[3];
sx q[3];
rz(0.46970126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8645653) q[0];
sx q[0];
rz(-1.9643354) q[0];
sx q[0];
rz(2.9789341) q[0];
rz(-2.3795369) q[1];
sx q[1];
rz(-2.9167852) q[1];
sx q[1];
rz(2.3077097) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5637109) q[0];
sx q[0];
rz(-2.1378008) q[0];
sx q[0];
rz(-1.7162697) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.44631) q[2];
sx q[2];
rz(-1.8911749) q[2];
sx q[2];
rz(1.2492391) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0625588) q[1];
sx q[1];
rz(-1.9818881) q[1];
sx q[1];
rz(-2.316019) q[1];
rz(3.0224946) q[3];
sx q[3];
rz(-0.4435215) q[3];
sx q[3];
rz(-1.8325072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.062408058) q[2];
sx q[2];
rz(-2.0024313) q[2];
sx q[2];
rz(2.8774234) q[2];
rz(-1.0583813) q[3];
sx q[3];
rz(-0.94401413) q[3];
sx q[3];
rz(-0.03037608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77085483) q[0];
sx q[0];
rz(-0.83468947) q[0];
sx q[0];
rz(-1.7133065) q[0];
rz(-2.3772073) q[1];
sx q[1];
rz(-1.9995707) q[1];
sx q[1];
rz(-1.8199325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8115615) q[0];
sx q[0];
rz(-1.223765) q[0];
sx q[0];
rz(-2.7156732) q[0];
x q[1];
rz(-3.0899449) q[2];
sx q[2];
rz(-1.902129) q[2];
sx q[2];
rz(0.29690642) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0881867) q[1];
sx q[1];
rz(-1.2866396) q[1];
sx q[1];
rz(-2.8233145) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43784053) q[3];
sx q[3];
rz(-0.90477123) q[3];
sx q[3];
rz(0.36208682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6617714) q[2];
sx q[2];
rz(-1.4134553) q[2];
sx q[2];
rz(-2.856355) q[2];
rz(1.2626922) q[3];
sx q[3];
rz(-1.217239) q[3];
sx q[3];
rz(-1.4234022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406463) q[0];
sx q[0];
rz(-0.74746376) q[0];
sx q[0];
rz(2.2593011) q[0];
rz(0.67059416) q[1];
sx q[1];
rz(-1.0958442) q[1];
sx q[1];
rz(2.1569599) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.479001) q[0];
sx q[0];
rz(-0.07018319) q[0];
sx q[0];
rz(-2.6778738) q[0];
rz(3.0987344) q[2];
sx q[2];
rz(-2.1120666) q[2];
sx q[2];
rz(-1.0554016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24422503) q[1];
sx q[1];
rz(-2.3363718) q[1];
sx q[1];
rz(0.77181384) q[1];
rz(0.71159848) q[3];
sx q[3];
rz(-0.67005537) q[3];
sx q[3];
rz(1.3942476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.275445) q[2];
sx q[2];
rz(-2.6068942) q[2];
sx q[2];
rz(-2.5246485) q[2];
rz(-2.3894737) q[3];
sx q[3];
rz(-1.3277466) q[3];
sx q[3];
rz(0.45929685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1521456) q[0];
sx q[0];
rz(-0.716827) q[0];
sx q[0];
rz(-0.77051198) q[0];
rz(2.41467) q[1];
sx q[1];
rz(-0.22218552) q[1];
sx q[1];
rz(2.4368584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.373823) q[0];
sx q[0];
rz(-2.2859068) q[0];
sx q[0];
rz(-2.3138347) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97874586) q[2];
sx q[2];
rz(-2.5368368) q[2];
sx q[2];
rz(2.3900685) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.875346) q[1];
sx q[1];
rz(-0.35142144) q[1];
sx q[1];
rz(0.09697341) q[1];
x q[2];
rz(1.1868531) q[3];
sx q[3];
rz(-0.60572165) q[3];
sx q[3];
rz(1.5662409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3720234) q[2];
sx q[2];
rz(-1.4906887) q[2];
sx q[2];
rz(-0.43323576) q[2];
rz(0.20089928) q[3];
sx q[3];
rz(-1.9379987) q[3];
sx q[3];
rz(0.5350298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83605003) q[0];
sx q[0];
rz(-1.1490281) q[0];
sx q[0];
rz(0.38247633) q[0];
rz(1.0410694) q[1];
sx q[1];
rz(-1.46547) q[1];
sx q[1];
rz(-0.44630757) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.324821) q[0];
sx q[0];
rz(-1.8292301) q[0];
sx q[0];
rz(1.517061) q[0];
rz(-pi) q[1];
rz(-0.52526229) q[2];
sx q[2];
rz(-0.43516152) q[2];
sx q[2];
rz(2.9829142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1416229) q[1];
sx q[1];
rz(-1.217662) q[1];
sx q[1];
rz(2.1976539) q[1];
rz(-pi) q[2];
rz(-2.2931523) q[3];
sx q[3];
rz(-2.1506243) q[3];
sx q[3];
rz(-2.8715796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3216766) q[2];
sx q[2];
rz(-2.0240462) q[2];
sx q[2];
rz(-0.72810158) q[2];
rz(2.8150832) q[3];
sx q[3];
rz(-1.6312586) q[3];
sx q[3];
rz(0.85730332) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2364872) q[0];
sx q[0];
rz(-0.40920722) q[0];
sx q[0];
rz(1.2976728) q[0];
rz(-1.0009276) q[1];
sx q[1];
rz(-2.3960254) q[1];
sx q[1];
rz(-0.40850684) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1649969) q[0];
sx q[0];
rz(-1.5566096) q[0];
sx q[0];
rz(0.00084695296) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0922008) q[2];
sx q[2];
rz(-1.9754306) q[2];
sx q[2];
rz(1.7737845) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6031987) q[1];
sx q[1];
rz(-1.4567175) q[1];
sx q[1];
rz(1.7054547) q[1];
rz(-pi) q[2];
rz(0.85902416) q[3];
sx q[3];
rz(-2.6005496) q[3];
sx q[3];
rz(0.67160254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0268176) q[2];
sx q[2];
rz(-1.6208181) q[2];
sx q[2];
rz(-0.82359037) q[2];
rz(-0.94281998) q[3];
sx q[3];
rz(-1.7190245) q[3];
sx q[3];
rz(1.5028809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5470062) q[0];
sx q[0];
rz(-2.6972045) q[0];
sx q[0];
rz(-1.8894926) q[0];
rz(-1.3638672) q[1];
sx q[1];
rz(-1.6412647) q[1];
sx q[1];
rz(1.3346671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9329488) q[0];
sx q[0];
rz(-1.874039) q[0];
sx q[0];
rz(-1.5327081) q[0];
x q[1];
rz(-2.4995915) q[2];
sx q[2];
rz(-1.4307012) q[2];
sx q[2];
rz(0.40792712) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.081055911) q[1];
sx q[1];
rz(-1.4976275) q[1];
sx q[1];
rz(-0.45159486) q[1];
rz(-pi) q[2];
rz(1.1030508) q[3];
sx q[3];
rz(-1.333101) q[3];
sx q[3];
rz(-0.89511824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.021726457) q[2];
sx q[2];
rz(-2.848337) q[2];
sx q[2];
rz(0.29652706) q[2];
rz(-1.4793388) q[3];
sx q[3];
rz(-1.5115503) q[3];
sx q[3];
rz(-0.1571981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.6327561) q[0];
sx q[0];
rz(-2.4907676) q[0];
sx q[0];
rz(2.346709) q[0];
rz(-2.5859313) q[1];
sx q[1];
rz(-1.2763005) q[1];
sx q[1];
rz(-1.6937675) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7894365) q[0];
sx q[0];
rz(-1.1934501) q[0];
sx q[0];
rz(-2.3579602) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5083639) q[2];
sx q[2];
rz(-1.0737906) q[2];
sx q[2];
rz(-1.7071498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9196133) q[1];
sx q[1];
rz(-1.3718425) q[1];
sx q[1];
rz(-2.0556021) q[1];
rz(-pi) q[2];
rz(2.6083412) q[3];
sx q[3];
rz(-2.0999895) q[3];
sx q[3];
rz(-1.8797678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8453703) q[2];
sx q[2];
rz(-1.6310548) q[2];
sx q[2];
rz(-1.7108062) q[2];
rz(0.75302643) q[3];
sx q[3];
rz(-0.81792653) q[3];
sx q[3];
rz(-0.16451612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13567781) q[0];
sx q[0];
rz(-2.4527241) q[0];
sx q[0];
rz(1.2149568) q[0];
rz(1.6750492) q[1];
sx q[1];
rz(-2.2129682) q[1];
sx q[1];
rz(0.60563544) q[1];
rz(2.1521592) q[2];
sx q[2];
rz(-2.0579915) q[2];
sx q[2];
rz(1.9111765) q[2];
rz(-3.0343227) q[3];
sx q[3];
rz(-1.7402049) q[3];
sx q[3];
rz(-2.3457478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
