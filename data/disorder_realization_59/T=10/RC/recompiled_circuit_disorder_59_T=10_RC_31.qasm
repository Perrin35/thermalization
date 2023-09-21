OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3988848) q[0];
sx q[0];
rz(-2.3595915) q[0];
sx q[0];
rz(-1.8703823) q[0];
rz(0.27702364) q[1];
sx q[1];
rz(-0.47200534) q[1];
sx q[1];
rz(0.0013874887) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3256677) q[0];
sx q[0];
rz(-2.7357833) q[0];
sx q[0];
rz(1.146846) q[0];
x q[1];
rz(1.6127365) q[2];
sx q[2];
rz(-1.1239927) q[2];
sx q[2];
rz(-0.97181335) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.672294) q[1];
sx q[1];
rz(-1.0743595) q[1];
sx q[1];
rz(1.4908355) q[1];
rz(1.9418342) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(2.1109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(0.74938613) q[2];
rz(2.1253712) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(-0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7063023) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(0.97066561) q[0];
rz(-2.1043815) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(0.81545365) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0464697) q[0];
sx q[0];
rz(-2.2182584) q[0];
sx q[0];
rz(1.4955273) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8380727) q[2];
sx q[2];
rz(-2.3887861) q[2];
sx q[2];
rz(-0.34740651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8784139) q[1];
sx q[1];
rz(-2.8385332) q[1];
sx q[1];
rz(-0.11942272) q[1];
x q[2];
rz(-1.3482413) q[3];
sx q[3];
rz(-0.9369623) q[3];
sx q[3];
rz(-2.8702877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4619535) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3011424) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(-0.87483037) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(2.1420746) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8529352) q[0];
sx q[0];
rz(-1.7203727) q[0];
sx q[0];
rz(1.268671) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61727662) q[2];
sx q[2];
rz(-2.1356138) q[2];
sx q[2];
rz(-1.9895983) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.42400186) q[1];
sx q[1];
rz(-2.058299) q[1];
sx q[1];
rz(1.2815777) q[1];
rz(-pi) q[2];
rz(1.0619034) q[3];
sx q[3];
rz(-1.6403927) q[3];
sx q[3];
rz(-3.0825465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6040566) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.375305) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(0.91039175) q[0];
rz(-2.6903649) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(-2.8667563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1823605) q[0];
sx q[0];
rz(-1.4920456) q[0];
sx q[0];
rz(1.6596646) q[0];
x q[1];
rz(-1.574013) q[2];
sx q[2];
rz(-0.46041691) q[2];
sx q[2];
rz(0.38052961) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.53510016) q[1];
sx q[1];
rz(-2.1537158) q[1];
sx q[1];
rz(0.54775723) q[1];
rz(-pi) q[2];
rz(1.4196017) q[3];
sx q[3];
rz(-0.70221838) q[3];
sx q[3];
rz(-0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.794902) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(-1.6332731) q[2];
rz(-1.9968962) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(0.16170734) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4836924) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(2.143798) q[0];
rz(0.18355852) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(1.516974) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51380101) q[0];
sx q[0];
rz(-0.9538981) q[0];
sx q[0];
rz(0.4517201) q[0];
x q[1];
rz(2.6078569) q[2];
sx q[2];
rz(-1.9118475) q[2];
sx q[2];
rz(-1.595572) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4167037) q[1];
sx q[1];
rz(-0.54425889) q[1];
sx q[1];
rz(-0.7087899) q[1];
x q[2];
rz(0.59668221) q[3];
sx q[3];
rz(-2.0697069) q[3];
sx q[3];
rz(-0.98017207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(2.8175763) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76535392) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(1.8776241) q[0];
rz(-0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(-0.34067672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6402123) q[0];
sx q[0];
rz(-3.0684154) q[0];
sx q[0];
rz(-1.120938) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2773877) q[2];
sx q[2];
rz(-2.3050606) q[2];
sx q[2];
rz(0.66213995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9601945) q[1];
sx q[1];
rz(-1.3409233) q[1];
sx q[1];
rz(0.48744907) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0601222) q[3];
sx q[3];
rz(-1.9483856) q[3];
sx q[3];
rz(-1.1946354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1910151) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(0.54780444) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1691549) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(0.46494928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60663215) q[0];
sx q[0];
rz(-0.50919845) q[0];
sx q[0];
rz(1.1688933) q[0];
rz(-pi) q[1];
rz(-0.99636997) q[2];
sx q[2];
rz(-1.1466221) q[2];
sx q[2];
rz(-3.0976354) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3073687) q[1];
sx q[1];
rz(-0.24337473) q[1];
sx q[1];
rz(0.4537531) q[1];
rz(-pi) q[2];
rz(3.1157324) q[3];
sx q[3];
rz(-1.5246632) q[3];
sx q[3];
rz(-2.9010454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1376301) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(-0.31759343) q[2];
rz(-0.57146227) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(2.8572594) q[0];
rz(2.590086) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(3.0632339) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4975472) q[0];
sx q[0];
rz(-0.80701485) q[0];
sx q[0];
rz(0.09597309) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8163082) q[2];
sx q[2];
rz(-1.3405372) q[2];
sx q[2];
rz(2.2895209) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9662712) q[1];
sx q[1];
rz(-0.92370719) q[1];
sx q[1];
rz(2.9233169) q[1];
rz(-1.4140698) q[3];
sx q[3];
rz(-1.5369475) q[3];
sx q[3];
rz(-2.1050997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7408961) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437929) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(-3.0902241) q[0];
rz(2.2180166) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-0.87402469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066814518) q[0];
sx q[0];
rz(-1.5621119) q[0];
sx q[0];
rz(2.3735223) q[0];
x q[1];
rz(-0.36058493) q[2];
sx q[2];
rz(-1.9546095) q[2];
sx q[2];
rz(-1.6894481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2487138) q[1];
sx q[1];
rz(-1.4993748) q[1];
sx q[1];
rz(1.0700657) q[1];
rz(-2.6690528) q[3];
sx q[3];
rz(-2.4774385) q[3];
sx q[3];
rz(-3.1411375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41436568) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(-1.9571346) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(-1.7782036) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0062362) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(2.4172879) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(1.1788517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2154685) q[0];
sx q[0];
rz(-1.7755531) q[0];
sx q[0];
rz(1.8490851) q[0];
rz(2.5095021) q[2];
sx q[2];
rz(-0.1569911) q[2];
sx q[2];
rz(1.1319515) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57943343) q[1];
sx q[1];
rz(-2.4907618) q[1];
sx q[1];
rz(1.8821554) q[1];
rz(1.3445271) q[3];
sx q[3];
rz(-2.0824277) q[3];
sx q[3];
rz(2.1567878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0835691) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(-0.89938346) q[2];
rz(-2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62109229) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(2.3836366) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(1.4964676) q[2];
sx q[2];
rz(-0.44114124) q[2];
sx q[2];
rz(0.91976358) q[2];
rz(-2.6562128) q[3];
sx q[3];
rz(-0.30967064) q[3];
sx q[3];
rz(-2.526024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
