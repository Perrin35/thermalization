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
rz(-1.7669825) q[0];
sx q[0];
rz(-0.97281015) q[0];
sx q[0];
rz(-1.2306124) q[0];
rz(-1.542701) q[1];
sx q[1];
rz(-0.35294947) q[1];
sx q[1];
rz(-2.7065839) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6957851) q[0];
sx q[0];
rz(-1.2381993) q[0];
sx q[0];
rz(0.10528721) q[0];
rz(1.4392008) q[2];
sx q[2];
rz(-1.3450661) q[2];
sx q[2];
rz(0.44039681) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3865996) q[1];
sx q[1];
rz(-2.3191929) q[1];
sx q[1];
rz(1.9564641) q[1];
rz(-1.1633662) q[3];
sx q[3];
rz(-3.1071783) q[3];
sx q[3];
rz(-0.84190166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.9445442) q[2];
sx q[2];
rz(-0.41434449) q[2];
sx q[2];
rz(-0.91828263) q[2];
rz(-0.36327547) q[3];
sx q[3];
rz(-1.8922292) q[3];
sx q[3];
rz(-1.3242807) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0488224) q[0];
sx q[0];
rz(-2.4019882) q[0];
sx q[0];
rz(-1.1356461) q[0];
rz(0.28226918) q[1];
sx q[1];
rz(-1.5156563) q[1];
sx q[1];
rz(0.20271066) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7736063) q[0];
sx q[0];
rz(-1.752089) q[0];
sx q[0];
rz(1.7530935) q[0];
x q[1];
rz(0.01303444) q[2];
sx q[2];
rz(-1.0891277) q[2];
sx q[2];
rz(2.3866057) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2308826) q[1];
sx q[1];
rz(-0.47180125) q[1];
sx q[1];
rz(-0.27179007) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1003157) q[3];
sx q[3];
rz(-0.011547877) q[3];
sx q[3];
rz(-1.601416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2477766) q[2];
sx q[2];
rz(-2.717369) q[2];
sx q[2];
rz(0.068537863) q[2];
rz(2.2872772) q[3];
sx q[3];
rz(-1.8770542) q[3];
sx q[3];
rz(-3.1356649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.051556) q[0];
sx q[0];
rz(-0.55110252) q[0];
sx q[0];
rz(1.0968444) q[0];
rz(-0.54496533) q[1];
sx q[1];
rz(-0.68383354) q[1];
sx q[1];
rz(0.1942689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504963) q[0];
sx q[0];
rz(-1.9197113) q[0];
sx q[0];
rz(2.0112841) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33325573) q[2];
sx q[2];
rz(-2.1001847) q[2];
sx q[2];
rz(-0.80798244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.0021469963) q[1];
sx q[1];
rz(-1.0771828) q[1];
sx q[1];
rz(-1.1612372) q[1];
rz(-pi) q[2];
rz(-1.9220334) q[3];
sx q[3];
rz(-2.7673169) q[3];
sx q[3];
rz(2.5881899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7668309) q[2];
sx q[2];
rz(-1.7685879) q[2];
sx q[2];
rz(2.4135446) q[2];
rz(2.9211365) q[3];
sx q[3];
rz(-0.15153344) q[3];
sx q[3];
rz(-0.63173405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8927638) q[0];
sx q[0];
rz(-1.4735824) q[0];
sx q[0];
rz(2.072075) q[0];
rz(2.2546774) q[1];
sx q[1];
rz(-2.6353757) q[1];
sx q[1];
rz(1.8198397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2489172) q[0];
sx q[0];
rz(-1.6098611) q[0];
sx q[0];
rz(2.7075753) q[0];
rz(-pi) q[1];
x q[1];
rz(0.063912674) q[2];
sx q[2];
rz(-0.23944975) q[2];
sx q[2];
rz(2.4452345) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2414714) q[1];
sx q[1];
rz(-1.6423499) q[1];
sx q[1];
rz(-1.1318737) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6245313) q[3];
sx q[3];
rz(-1.0365465) q[3];
sx q[3];
rz(-1.8815966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47361031) q[2];
sx q[2];
rz(-1.0429635) q[2];
sx q[2];
rz(-1.2874416) q[2];
rz(1.0188262) q[3];
sx q[3];
rz(-0.5580709) q[3];
sx q[3];
rz(-0.67600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92893112) q[0];
sx q[0];
rz(-0.23290817) q[0];
sx q[0];
rz(2.2610597) q[0];
rz(0.13678837) q[1];
sx q[1];
rz(-2.9158178) q[1];
sx q[1];
rz(-2.5313964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6132122) q[0];
sx q[0];
rz(-0.92862178) q[0];
sx q[0];
rz(-0.015676966) q[0];
rz(-1.7238823) q[2];
sx q[2];
rz(-1.9979852) q[2];
sx q[2];
rz(-1.2059905) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6366263) q[1];
sx q[1];
rz(-1.7646878) q[1];
sx q[1];
rz(2.8411179) q[1];
rz(-pi) q[2];
rz(1.4349157) q[3];
sx q[3];
rz(-2.0066094) q[3];
sx q[3];
rz(-2.6282981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5337164) q[2];
sx q[2];
rz(-2.4166962) q[2];
sx q[2];
rz(-1.2302715) q[2];
rz(-0.88165927) q[3];
sx q[3];
rz(-1.449838) q[3];
sx q[3];
rz(-1.8264044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-0.45516685) q[0];
sx q[0];
rz(-1.5979586) q[0];
sx q[0];
rz(1.0536739) q[0];
rz(-1.7651438) q[1];
sx q[1];
rz(-1.1076628) q[1];
sx q[1];
rz(2.460316) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69619715) q[0];
sx q[0];
rz(-0.70022178) q[0];
sx q[0];
rz(0.025431319) q[0];
rz(-1.2238962) q[2];
sx q[2];
rz(-0.85641801) q[2];
sx q[2];
rz(-2.6831085) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5743915) q[1];
sx q[1];
rz(-1.302845) q[1];
sx q[1];
rz(-2.0824357) q[1];
rz(-pi) q[2];
rz(-0.87879718) q[3];
sx q[3];
rz(-1.5996154) q[3];
sx q[3];
rz(-1.6520713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.072558746) q[2];
sx q[2];
rz(-2.7636038) q[2];
sx q[2];
rz(-2.6663713) q[2];
rz(1.9184939) q[3];
sx q[3];
rz(-2.1971072) q[3];
sx q[3];
rz(-0.16330115) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.908919) q[0];
sx q[0];
rz(-0.95218807) q[0];
sx q[0];
rz(-2.9222144) q[0];
rz(-2.1351922) q[1];
sx q[1];
rz(-1.4234797) q[1];
sx q[1];
rz(-1.9361608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74586263) q[0];
sx q[0];
rz(-0.59915483) q[0];
sx q[0];
rz(1.5320734) q[0];
x q[1];
rz(-2.7458887) q[2];
sx q[2];
rz(-1.8181516) q[2];
sx q[2];
rz(0.57952124) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6483375) q[1];
sx q[1];
rz(-1.3715136) q[1];
sx q[1];
rz(-1.553345) q[1];
rz(-1.1262283) q[3];
sx q[3];
rz(-1.4476814) q[3];
sx q[3];
rz(3.0266704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9203303) q[2];
sx q[2];
rz(-0.39322501) q[2];
sx q[2];
rz(-0.66366759) q[2];
rz(-0.71781939) q[3];
sx q[3];
rz(-0.65687537) q[3];
sx q[3];
rz(-2.9653911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32888907) q[0];
sx q[0];
rz(-0.95903522) q[0];
sx q[0];
rz(-1.9004199) q[0];
rz(2.7136956) q[1];
sx q[1];
rz(-1.5792081) q[1];
sx q[1];
rz(-1.7180299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0421906) q[0];
sx q[0];
rz(-1.4960327) q[0];
sx q[0];
rz(0.21850234) q[0];
rz(-pi) q[1];
rz(0.85278089) q[2];
sx q[2];
rz(-1.4508392) q[2];
sx q[2];
rz(-0.44242417) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9291722) q[1];
sx q[1];
rz(-0.66832405) q[1];
sx q[1];
rz(1.5364439) q[1];
rz(1.5345433) q[3];
sx q[3];
rz(-1.2360308) q[3];
sx q[3];
rz(-2.436155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9972035) q[2];
sx q[2];
rz(-0.68789566) q[2];
sx q[2];
rz(-2.8456861) q[2];
rz(2.0603254) q[3];
sx q[3];
rz(-1.6176977) q[3];
sx q[3];
rz(-2.9937939) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57701552) q[0];
sx q[0];
rz(-0.58203375) q[0];
sx q[0];
rz(-0.70593315) q[0];
rz(-0.18711807) q[1];
sx q[1];
rz(-0.83693224) q[1];
sx q[1];
rz(2.6954938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68413823) q[0];
sx q[0];
rz(-1.7230526) q[0];
sx q[0];
rz(1.5726552) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99775935) q[2];
sx q[2];
rz(-2.5964542) q[2];
sx q[2];
rz(1.0447431) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1343119) q[1];
sx q[1];
rz(-2.4680302) q[1];
sx q[1];
rz(2.1271655) q[1];
rz(-pi) q[2];
rz(0.57658617) q[3];
sx q[3];
rz(-1.714725) q[3];
sx q[3];
rz(-1.3670539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17911653) q[2];
sx q[2];
rz(-2.1160782) q[2];
sx q[2];
rz(-1.7402488) q[2];
rz(0.60351795) q[3];
sx q[3];
rz(-2.9602435) q[3];
sx q[3];
rz(0.13874273) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1097581) q[0];
sx q[0];
rz(-1.7367481) q[0];
sx q[0];
rz(0.11012125) q[0];
rz(-1.3212063) q[1];
sx q[1];
rz(-0.92480129) q[1];
sx q[1];
rz(0.22470156) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53407828) q[0];
sx q[0];
rz(-2.461643) q[0];
sx q[0];
rz(-2.2282766) q[0];
rz(-pi) q[1];
rz(2.9797033) q[2];
sx q[2];
rz(-1.8542552) q[2];
sx q[2];
rz(3.0155011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1096829) q[1];
sx q[1];
rz(-0.8485629) q[1];
sx q[1];
rz(-1.3501881) q[1];
rz(-pi) q[2];
rz(-1.3984001) q[3];
sx q[3];
rz(-1.5972923) q[3];
sx q[3];
rz(-2.4717028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3118185) q[2];
sx q[2];
rz(-0.31121397) q[2];
sx q[2];
rz(0.52660006) q[2];
rz(-0.38463587) q[3];
sx q[3];
rz(-1.2472109) q[3];
sx q[3];
rz(2.9068936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93152355) q[0];
sx q[0];
rz(-0.43545224) q[0];
sx q[0];
rz(-0.42238105) q[0];
rz(2.3474563) q[1];
sx q[1];
rz(-1.7105449) q[1];
sx q[1];
rz(1.7078043) q[1];
rz(0.64308249) q[2];
sx q[2];
rz(-0.98862917) q[2];
sx q[2];
rz(0.25884982) q[2];
rz(-1.9392813) q[3];
sx q[3];
rz(-1.9574584) q[3];
sx q[3];
rz(-3.0493469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
