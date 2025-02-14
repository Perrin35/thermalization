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
rz(1.9109803) q[0];
rz(1.5988916) q[1];
sx q[1];
rz(3.4945421) q[1];
sx q[1];
rz(8.9897692) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6957851) q[0];
sx q[0];
rz(-1.2381993) q[0];
sx q[0];
rz(3.0363054) q[0];
rz(-0.22763197) q[2];
sx q[2];
rz(-1.6990347) q[2];
sx q[2];
rz(-1.1007835) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8486464) q[1];
sx q[1];
rz(-0.82445626) q[1];
sx q[1];
rz(2.7566977) q[1];
x q[2];
rz(1.1633662) q[3];
sx q[3];
rz(-3.1071783) q[3];
sx q[3];
rz(0.84190166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9445442) q[2];
sx q[2];
rz(-2.7272482) q[2];
sx q[2];
rz(0.91828263) q[2];
rz(-2.7783172) q[3];
sx q[3];
rz(-1.8922292) q[3];
sx q[3];
rz(-1.8173119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0927703) q[0];
sx q[0];
rz(-0.73960441) q[0];
sx q[0];
rz(2.0059465) q[0];
rz(-0.28226918) q[1];
sx q[1];
rz(-1.6259364) q[1];
sx q[1];
rz(-2.938882) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36798635) q[0];
sx q[0];
rz(-1.3895036) q[0];
sx q[0];
rz(1.3884991) q[0];
rz(-1.0890929) q[2];
sx q[2];
rz(-1.559245) q[2];
sx q[2];
rz(2.3318219) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.41676471) q[1];
sx q[1];
rz(-1.4484805) q[1];
sx q[1];
rz(-0.45682795) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5810896) q[3];
sx q[3];
rz(-1.5760311) q[3];
sx q[3];
rz(-2.7017587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89381605) q[2];
sx q[2];
rz(-2.717369) q[2];
sx q[2];
rz(-0.068537863) q[2];
rz(2.2872772) q[3];
sx q[3];
rz(-1.8770542) q[3];
sx q[3];
rz(-3.1356649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.051556) q[0];
sx q[0];
rz(-0.55110252) q[0];
sx q[0];
rz(-2.0447482) q[0];
rz(0.54496533) q[1];
sx q[1];
rz(-0.68383354) q[1];
sx q[1];
rz(2.9473238) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30667728) q[0];
sx q[0];
rz(-2.5868572) q[0];
sx q[0];
rz(0.86443211) q[0];
rz(1.0164073) q[2];
sx q[2];
rz(-1.28456) q[2];
sx q[2];
rz(0.58974671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.771505) q[1];
sx q[1];
rz(-1.2125101) q[1];
sx q[1];
rz(-0.53046988) q[1];
x q[2];
rz(-1.9220334) q[3];
sx q[3];
rz(-0.3742758) q[3];
sx q[3];
rz(-2.5881899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37476173) q[2];
sx q[2];
rz(-1.7685879) q[2];
sx q[2];
rz(2.4135446) q[2];
rz(0.22045615) q[3];
sx q[3];
rz(-2.9900592) q[3];
sx q[3];
rz(-0.63173405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24882889) q[0];
sx q[0];
rz(-1.6680102) q[0];
sx q[0];
rz(-1.0695176) q[0];
rz(0.88691521) q[1];
sx q[1];
rz(-0.506217) q[1];
sx q[1];
rz(1.8198397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4815714) q[0];
sx q[0];
rz(-1.1371326) q[0];
sx q[0];
rz(1.5277442) q[0];
rz(-pi) q[1];
rz(0.063912674) q[2];
sx q[2];
rz(-2.9021429) q[2];
sx q[2];
rz(0.69635812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9001213) q[1];
sx q[1];
rz(-1.6423499) q[1];
sx q[1];
rz(-1.1318737) q[1];
rz(0.87476829) q[3];
sx q[3];
rz(-0.7255377) q[3];
sx q[3];
rz(1.0410795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47361031) q[2];
sx q[2];
rz(-1.0429635) q[2];
sx q[2];
rz(-1.854151) q[2];
rz(2.1227664) q[3];
sx q[3];
rz(-2.5835218) q[3];
sx q[3];
rz(-0.67600018) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92893112) q[0];
sx q[0];
rz(-2.9086845) q[0];
sx q[0];
rz(2.2610597) q[0];
rz(-0.13678837) q[1];
sx q[1];
rz(-2.9158178) q[1];
sx q[1];
rz(-0.61019623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1085668) q[0];
sx q[0];
rz(-1.5833502) q[0];
sx q[0];
rz(-0.92856284) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32322804) q[2];
sx q[2];
rz(-2.6893977) q[2];
sx q[2];
rz(1.5791073) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1353947) q[1];
sx q[1];
rz(-1.2761226) q[1];
sx q[1];
rz(-1.7735405) q[1];
rz(-pi) q[2];
rz(0.28308308) q[3];
sx q[3];
rz(-2.6863881) q[3];
sx q[3];
rz(-2.9415123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.60787624) q[2];
sx q[2];
rz(-0.72489649) q[2];
sx q[2];
rz(1.2302715) q[2];
rz(-0.88165927) q[3];
sx q[3];
rz(-1.6917546) q[3];
sx q[3];
rz(1.8264044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6864258) q[0];
sx q[0];
rz(-1.5979586) q[0];
sx q[0];
rz(-2.0879188) q[0];
rz(1.3764489) q[1];
sx q[1];
rz(-2.0339298) q[1];
sx q[1];
rz(-2.460316) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69619715) q[0];
sx q[0];
rz(-0.70022178) q[0];
sx q[0];
rz(3.1161613) q[0];
rz(-1.9176964) q[2];
sx q[2];
rz(-2.2851746) q[2];
sx q[2];
rz(-2.6831085) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5743915) q[1];
sx q[1];
rz(-1.302845) q[1];
sx q[1];
rz(2.0824357) q[1];
rz(-pi) q[2];
rz(2.2627955) q[3];
sx q[3];
rz(-1.5996154) q[3];
sx q[3];
rz(-1.6520713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0690339) q[2];
sx q[2];
rz(-0.37798887) q[2];
sx q[2];
rz(2.6663713) q[2];
rz(1.9184939) q[3];
sx q[3];
rz(-2.1971072) q[3];
sx q[3];
rz(2.9782915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.2326736) q[0];
sx q[0];
rz(-0.95218807) q[0];
sx q[0];
rz(-0.21937823) q[0];
rz(1.0064005) q[1];
sx q[1];
rz(-1.4234797) q[1];
sx q[1];
rz(1.2054319) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4426098) q[0];
sx q[0];
rz(-2.1694393) q[0];
sx q[0];
rz(-3.1151616) q[0];
x q[1];
rz(-0.57997616) q[2];
sx q[2];
rz(-2.6784228) q[2];
sx q[2];
rz(-1.620174) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58118966) q[1];
sx q[1];
rz(-0.20003527) q[1];
sx q[1];
rz(0.086189857) q[1];
x q[2];
rz(1.2906399) q[3];
sx q[3];
rz(-0.46020111) q[3];
sx q[3];
rz(1.4334219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22126234) q[2];
sx q[2];
rz(-2.7483676) q[2];
sx q[2];
rz(2.4779251) q[2];
rz(-2.4237733) q[3];
sx q[3];
rz(-0.65687537) q[3];
sx q[3];
rz(2.9653911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8127036) q[0];
sx q[0];
rz(-2.1825574) q[0];
sx q[0];
rz(-1.9004199) q[0];
rz(-2.7136956) q[1];
sx q[1];
rz(-1.5792081) q[1];
sx q[1];
rz(1.7180299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9374838) q[0];
sx q[0];
rz(-0.23074575) q[0];
sx q[0];
rz(2.8088914) q[0];
x q[1];
rz(1.752002) q[2];
sx q[2];
rz(-0.7262035) q[2];
sx q[2];
rz(2.1493634) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2124205) q[1];
sx q[1];
rz(-2.4732686) q[1];
sx q[1];
rz(1.5364439) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0377722) q[3];
sx q[3];
rz(-2.8049433) q[3];
sx q[3];
rz(2.3262084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9972035) q[2];
sx q[2];
rz(-0.68789566) q[2];
sx q[2];
rz(2.8456861) q[2];
rz(1.0812673) q[3];
sx q[3];
rz(-1.6176977) q[3];
sx q[3];
rz(-0.14779873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.57701552) q[0];
sx q[0];
rz(-0.58203375) q[0];
sx q[0];
rz(2.4356595) q[0];
rz(2.9544746) q[1];
sx q[1];
rz(-2.3046604) q[1];
sx q[1];
rz(0.44609889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68413823) q[0];
sx q[0];
rz(-1.7230526) q[0];
sx q[0];
rz(-1.5726552) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99775935) q[2];
sx q[2];
rz(-0.54513844) q[2];
sx q[2];
rz(1.0447431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0072807) q[1];
sx q[1];
rz(-2.4680302) q[1];
sx q[1];
rz(1.0144272) q[1];
rz(-1.7419849) q[3];
sx q[3];
rz(-1.0009196) q[3];
sx q[3];
rz(-0.2967473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9624761) q[2];
sx q[2];
rz(-2.1160782) q[2];
sx q[2];
rz(1.7402488) q[2];
rz(2.5380747) q[3];
sx q[3];
rz(-2.9602435) q[3];
sx q[3];
rz(-0.13874273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-3.1097581) q[0];
sx q[0];
rz(-1.4048445) q[0];
sx q[0];
rz(-3.0314714) q[0];
rz(-1.3212063) q[1];
sx q[1];
rz(-0.92480129) q[1];
sx q[1];
rz(0.22470156) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6455655) q[0];
sx q[0];
rz(-1.9651881) q[0];
sx q[0];
rz(1.0014707) q[0];
rz(-1.857809) q[2];
sx q[2];
rz(-1.415421) q[2];
sx q[2];
rz(-1.7425328) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4367392) q[1];
sx q[1];
rz(-0.74932832) q[1];
sx q[1];
rz(2.8981461) q[1];
rz(-pi) q[2];
rz(1.7240771) q[3];
sx q[3];
rz(-0.17440052) q[3];
sx q[3];
rz(-2.0896951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3118185) q[2];
sx q[2];
rz(-0.31121397) q[2];
sx q[2];
rz(2.6149926) q[2];
rz(-2.7569568) q[3];
sx q[3];
rz(-1.8943818) q[3];
sx q[3];
rz(2.9068936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.2100691) q[0];
sx q[0];
rz(-2.7061404) q[0];
sx q[0];
rz(2.7192116) q[0];
rz(2.3474563) q[1];
sx q[1];
rz(-1.7105449) q[1];
sx q[1];
rz(1.7078043) q[1];
rz(-2.3096397) q[2];
sx q[2];
rz(-0.83870287) q[2];
sx q[2];
rz(2.4626682) q[2];
rz(1.9392813) q[3];
sx q[3];
rz(-1.1841342) q[3];
sx q[3];
rz(0.092245734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
