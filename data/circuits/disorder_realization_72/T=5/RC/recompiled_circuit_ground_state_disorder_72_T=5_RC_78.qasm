OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.44719625) q[0];
sx q[0];
rz(-1.5153272) q[0];
sx q[0];
rz(0.84994999) q[0];
rz(1.3810459) q[1];
sx q[1];
rz(-0.84264207) q[1];
sx q[1];
rz(-3.091264) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.91431) q[0];
sx q[0];
rz(-1.3903872) q[0];
sx q[0];
rz(-1.6561899) q[0];
x q[1];
rz(-1.8036795) q[2];
sx q[2];
rz(-0.63265975) q[2];
sx q[2];
rz(-0.92756995) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.58568629) q[1];
sx q[1];
rz(-1.7650604) q[1];
sx q[1];
rz(1.6418841) q[1];
x q[2];
rz(-1.6292455) q[3];
sx q[3];
rz(-2.2104467) q[3];
sx q[3];
rz(-0.15098886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7488659) q[2];
sx q[2];
rz(-1.2458845) q[2];
sx q[2];
rz(2.5457814) q[2];
rz(-0.81123224) q[3];
sx q[3];
rz(-1.8758643) q[3];
sx q[3];
rz(2.5442512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8314464) q[0];
sx q[0];
rz(-0.50066384) q[0];
sx q[0];
rz(1.3432107) q[0];
rz(1.385618) q[1];
sx q[1];
rz(-2.0090397) q[1];
sx q[1];
rz(1.064942) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5249817) q[0];
sx q[0];
rz(-3.0654902) q[0];
sx q[0];
rz(0.2790926) q[0];
rz(-pi) q[1];
rz(-1.0491514) q[2];
sx q[2];
rz(-1.8966833) q[2];
sx q[2];
rz(2.2679595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2328223) q[1];
sx q[1];
rz(-1.576464) q[1];
sx q[1];
rz(-3.0012111) q[1];
rz(-1.4707758) q[3];
sx q[3];
rz(-1.6006114) q[3];
sx q[3];
rz(-0.884197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6407577) q[2];
sx q[2];
rz(-1.0229599) q[2];
sx q[2];
rz(0.41066059) q[2];
rz(-1.143035) q[3];
sx q[3];
rz(-0.7468907) q[3];
sx q[3];
rz(0.66543287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9329342) q[0];
sx q[0];
rz(-1.3422796) q[0];
sx q[0];
rz(2.4564273) q[0];
rz(2.4655828) q[1];
sx q[1];
rz(-1.490255) q[1];
sx q[1];
rz(2.4032059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9907407) q[0];
sx q[0];
rz(-2.5408077) q[0];
sx q[0];
rz(-0.59101253) q[0];
x q[1];
rz(-3.1072363) q[2];
sx q[2];
rz(-0.7329251) q[2];
sx q[2];
rz(0.71897163) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1967839) q[1];
sx q[1];
rz(-0.15722188) q[1];
sx q[1];
rz(-1.9095837) q[1];
rz(2.1687968) q[3];
sx q[3];
rz(-2.1913652) q[3];
sx q[3];
rz(0.52956328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43397778) q[2];
sx q[2];
rz(-1.4716163) q[2];
sx q[2];
rz(0.99226704) q[2];
rz(0.17539242) q[3];
sx q[3];
rz(-2.6447191) q[3];
sx q[3];
rz(2.8309256) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6076412) q[0];
sx q[0];
rz(-1.1247922) q[0];
sx q[0];
rz(-0.24096179) q[0];
rz(0.85722771) q[1];
sx q[1];
rz(-0.78397426) q[1];
sx q[1];
rz(1.2711752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60343307) q[0];
sx q[0];
rz(-2.6707244) q[0];
sx q[0];
rz(1.9290857) q[0];
rz(-pi) q[1];
rz(-0.52267646) q[2];
sx q[2];
rz(-1.1416404) q[2];
sx q[2];
rz(1.9945838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.037447151) q[1];
sx q[1];
rz(-2.0097783) q[1];
sx q[1];
rz(-0.46260117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.70636) q[3];
sx q[3];
rz(-0.80343825) q[3];
sx q[3];
rz(0.17882192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50095373) q[2];
sx q[2];
rz(-0.95879889) q[2];
sx q[2];
rz(2.1232088) q[2];
rz(1.9765249) q[3];
sx q[3];
rz(-2.1924721) q[3];
sx q[3];
rz(-2.693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.0388221) q[0];
sx q[0];
rz(-2.3508681) q[0];
sx q[0];
rz(0.08654174) q[0];
rz(1.6697007) q[1];
sx q[1];
rz(-2.4457928) q[1];
sx q[1];
rz(-2.5873628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51064483) q[0];
sx q[0];
rz(-1.3501337) q[0];
sx q[0];
rz(2.4359926) q[0];
rz(-pi) q[1];
rz(-2.9678095) q[2];
sx q[2];
rz(-0.67774583) q[2];
sx q[2];
rz(0.89231561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.94658417) q[1];
sx q[1];
rz(-1.1145419) q[1];
sx q[1];
rz(-0.66000451) q[1];
rz(-1.1301) q[3];
sx q[3];
rz(-1.8365068) q[3];
sx q[3];
rz(1.2540224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.294813) q[2];
sx q[2];
rz(-2.7816935) q[2];
sx q[2];
rz(2.2100718) q[2];
rz(-2.8281663) q[3];
sx q[3];
rz(-0.64160186) q[3];
sx q[3];
rz(0.48345598) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.038092) q[0];
sx q[0];
rz(-2.3374538) q[0];
sx q[0];
rz(-1.3793797) q[0];
rz(2.3992959) q[1];
sx q[1];
rz(-1.392044) q[1];
sx q[1];
rz(-2.0753863) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8277513) q[0];
sx q[0];
rz(-3.0004639) q[0];
sx q[0];
rz(-2.9147097) q[0];
rz(-pi) q[1];
rz(-3.0020038) q[2];
sx q[2];
rz(-2.0878279) q[2];
sx q[2];
rz(-1.4429744) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71519796) q[1];
sx q[1];
rz(-2.4860588) q[1];
sx q[1];
rz(-1.1525201) q[1];
x q[2];
rz(2.8703254) q[3];
sx q[3];
rz(-1.4888797) q[3];
sx q[3];
rz(1.0759491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.958287) q[2];
sx q[2];
rz(-0.97465193) q[2];
sx q[2];
rz(-0.24678123) q[2];
rz(1.0684446) q[3];
sx q[3];
rz(-1.1698134) q[3];
sx q[3];
rz(1.0698414) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5088155) q[0];
sx q[0];
rz(-2.1822073) q[0];
sx q[0];
rz(2.7799613) q[0];
rz(-1.5879141) q[1];
sx q[1];
rz(-1.5767153) q[1];
sx q[1];
rz(-1.9379001) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4483245) q[0];
sx q[0];
rz(-2.0965323) q[0];
sx q[0];
rz(2.3828302) q[0];
rz(-pi) q[1];
rz(-2.6404722) q[2];
sx q[2];
rz(-1.7351056) q[2];
sx q[2];
rz(-1.1962547) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3059235) q[1];
sx q[1];
rz(-0.9811372) q[1];
sx q[1];
rz(2.8980189) q[1];
x q[2];
rz(-2.5011115) q[3];
sx q[3];
rz(-0.52947097) q[3];
sx q[3];
rz(-2.01628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2162073) q[2];
sx q[2];
rz(-1.46571) q[2];
sx q[2];
rz(-2.1290131) q[2];
rz(0.25977627) q[3];
sx q[3];
rz(-0.24661073) q[3];
sx q[3];
rz(0.40464211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15553661) q[0];
sx q[0];
rz(-1.8590834) q[0];
sx q[0];
rz(0.78722659) q[0];
rz(-0.8440482) q[1];
sx q[1];
rz(-2.7828352) q[1];
sx q[1];
rz(1.2740096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.781141) q[0];
sx q[0];
rz(-1.5119189) q[0];
sx q[0];
rz(0.33712338) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2574848) q[2];
sx q[2];
rz(-0.85585218) q[2];
sx q[2];
rz(2.7471209) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7374271) q[1];
sx q[1];
rz(-2.3336172) q[1];
sx q[1];
rz(0.91973181) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6350217) q[3];
sx q[3];
rz(-1.1665441) q[3];
sx q[3];
rz(0.4413213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0446223) q[2];
sx q[2];
rz(-1.5771834) q[2];
sx q[2];
rz(0.38193646) q[2];
rz(-1.1218128) q[3];
sx q[3];
rz(-0.33676454) q[3];
sx q[3];
rz(-0.063966123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065351039) q[0];
sx q[0];
rz(-2.2168646) q[0];
sx q[0];
rz(1.4724154) q[0];
rz(-2.8352101) q[1];
sx q[1];
rz(-2.6173321) q[1];
sx q[1];
rz(-0.27149567) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24105212) q[0];
sx q[0];
rz(-1.1212305) q[0];
sx q[0];
rz(-2.7009214) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1708529) q[2];
sx q[2];
rz(-0.78151199) q[2];
sx q[2];
rz(1.140238) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9418464) q[1];
sx q[1];
rz(-1.0916657) q[1];
sx q[1];
rz(0.61750827) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90320194) q[3];
sx q[3];
rz(-2.4861504) q[3];
sx q[3];
rz(-3.1208861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86758119) q[2];
sx q[2];
rz(-2.0143955) q[2];
sx q[2];
rz(2.4809044) q[2];
rz(0.88179669) q[3];
sx q[3];
rz(-0.56395689) q[3];
sx q[3];
rz(-2.2753415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26447403) q[0];
sx q[0];
rz(-1.6082123) q[0];
sx q[0];
rz(-1.8033173) q[0];
rz(1.0461461) q[1];
sx q[1];
rz(-2.1144861) q[1];
sx q[1];
rz(-0.17790067) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7545181) q[0];
sx q[0];
rz(-0.99102586) q[0];
sx q[0];
rz(2.3405318) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17651813) q[2];
sx q[2];
rz(-2.1452139) q[2];
sx q[2];
rz(-0.18744577) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2819351) q[1];
sx q[1];
rz(-1.9023696) q[1];
sx q[1];
rz(2.8012707) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7081065) q[3];
sx q[3];
rz(-1.2009283) q[3];
sx q[3];
rz(-2.9183877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.88343128) q[2];
sx q[2];
rz(-1.0940172) q[2];
sx q[2];
rz(2.2056313) q[2];
rz(0.59518138) q[3];
sx q[3];
rz(-1.7116356) q[3];
sx q[3];
rz(-0.87966758) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7137939) q[0];
sx q[0];
rz(-1.4604764) q[0];
sx q[0];
rz(-1.9124734) q[0];
rz(-1.8145369) q[1];
sx q[1];
rz(-1.6358903) q[1];
sx q[1];
rz(-1.9930175) q[1];
rz(-3.0279866) q[2];
sx q[2];
rz(-1.3945711) q[2];
sx q[2];
rz(-0.61749189) q[2];
rz(-2.280664) q[3];
sx q[3];
rz(-0.39562514) q[3];
sx q[3];
rz(2.8277314) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
