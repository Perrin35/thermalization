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
rz(-1.1620698) q[0];
sx q[0];
rz(-0.23482366) q[0];
sx q[0];
rz(-1.809037) q[0];
rz(-1.8279583) q[1];
sx q[1];
rz(5.0143427) q[1];
sx q[1];
rz(11.799078) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7458079) q[0];
sx q[0];
rz(-1.1757478) q[0];
sx q[0];
rz(-1.5725738) q[0];
x q[1];
rz(-2.4440304) q[2];
sx q[2];
rz(-2.5346906) q[2];
sx q[2];
rz(1.6873492) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0815107) q[1];
sx q[1];
rz(-2.7742371) q[1];
sx q[1];
rz(1.0698331) q[1];
x q[2];
rz(0.49873036) q[3];
sx q[3];
rz(-1.937003) q[3];
sx q[3];
rz(1.6134451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2904498) q[2];
sx q[2];
rz(-0.81633154) q[2];
sx q[2];
rz(1.2510703) q[2];
rz(2.3276681) q[3];
sx q[3];
rz(-1.1000752) q[3];
sx q[3];
rz(-1.6479015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55452764) q[0];
sx q[0];
rz(-1.6870512) q[0];
sx q[0];
rz(-0.54334193) q[0];
rz(2.4045565) q[1];
sx q[1];
rz(-2.4529424) q[1];
sx q[1];
rz(0.063311689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2800326) q[0];
sx q[0];
rz(-0.82516042) q[0];
sx q[0];
rz(-0.78365032) q[0];
rz(-1.8149764) q[2];
sx q[2];
rz(-1.2968213) q[2];
sx q[2];
rz(-1.8802644) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0249512) q[1];
sx q[1];
rz(-2.4519741) q[1];
sx q[1];
rz(-2.078213) q[1];
rz(1.3593986) q[3];
sx q[3];
rz(-1.4305887) q[3];
sx q[3];
rz(-2.0003776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43725499) q[2];
sx q[2];
rz(-2.0180118) q[2];
sx q[2];
rz(-2.227318) q[2];
rz(2.8546913) q[3];
sx q[3];
rz(-0.56921452) q[3];
sx q[3];
rz(-0.60681075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8617926) q[0];
sx q[0];
rz(-1.8656116) q[0];
sx q[0];
rz(1.7472501) q[0];
rz(2.2718248) q[1];
sx q[1];
rz(-0.52647796) q[1];
sx q[1];
rz(0.24551749) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95076563) q[0];
sx q[0];
rz(-0.47266911) q[0];
sx q[0];
rz(2.5373775) q[0];
rz(1.0209544) q[2];
sx q[2];
rz(-1.6022953) q[2];
sx q[2];
rz(-1.3481639) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15233065) q[1];
sx q[1];
rz(-0.50571364) q[1];
sx q[1];
rz(2.5515517) q[1];
rz(0.80733733) q[3];
sx q[3];
rz(-1.5669979) q[3];
sx q[3];
rz(-1.8769052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8292344) q[2];
sx q[2];
rz(-1.979253) q[2];
sx q[2];
rz(-1.8355628) q[2];
rz(2.8570789) q[3];
sx q[3];
rz(-1.5008711) q[3];
sx q[3];
rz(-1.7623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5259842) q[0];
sx q[0];
rz(-0.59399501) q[0];
sx q[0];
rz(-2.3748412) q[0];
rz(-1.7150257) q[1];
sx q[1];
rz(-1.7568935) q[1];
sx q[1];
rz(-1.0624622) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7257197) q[0];
sx q[0];
rz(-1.5860646) q[0];
sx q[0];
rz(1.589561) q[0];
x q[1];
rz(-3.0059848) q[2];
sx q[2];
rz(-2.1580843) q[2];
sx q[2];
rz(2.7073017) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6973138) q[1];
sx q[1];
rz(-2.0525888) q[1];
sx q[1];
rz(-2.8413111) q[1];
x q[2];
rz(-2.9963909) q[3];
sx q[3];
rz(-1.6305974) q[3];
sx q[3];
rz(-0.54403323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95395025) q[2];
sx q[2];
rz(-3.1270449) q[2];
sx q[2];
rz(0.1178096) q[2];
rz(-2.5976962) q[3];
sx q[3];
rz(-1.0313279) q[3];
sx q[3];
rz(2.4081374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.2619005) q[0];
sx q[0];
rz(-1.6629135) q[0];
sx q[0];
rz(-0.53801584) q[0];
rz(-3.0777439) q[1];
sx q[1];
rz(-1.0212746) q[1];
sx q[1];
rz(-1.8003546) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8842476) q[0];
sx q[0];
rz(-1.0228511) q[0];
sx q[0];
rz(2.1764285) q[0];
x q[1];
rz(0.55300216) q[2];
sx q[2];
rz(-1.2931839) q[2];
sx q[2];
rz(2.0763408) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2423305) q[1];
sx q[1];
rz(-0.56277187) q[1];
sx q[1];
rz(0.11426781) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7085764) q[3];
sx q[3];
rz(-1.4524432) q[3];
sx q[3];
rz(2.1089212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4591879) q[2];
sx q[2];
rz(-0.29067278) q[2];
sx q[2];
rz(2.235137) q[2];
rz(-1.9393548) q[3];
sx q[3];
rz(-1.5463444) q[3];
sx q[3];
rz(-1.7513587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9850013) q[0];
sx q[0];
rz(-2.5522794) q[0];
sx q[0];
rz(-0.44664788) q[0];
rz(-2.7062972) q[1];
sx q[1];
rz(-0.63658249) q[1];
sx q[1];
rz(2.2925503) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0627611) q[0];
sx q[0];
rz(-2.9568892) q[0];
sx q[0];
rz(-1.9901426) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6094815) q[2];
sx q[2];
rz(-1.7052884) q[2];
sx q[2];
rz(-2.7464339) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0246692) q[1];
sx q[1];
rz(-1.797442) q[1];
sx q[1];
rz(-1.5818854) q[1];
rz(-pi) q[2];
rz(1.7432415) q[3];
sx q[3];
rz(-1.9151805) q[3];
sx q[3];
rz(1.1960967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1242096) q[2];
sx q[2];
rz(-2.1653192) q[2];
sx q[2];
rz(1.3503831) q[2];
rz(-0.24556686) q[3];
sx q[3];
rz(-0.96846247) q[3];
sx q[3];
rz(0.33026162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6733112) q[0];
sx q[0];
rz(-2.249233) q[0];
sx q[0];
rz(-2.6924676) q[0];
rz(-1.4324073) q[1];
sx q[1];
rz(-0.9811554) q[1];
sx q[1];
rz(1.2264576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0104727) q[0];
sx q[0];
rz(-0.95632271) q[0];
sx q[0];
rz(-3.0486186) q[0];
x q[1];
rz(-0.50275393) q[2];
sx q[2];
rz(-1.9079676) q[2];
sx q[2];
rz(1.3664811) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.439319) q[1];
sx q[1];
rz(-1.4838184) q[1];
sx q[1];
rz(-3.0371515) q[1];
rz(-0.65761538) q[3];
sx q[3];
rz(-2.6675993) q[3];
sx q[3];
rz(-2.4331987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2655098) q[2];
sx q[2];
rz(-0.86388695) q[2];
sx q[2];
rz(-0.49501219) q[2];
rz(-1.0692976) q[3];
sx q[3];
rz(-2.4102231) q[3];
sx q[3];
rz(-0.67720145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7345562) q[0];
sx q[0];
rz(-0.41663909) q[0];
sx q[0];
rz(0.46105841) q[0];
rz(3.0126493) q[1];
sx q[1];
rz(-2.415633) q[1];
sx q[1];
rz(-2.2393548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8811097) q[0];
sx q[0];
rz(-0.04180464) q[0];
sx q[0];
rz(0.16784541) q[0];
x q[1];
rz(3.0831218) q[2];
sx q[2];
rz(-1.2035511) q[2];
sx q[2];
rz(-2.1580687) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9736147) q[1];
sx q[1];
rz(-1.2149286) q[1];
sx q[1];
rz(2.0116429) q[1];
rz(-pi) q[2];
x q[2];
rz(0.095288988) q[3];
sx q[3];
rz(-2.2248292) q[3];
sx q[3];
rz(-0.44910329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43877131) q[2];
sx q[2];
rz(-2.2192025) q[2];
sx q[2];
rz(1.927446) q[2];
rz(-2.7625648) q[3];
sx q[3];
rz(-1.4641848) q[3];
sx q[3];
rz(-1.3488784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.2269065) q[0];
sx q[0];
rz(-1.560697) q[0];
sx q[0];
rz(-3.1353986) q[0];
rz(-0.50503039) q[1];
sx q[1];
rz(-2.6513702) q[1];
sx q[1];
rz(-0.44713155) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762334) q[0];
sx q[0];
rz(-2.0885944) q[0];
sx q[0];
rz(0.8014265) q[0];
x q[1];
rz(0.51905175) q[2];
sx q[2];
rz(-1.468892) q[2];
sx q[2];
rz(2.6554012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.3382692) q[1];
sx q[1];
rz(-1.4191966) q[1];
sx q[1];
rz(0.029241745) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7617486) q[3];
sx q[3];
rz(-1.1292283) q[3];
sx q[3];
rz(-0.34363765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7827683) q[2];
sx q[2];
rz(-1.1268104) q[2];
sx q[2];
rz(-0.44754851) q[2];
rz(-1.1131845) q[3];
sx q[3];
rz(-1.0380849) q[3];
sx q[3];
rz(0.32570496) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96414763) q[0];
sx q[0];
rz(-0.064346813) q[0];
sx q[0];
rz(-2.294975) q[0];
rz(2.8586491) q[1];
sx q[1];
rz(-1.6694262) q[1];
sx q[1];
rz(2.4683594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90712522) q[0];
sx q[0];
rz(-1.5801718) q[0];
sx q[0];
rz(-3.1251291) q[0];
rz(2.908082) q[2];
sx q[2];
rz(-2.3552009) q[2];
sx q[2];
rz(-0.76372432) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9615508) q[1];
sx q[1];
rz(-2.2030239) q[1];
sx q[1];
rz(0.2747196) q[1];
rz(-2.2825234) q[3];
sx q[3];
rz(-2.7112656) q[3];
sx q[3];
rz(2.8266065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8474569) q[2];
sx q[2];
rz(-2.2723276) q[2];
sx q[2];
rz(0.67226234) q[2];
rz(1.8654035) q[3];
sx q[3];
rz(-1.0381235) q[3];
sx q[3];
rz(2.7289895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1558253) q[0];
sx q[0];
rz(-1.1341996) q[0];
sx q[0];
rz(-1.2962935) q[0];
rz(-0.014421944) q[1];
sx q[1];
rz(-1.2851234) q[1];
sx q[1];
rz(-1.8869225) q[1];
rz(1.1070743) q[2];
sx q[2];
rz(-0.92841371) q[2];
sx q[2];
rz(-2.4910891) q[2];
rz(-2.9113967) q[3];
sx q[3];
rz(-1.7741813) q[3];
sx q[3];
rz(0.81188191) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
