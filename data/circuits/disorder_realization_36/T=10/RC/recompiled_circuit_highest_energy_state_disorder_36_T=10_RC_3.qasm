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
rz(-1.2688426) q[1];
sx q[1];
rz(2.3743) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7458079) q[0];
sx q[0];
rz(-1.1757478) q[0];
sx q[0];
rz(-1.5725738) q[0];
x q[1];
rz(0.69756223) q[2];
sx q[2];
rz(-2.5346906) q[2];
sx q[2];
rz(1.6873492) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5509294) q[1];
sx q[1];
rz(-1.8912705) q[1];
sx q[1];
rz(-2.9588353) q[1];
rz(-0.49873036) q[3];
sx q[3];
rz(-1.937003) q[3];
sx q[3];
rz(-1.6134451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2904498) q[2];
sx q[2];
rz(-2.3252611) q[2];
sx q[2];
rz(-1.2510703) q[2];
rz(0.81392455) q[3];
sx q[3];
rz(-1.1000752) q[3];
sx q[3];
rz(-1.4936911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.587065) q[0];
sx q[0];
rz(-1.4545414) q[0];
sx q[0];
rz(2.5982507) q[0];
rz(-0.73703611) q[1];
sx q[1];
rz(-0.68865028) q[1];
sx q[1];
rz(3.078281) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.30697) q[0];
sx q[0];
rz(-1.0234912) q[0];
sx q[0];
rz(-2.2234248) q[0];
rz(-2.43119) q[2];
sx q[2];
rz(-2.7766529) q[2];
sx q[2];
rz(1.1360363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1406537) q[1];
sx q[1];
rz(-1.88511) q[1];
sx q[1];
rz(-2.1953366) q[1];
rz(-2.1629093) q[3];
sx q[3];
rz(-0.25308713) q[3];
sx q[3];
rz(-1.0067948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7043377) q[2];
sx q[2];
rz(-2.0180118) q[2];
sx q[2];
rz(-0.91427461) q[2];
rz(-0.28690139) q[3];
sx q[3];
rz(-0.56921452) q[3];
sx q[3];
rz(-0.60681075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27980003) q[0];
sx q[0];
rz(-1.2759811) q[0];
sx q[0];
rz(1.7472501) q[0];
rz(-2.2718248) q[1];
sx q[1];
rz(-2.6151147) q[1];
sx q[1];
rz(-2.8960752) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29123339) q[0];
sx q[0];
rz(-1.1867674) q[0];
sx q[0];
rz(1.8535094) q[0];
rz(-1.5105702) q[2];
sx q[2];
rz(-0.55065075) q[2];
sx q[2];
rz(-2.8676195) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.88855381) q[1];
sx q[1];
rz(-1.2978857) q[1];
sx q[1];
rz(2.7103579) q[1];
rz(-pi) q[2];
rz(-2.3342553) q[3];
sx q[3];
rz(-1.5669979) q[3];
sx q[3];
rz(-1.8769052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3123582) q[2];
sx q[2];
rz(-1.1623397) q[2];
sx q[2];
rz(1.8355628) q[2];
rz(-0.28451377) q[3];
sx q[3];
rz(-1.6407216) q[3];
sx q[3];
rz(-1.3792926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5259842) q[0];
sx q[0];
rz(-2.5475976) q[0];
sx q[0];
rz(2.3748412) q[0];
rz(-1.7150257) q[1];
sx q[1];
rz(-1.7568935) q[1];
sx q[1];
rz(-1.0624622) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4158729) q[0];
sx q[0];
rz(-1.5860646) q[0];
sx q[0];
rz(-1.589561) q[0];
rz(1.771174) q[2];
sx q[2];
rz(-0.60094072) q[2];
sx q[2];
rz(-0.19285017) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8725682) q[1];
sx q[1];
rz(-1.3055798) q[1];
sx q[1];
rz(1.069963) q[1];
x q[2];
rz(-2.7492529) q[3];
sx q[3];
rz(-2.9846387) q[3];
sx q[3];
rz(0.63877393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95395025) q[2];
sx q[2];
rz(-0.01454777) q[2];
sx q[2];
rz(3.0237831) q[2];
rz(0.5438965) q[3];
sx q[3];
rz(-2.1102648) q[3];
sx q[3];
rz(-2.4081374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8796922) q[0];
sx q[0];
rz(-1.6629135) q[0];
sx q[0];
rz(-2.6035768) q[0];
rz(-0.063848786) q[1];
sx q[1];
rz(-1.0212746) q[1];
sx q[1];
rz(-1.341238) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1088561) q[0];
sx q[0];
rz(-2.0782315) q[0];
sx q[0];
rz(-2.5030337) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55300216) q[2];
sx q[2];
rz(-1.2931839) q[2];
sx q[2];
rz(-1.0652519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3772014) q[1];
sx q[1];
rz(-2.1294597) q[1];
sx q[1];
rz(-1.6425981) q[1];
rz(1.7010541) q[3];
sx q[3];
rz(-2.0005811) q[3];
sx q[3];
rz(0.48359475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4591879) q[2];
sx q[2];
rz(-0.29067278) q[2];
sx q[2];
rz(0.9064557) q[2];
rz(-1.9393548) q[3];
sx q[3];
rz(-1.5463444) q[3];
sx q[3];
rz(1.3902339) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9850013) q[0];
sx q[0];
rz(-2.5522794) q[0];
sx q[0];
rz(0.44664788) q[0];
rz(2.7062972) q[1];
sx q[1];
rz(-0.63658249) q[1];
sx q[1];
rz(-2.2925503) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.236625) q[0];
sx q[0];
rz(-1.4959488) q[0];
sx q[0];
rz(-1.4017795) q[0];
rz(2.6094815) q[2];
sx q[2];
rz(-1.7052884) q[2];
sx q[2];
rz(2.7464339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0246692) q[1];
sx q[1];
rz(-1.797442) q[1];
sx q[1];
rz(1.5818854) q[1];
rz(-pi) q[2];
rz(2.6953727) q[3];
sx q[3];
rz(-0.38360197) q[3];
sx q[3];
rz(-1.6723796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1242096) q[2];
sx q[2];
rz(-2.1653192) q[2];
sx q[2];
rz(-1.7912095) q[2];
rz(-2.8960258) q[3];
sx q[3];
rz(-2.1731302) q[3];
sx q[3];
rz(-2.811331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6733112) q[0];
sx q[0];
rz(-0.89235965) q[0];
sx q[0];
rz(-0.44912502) q[0];
rz(1.7091854) q[1];
sx q[1];
rz(-0.9811554) q[1];
sx q[1];
rz(1.2264576) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2914631) q[0];
sx q[0];
rz(-0.62056834) q[0];
sx q[0];
rz(-1.4399714) q[0];
x q[1];
rz(-1.9513556) q[2];
sx q[2];
rz(-2.042843) q[2];
sx q[2];
rz(-2.7573331) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82348541) q[1];
sx q[1];
rz(-0.13581443) q[1];
sx q[1];
rz(-0.69655711) q[1];
x q[2];
rz(1.2669446) q[3];
sx q[3];
rz(-1.2011853) q[3];
sx q[3];
rz(-0.0064484869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2655098) q[2];
sx q[2];
rz(-2.2777057) q[2];
sx q[2];
rz(0.49501219) q[2];
rz(2.072295) q[3];
sx q[3];
rz(-0.73136955) q[3];
sx q[3];
rz(0.67720145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7345562) q[0];
sx q[0];
rz(-0.41663909) q[0];
sx q[0];
rz(-0.46105841) q[0];
rz(-3.0126493) q[1];
sx q[1];
rz(-0.72595969) q[1];
sx q[1];
rz(-2.2393548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26048294) q[0];
sx q[0];
rz(-3.099788) q[0];
sx q[0];
rz(-2.9737472) q[0];
rz(-pi) q[1];
rz(1.2029776) q[2];
sx q[2];
rz(-1.6253643) q[2];
sx q[2];
rz(-0.56625783) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.56574351) q[1];
sx q[1];
rz(-1.9822803) q[1];
sx q[1];
rz(-0.38994148) q[1];
rz(-pi) q[2];
x q[2];
rz(0.095288988) q[3];
sx q[3];
rz(-2.2248292) q[3];
sx q[3];
rz(2.6924894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43877131) q[2];
sx q[2];
rz(-2.2192025) q[2];
sx q[2];
rz(-1.2141466) q[2];
rz(0.37902784) q[3];
sx q[3];
rz(-1.4641848) q[3];
sx q[3];
rz(-1.3488784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2269065) q[0];
sx q[0];
rz(-1.5808957) q[0];
sx q[0];
rz(3.1353986) q[0];
rz(2.6365623) q[1];
sx q[1];
rz(-0.49022245) q[1];
sx q[1];
rz(0.44713155) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4775708) q[0];
sx q[0];
rz(-2.24488) q[0];
sx q[0];
rz(0.88468219) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51905175) q[2];
sx q[2];
rz(-1.468892) q[2];
sx q[2];
rz(-2.6554012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8033235) q[1];
sx q[1];
rz(-1.4191966) q[1];
sx q[1];
rz(0.029241745) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7617486) q[3];
sx q[3];
rz(-2.0123643) q[3];
sx q[3];
rz(-2.797955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7827683) q[2];
sx q[2];
rz(-1.1268104) q[2];
sx q[2];
rz(2.6940441) q[2];
rz(-2.0284082) q[3];
sx q[3];
rz(-2.1035078) q[3];
sx q[3];
rz(-2.8158877) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.177445) q[0];
sx q[0];
rz(-0.064346813) q[0];
sx q[0];
rz(-0.84661761) q[0];
rz(-0.28294357) q[1];
sx q[1];
rz(-1.6694262) q[1];
sx q[1];
rz(2.4683594) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780759) q[0];
sx q[0];
rz(-1.5543335) q[0];
sx q[0];
rz(-1.580173) q[0];
rz(-pi) q[1];
rz(0.77263562) q[2];
sx q[2];
rz(-1.4062721) q[2];
sx q[2];
rz(0.64060894) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5557944) q[1];
sx q[1];
rz(-1.7914247) q[1];
sx q[1];
rz(0.92024721) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8503162) q[3];
sx q[3];
rz(-1.8921953) q[3];
sx q[3];
rz(1.074256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29413572) q[2];
sx q[2];
rz(-0.86926502) q[2];
sx q[2];
rz(-0.67226234) q[2];
rz(-1.8654035) q[3];
sx q[3];
rz(-1.0381235) q[3];
sx q[3];
rz(-2.7289895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1558253) q[0];
sx q[0];
rz(-2.0073931) q[0];
sx q[0];
rz(1.8452992) q[0];
rz(-3.1271707) q[1];
sx q[1];
rz(-1.8564693) q[1];
sx q[1];
rz(1.2546702) q[1];
rz(-0.6966656) q[2];
sx q[2];
rz(-1.9370542) q[2];
sx q[2];
rz(2.5123971) q[2];
rz(2.4066385) q[3];
sx q[3];
rz(-0.30597449) q[3];
sx q[3];
rz(1.6713175) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
