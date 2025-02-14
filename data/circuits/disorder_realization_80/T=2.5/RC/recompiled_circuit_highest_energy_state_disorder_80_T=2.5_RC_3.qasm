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
rz(0.12240527) q[0];
sx q[0];
rz(-0.91949099) q[0];
sx q[0];
rz(1.9424196) q[0];
rz(-2.9543258) q[1];
sx q[1];
rz(-0.54228243) q[1];
sx q[1];
rz(-1.6220925) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2080948) q[0];
sx q[0];
rz(-1.6483432) q[0];
sx q[0];
rz(-2.5258875) q[0];
x q[1];
rz(-2.7678732) q[2];
sx q[2];
rz(-1.9961832) q[2];
sx q[2];
rz(-0.56217867) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6295084) q[1];
sx q[1];
rz(-2.0633882) q[1];
sx q[1];
rz(-2.0297445) q[1];
rz(-pi) q[2];
rz(2.8974122) q[3];
sx q[3];
rz(-1.5236519) q[3];
sx q[3];
rz(-2.3650996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2970994) q[2];
sx q[2];
rz(-0.62701925) q[2];
sx q[2];
rz(1.373488) q[2];
rz(-2.3376236) q[3];
sx q[3];
rz(-1.4990025) q[3];
sx q[3];
rz(-1.1871626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3961769) q[0];
sx q[0];
rz(-0.71393037) q[0];
sx q[0];
rz(-2.7667238) q[0];
rz(-2.6238341) q[1];
sx q[1];
rz(-1.2034143) q[1];
sx q[1];
rz(1.7727324) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33297172) q[0];
sx q[0];
rz(-0.00095168984) q[0];
sx q[0];
rz(3.0406221) q[0];
rz(-pi) q[1];
rz(1.5438186) q[2];
sx q[2];
rz(-1.7671874) q[2];
sx q[2];
rz(-0.90919221) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79530935) q[1];
sx q[1];
rz(-2.4716271) q[1];
sx q[1];
rz(0.69366347) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1325466) q[3];
sx q[3];
rz(-0.25190946) q[3];
sx q[3];
rz(2.589165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.123473) q[2];
sx q[2];
rz(-0.70088434) q[2];
sx q[2];
rz(-2.4269721) q[2];
rz(-2.9361652) q[3];
sx q[3];
rz(-1.8193865) q[3];
sx q[3];
rz(-1.8429168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4648723) q[0];
sx q[0];
rz(-1.0370075) q[0];
sx q[0];
rz(2.6439457) q[0];
rz(2.5045555) q[1];
sx q[1];
rz(-0.62890816) q[1];
sx q[1];
rz(-0.10674891) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7926559) q[0];
sx q[0];
rz(-1.7829527) q[0];
sx q[0];
rz(1.2646219) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90107839) q[2];
sx q[2];
rz(-1.5656098) q[2];
sx q[2];
rz(2.5825495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.63697228) q[1];
sx q[1];
rz(-0.85831888) q[1];
sx q[1];
rz(-0.070657151) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40092881) q[3];
sx q[3];
rz(-1.2967921) q[3];
sx q[3];
rz(-1.2284492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5367624) q[2];
sx q[2];
rz(-0.75211516) q[2];
sx q[2];
rz(-0.19482782) q[2];
rz(0.005006494) q[3];
sx q[3];
rz(-2.5275793) q[3];
sx q[3];
rz(0.79202882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1057338) q[0];
sx q[0];
rz(-0.53426131) q[0];
sx q[0];
rz(2.1015097) q[0];
rz(0.19373521) q[1];
sx q[1];
rz(-1.6400784) q[1];
sx q[1];
rz(2.3949882) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23768856) q[0];
sx q[0];
rz(-2.6857566) q[0];
sx q[0];
rz(1.8346915) q[0];
rz(-pi) q[1];
rz(-0.14239133) q[2];
sx q[2];
rz(-2.2658544) q[2];
sx q[2];
rz(-0.54706878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6317217) q[1];
sx q[1];
rz(-1.6225909) q[1];
sx q[1];
rz(-2.447261) q[1];
x q[2];
rz(2.901731) q[3];
sx q[3];
rz(-1.4647533) q[3];
sx q[3];
rz(2.9221846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.99183434) q[2];
sx q[2];
rz(-1.6712302) q[2];
sx q[2];
rz(0.20450083) q[2];
rz(1.9483942) q[3];
sx q[3];
rz(-2.2843993) q[3];
sx q[3];
rz(0.96113718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0123154) q[0];
sx q[0];
rz(-2.9670872) q[0];
sx q[0];
rz(-2.0611064) q[0];
rz(0.37733817) q[1];
sx q[1];
rz(-1.6308547) q[1];
sx q[1];
rz(0.89471716) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6301439) q[0];
sx q[0];
rz(-2.7155657) q[0];
sx q[0];
rz(-0.57630868) q[0];
rz(2.3482125) q[2];
sx q[2];
rz(-1.0715683) q[2];
sx q[2];
rz(-2.3069059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1354358) q[1];
sx q[1];
rz(-0.73169152) q[1];
sx q[1];
rz(-0.77349551) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36578806) q[3];
sx q[3];
rz(-0.96082965) q[3];
sx q[3];
rz(0.47780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6696024) q[2];
sx q[2];
rz(-2.695638) q[2];
sx q[2];
rz(-3.0338244) q[2];
rz(2.9660411) q[3];
sx q[3];
rz(-1.8864417) q[3];
sx q[3];
rz(-0.56041437) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948697) q[0];
sx q[0];
rz(-2.6304498) q[0];
sx q[0];
rz(2.129659) q[0];
rz(1.5049505) q[1];
sx q[1];
rz(-2.4220146) q[1];
sx q[1];
rz(-1.0171657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6300679) q[0];
sx q[0];
rz(-1.2753133) q[0];
sx q[0];
rz(0.17952551) q[0];
rz(-0.58508137) q[2];
sx q[2];
rz(-0.85137109) q[2];
sx q[2];
rz(2.35308) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6222657) q[1];
sx q[1];
rz(-2.9510088) q[1];
sx q[1];
rz(-3.0960073) q[1];
rz(-pi) q[2];
rz(-1.1215707) q[3];
sx q[3];
rz(-1.4288846) q[3];
sx q[3];
rz(-0.52983879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4096058) q[2];
sx q[2];
rz(-0.66183949) q[2];
sx q[2];
rz(-0.99676639) q[2];
rz(-3.0766727) q[3];
sx q[3];
rz(-1.4109979) q[3];
sx q[3];
rz(-1.9916649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053059269) q[0];
sx q[0];
rz(-2.2011338) q[0];
sx q[0];
rz(-2.233182) q[0];
rz(2.9216595) q[1];
sx q[1];
rz(-1.5426153) q[1];
sx q[1];
rz(1.251108) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.096643) q[0];
sx q[0];
rz(-1.8306376) q[0];
sx q[0];
rz(1.0837062) q[0];
rz(-2.8562137) q[2];
sx q[2];
rz(-1.4754461) q[2];
sx q[2];
rz(2.2298663) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0798222) q[1];
sx q[1];
rz(-0.95937371) q[1];
sx q[1];
rz(0.76217125) q[1];
rz(-pi) q[2];
rz(1.7102107) q[3];
sx q[3];
rz(-2.406139) q[3];
sx q[3];
rz(-2.9289587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1088341) q[2];
sx q[2];
rz(-2.3075576) q[2];
sx q[2];
rz(1.6081107) q[2];
rz(-1.1533302) q[3];
sx q[3];
rz(-1.0554375) q[3];
sx q[3];
rz(1.6495033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.15025) q[0];
sx q[0];
rz(-0.38337502) q[0];
sx q[0];
rz(-0.94069329) q[0];
rz(0.13537814) q[1];
sx q[1];
rz(-1.7372513) q[1];
sx q[1];
rz(2.1167596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8405466) q[0];
sx q[0];
rz(-2.4468166) q[0];
sx q[0];
rz(-1.077681) q[0];
rz(2.5338855) q[2];
sx q[2];
rz(-1.6204699) q[2];
sx q[2];
rz(0.36591727) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8044628) q[1];
sx q[1];
rz(-1.4373597) q[1];
sx q[1];
rz(0.19652469) q[1];
x q[2];
rz(1.5190184) q[3];
sx q[3];
rz(-1.9815784) q[3];
sx q[3];
rz(0.67922986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3465053) q[2];
sx q[2];
rz(-2.2515209) q[2];
sx q[2];
rz(0.65004641) q[2];
rz(2.5551689) q[3];
sx q[3];
rz(-1.4805877) q[3];
sx q[3];
rz(2.39095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8745678) q[0];
sx q[0];
rz(-0.50849193) q[0];
sx q[0];
rz(1.1453999) q[0];
rz(-1.3151431) q[1];
sx q[1];
rz(-1.2329085) q[1];
sx q[1];
rz(-2.4536536) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8696339) q[0];
sx q[0];
rz(-0.42068538) q[0];
sx q[0];
rz(-1.0933594) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0264977) q[2];
sx q[2];
rz(-0.67152464) q[2];
sx q[2];
rz(-1.4948927) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.1286182) q[1];
sx q[1];
rz(-0.77207295) q[1];
sx q[1];
rz(2.9086604) q[1];
rz(2.7642058) q[3];
sx q[3];
rz(-2.0321369) q[3];
sx q[3];
rz(2.2113706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73100662) q[2];
sx q[2];
rz(-1.1062016) q[2];
sx q[2];
rz(0.86622396) q[2];
rz(-2.2335562) q[3];
sx q[3];
rz(-2.3767545) q[3];
sx q[3];
rz(-1.0959371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9662125) q[0];
sx q[0];
rz(-1.1840273) q[0];
sx q[0];
rz(0.41561919) q[0];
rz(-0.79187727) q[1];
sx q[1];
rz(-1.917058) q[1];
sx q[1];
rz(0.22722879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2450175) q[0];
sx q[0];
rz(-1.554721) q[0];
sx q[0];
rz(1.5825558) q[0];
rz(-1.6516901) q[2];
sx q[2];
rz(-2.3757907) q[2];
sx q[2];
rz(1.7532938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1243678) q[1];
sx q[1];
rz(-1.3137046) q[1];
sx q[1];
rz(0.59558792) q[1];
rz(-pi) q[2];
rz(-2.6260904) q[3];
sx q[3];
rz(-1.1519264) q[3];
sx q[3];
rz(2.2282698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8772584) q[2];
sx q[2];
rz(-2.4206968) q[2];
sx q[2];
rz(-0.43992511) q[2];
rz(2.544493) q[3];
sx q[3];
rz(-0.66949451) q[3];
sx q[3];
rz(1.7841024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6482342) q[0];
sx q[0];
rz(-1.5879205) q[0];
sx q[0];
rz(1.2863202) q[0];
rz(-1.413912) q[1];
sx q[1];
rz(-1.021011) q[1];
sx q[1];
rz(0.11722142) q[1];
rz(-1.6261423) q[2];
sx q[2];
rz(-1.9774441) q[2];
sx q[2];
rz(-1.0033506) q[2];
rz(-1.8341919) q[3];
sx q[3];
rz(-0.7701244) q[3];
sx q[3];
rz(-0.27944892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
