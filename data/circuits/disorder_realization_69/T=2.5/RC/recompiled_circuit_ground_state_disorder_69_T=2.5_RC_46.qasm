OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90836877) q[0];
sx q[0];
rz(4.6908661) q[0];
sx q[0];
rz(9.6022845) q[0];
rz(-1.5818051) q[1];
sx q[1];
rz(-1.3270562) q[1];
sx q[1];
rz(1.1787193) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4807199) q[0];
sx q[0];
rz(-0.62107039) q[0];
sx q[0];
rz(1.8939411) q[0];
rz(-pi) q[1];
rz(0.29961961) q[2];
sx q[2];
rz(-2.1914406) q[2];
sx q[2];
rz(-0.52052467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.385966) q[1];
sx q[1];
rz(-0.75191359) q[1];
sx q[1];
rz(0.90689838) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3050582) q[3];
sx q[3];
rz(-1.4713376) q[3];
sx q[3];
rz(0.11900706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.061964758) q[2];
sx q[2];
rz(-0.92755932) q[2];
sx q[2];
rz(0.24688841) q[2];
rz(-2.7039792) q[3];
sx q[3];
rz(-2.103431) q[3];
sx q[3];
rz(-0.26052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0831809) q[0];
sx q[0];
rz(-2.658168) q[0];
sx q[0];
rz(-1.1046945) q[0];
rz(2.3339234) q[1];
sx q[1];
rz(-0.75391114) q[1];
sx q[1];
rz(-2.5025936) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4527567) q[0];
sx q[0];
rz(-2.8114744) q[0];
sx q[0];
rz(-2.5725276) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17484395) q[2];
sx q[2];
rz(-2.8626466) q[2];
sx q[2];
rz(2.1871559) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39589992) q[1];
sx q[1];
rz(-1.9338765) q[1];
sx q[1];
rz(2.3895742) q[1];
rz(2.8547332) q[3];
sx q[3];
rz(-1.8897227) q[3];
sx q[3];
rz(-0.87251679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4421926) q[2];
sx q[2];
rz(-1.7510119) q[2];
sx q[2];
rz(2.3244582) q[2];
rz(-2.6160431) q[3];
sx q[3];
rz(-0.81997973) q[3];
sx q[3];
rz(-1.1901963) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84080559) q[0];
sx q[0];
rz(-1.5504993) q[0];
sx q[0];
rz(0.15980414) q[0];
rz(0.5197168) q[1];
sx q[1];
rz(-1.0451885) q[1];
sx q[1];
rz(-2.2249178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778939) q[0];
sx q[0];
rz(-0.5710154) q[0];
sx q[0];
rz(-3.0932941) q[0];
x q[1];
rz(2.8755912) q[2];
sx q[2];
rz(-2.383138) q[2];
sx q[2];
rz(-0.11125362) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0708051) q[1];
sx q[1];
rz(-0.19529058) q[1];
sx q[1];
rz(-1.3210123) q[1];
rz(-pi) q[2];
rz(2.9514922) q[3];
sx q[3];
rz(-2.5203966) q[3];
sx q[3];
rz(-2.7060696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2436169) q[2];
sx q[2];
rz(-2.6349082) q[2];
sx q[2];
rz(1.766073) q[2];
rz(-3.1084133) q[3];
sx q[3];
rz(-1.5539919) q[3];
sx q[3];
rz(-0.83056617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0631977) q[0];
sx q[0];
rz(-0.66135186) q[0];
sx q[0];
rz(0.98264328) q[0];
rz(-2.4962418) q[1];
sx q[1];
rz(-1.2256365) q[1];
sx q[1];
rz(-0.4172079) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44051927) q[0];
sx q[0];
rz(-0.60581453) q[0];
sx q[0];
rz(-2.5941501) q[0];
rz(1.190416) q[2];
sx q[2];
rz(-1.5065168) q[2];
sx q[2];
rz(0.96997875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6748089) q[1];
sx q[1];
rz(-0.67608139) q[1];
sx q[1];
rz(-0.095211857) q[1];
rz(-0.54822918) q[3];
sx q[3];
rz(-1.8129725) q[3];
sx q[3];
rz(3.0088923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8365606) q[2];
sx q[2];
rz(-2.8053668) q[2];
sx q[2];
rz(1.9545371) q[2];
rz(0.120397) q[3];
sx q[3];
rz(-1.7359366) q[3];
sx q[3];
rz(-2.9482237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3099986) q[0];
sx q[0];
rz(-0.11045063) q[0];
sx q[0];
rz(2.4647392) q[0];
rz(1.5325158) q[1];
sx q[1];
rz(-1.1108421) q[1];
sx q[1];
rz(-0.19930509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88955697) q[0];
sx q[0];
rz(-1.8235208) q[0];
sx q[0];
rz(-3.0480372) q[0];
x q[1];
rz(-1.4485613) q[2];
sx q[2];
rz(-1.9150371) q[2];
sx q[2];
rz(2.2064672) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78969722) q[1];
sx q[1];
rz(-1.1994065) q[1];
sx q[1];
rz(1.2523734) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0065157117) q[3];
sx q[3];
rz(-1.2938936) q[3];
sx q[3];
rz(-2.2369825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6469595) q[2];
sx q[2];
rz(-2.5854752) q[2];
sx q[2];
rz(2.6456918) q[2];
rz(2.2206709) q[3];
sx q[3];
rz(-2.0995188) q[3];
sx q[3];
rz(0.31057772) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344236) q[0];
sx q[0];
rz(-1.9155707) q[0];
sx q[0];
rz(2.7219211) q[0];
rz(-0.93027973) q[1];
sx q[1];
rz(-1.0134965) q[1];
sx q[1];
rz(1.0118265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8735805) q[0];
sx q[0];
rz(-2.4873864) q[0];
sx q[0];
rz(3.0604721) q[0];
rz(-0.65134766) q[2];
sx q[2];
rz(-2.911951) q[2];
sx q[2];
rz(2.1852213) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4395482) q[1];
sx q[1];
rz(-1.3633894) q[1];
sx q[1];
rz(-1.009936) q[1];
x q[2];
rz(-2.9568372) q[3];
sx q[3];
rz(-2.1596381) q[3];
sx q[3];
rz(-2.7240804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4100627) q[2];
sx q[2];
rz(-2.5750934) q[2];
sx q[2];
rz(-1.1140964) q[2];
rz(-1.3725494) q[3];
sx q[3];
rz(-2.1362344) q[3];
sx q[3];
rz(-0.70926386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0396742) q[0];
sx q[0];
rz(-1.2693951) q[0];
sx q[0];
rz(1.6495548) q[0];
rz(2.6986625) q[1];
sx q[1];
rz(-1.8699173) q[1];
sx q[1];
rz(2.2992004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32552606) q[0];
sx q[0];
rz(-2.2057945) q[0];
sx q[0];
rz(2.2731337) q[0];
x q[1];
rz(1.2356495) q[2];
sx q[2];
rz(-2.4494) q[2];
sx q[2];
rz(0.38470925) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29553005) q[1];
sx q[1];
rz(-1.1605442) q[1];
sx q[1];
rz(-2.2890291) q[1];
rz(1.674323) q[3];
sx q[3];
rz(-1.5281505) q[3];
sx q[3];
rz(-1.9089684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3288021) q[2];
sx q[2];
rz(-2.2219358) q[2];
sx q[2];
rz(3.1205422) q[2];
rz(2.1444495) q[3];
sx q[3];
rz(-2.5992584) q[3];
sx q[3];
rz(-1.5890582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.3461304) q[0];
sx q[0];
rz(-1.2393476) q[0];
sx q[0];
rz(-1.963266) q[0];
rz(1.0109674) q[1];
sx q[1];
rz(-1.4821472) q[1];
sx q[1];
rz(2.7706026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31974935) q[0];
sx q[0];
rz(-2.8822063) q[0];
sx q[0];
rz(0.46617561) q[0];
x q[1];
rz(1.0576789) q[2];
sx q[2];
rz(-2.9484595) q[2];
sx q[2];
rz(2.4095175) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.64196649) q[1];
sx q[1];
rz(-1.2807944) q[1];
sx q[1];
rz(0.15695842) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1493672) q[3];
sx q[3];
rz(-2.605509) q[3];
sx q[3];
rz(1.9158196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79169881) q[2];
sx q[2];
rz(-2.0588304) q[2];
sx q[2];
rz(0.91378158) q[2];
rz(0.53798109) q[3];
sx q[3];
rz(-0.83298433) q[3];
sx q[3];
rz(-0.37916455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0898042) q[0];
sx q[0];
rz(-1.3642949) q[0];
sx q[0];
rz(-2.9096933) q[0];
rz(1.6357577) q[1];
sx q[1];
rz(-1.6927405) q[1];
sx q[1];
rz(-0.64461446) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8763257) q[0];
sx q[0];
rz(-0.4943706) q[0];
sx q[0];
rz(-2.3571924) q[0];
rz(3.1050131) q[2];
sx q[2];
rz(-2.5365005) q[2];
sx q[2];
rz(0.99672752) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7773103) q[1];
sx q[1];
rz(-2.9743618) q[1];
sx q[1];
rz(-1.4521397) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0057974) q[3];
sx q[3];
rz(-0.91640546) q[3];
sx q[3];
rz(-1.9616579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62442786) q[2];
sx q[2];
rz(-1.6021873) q[2];
sx q[2];
rz(-1.8713162) q[2];
rz(-0.43954784) q[3];
sx q[3];
rz(-0.63616532) q[3];
sx q[3];
rz(-0.51072454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21580932) q[0];
sx q[0];
rz(-1.9010811) q[0];
sx q[0];
rz(1.0260169) q[0];
rz(-0.45694524) q[1];
sx q[1];
rz(-0.89070717) q[1];
sx q[1];
rz(0.93422008) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4480058) q[0];
sx q[0];
rz(-1.7035055) q[0];
sx q[0];
rz(-0.42610077) q[0];
rz(-pi) q[1];
x q[1];
rz(1.825594) q[2];
sx q[2];
rz(-1.464286) q[2];
sx q[2];
rz(-0.62588706) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10623744) q[1];
sx q[1];
rz(-0.83939531) q[1];
sx q[1];
rz(1.7994355) q[1];
x q[2];
rz(2.5601897) q[3];
sx q[3];
rz(-2.1384967) q[3];
sx q[3];
rz(-2.3025177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43941867) q[2];
sx q[2];
rz(-0.69434035) q[2];
sx q[2];
rz(1.0609421) q[2];
rz(2.151978) q[3];
sx q[3];
rz(-1.3720392) q[3];
sx q[3];
rz(2.8770679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92985741) q[0];
sx q[0];
rz(-1.095357) q[0];
sx q[0];
rz(2.3768421) q[0];
rz(1.1872956) q[1];
sx q[1];
rz(-1.7056414) q[1];
sx q[1];
rz(2.5797226) q[1];
rz(1.622772) q[2];
sx q[2];
rz(-1.5904273) q[2];
sx q[2];
rz(-1.2640819) q[2];
rz(1.1239617) q[3];
sx q[3];
rz(-1.1458967) q[3];
sx q[3];
rz(0.55895373) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
