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
rz(-1.1779397) q[0];
sx q[0];
rz(-3.0484634) q[0];
sx q[0];
rz(2.4094474) q[0];
rz(-1.5692476) q[1];
sx q[1];
rz(-0.88580004) q[1];
sx q[1];
rz(-2.7028309) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6308545) q[0];
sx q[0];
rz(-2.8048212) q[0];
sx q[0];
rz(1.4022409) q[0];
x q[1];
rz(-2.311539) q[2];
sx q[2];
rz(-2.7657653) q[2];
sx q[2];
rz(-2.2245882) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0369934) q[1];
sx q[1];
rz(-1.5244836) q[1];
sx q[1];
rz(0.258567) q[1];
rz(3.0436727) q[3];
sx q[3];
rz(-0.46028462) q[3];
sx q[3];
rz(-1.8079881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0438805) q[2];
sx q[2];
rz(-0.040412929) q[2];
sx q[2];
rz(0.91862339) q[2];
rz(2.0888603) q[3];
sx q[3];
rz(-2.2248) q[3];
sx q[3];
rz(-2.2241433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265428) q[0];
sx q[0];
rz(-0.81887236) q[0];
sx q[0];
rz(-0.87232605) q[0];
rz(1.0126975) q[1];
sx q[1];
rz(-1.5112977) q[1];
sx q[1];
rz(0.33908078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0377308) q[0];
sx q[0];
rz(-1.6997689) q[0];
sx q[0];
rz(0.43855389) q[0];
rz(-pi) q[1];
rz(-2.0929957) q[2];
sx q[2];
rz(-1.1638008) q[2];
sx q[2];
rz(-2.5841449) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87552947) q[1];
sx q[1];
rz(-1.5782981) q[1];
sx q[1];
rz(-0.041635978) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9267155) q[3];
sx q[3];
rz(-0.3444852) q[3];
sx q[3];
rz(0.075862715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8104711) q[2];
sx q[2];
rz(-1.1473742) q[2];
sx q[2];
rz(2.5362711) q[2];
rz(0.33453861) q[3];
sx q[3];
rz(-2.7693222) q[3];
sx q[3];
rz(-0.076586671) q[3];
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
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89638585) q[0];
sx q[0];
rz(-0.67436445) q[0];
sx q[0];
rz(1.4728004) q[0];
rz(0.32707602) q[1];
sx q[1];
rz(-1.8388803) q[1];
sx q[1];
rz(2.8715141) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.208858) q[0];
sx q[0];
rz(-1.3737824) q[0];
sx q[0];
rz(2.6126562) q[0];
x q[1];
rz(1.4065341) q[2];
sx q[2];
rz(-2.7077419) q[2];
sx q[2];
rz(-2.1151154) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9880319) q[1];
sx q[1];
rz(-0.84756339) q[1];
sx q[1];
rz(-0.60232298) q[1];
rz(-1.5746501) q[3];
sx q[3];
rz(-1.650963) q[3];
sx q[3];
rz(1.4269258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0373056) q[2];
sx q[2];
rz(-1.1423926) q[2];
sx q[2];
rz(1.2874862) q[2];
rz(1.2152952) q[3];
sx q[3];
rz(-1.7561965) q[3];
sx q[3];
rz(-2.9638885) q[3];
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
rz(0.52008587) q[0];
sx q[0];
rz(-0.44317133) q[0];
sx q[0];
rz(-0.31975123) q[0];
rz(-2.7301835) q[1];
sx q[1];
rz(-1.5450059) q[1];
sx q[1];
rz(-3.06126) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4302931) q[0];
sx q[0];
rz(-1.3724494) q[0];
sx q[0];
rz(-1.2567149) q[0];
rz(0.80647237) q[2];
sx q[2];
rz(-1.4676759) q[2];
sx q[2];
rz(-1.1191302) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6957533) q[1];
sx q[1];
rz(-1.8421122) q[1];
sx q[1];
rz(-0.81636565) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37253054) q[3];
sx q[3];
rz(-0.64465678) q[3];
sx q[3];
rz(-2.9304867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2692928) q[2];
sx q[2];
rz(-0.23739693) q[2];
sx q[2];
rz(0.044895127) q[2];
rz(0.53523713) q[3];
sx q[3];
rz(-1.5072631) q[3];
sx q[3];
rz(0.11307344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.991268) q[0];
sx q[0];
rz(-0.70882216) q[0];
sx q[0];
rz(-0.87565652) q[0];
rz(2.9278897) q[1];
sx q[1];
rz(-0.49476606) q[1];
sx q[1];
rz(1.5692086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.101709) q[0];
sx q[0];
rz(-1.3613762) q[0];
sx q[0];
rz(1.2509734) q[0];
rz(-pi) q[1];
rz(-0.46342586) q[2];
sx q[2];
rz(-2.6604386) q[2];
sx q[2];
rz(-2.0068741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24092785) q[1];
sx q[1];
rz(-2.2655409) q[1];
sx q[1];
rz(-1.7068295) q[1];
rz(-2.8646889) q[3];
sx q[3];
rz(-0.33280643) q[3];
sx q[3];
rz(1.8148418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.18465061) q[2];
sx q[2];
rz(-2.9968379) q[2];
sx q[2];
rz(-2.0733898) q[2];
rz(0.60643658) q[3];
sx q[3];
rz(-1.2111827) q[3];
sx q[3];
rz(-0.010312168) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73685) q[0];
sx q[0];
rz(-2.8917942) q[0];
sx q[0];
rz(-1.3909719) q[0];
rz(-0.5237611) q[1];
sx q[1];
rz(-0.57558376) q[1];
sx q[1];
rz(-1.1837122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9212657) q[0];
sx q[0];
rz(-1.282721) q[0];
sx q[0];
rz(0.057805268) q[0];
rz(-pi) q[1];
rz(-1.4951811) q[2];
sx q[2];
rz(-1.6241777) q[2];
sx q[2];
rz(0.23301197) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1866926) q[1];
sx q[1];
rz(-0.96265332) q[1];
sx q[1];
rz(-0.18540877) q[1];
rz(-pi) q[2];
rz(3.0314702) q[3];
sx q[3];
rz(-1.4372794) q[3];
sx q[3];
rz(0.27329521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4411053) q[2];
sx q[2];
rz(-0.59868559) q[2];
sx q[2];
rz(1.0318476) q[2];
rz(-1.6796238) q[3];
sx q[3];
rz(-0.69457355) q[3];
sx q[3];
rz(1.3612548) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66330355) q[0];
sx q[0];
rz(-2.7155868) q[0];
sx q[0];
rz(1.7167094) q[0];
rz(-2.5583963) q[1];
sx q[1];
rz(-1.2251264) q[1];
sx q[1];
rz(0.057417631) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9571788) q[0];
sx q[0];
rz(-2.3382362) q[0];
sx q[0];
rz(-3.1367541) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67441794) q[2];
sx q[2];
rz(-0.73235529) q[2];
sx q[2];
rz(-0.44841097) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4392885) q[1];
sx q[1];
rz(-1.3953475) q[1];
sx q[1];
rz(2.5006341) q[1];
rz(-pi) q[2];
rz(1.1693277) q[3];
sx q[3];
rz(-0.9801995) q[3];
sx q[3];
rz(-0.99311738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3164369) q[2];
sx q[2];
rz(-1.8437443) q[2];
sx q[2];
rz(1.8334897) q[2];
rz(1.3206652) q[3];
sx q[3];
rz(-2.2949009) q[3];
sx q[3];
rz(-2.2804885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2945781) q[0];
sx q[0];
rz(-0.68462831) q[0];
sx q[0];
rz(1.2680049) q[0];
rz(-1.1216724) q[1];
sx q[1];
rz(-1.2753762) q[1];
sx q[1];
rz(-2.4915288) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0065816) q[0];
sx q[0];
rz(-2.0913634) q[0];
sx q[0];
rz(2.6755946) q[0];
rz(-pi) q[1];
rz(-0.6261601) q[2];
sx q[2];
rz(-0.78842064) q[2];
sx q[2];
rz(1.5962708) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6396811) q[1];
sx q[1];
rz(-2.4458652) q[1];
sx q[1];
rz(0.63668107) q[1];
rz(-pi) q[2];
x q[2];
rz(0.995618) q[3];
sx q[3];
rz(-1.9891959) q[3];
sx q[3];
rz(-2.0783238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0012297) q[2];
sx q[2];
rz(-1.659212) q[2];
sx q[2];
rz(-2.7492827) q[2];
rz(0.20491925) q[3];
sx q[3];
rz(-2.6409918) q[3];
sx q[3];
rz(-2.4216381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1392764) q[0];
sx q[0];
rz(-1.5895695) q[0];
sx q[0];
rz(-1.4578777) q[0];
rz(2.814759) q[1];
sx q[1];
rz(-1.146233) q[1];
sx q[1];
rz(-2.0976417) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21699587) q[0];
sx q[0];
rz(-1.1490273) q[0];
sx q[0];
rz(0.20713465) q[0];
rz(-pi) q[1];
rz(-2.86356) q[2];
sx q[2];
rz(-2.6671609) q[2];
sx q[2];
rz(2.5086046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1684224) q[1];
sx q[1];
rz(-2.7379824) q[1];
sx q[1];
rz(1.7803867) q[1];
x q[2];
rz(-0.92071988) q[3];
sx q[3];
rz(-1.8333922) q[3];
sx q[3];
rz(0.85185862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8331208) q[2];
sx q[2];
rz(-1.0288419) q[2];
sx q[2];
rz(2.2157045) q[2];
rz(-1.067767) q[3];
sx q[3];
rz(-1.1238778) q[3];
sx q[3];
rz(-1.3656176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2866216) q[0];
sx q[0];
rz(-1.6678896) q[0];
sx q[0];
rz(2.5545252) q[0];
rz(-0.58285561) q[1];
sx q[1];
rz(-2.8522377) q[1];
sx q[1];
rz(-0.48097441) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60763393) q[0];
sx q[0];
rz(-1.9780428) q[0];
sx q[0];
rz(-1.3429705) q[0];
x q[1];
rz(-1.2103542) q[2];
sx q[2];
rz(-1.3925465) q[2];
sx q[2];
rz(1.0457863) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9900283) q[1];
sx q[1];
rz(-1.0326325) q[1];
sx q[1];
rz(1.0721704) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.770788) q[3];
sx q[3];
rz(-1.963435) q[3];
sx q[3];
rz(-2.5225366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9440072) q[2];
sx q[2];
rz(-2.6052167) q[2];
sx q[2];
rz(-1.3295004) q[2];
rz(2.3594989) q[3];
sx q[3];
rz(-0.32556459) q[3];
sx q[3];
rz(1.1480924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(2.6999577) q[0];
sx q[0];
rz(-1.4186207) q[0];
sx q[0];
rz(1.6765208) q[0];
rz(1.5326473) q[1];
sx q[1];
rz(-1.9736704) q[1];
sx q[1];
rz(-0.71221487) q[1];
rz(-0.92171348) q[2];
sx q[2];
rz(-0.94759181) q[2];
sx q[2];
rz(0.0061719051) q[2];
rz(0.28859517) q[3];
sx q[3];
rz(-2.6251818) q[3];
sx q[3];
rz(2.8040721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
