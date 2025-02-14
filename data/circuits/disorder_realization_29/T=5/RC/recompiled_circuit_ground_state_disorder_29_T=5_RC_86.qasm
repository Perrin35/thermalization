OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.090341181) q[0];
sx q[0];
rz(2.0988965) q[0];
sx q[0];
rz(12.287696) q[0];
rz(-1.9947808) q[1];
sx q[1];
rz(-2.6169701) q[1];
sx q[1];
rz(-1.8705179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4500101) q[0];
sx q[0];
rz(-1.5890536) q[0];
sx q[0];
rz(0.26055793) q[0];
rz(-pi) q[1];
rz(-2.991011) q[2];
sx q[2];
rz(-1.7007593) q[2];
sx q[2];
rz(-0.58565631) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6587352) q[1];
sx q[1];
rz(-2.4244011) q[1];
sx q[1];
rz(1.4499922) q[1];
rz(-pi) q[2];
rz(2.5335629) q[3];
sx q[3];
rz(-1.078106) q[3];
sx q[3];
rz(1.0580491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0311387) q[2];
sx q[2];
rz(-1.9240856) q[2];
sx q[2];
rz(-0.26220599) q[2];
rz(-1.0555222) q[3];
sx q[3];
rz(-1.8942098) q[3];
sx q[3];
rz(-0.087337703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6156886) q[0];
sx q[0];
rz(-0.50351244) q[0];
sx q[0];
rz(-0.18158922) q[0];
rz(1.7407821) q[1];
sx q[1];
rz(-1.1927651) q[1];
sx q[1];
rz(-2.3006732) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940482) q[0];
sx q[0];
rz(-1.2508014) q[0];
sx q[0];
rz(-1.589369) q[0];
rz(-2.1757751) q[2];
sx q[2];
rz(-2.4060095) q[2];
sx q[2];
rz(0.090026131) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.629417) q[1];
sx q[1];
rz(-1.4025153) q[1];
sx q[1];
rz(3.1049033) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1827994) q[3];
sx q[3];
rz(-1.8763767) q[3];
sx q[3];
rz(-0.45388729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38222868) q[2];
sx q[2];
rz(-0.21491924) q[2];
sx q[2];
rz(-2.0767029) q[2];
rz(3.0801638) q[3];
sx q[3];
rz(-1.8062785) q[3];
sx q[3];
rz(1.6399062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29838022) q[0];
sx q[0];
rz(-1.9623663) q[0];
sx q[0];
rz(-0.19905736) q[0];
rz(0.82636034) q[1];
sx q[1];
rz(-0.91824707) q[1];
sx q[1];
rz(0.77702776) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26829241) q[0];
sx q[0];
rz(-1.4140633) q[0];
sx q[0];
rz(0.066989338) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42962809) q[2];
sx q[2];
rz(-1.0681392) q[2];
sx q[2];
rz(2.4328897) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72973452) q[1];
sx q[1];
rz(-2.5246906) q[1];
sx q[1];
rz(-2.0286948) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8641508) q[3];
sx q[3];
rz(-2.146477) q[3];
sx q[3];
rz(0.68439181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1021154) q[2];
sx q[2];
rz(-1.9934883) q[2];
sx q[2];
rz(2.6285505) q[2];
rz(-0.78220621) q[3];
sx q[3];
rz(-1.2057883) q[3];
sx q[3];
rz(-1.7647083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.9272598) q[0];
sx q[0];
rz(-1.5981277) q[0];
sx q[0];
rz(-1.4904892) q[0];
rz(-0.53228846) q[1];
sx q[1];
rz(-0.63096255) q[1];
sx q[1];
rz(1.6040364) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69025264) q[0];
sx q[0];
rz(-1.7336539) q[0];
sx q[0];
rz(-1.4630408) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093393383) q[2];
sx q[2];
rz(-1.0265145) q[2];
sx q[2];
rz(-1.5285847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73276943) q[1];
sx q[1];
rz(-0.63946001) q[1];
sx q[1];
rz(-1.4519522) q[1];
x q[2];
rz(0.52738366) q[3];
sx q[3];
rz(-1.77447) q[3];
sx q[3];
rz(-2.7319997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6060467) q[2];
sx q[2];
rz(-2.594785) q[2];
sx q[2];
rz(0.01288506) q[2];
rz(1.9384725) q[3];
sx q[3];
rz(-1.4667526) q[3];
sx q[3];
rz(-2.119868) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0032229) q[0];
sx q[0];
rz(-2.5165181) q[0];
sx q[0];
rz(-1.3997929) q[0];
rz(-2.7895582) q[1];
sx q[1];
rz(-2.0875918) q[1];
sx q[1];
rz(-2.3037516) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8796757) q[0];
sx q[0];
rz(-2.3220372) q[0];
sx q[0];
rz(-2.9330334) q[0];
rz(2.1878686) q[2];
sx q[2];
rz(-1.0487818) q[2];
sx q[2];
rz(-1.2921289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4427691) q[1];
sx q[1];
rz(-1.6976446) q[1];
sx q[1];
rz(-2.5731099) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9772163) q[3];
sx q[3];
rz(-1.419074) q[3];
sx q[3];
rz(-2.3325932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0248523) q[2];
sx q[2];
rz(-0.96472538) q[2];
sx q[2];
rz(-1.1836729) q[2];
rz(-0.27070326) q[3];
sx q[3];
rz(-0.94047061) q[3];
sx q[3];
rz(-1.6089599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026641332) q[0];
sx q[0];
rz(-2.5178435) q[0];
sx q[0];
rz(2.5624516) q[0];
rz(-1.9193513) q[1];
sx q[1];
rz(-1.3074343) q[1];
sx q[1];
rz(2.0319669) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1117119) q[0];
sx q[0];
rz(-1.7222026) q[0];
sx q[0];
rz(-0.16686186) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8366361) q[2];
sx q[2];
rz(-1.6246319) q[2];
sx q[2];
rz(0.75537813) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5803924) q[1];
sx q[1];
rz(-1.8186186) q[1];
sx q[1];
rz(2.7996254) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3253981) q[3];
sx q[3];
rz(-2.4494054) q[3];
sx q[3];
rz(-0.029267197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1960725) q[2];
sx q[2];
rz(-1.2976982) q[2];
sx q[2];
rz(-0.35745364) q[2];
rz(-1.9516021) q[3];
sx q[3];
rz(-1.9414732) q[3];
sx q[3];
rz(-1.801871) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59915197) q[0];
sx q[0];
rz(-1.0170794) q[0];
sx q[0];
rz(-2.3495112) q[0];
rz(2.4374841) q[1];
sx q[1];
rz(-0.45583615) q[1];
sx q[1];
rz(-1.5584996) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38403758) q[0];
sx q[0];
rz(-2.677394) q[0];
sx q[0];
rz(3.1350101) q[0];
rz(-pi) q[1];
rz(-0.94740811) q[2];
sx q[2];
rz(-2.0724503) q[2];
sx q[2];
rz(-2.9497408) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.577282) q[1];
sx q[1];
rz(-2.7140712) q[1];
sx q[1];
rz(-0.083210817) q[1];
rz(-pi) q[2];
rz(-1.3450479) q[3];
sx q[3];
rz(-2.3093947) q[3];
sx q[3];
rz(-2.3647472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11129561) q[2];
sx q[2];
rz(-2.3534677) q[2];
sx q[2];
rz(0.0041858717) q[2];
rz(0.26675102) q[3];
sx q[3];
rz(-1.5240074) q[3];
sx q[3];
rz(-0.72223103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5608212) q[0];
sx q[0];
rz(-2.4120791) q[0];
sx q[0];
rz(0.15604493) q[0];
rz(1.5566298) q[1];
sx q[1];
rz(-2.2792351) q[1];
sx q[1];
rz(2.5468266) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9879919) q[0];
sx q[0];
rz(-0.56454851) q[0];
sx q[0];
rz(-0.49086824) q[0];
rz(2.4492599) q[2];
sx q[2];
rz(-1.7966447) q[2];
sx q[2];
rz(-0.69162265) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.22560355) q[1];
sx q[1];
rz(-2.2609613) q[1];
sx q[1];
rz(0.76721104) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72130008) q[3];
sx q[3];
rz(-2.2709322) q[3];
sx q[3];
rz(-0.68381468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5742089) q[2];
sx q[2];
rz(-1.5625861) q[2];
sx q[2];
rz(2.1136005) q[2];
rz(-1.3422286) q[3];
sx q[3];
rz(-2.4864311) q[3];
sx q[3];
rz(1.8067693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26355711) q[0];
sx q[0];
rz(-1.4565383) q[0];
sx q[0];
rz(-0.79826075) q[0];
rz(2.3410666) q[1];
sx q[1];
rz(-2.6526178) q[1];
sx q[1];
rz(-0.54623234) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52926621) q[0];
sx q[0];
rz(-0.61175013) q[0];
sx q[0];
rz(1.7609558) q[0];
rz(1.1425708) q[2];
sx q[2];
rz(-1.0248803) q[2];
sx q[2];
rz(0.8949309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7475976) q[1];
sx q[1];
rz(-1.0859274) q[1];
sx q[1];
rz(-2.5970244) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22643666) q[3];
sx q[3];
rz(-1.4996935) q[3];
sx q[3];
rz(-1.16217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3937248) q[2];
sx q[2];
rz(-2.0515029) q[2];
sx q[2];
rz(-3.0785576) q[2];
rz(2.4437599) q[3];
sx q[3];
rz(-2.883426) q[3];
sx q[3];
rz(-0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49552712) q[0];
sx q[0];
rz(-0.87764469) q[0];
sx q[0];
rz(-0.416042) q[0];
rz(-1.3620194) q[1];
sx q[1];
rz(-2.0174618) q[1];
sx q[1];
rz(-1.8488041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1071615) q[0];
sx q[0];
rz(-0.35625544) q[0];
sx q[0];
rz(0.70888575) q[0];
rz(-pi) q[1];
rz(0.52915574) q[2];
sx q[2];
rz(-2.6247027) q[2];
sx q[2];
rz(-0.86597499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19504612) q[1];
sx q[1];
rz(-2.1113696) q[1];
sx q[1];
rz(-0.54624301) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1077652) q[3];
sx q[3];
rz(-2.6601048) q[3];
sx q[3];
rz(-0.010537174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4585939) q[2];
sx q[2];
rz(-1.3364044) q[2];
sx q[2];
rz(2.2340753) q[2];
rz(1.2609743) q[3];
sx q[3];
rz(-2.5670299) q[3];
sx q[3];
rz(-1.6921836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4954729) q[0];
sx q[0];
rz(-1.4541805) q[0];
sx q[0];
rz(-2.4021586) q[0];
rz(-0.44838913) q[1];
sx q[1];
rz(-1.3403475) q[1];
sx q[1];
rz(1.1833804) q[1];
rz(-0.72095966) q[2];
sx q[2];
rz(-2.6868281) q[2];
sx q[2];
rz(1.6285298) q[2];
rz(-0.84239324) q[3];
sx q[3];
rz(-1.471023) q[3];
sx q[3];
rz(-1.4820549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
