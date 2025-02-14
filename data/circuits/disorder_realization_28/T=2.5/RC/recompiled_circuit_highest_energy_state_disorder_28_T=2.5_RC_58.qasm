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
rz(-2.3856491) q[0];
sx q[0];
rz(-0.33060253) q[0];
sx q[0];
rz(0.27275738) q[0];
rz(2.7859712) q[1];
sx q[1];
rz(1.0300809) q[1];
sx q[1];
rz(7.8539943) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2897252) q[0];
sx q[0];
rz(-0.12061943) q[0];
sx q[0];
rz(-0.9319181) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9000969) q[2];
sx q[2];
rz(-2.0138171) q[2];
sx q[2];
rz(-1.4141413) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.570191) q[1];
sx q[1];
rz(-0.68059151) q[1];
sx q[1];
rz(-0.067318214) q[1];
rz(-pi) q[2];
rz(1.4727888) q[3];
sx q[3];
rz(-1.2036754) q[3];
sx q[3];
rz(-2.844732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1640799) q[2];
sx q[2];
rz(-1.432632) q[2];
sx q[2];
rz(2.254503) q[2];
rz(1.8767493) q[3];
sx q[3];
rz(-2.0829468) q[3];
sx q[3];
rz(-3.1177055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7510659) q[0];
sx q[0];
rz(-0.85705119) q[0];
sx q[0];
rz(-2.8476025) q[0];
rz(1.7901621) q[1];
sx q[1];
rz(-1.0963115) q[1];
sx q[1];
rz(1.9757804) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6751926) q[0];
sx q[0];
rz(-1.2443911) q[0];
sx q[0];
rz(-2.9221888) q[0];
rz(-0.027539081) q[2];
sx q[2];
rz(-1.6438494) q[2];
sx q[2];
rz(0.44833427) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6752515) q[1];
sx q[1];
rz(-2.0077809) q[1];
sx q[1];
rz(-2.5732424) q[1];
x q[2];
rz(-0.56880672) q[3];
sx q[3];
rz(-1.8173462) q[3];
sx q[3];
rz(2.5914471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4531648) q[2];
sx q[2];
rz(-2.5025949) q[2];
sx q[2];
rz(-1.9387908) q[2];
rz(1.5419675) q[3];
sx q[3];
rz(-1.8034233) q[3];
sx q[3];
rz(-2.8560824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11005345) q[0];
sx q[0];
rz(-3.0634614) q[0];
sx q[0];
rz(1.9670638) q[0];
rz(-0.9749167) q[1];
sx q[1];
rz(-0.87923032) q[1];
sx q[1];
rz(2.6851795) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74921319) q[0];
sx q[0];
rz(-1.8421041) q[0];
sx q[0];
rz(3.1086568) q[0];
x q[1];
rz(0.24791294) q[2];
sx q[2];
rz(-1.1968024) q[2];
sx q[2];
rz(-0.1318814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2043201) q[1];
sx q[1];
rz(-1.0144925) q[1];
sx q[1];
rz(2.4672303) q[1];
rz(-pi) q[2];
rz(-1.0393082) q[3];
sx q[3];
rz(-2.3520497) q[3];
sx q[3];
rz(1.9122151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5259033) q[2];
sx q[2];
rz(-0.38280767) q[2];
sx q[2];
rz(0.75695401) q[2];
rz(-3.1225539) q[3];
sx q[3];
rz(-1.2913387) q[3];
sx q[3];
rz(-0.835787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86209908) q[0];
sx q[0];
rz(-2.4308496) q[0];
sx q[0];
rz(0.94974649) q[0];
rz(0.21385916) q[1];
sx q[1];
rz(-2.2016134) q[1];
sx q[1];
rz(-0.59923879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1480162) q[0];
sx q[0];
rz(-2.2207119) q[0];
sx q[0];
rz(-0.17493576) q[0];
rz(-1.0548575) q[2];
sx q[2];
rz(-1.4241059) q[2];
sx q[2];
rz(-0.88146081) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.039700198) q[1];
sx q[1];
rz(-1.779161) q[1];
sx q[1];
rz(-2.1096257) q[1];
rz(-2.9327389) q[3];
sx q[3];
rz(-2.5503272) q[3];
sx q[3];
rz(0.2805191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5102101) q[2];
sx q[2];
rz(-1.6333132) q[2];
sx q[2];
rz(-0.82376662) q[2];
rz(2.0508164) q[3];
sx q[3];
rz(-2.3407276) q[3];
sx q[3];
rz(2.4763079) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7432231) q[0];
sx q[0];
rz(-0.38341612) q[0];
sx q[0];
rz(-0.98633352) q[0];
rz(-2.8979454) q[1];
sx q[1];
rz(-1.2657093) q[1];
sx q[1];
rz(2.9820014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3004813) q[0];
sx q[0];
rz(-2.2694025) q[0];
sx q[0];
rz(-2.7270856) q[0];
rz(-0.78686495) q[2];
sx q[2];
rz(-1.9965916) q[2];
sx q[2];
rz(2.6540861) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21503285) q[1];
sx q[1];
rz(-1.6590282) q[1];
sx q[1];
rz(2.0162321) q[1];
rz(-1.9165841) q[3];
sx q[3];
rz(-2.4107614) q[3];
sx q[3];
rz(-2.2737434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7470982) q[2];
sx q[2];
rz(-1.8793841) q[2];
sx q[2];
rz(-2.0060284) q[2];
rz(0.46843946) q[3];
sx q[3];
rz(-1.9197542) q[3];
sx q[3];
rz(-0.92456094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.121345) q[0];
sx q[0];
rz(-2.0609042) q[0];
sx q[0];
rz(-2.8998846) q[0];
rz(-1.3555591) q[1];
sx q[1];
rz(-2.2691998) q[1];
sx q[1];
rz(0.88368574) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.396326) q[0];
sx q[0];
rz(-0.81423212) q[0];
sx q[0];
rz(-7*pi/9) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4587347) q[2];
sx q[2];
rz(-0.39432967) q[2];
sx q[2];
rz(0.71428052) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8203796) q[1];
sx q[1];
rz(-1.1324778) q[1];
sx q[1];
rz(2.9258116) q[1];
rz(-0.00022907654) q[3];
sx q[3];
rz(-1.159824) q[3];
sx q[3];
rz(0.64490333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22150618) q[2];
sx q[2];
rz(-2.2111427) q[2];
sx q[2];
rz(2.7152756) q[2];
rz(0.91655556) q[3];
sx q[3];
rz(-1.5896268) q[3];
sx q[3];
rz(2.0385273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.351848) q[0];
sx q[0];
rz(-1.9639356) q[0];
sx q[0];
rz(-2.015693) q[0];
rz(0.060983505) q[1];
sx q[1];
rz(-2.1911502) q[1];
sx q[1];
rz(1.0152063) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7506152) q[0];
sx q[0];
rz(-0.96834194) q[0];
sx q[0];
rz(1.0573122) q[0];
x q[1];
rz(1.7776958) q[2];
sx q[2];
rz(-2.7963429) q[2];
sx q[2];
rz(2.783684) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.12830827) q[1];
sx q[1];
rz(-2.7375712) q[1];
sx q[1];
rz(-1.4133417) q[1];
rz(1.531008) q[3];
sx q[3];
rz(-1.6841751) q[3];
sx q[3];
rz(1.461352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.31275493) q[2];
sx q[2];
rz(-1.9373031) q[2];
sx q[2];
rz(-0.22605669) q[2];
rz(2.1894646) q[3];
sx q[3];
rz(-0.13855562) q[3];
sx q[3];
rz(-0.8086732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011618135) q[0];
sx q[0];
rz(-1.8267153) q[0];
sx q[0];
rz(-2.8011978) q[0];
rz(3.0420692) q[1];
sx q[1];
rz(-2.0390022) q[1];
sx q[1];
rz(1.5955101) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52852825) q[0];
sx q[0];
rz(-2.1648667) q[0];
sx q[0];
rz(1.1737002) q[0];
rz(-1.7695707) q[2];
sx q[2];
rz(-0.27789657) q[2];
sx q[2];
rz(3.1350435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73879646) q[1];
sx q[1];
rz(-1.5288884) q[1];
sx q[1];
rz(-3.0422076) q[1];
rz(-pi) q[2];
rz(-1.4813878) q[3];
sx q[3];
rz(-1.1173038) q[3];
sx q[3];
rz(-2.2553473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16326363) q[2];
sx q[2];
rz(-1.1922057) q[2];
sx q[2];
rz(-2.1400129) q[2];
rz(1.352365) q[3];
sx q[3];
rz(-2.6632301) q[3];
sx q[3];
rz(-1.9549687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773862) q[0];
sx q[0];
rz(-2.0851676) q[0];
sx q[0];
rz(0.00038432234) q[0];
rz(-2.7035825) q[1];
sx q[1];
rz(-1.6338467) q[1];
sx q[1];
rz(2.6843574) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9033636) q[0];
sx q[0];
rz(-1.6203124) q[0];
sx q[0];
rz(-1.7286506) q[0];
rz(-0.57827823) q[2];
sx q[2];
rz(-1.0338898) q[2];
sx q[2];
rz(-2.0202877) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8548566) q[1];
sx q[1];
rz(-1.7932307) q[1];
sx q[1];
rz(2.7871034) q[1];
x q[2];
rz(-2.6351686) q[3];
sx q[3];
rz(-1.862414) q[3];
sx q[3];
rz(-0.17403655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75927258) q[2];
sx q[2];
rz(-0.98733941) q[2];
sx q[2];
rz(-1.6299204) q[2];
rz(-3.0509389) q[3];
sx q[3];
rz(-2.1001308) q[3];
sx q[3];
rz(-2.4618728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35336211) q[0];
sx q[0];
rz(-1.7857977) q[0];
sx q[0];
rz(-2.8613388) q[0];
rz(-1.7578112) q[1];
sx q[1];
rz(-0.46418142) q[1];
sx q[1];
rz(1.2253449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68968407) q[0];
sx q[0];
rz(-1.6674918) q[0];
sx q[0];
rz(2.9908331) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7011048) q[2];
sx q[2];
rz(-2.2058528) q[2];
sx q[2];
rz(-1.6556513) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8333421) q[1];
sx q[1];
rz(-2.3099515) q[1];
sx q[1];
rz(1.9019446) q[1];
rz(-pi) q[2];
rz(-0.7627714) q[3];
sx q[3];
rz(-2.7342334) q[3];
sx q[3];
rz(2.6727904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19259351) q[2];
sx q[2];
rz(-1.5263564) q[2];
sx q[2];
rz(-2.4601649) q[2];
rz(-1.1000018) q[3];
sx q[3];
rz(-0.93183485) q[3];
sx q[3];
rz(-1.3070235) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6901907) q[0];
sx q[0];
rz(-2.207353) q[0];
sx q[0];
rz(1.9718715) q[0];
rz(2.7611217) q[1];
sx q[1];
rz(-0.66743757) q[1];
sx q[1];
rz(0.22421511) q[1];
rz(-1.2606332) q[2];
sx q[2];
rz(-2.416353) q[2];
sx q[2];
rz(-2.261434) q[2];
rz(-0.43222618) q[3];
sx q[3];
rz(-1.0226915) q[3];
sx q[3];
rz(2.6998479) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
