OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1563675) q[0];
sx q[0];
rz(-1.2824143) q[0];
sx q[0];
rz(3.0106944) q[0];
rz(-2.6137597) q[1];
sx q[1];
rz(-0.37017828) q[1];
sx q[1];
rz(0.17732492) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2180358) q[0];
sx q[0];
rz(-2.3157488) q[0];
sx q[0];
rz(2.7808583) q[0];
rz(-pi) q[1];
rz(-1.8715636) q[2];
sx q[2];
rz(-1.7890777) q[2];
sx q[2];
rz(1.0839562) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0874112) q[1];
sx q[1];
rz(-1.0980532) q[1];
sx q[1];
rz(-0.26473882) q[1];
rz(-pi) q[2];
rz(1.9779284) q[3];
sx q[3];
rz(-0.61591776) q[3];
sx q[3];
rz(-0.27995279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6191972) q[2];
sx q[2];
rz(-1.088524) q[2];
sx q[2];
rz(-1.3414475) q[2];
rz(-1.3522735) q[3];
sx q[3];
rz(-1.5315346) q[3];
sx q[3];
rz(1.8488098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.751048) q[0];
sx q[0];
rz(-1.9123257) q[0];
sx q[0];
rz(2.6265662) q[0];
rz(-0.36188778) q[1];
sx q[1];
rz(-1.3970951) q[1];
sx q[1];
rz(-2.1451758) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5866885) q[0];
sx q[0];
rz(-1.3350272) q[0];
sx q[0];
rz(-2.9352208) q[0];
rz(-0.17536945) q[2];
sx q[2];
rz(-1.9992454) q[2];
sx q[2];
rz(1.0340978) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0139549) q[1];
sx q[1];
rz(-1.6628237) q[1];
sx q[1];
rz(-2.0929232) q[1];
x q[2];
rz(0.44307905) q[3];
sx q[3];
rz(-1.0765809) q[3];
sx q[3];
rz(0.48393238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4126052) q[2];
sx q[2];
rz(-0.94634405) q[2];
sx q[2];
rz(1.4808572) q[2];
rz(-0.49731538) q[3];
sx q[3];
rz(-0.9484843) q[3];
sx q[3];
rz(-2.4174387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6052674) q[0];
sx q[0];
rz(-1.8999506) q[0];
sx q[0];
rz(-1.190825) q[0];
rz(0.39769998) q[1];
sx q[1];
rz(-1.5813446) q[1];
sx q[1];
rz(-2.2427028) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54972285) q[0];
sx q[0];
rz(-1.6373349) q[0];
sx q[0];
rz(-1.9107242) q[0];
rz(-pi) q[1];
rz(1.4853046) q[2];
sx q[2];
rz(-1.3807266) q[2];
sx q[2];
rz(-2.9051733) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0278416) q[1];
sx q[1];
rz(-1.6154624) q[1];
sx q[1];
rz(-1.0381446) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8679138) q[3];
sx q[3];
rz(-2.6032902) q[3];
sx q[3];
rz(-2.3915714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1206104) q[2];
sx q[2];
rz(-1.9526498) q[2];
sx q[2];
rz(0.19213842) q[2];
rz(2.4118679) q[3];
sx q[3];
rz(-0.14484043) q[3];
sx q[3];
rz(0.63989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9925053) q[0];
sx q[0];
rz(-2.3893116) q[0];
sx q[0];
rz(1.8335023) q[0];
rz(-0.84450841) q[1];
sx q[1];
rz(-1.4627855) q[1];
sx q[1];
rz(-0.62215296) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62646851) q[0];
sx q[0];
rz(-2.4822786) q[0];
sx q[0];
rz(2.4298682) q[0];
rz(-1.6255369) q[2];
sx q[2];
rz(-1.7559045) q[2];
sx q[2];
rz(-0.4601882) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.51650652) q[1];
sx q[1];
rz(-1.4072021) q[1];
sx q[1];
rz(-1.518599) q[1];
rz(0.34028168) q[3];
sx q[3];
rz(-2.8064686) q[3];
sx q[3];
rz(-1.7652896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1239803) q[2];
sx q[2];
rz(-1.7698741) q[2];
sx q[2];
rz(-2.2464216) q[2];
rz(0.34919843) q[3];
sx q[3];
rz(-0.57217351) q[3];
sx q[3];
rz(-2.9077742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2772086) q[0];
sx q[0];
rz(-1.207749) q[0];
sx q[0];
rz(2.5982017) q[0];
rz(3.0282989) q[1];
sx q[1];
rz(-1.7430867) q[1];
sx q[1];
rz(1.8494122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6627412) q[0];
sx q[0];
rz(-0.9060735) q[0];
sx q[0];
rz(0.49481884) q[0];
rz(-pi) q[1];
rz(0.92083709) q[2];
sx q[2];
rz(-1.2652768) q[2];
sx q[2];
rz(2.4561575) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0185011) q[1];
sx q[1];
rz(-0.95884174) q[1];
sx q[1];
rz(2.3101644) q[1];
x q[2];
rz(-0.75877996) q[3];
sx q[3];
rz(-0.89836392) q[3];
sx q[3];
rz(-2.1684949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.31305227) q[2];
sx q[2];
rz(-1.3581759) q[2];
sx q[2];
rz(-0.48946112) q[2];
rz(-0.66568565) q[3];
sx q[3];
rz(-0.61167115) q[3];
sx q[3];
rz(2.3556975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8996443) q[0];
sx q[0];
rz(-0.90270942) q[0];
sx q[0];
rz(3.0808501) q[0];
rz(-1.863106) q[1];
sx q[1];
rz(-2.2272031) q[1];
sx q[1];
rz(-2.1155105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2506977) q[0];
sx q[0];
rz(-2.1259667) q[0];
sx q[0];
rz(-2.5722136) q[0];
x q[1];
rz(2.0572971) q[2];
sx q[2];
rz(-2.5160273) q[2];
sx q[2];
rz(1.4826258) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85426165) q[1];
sx q[1];
rz(-1.8214487) q[1];
sx q[1];
rz(-1.4297559) q[1];
rz(-pi) q[2];
rz(-2.6137976) q[3];
sx q[3];
rz(-1.286881) q[3];
sx q[3];
rz(-0.75263587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8817899) q[2];
sx q[2];
rz(-1.5851574) q[2];
sx q[2];
rz(0.078710236) q[2];
rz(-1.6591266) q[3];
sx q[3];
rz(-2.8251298) q[3];
sx q[3];
rz(2.6993774) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.215613) q[0];
sx q[0];
rz(-0.38023606) q[0];
sx q[0];
rz(1.444814) q[0];
rz(-0.80698693) q[1];
sx q[1];
rz(-1.746256) q[1];
sx q[1];
rz(0.011946202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5552607) q[0];
sx q[0];
rz(-1.5440479) q[0];
sx q[0];
rz(2.2097951) q[0];
x q[1];
rz(1.4623227) q[2];
sx q[2];
rz(-2.0758934) q[2];
sx q[2];
rz(-2.6384356) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0199737) q[1];
sx q[1];
rz(-1.498292) q[1];
sx q[1];
rz(2.1825779) q[1];
x q[2];
rz(1.9254382) q[3];
sx q[3];
rz(-2.7193448) q[3];
sx q[3];
rz(2.1481882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5936467) q[2];
sx q[2];
rz(-0.30210364) q[2];
sx q[2];
rz(-2.3054403) q[2];
rz(3.0892843) q[3];
sx q[3];
rz(-1.6736504) q[3];
sx q[3];
rz(-2.2056495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9713822) q[0];
sx q[0];
rz(-0.8571856) q[0];
sx q[0];
rz(-0.49825391) q[0];
rz(1.0055297) q[1];
sx q[1];
rz(-1.6906831) q[1];
sx q[1];
rz(3.0301869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94994369) q[0];
sx q[0];
rz(-2.7102817) q[0];
sx q[0];
rz(1.3383207) q[0];
x q[1];
rz(-2.6314503) q[2];
sx q[2];
rz(-2.6594793) q[2];
sx q[2];
rz(0.040566074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5222539) q[1];
sx q[1];
rz(-2.2845387) q[1];
sx q[1];
rz(0.38222169) q[1];
rz(1.697568) q[3];
sx q[3];
rz(-1.2283192) q[3];
sx q[3];
rz(0.7966744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9446543) q[2];
sx q[2];
rz(-1.780218) q[2];
sx q[2];
rz(1.7353752) q[2];
rz(1.1574636) q[3];
sx q[3];
rz(-0.99298733) q[3];
sx q[3];
rz(1.5895313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.120753) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(-1.7027759) q[0];
rz(2.2976047) q[1];
sx q[1];
rz(-1.3071209) q[1];
sx q[1];
rz(-2.9356975) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5563894) q[0];
sx q[0];
rz(-2.2015338) q[0];
sx q[0];
rz(-1.5205193) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6240142) q[2];
sx q[2];
rz(-0.88376617) q[2];
sx q[2];
rz(-1.991697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8329551) q[1];
sx q[1];
rz(-1.7058733) q[1];
sx q[1];
rz(-1.6494688) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.467942) q[3];
sx q[3];
rz(-0.60551548) q[3];
sx q[3];
rz(-0.57693627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0042808) q[2];
sx q[2];
rz(-1.0575123) q[2];
sx q[2];
rz(-1.6005969) q[2];
rz(2.6565523) q[3];
sx q[3];
rz(-1.78616) q[3];
sx q[3];
rz(-3.0276827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4892905) q[0];
sx q[0];
rz(-0.052209608) q[0];
sx q[0];
rz(1.8452277) q[0];
rz(1.0507874) q[1];
sx q[1];
rz(-1.4605582) q[1];
sx q[1];
rz(2.5349862) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0832336) q[0];
sx q[0];
rz(-2.9748821) q[0];
sx q[0];
rz(1.7687083) q[0];
x q[1];
rz(1.6982618) q[2];
sx q[2];
rz(-2.1268788) q[2];
sx q[2];
rz(3.0024204) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1639599) q[1];
sx q[1];
rz(-0.31144413) q[1];
sx q[1];
rz(-2.9410579) q[1];
rz(-1.6876771) q[3];
sx q[3];
rz(-1.28834) q[3];
sx q[3];
rz(2.4931537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7217094) q[2];
sx q[2];
rz(-1.1315283) q[2];
sx q[2];
rz(-2.4427872) q[2];
rz(1.6803668) q[3];
sx q[3];
rz(-2.1278087) q[3];
sx q[3];
rz(-0.47718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2496495) q[0];
sx q[0];
rz(-1.5105381) q[0];
sx q[0];
rz(1.426209) q[0];
rz(2.7632948) q[1];
sx q[1];
rz(-0.62722798) q[1];
sx q[1];
rz(2.41082) q[1];
rz(-1.0155914) q[2];
sx q[2];
rz(-2.5980661) q[2];
sx q[2];
rz(-0.51474434) q[2];
rz(0.049765718) q[3];
sx q[3];
rz(-0.55703041) q[3];
sx q[3];
rz(1.218474) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
