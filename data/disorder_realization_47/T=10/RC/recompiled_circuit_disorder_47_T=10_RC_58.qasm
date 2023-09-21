OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(2.3072825) q[0];
sx q[0];
rz(6.351525) q[0];
rz(2.4812658) q[1];
sx q[1];
rz(3.9897444) q[1];
sx q[1];
rz(6.2453649) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4662287) q[0];
sx q[0];
rz(-0.20875202) q[0];
sx q[0];
rz(-1.9573613) q[0];
x q[1];
rz(0.1184138) q[2];
sx q[2];
rz(-2.2288433) q[2];
sx q[2];
rz(-2.7409035) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1483128) q[1];
sx q[1];
rz(-0.82511745) q[1];
sx q[1];
rz(-1.5013904) q[1];
rz(1.2245523) q[3];
sx q[3];
rz(-1.7131221) q[3];
sx q[3];
rz(3.0179253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7906856) q[2];
sx q[2];
rz(-0.51351341) q[2];
sx q[2];
rz(-1.2834056) q[2];
rz(-0.12617271) q[3];
sx q[3];
rz(-1.7254555) q[3];
sx q[3];
rz(-3.0509907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.693817) q[0];
sx q[0];
rz(-2.2651146) q[0];
sx q[0];
rz(0.59666657) q[0];
rz(1.5555351) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(-1.7780875) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7994021) q[0];
sx q[0];
rz(-2.0997093) q[0];
sx q[0];
rz(1.0728114) q[0];
rz(-1.5200137) q[2];
sx q[2];
rz(-1.0154468) q[2];
sx q[2];
rz(-2.3107049) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9047782) q[1];
sx q[1];
rz(-1.8493376) q[1];
sx q[1];
rz(-0.65794557) q[1];
rz(-pi) q[2];
rz(-2.2927106) q[3];
sx q[3];
rz(-2.1697681) q[3];
sx q[3];
rz(2.2357383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3229225) q[2];
sx q[2];
rz(-1.0007891) q[2];
sx q[2];
rz(-2.4070516) q[2];
rz(2.5143886) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4194141) q[0];
sx q[0];
rz(-2.3086771) q[0];
sx q[0];
rz(-2.869379) q[0];
rz(0.84725562) q[1];
sx q[1];
rz(-1.3191185) q[1];
sx q[1];
rz(-0.31262696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0659098) q[0];
sx q[0];
rz(-1.5172177) q[0];
sx q[0];
rz(-2.9220394) q[0];
rz(-2.0414511) q[2];
sx q[2];
rz(-1.0725642) q[2];
sx q[2];
rz(0.91980308) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0215175) q[1];
sx q[1];
rz(-2.8719963) q[1];
sx q[1];
rz(1.8345548) q[1];
rz(-pi) q[2];
rz(-2.3787141) q[3];
sx q[3];
rz(-1.2925914) q[3];
sx q[3];
rz(1.1080081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.12038885) q[2];
sx q[2];
rz(-0.52336064) q[2];
sx q[2];
rz(-2.3266501) q[2];
rz(0.96873823) q[3];
sx q[3];
rz(-2.0961943) q[3];
sx q[3];
rz(2.708784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.75354904) q[0];
sx q[0];
rz(-0.20064813) q[0];
sx q[0];
rz(2.5182305) q[0];
rz(-2.3240044) q[1];
sx q[1];
rz(-2.8249884) q[1];
sx q[1];
rz(-3.0923016) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.374644) q[0];
sx q[0];
rz(-0.76061237) q[0];
sx q[0];
rz(-1.5006554) q[0];
x q[1];
rz(-1.9360147) q[2];
sx q[2];
rz(-2.3962002) q[2];
sx q[2];
rz(0.12860194) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4781487) q[1];
sx q[1];
rz(-1.899292) q[1];
sx q[1];
rz(-0.25681396) q[1];
x q[2];
rz(-0.44932511) q[3];
sx q[3];
rz(-1.0885363) q[3];
sx q[3];
rz(-2.0039909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7307044) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(0.98497406) q[2];
rz(-0.90304053) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(2.8682017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7966998) q[0];
sx q[0];
rz(-1.6051689) q[0];
sx q[0];
rz(3.0539736) q[0];
rz(2.9852729) q[1];
sx q[1];
rz(-0.55115288) q[1];
sx q[1];
rz(0.87096754) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5694215) q[0];
sx q[0];
rz(-1.2835802) q[0];
sx q[0];
rz(-1.5409018) q[0];
x q[1];
rz(2.9882177) q[2];
sx q[2];
rz(-1.9605637) q[2];
sx q[2];
rz(2.6467121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.85775447) q[1];
sx q[1];
rz(-2.4205967) q[1];
sx q[1];
rz(2.7096219) q[1];
x q[2];
rz(-0.80645873) q[3];
sx q[3];
rz(-1.4440923) q[3];
sx q[3];
rz(0.63092953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0481723) q[2];
sx q[2];
rz(-2.293736) q[2];
sx q[2];
rz(2.7887153) q[2];
rz(-0.86587632) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(-2.849259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46309328) q[0];
sx q[0];
rz(-2.0149639) q[0];
sx q[0];
rz(-0.10678664) q[0];
rz(1.9550025) q[1];
sx q[1];
rz(-1.0266961) q[1];
sx q[1];
rz(2.1170763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5750676) q[0];
sx q[0];
rz(-2.0216414) q[0];
sx q[0];
rz(1.8917811) q[0];
rz(-0.60501955) q[2];
sx q[2];
rz(-1.9540678) q[2];
sx q[2];
rz(0.59125102) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34895996) q[1];
sx q[1];
rz(-2.0434725) q[1];
sx q[1];
rz(2.2036168) q[1];
x q[2];
rz(-0.6242574) q[3];
sx q[3];
rz(-0.58371021) q[3];
sx q[3];
rz(-2.7492085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7496877) q[2];
sx q[2];
rz(-1.9275894) q[2];
sx q[2];
rz(0.8708896) q[2];
rz(1.332256) q[3];
sx q[3];
rz(-0.94227666) q[3];
sx q[3];
rz(1.7308621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1720599) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(-2.5906738) q[0];
rz(0.46539601) q[1];
sx q[1];
rz(-0.67133633) q[1];
sx q[1];
rz(-0.23682061) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.119334) q[0];
sx q[0];
rz(-0.10611457) q[0];
sx q[0];
rz(-1.2073713) q[0];
rz(-pi) q[1];
rz(0.74117629) q[2];
sx q[2];
rz(-2.9393662) q[2];
sx q[2];
rz(-1.2349595) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.782978) q[1];
sx q[1];
rz(-1.8404984) q[1];
sx q[1];
rz(-0.077809212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2348433) q[3];
sx q[3];
rz(-0.7630322) q[3];
sx q[3];
rz(-1.2515765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6992496) q[2];
sx q[2];
rz(-0.68168679) q[2];
sx q[2];
rz(2.701475) q[2];
rz(-2.0464499) q[3];
sx q[3];
rz(-1.4705855) q[3];
sx q[3];
rz(-1.346689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31297627) q[0];
sx q[0];
rz(-1.3422817) q[0];
sx q[0];
rz(3.112088) q[0];
rz(0.94738952) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(3.0227919) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9751562) q[0];
sx q[0];
rz(-1.6140037) q[0];
sx q[0];
rz(2.8309566) q[0];
rz(-2.115909) q[2];
sx q[2];
rz(-1.9695749) q[2];
sx q[2];
rz(-1.1023956) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.006281) q[1];
sx q[1];
rz(-2.0326008) q[1];
sx q[1];
rz(-1.9001574) q[1];
x q[2];
rz(0.16468594) q[3];
sx q[3];
rz(-2.0593004) q[3];
sx q[3];
rz(0.96795852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3999195) q[2];
sx q[2];
rz(-0.46367773) q[2];
sx q[2];
rz(1.5734394) q[2];
rz(1.167477) q[3];
sx q[3];
rz(-1.4195331) q[3];
sx q[3];
rz(-1.8673816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90299273) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(2.9883244) q[0];
rz(1.0614456) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(-1.1538039) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35837367) q[0];
sx q[0];
rz(-1.8405387) q[0];
sx q[0];
rz(-0.37958919) q[0];
x q[1];
rz(0.93038721) q[2];
sx q[2];
rz(-2.3500035) q[2];
sx q[2];
rz(1.9554348) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3378539) q[1];
sx q[1];
rz(-1.6263464) q[1];
sx q[1];
rz(-1.3615723) q[1];
rz(2.2581402) q[3];
sx q[3];
rz(-1.6646107) q[3];
sx q[3];
rz(2.1283619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29685059) q[2];
sx q[2];
rz(-1.9694318) q[2];
sx q[2];
rz(-1.5273757) q[2];
rz(2.1447003) q[3];
sx q[3];
rz(-1.6485873) q[3];
sx q[3];
rz(0.87866384) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8940354) q[0];
sx q[0];
rz(-1.0972801) q[0];
sx q[0];
rz(-0.19009185) q[0];
rz(2.4709573) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(1.488283) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7877486) q[0];
sx q[0];
rz(-1.6962546) q[0];
sx q[0];
rz(1.6018484) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3004296) q[2];
sx q[2];
rz(-1.1425945) q[2];
sx q[2];
rz(2.2141475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.59939811) q[1];
sx q[1];
rz(-2.290526) q[1];
sx q[1];
rz(-2.4464314) q[1];
rz(-pi) q[2];
rz(2.2641803) q[3];
sx q[3];
rz(-2.0392232) q[3];
sx q[3];
rz(-1.3437769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0037447475) q[2];
sx q[2];
rz(-2.2403084) q[2];
sx q[2];
rz(2.2424662) q[2];
rz(2.6265465) q[3];
sx q[3];
rz(-2.6656272) q[3];
sx q[3];
rz(-0.17764828) q[3];
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
rz(-1.0777733) q[0];
sx q[0];
rz(-1.780092) q[0];
sx q[0];
rz(-0.5583981) q[0];
rz(-2.8521815) q[1];
sx q[1];
rz(-0.90129539) q[1];
sx q[1];
rz(1.7064066) q[1];
rz(-2.4317447) q[2];
sx q[2];
rz(-1.3550497) q[2];
sx q[2];
rz(1.4646127) q[2];
rz(-1.297847) q[3];
sx q[3];
rz(-1.6932586) q[3];
sx q[3];
rz(0.3441588) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];