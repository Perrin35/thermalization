OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9632602) q[0];
sx q[0];
rz(-1.652521) q[0];
sx q[0];
rz(0.89515495) q[0];
rz(-0.31495467) q[1];
sx q[1];
rz(-2.0575674) q[1];
sx q[1];
rz(-1.6853583) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8946978) q[0];
sx q[0];
rz(-0.1323192) q[0];
sx q[0];
rz(1.4207065) q[0];
rz(-pi) q[1];
rz(-0.37402447) q[2];
sx q[2];
rz(-1.0569388) q[2];
sx q[2];
rz(-2.6467269) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.734182) q[1];
sx q[1];
rz(-1.721518) q[1];
sx q[1];
rz(0.68024866) q[1];
rz(-2.9772894) q[3];
sx q[3];
rz(-2.8176753) q[3];
sx q[3];
rz(1.8620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75498092) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(-2.6888729) q[2];
rz(2.9833941) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92900705) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(-1.989495) q[0];
rz(1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(-0.47098413) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4244714) q[0];
sx q[0];
rz(-1.0891799) q[0];
sx q[0];
rz(-0.061901285) q[0];
x q[1];
rz(0.23439622) q[2];
sx q[2];
rz(-0.98653754) q[2];
sx q[2];
rz(2.5258979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4091332) q[1];
sx q[1];
rz(-2.0583378) q[1];
sx q[1];
rz(3.0797466) q[1];
x q[2];
rz(-1.6509634) q[3];
sx q[3];
rz(-0.87682322) q[3];
sx q[3];
rz(0.44192867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.058078893) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(-0.2557959) q[2];
rz(-1.6563709) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(-2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(2.8783022) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(-0.27134744) q[0];
rz(2.4052606) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(2.7022865) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8894316) q[0];
sx q[0];
rz(-1.6003506) q[0];
sx q[0];
rz(-1.6512524) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13621026) q[2];
sx q[2];
rz(-0.36362193) q[2];
sx q[2];
rz(-1.1656392) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1767133) q[1];
sx q[1];
rz(-0.63767725) q[1];
sx q[1];
rz(1.8646851) q[1];
x q[2];
rz(-0.93033354) q[3];
sx q[3];
rz(-2.4822682) q[3];
sx q[3];
rz(-1.7742771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1674041) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(-0.75508562) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(-2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(-1.8970998) q[0];
rz(2.3311133) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(-0.8746075) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0229605) q[0];
sx q[0];
rz(-1.7264139) q[0];
sx q[0];
rz(-2.1778584) q[0];
rz(-pi) q[1];
rz(1.7165519) q[2];
sx q[2];
rz(-1.1894023) q[2];
sx q[2];
rz(-0.33304735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3687467) q[1];
sx q[1];
rz(-0.44279848) q[1];
sx q[1];
rz(2.5889791) q[1];
rz(1.6455669) q[3];
sx q[3];
rz(-1.6755591) q[3];
sx q[3];
rz(-1.2391702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0041634) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(-1.4271663) q[2];
rz(-0.066120474) q[3];
sx q[3];
rz(-2.7757006) q[3];
sx q[3];
rz(1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(0.57058913) q[0];
rz(-2.5866306) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(-0.74329174) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76737228) q[0];
sx q[0];
rz(-1.351007) q[0];
sx q[0];
rz(2.2348997) q[0];
rz(-pi) q[1];
rz(0.18851738) q[2];
sx q[2];
rz(-0.50893116) q[2];
sx q[2];
rz(-2.237042) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9651282) q[1];
sx q[1];
rz(-1.3511409) q[1];
sx q[1];
rz(-1.7996644) q[1];
x q[2];
rz(-1.8573496) q[3];
sx q[3];
rz(-0.78714579) q[3];
sx q[3];
rz(-0.033165008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(0.26838475) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(-2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140117) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-2.6348689) q[0];
rz(2.8880033) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(-1.0553029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6049833) q[0];
sx q[0];
rz(-2.4924926) q[0];
sx q[0];
rz(1.2654632) q[0];
rz(-pi) q[1];
rz(-0.92743404) q[2];
sx q[2];
rz(-2.3775527) q[2];
sx q[2];
rz(1.6709136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9678952) q[1];
sx q[1];
rz(-2.3579862) q[1];
sx q[1];
rz(-2.150279) q[1];
x q[2];
rz(0.25552337) q[3];
sx q[3];
rz(-0.75948411) q[3];
sx q[3];
rz(-1.7373191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6598597) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(2.4957538) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6362474) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(-2.3690467) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(-2.5678182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0550802) q[0];
sx q[0];
rz(-1.3801314) q[0];
sx q[0];
rz(3.029682) q[0];
rz(-pi) q[1];
rz(3.119486) q[2];
sx q[2];
rz(-1.7709641) q[2];
sx q[2];
rz(-2.9526763) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7635203) q[1];
sx q[1];
rz(-1.0135279) q[1];
sx q[1];
rz(0.65655638) q[1];
rz(2.1366828) q[3];
sx q[3];
rz(-2.8661869) q[3];
sx q[3];
rz(-2.3435081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1022169) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(-2.3708564) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.27281) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(-1.1707206) q[0];
rz(2.6314578) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(1.8458813) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9840178) q[0];
sx q[0];
rz(-0.57849738) q[0];
sx q[0];
rz(2.5111141) q[0];
rz(-pi) q[1];
rz(1.2175351) q[2];
sx q[2];
rz(-1.8262987) q[2];
sx q[2];
rz(-3.0614292) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1295373) q[1];
sx q[1];
rz(-0.66857282) q[1];
sx q[1];
rz(1.2052016) q[1];
x q[2];
rz(1.6298953) q[3];
sx q[3];
rz(-2.700138) q[3];
sx q[3];
rz(0.7837226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7065113) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(-0.56813017) q[2];
rz(2.1614697) q[3];
sx q[3];
rz(-1.0390493) q[3];
sx q[3];
rz(-0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66529626) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(2.4639159) q[0];
rz(-0.19605818) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(-2.303404) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8294551) q[0];
sx q[0];
rz(-1.9337618) q[0];
sx q[0];
rz(1.3296933) q[0];
rz(-pi) q[1];
rz(-1.499275) q[2];
sx q[2];
rz(-0.75746775) q[2];
sx q[2];
rz(1.0345936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6088241) q[1];
sx q[1];
rz(-0.84034398) q[1];
sx q[1];
rz(2.8735012) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6543051) q[3];
sx q[3];
rz(-1.3132846) q[3];
sx q[3];
rz(-0.53572342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80205408) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(-2.1789815) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4203913) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(2.904073) q[0];
rz(-1.0182084) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(0.231803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0203665) q[0];
sx q[0];
rz(-1.0640642) q[0];
sx q[0];
rz(-0.10911848) q[0];
x q[1];
rz(2.3583057) q[2];
sx q[2];
rz(-2.8056393) q[2];
sx q[2];
rz(2.96539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.17813645) q[1];
sx q[1];
rz(-1.3061211) q[1];
sx q[1];
rz(-1.4241649) q[1];
x q[2];
rz(-2.9276804) q[3];
sx q[3];
rz(-2.247346) q[3];
sx q[3];
rz(-3.1090581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(0.38816372) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647472) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(-0.96881962) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(1.8112524) q[2];
sx q[2];
rz(-1.8869055) q[2];
sx q[2];
rz(-1.4958924) q[2];
rz(-0.10235056) q[3];
sx q[3];
rz(-0.1889189) q[3];
sx q[3];
rz(1.2253996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
