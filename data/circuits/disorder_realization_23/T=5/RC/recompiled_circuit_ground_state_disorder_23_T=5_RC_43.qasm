OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57103676) q[0];
sx q[0];
rz(3.8138226) q[0];
sx q[0];
rz(11.091118) q[0];
rz(0.9530468) q[1];
sx q[1];
rz(3.3068125) q[1];
sx q[1];
rz(10.693476) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3606963) q[0];
sx q[0];
rz(-1.3421881) q[0];
sx q[0];
rz(0.40970384) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95999865) q[2];
sx q[2];
rz(-0.75163254) q[2];
sx q[2];
rz(-3.1267191) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7172473) q[1];
sx q[1];
rz(-1.0043036) q[1];
sx q[1];
rz(-0.15498181) q[1];
x q[2];
rz(-1.8016863) q[3];
sx q[3];
rz(-0.82726631) q[3];
sx q[3];
rz(-0.32318599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1573726) q[2];
sx q[2];
rz(-2.7998388) q[2];
sx q[2];
rz(-2.3483707) q[2];
rz(-0.67304099) q[3];
sx q[3];
rz(-1.0854191) q[3];
sx q[3];
rz(-1.8719155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15758812) q[0];
sx q[0];
rz(-2.3389811) q[0];
sx q[0];
rz(-0.60930914) q[0];
rz(2.9318103) q[1];
sx q[1];
rz(-0.79133004) q[1];
sx q[1];
rz(-0.9598859) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083733746) q[0];
sx q[0];
rz(-1.4918033) q[0];
sx q[0];
rz(1.5261488) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0591956) q[2];
sx q[2];
rz(-1.5090183) q[2];
sx q[2];
rz(-0.90324963) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4757735) q[1];
sx q[1];
rz(-0.83516652) q[1];
sx q[1];
rz(-0.084666157) q[1];
rz(-pi) q[2];
rz(-2.2678492) q[3];
sx q[3];
rz(-0.62752073) q[3];
sx q[3];
rz(2.2604159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7858872) q[2];
sx q[2];
rz(-2.3023119) q[2];
sx q[2];
rz(-0.4915702) q[2];
rz(0.0056886557) q[3];
sx q[3];
rz(-2.0523043) q[3];
sx q[3];
rz(2.6277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.046535) q[0];
sx q[0];
rz(-2.6061366) q[0];
sx q[0];
rz(2.3989578) q[0];
rz(-0.62478089) q[1];
sx q[1];
rz(-2.2014328) q[1];
sx q[1];
rz(2.0754441) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4542592) q[0];
sx q[0];
rz(-1.3479184) q[0];
sx q[0];
rz(-2.9965626) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9511388) q[2];
sx q[2];
rz(-1.2553658) q[2];
sx q[2];
rz(-0.73376943) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6934476) q[1];
sx q[1];
rz(-1.9733917) q[1];
sx q[1];
rz(-0.13685302) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0622018) q[3];
sx q[3];
rz(-1.3771476) q[3];
sx q[3];
rz(-1.0102864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2523969) q[2];
sx q[2];
rz(-2.9434581) q[2];
sx q[2];
rz(-2.5890217) q[2];
rz(-2.6342454) q[3];
sx q[3];
rz(-1.0891958) q[3];
sx q[3];
rz(2.5721917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6602537) q[0];
sx q[0];
rz(-2.5642671) q[0];
sx q[0];
rz(-1.2465771) q[0];
rz(-1.4056816) q[1];
sx q[1];
rz(-0.77189267) q[1];
sx q[1];
rz(-0.71800047) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3858741) q[0];
sx q[0];
rz(-2.0760975) q[0];
sx q[0];
rz(0.85407172) q[0];
rz(2.1338855) q[2];
sx q[2];
rz(-1.3674842) q[2];
sx q[2];
rz(2.2846534) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2110465) q[1];
sx q[1];
rz(-0.77096838) q[1];
sx q[1];
rz(0.20937415) q[1];
x q[2];
rz(-0.15851373) q[3];
sx q[3];
rz(-1.6974546) q[3];
sx q[3];
rz(1.4617006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8743073) q[2];
sx q[2];
rz(-2.3882046) q[2];
sx q[2];
rz(2.4370952) q[2];
rz(0.38257515) q[3];
sx q[3];
rz(-1.4380598) q[3];
sx q[3];
rz(0.85324919) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.660897) q[0];
sx q[0];
rz(-0.040204164) q[0];
sx q[0];
rz(-0.30886343) q[0];
rz(2.1924696) q[1];
sx q[1];
rz(-2.821065) q[1];
sx q[1];
rz(2.5905051) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68252045) q[0];
sx q[0];
rz(-0.42051747) q[0];
sx q[0];
rz(0.38078749) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4818207) q[2];
sx q[2];
rz(-0.65976652) q[2];
sx q[2];
rz(-0.24487409) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6092564) q[1];
sx q[1];
rz(-0.21363959) q[1];
sx q[1];
rz(-2.3268301) q[1];
rz(-pi) q[2];
rz(-1.7150319) q[3];
sx q[3];
rz(-2.3420685) q[3];
sx q[3];
rz(2.1790847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.71517313) q[2];
sx q[2];
rz(-2.3931563) q[2];
sx q[2];
rz(0.61881649) q[2];
rz(0.62301451) q[3];
sx q[3];
rz(-2.2996733) q[3];
sx q[3];
rz(0.4272517) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039336786) q[0];
sx q[0];
rz(-2.4781041) q[0];
sx q[0];
rz(-0.46184194) q[0];
rz(2.6448008) q[1];
sx q[1];
rz(-1.8180314) q[1];
sx q[1];
rz(1.7740446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0858618) q[0];
sx q[0];
rz(-1.6501969) q[0];
sx q[0];
rz(-2.2532796) q[0];
rz(-0.15256792) q[2];
sx q[2];
rz(-1.0786453) q[2];
sx q[2];
rz(-0.14482611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.313844) q[1];
sx q[1];
rz(-1.5987458) q[1];
sx q[1];
rz(-1.8716783) q[1];
rz(0.9487834) q[3];
sx q[3];
rz(-1.7025196) q[3];
sx q[3];
rz(1.6831116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0355012) q[2];
sx q[2];
rz(-1.1320628) q[2];
sx q[2];
rz(-2.3512225) q[2];
rz(-2.5370989) q[3];
sx q[3];
rz(-2.1166182) q[3];
sx q[3];
rz(0.42009556) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6919493) q[0];
sx q[0];
rz(-0.88570166) q[0];
sx q[0];
rz(2.7272136) q[0];
rz(-0.62037933) q[1];
sx q[1];
rz(-1.9199771) q[1];
sx q[1];
rz(1.5827804) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5511502) q[0];
sx q[0];
rz(-1.4865459) q[0];
sx q[0];
rz(0.1131375) q[0];
rz(-pi) q[1];
rz(-1.5887268) q[2];
sx q[2];
rz(-1.6758783) q[2];
sx q[2];
rz(-1.1370657) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1852946) q[1];
sx q[1];
rz(-1.5130338) q[1];
sx q[1];
rz(0.0094780427) q[1];
x q[2];
rz(-1.9314194) q[3];
sx q[3];
rz(-0.40935959) q[3];
sx q[3];
rz(-1.6334074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15535007) q[2];
sx q[2];
rz(-0.68779951) q[2];
sx q[2];
rz(-0.17779329) q[2];
rz(-2.6627461) q[3];
sx q[3];
rz(-2.0279453) q[3];
sx q[3];
rz(3.0732885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9947522) q[0];
sx q[0];
rz(-0.9786334) q[0];
sx q[0];
rz(-0.11257182) q[0];
rz(0.62537891) q[1];
sx q[1];
rz(-2.0140779) q[1];
sx q[1];
rz(-0.88923997) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.205394) q[0];
sx q[0];
rz(-1.5279211) q[0];
sx q[0];
rz(-1.491719) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26548141) q[2];
sx q[2];
rz(-1.4944139) q[2];
sx q[2];
rz(2.7838109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4496838) q[1];
sx q[1];
rz(-1.050321) q[1];
sx q[1];
rz(0.27218216) q[1];
rz(-pi) q[2];
rz(0.067178848) q[3];
sx q[3];
rz(-1.805873) q[3];
sx q[3];
rz(2.377411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2832977) q[2];
sx q[2];
rz(-3.0136643) q[2];
sx q[2];
rz(-2.8683786) q[2];
rz(-2.9122399) q[3];
sx q[3];
rz(-1.554824) q[3];
sx q[3];
rz(1.2685512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-0.44723085) q[0];
sx q[0];
rz(-0.87089592) q[0];
sx q[0];
rz(0.91621512) q[0];
rz(-2.9592196) q[1];
sx q[1];
rz(-0.20621754) q[1];
sx q[1];
rz(1.1590385) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0908546) q[0];
sx q[0];
rz(-1.3373475) q[0];
sx q[0];
rz(1.3755193) q[0];
rz(0.49075895) q[2];
sx q[2];
rz(-0.10819866) q[2];
sx q[2];
rz(0.27709093) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68284833) q[1];
sx q[1];
rz(-0.35431752) q[1];
sx q[1];
rz(-1.09642) q[1];
rz(-pi) q[2];
rz(-0.69144122) q[3];
sx q[3];
rz(-0.85606282) q[3];
sx q[3];
rz(0.5428398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75659043) q[2];
sx q[2];
rz(-0.80005163) q[2];
sx q[2];
rz(0.315256) q[2];
rz(-0.59424019) q[3];
sx q[3];
rz(-2.2599594) q[3];
sx q[3];
rz(-0.63410223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88157982) q[0];
sx q[0];
rz(-0.6686815) q[0];
sx q[0];
rz(-0.15923937) q[0];
rz(2.0533994) q[1];
sx q[1];
rz(-2.2290778) q[1];
sx q[1];
rz(2.870627) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3581287) q[0];
sx q[0];
rz(-0.53345752) q[0];
sx q[0];
rz(-1.0607884) q[0];
x q[1];
rz(0.6728716) q[2];
sx q[2];
rz(-1.6511763) q[2];
sx q[2];
rz(2.5774235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9393255) q[1];
sx q[1];
rz(-1.9069873) q[1];
sx q[1];
rz(2.4656117) q[1];
rz(1.9064398) q[3];
sx q[3];
rz(-2.0291532) q[3];
sx q[3];
rz(3.1202849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8012041) q[2];
sx q[2];
rz(-0.77216721) q[2];
sx q[2];
rz(-0.32269746) q[2];
rz(0.05803756) q[3];
sx q[3];
rz(-0.80153424) q[3];
sx q[3];
rz(-0.29840741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67644453) q[0];
sx q[0];
rz(-1.5100751) q[0];
sx q[0];
rz(-0.81612192) q[0];
rz(2.8836518) q[1];
sx q[1];
rz(-1.1358658) q[1];
sx q[1];
rz(1.6075016) q[1];
rz(1.285124) q[2];
sx q[2];
rz(-1.2930047) q[2];
sx q[2];
rz(1.5563896) q[2];
rz(-0.54775379) q[3];
sx q[3];
rz(-2.4079101) q[3];
sx q[3];
rz(-1.1338601) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
