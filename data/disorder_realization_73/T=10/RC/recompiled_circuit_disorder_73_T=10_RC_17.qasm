OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(-2.6383658) q[0];
sx q[0];
rz(0.72416645) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(0.78483265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1278348) q[0];
sx q[0];
rz(-2.7779967) q[0];
sx q[0];
rz(-2.512393) q[0];
x q[1];
rz(-0.22462331) q[2];
sx q[2];
rz(-2.7135239) q[2];
sx q[2];
rz(0.12878865) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29014978) q[1];
sx q[1];
rz(-2.5622497) q[1];
sx q[1];
rz(1.1555374) q[1];
x q[2];
rz(-1.9853398) q[3];
sx q[3];
rz(-1.6016377) q[3];
sx q[3];
rz(2.8950092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-0.067967728) q[2];
rz(-3.0170278) q[3];
sx q[3];
rz(-2.8187276) q[3];
sx q[3];
rz(-1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2215866) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(2.8979229) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(1.3557281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97665635) q[0];
sx q[0];
rz(-1.5979842) q[0];
sx q[0];
rz(-1.8261766) q[0];
rz(-pi) q[1];
rz(2.3089356) q[2];
sx q[2];
rz(-2.6359733) q[2];
sx q[2];
rz(-0.51072272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1214952) q[1];
sx q[1];
rz(-1.1250245) q[1];
sx q[1];
rz(0.14264588) q[1];
rz(-pi) q[2];
rz(-0.59216604) q[3];
sx q[3];
rz(-1.5321931) q[3];
sx q[3];
rz(-1.5241227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(2.8919354) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8963985) q[0];
sx q[0];
rz(-1.9165374) q[0];
sx q[0];
rz(2.2431592) q[0];
rz(1.8067182) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(-1.2737087) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0262336) q[0];
sx q[0];
rz(-2.6991978) q[0];
sx q[0];
rz(2.1917403) q[0];
rz(-2.8286335) q[2];
sx q[2];
rz(-2.8319781) q[2];
sx q[2];
rz(1.5086053) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.102036) q[1];
sx q[1];
rz(-1.9196871) q[1];
sx q[1];
rz(0.7015014) q[1];
x q[2];
rz(2.418163) q[3];
sx q[3];
rz(-1.1838786) q[3];
sx q[3];
rz(0.5667516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0597824) q[2];
sx q[2];
rz(-2.5368097) q[2];
sx q[2];
rz(-2.188142) q[2];
rz(0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(-2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.8614486) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(2.8934073) q[0];
rz(-1.0379627) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(0.074137069) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0291245) q[0];
sx q[0];
rz(-1.0974786) q[0];
sx q[0];
rz(0.28664696) q[0];
x q[1];
rz(-0.89216994) q[2];
sx q[2];
rz(-1.2954419) q[2];
sx q[2];
rz(-2.8868669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3809507) q[1];
sx q[1];
rz(-1.6886854) q[1];
sx q[1];
rz(-1.8878493) q[1];
x q[2];
rz(0.24012633) q[3];
sx q[3];
rz(-1.4014981) q[3];
sx q[3];
rz(-1.5907767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-0.038643535) q[2];
rz(2.1679227) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(-2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(1.3624396) q[0];
rz(0.81659395) q[1];
sx q[1];
rz(-1.8530308) q[1];
sx q[1];
rz(1.978925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1122738) q[0];
sx q[0];
rz(-0.11538878) q[0];
sx q[0];
rz(-1.3089887) q[0];
x q[1];
rz(1.2478998) q[2];
sx q[2];
rz(-1.731551) q[2];
sx q[2];
rz(1.2231183) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.34449023) q[1];
sx q[1];
rz(-1.7624197) q[1];
sx q[1];
rz(-1.2002798) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88697042) q[3];
sx q[3];
rz(-2.5081222) q[3];
sx q[3];
rz(-0.57852832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(0.14222063) q[2];
rz(-2.2375315) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(2.952125) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11480039) q[0];
sx q[0];
rz(-0.27652201) q[0];
sx q[0];
rz(-1.4676771) q[0];
rz(-2.5698075) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(-0.30803672) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80297563) q[0];
sx q[0];
rz(-2.9836285) q[0];
sx q[0];
rz(-0.64361848) q[0];
rz(2.5295837) q[2];
sx q[2];
rz(-2.0888121) q[2];
sx q[2];
rz(2.3305364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.249237) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(-1.874079) q[1];
rz(2.7932348) q[3];
sx q[3];
rz(-1.6568686) q[3];
sx q[3];
rz(-1.100988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77928153) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(1.1479088) q[2];
rz(2.4273196) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(2.4961297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797453) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(-3.0116144) q[0];
rz(3.1107483) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(-0.67108363) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13043159) q[0];
sx q[0];
rz(-1.6058308) q[0];
sx q[0];
rz(-3.0879211) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7332156) q[2];
sx q[2];
rz(-0.186609) q[2];
sx q[2];
rz(-2.9090372) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.259686) q[1];
sx q[1];
rz(-0.38609186) q[1];
sx q[1];
rz(-1.2950456) q[1];
rz(0.89550771) q[3];
sx q[3];
rz(-1.5338147) q[3];
sx q[3];
rz(2.0946338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0188296) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(-3.1130062) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9776483) q[0];
sx q[0];
rz(-2.2387235) q[0];
sx q[0];
rz(-0.41241616) q[0];
rz(1.6917797) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(1.9746045) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8186504) q[0];
sx q[0];
rz(-2.4837821) q[0];
sx q[0];
rz(-1.2176745) q[0];
rz(-1.251986) q[2];
sx q[2];
rz(-2.406771) q[2];
sx q[2];
rz(0.97359818) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1046909) q[1];
sx q[1];
rz(-0.19100405) q[1];
sx q[1];
rz(0.16049338) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9917631) q[3];
sx q[3];
rz(-1.0458046) q[3];
sx q[3];
rz(-2.6168952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-2.9525625) q[2];
rz(-0.14686251) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(1.3930901) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1259595) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(2.2145859) q[0];
rz(-1.3828297) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(1.4896726) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4179174) q[0];
sx q[0];
rz(-2.0360887) q[0];
sx q[0];
rz(-2.377541) q[0];
x q[1];
rz(1.4719047) q[2];
sx q[2];
rz(-0.47404587) q[2];
sx q[2];
rz(-0.74741077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5786963) q[1];
sx q[1];
rz(-1.3778731) q[1];
sx q[1];
rz(-1.3955411) q[1];
rz(0.82724039) q[3];
sx q[3];
rz(-2.4628371) q[3];
sx q[3];
rz(-0.37952207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5902517) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(1.7112973) q[2];
rz(2.5643505) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(-1.757471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5532613) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(-0.28840315) q[0];
rz(2.6092031) q[1];
sx q[1];
rz(-2.6817697) q[1];
sx q[1];
rz(0.14702252) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64919103) q[0];
sx q[0];
rz(-2.1426139) q[0];
sx q[0];
rz(-2.6834821) q[0];
rz(0.23994259) q[2];
sx q[2];
rz(-1.502617) q[2];
sx q[2];
rz(-0.49382526) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.22132561) q[1];
sx q[1];
rz(-2.1888071) q[1];
sx q[1];
rz(-0.048939181) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5649928) q[3];
sx q[3];
rz(-1.4398265) q[3];
sx q[3];
rz(0.031158202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0599351) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(-2.3975513) q[2];
rz(2.384281) q[3];
sx q[3];
rz(-1.363874) q[3];
sx q[3];
rz(1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1159146) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(-2.3241282) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(-0.4521162) q[2];
sx q[2];
rz(-1.5099667) q[2];
sx q[2];
rz(-2.920334) q[2];
rz(-0.13027262) q[3];
sx q[3];
rz(-1.0126922) q[3];
sx q[3];
rz(-1.0425413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];