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
rz(0.0042301099) q[0];
sx q[0];
rz(-0.63679305) q[0];
sx q[0];
rz(-1.1871583) q[0];
rz(4.123426) q[1];
sx q[1];
rz(0.47458664) q[1];
sx q[1];
rz(8.4835806) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.610674) q[0];
sx q[0];
rz(-2.0661021) q[0];
sx q[0];
rz(-1.5688867) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12856483) q[2];
sx q[2];
rz(-1.7714273) q[2];
sx q[2];
rz(3.0561997) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8601712) q[1];
sx q[1];
rz(-1.0152536) q[1];
sx q[1];
rz(-1.9940328) q[1];
x q[2];
rz(-1.7205129) q[3];
sx q[3];
rz(-0.46105534) q[3];
sx q[3];
rz(2.113101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59419751) q[2];
sx q[2];
rz(-1.6387458) q[2];
sx q[2];
rz(2.1217864) q[2];
rz(2.8504573) q[3];
sx q[3];
rz(-2.9499493) q[3];
sx q[3];
rz(-1.8118793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52784598) q[0];
sx q[0];
rz(-2.3044523) q[0];
sx q[0];
rz(-1.5035195) q[0];
rz(1.232052) q[1];
sx q[1];
rz(-0.82545311) q[1];
sx q[1];
rz(1.253461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37663662) q[0];
sx q[0];
rz(-0.47237637) q[0];
sx q[0];
rz(1.0563904) q[0];
x q[1];
rz(0.49053662) q[2];
sx q[2];
rz(-0.56714688) q[2];
sx q[2];
rz(-0.74037742) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4644891) q[1];
sx q[1];
rz(-1.3687412) q[1];
sx q[1];
rz(1.7176164) q[1];
x q[2];
rz(-2.3971629) q[3];
sx q[3];
rz(-1.9081313) q[3];
sx q[3];
rz(2.6312089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2848009) q[2];
sx q[2];
rz(-2.5992726) q[2];
sx q[2];
rz(0.33033237) q[2];
rz(3.0387943) q[3];
sx q[3];
rz(-1.8967352) q[3];
sx q[3];
rz(-2.5225294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61515808) q[0];
sx q[0];
rz(-2.0005083) q[0];
sx q[0];
rz(-1.0636299) q[0];
rz(0.10387736) q[1];
sx q[1];
rz(-2.3492298) q[1];
sx q[1];
rz(1.1054543) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55872592) q[0];
sx q[0];
rz(-1.6028645) q[0];
sx q[0];
rz(1.523287) q[0];
rz(-1.8806522) q[2];
sx q[2];
rz(-0.83854691) q[2];
sx q[2];
rz(-2.1149879) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7522831) q[1];
sx q[1];
rz(-2.0842617) q[1];
sx q[1];
rz(-0.13996831) q[1];
rz(1.5289588) q[3];
sx q[3];
rz(-1.5825019) q[3];
sx q[3];
rz(-0.042543471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1198472) q[2];
sx q[2];
rz(-2.3357119) q[2];
sx q[2];
rz(0.22551192) q[2];
rz(-1.3269904) q[3];
sx q[3];
rz(-0.83695379) q[3];
sx q[3];
rz(2.5007611) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2284112) q[0];
sx q[0];
rz(-1.5786194) q[0];
sx q[0];
rz(2.5033503) q[0];
rz(-1.111521) q[1];
sx q[1];
rz(-0.94739729) q[1];
sx q[1];
rz(-0.18724719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3950662) q[0];
sx q[0];
rz(-1.2842261) q[0];
sx q[0];
rz(2.2118083) q[0];
rz(-pi) q[1];
rz(-0.78088607) q[2];
sx q[2];
rz(-1.5263288) q[2];
sx q[2];
rz(-1.5574297) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.79689246) q[1];
sx q[1];
rz(-2.6019367) q[1];
sx q[1];
rz(-0.15426387) q[1];
x q[2];
rz(1.4832014) q[3];
sx q[3];
rz(-2.2822652) q[3];
sx q[3];
rz(0.30184823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1049261) q[2];
sx q[2];
rz(-1.5035524) q[2];
sx q[2];
rz(0.46910134) q[2];
rz(-2.8406298) q[3];
sx q[3];
rz(-2.8807237) q[3];
sx q[3];
rz(-1.652396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24333532) q[0];
sx q[0];
rz(-1.5168334) q[0];
sx q[0];
rz(0.54018706) q[0];
rz(0.46145269) q[1];
sx q[1];
rz(-1.6199473) q[1];
sx q[1];
rz(2.8823421) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40045462) q[0];
sx q[0];
rz(-1.8807066) q[0];
sx q[0];
rz(-1.8551338) q[0];
x q[1];
rz(-2.2770693) q[2];
sx q[2];
rz(-2.5062423) q[2];
sx q[2];
rz(-2.8453448) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.66189761) q[1];
sx q[1];
rz(-0.59407083) q[1];
sx q[1];
rz(-1.8516783) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9018039) q[3];
sx q[3];
rz(-1.3223303) q[3];
sx q[3];
rz(2.1295761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3429012) q[2];
sx q[2];
rz(-2.4553802) q[2];
sx q[2];
rz(0.9101103) q[2];
rz(-2.5765007) q[3];
sx q[3];
rz(-1.3684973) q[3];
sx q[3];
rz(2.4250987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27123505) q[0];
sx q[0];
rz(-0.39327455) q[0];
sx q[0];
rz(-1.2275335) q[0];
rz(-0.50499376) q[1];
sx q[1];
rz(-2.125449) q[1];
sx q[1];
rz(1.9379157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6617463) q[0];
sx q[0];
rz(-1.8176314) q[0];
sx q[0];
rz(-1.6936452) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2199322) q[2];
sx q[2];
rz(-1.9699012) q[2];
sx q[2];
rz(1.6726577) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5745637) q[1];
sx q[1];
rz(-2.1334799) q[1];
sx q[1];
rz(2.3966054) q[1];
rz(-1.1111383) q[3];
sx q[3];
rz(-1.3317923) q[3];
sx q[3];
rz(0.29005946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0248854) q[2];
sx q[2];
rz(-2.1163157) q[2];
sx q[2];
rz(-0.95899478) q[2];
rz(2.8420908) q[3];
sx q[3];
rz(-3.0193269) q[3];
sx q[3];
rz(-1.6925156) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9129979) q[0];
sx q[0];
rz(-1.2499502) q[0];
sx q[0];
rz(0.49384299) q[0];
rz(0.85810581) q[1];
sx q[1];
rz(-1.162642) q[1];
sx q[1];
rz(1.3546622) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1345024) q[0];
sx q[0];
rz(-2.1190152) q[0];
sx q[0];
rz(-2.8584418) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0042393) q[2];
sx q[2];
rz(-1.0826515) q[2];
sx q[2];
rz(3.1078095) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63039102) q[1];
sx q[1];
rz(-1.5005209) q[1];
sx q[1];
rz(-2.3831638) q[1];
x q[2];
rz(2.2432547) q[3];
sx q[3];
rz(-1.1500937) q[3];
sx q[3];
rz(1.4843693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4570423) q[2];
sx q[2];
rz(-0.6251308) q[2];
sx q[2];
rz(-0.54614145) q[2];
rz(-2.9438733) q[3];
sx q[3];
rz(-1.0309912) q[3];
sx q[3];
rz(-0.24699591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2208743) q[0];
sx q[0];
rz(-1.1552759) q[0];
sx q[0];
rz(-2.0350631) q[0];
rz(0.58440343) q[1];
sx q[1];
rz(-1.1738651) q[1];
sx q[1];
rz(2.1388785) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2276777) q[0];
sx q[0];
rz(-1.4105689) q[0];
sx q[0];
rz(-2.0100309) q[0];
rz(-pi) q[1];
rz(-0.8876351) q[2];
sx q[2];
rz(-2.7681367) q[2];
sx q[2];
rz(1.68242) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4312736) q[1];
sx q[1];
rz(-2.6280676) q[1];
sx q[1];
rz(-3.0679614) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9225305) q[3];
sx q[3];
rz(-0.98120839) q[3];
sx q[3];
rz(-2.0812579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5661085) q[2];
sx q[2];
rz(-0.9245975) q[2];
sx q[2];
rz(-1.7522579) q[2];
rz(-0.53650457) q[3];
sx q[3];
rz(-1.5080695) q[3];
sx q[3];
rz(2.2947521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.093512) q[0];
sx q[0];
rz(-0.41963136) q[0];
sx q[0];
rz(2.804948) q[0];
rz(-0.078314217) q[1];
sx q[1];
rz(-1.2093733) q[1];
sx q[1];
rz(1.9858817) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5040249) q[0];
sx q[0];
rz(-0.89533606) q[0];
sx q[0];
rz(-0.43677377) q[0];
x q[1];
rz(1.1572273) q[2];
sx q[2];
rz(-1.6317657) q[2];
sx q[2];
rz(-2.3496036) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3382321) q[1];
sx q[1];
rz(-0.50932099) q[1];
sx q[1];
rz(-1.2912911) q[1];
rz(-pi) q[2];
x q[2];
rz(0.09241032) q[3];
sx q[3];
rz(-0.31805491) q[3];
sx q[3];
rz(0.89124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.081850514) q[2];
sx q[2];
rz(-2.0589224) q[2];
sx q[2];
rz(-0.82320631) q[2];
rz(-2.1246223) q[3];
sx q[3];
rz(-2.7401994) q[3];
sx q[3];
rz(3.1034191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80242771) q[0];
sx q[0];
rz(-2.9186354) q[0];
sx q[0];
rz(1.7145351) q[0];
rz(1.4786221) q[1];
sx q[1];
rz(-1.0717816) q[1];
sx q[1];
rz(-2.2907385) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1647536) q[0];
sx q[0];
rz(-1.098159) q[0];
sx q[0];
rz(-3.0125299) q[0];
x q[1];
rz(3.0548373) q[2];
sx q[2];
rz(-1.6977001) q[2];
sx q[2];
rz(-2.639304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6809594) q[1];
sx q[1];
rz(-0.52648989) q[1];
sx q[1];
rz(0.10641702) q[1];
rz(0.14031841) q[3];
sx q[3];
rz(-1.6346953) q[3];
sx q[3];
rz(1.7238703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.458638) q[2];
sx q[2];
rz(-1.3502324) q[2];
sx q[2];
rz(0.82009912) q[2];
rz(-1.5470777) q[3];
sx q[3];
rz(-1.4328512) q[3];
sx q[3];
rz(-1.1944176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9783258) q[0];
sx q[0];
rz(-1.5112725) q[0];
sx q[0];
rz(1.5164966) q[0];
rz(-0.77240472) q[1];
sx q[1];
rz(-2.518387) q[1];
sx q[1];
rz(1.0060681) q[1];
rz(-0.27323451) q[2];
sx q[2];
rz(-2.4315576) q[2];
sx q[2];
rz(-0.72080056) q[2];
rz(-2.3631572) q[3];
sx q[3];
rz(-0.62494559) q[3];
sx q[3];
rz(0.86452053) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
