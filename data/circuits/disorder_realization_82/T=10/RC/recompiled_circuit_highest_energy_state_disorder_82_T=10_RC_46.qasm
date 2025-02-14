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
rz(3.1458228) q[0];
sx q[0];
rz(0.63679305) q[0];
sx q[0];
rz(11.379212) q[0];
rz(-2.1597593) q[1];
sx q[1];
rz(-2.667006) q[1];
sx q[1];
rz(0.9411974) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1008074) q[0];
sx q[0];
rz(-1.5724764) q[0];
sx q[0];
rz(0.4953065) q[0];
x q[1];
rz(3.0130278) q[2];
sx q[2];
rz(-1.7714273) q[2];
sx q[2];
rz(-3.0561997) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8601712) q[1];
sx q[1];
rz(-1.0152536) q[1];
sx q[1];
rz(1.1475599) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7205129) q[3];
sx q[3];
rz(-0.46105534) q[3];
sx q[3];
rz(1.0284916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59419751) q[2];
sx q[2];
rz(-1.6387458) q[2];
sx q[2];
rz(2.1217864) q[2];
rz(-2.8504573) q[3];
sx q[3];
rz(-0.19164339) q[3];
sx q[3];
rz(1.3297133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52784598) q[0];
sx q[0];
rz(-0.83714038) q[0];
sx q[0];
rz(1.5035195) q[0];
rz(1.9095406) q[1];
sx q[1];
rz(-0.82545311) q[1];
sx q[1];
rz(-1.253461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37663662) q[0];
sx q[0];
rz(-0.47237637) q[0];
sx q[0];
rz(-2.0852023) q[0];
rz(-0.5118891) q[2];
sx q[2];
rz(-1.3149258) q[2];
sx q[2];
rz(2.7344109) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67710357) q[1];
sx q[1];
rz(-1.7728514) q[1];
sx q[1];
rz(1.7176164) q[1];
rz(-2.6639179) q[3];
sx q[3];
rz(-0.80373418) q[3];
sx q[3];
rz(-0.71550575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2848009) q[2];
sx q[2];
rz(-0.54232001) q[2];
sx q[2];
rz(-2.8112603) q[2];
rz(0.10279837) q[3];
sx q[3];
rz(-1.8967352) q[3];
sx q[3];
rz(-0.61906329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61515808) q[0];
sx q[0];
rz(-1.1410843) q[0];
sx q[0];
rz(2.0779628) q[0];
rz(-0.10387736) q[1];
sx q[1];
rz(-2.3492298) q[1];
sx q[1];
rz(2.0361384) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5362107) q[0];
sx q[0];
rz(-3.0842801) q[0];
sx q[0];
rz(2.1648699) q[0];
rz(0.75656192) q[2];
sx q[2];
rz(-1.7995477) q[2];
sx q[2];
rz(-2.808266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.11239219) q[1];
sx q[1];
rz(-1.6926188) q[1];
sx q[1];
rz(-1.0531154) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8436935) q[3];
sx q[3];
rz(-0.04344329) q[3];
sx q[3];
rz(1.3406875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1198472) q[2];
sx q[2];
rz(-0.80588078) q[2];
sx q[2];
rz(0.22551192) q[2];
rz(1.3269904) q[3];
sx q[3];
rz(-0.83695379) q[3];
sx q[3];
rz(-2.5007611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2284112) q[0];
sx q[0];
rz(-1.5629733) q[0];
sx q[0];
rz(2.5033503) q[0];
rz(-1.111521) q[1];
sx q[1];
rz(-2.1941954) q[1];
sx q[1];
rz(-2.9543455) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1094681) q[0];
sx q[0];
rz(-0.95987849) q[0];
sx q[0];
rz(0.35232589) q[0];
x q[1];
rz(-0.78088607) q[2];
sx q[2];
rz(-1.5263288) q[2];
sx q[2];
rz(-1.5574297) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3447002) q[1];
sx q[1];
rz(-2.6019367) q[1];
sx q[1];
rz(-0.15426387) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6583913) q[3];
sx q[3];
rz(-2.2822652) q[3];
sx q[3];
rz(2.8397444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1049261) q[2];
sx q[2];
rz(-1.6380402) q[2];
sx q[2];
rz(-2.6724913) q[2];
rz(2.8406298) q[3];
sx q[3];
rz(-0.260869) q[3];
sx q[3];
rz(-1.652396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24333532) q[0];
sx q[0];
rz(-1.5168334) q[0];
sx q[0];
rz(-2.6014056) q[0];
rz(0.46145269) q[1];
sx q[1];
rz(-1.6199473) q[1];
sx q[1];
rz(2.8823421) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8823555) q[0];
sx q[0];
rz(-1.8412388) q[0];
sx q[0];
rz(-2.8195802) q[0];
rz(-0.86452338) q[2];
sx q[2];
rz(-2.5062423) q[2];
sx q[2];
rz(2.8453448) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67422489) q[1];
sx q[1];
rz(-1.4150054) q[1];
sx q[1];
rz(0.99512659) q[1];
x q[2];
rz(-0.90796538) q[3];
sx q[3];
rz(-0.41112152) q[3];
sx q[3];
rz(-3.079252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.79869142) q[2];
sx q[2];
rz(-2.4553802) q[2];
sx q[2];
rz(-2.2314824) q[2];
rz(2.5765007) q[3];
sx q[3];
rz(-1.3684973) q[3];
sx q[3];
rz(0.71649396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703576) q[0];
sx q[0];
rz(-2.7483181) q[0];
sx q[0];
rz(-1.9140592) q[0];
rz(-2.6365989) q[1];
sx q[1];
rz(-1.0161437) q[1];
sx q[1];
rz(1.9379157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01194309) q[0];
sx q[0];
rz(-2.8664358) q[0];
sx q[0];
rz(2.6889474) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9616379) q[2];
sx q[2];
rz(-0.74660245) q[2];
sx q[2];
rz(2.7701829) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5281637) q[1];
sx q[1];
rz(-0.89980308) q[1];
sx q[1];
rz(0.74929418) q[1];
rz(-pi) q[2];
rz(2.0730482) q[3];
sx q[3];
rz(-2.6274791) q[3];
sx q[3];
rz(0.83465605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0248854) q[2];
sx q[2];
rz(-2.1163157) q[2];
sx q[2];
rz(-2.1825979) q[2];
rz(-2.8420908) q[3];
sx q[3];
rz(-0.12226573) q[3];
sx q[3];
rz(-1.6925156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9129979) q[0];
sx q[0];
rz(-1.8916425) q[0];
sx q[0];
rz(2.6477497) q[0];
rz(-2.2834868) q[1];
sx q[1];
rz(-1.162642) q[1];
sx q[1];
rz(1.3546622) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49789594) q[0];
sx q[0];
rz(-2.5312811) q[0];
sx q[0];
rz(-1.1417139) q[0];
x q[1];
rz(-0.13735339) q[2];
sx q[2];
rz(-2.0589411) q[2];
sx q[2];
rz(3.1078095) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2676182) q[1];
sx q[1];
rz(-0.81470352) q[1];
sx q[1];
rz(-1.6674629) q[1];
x q[2];
rz(-2.6220697) q[3];
sx q[3];
rz(-2.1755927) q[3];
sx q[3];
rz(2.7407569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.4570423) q[2];
sx q[2];
rz(-2.5164618) q[2];
sx q[2];
rz(0.54614145) q[2];
rz(0.1977194) q[3];
sx q[3];
rz(-1.0309912) q[3];
sx q[3];
rz(2.8945967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2208743) q[0];
sx q[0];
rz(-1.9863167) q[0];
sx q[0];
rz(-1.1065296) q[0];
rz(-0.58440343) q[1];
sx q[1];
rz(-1.9677275) q[1];
sx q[1];
rz(2.1388785) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015771078) q[0];
sx q[0];
rz(-2.6758411) q[0];
sx q[0];
rz(-1.9339824) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8658337) q[2];
sx q[2];
rz(-1.8031839) q[2];
sx q[2];
rz(0.7600998) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9379077) q[1];
sx q[1];
rz(-1.6069429) q[1];
sx q[1];
rz(2.6292277) q[1];
rz(-pi) q[2];
rz(0.47559019) q[3];
sx q[3];
rz(-0.67567837) q[3];
sx q[3];
rz(-1.4978564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57548412) q[2];
sx q[2];
rz(-2.2169952) q[2];
sx q[2];
rz(1.3893348) q[2];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0480807) q[0];
sx q[0];
rz(-0.41963136) q[0];
sx q[0];
rz(-0.33664465) q[0];
rz(3.0632784) q[1];
sx q[1];
rz(-1.2093733) q[1];
sx q[1];
rz(1.9858817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2172359) q[0];
sx q[0];
rz(-1.2343533) q[0];
sx q[0];
rz(2.2948059) q[0];
x q[1];
rz(-1.9843654) q[2];
sx q[2];
rz(-1.509827) q[2];
sx q[2];
rz(2.3496036) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1286238) q[1];
sx q[1];
rz(-1.4358724) q[1];
sx q[1];
rz(1.0781481) q[1];
rz(-pi) q[2];
rz(-1.6011681) q[3];
sx q[3];
rz(-1.8874468) q[3];
sx q[3];
rz(-0.79398549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.081850514) q[2];
sx q[2];
rz(-2.0589224) q[2];
sx q[2];
rz(-2.3183863) q[2];
rz(-1.0169704) q[3];
sx q[3];
rz(-2.7401994) q[3];
sx q[3];
rz(0.038173525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3391649) q[0];
sx q[0];
rz(-0.22295727) q[0];
sx q[0];
rz(-1.7145351) q[0];
rz(1.4786221) q[1];
sx q[1];
rz(-1.0717816) q[1];
sx q[1];
rz(0.85085416) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9768391) q[0];
sx q[0];
rz(-1.098159) q[0];
sx q[0];
rz(-0.12906277) q[0];
x q[1];
rz(-2.1673604) q[2];
sx q[2];
rz(-0.15359226) q[2];
sx q[2];
rz(-0.099791137) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58356279) q[1];
sx q[1];
rz(-2.0940015) q[1];
sx q[1];
rz(-1.5091404) q[1];
rz(-pi) q[2];
rz(0.42907866) q[3];
sx q[3];
rz(-2.9874969) q[3];
sx q[3];
rz(-0.27151268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.458638) q[2];
sx q[2];
rz(-1.7913603) q[2];
sx q[2];
rz(-0.82009912) q[2];
rz(1.594515) q[3];
sx q[3];
rz(-1.4328512) q[3];
sx q[3];
rz(1.9471751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-2.3691879) q[1];
sx q[1];
rz(-0.62320566) q[1];
sx q[1];
rz(-2.1355245) q[1];
rz(0.27323451) q[2];
sx q[2];
rz(-0.71003503) q[2];
sx q[2];
rz(2.4207921) q[2];
rz(-1.1019273) q[3];
sx q[3];
rz(-2.0004604) q[3];
sx q[3];
rz(-0.018044005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
