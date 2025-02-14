OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3352873) q[0];
sx q[0];
rz(-1.0399613) q[0];
sx q[0];
rz(2.7916173) q[0];
rz(1.5179874) q[1];
sx q[1];
rz(-2.1972158) q[1];
sx q[1];
rz(1.1854393) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0703263) q[0];
sx q[0];
rz(-1.2375273) q[0];
sx q[0];
rz(-0.44227483) q[0];
rz(-pi) q[1];
rz(-1.3385504) q[2];
sx q[2];
rz(-1.9793538) q[2];
sx q[2];
rz(-2.3790155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9361286) q[1];
sx q[1];
rz(-1.3359114) q[1];
sx q[1];
rz(-1.0297736) q[1];
rz(2.9583232) q[3];
sx q[3];
rz(-2.1784867) q[3];
sx q[3];
rz(0.83848467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27652937) q[2];
sx q[2];
rz(-2.7294366) q[2];
sx q[2];
rz(2.7318447) q[2];
rz(-2.3661738) q[3];
sx q[3];
rz(-1.6131468) q[3];
sx q[3];
rz(0.30737901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.2839117) q[0];
sx q[0];
rz(-2.6429521) q[0];
sx q[0];
rz(-2.6890802) q[0];
rz(-2.9127938) q[1];
sx q[1];
rz(-0.50140536) q[1];
sx q[1];
rz(-0.14678821) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7797557) q[0];
sx q[0];
rz(-2.728403) q[0];
sx q[0];
rz(-1.0784515) q[0];
rz(-pi) q[1];
rz(-0.37756672) q[2];
sx q[2];
rz(-1.7930328) q[2];
sx q[2];
rz(1.8309895) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9917069) q[1];
sx q[1];
rz(-0.14016166) q[1];
sx q[1];
rz(-0.28729846) q[1];
x q[2];
rz(0.5654752) q[3];
sx q[3];
rz(-2.3982022) q[3];
sx q[3];
rz(-0.88228031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16810922) q[2];
sx q[2];
rz(-2.8066469) q[2];
sx q[2];
rz(-2.2701263) q[2];
rz(0.35428366) q[3];
sx q[3];
rz(-2.2058217) q[3];
sx q[3];
rz(-1.3080477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63534921) q[0];
sx q[0];
rz(-0.62017089) q[0];
sx q[0];
rz(-2.035602) q[0];
rz(-2.1982819) q[1];
sx q[1];
rz(-2.3425075) q[1];
sx q[1];
rz(1.4897289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8603926) q[0];
sx q[0];
rz(-1.6985122) q[0];
sx q[0];
rz(0.010334947) q[0];
x q[1];
rz(2.9331384) q[2];
sx q[2];
rz(-1.4770419) q[2];
sx q[2];
rz(-1.9970279) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5691487) q[1];
sx q[1];
rz(-0.97578543) q[1];
sx q[1];
rz(-1.6165363) q[1];
rz(2.843804) q[3];
sx q[3];
rz(-0.2486767) q[3];
sx q[3];
rz(2.3657048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6279383) q[2];
sx q[2];
rz(-1.1955806) q[2];
sx q[2];
rz(0.99620831) q[2];
rz(-0.57322383) q[3];
sx q[3];
rz(-2.6300391) q[3];
sx q[3];
rz(-0.63428026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-1.2117677) q[0];
sx q[0];
rz(-1.8695762) q[0];
sx q[0];
rz(1.5856532) q[0];
rz(2.8171483) q[1];
sx q[1];
rz(-0.8322081) q[1];
sx q[1];
rz(-1.197804) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65069136) q[0];
sx q[0];
rz(-1.317121) q[0];
sx q[0];
rz(0.59532292) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3487071) q[2];
sx q[2];
rz(-2.54455) q[2];
sx q[2];
rz(-0.82840234) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9442695) q[1];
sx q[1];
rz(-1.0900647) q[1];
sx q[1];
rz(-2.9949147) q[1];
rz(-1.3088552) q[3];
sx q[3];
rz(-1.4125694) q[3];
sx q[3];
rz(-2.5867274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0962301) q[2];
sx q[2];
rz(-0.54916334) q[2];
sx q[2];
rz(1.9722923) q[2];
rz(-0.62396389) q[3];
sx q[3];
rz(-1.4493161) q[3];
sx q[3];
rz(0.72871488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28155926) q[0];
sx q[0];
rz(-2.7427854) q[0];
sx q[0];
rz(-1.1997724) q[0];
rz(-1.3072321) q[1];
sx q[1];
rz(-0.93991005) q[1];
sx q[1];
rz(1.4365546) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4685637) q[0];
sx q[0];
rz(-1.1832245) q[0];
sx q[0];
rz(0.065385575) q[0];
rz(-2.0827411) q[2];
sx q[2];
rz(-0.26341715) q[2];
sx q[2];
rz(-2.9204766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87203997) q[1];
sx q[1];
rz(-1.2721354) q[1];
sx q[1];
rz(2.7725817) q[1];
rz(-0.41884274) q[3];
sx q[3];
rz(-2.1757033) q[3];
sx q[3];
rz(2.1640282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78937078) q[2];
sx q[2];
rz(-1.7033966) q[2];
sx q[2];
rz(0.5729202) q[2];
rz(0.88417435) q[3];
sx q[3];
rz(-0.97777647) q[3];
sx q[3];
rz(-0.49527112) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5837412) q[0];
sx q[0];
rz(-1.8889677) q[0];
sx q[0];
rz(-1.0766693) q[0];
rz(-0.50648266) q[1];
sx q[1];
rz(-0.61512893) q[1];
sx q[1];
rz(1.3004251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69715766) q[0];
sx q[0];
rz(-1.5570672) q[0];
sx q[0];
rz(1.619471) q[0];
rz(-2.9175148) q[2];
sx q[2];
rz(-2.5505318) q[2];
sx q[2];
rz(3.0853809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9058076) q[1];
sx q[1];
rz(-1.7828568) q[1];
sx q[1];
rz(-3.0074988) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9415116) q[3];
sx q[3];
rz(-1.225138) q[3];
sx q[3];
rz(2.8576771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41636813) q[2];
sx q[2];
rz(-2.29795) q[2];
sx q[2];
rz(1.9492524) q[2];
rz(-0.013817712) q[3];
sx q[3];
rz(-2.4255224) q[3];
sx q[3];
rz(2.1684087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8927807) q[0];
sx q[0];
rz(-0.92614323) q[0];
sx q[0];
rz(-2.7698351) q[0];
rz(0.88428503) q[1];
sx q[1];
rz(-2.6339032) q[1];
sx q[1];
rz(2.5977792) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.666558) q[0];
sx q[0];
rz(-1.0109996) q[0];
sx q[0];
rz(2.0535357) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0804266) q[2];
sx q[2];
rz(-0.78713929) q[2];
sx q[2];
rz(2.1619881) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71173271) q[1];
sx q[1];
rz(-1.5641583) q[1];
sx q[1];
rz(0.37488519) q[1];
rz(-pi) q[2];
rz(-1.722808) q[3];
sx q[3];
rz(-1.6789888) q[3];
sx q[3];
rz(-3.0967847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1112993) q[2];
sx q[2];
rz(-1.0617504) q[2];
sx q[2];
rz(-0.57478762) q[2];
rz(2.287367) q[3];
sx q[3];
rz(-1.4762907) q[3];
sx q[3];
rz(-1.1150208) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2146571) q[0];
sx q[0];
rz(-0.68100005) q[0];
sx q[0];
rz(-2.8171203) q[0];
rz(0.60943162) q[1];
sx q[1];
rz(-2.29988) q[1];
sx q[1];
rz(0.51469523) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99583944) q[0];
sx q[0];
rz(-2.0059364) q[0];
sx q[0];
rz(3.1313722) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35867547) q[2];
sx q[2];
rz(-2.2513362) q[2];
sx q[2];
rz(0.84748617) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1140442) q[1];
sx q[1];
rz(-1.0107155) q[1];
sx q[1];
rz(1.8482988) q[1];
rz(1.8466161) q[3];
sx q[3];
rz(-1.4657186) q[3];
sx q[3];
rz(0.17227473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9784166) q[2];
sx q[2];
rz(-1.6360444) q[2];
sx q[2];
rz(-2.8669538) q[2];
rz(-0.93242532) q[3];
sx q[3];
rz(-2.7115188) q[3];
sx q[3];
rz(2.8453804) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68359971) q[0];
sx q[0];
rz(-1.4608811) q[0];
sx q[0];
rz(3.006111) q[0];
rz(0.97700351) q[1];
sx q[1];
rz(-0.75825399) q[1];
sx q[1];
rz(-0.0010842222) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1019395) q[0];
sx q[0];
rz(-2.6216767) q[0];
sx q[0];
rz(2.824409) q[0];
x q[1];
rz(3.0227109) q[2];
sx q[2];
rz(-2.0914051) q[2];
sx q[2];
rz(-0.05725459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75669589) q[1];
sx q[1];
rz(-0.86844173) q[1];
sx q[1];
rz(-2.3780502) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67768456) q[3];
sx q[3];
rz(-2.2608058) q[3];
sx q[3];
rz(-3.1180988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5486501) q[2];
sx q[2];
rz(-1.0385916) q[2];
sx q[2];
rz(-0.20498094) q[2];
rz(0.18260469) q[3];
sx q[3];
rz(-0.20668106) q[3];
sx q[3];
rz(-2.4038147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36271998) q[0];
sx q[0];
rz(-2.7474779) q[0];
sx q[0];
rz(0.47738281) q[0];
rz(3.0248771) q[1];
sx q[1];
rz(-1.621403) q[1];
sx q[1];
rz(2.3027072) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8461743) q[0];
sx q[0];
rz(-1.3059761) q[0];
sx q[0];
rz(2.2841262) q[0];
rz(2.7493189) q[2];
sx q[2];
rz(-0.76425154) q[2];
sx q[2];
rz(-1.6583673) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.15167639) q[1];
sx q[1];
rz(-0.71571678) q[1];
sx q[1];
rz(0.28629656) q[1];
rz(-3.1118591) q[3];
sx q[3];
rz(-1.8095867) q[3];
sx q[3];
rz(-0.58837552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78187752) q[2];
sx q[2];
rz(-2.8350267) q[2];
sx q[2];
rz(0.26620418) q[2];
rz(2.6340458) q[3];
sx q[3];
rz(-0.72269732) q[3];
sx q[3];
rz(2.7005196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97892852) q[0];
sx q[0];
rz(-1.6829818) q[0];
sx q[0];
rz(-0.80243954) q[0];
rz(-0.91213999) q[1];
sx q[1];
rz(-1.1911387) q[1];
sx q[1];
rz(-2.3008507) q[1];
rz(0.56923127) q[2];
sx q[2];
rz(-0.29217671) q[2];
sx q[2];
rz(3.0859699) q[2];
rz(1.1900564) q[3];
sx q[3];
rz(-1.1968812) q[3];
sx q[3];
rz(1.4637917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
