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
rz(0.98183331) q[1];
sx q[1];
rz(-0.47458664) q[1];
sx q[1];
rz(2.2003953) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5309187) q[0];
sx q[0];
rz(-1.0754906) q[0];
sx q[0];
rz(-1.5727059) q[0];
rz(-pi) q[1];
rz(-1.0082939) q[2];
sx q[2];
rz(-0.23782158) q[2];
sx q[2];
rz(0.49006762) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.57454771) q[1];
sx q[1];
rz(-0.68463782) q[1];
sx q[1];
rz(2.5570832) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4210798) q[3];
sx q[3];
rz(-0.46105534) q[3];
sx q[3];
rz(1.0284916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59419751) q[2];
sx q[2];
rz(-1.6387458) q[2];
sx q[2];
rz(-2.1217864) q[2];
rz(-2.8504573) q[3];
sx q[3];
rz(-2.9499493) q[3];
sx q[3];
rz(1.8118793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137467) q[0];
sx q[0];
rz(-0.83714038) q[0];
sx q[0];
rz(-1.6380731) q[0];
rz(-1.9095406) q[1];
sx q[1];
rz(-2.3161395) q[1];
sx q[1];
rz(-1.253461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4811752) q[0];
sx q[0];
rz(-1.345013) q[0];
sx q[0];
rz(-1.9893454) q[0];
rz(-2.651056) q[2];
sx q[2];
rz(-0.56714688) q[2];
sx q[2];
rz(-0.74037742) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3121705) q[1];
sx q[1];
rz(-2.8924184) q[1];
sx q[1];
rz(0.62015688) q[1];
rz(-0.4776748) q[3];
sx q[3];
rz(-0.80373418) q[3];
sx q[3];
rz(0.71550575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85679179) q[2];
sx q[2];
rz(-0.54232001) q[2];
sx q[2];
rz(0.33033237) q[2];
rz(-0.10279837) q[3];
sx q[3];
rz(-1.2448575) q[3];
sx q[3];
rz(-0.61906329) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61515808) q[0];
sx q[0];
rz(-1.1410843) q[0];
sx q[0];
rz(2.0779628) q[0];
rz(3.0377153) q[1];
sx q[1];
rz(-2.3492298) q[1];
sx q[1];
rz(-1.1054543) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.010546) q[0];
sx q[0];
rz(-1.6182812) q[0];
sx q[0];
rz(-0.032104339) q[0];
rz(-pi) q[1];
rz(-1.2609405) q[2];
sx q[2];
rz(-2.3030457) q[2];
sx q[2];
rz(-2.1149879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0292005) q[1];
sx q[1];
rz(-1.6926188) q[1];
sx q[1];
rz(-2.0884772) q[1];
rz(3.1298769) q[3];
sx q[3];
rz(-1.612631) q[3];
sx q[3];
rz(-1.5277629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1198472) q[2];
sx q[2];
rz(-2.3357119) q[2];
sx q[2];
rz(-0.22551192) q[2];
rz(-1.8146023) q[3];
sx q[3];
rz(-2.3046389) q[3];
sx q[3];
rz(2.5007611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2284112) q[0];
sx q[0];
rz(-1.5786194) q[0];
sx q[0];
rz(-2.5033503) q[0];
rz(-2.0300716) q[1];
sx q[1];
rz(-2.1941954) q[1];
sx q[1];
rz(2.9543455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7465265) q[0];
sx q[0];
rz(-1.8573666) q[0];
sx q[0];
rz(2.2118083) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0784628) q[2];
sx q[2];
rz(-2.3597096) q[2];
sx q[2];
rz(3.1101335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.9065215) q[1];
sx q[1];
rz(-1.6498317) q[1];
sx q[1];
rz(-2.6071878) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10113206) q[3];
sx q[3];
rz(-0.71590483) q[3];
sx q[3];
rz(0.16815312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1049261) q[2];
sx q[2];
rz(-1.6380402) q[2];
sx q[2];
rz(-0.46910134) q[2];
rz(0.30096287) q[3];
sx q[3];
rz(-2.8807237) q[3];
sx q[3];
rz(1.4891967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24333532) q[0];
sx q[0];
rz(-1.5168334) q[0];
sx q[0];
rz(-2.6014056) q[0];
rz(-2.68014) q[1];
sx q[1];
rz(-1.5216454) q[1];
sx q[1];
rz(-2.8823421) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2592372) q[0];
sx q[0];
rz(-1.8412388) q[0];
sx q[0];
rz(-0.32201249) q[0];
rz(0.44633003) q[2];
sx q[2];
rz(-1.1023554) q[2];
sx q[2];
rz(0.51821741) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.67422489) q[1];
sx q[1];
rz(-1.7265872) q[1];
sx q[1];
rz(0.99512659) q[1];
rz(-pi) q[2];
rz(-1.2397887) q[3];
sx q[3];
rz(-1.8192623) q[3];
sx q[3];
rz(2.1295761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3429012) q[2];
sx q[2];
rz(-0.68621245) q[2];
sx q[2];
rz(-0.9101103) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27123505) q[0];
sx q[0];
rz(-2.7483181) q[0];
sx q[0];
rz(-1.9140592) q[0];
rz(-0.50499376) q[1];
sx q[1];
rz(-1.0161437) q[1];
sx q[1];
rz(1.203677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1211091) q[0];
sx q[0];
rz(-1.6899037) q[0];
sx q[0];
rz(0.24863338) q[0];
x q[1];
rz(0.9616379) q[2];
sx q[2];
rz(-0.74660245) q[2];
sx q[2];
rz(-0.37140977) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54654361) q[1];
sx q[1];
rz(-2.1814934) q[1];
sx q[1];
rz(-2.2799019) q[1];
x q[2];
rz(2.0304544) q[3];
sx q[3];
rz(-1.8098003) q[3];
sx q[3];
rz(-0.29005946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0248854) q[2];
sx q[2];
rz(-2.1163157) q[2];
sx q[2];
rz(2.1825979) q[2];
rz(0.2995019) q[3];
sx q[3];
rz(-3.0193269) q[3];
sx q[3];
rz(1.6925156) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9129979) q[0];
sx q[0];
rz(-1.2499502) q[0];
sx q[0];
rz(-0.49384299) q[0];
rz(2.2834868) q[1];
sx q[1];
rz(-1.9789507) q[1];
sx q[1];
rz(-1.7869305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49789594) q[0];
sx q[0];
rz(-2.5312811) q[0];
sx q[0];
rz(1.9998788) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0628711) q[2];
sx q[2];
rz(-1.6920231) q[2];
sx q[2];
rz(1.4722784) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1272116) q[1];
sx q[1];
rz(-0.76103044) q[1];
sx q[1];
rz(-3.0396023) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89833791) q[3];
sx q[3];
rz(-1.1500937) q[3];
sx q[3];
rz(1.4843693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6845503) q[2];
sx q[2];
rz(-0.6251308) q[2];
sx q[2];
rz(2.5954512) q[2];
rz(-2.9438733) q[3];
sx q[3];
rz(-2.1106014) q[3];
sx q[3];
rz(0.24699591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9207183) q[0];
sx q[0];
rz(-1.9863167) q[0];
sx q[0];
rz(2.0350631) q[0];
rz(0.58440343) q[1];
sx q[1];
rz(-1.1738651) q[1];
sx q[1];
rz(-1.0027142) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4179392) q[0];
sx q[0];
rz(-2.0040211) q[0];
sx q[0];
rz(0.17670011) q[0];
rz(-pi) q[1];
rz(-1.2757589) q[2];
sx q[2];
rz(-1.8031839) q[2];
sx q[2];
rz(0.7600998) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3467873) q[1];
sx q[1];
rz(-2.0827939) q[1];
sx q[1];
rz(1.5293299) q[1];
rz(-pi) q[2];
rz(-0.47559019) q[3];
sx q[3];
rz(-2.4659143) q[3];
sx q[3];
rz(1.6437363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.57548412) q[2];
sx q[2];
rz(-0.9245975) q[2];
sx q[2];
rz(1.7522579) q[2];
rz(-0.53650457) q[3];
sx q[3];
rz(-1.6335231) q[3];
sx q[3];
rz(0.84684053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.093512) q[0];
sx q[0];
rz(-2.7219613) q[0];
sx q[0];
rz(-2.804948) q[0];
rz(0.078314217) q[1];
sx q[1];
rz(-1.2093733) q[1];
sx q[1];
rz(-1.9858817) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9243568) q[0];
sx q[0];
rz(-1.2343533) q[0];
sx q[0];
rz(2.2948059) q[0];
x q[1];
rz(-0.066566932) q[2];
sx q[2];
rz(-1.9835501) q[2];
sx q[2];
rz(-2.33605) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0129688) q[1];
sx q[1];
rz(-1.4358724) q[1];
sx q[1];
rz(2.0634445) q[1];
rz(-2.8248057) q[3];
sx q[3];
rz(-1.5419349) q[3];
sx q[3];
rz(2.3742418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.081850514) q[2];
sx q[2];
rz(-1.0826702) q[2];
sx q[2];
rz(0.82320631) q[2];
rz(1.0169704) q[3];
sx q[3];
rz(-0.40139324) q[3];
sx q[3];
rz(0.038173525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3391649) q[0];
sx q[0];
rz(-0.22295727) q[0];
sx q[0];
rz(1.4270576) q[0];
rz(1.4786221) q[1];
sx q[1];
rz(-2.069811) q[1];
sx q[1];
rz(-0.85085416) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6991147) q[0];
sx q[0];
rz(-2.6529387) q[0];
sx q[0];
rz(1.8173993) q[0];
rz(-pi) q[1];
rz(3.0548373) q[2];
sx q[2];
rz(-1.6977001) q[2];
sx q[2];
rz(0.50228861) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5580299) q[1];
sx q[1];
rz(-1.0475912) q[1];
sx q[1];
rz(-1.6324522) q[1];
rz(-0.14031841) q[3];
sx q[3];
rz(-1.5068973) q[3];
sx q[3];
rz(-1.4177223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9783258) q[0];
sx q[0];
rz(-1.6303202) q[0];
sx q[0];
rz(-1.6250961) q[0];
rz(2.3691879) q[1];
sx q[1];
rz(-2.518387) q[1];
sx q[1];
rz(1.0060681) q[1];
rz(-1.7987235) q[2];
sx q[2];
rz(-0.89222903) q[2];
sx q[2];
rz(-1.0747838) q[2];
rz(1.1019273) q[3];
sx q[3];
rz(-1.1411323) q[3];
sx q[3];
rz(3.1235486) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
