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
rz(-1.0919548) q[0];
sx q[0];
rz(-1.0821992) q[0];
sx q[0];
rz(-1.9424633) q[0];
rz(0.24425976) q[1];
sx q[1];
rz(1.8241939) q[1];
sx q[1];
rz(7.1526935) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0713214) q[0];
sx q[0];
rz(-2.3247221) q[0];
sx q[0];
rz(-0.29036291) q[0];
rz(-1.0788953) q[2];
sx q[2];
rz(-1.190335) q[2];
sx q[2];
rz(1.5306533) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8972733) q[1];
sx q[1];
rz(-0.88855511) q[1];
sx q[1];
rz(1.2644493) q[1];
x q[2];
rz(-2.364089) q[3];
sx q[3];
rz(-0.99626675) q[3];
sx q[3];
rz(-3.0878909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96183744) q[2];
sx q[2];
rz(-1.1575674) q[2];
sx q[2];
rz(1.8400787) q[2];
rz(2.2573722) q[3];
sx q[3];
rz(-1.0853465) q[3];
sx q[3];
rz(1.6940544) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41534153) q[0];
sx q[0];
rz(-1.911442) q[0];
sx q[0];
rz(3.0905261) q[0];
rz(2.6452433) q[1];
sx q[1];
rz(-1.8615078) q[1];
sx q[1];
rz(3.0275717) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89897777) q[0];
sx q[0];
rz(-1.0568084) q[0];
sx q[0];
rz(-1.8144426) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.013234303) q[2];
sx q[2];
rz(-2.2639788) q[2];
sx q[2];
rz(2.7342094) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1920701) q[1];
sx q[1];
rz(-1.1382427) q[1];
sx q[1];
rz(-1.0640952) q[1];
x q[2];
rz(0.9468803) q[3];
sx q[3];
rz(-1.0495473) q[3];
sx q[3];
rz(0.60714591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30782917) q[2];
sx q[2];
rz(-1.2276063) q[2];
sx q[2];
rz(-2.7755136) q[2];
rz(2.6858373) q[3];
sx q[3];
rz(-2.3594806) q[3];
sx q[3];
rz(-0.2987878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5251821) q[0];
sx q[0];
rz(-2.1822378) q[0];
sx q[0];
rz(-0.33732238) q[0];
rz(1.4704544) q[1];
sx q[1];
rz(-0.93177876) q[1];
sx q[1];
rz(0.09387389) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81525345) q[0];
sx q[0];
rz(-1.4960644) q[0];
sx q[0];
rz(2.4228404) q[0];
rz(-pi) q[1];
rz(-2.502021) q[2];
sx q[2];
rz(-0.54023114) q[2];
sx q[2];
rz(-2.8566) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.962576) q[1];
sx q[1];
rz(-1.1454795) q[1];
sx q[1];
rz(3.1326181) q[1];
rz(-1.2971506) q[3];
sx q[3];
rz(-0.85262596) q[3];
sx q[3];
rz(2.9349365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1845392) q[2];
sx q[2];
rz(-2.5856555) q[2];
sx q[2];
rz(-0.24587336) q[2];
rz(-0.18694123) q[3];
sx q[3];
rz(-1.4176466) q[3];
sx q[3];
rz(0.42417446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251855) q[0];
sx q[0];
rz(-0.67916361) q[0];
sx q[0];
rz(-2.9492522) q[0];
rz(-1.964407) q[1];
sx q[1];
rz(-2.3277551) q[1];
sx q[1];
rz(-2.7331533) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.724736) q[0];
sx q[0];
rz(-1.724986) q[0];
sx q[0];
rz(1.6282999) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2871509) q[2];
sx q[2];
rz(-0.39098323) q[2];
sx q[2];
rz(1.6765082) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.08800297) q[1];
sx q[1];
rz(-1.9710959) q[1];
sx q[1];
rz(2.4884239) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3173265) q[3];
sx q[3];
rz(-2.4400024) q[3];
sx q[3];
rz(-1.4647558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.283215) q[2];
sx q[2];
rz(-0.89906162) q[2];
sx q[2];
rz(-0.176972) q[2];
rz(2.5519798) q[3];
sx q[3];
rz(-1.4859345) q[3];
sx q[3];
rz(-2.1069215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.706834) q[0];
sx q[0];
rz(-2.28573) q[0];
sx q[0];
rz(-0.51704299) q[0];
rz(2.928858) q[1];
sx q[1];
rz(-1.0849625) q[1];
sx q[1];
rz(-2.3092666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2281003) q[0];
sx q[0];
rz(-1.5287762) q[0];
sx q[0];
rz(1.5629004) q[0];
rz(-pi) q[1];
rz(2.2465785) q[2];
sx q[2];
rz(-2.7594406) q[2];
sx q[2];
rz(-2.7016751) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4246108) q[1];
sx q[1];
rz(-1.2094133) q[1];
sx q[1];
rz(-0.34504621) q[1];
x q[2];
rz(-1.4067235) q[3];
sx q[3];
rz(-2.0436156) q[3];
sx q[3];
rz(-2.7834322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0011255) q[2];
sx q[2];
rz(-2.3095864) q[2];
sx q[2];
rz(-2.7531085) q[2];
rz(-1.8899941) q[3];
sx q[3];
rz(-1.783071) q[3];
sx q[3];
rz(2.0775011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0851704) q[0];
sx q[0];
rz(-0.52260411) q[0];
sx q[0];
rz(1.1473468) q[0];
rz(-1.3832248) q[1];
sx q[1];
rz(-1.5937832) q[1];
sx q[1];
rz(0.17471084) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4736814) q[0];
sx q[0];
rz(-1.5725878) q[0];
sx q[0];
rz(-2.4666193) q[0];
rz(-1.5632929) q[2];
sx q[2];
rz(-0.96133366) q[2];
sx q[2];
rz(-2.1799019) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6683674) q[1];
sx q[1];
rz(-2.7764634) q[1];
sx q[1];
rz(-0.58546807) q[1];
rz(1.6047603) q[3];
sx q[3];
rz(-1.782825) q[3];
sx q[3];
rz(0.20454031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9354348) q[2];
sx q[2];
rz(-1.9400699) q[2];
sx q[2];
rz(-2.5344482) q[2];
rz(-0.28572765) q[3];
sx q[3];
rz(-2.5305179) q[3];
sx q[3];
rz(0.40209517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5173335) q[0];
sx q[0];
rz(-1.7835971) q[0];
sx q[0];
rz(0.72878033) q[0];
rz(1.1392611) q[1];
sx q[1];
rz(-2.0323205) q[1];
sx q[1];
rz(1.2713825) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8761638) q[0];
sx q[0];
rz(-3.1294397) q[0];
sx q[0];
rz(1.7677714) q[0];
rz(-pi) q[1];
rz(1.8461349) q[2];
sx q[2];
rz(-1.0022114) q[2];
sx q[2];
rz(-2.8953573) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67290689) q[1];
sx q[1];
rz(-1.7727643) q[1];
sx q[1];
rz(0.57909261) q[1];
rz(2.0470601) q[3];
sx q[3];
rz(-1.4753398) q[3];
sx q[3];
rz(0.80194471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3108958) q[2];
sx q[2];
rz(-1.7252012) q[2];
sx q[2];
rz(-0.26965109) q[2];
rz(2.0693178) q[3];
sx q[3];
rz(-2.8282073) q[3];
sx q[3];
rz(2.0872033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6894182) q[0];
sx q[0];
rz(-1.1913238) q[0];
sx q[0];
rz(-2.7574975) q[0];
rz(2.0317888) q[1];
sx q[1];
rz(-2.2566819) q[1];
sx q[1];
rz(-1.6222515) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9948331) q[0];
sx q[0];
rz(-1.681856) q[0];
sx q[0];
rz(1.5572118) q[0];
rz(1.0949959) q[2];
sx q[2];
rz(-0.738022) q[2];
sx q[2];
rz(1.2890146) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1422524) q[1];
sx q[1];
rz(-2.9007962) q[1];
sx q[1];
rz(-1.1204512) q[1];
rz(-0.31838454) q[3];
sx q[3];
rz(-0.51162135) q[3];
sx q[3];
rz(-1.4871396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9263837) q[2];
sx q[2];
rz(-1.0101725) q[2];
sx q[2];
rz(2.9877648) q[2];
rz(0.93872968) q[3];
sx q[3];
rz(-1.2490844) q[3];
sx q[3];
rz(1.4332829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.4140103) q[0];
sx q[0];
rz(-1.753597) q[0];
sx q[0];
rz(-2.599732) q[0];
rz(-1.9210723) q[1];
sx q[1];
rz(-0.6747171) q[1];
sx q[1];
rz(-0.064124785) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42196938) q[0];
sx q[0];
rz(-3.0753284) q[0];
sx q[0];
rz(0.55304773) q[0];
rz(-0.86875963) q[2];
sx q[2];
rz(-1.6693309) q[2];
sx q[2];
rz(1.2159525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2623973) q[1];
sx q[1];
rz(-0.91934312) q[1];
sx q[1];
rz(-1.5495954) q[1];
rz(-pi) q[2];
rz(-0.4145741) q[3];
sx q[3];
rz(-1.7547973) q[3];
sx q[3];
rz(2.6340226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8378143) q[2];
sx q[2];
rz(-2.1393675) q[2];
sx q[2];
rz(0.079806002) q[2];
rz(-0.59085733) q[3];
sx q[3];
rz(-0.30954027) q[3];
sx q[3];
rz(-3.0965366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87126842) q[0];
sx q[0];
rz(-2.7239983) q[0];
sx q[0];
rz(-2.8840892) q[0];
rz(-0.90266699) q[1];
sx q[1];
rz(-2.3011484) q[1];
sx q[1];
rz(2.2534175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1290263) q[0];
sx q[0];
rz(-1.9309224) q[0];
sx q[0];
rz(-8/(11*pi)) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3357847) q[2];
sx q[2];
rz(-2.1099612) q[2];
sx q[2];
rz(-1.6279814) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1969704) q[1];
sx q[1];
rz(-2.053481) q[1];
sx q[1];
rz(-0.75782366) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4367797) q[3];
sx q[3];
rz(-2.146943) q[3];
sx q[3];
rz(1.5004053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4070134) q[2];
sx q[2];
rz(-2.1740156) q[2];
sx q[2];
rz(0.61545294) q[2];
rz(2.9493799) q[3];
sx q[3];
rz(-0.8513611) q[3];
sx q[3];
rz(1.6063469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2007582) q[0];
sx q[0];
rz(-1.8466908) q[0];
sx q[0];
rz(1.2373663) q[0];
rz(1.0279961) q[1];
sx q[1];
rz(-0.68872394) q[1];
sx q[1];
rz(0.56180305) q[1];
rz(3.0801283) q[2];
sx q[2];
rz(-2.2066084) q[2];
sx q[2];
rz(1.6359272) q[2];
rz(-1.2067704) q[3];
sx q[3];
rz(-1.0586998) q[3];
sx q[3];
rz(2.8431526) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
