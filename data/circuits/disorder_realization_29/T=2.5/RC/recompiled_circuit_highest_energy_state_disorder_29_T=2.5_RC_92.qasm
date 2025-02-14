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
rz(1.2443378) q[0];
sx q[0];
rz(2.3490348) q[0];
sx q[0];
rz(9.6777182) q[0];
rz(2.3674372) q[1];
sx q[1];
rz(-2.7634662) q[1];
sx q[1];
rz(2.7546496) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9772676) q[0];
sx q[0];
rz(-1.2480134) q[0];
sx q[0];
rz(2.8708145) q[0];
x q[1];
rz(0.57780452) q[2];
sx q[2];
rz(-1.6434323) q[2];
sx q[2];
rz(-2.5355123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2110212) q[1];
sx q[1];
rz(-1.6503578) q[1];
sx q[1];
rz(-3.0573175) q[1];
rz(-pi) q[2];
rz(2.4693842) q[3];
sx q[3];
rz(-0.70939964) q[3];
sx q[3];
rz(-1.0649452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0953377) q[2];
sx q[2];
rz(-0.75901186) q[2];
sx q[2];
rz(1.7389899) q[2];
rz(1.7272353) q[3];
sx q[3];
rz(-0.71310133) q[3];
sx q[3];
rz(0.49899092) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6801179) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(2.4312191) q[0];
rz(1.0164227) q[1];
sx q[1];
rz(-1.2998394) q[1];
sx q[1];
rz(0.76510915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7920096) q[0];
sx q[0];
rz(-0.63711626) q[0];
sx q[0];
rz(0.55320112) q[0];
x q[1];
rz(2.4508453) q[2];
sx q[2];
rz(-1.0260858) q[2];
sx q[2];
rz(2.4360457) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4087832) q[1];
sx q[1];
rz(-1.8955232) q[1];
sx q[1];
rz(0.88692437) q[1];
x q[2];
rz(1.2364976) q[3];
sx q[3];
rz(-1.7124885) q[3];
sx q[3];
rz(-0.96135512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0244828) q[2];
sx q[2];
rz(-0.64579248) q[2];
sx q[2];
rz(0.70466858) q[2];
rz(-2.9674528) q[3];
sx q[3];
rz(-2.2691059) q[3];
sx q[3];
rz(-0.38159889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0891377) q[0];
sx q[0];
rz(-1.8280886) q[0];
sx q[0];
rz(-0.88743368) q[0];
rz(1.2499836) q[1];
sx q[1];
rz(-1.2108112) q[1];
sx q[1];
rz(-1.7807622) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3850573) q[0];
sx q[0];
rz(-2.4891315) q[0];
sx q[0];
rz(-0.65076179) q[0];
x q[1];
rz(1.046692) q[2];
sx q[2];
rz(-2.9827318) q[2];
sx q[2];
rz(2.9422974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20818921) q[1];
sx q[1];
rz(-1.4706685) q[1];
sx q[1];
rz(-0.80762819) q[1];
rz(-3.0988337) q[3];
sx q[3];
rz(-1.1819289) q[3];
sx q[3];
rz(-1.0145018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36797324) q[2];
sx q[2];
rz(-0.36586389) q[2];
sx q[2];
rz(-1.492929) q[2];
rz(0.44935539) q[3];
sx q[3];
rz(-1.8466693) q[3];
sx q[3];
rz(1.9213283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.087695) q[0];
sx q[0];
rz(-2.5515285) q[0];
sx q[0];
rz(1.7328523) q[0];
rz(2.159481) q[1];
sx q[1];
rz(-1.5928007) q[1];
sx q[1];
rz(-1.4062448) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2236693) q[0];
sx q[0];
rz(-1.6817998) q[0];
sx q[0];
rz(1.7864321) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6236308) q[2];
sx q[2];
rz(-1.6103193) q[2];
sx q[2];
rz(-1.7588577) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6342369) q[1];
sx q[1];
rz(-1.0546396) q[1];
sx q[1];
rz(-2.8576351) q[1];
rz(-pi) q[2];
rz(0.26867433) q[3];
sx q[3];
rz(-0.93963913) q[3];
sx q[3];
rz(2.1674066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9881607) q[2];
sx q[2];
rz(-1.016523) q[2];
sx q[2];
rz(-0.15596685) q[2];
rz(0.65034741) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(0.11421886) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0735737) q[0];
sx q[0];
rz(-1.8719712) q[0];
sx q[0];
rz(-0.20075783) q[0];
rz(-1.1772032) q[1];
sx q[1];
rz(-0.8526082) q[1];
sx q[1];
rz(-1.9815365) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2748799) q[0];
sx q[0];
rz(-1.794941) q[0];
sx q[0];
rz(-0.17205162) q[0];
rz(-0.54520901) q[2];
sx q[2];
rz(-0.13449796) q[2];
sx q[2];
rz(2.4328977) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0811942) q[1];
sx q[1];
rz(-1.2676928) q[1];
sx q[1];
rz(0.36784192) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9405064) q[3];
sx q[3];
rz(-2.0917836) q[3];
sx q[3];
rz(-2.9231284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7572299) q[2];
sx q[2];
rz(-2.195916) q[2];
sx q[2];
rz(2.6833656) q[2];
rz(2.5107757) q[3];
sx q[3];
rz(-1.9061371) q[3];
sx q[3];
rz(2.6423776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.683627) q[0];
sx q[0];
rz(-2.2891335) q[0];
sx q[0];
rz(-0.0068579554) q[0];
rz(-2.9099756) q[1];
sx q[1];
rz(-1.3791142) q[1];
sx q[1];
rz(1.0622271) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7201014) q[0];
sx q[0];
rz(-1.694223) q[0];
sx q[0];
rz(-2.8272259) q[0];
rz(-2.7659645) q[2];
sx q[2];
rz(-2.6115719) q[2];
sx q[2];
rz(-0.14800581) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66328996) q[1];
sx q[1];
rz(-1.3172825) q[1];
sx q[1];
rz(0.67914825) q[1];
rz(-pi) q[2];
rz(0.8616448) q[3];
sx q[3];
rz(-1.6966736) q[3];
sx q[3];
rz(1.5365841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14083938) q[2];
sx q[2];
rz(-1.2522298) q[2];
sx q[2];
rz(-1.8737277) q[2];
rz(2.2800692) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2086901) q[0];
sx q[0];
rz(-2.8080495) q[0];
sx q[0];
rz(2.9260337) q[0];
rz(-0.49628273) q[1];
sx q[1];
rz(-1.1880621) q[1];
sx q[1];
rz(0.15942474) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2886141) q[0];
sx q[0];
rz(-0.63487494) q[0];
sx q[0];
rz(-1.5663516) q[0];
x q[1];
rz(2.6706598) q[2];
sx q[2];
rz(-1.3976946) q[2];
sx q[2];
rz(1.3927695) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7466399) q[1];
sx q[1];
rz(-1.3356908) q[1];
sx q[1];
rz(3.0993153) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46298085) q[3];
sx q[3];
rz(-2.0714875) q[3];
sx q[3];
rz(0.41741308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7364007) q[2];
sx q[2];
rz(-1.1242194) q[2];
sx q[2];
rz(2.6250725) q[2];
rz(-1.2107595) q[3];
sx q[3];
rz(-1.4529994) q[3];
sx q[3];
rz(1.7621382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4844168) q[0];
sx q[0];
rz(-0.076198904) q[0];
sx q[0];
rz(2.6013689) q[0];
rz(1.5243439) q[1];
sx q[1];
rz(-2.1397619) q[1];
sx q[1];
rz(1.8744972) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57813533) q[0];
sx q[0];
rz(-1.7892013) q[0];
sx q[0];
rz(-0.51213034) q[0];
rz(-pi) q[1];
rz(-2.3643963) q[2];
sx q[2];
rz(-1.0269477) q[2];
sx q[2];
rz(1.0842619) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1649176) q[1];
sx q[1];
rz(-1.7150214) q[1];
sx q[1];
rz(-0.37101908) q[1];
rz(-0.35532002) q[3];
sx q[3];
rz(-2.1815724) q[3];
sx q[3];
rz(-0.54278096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8871062) q[2];
sx q[2];
rz(-1.9820513) q[2];
sx q[2];
rz(0.17733388) q[2];
rz(-1.0698498) q[3];
sx q[3];
rz(-1.0515352) q[3];
sx q[3];
rz(-0.18226084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808481) q[0];
sx q[0];
rz(-1.1540664) q[0];
sx q[0];
rz(-3.0408707) q[0];
rz(1.1489457) q[1];
sx q[1];
rz(-1.9457685) q[1];
sx q[1];
rz(1.1740059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3927395) q[0];
sx q[0];
rz(-2.1707188) q[0];
sx q[0];
rz(1.8003182) q[0];
x q[1];
rz(-3.0413925) q[2];
sx q[2];
rz(-0.76876193) q[2];
sx q[2];
rz(-2.5454846) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2795978) q[1];
sx q[1];
rz(-1.5902385) q[1];
sx q[1];
rz(1.6085546) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4467008) q[3];
sx q[3];
rz(-0.69487232) q[3];
sx q[3];
rz(1.2563039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.164087) q[2];
sx q[2];
rz(-1.4848494) q[2];
sx q[2];
rz(2.738415) q[2];
rz(2.4781503) q[3];
sx q[3];
rz(-0.57938975) q[3];
sx q[3];
rz(-3.1373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73096257) q[0];
sx q[0];
rz(-1.0799438) q[0];
sx q[0];
rz(-0.078911111) q[0];
rz(-0.90947378) q[1];
sx q[1];
rz(-1.0458922) q[1];
sx q[1];
rz(-0.39631072) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3672597) q[0];
sx q[0];
rz(-3.1235187) q[0];
sx q[0];
rz(-2.0967183) q[0];
x q[1];
rz(-0.79371039) q[2];
sx q[2];
rz(-1.5657445) q[2];
sx q[2];
rz(1.2271301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4034319) q[1];
sx q[1];
rz(-2.6034174) q[1];
sx q[1];
rz(-2.6713283) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86045281) q[3];
sx q[3];
rz(-0.75569587) q[3];
sx q[3];
rz(0.44482728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69184715) q[2];
sx q[2];
rz(-2.2578466) q[2];
sx q[2];
rz(0.79908243) q[2];
rz(0.07829047) q[3];
sx q[3];
rz(-1.2543863) q[3];
sx q[3];
rz(2.4962795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.66452022) q[0];
sx q[0];
rz(-1.7417396) q[0];
sx q[0];
rz(-2.59792) q[0];
rz(1.1935344) q[1];
sx q[1];
rz(-1.8652893) q[1];
sx q[1];
rz(0.51020772) q[1];
rz(1.2340679) q[2];
sx q[2];
rz(-2.6609145) q[2];
sx q[2];
rz(2.547472) q[2];
rz(-0.54936784) q[3];
sx q[3];
rz(-2.3372169) q[3];
sx q[3];
rz(-1.611471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
