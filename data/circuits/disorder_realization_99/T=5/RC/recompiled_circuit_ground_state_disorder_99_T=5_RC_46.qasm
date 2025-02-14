OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91166624) q[0];
sx q[0];
rz(-0.32719964) q[0];
sx q[0];
rz(-1.8157475) q[0];
rz(-1.7757379) q[1];
sx q[1];
rz(-2.7471625) q[1];
sx q[1];
rz(2.3529513) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74371007) q[0];
sx q[0];
rz(-1.8787483) q[0];
sx q[0];
rz(-0.90713199) q[0];
x q[1];
rz(2.4633258) q[2];
sx q[2];
rz(-2.1951198) q[2];
sx q[2];
rz(2.2140142) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4439075) q[1];
sx q[1];
rz(-0.34516066) q[1];
sx q[1];
rz(-2.5408695) q[1];
rz(-1.7472505) q[3];
sx q[3];
rz(-2.5061228) q[3];
sx q[3];
rz(0.42144767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.49393645) q[2];
sx q[2];
rz(-1.6948573) q[2];
sx q[2];
rz(-0.20230618) q[2];
rz(3.0322187) q[3];
sx q[3];
rz(-2.5325363) q[3];
sx q[3];
rz(-0.085453184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0789455) q[0];
sx q[0];
rz(-0.92342347) q[0];
sx q[0];
rz(-1.1751291) q[0];
rz(-1.6773978) q[1];
sx q[1];
rz(-1.0025832) q[1];
sx q[1];
rz(2.7820803) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4660366) q[0];
sx q[0];
rz(-0.17854843) q[0];
sx q[0];
rz(-1.9738865) q[0];
rz(-1.9486675) q[2];
sx q[2];
rz(-1.0456351) q[2];
sx q[2];
rz(1.1532056) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4103315) q[1];
sx q[1];
rz(-0.71440694) q[1];
sx q[1];
rz(-0.53602119) q[1];
rz(-2.5595846) q[3];
sx q[3];
rz(-2.0872384) q[3];
sx q[3];
rz(-0.88508115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29490617) q[2];
sx q[2];
rz(-1.0470231) q[2];
sx q[2];
rz(0.27493757) q[2];
rz(-1.8168195) q[3];
sx q[3];
rz(-1.6144269) q[3];
sx q[3];
rz(-1.8435318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0777968) q[0];
sx q[0];
rz(-2.9834788) q[0];
sx q[0];
rz(2.7445444) q[0];
rz(2.5322757) q[1];
sx q[1];
rz(-1.3343697) q[1];
sx q[1];
rz(-1.5416001) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3491495) q[0];
sx q[0];
rz(-2.8314548) q[0];
sx q[0];
rz(-1.9196426) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6360388) q[2];
sx q[2];
rz(-1.5833304) q[2];
sx q[2];
rz(-1.6103075) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1101602) q[1];
sx q[1];
rz(-2.2972887) q[1];
sx q[1];
rz(-3.1189819) q[1];
rz(-2.479781) q[3];
sx q[3];
rz(-0.99259171) q[3];
sx q[3];
rz(2.6933262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93727532) q[2];
sx q[2];
rz(-1.5466362) q[2];
sx q[2];
rz(-2.9079962) q[2];
rz(1.2967845) q[3];
sx q[3];
rz(-1.9498884) q[3];
sx q[3];
rz(-0.72062033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5066756) q[0];
sx q[0];
rz(-1.9838061) q[0];
sx q[0];
rz(1.7764212) q[0];
rz(-1.2719951) q[1];
sx q[1];
rz(-1.1993661) q[1];
sx q[1];
rz(-1.2619527) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.303971) q[0];
sx q[0];
rz(-0.6084992) q[0];
sx q[0];
rz(-2.801218) q[0];
rz(-pi) q[1];
rz(1.8795525) q[2];
sx q[2];
rz(-1.595904) q[2];
sx q[2];
rz(2.0659735) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.32324821) q[1];
sx q[1];
rz(-1.3882611) q[1];
sx q[1];
rz(1.6365461) q[1];
rz(-0.82437317) q[3];
sx q[3];
rz(-1.715987) q[3];
sx q[3];
rz(-0.33156013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23123732) q[2];
sx q[2];
rz(-1.0301901) q[2];
sx q[2];
rz(-0.78872284) q[2];
rz(-1.2809058) q[3];
sx q[3];
rz(-1.7773881) q[3];
sx q[3];
rz(1.0219215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(1.7541517) q[0];
sx q[0];
rz(-1.0197637) q[0];
sx q[0];
rz(-2.4805241) q[0];
rz(1.1152274) q[1];
sx q[1];
rz(-1.8293569) q[1];
sx q[1];
rz(0.95975319) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1486014) q[0];
sx q[0];
rz(-2.0856983) q[0];
sx q[0];
rz(1.2599676) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8371253) q[2];
sx q[2];
rz(-2.5779471) q[2];
sx q[2];
rz(-3.0802022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5670523) q[1];
sx q[1];
rz(-2.8216719) q[1];
sx q[1];
rz(-2.0053276) q[1];
rz(0.92613756) q[3];
sx q[3];
rz(-1.7297267) q[3];
sx q[3];
rz(1.8862806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58012086) q[2];
sx q[2];
rz(-0.46130195) q[2];
sx q[2];
rz(-2.0816154) q[2];
rz(2.3914242) q[3];
sx q[3];
rz(-1.7629905) q[3];
sx q[3];
rz(1.9157971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2724514) q[0];
sx q[0];
rz(-1.5584109) q[0];
sx q[0];
rz(-0.22431746) q[0];
rz(1.6250826) q[1];
sx q[1];
rz(-2.5254011) q[1];
sx q[1];
rz(0.075627653) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13091732) q[0];
sx q[0];
rz(-2.11158) q[0];
sx q[0];
rz(2.6877139) q[0];
x q[1];
rz(-0.096616726) q[2];
sx q[2];
rz(-1.268309) q[2];
sx q[2];
rz(-2.5563478) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5420679) q[1];
sx q[1];
rz(-2.0318527) q[1];
sx q[1];
rz(2.6510524) q[1];
rz(-0.84556112) q[3];
sx q[3];
rz(-1.2367833) q[3];
sx q[3];
rz(0.090868252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.533796) q[2];
sx q[2];
rz(-2.3888612) q[2];
sx q[2];
rz(0.047018615) q[2];
rz(0.67695824) q[3];
sx q[3];
rz(-1.8432063) q[3];
sx q[3];
rz(0.025378749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15636477) q[0];
sx q[0];
rz(-1.4893091) q[0];
sx q[0];
rz(-1.4439247) q[0];
rz(0.9230744) q[1];
sx q[1];
rz(-2.1363027) q[1];
sx q[1];
rz(2.9583171) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3175439) q[0];
sx q[0];
rz(-0.62301659) q[0];
sx q[0];
rz(-0.43806324) q[0];
rz(1.9940111) q[2];
sx q[2];
rz(-1.1626889) q[2];
sx q[2];
rz(-1.0427047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0959151) q[1];
sx q[1];
rz(-2.5325477) q[1];
sx q[1];
rz(-0.97859971) q[1];
x q[2];
rz(-2.1499726) q[3];
sx q[3];
rz(-1.8133111) q[3];
sx q[3];
rz(2.8057867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.034059374) q[2];
sx q[2];
rz(-1.4564161) q[2];
sx q[2];
rz(-3.122186) q[2];
rz(-0.16128811) q[3];
sx q[3];
rz(-2.1262157) q[3];
sx q[3];
rz(-1.9136782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1041383) q[0];
sx q[0];
rz(-2.7327974) q[0];
sx q[0];
rz(-3.0492875) q[0];
rz(1.1760938) q[1];
sx q[1];
rz(-2.8832925) q[1];
sx q[1];
rz(-1.4091122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21299674) q[0];
sx q[0];
rz(-0.11830506) q[0];
sx q[0];
rz(-2.0483093) q[0];
x q[1];
rz(0.94495602) q[2];
sx q[2];
rz(-0.71364489) q[2];
sx q[2];
rz(1.8533486) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84687877) q[1];
sx q[1];
rz(-2.6637609) q[1];
sx q[1];
rz(-0.61378693) q[1];
rz(0.73958379) q[3];
sx q[3];
rz(-0.7639262) q[3];
sx q[3];
rz(-0.52810625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1000503) q[2];
sx q[2];
rz(-1.5771733) q[2];
sx q[2];
rz(-1.1735631) q[2];
rz(-2.7581577) q[3];
sx q[3];
rz(-1.8337367) q[3];
sx q[3];
rz(1.0907382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1293056) q[0];
sx q[0];
rz(-1.4947991) q[0];
sx q[0];
rz(-2.9360085) q[0];
rz(2.3221817) q[1];
sx q[1];
rz(-1.2148379) q[1];
sx q[1];
rz(0.49088556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48614472) q[0];
sx q[0];
rz(-1.9441588) q[0];
sx q[0];
rz(2.4400737) q[0];
rz(2.728168) q[2];
sx q[2];
rz(-0.67081645) q[2];
sx q[2];
rz(2.6476423) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.017990503) q[1];
sx q[1];
rz(-1.6185056) q[1];
sx q[1];
rz(1.0022267) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.396046) q[3];
sx q[3];
rz(-2.6239003) q[3];
sx q[3];
rz(-0.92624901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4129591) q[2];
sx q[2];
rz(-2.1042991) q[2];
sx q[2];
rz(-1.884985) q[2];
rz(2.7058153) q[3];
sx q[3];
rz(-2.0332917) q[3];
sx q[3];
rz(0.88424879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4202145) q[0];
sx q[0];
rz(-0.90737897) q[0];
sx q[0];
rz(-0.23289982) q[0];
rz(0.13889343) q[1];
sx q[1];
rz(-1.1537617) q[1];
sx q[1];
rz(0.5084261) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.229436) q[0];
sx q[0];
rz(-2.5340853) q[0];
sx q[0];
rz(2.0338533) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0632238) q[2];
sx q[2];
rz(-1.3295577) q[2];
sx q[2];
rz(0.2080179) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0838937) q[1];
sx q[1];
rz(-3.0073799) q[1];
sx q[1];
rz(-1.5115518) q[1];
rz(2.2670161) q[3];
sx q[3];
rz(-1.1406058) q[3];
sx q[3];
rz(-3.0759574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4109219) q[2];
sx q[2];
rz(-1.6795009) q[2];
sx q[2];
rz(0.81653583) q[2];
rz(2.7517892) q[3];
sx q[3];
rz(-2.1392418) q[3];
sx q[3];
rz(1.3780814) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.447406) q[0];
sx q[0];
rz(-2.2185855) q[0];
sx q[0];
rz(2.9119281) q[0];
rz(1.5028839) q[1];
sx q[1];
rz(-1.8062183) q[1];
sx q[1];
rz(0.91581215) q[1];
rz(-2.7583078) q[2];
sx q[2];
rz(-1.7310168) q[2];
sx q[2];
rz(1.0723266) q[2];
rz(-2.2727454) q[3];
sx q[3];
rz(-1.6985536) q[3];
sx q[3];
rz(2.539195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
