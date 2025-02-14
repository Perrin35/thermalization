OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2158382) q[0];
sx q[0];
rz(-2.821142) q[0];
sx q[0];
rz(-3.0132063) q[0];
rz(-2.1865891) q[1];
sx q[1];
rz(-0.65989143) q[1];
sx q[1];
rz(-0.58712062) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40872657) q[0];
sx q[0];
rz(-0.48983296) q[0];
sx q[0];
rz(0.31087713) q[0];
rz(-pi) q[1];
rz(0.70694114) q[2];
sx q[2];
rz(-0.64919186) q[2];
sx q[2];
rz(-2.6892218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.93361357) q[1];
sx q[1];
rz(-0.99665239) q[1];
sx q[1];
rz(2.0962534) q[1];
x q[2];
rz(-0.0026851459) q[3];
sx q[3];
rz(-1.955282) q[3];
sx q[3];
rz(0.17385305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2653653) q[2];
sx q[2];
rz(-2.9076125) q[2];
sx q[2];
rz(-2.2205676) q[2];
rz(-1.6169351) q[3];
sx q[3];
rz(-1.8504668) q[3];
sx q[3];
rz(0.18572148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28144535) q[0];
sx q[0];
rz(-2.7466819) q[0];
sx q[0];
rz(1.7922147) q[0];
rz(2.7941864) q[1];
sx q[1];
rz(-1.9527304) q[1];
sx q[1];
rz(1.7385534) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31579298) q[0];
sx q[0];
rz(-0.99195882) q[0];
sx q[0];
rz(-1.3114503) q[0];
rz(-pi) q[1];
rz(2.1881869) q[2];
sx q[2];
rz(-2.2917903) q[2];
sx q[2];
rz(-1.7252201) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7561556) q[1];
sx q[1];
rz(-1.4175347) q[1];
sx q[1];
rz(0.71716829) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2590257) q[3];
sx q[3];
rz(-1.360713) q[3];
sx q[3];
rz(-0.20245598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8189524) q[2];
sx q[2];
rz(-1.3078657) q[2];
sx q[2];
rz(2.9827706) q[2];
rz(-2.9239376) q[3];
sx q[3];
rz(-1.7526151) q[3];
sx q[3];
rz(1.7910819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36667103) q[0];
sx q[0];
rz(-1.8407624) q[0];
sx q[0];
rz(1.8721254) q[0];
rz(0.41995755) q[1];
sx q[1];
rz(-2.1407514) q[1];
sx q[1];
rz(-1.5392083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8458662) q[0];
sx q[0];
rz(-0.55385607) q[0];
sx q[0];
rz(0.94828301) q[0];
rz(-pi) q[1];
rz(0.69373083) q[2];
sx q[2];
rz(-0.28273928) q[2];
sx q[2];
rz(-3.0764584) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0809787) q[1];
sx q[1];
rz(-1.1757869) q[1];
sx q[1];
rz(0.057821349) q[1];
rz(2.5747803) q[3];
sx q[3];
rz(-1.8517963) q[3];
sx q[3];
rz(2.5917201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0594844) q[2];
sx q[2];
rz(-2.2993645) q[2];
sx q[2];
rz(2.9498937) q[2];
rz(2.5890403) q[3];
sx q[3];
rz(-1.6165761) q[3];
sx q[3];
rz(-1.0244757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4100274) q[0];
sx q[0];
rz(-0.64842328) q[0];
sx q[0];
rz(2.5372274) q[0];
rz(-1.2454237) q[1];
sx q[1];
rz(-2.3912997) q[1];
sx q[1];
rz(2.2818458) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3533136) q[0];
sx q[0];
rz(-1.2333567) q[0];
sx q[0];
rz(2.2252625) q[0];
x q[1];
rz(1.4330788) q[2];
sx q[2];
rz(-3.0385512) q[2];
sx q[2];
rz(2.2400165) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27450505) q[1];
sx q[1];
rz(-2.9531859) q[1];
sx q[1];
rz(-2.8727358) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8243916) q[3];
sx q[3];
rz(-2.275064) q[3];
sx q[3];
rz(-0.21237838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8613646) q[2];
sx q[2];
rz(-2.1275529) q[2];
sx q[2];
rz(1.4385983) q[2];
rz(1.8493308) q[3];
sx q[3];
rz(-2.3137213) q[3];
sx q[3];
rz(-1.0546225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0786521) q[0];
sx q[0];
rz(-2.8684454) q[0];
sx q[0];
rz(0.31785059) q[0];
rz(-1.802313) q[1];
sx q[1];
rz(-2.7236718) q[1];
sx q[1];
rz(-0.97871614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1237549) q[0];
sx q[0];
rz(-2.7394501) q[0];
sx q[0];
rz(0.077362424) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42338223) q[2];
sx q[2];
rz(-2.9278594) q[2];
sx q[2];
rz(2.9345484) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3826344) q[1];
sx q[1];
rz(-2.3178604) q[1];
sx q[1];
rz(1.7428223) q[1];
rz(2.948165) q[3];
sx q[3];
rz(-1.3816091) q[3];
sx q[3];
rz(-2.6003169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9432482) q[2];
sx q[2];
rz(-2.1965616) q[2];
sx q[2];
rz(2.0242958) q[2];
rz(0.18103389) q[3];
sx q[3];
rz(-1.2369316) q[3];
sx q[3];
rz(1.9443289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4365874) q[0];
sx q[0];
rz(-0.26708189) q[0];
sx q[0];
rz(1.3096814) q[0];
rz(2.1197223) q[1];
sx q[1];
rz(-0.836687) q[1];
sx q[1];
rz(1.4885363) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7393417) q[0];
sx q[0];
rz(-1.2381123) q[0];
sx q[0];
rz(2.2646623) q[0];
rz(-pi) q[1];
rz(1.3371633) q[2];
sx q[2];
rz(-0.90590817) q[2];
sx q[2];
rz(-0.49371749) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52681953) q[1];
sx q[1];
rz(-1.0788016) q[1];
sx q[1];
rz(-1.7326049) q[1];
rz(-0.11373489) q[3];
sx q[3];
rz(-1.0128504) q[3];
sx q[3];
rz(-0.21765366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66270193) q[2];
sx q[2];
rz(-1.3871437) q[2];
sx q[2];
rz(-0.26408163) q[2];
rz(-1.665202) q[3];
sx q[3];
rz(-2.4614406) q[3];
sx q[3];
rz(0.41980729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5298117) q[0];
sx q[0];
rz(-1.6407069) q[0];
sx q[0];
rz(1.06426) q[0];
rz(-3.0423959) q[1];
sx q[1];
rz(-1.3490889) q[1];
sx q[1];
rz(1.4605716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6457466) q[0];
sx q[0];
rz(-0.48862132) q[0];
sx q[0];
rz(0.13617985) q[0];
rz(1.6978092) q[2];
sx q[2];
rz(-0.30108967) q[2];
sx q[2];
rz(1.3236698) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.34014186) q[1];
sx q[1];
rz(-0.51348084) q[1];
sx q[1];
rz(-2.5380847) q[1];
rz(-2.0123294) q[3];
sx q[3];
rz(-2.0534424) q[3];
sx q[3];
rz(-3.1114374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3508241) q[2];
sx q[2];
rz(-1.1188353) q[2];
sx q[2];
rz(-2.5420945) q[2];
rz(2.2439469) q[3];
sx q[3];
rz(-1.41956) q[3];
sx q[3];
rz(0.038014855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3277603) q[0];
sx q[0];
rz(-1.0603511) q[0];
sx q[0];
rz(1.7401485) q[0];
rz(3.0701045) q[1];
sx q[1];
rz(-1.46336) q[1];
sx q[1];
rz(-2.1404526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6342696) q[0];
sx q[0];
rz(-1.8232913) q[0];
sx q[0];
rz(2.6802313) q[0];
x q[1];
rz(-0.48076081) q[2];
sx q[2];
rz(-2.4188453) q[2];
sx q[2];
rz(1.9113845) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0948755) q[1];
sx q[1];
rz(-0.64697274) q[1];
sx q[1];
rz(0.21979522) q[1];
rz(-0.66949943) q[3];
sx q[3];
rz(-2.0560798) q[3];
sx q[3];
rz(0.2822596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1772168) q[2];
sx q[2];
rz(-1.588725) q[2];
sx q[2];
rz(-0.71844086) q[2];
rz(0.046772379) q[3];
sx q[3];
rz(-0.7074357) q[3];
sx q[3];
rz(2.5717946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30963323) q[0];
sx q[0];
rz(-0.15441144) q[0];
sx q[0];
rz(-0.69044789) q[0];
rz(0.18028232) q[1];
sx q[1];
rz(-2.0803662) q[1];
sx q[1];
rz(-0.32275018) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13234102) q[0];
sx q[0];
rz(-1.4173495) q[0];
sx q[0];
rz(1.6770393) q[0];
rz(-1.5602925) q[2];
sx q[2];
rz(-0.5521419) q[2];
sx q[2];
rz(0.10607468) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1806644) q[1];
sx q[1];
rz(-2.3276648) q[1];
sx q[1];
rz(-0.77427534) q[1];
rz(-2.2261103) q[3];
sx q[3];
rz(-2.142066) q[3];
sx q[3];
rz(2.9397428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6705769) q[2];
sx q[2];
rz(-1.7936423) q[2];
sx q[2];
rz(1.0178817) q[2];
rz(0.0035303591) q[3];
sx q[3];
rz(-2.1719833) q[3];
sx q[3];
rz(1.2118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3280535) q[0];
sx q[0];
rz(-0.73754755) q[0];
sx q[0];
rz(-0.92593431) q[0];
rz(3.0042341) q[1];
sx q[1];
rz(-2.4691212) q[1];
sx q[1];
rz(0.59991178) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023742763) q[0];
sx q[0];
rz(-1.5735007) q[0];
sx q[0];
rz(2.5976624) q[0];
rz(-pi) q[1];
rz(-0.62905797) q[2];
sx q[2];
rz(-1.1060083) q[2];
sx q[2];
rz(-2.128922) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2221189) q[1];
sx q[1];
rz(-1.6012234) q[1];
sx q[1];
rz(-0.98201507) q[1];
rz(-pi) q[2];
rz(-2.0447015) q[3];
sx q[3];
rz(-1.2792493) q[3];
sx q[3];
rz(0.64563676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4124734) q[2];
sx q[2];
rz(-1.0189265) q[2];
sx q[2];
rz(-2.4998383) q[2];
rz(1.1563835) q[3];
sx q[3];
rz(-0.76120794) q[3];
sx q[3];
rz(1.523264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625576) q[0];
sx q[0];
rz(-2.2618444) q[0];
sx q[0];
rz(-2.3656144) q[0];
rz(-1.7017801) q[1];
sx q[1];
rz(-1.4714614) q[1];
sx q[1];
rz(0.068838483) q[1];
rz(0.25664888) q[2];
sx q[2];
rz(-1.3779852) q[2];
sx q[2];
rz(0.040711395) q[2];
rz(-0.49700789) q[3];
sx q[3];
rz(-1.0773226) q[3];
sx q[3];
rz(1.1763431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
