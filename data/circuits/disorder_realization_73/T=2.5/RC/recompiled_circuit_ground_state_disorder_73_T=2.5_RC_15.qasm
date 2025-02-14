OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.29326987) q[0];
sx q[0];
rz(3.3942437) q[0];
sx q[0];
rz(10.630339) q[0];
rz(2.0134917) q[1];
sx q[1];
rz(-1.3265346) q[1];
sx q[1];
rz(0.74572745) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4843895) q[0];
sx q[0];
rz(-1.6307388) q[0];
sx q[0];
rz(-1.7933229) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9171997) q[2];
sx q[2];
rz(-1.577415) q[2];
sx q[2];
rz(1.4057019) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2066951) q[1];
sx q[1];
rz(-1.3612011) q[1];
sx q[1];
rz(-1.9137726) q[1];
x q[2];
rz(2.0449355) q[3];
sx q[3];
rz(-1.385194) q[3];
sx q[3];
rz(-1.4193648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2241609) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(1.6708299) q[2];
rz(1.6338232) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(-2.2182218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5690174) q[0];
sx q[0];
rz(-1.9439531) q[0];
sx q[0];
rz(1.5731328) q[0];
rz(-2.9720427) q[1];
sx q[1];
rz(-3.0270271) q[1];
sx q[1];
rz(-0.13470185) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7260774) q[0];
sx q[0];
rz(-1.9163791) q[0];
sx q[0];
rz(-2.6879361) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6394272) q[2];
sx q[2];
rz(-1.5830212) q[2];
sx q[2];
rz(-3.1032012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33011064) q[1];
sx q[1];
rz(-0.54381424) q[1];
sx q[1];
rz(0.61119975) q[1];
rz(-pi) q[2];
x q[2];
rz(0.080218519) q[3];
sx q[3];
rz(-1.3933946) q[3];
sx q[3];
rz(-0.77840786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0160825) q[2];
sx q[2];
rz(-1.4849911) q[2];
sx q[2];
rz(2.9888195) q[2];
rz(-1.3656535) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(2.9735145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.633054) q[0];
sx q[0];
rz(-0.80798739) q[0];
sx q[0];
rz(0.48164865) q[0];
rz(-2.956849) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(-0.99536037) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99670519) q[0];
sx q[0];
rz(-1.2249682) q[0];
sx q[0];
rz(0.24390999) q[0];
rz(3.0929933) q[2];
sx q[2];
rz(-1.6311797) q[2];
sx q[2];
rz(-0.27602613) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0479483) q[1];
sx q[1];
rz(-1.8694832) q[1];
sx q[1];
rz(-0.68409749) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97484346) q[3];
sx q[3];
rz(-1.5961093) q[3];
sx q[3];
rz(1.5060177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2153726) q[2];
sx q[2];
rz(-0.054457713) q[2];
sx q[2];
rz(-3.0769297) q[2];
rz(-2.0926545) q[3];
sx q[3];
rz(-3.1147396) q[3];
sx q[3];
rz(-1.8612727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.8365086) q[0];
sx q[0];
rz(-0.11798141) q[0];
sx q[0];
rz(2.3205561) q[0];
rz(-3.0631284) q[1];
sx q[1];
rz(-1.6946225) q[1];
sx q[1];
rz(0.99748126) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0745463) q[0];
sx q[0];
rz(-0.95854488) q[0];
sx q[0];
rz(-1.1825829) q[0];
rz(-pi) q[1];
rz(1.6299963) q[2];
sx q[2];
rz(-1.6069876) q[2];
sx q[2];
rz(2.9179887) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23826829) q[1];
sx q[1];
rz(-1.0449755) q[1];
sx q[1];
rz(0.30491288) q[1];
x q[2];
rz(1.6281582) q[3];
sx q[3];
rz(-1.093075) q[3];
sx q[3];
rz(0.78335947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0923826) q[2];
sx q[2];
rz(-3.1123078) q[2];
sx q[2];
rz(1.1537665) q[2];
rz(-2.9554101) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(0.33128273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.004772923) q[0];
sx q[0];
rz(-0.81757075) q[0];
sx q[0];
rz(1.1432884) q[0];
rz(-2.000287) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(-0.51796651) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1521027) q[0];
sx q[0];
rz(-2.1955262) q[0];
sx q[0];
rz(-2.2476907) q[0];
x q[1];
rz(1.5557655) q[2];
sx q[2];
rz(-1.5717197) q[2];
sx q[2];
rz(-0.24713384) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6828047) q[1];
sx q[1];
rz(-2.7066351) q[1];
sx q[1];
rz(-0.87587237) q[1];
rz(-pi) q[2];
rz(-2.84359) q[3];
sx q[3];
rz(-1.7116705) q[3];
sx q[3];
rz(-0.44671392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(0.32773584) q[2];
rz(1.1945126) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(2.2849042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1389393) q[0];
sx q[0];
rz(-0.26873538) q[0];
sx q[0];
rz(0.56185454) q[0];
rz(-1.6506763) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(0.098310016) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54360169) q[0];
sx q[0];
rz(-2.7181135) q[0];
sx q[0];
rz(1.6571568) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.570669) q[2];
sx q[2];
rz(-1.5701558) q[2];
sx q[2];
rz(-2.8252779) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6209539) q[1];
sx q[1];
rz(-2.6700182) q[1];
sx q[1];
rz(-0.10271272) q[1];
x q[2];
rz(-0.50054153) q[3];
sx q[3];
rz(-1.9841474) q[3];
sx q[3];
rz(2.7374637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0564698) q[2];
sx q[2];
rz(-0.19530185) q[2];
sx q[2];
rz(1.0657715) q[2];
rz(-2.7417475) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(1.8736418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14168508) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(-1.6837233) q[0];
rz(2.0211925) q[1];
sx q[1];
rz(-3.0034062) q[1];
sx q[1];
rz(-0.33946005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0908112) q[0];
sx q[0];
rz(-3.0635186) q[0];
sx q[0];
rz(0.17103057) q[0];
x q[1];
rz(-1.5846662) q[2];
sx q[2];
rz(-1.0161581) q[2];
sx q[2];
rz(-0.0019794606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.554018) q[1];
sx q[1];
rz(-1.6563935) q[1];
sx q[1];
rz(0.0016335131) q[1];
rz(2.8651627) q[3];
sx q[3];
rz(-0.35110858) q[3];
sx q[3];
rz(1.029226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3642984) q[2];
sx q[2];
rz(-0.0019145049) q[2];
sx q[2];
rz(0.36336362) q[2];
rz(-2.0511138) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(1.1183848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5217487) q[0];
sx q[0];
rz(-0.32829568) q[0];
sx q[0];
rz(2.2186665) q[0];
rz(1.482831) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(0.0042075687) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.286233) q[0];
sx q[0];
rz(-0.84001937) q[0];
sx q[0];
rz(-0.71536163) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9816859) q[2];
sx q[2];
rz(-1.5762323) q[2];
sx q[2];
rz(1.5856575) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1394704) q[1];
sx q[1];
rz(-1.0957901) q[1];
sx q[1];
rz(0.072474555) q[1];
x q[2];
rz(0.065807314) q[3];
sx q[3];
rz(-2.5281457) q[3];
sx q[3];
rz(-0.25019161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5620455) q[2];
sx q[2];
rz(-1.5374708) q[2];
sx q[2];
rz(1.9407678) q[2];
rz(1.7480525) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(2.4414731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30018184) q[0];
sx q[0];
rz(-0.44898471) q[0];
sx q[0];
rz(-1.2196983) q[0];
rz(1.3297184) q[1];
sx q[1];
rz(-1.1390353) q[1];
sx q[1];
rz(2.9776998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9326646) q[0];
sx q[0];
rz(-0.39723662) q[0];
sx q[0];
rz(-2.0276638) q[0];
rz(-pi) q[1];
rz(-2.5570325) q[2];
sx q[2];
rz(-1.8001846) q[2];
sx q[2];
rz(0.10525119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1023952) q[1];
sx q[1];
rz(-1.6886687) q[1];
sx q[1];
rz(3.0815691) q[1];
x q[2];
rz(0.90198646) q[3];
sx q[3];
rz(-2.1928582) q[3];
sx q[3];
rz(-2.7354776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4280052) q[2];
sx q[2];
rz(-0.089944936) q[2];
sx q[2];
rz(-0.62291992) q[2];
rz(-0.04341393) q[3];
sx q[3];
rz(-0.89689887) q[3];
sx q[3];
rz(-0.77891946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.49122214) q[0];
sx q[0];
rz(-0.0064042052) q[0];
sx q[0];
rz(2.6553335) q[0];
rz(2.4468415) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(2.6659226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93906389) q[0];
sx q[0];
rz(-0.049204218) q[0];
sx q[0];
rz(-2.3243107) q[0];
rz(-pi) q[1];
x q[1];
rz(0.034157201) q[2];
sx q[2];
rz(-2.2766621) q[2];
sx q[2];
rz(-1.5547084) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.116058) q[1];
sx q[1];
rz(-2.5339911) q[1];
sx q[1];
rz(-3.1302384) q[1];
rz(-pi) q[2];
rz(0.52880295) q[3];
sx q[3];
rz(-0.35548726) q[3];
sx q[3];
rz(-2.5606009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67705578) q[2];
sx q[2];
rz(-0.060364351) q[2];
sx q[2];
rz(-1.8439058) q[2];
rz(-1.8929947) q[3];
sx q[3];
rz(-0.55515754) q[3];
sx q[3];
rz(-2.7738074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1260592) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(-2.2592648) q[1];
sx q[1];
rz(-3.075141) q[1];
sx q[1];
rz(-2.3022423) q[1];
rz(-2.7693314) q[2];
sx q[2];
rz(-3.0426171) q[2];
sx q[2];
rz(-0.097886861) q[2];
rz(1.5757379) q[3];
sx q[3];
rz(-1.3316122) q[3];
sx q[3];
rz(0.022461654) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
