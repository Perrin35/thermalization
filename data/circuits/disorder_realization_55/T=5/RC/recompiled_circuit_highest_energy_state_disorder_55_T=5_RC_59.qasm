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
rz(-3.0814085) q[0];
sx q[0];
rz(-1.0589851) q[0];
sx q[0];
rz(-2.1314148) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(-0.29616907) q[1];
sx q[1];
rz(-2.8316166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.298807) q[0];
sx q[0];
rz(-1.1514542) q[0];
sx q[0];
rz(1.4489277) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9363382) q[2];
sx q[2];
rz(-1.5927093) q[2];
sx q[2];
rz(-2.2147629) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12116005) q[1];
sx q[1];
rz(-1.766664) q[1];
sx q[1];
rz(-1.905026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3499851) q[3];
sx q[3];
rz(-2.3901) q[3];
sx q[3];
rz(1.3205075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6628722) q[2];
sx q[2];
rz(-2.8933849) q[2];
sx q[2];
rz(-0.09566801) q[2];
rz(-3.0293448) q[3];
sx q[3];
rz(-0.87533689) q[3];
sx q[3];
rz(1.4314502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2410759) q[0];
sx q[0];
rz(-0.88923419) q[0];
sx q[0];
rz(-2.526793) q[0];
rz(-0.88848937) q[1];
sx q[1];
rz(-1.588984) q[1];
sx q[1];
rz(-0.49957553) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1785806) q[0];
sx q[0];
rz(-2.7546429) q[0];
sx q[0];
rz(-2.1301756) q[0];
x q[1];
rz(1.7287298) q[2];
sx q[2];
rz(-1.0519895) q[2];
sx q[2];
rz(3.0268596) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.080729) q[1];
sx q[1];
rz(-1.398096) q[1];
sx q[1];
rz(-2.2279146) q[1];
x q[2];
rz(-2.0604443) q[3];
sx q[3];
rz(-2.5003365) q[3];
sx q[3];
rz(-1.7288417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8932314) q[2];
sx q[2];
rz(-1.9042559) q[2];
sx q[2];
rz(0.020817967) q[2];
rz(-1.2402395) q[3];
sx q[3];
rz(-1.6720684) q[3];
sx q[3];
rz(1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5027387) q[0];
sx q[0];
rz(-2.8732712) q[0];
sx q[0];
rz(-2.8079206) q[0];
rz(-2.4036713) q[1];
sx q[1];
rz(-1.9155904) q[1];
sx q[1];
rz(-3.0121682) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5564011) q[0];
sx q[0];
rz(-1.7298152) q[0];
sx q[0];
rz(-0.88944737) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1652075) q[2];
sx q[2];
rz(-0.98862851) q[2];
sx q[2];
rz(-0.79603031) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1947965) q[1];
sx q[1];
rz(-1.4490286) q[1];
sx q[1];
rz(3.0508079) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78227104) q[3];
sx q[3];
rz(-1.2201628) q[3];
sx q[3];
rz(-2.3334437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1239329) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(-1.370149) q[2];
rz(-0.74639368) q[3];
sx q[3];
rz(-1.3179702) q[3];
sx q[3];
rz(-3.0588176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0081886) q[0];
sx q[0];
rz(-1.9318102) q[0];
sx q[0];
rz(-2.5764537) q[0];
rz(2.7124229) q[1];
sx q[1];
rz(-2.4950835) q[1];
sx q[1];
rz(0.82829222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9681184) q[0];
sx q[0];
rz(-0.83195639) q[0];
sx q[0];
rz(1.4627187) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0617375) q[2];
sx q[2];
rz(-0.54414302) q[2];
sx q[2];
rz(2.6119815) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1431378) q[1];
sx q[1];
rz(-1.3052485) q[1];
sx q[1];
rz(2.6396855) q[1];
x q[2];
rz(-2.8935562) q[3];
sx q[3];
rz(-1.0674849) q[3];
sx q[3];
rz(2.4865884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4565178) q[2];
sx q[2];
rz(-2.0802616) q[2];
sx q[2];
rz(-1.7100517) q[2];
rz(-0.93959129) q[3];
sx q[3];
rz(-2.6288433) q[3];
sx q[3];
rz(0.00069869839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13570304) q[0];
sx q[0];
rz(-0.26517427) q[0];
sx q[0];
rz(1.8121207) q[0];
rz(-1.9895408) q[1];
sx q[1];
rz(-1.0877129) q[1];
sx q[1];
rz(-0.00024814127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2787001) q[0];
sx q[0];
rz(-0.30847806) q[0];
sx q[0];
rz(-2.2174979) q[0];
rz(-pi) q[1];
rz(1.7176487) q[2];
sx q[2];
rz(-0.53895437) q[2];
sx q[2];
rz(2.1951064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8793638) q[1];
sx q[1];
rz(-1.2254224) q[1];
sx q[1];
rz(2.097258) q[1];
x q[2];
rz(2.4100634) q[3];
sx q[3];
rz(-1.5168744) q[3];
sx q[3];
rz(1.1889072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1382711) q[2];
sx q[2];
rz(-1.2500117) q[2];
sx q[2];
rz(-0.0058343466) q[2];
rz(2.2516069) q[3];
sx q[3];
rz(-2.4599288) q[3];
sx q[3];
rz(-1.1267004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936845) q[0];
sx q[0];
rz(-1.6981145) q[0];
sx q[0];
rz(-1.4877315) q[0];
rz(3.1112025) q[1];
sx q[1];
rz(-1.1266212) q[1];
sx q[1];
rz(0.94246513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2104032) q[0];
sx q[0];
rz(-1.382834) q[0];
sx q[0];
rz(-2.8920449) q[0];
x q[1];
rz(-1.312227) q[2];
sx q[2];
rz(-1.2814008) q[2];
sx q[2];
rz(-2.0813775) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2487464) q[1];
sx q[1];
rz(-2.5809118) q[1];
sx q[1];
rz(2.1302724) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0464155) q[3];
sx q[3];
rz(-1.3381357) q[3];
sx q[3];
rz(1.7847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59948644) q[2];
sx q[2];
rz(-0.78080559) q[2];
sx q[2];
rz(1.8360809) q[2];
rz(-0.54715884) q[3];
sx q[3];
rz(-1.2055509) q[3];
sx q[3];
rz(0.80593306) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18043537) q[0];
sx q[0];
rz(-2.1077709) q[0];
sx q[0];
rz(-3.0410774) q[0];
rz(0.48249498) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(-0.85711342) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6109638) q[0];
sx q[0];
rz(-0.93994323) q[0];
sx q[0];
rz(-1.1683589) q[0];
x q[1];
rz(3.0761511) q[2];
sx q[2];
rz(-0.64787728) q[2];
sx q[2];
rz(2.6104948) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.313011) q[1];
sx q[1];
rz(-1.9691113) q[1];
sx q[1];
rz(-1.3446773) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1825652) q[3];
sx q[3];
rz(-2.534158) q[3];
sx q[3];
rz(1.5807815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7414005) q[2];
sx q[2];
rz(-0.97101784) q[2];
sx q[2];
rz(2.8391489) q[2];
rz(-0.97964573) q[3];
sx q[3];
rz(-1.0023578) q[3];
sx q[3];
rz(1.2790595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7804467) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(0.030990344) q[0];
rz(0.57394761) q[1];
sx q[1];
rz(-1.6684883) q[1];
sx q[1];
rz(-1.8642289) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.325186) q[0];
sx q[0];
rz(-1.2379097) q[0];
sx q[0];
rz(-1.216796) q[0];
rz(2.8988016) q[2];
sx q[2];
rz(-2.0567679) q[2];
sx q[2];
rz(-2.9832341) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9121418) q[1];
sx q[1];
rz(-1.0849285) q[1];
sx q[1];
rz(-3.058601) q[1];
rz(1.8559009) q[3];
sx q[3];
rz(-1.630097) q[3];
sx q[3];
rz(1.3467195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.10451) q[2];
sx q[2];
rz(-0.13585486) q[2];
sx q[2];
rz(2.6541397) q[2];
rz(-2.4449352) q[3];
sx q[3];
rz(-2.2127377) q[3];
sx q[3];
rz(-0.063974403) q[3];
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
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697407) q[0];
sx q[0];
rz(-2.6672279) q[0];
sx q[0];
rz(-0.35097861) q[0];
rz(0.17414302) q[1];
sx q[1];
rz(-1.5085647) q[1];
sx q[1];
rz(-1.7399656) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49920666) q[0];
sx q[0];
rz(-2.8642352) q[0];
sx q[0];
rz(-2.6794898) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0171595) q[2];
sx q[2];
rz(-2.2702262) q[2];
sx q[2];
rz(2.1289189) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8216711) q[1];
sx q[1];
rz(-0.18583365) q[1];
sx q[1];
rz(-1.4129078) q[1];
x q[2];
rz(-0.74123592) q[3];
sx q[3];
rz(-2.0307856) q[3];
sx q[3];
rz(0.408084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1104687) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(-2.5049211) q[2];
rz(2.0311671) q[3];
sx q[3];
rz(-1.597155) q[3];
sx q[3];
rz(-0.5184263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3670032) q[0];
sx q[0];
rz(-2.84802) q[0];
sx q[0];
rz(2.0830182) q[0];
rz(0.038657945) q[1];
sx q[1];
rz(-1.6148753) q[1];
sx q[1];
rz(-2.0775332) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2702613) q[0];
sx q[0];
rz(-0.0055905213) q[0];
sx q[0];
rz(-2.5020775) q[0];
x q[1];
rz(1.2821141) q[2];
sx q[2];
rz(-1.836986) q[2];
sx q[2];
rz(2.4634354) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1326931) q[1];
sx q[1];
rz(-1.6311967) q[1];
sx q[1];
rz(-1.5302883) q[1];
rz(-1.6056772) q[3];
sx q[3];
rz(-2.0300755) q[3];
sx q[3];
rz(-1.1920795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71558636) q[2];
sx q[2];
rz(-0.49804372) q[2];
sx q[2];
rz(3.0675724) q[2];
rz(0.38153875) q[3];
sx q[3];
rz(-1.3149202) q[3];
sx q[3];
rz(0.35416245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.030180177) q[0];
sx q[0];
rz(-1.8288061) q[0];
sx q[0];
rz(-0.70963138) q[0];
rz(-2.5297655) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(-2.5135573) q[2];
sx q[2];
rz(-0.59878329) q[2];
sx q[2];
rz(-1.5046635) q[2];
rz(-1.0985804) q[3];
sx q[3];
rz(-2.2618812) q[3];
sx q[3];
rz(-1.19899) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
