OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3684664) q[0];
sx q[0];
rz(-2.3946895) q[0];
sx q[0];
rz(-0.84063831) q[0];
rz(-3.0199938) q[1];
sx q[1];
rz(-1.8687948) q[1];
sx q[1];
rz(2.8425541) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75343412) q[0];
sx q[0];
rz(-0.82601604) q[0];
sx q[0];
rz(2.1623934) q[0];
rz(1.7173355) q[2];
sx q[2];
rz(-0.49078178) q[2];
sx q[2];
rz(-2.7637568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4579276) q[1];
sx q[1];
rz(-1.8103765) q[1];
sx q[1];
rz(1.423382) q[1];
rz(-pi) q[2];
rz(0.14278966) q[3];
sx q[3];
rz(-1.1200818) q[3];
sx q[3];
rz(0.55270178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3806939) q[2];
sx q[2];
rz(-2.0913405) q[2];
sx q[2];
rz(-0.32763457) q[2];
rz(-1.7662175) q[3];
sx q[3];
rz(-1.7211434) q[3];
sx q[3];
rz(-2.6473141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7647917) q[0];
sx q[0];
rz(-1.6015653) q[0];
sx q[0];
rz(-0.85533992) q[0];
rz(-0.0080464706) q[1];
sx q[1];
rz(-1.2866373) q[1];
sx q[1];
rz(2.6904552) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6159003) q[0];
sx q[0];
rz(-2.2886758) q[0];
sx q[0];
rz(0.20882102) q[0];
rz(2.3694384) q[2];
sx q[2];
rz(-1.7092961) q[2];
sx q[2];
rz(-2.0318299) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0476183) q[1];
sx q[1];
rz(-2.6674358) q[1];
sx q[1];
rz(-0.38800254) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9782009) q[3];
sx q[3];
rz(-1.1278858) q[3];
sx q[3];
rz(0.16678424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1796639) q[2];
sx q[2];
rz(-2.0018061) q[2];
sx q[2];
rz(-3.0120604) q[2];
rz(0.1772964) q[3];
sx q[3];
rz(-2.5640021) q[3];
sx q[3];
rz(3.0887443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.391908) q[0];
sx q[0];
rz(-2.7295697) q[0];
sx q[0];
rz(2.5823197) q[0];
rz(3.0468805) q[1];
sx q[1];
rz(-1.5385224) q[1];
sx q[1];
rz(-0.47725484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4995296) q[0];
sx q[0];
rz(-1.1746049) q[0];
sx q[0];
rz(2.8642333) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7202708) q[2];
sx q[2];
rz(-1.9286619) q[2];
sx q[2];
rz(-1.1015111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6544158) q[1];
sx q[1];
rz(-2.3337939) q[1];
sx q[1];
rz(-2.828965) q[1];
x q[2];
rz(-0.52366728) q[3];
sx q[3];
rz(-0.91991495) q[3];
sx q[3];
rz(-1.6802579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39530784) q[2];
sx q[2];
rz(-1.2664653) q[2];
sx q[2];
rz(-0.9300119) q[2];
rz(1.2578472) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(3.0959082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78087085) q[0];
sx q[0];
rz(-1.5732795) q[0];
sx q[0];
rz(-2.1029396) q[0];
rz(2.5301798) q[1];
sx q[1];
rz(-0.88714209) q[1];
sx q[1];
rz(-2.2183653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8002565) q[0];
sx q[0];
rz(-2.1390599) q[0];
sx q[0];
rz(-2.7305014) q[0];
rz(-pi) q[1];
rz(1.5973041) q[2];
sx q[2];
rz(-0.4178646) q[2];
sx q[2];
rz(-2.1919427) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9070396) q[1];
sx q[1];
rz(-0.72826339) q[1];
sx q[1];
rz(-1.4437463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8037854) q[3];
sx q[3];
rz(-0.67854133) q[3];
sx q[3];
rz(-1.9935009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5103147) q[2];
sx q[2];
rz(-0.46773657) q[2];
sx q[2];
rz(-2.5941217) q[2];
rz(-0.064420961) q[3];
sx q[3];
rz(-1.8378601) q[3];
sx q[3];
rz(2.2678383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.9861458) q[0];
sx q[0];
rz(-2.2387945) q[0];
sx q[0];
rz(1.6814394) q[0];
rz(2.1414781) q[1];
sx q[1];
rz(-0.8668879) q[1];
sx q[1];
rz(0.11988457) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6836707) q[0];
sx q[0];
rz(-0.73212762) q[0];
sx q[0];
rz(1.2326272) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91953711) q[2];
sx q[2];
rz(-2.0555758) q[2];
sx q[2];
rz(-2.9491021) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.29741187) q[1];
sx q[1];
rz(-1.4919123) q[1];
sx q[1];
rz(-2.9260738) q[1];
x q[2];
rz(-2.5402903) q[3];
sx q[3];
rz(-1.9996694) q[3];
sx q[3];
rz(1.2473904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.039006058) q[2];
sx q[2];
rz(-0.90709364) q[2];
sx q[2];
rz(-2.8809123) q[2];
rz(-0.87812224) q[3];
sx q[3];
rz(-1.9120646) q[3];
sx q[3];
rz(-2.3584283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53073019) q[0];
sx q[0];
rz(-0.090228883) q[0];
sx q[0];
rz(-2.1395785) q[0];
rz(2.7643381) q[1];
sx q[1];
rz(-0.95163029) q[1];
sx q[1];
rz(2.196905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024600895) q[0];
sx q[0];
rz(-0.62378609) q[0];
sx q[0];
rz(2.7626415) q[0];
rz(-pi) q[1];
rz(2.5414921) q[2];
sx q[2];
rz(-2.1515111) q[2];
sx q[2];
rz(-2.3490259) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2608883) q[1];
sx q[1];
rz(-0.35683888) q[1];
sx q[1];
rz(0.44512213) q[1];
rz(-pi) q[2];
rz(-0.11774534) q[3];
sx q[3];
rz(-2.1097933) q[3];
sx q[3];
rz(0.6955516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79823309) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(-0.043924335) q[2];
rz(-1.515306) q[3];
sx q[3];
rz(-2.01912) q[3];
sx q[3];
rz(1.4656434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2783022) q[0];
sx q[0];
rz(-0.55178061) q[0];
sx q[0];
rz(0.58468753) q[0];
rz(-2.0901285) q[1];
sx q[1];
rz(-0.81948558) q[1];
sx q[1];
rz(-2.6712766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4620275) q[0];
sx q[0];
rz(-0.60067486) q[0];
sx q[0];
rz(0.26058773) q[0];
rz(-pi) q[1];
rz(0.83018556) q[2];
sx q[2];
rz(-2.6174394) q[2];
sx q[2];
rz(2.1955127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0329758) q[1];
sx q[1];
rz(-1.2101296) q[1];
sx q[1];
rz(-2.3920139) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.028585) q[3];
sx q[3];
rz(-1.5829493) q[3];
sx q[3];
rz(3.0719824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8768002) q[2];
sx q[2];
rz(-0.53692836) q[2];
sx q[2];
rz(0.31141591) q[2];
rz(0.16820678) q[3];
sx q[3];
rz(-1.5663389) q[3];
sx q[3];
rz(0.28347191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6300221) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(-0.93609634) q[0];
rz(0.49631897) q[1];
sx q[1];
rz(-2.4744108) q[1];
sx q[1];
rz(-0.47028968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054809991) q[0];
sx q[0];
rz(-2.5544689) q[0];
sx q[0];
rz(2.5820288) q[0];
x q[1];
rz(-1.7685212) q[2];
sx q[2];
rz(-1.573296) q[2];
sx q[2];
rz(-0.041415215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.072402231) q[1];
sx q[1];
rz(-1.6399929) q[1];
sx q[1];
rz(1.9131768) q[1];
rz(-pi) q[2];
rz(1.3948729) q[3];
sx q[3];
rz(-1.8430437) q[3];
sx q[3];
rz(1.7737089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63142598) q[2];
sx q[2];
rz(-1.7743856) q[2];
sx q[2];
rz(-1.6660956) q[2];
rz(-0.79536074) q[3];
sx q[3];
rz(-0.16203351) q[3];
sx q[3];
rz(1.2696666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0610166) q[0];
sx q[0];
rz(-0.59122714) q[0];
sx q[0];
rz(0.06037816) q[0];
rz(0.16054842) q[1];
sx q[1];
rz(-1.5269273) q[1];
sx q[1];
rz(-2.1626332) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6290603) q[0];
sx q[0];
rz(-2.1950025) q[0];
sx q[0];
rz(-1.4275622) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8674815) q[2];
sx q[2];
rz(-2.2345671) q[2];
sx q[2];
rz(1.889515) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.97040788) q[1];
sx q[1];
rz(-1.1029585) q[1];
sx q[1];
rz(-2.4572608) q[1];
rz(0.03421182) q[3];
sx q[3];
rz(-2.1497823) q[3];
sx q[3];
rz(1.8949231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0673361) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(-3.0687029) q[2];
rz(-0.59761754) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(-3.1387175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6947967) q[0];
sx q[0];
rz(-1.0121166) q[0];
sx q[0];
rz(0.52892518) q[0];
rz(2.9534598) q[1];
sx q[1];
rz(-0.70536047) q[1];
sx q[1];
rz(2.5573152) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0952632) q[0];
sx q[0];
rz(-1.8285311) q[0];
sx q[0];
rz(1.8247752) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4227117) q[2];
sx q[2];
rz(-2.5846057) q[2];
sx q[2];
rz(0.62918909) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6021298) q[1];
sx q[1];
rz(-1.5279084) q[1];
sx q[1];
rz(-0.40956386) q[1];
rz(0.95240611) q[3];
sx q[3];
rz(-1.3460396) q[3];
sx q[3];
rz(-2.145203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75954413) q[2];
sx q[2];
rz(-2.1376762) q[2];
sx q[2];
rz(3.0697401) q[2];
rz(2.1182649) q[3];
sx q[3];
rz(-1.5724678) q[3];
sx q[3];
rz(2.6527827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.95579424) q[0];
sx q[0];
rz(-2.4837942) q[0];
sx q[0];
rz(2.876045) q[0];
rz(0.67509782) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(3.0037389) q[2];
sx q[2];
rz(-2.6225435) q[2];
sx q[2];
rz(-0.23487716) q[2];
rz(1.731338) q[3];
sx q[3];
rz(-1.4057126) q[3];
sx q[3];
rz(-2.2688903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
