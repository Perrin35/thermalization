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
rz(1.2441128) q[0];
sx q[0];
rz(3.4611995) q[0];
sx q[0];
rz(9.925579) q[0];
rz(-2.8473941) q[1];
sx q[1];
rz(-0.74217141) q[1];
sx q[1];
rz(2.5941526) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476737) q[0];
sx q[0];
rz(-1.8800056) q[0];
sx q[0];
rz(1.7237612) q[0];
x q[1];
rz(0.4644248) q[2];
sx q[2];
rz(-1.8193442) q[2];
sx q[2];
rz(-2.3119702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96066957) q[1];
sx q[1];
rz(-1.4804739) q[1];
sx q[1];
rz(2.6336843) q[1];
x q[2];
rz(0.35695255) q[3];
sx q[3];
rz(-2.0805177) q[3];
sx q[3];
rz(1.3689643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9886292) q[2];
sx q[2];
rz(-2.2965778) q[2];
sx q[2];
rz(-0.60626924) q[2];
rz(-1.2074977) q[3];
sx q[3];
rz(-1.8469801) q[3];
sx q[3];
rz(-1.5906364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2400874) q[0];
sx q[0];
rz(-1.5272239) q[0];
sx q[0];
rz(0.69183451) q[0];
rz(1.2884619) q[1];
sx q[1];
rz(-1.9970147) q[1];
sx q[1];
rz(0.28786927) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8765204) q[0];
sx q[0];
rz(-1.2848043) q[0];
sx q[0];
rz(-3.0080481) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1580196) q[2];
sx q[2];
rz(-2.1701239) q[2];
sx q[2];
rz(2.9119208) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46241847) q[1];
sx q[1];
rz(-0.44826213) q[1];
sx q[1];
rz(-3.03469) q[1];
rz(1.4631264) q[3];
sx q[3];
rz(-1.7070908) q[3];
sx q[3];
rz(1.8267269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3466779) q[2];
sx q[2];
rz(-1.6592462) q[2];
sx q[2];
rz(-2.1686926) q[2];
rz(1.3341303) q[3];
sx q[3];
rz(-0.79136807) q[3];
sx q[3];
rz(2.5202675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.36488229) q[0];
sx q[0];
rz(-2.8948687) q[0];
sx q[0];
rz(-0.2488939) q[0];
rz(1.6913255) q[1];
sx q[1];
rz(-1.1509117) q[1];
sx q[1];
rz(2.2379564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2118857) q[0];
sx q[0];
rz(-2.531966) q[0];
sx q[0];
rz(-3.097405) q[0];
x q[1];
rz(-1.9157639) q[2];
sx q[2];
rz(-1.6716701) q[2];
sx q[2];
rz(0.28988923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98164046) q[1];
sx q[1];
rz(-0.86152309) q[1];
sx q[1];
rz(2.0877128) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1670985) q[3];
sx q[3];
rz(-2.1362602) q[3];
sx q[3];
rz(1.3270449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16477188) q[2];
sx q[2];
rz(-0.73126078) q[2];
sx q[2];
rz(2.8110647) q[2];
rz(0.22322379) q[3];
sx q[3];
rz(-1.1130755) q[3];
sx q[3];
rz(-2.7631675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44499236) q[0];
sx q[0];
rz(-1.7855423) q[0];
sx q[0];
rz(-0.14183216) q[0];
rz(-3.0524571) q[1];
sx q[1];
rz(-0.24239692) q[1];
sx q[1];
rz(2.1884122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0702111) q[0];
sx q[0];
rz(-1.1888904) q[0];
sx q[0];
rz(0.75526636) q[0];
rz(-pi) q[1];
rz(-0.29229477) q[2];
sx q[2];
rz(-0.39543786) q[2];
sx q[2];
rz(-1.2751477) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.61342) q[1];
sx q[1];
rz(-2.0866359) q[1];
sx q[1];
rz(0.9524166) q[1];
x q[2];
rz(-1.370055) q[3];
sx q[3];
rz(-1.0827409) q[3];
sx q[3];
rz(1.0112011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7715093) q[2];
sx q[2];
rz(-0.13088317) q[2];
sx q[2];
rz(0.5292325) q[2];
rz(2.0970763) q[3];
sx q[3];
rz(-1.7525201) q[3];
sx q[3];
rz(2.9862826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20045497) q[0];
sx q[0];
rz(-1.2417685) q[0];
sx q[0];
rz(0.45503765) q[0];
rz(0.014668839) q[1];
sx q[1];
rz(-0.55611062) q[1];
sx q[1];
rz(-2.1591149) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1896613) q[0];
sx q[0];
rz(-1.2711382) q[0];
sx q[0];
rz(-1.4656187) q[0];
rz(2.234794) q[2];
sx q[2];
rz(-2.325146) q[2];
sx q[2];
rz(2.0891669) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8019218) q[1];
sx q[1];
rz(-1.4661487) q[1];
sx q[1];
rz(-0.27709477) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69070821) q[3];
sx q[3];
rz(-2.7793607) q[3];
sx q[3];
rz(-1.2129403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.602953) q[2];
sx q[2];
rz(-2.2621138) q[2];
sx q[2];
rz(-2.5243916) q[2];
rz(1.7547102) q[3];
sx q[3];
rz(-0.91104561) q[3];
sx q[3];
rz(2.9673747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0571601) q[0];
sx q[0];
rz(-0.76157695) q[0];
sx q[0];
rz(2.1630951) q[0];
rz(2.1242566) q[1];
sx q[1];
rz(-0.2778191) q[1];
sx q[1];
rz(0.4496347) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.250424) q[0];
sx q[0];
rz(-0.40968597) q[0];
sx q[0];
rz(-0.72575672) q[0];
rz(-2.6891031) q[2];
sx q[2];
rz(-0.86039484) q[2];
sx q[2];
rz(1.2373287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3842135) q[1];
sx q[1];
rz(-1.5259201) q[1];
sx q[1];
rz(-1.0991389) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71549039) q[3];
sx q[3];
rz(-2.3187197) q[3];
sx q[3];
rz(0.70586849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.52266878) q[2];
sx q[2];
rz(-1.7974682) q[2];
sx q[2];
rz(1.7632923) q[2];
rz(1.9887234) q[3];
sx q[3];
rz(-1.1533573) q[3];
sx q[3];
rz(2.8809179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5744837) q[0];
sx q[0];
rz(-2.4686047) q[0];
sx q[0];
rz(2.8048977) q[0];
rz(-0.7252655) q[1];
sx q[1];
rz(-1.5870353) q[1];
sx q[1];
rz(-2.7560962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8367675) q[0];
sx q[0];
rz(-2.4186717) q[0];
sx q[0];
rz(-0.54534955) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3213089) q[2];
sx q[2];
rz(-1.4324585) q[2];
sx q[2];
rz(-2.4398566) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32371569) q[1];
sx q[1];
rz(-2.1855832) q[1];
sx q[1];
rz(2.7842658) q[1];
rz(-0.5637873) q[3];
sx q[3];
rz(-0.70433378) q[3];
sx q[3];
rz(-0.35455656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5943299) q[2];
sx q[2];
rz(-1.9219834) q[2];
sx q[2];
rz(2.8426113) q[2];
rz(2.3494521) q[3];
sx q[3];
rz(-2.7463089) q[3];
sx q[3];
rz(-1.3778752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91656536) q[0];
sx q[0];
rz(-1.7074317) q[0];
sx q[0];
rz(0.7780956) q[0];
rz(-1.3968702) q[1];
sx q[1];
rz(-2.2010746) q[1];
sx q[1];
rz(3.0527589) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0087826) q[0];
sx q[0];
rz(-0.36160053) q[0];
sx q[0];
rz(2.8426991) q[0];
x q[1];
rz(-0.5702604) q[2];
sx q[2];
rz(-2.2024254) q[2];
sx q[2];
rz(-0.94550242) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0263097) q[1];
sx q[1];
rz(-2.6404318) q[1];
sx q[1];
rz(2.1351241) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2031112) q[3];
sx q[3];
rz(-1.5422149) q[3];
sx q[3];
rz(2.1052484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5953411) q[2];
sx q[2];
rz(-0.83117861) q[2];
sx q[2];
rz(-2.0383932) q[2];
rz(-2.5698419) q[3];
sx q[3];
rz(-0.56643707) q[3];
sx q[3];
rz(-1.5248689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53200805) q[0];
sx q[0];
rz(-0.12140618) q[0];
sx q[0];
rz(2.4914361) q[0];
rz(-2.0434168) q[1];
sx q[1];
rz(-1.8868586) q[1];
sx q[1];
rz(3.0756899) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9377334) q[0];
sx q[0];
rz(-2.2324076) q[0];
sx q[0];
rz(-0.96142049) q[0];
rz(-pi) q[1];
rz(0.19005929) q[2];
sx q[2];
rz(-1.5134317) q[2];
sx q[2];
rz(-0.8105958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84012953) q[1];
sx q[1];
rz(-0.80706978) q[1];
sx q[1];
rz(0.21424861) q[1];
rz(0.32089837) q[3];
sx q[3];
rz(-2.4545752) q[3];
sx q[3];
rz(-1.7248169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1419475) q[2];
sx q[2];
rz(-1.9530752) q[2];
sx q[2];
rz(-0.56335062) q[2];
rz(-0.79936409) q[3];
sx q[3];
rz(-2.1269709) q[3];
sx q[3];
rz(0.25580251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946741) q[0];
sx q[0];
rz(-0.020337157) q[0];
sx q[0];
rz(0.49816966) q[0];
rz(0.26214552) q[1];
sx q[1];
rz(-1.7356977) q[1];
sx q[1];
rz(-1.795305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7243626) q[0];
sx q[0];
rz(-0.87275332) q[0];
sx q[0];
rz(-2.5265187) q[0];
x q[1];
rz(-1.5602697) q[2];
sx q[2];
rz(-1.5684279) q[2];
sx q[2];
rz(-0.51711581) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0586068) q[1];
sx q[1];
rz(-0.78205651) q[1];
sx q[1];
rz(-2.6460827) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7896762) q[3];
sx q[3];
rz(-0.8026826) q[3];
sx q[3];
rz(-1.7520394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4761618) q[2];
sx q[2];
rz(-2.4527145) q[2];
sx q[2];
rz(-0.8821744) q[2];
rz(-2.4836508) q[3];
sx q[3];
rz(-2.0177757) q[3];
sx q[3];
rz(-2.1420124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5049725) q[0];
sx q[0];
rz(-1.5047147) q[0];
sx q[0];
rz(2.104105) q[0];
rz(0.50719117) q[1];
sx q[1];
rz(-0.7829983) q[1];
sx q[1];
rz(-1.0601039) q[1];
rz(-2.5835434) q[2];
sx q[2];
rz(-0.34172575) q[2];
sx q[2];
rz(1.5321129) q[2];
rz(-0.92978033) q[3];
sx q[3];
rz(-0.92757436) q[3];
sx q[3];
rz(-1.0074914) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
