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
rz(-2.8219858) q[0];
sx q[0];
rz(0.50080103) q[0];
rz(-2.8473941) q[1];
sx q[1];
rz(-0.74217141) q[1];
sx q[1];
rz(-0.54744005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8590605) q[0];
sx q[0];
rz(-1.4251389) q[0];
sx q[0];
rz(-2.8289625) q[0];
rz(1.8473832) q[2];
sx q[2];
rz(-1.1217077) q[2];
sx q[2];
rz(-0.61855471) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3708122) q[1];
sx q[1];
rz(-0.51518422) q[1];
sx q[1];
rz(-2.9574802) q[1];
rz(-pi) q[2];
rz(0.35695255) q[3];
sx q[3];
rz(-2.0805177) q[3];
sx q[3];
rz(-1.7726284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9886292) q[2];
sx q[2];
rz(-0.84501481) q[2];
sx q[2];
rz(2.5353234) q[2];
rz(1.9340949) q[3];
sx q[3];
rz(-1.8469801) q[3];
sx q[3];
rz(1.5509563) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2400874) q[0];
sx q[0];
rz(-1.6143687) q[0];
sx q[0];
rz(2.4497581) q[0];
rz(1.2884619) q[1];
sx q[1];
rz(-1.9970147) q[1];
sx q[1];
rz(-2.8537234) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26507227) q[0];
sx q[0];
rz(-1.8567883) q[0];
sx q[0];
rz(-3.0080481) q[0];
rz(-2.4601654) q[2];
sx q[2];
rz(-2.3286901) q[2];
sx q[2];
rz(-0.63804783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1296245) q[1];
sx q[1];
rz(-1.6170562) q[1];
sx q[1];
rz(2.6955626) q[1];
rz(-pi) q[2];
rz(2.4769529) q[3];
sx q[3];
rz(-0.17348504) q[3];
sx q[3];
rz(-1.9868614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3466779) q[2];
sx q[2];
rz(-1.6592462) q[2];
sx q[2];
rz(2.1686926) q[2];
rz(-1.8074624) q[3];
sx q[3];
rz(-2.3502246) q[3];
sx q[3];
rz(-2.5202675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36488229) q[0];
sx q[0];
rz(-0.24672395) q[0];
sx q[0];
rz(2.8926988) q[0];
rz(-1.6913255) q[1];
sx q[1];
rz(-1.9906809) q[1];
sx q[1];
rz(2.2379564) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8758276) q[0];
sx q[0];
rz(-2.1797414) q[0];
sx q[0];
rz(1.5399571) q[0];
x q[1];
rz(-0.10714178) q[2];
sx q[2];
rz(-1.2276548) q[2];
sx q[2];
rz(1.8245152) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.98164046) q[1];
sx q[1];
rz(-2.2800696) q[1];
sx q[1];
rz(-2.0877128) q[1];
x q[2];
rz(1.9744942) q[3];
sx q[3];
rz(-2.1362602) q[3];
sx q[3];
rz(1.3270449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9768208) q[2];
sx q[2];
rz(-2.4103319) q[2];
sx q[2];
rz(0.33052793) q[2];
rz(-2.9183689) q[3];
sx q[3];
rz(-2.0285172) q[3];
sx q[3];
rz(2.7631675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44499236) q[0];
sx q[0];
rz(-1.3560504) q[0];
sx q[0];
rz(0.14183216) q[0];
rz(0.089135535) q[1];
sx q[1];
rz(-0.24239692) q[1];
sx q[1];
rz(2.1884122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16312604) q[0];
sx q[0];
rz(-0.88136601) q[0];
sx q[0];
rz(-1.0667144) q[0];
rz(-pi) q[1];
rz(0.3802657) q[2];
sx q[2];
rz(-1.6820246) q[2];
sx q[2];
rz(-2.5750772) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2948896) q[1];
sx q[1];
rz(-2.0993471) q[1];
sx q[1];
rz(-2.5336086) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.370055) q[3];
sx q[3];
rz(-1.0827409) q[3];
sx q[3];
rz(-2.1303915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7715093) q[2];
sx q[2];
rz(-3.0107095) q[2];
sx q[2];
rz(-0.5292325) q[2];
rz(2.0970763) q[3];
sx q[3];
rz(-1.7525201) q[3];
sx q[3];
rz(2.9862826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20045497) q[0];
sx q[0];
rz(-1.8998242) q[0];
sx q[0];
rz(-0.45503765) q[0];
rz(0.014668839) q[1];
sx q[1];
rz(-0.55611062) q[1];
sx q[1];
rz(0.98247772) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1896613) q[0];
sx q[0];
rz(-1.8704544) q[0];
sx q[0];
rz(1.675974) q[0];
rz(-0.90679866) q[2];
sx q[2];
rz(-0.81644662) q[2];
sx q[2];
rz(-2.0891669) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0207386) q[1];
sx q[1];
rz(-2.845872) q[1];
sx q[1];
rz(-0.3665847) q[1];
rz(1.807688) q[3];
sx q[3];
rz(-1.29414) q[3];
sx q[3];
rz(-1.2048126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.602953) q[2];
sx q[2];
rz(-0.87947881) q[2];
sx q[2];
rz(0.61720103) q[2];
rz(1.3868825) q[3];
sx q[3];
rz(-0.91104561) q[3];
sx q[3];
rz(0.174218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0571601) q[0];
sx q[0];
rz(-2.3800157) q[0];
sx q[0];
rz(2.1630951) q[0];
rz(-1.017336) q[1];
sx q[1];
rz(-0.2778191) q[1];
sx q[1];
rz(0.4496347) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7788197) q[0];
sx q[0];
rz(-1.8383433) q[0];
sx q[0];
rz(-2.8275201) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6891031) q[2];
sx q[2];
rz(-0.86039484) q[2];
sx q[2];
rz(1.2373287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3510531) q[1];
sx q[1];
rz(-1.0996523) q[1];
sx q[1];
rz(0.050367722) q[1];
rz(-0.95532598) q[3];
sx q[3];
rz(-2.1571473) q[3];
sx q[3];
rz(2.9406656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6189239) q[2];
sx q[2];
rz(-1.3441244) q[2];
sx q[2];
rz(-1.3783003) q[2];
rz(1.9887234) q[3];
sx q[3];
rz(-1.9882354) q[3];
sx q[3];
rz(0.26067477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5744837) q[0];
sx q[0];
rz(-2.4686047) q[0];
sx q[0];
rz(-2.8048977) q[0];
rz(-2.4163272) q[1];
sx q[1];
rz(-1.5545574) q[1];
sx q[1];
rz(0.38549647) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9850901) q[0];
sx q[0];
rz(-2.1719731) q[0];
sx q[0];
rz(1.1416092) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0573559) q[2];
sx q[2];
rz(-0.2845736) q[2];
sx q[2];
rz(-1.3651266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.817877) q[1];
sx q[1];
rz(-2.1855832) q[1];
sx q[1];
rz(2.7842658) q[1];
rz(1.9970421) q[3];
sx q[3];
rz(-2.1499471) q[3];
sx q[3];
rz(-0.338011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5943299) q[2];
sx q[2];
rz(-1.2196093) q[2];
sx q[2];
rz(2.8426113) q[2];
rz(0.79214054) q[3];
sx q[3];
rz(-0.39528379) q[3];
sx q[3];
rz(-1.3778752) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2250273) q[0];
sx q[0];
rz(-1.4341609) q[0];
sx q[0];
rz(-0.7780956) q[0];
rz(1.3968702) q[1];
sx q[1];
rz(-0.94051802) q[1];
sx q[1];
rz(3.0527589) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7185811) q[0];
sx q[0];
rz(-1.4664343) q[0];
sx q[0];
rz(-2.7947438) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5702604) q[2];
sx q[2];
rz(-2.2024254) q[2];
sx q[2];
rz(-2.1960902) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.962304) q[1];
sx q[1];
rz(-1.8306762) q[1];
sx q[1];
rz(-2.0043027) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.030627294) q[3];
sx q[3];
rz(-1.9383241) q[3];
sx q[3];
rz(-2.5961329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54625154) q[2];
sx q[2];
rz(-0.83117861) q[2];
sx q[2];
rz(-1.1031995) q[2];
rz(0.57175076) q[3];
sx q[3];
rz(-0.56643707) q[3];
sx q[3];
rz(1.6167238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53200805) q[0];
sx q[0];
rz(-3.0201865) q[0];
sx q[0];
rz(-2.4914361) q[0];
rz(1.0981759) q[1];
sx q[1];
rz(-1.8868586) q[1];
sx q[1];
rz(-0.065902725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1797721) q[0];
sx q[0];
rz(-2.0393436) q[0];
sx q[0];
rz(0.75956068) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6292105) q[2];
sx q[2];
rz(-1.7605392) q[2];
sx q[2];
rz(-0.74917114) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.84012953) q[1];
sx q[1];
rz(-2.3345229) q[1];
sx q[1];
rz(-0.21424861) q[1];
rz(-pi) q[2];
rz(-1.8239924) q[3];
sx q[3];
rz(-2.2166219) q[3];
sx q[3];
rz(1.822804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.99964511) q[2];
sx q[2];
rz(-1.1885175) q[2];
sx q[2];
rz(2.578242) q[2];
rz(-0.79936409) q[3];
sx q[3];
rz(-1.0146217) q[3];
sx q[3];
rz(-0.25580251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.946741) q[0];
sx q[0];
rz(-0.020337157) q[0];
sx q[0];
rz(-0.49816966) q[0];
rz(0.26214552) q[1];
sx q[1];
rz(-1.7356977) q[1];
sx q[1];
rz(-1.795305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58459613) q[0];
sx q[0];
rz(-0.89476953) q[0];
sx q[0];
rz(-0.96831324) q[0];
x q[1];
rz(1.5813229) q[2];
sx q[2];
rz(-1.5731647) q[2];
sx q[2];
rz(0.51711581) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0586068) q[1];
sx q[1];
rz(-0.78205651) q[1];
sx q[1];
rz(2.6460827) q[1];
rz(-2.9204921) q[3];
sx q[3];
rz(-2.3490864) q[3];
sx q[3];
rz(-1.6994051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4761618) q[2];
sx q[2];
rz(-0.68887812) q[2];
sx q[2];
rz(-0.8821744) q[2];
rz(2.4836508) q[3];
sx q[3];
rz(-1.123817) q[3];
sx q[3];
rz(0.9995802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.6366202) q[0];
sx q[0];
rz(-1.5047147) q[0];
sx q[0];
rz(2.104105) q[0];
rz(0.50719117) q[1];
sx q[1];
rz(-0.7829983) q[1];
sx q[1];
rz(-1.0601039) q[1];
rz(1.7569594) q[2];
sx q[2];
rz(-1.8590448) q[2];
sx q[2];
rz(0.94696905) q[2];
rz(0.92978033) q[3];
sx q[3];
rz(-2.2140183) q[3];
sx q[3];
rz(2.1341013) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
