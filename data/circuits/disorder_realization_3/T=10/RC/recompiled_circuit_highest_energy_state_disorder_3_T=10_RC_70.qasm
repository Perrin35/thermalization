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
rz(-2.6920707) q[0];
sx q[0];
rz(-1.7188526) q[0];
sx q[0];
rz(2.0916405) q[0];
rz(-0.51896754) q[1];
sx q[1];
rz(-1.5331886) q[1];
sx q[1];
rz(-3.0906711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1919132) q[0];
sx q[0];
rz(-0.55242071) q[0];
sx q[0];
rz(-0.64903736) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21587431) q[2];
sx q[2];
rz(-1.2435438) q[2];
sx q[2];
rz(2.808771) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4733153) q[1];
sx q[1];
rz(-0.44679579) q[1];
sx q[1];
rz(2.3832641) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9031119) q[3];
sx q[3];
rz(-0.89934228) q[3];
sx q[3];
rz(0.32794288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.27429399) q[2];
sx q[2];
rz(-0.87624246) q[2];
sx q[2];
rz(-1.743861) q[2];
rz(2.6869669) q[3];
sx q[3];
rz(-0.8780829) q[3];
sx q[3];
rz(-1.7129869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.4947263) q[0];
sx q[0];
rz(-0.88748256) q[0];
sx q[0];
rz(0.60321641) q[0];
rz(-1.2546722) q[1];
sx q[1];
rz(-0.61379495) q[1];
sx q[1];
rz(0.57993728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80548426) q[0];
sx q[0];
rz(-1.0382129) q[0];
sx q[0];
rz(1.7642154) q[0];
rz(-pi) q[1];
rz(-2.0236778) q[2];
sx q[2];
rz(-1.0476026) q[2];
sx q[2];
rz(2.4459237) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7407826) q[1];
sx q[1];
rz(-2.1745958) q[1];
sx q[1];
rz(-0.048193805) q[1];
rz(-pi) q[2];
rz(-1.9469684) q[3];
sx q[3];
rz(-2.1056089) q[3];
sx q[3];
rz(-2.4589698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6569528) q[2];
sx q[2];
rz(-0.75289774) q[2];
sx q[2];
rz(0.37772712) q[2];
rz(1.4332625) q[3];
sx q[3];
rz(-1.5092311) q[3];
sx q[3];
rz(1.9506955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79320532) q[0];
sx q[0];
rz(-0.054940104) q[0];
sx q[0];
rz(1.4980263) q[0];
rz(1.161423) q[1];
sx q[1];
rz(-1.8893416) q[1];
sx q[1];
rz(-2.2983671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0204165) q[0];
sx q[0];
rz(-1.5112229) q[0];
sx q[0];
rz(-2.1942433) q[0];
x q[1];
rz(0.21768985) q[2];
sx q[2];
rz(-1.161631) q[2];
sx q[2];
rz(-0.37693757) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9595782) q[1];
sx q[1];
rz(-1.6977662) q[1];
sx q[1];
rz(-2.2499529) q[1];
x q[2];
rz(-2.7494861) q[3];
sx q[3];
rz(-0.74015731) q[3];
sx q[3];
rz(1.8405434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8680385) q[2];
sx q[2];
rz(-0.66323438) q[2];
sx q[2];
rz(-2.3466477) q[2];
rz(-2.3718209) q[3];
sx q[3];
rz(-2.218518) q[3];
sx q[3];
rz(0.15708378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5539826) q[0];
sx q[0];
rz(-1.8164604) q[0];
sx q[0];
rz(2.8532568) q[0];
rz(0.68710697) q[1];
sx q[1];
rz(-1.4930875) q[1];
sx q[1];
rz(1.4289325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7369722) q[0];
sx q[0];
rz(-1.5828743) q[0];
sx q[0];
rz(1.5642479) q[0];
rz(0.32055118) q[2];
sx q[2];
rz(-0.42141576) q[2];
sx q[2];
rz(-1.5326064) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7894788) q[1];
sx q[1];
rz(-1.6360456) q[1];
sx q[1];
rz(-2.8852374) q[1];
x q[2];
rz(2.3219548) q[3];
sx q[3];
rz(-2.313531) q[3];
sx q[3];
rz(2.6520906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5700506) q[2];
sx q[2];
rz(-2.1752581) q[2];
sx q[2];
rz(0.16560444) q[2];
rz(-3.0770732) q[3];
sx q[3];
rz(-0.12972984) q[3];
sx q[3];
rz(-1.8701514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9128543) q[0];
sx q[0];
rz(-1.9485291) q[0];
sx q[0];
rz(-0.91113973) q[0];
rz(0.97081026) q[1];
sx q[1];
rz(-1.3263005) q[1];
sx q[1];
rz(0.86311805) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1939094) q[0];
sx q[0];
rz(-1.2811617) q[0];
sx q[0];
rz(2.0783483) q[0];
rz(-1.4865033) q[2];
sx q[2];
rz(-2.8721923) q[2];
sx q[2];
rz(1.0625372) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.48820254) q[1];
sx q[1];
rz(-1.0428671) q[1];
sx q[1];
rz(1.8406244) q[1];
rz(-1.4310525) q[3];
sx q[3];
rz(-0.31523963) q[3];
sx q[3];
rz(1.8342575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0866278) q[2];
sx q[2];
rz(-1.8885771) q[2];
sx q[2];
rz(-2.6524554) q[2];
rz(2.1431811) q[3];
sx q[3];
rz(-2.9676134) q[3];
sx q[3];
rz(-3.0424931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84422207) q[0];
sx q[0];
rz(-2.4956644) q[0];
sx q[0];
rz(-2.5883664) q[0];
rz(-0.79752254) q[1];
sx q[1];
rz(-2.1452466) q[1];
sx q[1];
rz(-1.9901989) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11565514) q[0];
sx q[0];
rz(-2.135072) q[0];
sx q[0];
rz(-1.4141809) q[0];
rz(-pi) q[1];
rz(1.3998447) q[2];
sx q[2];
rz(-2.1151849) q[2];
sx q[2];
rz(2.483727) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4321255) q[1];
sx q[1];
rz(-1.3015987) q[1];
sx q[1];
rz(1.7749191) q[1];
x q[2];
rz(1.0955515) q[3];
sx q[3];
rz(-2.4184347) q[3];
sx q[3];
rz(2.1840546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79661757) q[2];
sx q[2];
rz(-1.4130219) q[2];
sx q[2];
rz(0.83149347) q[2];
rz(-1.8031395) q[3];
sx q[3];
rz(-2.0667388) q[3];
sx q[3];
rz(-0.89053806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7688585) q[0];
sx q[0];
rz(-1.2059728) q[0];
sx q[0];
rz(-0.15705577) q[0];
rz(1.8231237) q[1];
sx q[1];
rz(-0.942197) q[1];
sx q[1];
rz(0.35776055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6899932) q[0];
sx q[0];
rz(-2.4781961) q[0];
sx q[0];
rz(-0.76865743) q[0];
x q[1];
rz(0.90334185) q[2];
sx q[2];
rz(-1.9304262) q[2];
sx q[2];
rz(-1.0718759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.58531144) q[1];
sx q[1];
rz(-1.0204633) q[1];
sx q[1];
rz(-2.354391) q[1];
rz(-0.57447432) q[3];
sx q[3];
rz(-1.7041429) q[3];
sx q[3];
rz(2.3337618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1498108) q[2];
sx q[2];
rz(-1.2573743) q[2];
sx q[2];
rz(-0.020616654) q[2];
rz(-1.8990382) q[3];
sx q[3];
rz(-2.3239682) q[3];
sx q[3];
rz(0.29153618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.6108516) q[0];
sx q[0];
rz(-1.956097) q[0];
sx q[0];
rz(2.8148742) q[0];
rz(2.2945981) q[1];
sx q[1];
rz(-2.0920483) q[1];
sx q[1];
rz(1.0583896) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7493454) q[0];
sx q[0];
rz(-1.6405182) q[0];
sx q[0];
rz(2.0850943) q[0];
x q[1];
rz(-2.5024274) q[2];
sx q[2];
rz(-1.4286388) q[2];
sx q[2];
rz(-2.7246812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0363473) q[1];
sx q[1];
rz(-0.39064841) q[1];
sx q[1];
rz(1.3643144) q[1];
x q[2];
rz(0.31553573) q[3];
sx q[3];
rz(-1.0275176) q[3];
sx q[3];
rz(-1.4062987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1191001) q[2];
sx q[2];
rz(-2.7969226) q[2];
sx q[2];
rz(1.9692839) q[2];
rz(-0.0017702866) q[3];
sx q[3];
rz(-2.6383196) q[3];
sx q[3];
rz(-2.1625904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2718647) q[0];
sx q[0];
rz(-0.80369049) q[0];
sx q[0];
rz(-1.8029689) q[0];
rz(-1.1805234) q[1];
sx q[1];
rz(-1.0931284) q[1];
sx q[1];
rz(0.48042935) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1378616) q[0];
sx q[0];
rz(-2.6380153) q[0];
sx q[0];
rz(-1.1409111) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8992041) q[2];
sx q[2];
rz(-0.17504642) q[2];
sx q[2];
rz(0.72358658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5753382) q[1];
sx q[1];
rz(-1.4713788) q[1];
sx q[1];
rz(1.0097617) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3442743) q[3];
sx q[3];
rz(-1.9206646) q[3];
sx q[3];
rz(-2.1893152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43634811) q[2];
sx q[2];
rz(-0.99506012) q[2];
sx q[2];
rz(1.0850517) q[2];
rz(2.0294225) q[3];
sx q[3];
rz(-2.2472436) q[3];
sx q[3];
rz(-2.3280242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28806624) q[0];
sx q[0];
rz(-2.1948094) q[0];
sx q[0];
rz(-2.1296401) q[0];
rz(-2.2400253) q[1];
sx q[1];
rz(-0.26600599) q[1];
sx q[1];
rz(2.576135) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0386977) q[0];
sx q[0];
rz(-2.2129472) q[0];
sx q[0];
rz(-2.8805634) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1419292) q[2];
sx q[2];
rz(-1.8151917) q[2];
sx q[2];
rz(-1.7606869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3348166) q[1];
sx q[1];
rz(-2.2851181) q[1];
sx q[1];
rz(-3.0187155) q[1];
rz(-3.0363085) q[3];
sx q[3];
rz(-1.255688) q[3];
sx q[3];
rz(-1.1694825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7908343) q[2];
sx q[2];
rz(-1.5956722) q[2];
sx q[2];
rz(1.8809543) q[2];
rz(-2.6993921) q[3];
sx q[3];
rz(-2.111777) q[3];
sx q[3];
rz(-1.7650167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5653391) q[0];
sx q[0];
rz(-0.88247846) q[0];
sx q[0];
rz(-0.89378617) q[0];
rz(-1.5743938) q[1];
sx q[1];
rz(-1.6849453) q[1];
sx q[1];
rz(-1.4243855) q[1];
rz(-0.42548634) q[2];
sx q[2];
rz(-1.31447) q[2];
sx q[2];
rz(-0.18438495) q[2];
rz(1.1461729) q[3];
sx q[3];
rz(-2.7594447) q[3];
sx q[3];
rz(-1.2654163) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
