OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.469874) q[0];
sx q[0];
rz(-1.9865541) q[0];
sx q[0];
rz(1.7116829) q[0];
rz(0.45541304) q[1];
sx q[1];
rz(5.6918511) q[1];
sx q[1];
rz(11.439352) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9874632) q[0];
sx q[0];
rz(-0.23079458) q[0];
sx q[0];
rz(-1.6213645) q[0];
rz(-pi) q[1];
rz(-0.65140407) q[2];
sx q[2];
rz(-2.7198875) q[2];
sx q[2];
rz(2.5159474) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0309567) q[1];
sx q[1];
rz(-1.0070166) q[1];
sx q[1];
rz(-2.8085676) q[1];
rz(-pi) q[2];
x q[2];
rz(0.092394775) q[3];
sx q[3];
rz(-0.1915598) q[3];
sx q[3];
rz(2.5026623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0226125) q[2];
sx q[2];
rz(-0.70465124) q[2];
sx q[2];
rz(-1.653778) q[2];
rz(-2.7704499) q[3];
sx q[3];
rz(-1.1266339) q[3];
sx q[3];
rz(0.77905542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30598518) q[0];
sx q[0];
rz(-1.7636517) q[0];
sx q[0];
rz(-2.4775179) q[0];
rz(-0.56655073) q[1];
sx q[1];
rz(-1.5357176) q[1];
sx q[1];
rz(0.18633349) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9772352) q[0];
sx q[0];
rz(-2.1002227) q[0];
sx q[0];
rz(0.25821547) q[0];
x q[1];
rz(1.1900224) q[2];
sx q[2];
rz(-1.8203893) q[2];
sx q[2];
rz(-1.4480424) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1809606) q[1];
sx q[1];
rz(-2.4419028) q[1];
sx q[1];
rz(-1.8262499) q[1];
x q[2];
rz(1.6870895) q[3];
sx q[3];
rz(-0.52403677) q[3];
sx q[3];
rz(-0.91479036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8191007) q[2];
sx q[2];
rz(-1.842061) q[2];
sx q[2];
rz(-1.5619649) q[2];
rz(-1.1473848) q[3];
sx q[3];
rz(-0.29427823) q[3];
sx q[3];
rz(-1.3636205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550302) q[0];
sx q[0];
rz(-1.6494305) q[0];
sx q[0];
rz(0.30558875) q[0];
rz(1.2809523) q[1];
sx q[1];
rz(-0.31550229) q[1];
sx q[1];
rz(-1.5234647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8620096) q[0];
sx q[0];
rz(-1.9317937) q[0];
sx q[0];
rz(1.8185413) q[0];
x q[1];
rz(-0.045133682) q[2];
sx q[2];
rz(-1.5861057) q[2];
sx q[2];
rz(1.65427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9429607) q[1];
sx q[1];
rz(-2.3395798) q[1];
sx q[1];
rz(-0.50870163) q[1];
rz(-1.9033857) q[3];
sx q[3];
rz(-2.5032515) q[3];
sx q[3];
rz(-2.7901543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4898701) q[2];
sx q[2];
rz(-1.8946596) q[2];
sx q[2];
rz(-2.9215802) q[2];
rz(-3.0618727) q[3];
sx q[3];
rz(-0.97697512) q[3];
sx q[3];
rz(2.2102833) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43059573) q[0];
sx q[0];
rz(-1.7544704) q[0];
sx q[0];
rz(-3.1331449) q[0];
rz(-1.9795817) q[1];
sx q[1];
rz(-1.153667) q[1];
sx q[1];
rz(-2.0268424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8793068) q[0];
sx q[0];
rz(-1.7588076) q[0];
sx q[0];
rz(1.2557058) q[0];
x q[1];
rz(-2.3643119) q[2];
sx q[2];
rz(-0.72597625) q[2];
sx q[2];
rz(2.7216146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9544008) q[1];
sx q[1];
rz(-1.9241417) q[1];
sx q[1];
rz(-0.29275972) q[1];
rz(1.1590121) q[3];
sx q[3];
rz(-1.1819541) q[3];
sx q[3];
rz(-1.4847311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4343425) q[2];
sx q[2];
rz(-1.0184526) q[2];
sx q[2];
rz(0.84558359) q[2];
rz(2.8905458) q[3];
sx q[3];
rz(-1.8699402) q[3];
sx q[3];
rz(2.3959851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4276328) q[0];
sx q[0];
rz(-0.41249713) q[0];
sx q[0];
rz(1.0300256) q[0];
rz(0.59459844) q[1];
sx q[1];
rz(-1.7183036) q[1];
sx q[1];
rz(-1.7880218) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0145688) q[0];
sx q[0];
rz(-1.5480124) q[0];
sx q[0];
rz(1.5915074) q[0];
x q[1];
rz(0.51026235) q[2];
sx q[2];
rz(-2.7760421) q[2];
sx q[2];
rz(-2.9175959) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8882519) q[1];
sx q[1];
rz(-2.2512813) q[1];
sx q[1];
rz(-0.91618211) q[1];
rz(-0.59242512) q[3];
sx q[3];
rz(-1.1325784) q[3];
sx q[3];
rz(0.93417227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9109965) q[2];
sx q[2];
rz(-1.6350919) q[2];
sx q[2];
rz(0.050749151) q[2];
rz(1.1132318) q[3];
sx q[3];
rz(-1.166393) q[3];
sx q[3];
rz(0.39065233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0942866) q[0];
sx q[0];
rz(-1.0599437) q[0];
sx q[0];
rz(-0.37288368) q[0];
rz(-1.8376384) q[1];
sx q[1];
rz(-1.8770437) q[1];
sx q[1];
rz(-1.0252999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7367497) q[0];
sx q[0];
rz(-1.4944634) q[0];
sx q[0];
rz(-0.45053225) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16496678) q[2];
sx q[2];
rz(-1.8968762) q[2];
sx q[2];
rz(-0.13298377) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91955012) q[1];
sx q[1];
rz(-1.371438) q[1];
sx q[1];
rz(-2.2712703) q[1];
rz(-pi) q[2];
rz(-2.0399569) q[3];
sx q[3];
rz(-2.6734452) q[3];
sx q[3];
rz(-0.461138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.23124) q[2];
sx q[2];
rz(-2.9515036) q[2];
sx q[2];
rz(1.4300038) q[2];
rz(1.6010239) q[3];
sx q[3];
rz(-2.4547596) q[3];
sx q[3];
rz(-1.6704667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1533399) q[0];
sx q[0];
rz(-1.8415469) q[0];
sx q[0];
rz(-2.7100995) q[0];
rz(1.2639812) q[1];
sx q[1];
rz(-1.6004205) q[1];
sx q[1];
rz(-0.65013179) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1558485) q[0];
sx q[0];
rz(-1.2835644) q[0];
sx q[0];
rz(1.656507) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2540497) q[2];
sx q[2];
rz(-1.6517963) q[2];
sx q[2];
rz(2.0320867) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.35061478) q[1];
sx q[1];
rz(-0.87360604) q[1];
sx q[1];
rz(0.15048262) q[1];
rz(-pi) q[2];
rz(1.1129208) q[3];
sx q[3];
rz(-1.5639389) q[3];
sx q[3];
rz(-2.6049861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5905137) q[2];
sx q[2];
rz(-1.3912018) q[2];
sx q[2];
rz(-0.004465731) q[2];
rz(-3.1320069) q[3];
sx q[3];
rz(-1.3349345) q[3];
sx q[3];
rz(-2.7981304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15653217) q[0];
sx q[0];
rz(-2.8592906) q[0];
sx q[0];
rz(-2.9811133) q[0];
rz(-1.7566682) q[1];
sx q[1];
rz(-2.3424708) q[1];
sx q[1];
rz(-2.8327732) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35058103) q[0];
sx q[0];
rz(-1.6948943) q[0];
sx q[0];
rz(1.3826293) q[0];
x q[1];
rz(2.7032531) q[2];
sx q[2];
rz(-0.76073863) q[2];
sx q[2];
rz(2.1957514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0817326) q[1];
sx q[1];
rz(-2.5860021) q[1];
sx q[1];
rz(-0.76738552) q[1];
rz(-pi) q[2];
rz(0.71864031) q[3];
sx q[3];
rz(-2.6860533) q[3];
sx q[3];
rz(2.9589341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7979692) q[2];
sx q[2];
rz(-0.60911959) q[2];
sx q[2];
rz(1.5416175) q[2];
rz(-0.68495098) q[3];
sx q[3];
rz(-1.5545605) q[3];
sx q[3];
rz(-1.6182914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5110382) q[0];
sx q[0];
rz(-1.7031952) q[0];
sx q[0];
rz(2.5901929) q[0];
rz(-0.4772056) q[1];
sx q[1];
rz(-0.91608945) q[1];
sx q[1];
rz(-0.83470693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.793551) q[0];
sx q[0];
rz(-2.0299596) q[0];
sx q[0];
rz(-0.55450098) q[0];
x q[1];
rz(0.040933523) q[2];
sx q[2];
rz(-2.4403095) q[2];
sx q[2];
rz(-2.3778039) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1795802) q[1];
sx q[1];
rz(-2.5287345) q[1];
sx q[1];
rz(1.5108766) q[1];
x q[2];
rz(2.1797997) q[3];
sx q[3];
rz(-2.0941705) q[3];
sx q[3];
rz(-0.23853049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7029552) q[2];
sx q[2];
rz(-0.96200395) q[2];
sx q[2];
rz(2.9023602) q[2];
rz(-2.8171825) q[3];
sx q[3];
rz(-1.4374461) q[3];
sx q[3];
rz(1.602406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6259554) q[0];
sx q[0];
rz(-0.13787585) q[0];
sx q[0];
rz(-0.71664083) q[0];
rz(-0.58249885) q[1];
sx q[1];
rz(-1.7889675) q[1];
sx q[1];
rz(-2.8315721) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17663613) q[0];
sx q[0];
rz(-0.88439098) q[0];
sx q[0];
rz(2.6841037) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7298195) q[2];
sx q[2];
rz(-1.207282) q[2];
sx q[2];
rz(-0.95136729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2523732) q[1];
sx q[1];
rz(-0.82102005) q[1];
sx q[1];
rz(1.0733502) q[1];
rz(0.54599793) q[3];
sx q[3];
rz(-1.1359478) q[3];
sx q[3];
rz(-0.42213371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6325355) q[2];
sx q[2];
rz(-2.5639503) q[2];
sx q[2];
rz(-1.4198111) q[2];
rz(1.696473) q[3];
sx q[3];
rz(-2.4308379) q[3];
sx q[3];
rz(-0.76326171) q[3];
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
rz(1.9025018) q[0];
sx q[0];
rz(-0.9466753) q[0];
sx q[0];
rz(-1.7539903) q[0];
rz(1.4801964) q[1];
sx q[1];
rz(-1.4085242) q[1];
sx q[1];
rz(-2.9396802) q[1];
rz(1.4954064) q[2];
sx q[2];
rz(-1.8373377) q[2];
sx q[2];
rz(-1.1129825) q[2];
rz(-2.1275413) q[3];
sx q[3];
rz(-2.378856) q[3];
sx q[3];
rz(3.0434276) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
