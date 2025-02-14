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
rz(0.12198099) q[0];
sx q[0];
rz(3.0966336) q[0];
sx q[0];
rz(9.3293204) q[0];
rz(1.8010315) q[1];
sx q[1];
rz(3.0756693) q[1];
sx q[1];
rz(8.8311721) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55768591) q[0];
sx q[0];
rz(-0.61827165) q[0];
sx q[0];
rz(0.64557155) q[0];
rz(-pi) q[1];
rz(-2.4727137) q[2];
sx q[2];
rz(-2.9496585) q[2];
sx q[2];
rz(2.8366983) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0029353142) q[1];
sx q[1];
rz(-1.3627009) q[1];
sx q[1];
rz(-0.63677019) q[1];
x q[2];
rz(1.0450706) q[3];
sx q[3];
rz(-1.330867) q[3];
sx q[3];
rz(-0.6434427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0815214) q[2];
sx q[2];
rz(-1.8895443) q[2];
sx q[2];
rz(-0.56212765) q[2];
rz(0.33956042) q[3];
sx q[3];
rz(-1.5226676) q[3];
sx q[3];
rz(2.0558527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18315166) q[0];
sx q[0];
rz(-2.9957132) q[0];
sx q[0];
rz(-1.1613783) q[0];
rz(-0.24102744) q[1];
sx q[1];
rz(-0.82254219) q[1];
sx q[1];
rz(0.30311146) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89592349) q[0];
sx q[0];
rz(-1.637994) q[0];
sx q[0];
rz(-1.7691866) q[0];
rz(-pi) q[1];
rz(0.14601082) q[2];
sx q[2];
rz(-1.4476629) q[2];
sx q[2];
rz(1.5068123) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1352586) q[1];
sx q[1];
rz(-2.5016682) q[1];
sx q[1];
rz(2.695822) q[1];
x q[2];
rz(1.0996898) q[3];
sx q[3];
rz(-2.0031824) q[3];
sx q[3];
rz(1.3392836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8384398) q[2];
sx q[2];
rz(-1.2829605) q[2];
sx q[2];
rz(0.34070936) q[2];
rz(-2.7029612) q[3];
sx q[3];
rz(-2.3352968) q[3];
sx q[3];
rz(0.058569245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4230147) q[0];
sx q[0];
rz(-0.61841643) q[0];
sx q[0];
rz(-0.83786905) q[0];
rz(0.92309976) q[1];
sx q[1];
rz(-1.2078614) q[1];
sx q[1];
rz(-2.7941678) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7131963) q[0];
sx q[0];
rz(-1.5335463) q[0];
sx q[0];
rz(1.6646882) q[0];
x q[1];
rz(-1.3143888) q[2];
sx q[2];
rz(-1.6056345) q[2];
sx q[2];
rz(-0.23171356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3795137) q[1];
sx q[1];
rz(-2.3752932) q[1];
sx q[1];
rz(2.422782) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33073552) q[3];
sx q[3];
rz(-1.3130672) q[3];
sx q[3];
rz(2.7455519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.122494) q[2];
sx q[2];
rz(-0.118003) q[2];
sx q[2];
rz(-2.4180491) q[2];
rz(-1.8540234) q[3];
sx q[3];
rz(-2.5836594) q[3];
sx q[3];
rz(-3.0961228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36640722) q[0];
sx q[0];
rz(-2.7025096) q[0];
sx q[0];
rz(-2.0570237) q[0];
rz(-1.9547801) q[1];
sx q[1];
rz(-1.2326406) q[1];
sx q[1];
rz(1.67217) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0908244) q[0];
sx q[0];
rz(-1.6009129) q[0];
sx q[0];
rz(1.5539021) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3575063) q[2];
sx q[2];
rz(-2.6166281) q[2];
sx q[2];
rz(0.80277473) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9778825) q[1];
sx q[1];
rz(-2.5228466) q[1];
sx q[1];
rz(2.2676022) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3923548) q[3];
sx q[3];
rz(-0.47296528) q[3];
sx q[3];
rz(1.3441031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12545086) q[2];
sx q[2];
rz(-2.2405388) q[2];
sx q[2];
rz(0.93552843) q[2];
rz(0.15141307) q[3];
sx q[3];
rz(-2.1617523) q[3];
sx q[3];
rz(-2.6305731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22684111) q[0];
sx q[0];
rz(-1.5837357) q[0];
sx q[0];
rz(0.75411183) q[0];
rz(-0.61112815) q[1];
sx q[1];
rz(-1.9978304) q[1];
sx q[1];
rz(-2.4653844) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.138984) q[0];
sx q[0];
rz(-1.4088641) q[0];
sx q[0];
rz(1.2404299) q[0];
rz(-pi) q[1];
rz(-1.0180255) q[2];
sx q[2];
rz(-1.2571196) q[2];
sx q[2];
rz(-0.78757912) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76393781) q[1];
sx q[1];
rz(-1.7309524) q[1];
sx q[1];
rz(0.17367823) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4884454) q[3];
sx q[3];
rz(-2.1550278) q[3];
sx q[3];
rz(-2.1305934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.0040697441) q[2];
sx q[2];
rz(-0.71228945) q[2];
sx q[2];
rz(-0.75312692) q[2];
rz(2.5774041) q[3];
sx q[3];
rz(-2.0337532) q[3];
sx q[3];
rz(-0.75642014) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9516528) q[0];
sx q[0];
rz(-2.6615182) q[0];
sx q[0];
rz(-0.22115627) q[0];
rz(3.0033374) q[1];
sx q[1];
rz(-1.6860551) q[1];
sx q[1];
rz(2.5855605) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.611803) q[0];
sx q[0];
rz(-1.3767813) q[0];
sx q[0];
rz(2.89286) q[0];
rz(-pi) q[1];
rz(0.80433455) q[2];
sx q[2];
rz(-0.51689878) q[2];
sx q[2];
rz(2.3625952) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1631524) q[1];
sx q[1];
rz(-1.5500871) q[1];
sx q[1];
rz(-2.9881575) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8121427) q[3];
sx q[3];
rz(-2.6785319) q[3];
sx q[3];
rz(-1.1487414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.14744645) q[2];
sx q[2];
rz(-2.3799956) q[2];
sx q[2];
rz(-2.2046294) q[2];
rz(-0.41038904) q[3];
sx q[3];
rz(-0.97987163) q[3];
sx q[3];
rz(2.4858937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74528247) q[0];
sx q[0];
rz(-0.87438011) q[0];
sx q[0];
rz(-2.1790047) q[0];
rz(0.67086041) q[1];
sx q[1];
rz(-0.52898359) q[1];
sx q[1];
rz(-0.76505351) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017185171) q[0];
sx q[0];
rz(-1.1721503) q[0];
sx q[0];
rz(0.17332698) q[0];
rz(-pi) q[1];
rz(0.97114886) q[2];
sx q[2];
rz(-3.0191688) q[2];
sx q[2];
rz(0.79736191) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9967186) q[1];
sx q[1];
rz(-1.178374) q[1];
sx q[1];
rz(2.1226354) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9387127) q[3];
sx q[3];
rz(-1.6867609) q[3];
sx q[3];
rz(-2.396317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44136167) q[2];
sx q[2];
rz(-0.92741489) q[2];
sx q[2];
rz(-0.78511635) q[2];
rz(-2.650812) q[3];
sx q[3];
rz(-1.9793341) q[3];
sx q[3];
rz(-0.47500113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15071507) q[0];
sx q[0];
rz(-3.0325723) q[0];
sx q[0];
rz(0.14608598) q[0];
rz(-0.27169216) q[1];
sx q[1];
rz(-2.3482595) q[1];
sx q[1];
rz(-0.21154107) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78286029) q[0];
sx q[0];
rz(-1.4050589) q[0];
sx q[0];
rz(-0.77605794) q[0];
rz(-pi) q[1];
rz(-0.11577932) q[2];
sx q[2];
rz(-1.2010788) q[2];
sx q[2];
rz(-3.0860975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.52741963) q[1];
sx q[1];
rz(-2.895311) q[1];
sx q[1];
rz(2.8524532) q[1];
rz(-pi) q[2];
rz(1.7320427) q[3];
sx q[3];
rz(-0.60617709) q[3];
sx q[3];
rz(-2.2796749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3286256) q[2];
sx q[2];
rz(-2.2058637) q[2];
sx q[2];
rz(-2.3198371) q[2];
rz(1.3385319) q[3];
sx q[3];
rz(-2.0526363) q[3];
sx q[3];
rz(0.72430044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61447918) q[0];
sx q[0];
rz(-2.4452657) q[0];
sx q[0];
rz(1.7823727) q[0];
rz(2.1837153) q[1];
sx q[1];
rz(-0.33708894) q[1];
sx q[1];
rz(-2.8639796) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.255186) q[0];
sx q[0];
rz(-2.1069114) q[0];
sx q[0];
rz(-0.84250959) q[0];
x q[1];
rz(-2.3657655) q[2];
sx q[2];
rz(-1.3444573) q[2];
sx q[2];
rz(-0.090858484) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.094885849) q[1];
sx q[1];
rz(-1.1562655) q[1];
sx q[1];
rz(-1.4317896) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7630799) q[3];
sx q[3];
rz(-1.4106102) q[3];
sx q[3];
rz(1.9122363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9684888) q[2];
sx q[2];
rz(-2.5471881) q[2];
sx q[2];
rz(-1.8572726) q[2];
rz(-2.3641018) q[3];
sx q[3];
rz(-2.5992664) q[3];
sx q[3];
rz(-0.29977453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023329968) q[0];
sx q[0];
rz(-1.6602004) q[0];
sx q[0];
rz(-0.0059286038) q[0];
rz(1.802035) q[1];
sx q[1];
rz(-2.1052994) q[1];
sx q[1];
rz(2.6609227) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6193793) q[0];
sx q[0];
rz(-1.5565158) q[0];
sx q[0];
rz(-3.1254053) q[0];
rz(-0.43227682) q[2];
sx q[2];
rz(-2.1735809) q[2];
sx q[2];
rz(-0.31291459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1263894) q[1];
sx q[1];
rz(-1.1600094) q[1];
sx q[1];
rz(0.65151626) q[1];
x q[2];
rz(1.3044429) q[3];
sx q[3];
rz(-1.7031809) q[3];
sx q[3];
rz(1.799982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19114755) q[2];
sx q[2];
rz(-0.55926776) q[2];
sx q[2];
rz(2.9648103) q[2];
rz(-2.2117129) q[3];
sx q[3];
rz(-1.1888489) q[3];
sx q[3];
rz(-1.0298347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69242351) q[0];
sx q[0];
rz(-1.7208736) q[0];
sx q[0];
rz(2.0509913) q[0];
rz(1.0724267) q[1];
sx q[1];
rz(-1.4785531) q[1];
sx q[1];
rz(-1.4428152) q[1];
rz(0.011913997) q[2];
sx q[2];
rz(-0.49187751) q[2];
sx q[2];
rz(-0.054452489) q[2];
rz(0.60521331) q[3];
sx q[3];
rz(-1.6463668) q[3];
sx q[3];
rz(2.2201408) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
