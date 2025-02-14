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
rz(2.9003484) q[0];
sx q[0];
rz(3.3527346) q[0];
sx q[0];
rz(9.0169173) q[0];
rz(-1.7653699) q[1];
sx q[1];
rz(5.0566109) q[1];
sx q[1];
rz(15.772718) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61633793) q[0];
sx q[0];
rz(-3.0260575) q[0];
sx q[0];
rz(-1.6166685) q[0];
rz(-0.98528905) q[2];
sx q[2];
rz(-1.3006214) q[2];
sx q[2];
rz(-0.5952685) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.237707) q[1];
sx q[1];
rz(-0.14169417) q[1];
sx q[1];
rz(-0.078448729) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14459855) q[3];
sx q[3];
rz(-0.62718348) q[3];
sx q[3];
rz(2.7759564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.18225081) q[2];
sx q[2];
rz(-2.0472517) q[2];
sx q[2];
rz(1.135929) q[2];
rz(0.7302537) q[3];
sx q[3];
rz(-1.3948995) q[3];
sx q[3];
rz(-1.6212538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092904329) q[0];
sx q[0];
rz(-2.1850259) q[0];
sx q[0];
rz(-0.055305716) q[0];
rz(-0.37120184) q[1];
sx q[1];
rz(-0.35938811) q[1];
sx q[1];
rz(-2.7296383) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7474351) q[0];
sx q[0];
rz(-2.0965946) q[0];
sx q[0];
rz(1.6461685) q[0];
rz(-pi) q[1];
rz(2.8438004) q[2];
sx q[2];
rz(-1.2504559) q[2];
sx q[2];
rz(-3.0646119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4186379) q[1];
sx q[1];
rz(-1.4297863) q[1];
sx q[1];
rz(-0.73608558) q[1];
rz(-pi) q[2];
rz(2.2188236) q[3];
sx q[3];
rz(-1.8776181) q[3];
sx q[3];
rz(2.0056779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5281795) q[2];
sx q[2];
rz(-1.4643022) q[2];
sx q[2];
rz(0.8826274) q[2];
rz(-1.4334076) q[3];
sx q[3];
rz(-1.9288758) q[3];
sx q[3];
rz(2.9679756) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16647896) q[0];
sx q[0];
rz(-3.1298895) q[0];
sx q[0];
rz(0.56667462) q[0];
rz(-0.4176248) q[1];
sx q[1];
rz(-1.6028701) q[1];
sx q[1];
rz(-1.8142987) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6782129) q[0];
sx q[0];
rz(-0.87317077) q[0];
sx q[0];
rz(-0.11059983) q[0];
x q[1];
rz(-1.3603052) q[2];
sx q[2];
rz(-2.5062525) q[2];
sx q[2];
rz(0.47266211) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9108563) q[1];
sx q[1];
rz(-1.5720782) q[1];
sx q[1];
rz(-1.8788473) q[1];
x q[2];
rz(2.684115) q[3];
sx q[3];
rz(-0.8972392) q[3];
sx q[3];
rz(-0.20041538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7495482) q[2];
sx q[2];
rz(-1.0907402) q[2];
sx q[2];
rz(0.90847477) q[2];
rz(2.3467017) q[3];
sx q[3];
rz(-2.0386233) q[3];
sx q[3];
rz(-3.095043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4817151) q[0];
sx q[0];
rz(-0.91232038) q[0];
sx q[0];
rz(2.4758441) q[0];
rz(2.2944229) q[1];
sx q[1];
rz(-2.8916292) q[1];
sx q[1];
rz(0.4096823) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2100135) q[0];
sx q[0];
rz(-1.0780452) q[0];
sx q[0];
rz(-1.7442305) q[0];
rz(-pi) q[1];
rz(0.3607765) q[2];
sx q[2];
rz(-0.84178001) q[2];
sx q[2];
rz(-1.3363802) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.078881513) q[1];
sx q[1];
rz(-2.1797175) q[1];
sx q[1];
rz(-2.8669861) q[1];
x q[2];
rz(-1.631103) q[3];
sx q[3];
rz(-0.44476393) q[3];
sx q[3];
rz(-0.30032762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72982558) q[2];
sx q[2];
rz(-2.4335786) q[2];
sx q[2];
rz(-1.3789619) q[2];
rz(1.1369368) q[3];
sx q[3];
rz(-1.7863019) q[3];
sx q[3];
rz(2.304145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9627422) q[0];
sx q[0];
rz(-1.0609635) q[0];
sx q[0];
rz(2.7233997) q[0];
rz(1.7182982) q[1];
sx q[1];
rz(-1.9530674) q[1];
sx q[1];
rz(1.6421912) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20883501) q[0];
sx q[0];
rz(-1.3262188) q[0];
sx q[0];
rz(-2.6840058) q[0];
rz(-pi) q[1];
rz(2.4046481) q[2];
sx q[2];
rz(-0.43790753) q[2];
sx q[2];
rz(-2.2543668) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.114257) q[1];
sx q[1];
rz(-1.789195) q[1];
sx q[1];
rz(-3.0824003) q[1];
x q[2];
rz(1.1058331) q[3];
sx q[3];
rz(-0.14497193) q[3];
sx q[3];
rz(-1.5752058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2749918) q[2];
sx q[2];
rz(-1.9860705) q[2];
sx q[2];
rz(0.17821136) q[2];
rz(-2.7675659) q[3];
sx q[3];
rz(-2.1686797) q[3];
sx q[3];
rz(-1.1427574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48443925) q[0];
sx q[0];
rz(-1.3244119) q[0];
sx q[0];
rz(1.665218) q[0];
rz(-0.58552512) q[1];
sx q[1];
rz(-0.89591566) q[1];
sx q[1];
rz(2.4077328) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6505665) q[0];
sx q[0];
rz(-1.125245) q[0];
sx q[0];
rz(-2.4413021) q[0];
x q[1];
rz(1.2452872) q[2];
sx q[2];
rz(-1.8625038) q[2];
sx q[2];
rz(-2.0355647) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31325737) q[1];
sx q[1];
rz(-2.0861107) q[1];
sx q[1];
rz(0.22586598) q[1];
rz(-pi) q[2];
rz(1.9820342) q[3];
sx q[3];
rz(-0.48732584) q[3];
sx q[3];
rz(-3.0912378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0840941) q[2];
sx q[2];
rz(-1.7666662) q[2];
sx q[2];
rz(1.6424087) q[2];
rz(1.5205787) q[3];
sx q[3];
rz(-2.4937544) q[3];
sx q[3];
rz(-2.9852941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75477377) q[0];
sx q[0];
rz(-0.91803011) q[0];
sx q[0];
rz(1.7479489) q[0];
rz(0.58760324) q[1];
sx q[1];
rz(-1.6410476) q[1];
sx q[1];
rz(2.844152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7845532) q[0];
sx q[0];
rz(-1.2149356) q[0];
sx q[0];
rz(-0.83609348) q[0];
rz(3.1226899) q[2];
sx q[2];
rz(-1.3002965) q[2];
sx q[2];
rz(2.8342385) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.12402) q[1];
sx q[1];
rz(-2.4618853) q[1];
sx q[1];
rz(1.5862096) q[1];
x q[2];
rz(0.92771156) q[3];
sx q[3];
rz(-1.7684019) q[3];
sx q[3];
rz(1.2090982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8224767) q[2];
sx q[2];
rz(-1.3603223) q[2];
sx q[2];
rz(-0.3817257) q[2];
rz(2.5413359) q[3];
sx q[3];
rz(-2.468942) q[3];
sx q[3];
rz(2.8762347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6622019) q[0];
sx q[0];
rz(-1.5150161) q[0];
sx q[0];
rz(-1.3375244) q[0];
rz(0.0094825347) q[1];
sx q[1];
rz(-0.9597221) q[1];
sx q[1];
rz(1.3710075) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0616546) q[0];
sx q[0];
rz(-1.9091055) q[0];
sx q[0];
rz(-1.7163971) q[0];
x q[1];
rz(2.2393537) q[2];
sx q[2];
rz(-1.1663365) q[2];
sx q[2];
rz(2.8079833) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1320033) q[1];
sx q[1];
rz(-1.4456962) q[1];
sx q[1];
rz(-0.78937197) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6608606) q[3];
sx q[3];
rz(-1.9833296) q[3];
sx q[3];
rz(-1.6270669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.81801549) q[2];
sx q[2];
rz(-0.90138268) q[2];
sx q[2];
rz(0.14249194) q[2];
rz(0.88322181) q[3];
sx q[3];
rz(-1.764069) q[3];
sx q[3];
rz(1.6787329) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.091081) q[0];
sx q[0];
rz(-0.091878042) q[0];
sx q[0];
rz(-1.8930513) q[0];
rz(-2.7342791) q[1];
sx q[1];
rz(-1.7946323) q[1];
sx q[1];
rz(-1.9479082) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3996959) q[0];
sx q[0];
rz(-0.4238766) q[0];
sx q[0];
rz(-0.90452832) q[0];
rz(2.2676554) q[2];
sx q[2];
rz(-2.1997582) q[2];
sx q[2];
rz(-3.0203117) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5919009) q[1];
sx q[1];
rz(-0.61335269) q[1];
sx q[1];
rz(2.5332622) q[1];
rz(-pi) q[2];
x q[2];
rz(3.3832559e-05) q[3];
sx q[3];
rz(-1.8305998) q[3];
sx q[3];
rz(0.75601789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.074389) q[2];
sx q[2];
rz(-2.1838102) q[2];
sx q[2];
rz(2.0286782) q[2];
rz(-0.042996081) q[3];
sx q[3];
rz(-1.4837416) q[3];
sx q[3];
rz(-2.903741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064903108) q[0];
sx q[0];
rz(-1.4142798) q[0];
sx q[0];
rz(-0.195737) q[0];
rz(1.2204569) q[1];
sx q[1];
rz(-2.7595322) q[1];
sx q[1];
rz(-0.97741309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.138151) q[0];
sx q[0];
rz(-1.8288946) q[0];
sx q[0];
rz(-1.6736223) q[0];
rz(-pi) q[1];
rz(2.9896792) q[2];
sx q[2];
rz(-1.8679233) q[2];
sx q[2];
rz(-3.0346485) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4949242) q[1];
sx q[1];
rz(-0.99993247) q[1];
sx q[1];
rz(-2.5350556) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5217358) q[3];
sx q[3];
rz(-2.6708948) q[3];
sx q[3];
rz(2.6420267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1310554) q[2];
sx q[2];
rz(-0.26630339) q[2];
sx q[2];
rz(1.4492501) q[2];
rz(2.2629755) q[3];
sx q[3];
rz(-2.0802458) q[3];
sx q[3];
rz(1.7884458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
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
rz(1.5129678) q[0];
sx q[0];
rz(-1.9504564) q[0];
sx q[0];
rz(1.7497028) q[0];
rz(-1.7616918) q[1];
sx q[1];
rz(-1.6086144) q[1];
sx q[1];
rz(-0.10133941) q[1];
rz(-2.9607282) q[2];
sx q[2];
rz(-1.7500063) q[2];
sx q[2];
rz(1.0371006) q[2];
rz(0.64639277) q[3];
sx q[3];
rz(-1.7714995) q[3];
sx q[3];
rz(-2.5378791) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
