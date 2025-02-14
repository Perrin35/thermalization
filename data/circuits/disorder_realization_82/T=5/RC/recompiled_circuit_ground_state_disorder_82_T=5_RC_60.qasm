OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6717186) q[0];
sx q[0];
rz(-1.1550386) q[0];
sx q[0];
rz(-1.7116829) q[0];
rz(-2.6861796) q[1];
sx q[1];
rz(-2.5502584) q[1];
sx q[1];
rz(1.1270181) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4658964) q[0];
sx q[0];
rz(-1.5823592) q[0];
sx q[0];
rz(1.8013062) q[0];
rz(1.8363709) q[2];
sx q[2];
rz(-1.2392534) q[2];
sx q[2];
rz(-1.8198554) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8642173) q[1];
sx q[1];
rz(-1.2908415) q[1];
sx q[1];
rz(0.98119189) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0491979) q[3];
sx q[3];
rz(-2.9500329) q[3];
sx q[3];
rz(2.5026623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.11898018) q[2];
sx q[2];
rz(-0.70465124) q[2];
sx q[2];
rz(1.653778) q[2];
rz(0.37114272) q[3];
sx q[3];
rz(-1.1266339) q[3];
sx q[3];
rz(-2.3625372) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8356075) q[0];
sx q[0];
rz(-1.7636517) q[0];
sx q[0];
rz(-2.4775179) q[0];
rz(-2.5750419) q[1];
sx q[1];
rz(-1.5357176) q[1];
sx q[1];
rz(-0.18633349) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16435745) q[0];
sx q[0];
rz(-2.1002227) q[0];
sx q[0];
rz(-0.25821547) q[0];
x q[1];
rz(-1.1900224) q[2];
sx q[2];
rz(-1.3212034) q[2];
sx q[2];
rz(1.6935503) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6316698) q[1];
sx q[1];
rz(-2.2434739) q[1];
sx q[1];
rz(2.9320168) q[1];
rz(-pi) q[2];
rz(1.6870895) q[3];
sx q[3];
rz(-2.6175559) q[3];
sx q[3];
rz(0.91479036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8191007) q[2];
sx q[2];
rz(-1.842061) q[2];
sx q[2];
rz(1.5796278) q[2];
rz(1.1473848) q[3];
sx q[3];
rz(-2.8473144) q[3];
sx q[3];
rz(-1.3636205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086562432) q[0];
sx q[0];
rz(-1.4921621) q[0];
sx q[0];
rz(-2.8360039) q[0];
rz(1.2809523) q[1];
sx q[1];
rz(-2.8260904) q[1];
sx q[1];
rz(1.5234647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9394835) q[0];
sx q[0];
rz(-1.8022707) q[0];
sx q[0];
rz(0.37136308) q[0];
rz(-pi) q[1];
rz(1.5554713) q[2];
sx q[2];
rz(-1.6159247) q[2];
sx q[2];
rz(3.0588104) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1393006) q[1];
sx q[1];
rz(-1.2131547) q[1];
sx q[1];
rz(2.4071818) q[1];
rz(0.23766808) q[3];
sx q[3];
rz(-0.97252995) q[3];
sx q[3];
rz(0.054758398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.65172255) q[2];
sx q[2];
rz(-1.8946596) q[2];
sx q[2];
rz(0.22001246) q[2];
rz(3.0618727) q[3];
sx q[3];
rz(-0.97697512) q[3];
sx q[3];
rz(0.93130934) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43059573) q[0];
sx q[0];
rz(-1.3871223) q[0];
sx q[0];
rz(-0.0084477607) q[0];
rz(-1.9795817) q[1];
sx q[1];
rz(-1.153667) q[1];
sx q[1];
rz(-2.0268424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8793068) q[0];
sx q[0];
rz(-1.7588076) q[0];
sx q[0];
rz(-1.2557058) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1276541) q[2];
sx q[2];
rz(-1.077855) q[2];
sx q[2];
rz(1.8007474) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1871918) q[1];
sx q[1];
rz(-1.2174509) q[1];
sx q[1];
rz(0.29275972) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3678701) q[3];
sx q[3];
rz(-0.55860177) q[3];
sx q[3];
rz(-0.62880558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4343425) q[2];
sx q[2];
rz(-2.12314) q[2];
sx q[2];
rz(0.84558359) q[2];
rz(2.8905458) q[3];
sx q[3];
rz(-1.2716525) q[3];
sx q[3];
rz(0.74560753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4276328) q[0];
sx q[0];
rz(-2.7290955) q[0];
sx q[0];
rz(-1.0300256) q[0];
rz(2.5469942) q[1];
sx q[1];
rz(-1.7183036) q[1];
sx q[1];
rz(1.7880218) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524588) q[0];
sx q[0];
rz(-3.1108034) q[0];
sx q[0];
rz(-2.4039387) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7556023) q[2];
sx q[2];
rz(-1.2535742) q[2];
sx q[2];
rz(-0.76390195) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.7673257) q[1];
sx q[1];
rz(-1.077829) q[1];
sx q[1];
rz(-2.3460991) q[1];
rz(-0.59242512) q[3];
sx q[3];
rz(-1.1325784) q[3];
sx q[3];
rz(-2.2074204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2305962) q[2];
sx q[2];
rz(-1.6350919) q[2];
sx q[2];
rz(3.0908435) q[2];
rz(-2.0283608) q[3];
sx q[3];
rz(-1.166393) q[3];
sx q[3];
rz(-2.7509403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0942866) q[0];
sx q[0];
rz(-1.0599437) q[0];
sx q[0];
rz(0.37288368) q[0];
rz(-1.3039543) q[1];
sx q[1];
rz(-1.8770437) q[1];
sx q[1];
rz(-2.1162927) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4639125) q[0];
sx q[0];
rz(-2.6850774) q[0];
sx q[0];
rz(-2.9677261) q[0];
rz(-2.0228902) q[2];
sx q[2];
rz(-0.36410022) q[2];
sx q[2];
rz(-0.61227476) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91955012) q[1];
sx q[1];
rz(-1.371438) q[1];
sx q[1];
rz(-2.2712703) q[1];
x q[2];
rz(-2.9168374) q[3];
sx q[3];
rz(-1.9850176) q[3];
sx q[3];
rz(-0.97771588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9103526) q[2];
sx q[2];
rz(-0.19008907) q[2];
sx q[2];
rz(-1.7115889) q[2];
rz(1.6010239) q[3];
sx q[3];
rz(-2.4547596) q[3];
sx q[3];
rz(1.471126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1533399) q[0];
sx q[0];
rz(-1.3000458) q[0];
sx q[0];
rz(-2.7100995) q[0];
rz(1.2639812) q[1];
sx q[1];
rz(-1.6004205) q[1];
sx q[1];
rz(-0.65013179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1558485) q[0];
sx q[0];
rz(-1.8580282) q[0];
sx q[0];
rz(-1.656507) q[0];
rz(-pi) q[1];
rz(-2.2540497) q[2];
sx q[2];
rz(-1.4897963) q[2];
sx q[2];
rz(-2.0320867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35061478) q[1];
sx q[1];
rz(-2.2679866) q[1];
sx q[1];
rz(2.99111) q[1];
rz(-pi) q[2];
rz(-1.5552844) q[3];
sx q[3];
rz(-0.45792327) q[3];
sx q[3];
rz(1.0481038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5905137) q[2];
sx q[2];
rz(-1.7503909) q[2];
sx q[2];
rz(0.004465731) q[2];
rz(3.1320069) q[3];
sx q[3];
rz(-1.3349345) q[3];
sx q[3];
rz(-0.34346223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9850605) q[0];
sx q[0];
rz(-2.8592906) q[0];
sx q[0];
rz(-0.16047934) q[0];
rz(-1.7566682) q[1];
sx q[1];
rz(-0.79912186) q[1];
sx q[1];
rz(2.8327732) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.243781) q[0];
sx q[0];
rz(-1.7574991) q[0];
sx q[0];
rz(-3.0152882) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1868434) q[2];
sx q[2];
rz(-0.89659474) q[2];
sx q[2];
rz(0.37146682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9659979) q[1];
sx q[1];
rz(-1.9456989) q[1];
sx q[1];
rz(2.7213827) q[1];
rz(-0.35328226) q[3];
sx q[3];
rz(-1.2769413) q[3];
sx q[3];
rz(2.053956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7979692) q[2];
sx q[2];
rz(-0.60911959) q[2];
sx q[2];
rz(-1.5416175) q[2];
rz(2.4566417) q[3];
sx q[3];
rz(-1.5545605) q[3];
sx q[3];
rz(1.5233013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5110382) q[0];
sx q[0];
rz(-1.7031952) q[0];
sx q[0];
rz(-2.5901929) q[0];
rz(0.4772056) q[1];
sx q[1];
rz(-0.91608945) q[1];
sx q[1];
rz(-2.3068857) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3480417) q[0];
sx q[0];
rz(-2.0299596) q[0];
sx q[0];
rz(2.5870917) q[0];
x q[1];
rz(-2.4407225) q[2];
sx q[2];
rz(-1.5972023) q[2];
sx q[2];
rz(2.3658662) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0352201) q[1];
sx q[1];
rz(-2.1823931) q[1];
sx q[1];
rz(3.0995083) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5284994) q[3];
sx q[3];
rz(-1.0524155) q[3];
sx q[3];
rz(1.4739153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7029552) q[2];
sx q[2];
rz(-0.96200395) q[2];
sx q[2];
rz(0.23923242) q[2];
rz(2.8171825) q[3];
sx q[3];
rz(-1.7041465) q[3];
sx q[3];
rz(-1.5391866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5156373) q[0];
sx q[0];
rz(-0.13787585) q[0];
sx q[0];
rz(-0.71664083) q[0];
rz(0.58249885) q[1];
sx q[1];
rz(-1.3526252) q[1];
sx q[1];
rz(-2.8315721) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6577756) q[0];
sx q[0];
rz(-0.80379009) q[0];
sx q[0];
rz(-1.0763542) q[0];
rz(1.9642682) q[2];
sx q[2];
rz(-1.1873909) q[2];
sx q[2];
rz(-2.368106) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0359874) q[1];
sx q[1];
rz(-1.2140555) q[1];
sx q[1];
rz(2.3272661) q[1];
x q[2];
rz(2.0686555) q[3];
sx q[3];
rz(-1.080435) q[3];
sx q[3];
rz(2.2435202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.50905716) q[2];
sx q[2];
rz(-2.5639503) q[2];
sx q[2];
rz(1.7217815) q[2];
rz(-1.696473) q[3];
sx q[3];
rz(-0.71075478) q[3];
sx q[3];
rz(-0.76326171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9025018) q[0];
sx q[0];
rz(-2.1949174) q[0];
sx q[0];
rz(1.3876023) q[0];
rz(1.4801964) q[1];
sx q[1];
rz(-1.4085242) q[1];
sx q[1];
rz(-2.9396802) q[1];
rz(2.8724332) q[2];
sx q[2];
rz(-0.27675376) q[2];
sx q[2];
rz(-1.3922538) q[2];
rz(1.0140513) q[3];
sx q[3];
rz(-2.378856) q[3];
sx q[3];
rz(3.0434276) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
