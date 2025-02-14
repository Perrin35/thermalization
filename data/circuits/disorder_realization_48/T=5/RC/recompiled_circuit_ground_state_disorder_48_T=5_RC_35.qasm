OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(-1.1416924) q[0];
sx q[0];
rz(2.8853631) q[0];
rz(3.5186634) q[1];
sx q[1];
rz(4.9405603) q[1];
sx q[1];
rz(8.3648051) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1621583) q[0];
sx q[0];
rz(-0.16955626) q[0];
sx q[0];
rz(-2.3441699) q[0];
rz(1.5358309) q[2];
sx q[2];
rz(-0.87858534) q[2];
sx q[2];
rz(0.12223211) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1187073) q[1];
sx q[1];
rz(-1.361025) q[1];
sx q[1];
rz(-0.24521141) q[1];
x q[2];
rz(-0.14446958) q[3];
sx q[3];
rz(-1.472755) q[3];
sx q[3];
rz(-1.4812183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9609191) q[2];
sx q[2];
rz(-2.2283165) q[2];
sx q[2];
rz(1.712435) q[2];
rz(-0.6116496) q[3];
sx q[3];
rz(-1.4432171) q[3];
sx q[3];
rz(-0.071368607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.3058158) q[0];
sx q[0];
rz(-2.340305) q[0];
sx q[0];
rz(0.74547705) q[0];
rz(1.2019134) q[1];
sx q[1];
rz(-1.8169836) q[1];
sx q[1];
rz(2.1048996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3681524) q[0];
sx q[0];
rz(-0.37640171) q[0];
sx q[0];
rz(-1.8252116) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3747524) q[2];
sx q[2];
rz(-1.5897703) q[2];
sx q[2];
rz(-0.60650596) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.503232) q[1];
sx q[1];
rz(-0.88729862) q[1];
sx q[1];
rz(1.2018426) q[1];
rz(2.8193982) q[3];
sx q[3];
rz(-2.5327842) q[3];
sx q[3];
rz(2.6170066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1408954) q[2];
sx q[2];
rz(-1.8706198) q[2];
sx q[2];
rz(-0.97266436) q[2];
rz(-2.2843212) q[3];
sx q[3];
rz(-1.0664777) q[3];
sx q[3];
rz(-1.2681786) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.641441) q[0];
sx q[0];
rz(-0.86857906) q[0];
sx q[0];
rz(-0.80297536) q[0];
rz(2.4909486) q[1];
sx q[1];
rz(-2.3425808) q[1];
sx q[1];
rz(-0.21779901) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2265711) q[0];
sx q[0];
rz(-0.75998291) q[0];
sx q[0];
rz(1.0081916) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3399383) q[2];
sx q[2];
rz(-2.1862324) q[2];
sx q[2];
rz(0.99986693) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.21325363) q[1];
sx q[1];
rz(-1.4782741) q[1];
sx q[1];
rz(1.0800581) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59813114) q[3];
sx q[3];
rz(-2.7418828) q[3];
sx q[3];
rz(2.3593115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.190072) q[2];
sx q[2];
rz(-1.051544) q[2];
sx q[2];
rz(1.5104843) q[2];
rz(1.914628) q[3];
sx q[3];
rz(-1.0534143) q[3];
sx q[3];
rz(2.1372883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25201061) q[0];
sx q[0];
rz(-2.4376106) q[0];
sx q[0];
rz(-0.66396436) q[0];
rz(-0.12531677) q[1];
sx q[1];
rz(-1.7071416) q[1];
sx q[1];
rz(1.0391611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1215724) q[0];
sx q[0];
rz(-1.5680321) q[0];
sx q[0];
rz(-0.0098258709) q[0];
rz(2.8739615) q[2];
sx q[2];
rz(-1.1881141) q[2];
sx q[2];
rz(-3.1185993) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8092371) q[1];
sx q[1];
rz(-1.0910104) q[1];
sx q[1];
rz(-2.9055041) q[1];
rz(-pi) q[2];
rz(-0.44289808) q[3];
sx q[3];
rz(-1.4848466) q[3];
sx q[3];
rz(1.7922873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30109721) q[2];
sx q[2];
rz(-1.4157462) q[2];
sx q[2];
rz(-3.0636129) q[2];
rz(1.1178499) q[3];
sx q[3];
rz(-1.0406787) q[3];
sx q[3];
rz(2.1054721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9341105) q[0];
sx q[0];
rz(-2.9353862) q[0];
sx q[0];
rz(-2.7284486) q[0];
rz(-1.9059937) q[1];
sx q[1];
rz(-1.3256336) q[1];
sx q[1];
rz(1.3135757) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.006047) q[0];
sx q[0];
rz(-0.94903799) q[0];
sx q[0];
rz(-0.22613392) q[0];
rz(-pi) q[1];
rz(-2.2143873) q[2];
sx q[2];
rz(-2.5816133) q[2];
sx q[2];
rz(0.72712979) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47258618) q[1];
sx q[1];
rz(-0.60429497) q[1];
sx q[1];
rz(0.86493203) q[1];
rz(-pi) q[2];
rz(2.6875273) q[3];
sx q[3];
rz(-0.84058529) q[3];
sx q[3];
rz(2.8968352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95973394) q[2];
sx q[2];
rz(-1.9071969) q[2];
sx q[2];
rz(2.6971297) q[2];
rz(-1.2005165) q[3];
sx q[3];
rz(-1.8903172) q[3];
sx q[3];
rz(-3.0357231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926517) q[0];
sx q[0];
rz(-2.7466725) q[0];
sx q[0];
rz(-0.077202395) q[0];
rz(-0.96241799) q[1];
sx q[1];
rz(-1.187477) q[1];
sx q[1];
rz(3.107792) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156436) q[0];
sx q[0];
rz(-1.2275877) q[0];
sx q[0];
rz(2.6604466) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.097418) q[2];
sx q[2];
rz(-0.54810537) q[2];
sx q[2];
rz(-1.4352611) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.56735605) q[1];
sx q[1];
rz(-2.1891382) q[1];
sx q[1];
rz(-2.7537557) q[1];
rz(-pi) q[2];
rz(0.77671364) q[3];
sx q[3];
rz(-1.6459822) q[3];
sx q[3];
rz(3.1326339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9683711) q[2];
sx q[2];
rz(-1.6049478) q[2];
sx q[2];
rz(0.25924337) q[2];
rz(0.94414583) q[3];
sx q[3];
rz(-2.9278432) q[3];
sx q[3];
rz(1.3313782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060870085) q[0];
sx q[0];
rz(-0.20651564) q[0];
sx q[0];
rz(0.61332214) q[0];
rz(2.2382286) q[1];
sx q[1];
rz(-0.84066835) q[1];
sx q[1];
rz(2.9936252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4658303) q[0];
sx q[0];
rz(-2.0480541) q[0];
sx q[0];
rz(2.3289519) q[0];
rz(3.0736835) q[2];
sx q[2];
rz(-2.7017044) q[2];
sx q[2];
rz(-1.875017) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1175673) q[1];
sx q[1];
rz(-0.88741517) q[1];
sx q[1];
rz(-2.1272117) q[1];
x q[2];
rz(2.1359753) q[3];
sx q[3];
rz(-1.5907173) q[3];
sx q[3];
rz(-1.3764149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2030187) q[2];
sx q[2];
rz(-1.729894) q[2];
sx q[2];
rz(-2.1417248) q[2];
rz(-1.3153007) q[3];
sx q[3];
rz(-0.2838997) q[3];
sx q[3];
rz(2.4711173) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1891747) q[0];
sx q[0];
rz(-2.2947831) q[0];
sx q[0];
rz(-2.1536105) q[0];
rz(1.8723764) q[1];
sx q[1];
rz(-2.3321584) q[1];
sx q[1];
rz(-0.16407897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7994499) q[0];
sx q[0];
rz(-2.9963494) q[0];
sx q[0];
rz(2.5119315) q[0];
rz(-pi) q[1];
rz(2.0664178) q[2];
sx q[2];
rz(-1.7640503) q[2];
sx q[2];
rz(2.9444864) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8089701) q[1];
sx q[1];
rz(-1.5738942) q[1];
sx q[1];
rz(-1.078152) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3173728) q[3];
sx q[3];
rz(-2.4755423) q[3];
sx q[3];
rz(-1.0354774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24171955) q[2];
sx q[2];
rz(-2.6417929) q[2];
sx q[2];
rz(-2.0797753) q[2];
rz(0.1344943) q[3];
sx q[3];
rz(-0.89956346) q[3];
sx q[3];
rz(0.053248052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4835085) q[0];
sx q[0];
rz(-0.71165076) q[0];
sx q[0];
rz(2.9363976) q[0];
rz(2.9070053) q[1];
sx q[1];
rz(-1.7154452) q[1];
sx q[1];
rz(0.1964143) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.043606) q[0];
sx q[0];
rz(-0.65904407) q[0];
sx q[0];
rz(1.0511257) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8750849) q[2];
sx q[2];
rz(-0.54000137) q[2];
sx q[2];
rz(2.2842479) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.29039055) q[1];
sx q[1];
rz(-2.7702077) q[1];
sx q[1];
rz(-2.0773649) q[1];
rz(2.100575) q[3];
sx q[3];
rz(-2.2781721) q[3];
sx q[3];
rz(-2.012697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74497574) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(1.9632001) q[2];
rz(-0.34934238) q[3];
sx q[3];
rz(-2.0125466) q[3];
sx q[3];
rz(1.0497302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1402533) q[0];
sx q[0];
rz(-2.0297191) q[0];
sx q[0];
rz(-2.4358791) q[0];
rz(-2.4440675) q[1];
sx q[1];
rz(-1.3730647) q[1];
sx q[1];
rz(-0.90550214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4587935) q[0];
sx q[0];
rz(-2.0044649) q[0];
sx q[0];
rz(-1.5521897) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36351136) q[2];
sx q[2];
rz(-1.8773407) q[2];
sx q[2];
rz(-0.53378045) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9230629) q[1];
sx q[1];
rz(-0.31421146) q[1];
sx q[1];
rz(2.988222) q[1];
rz(-pi) q[2];
rz(-2.3042514) q[3];
sx q[3];
rz(-0.47347927) q[3];
sx q[3];
rz(0.66094962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2424348) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(2.4230797) q[2];
rz(0.68273035) q[3];
sx q[3];
rz(-1.1444789) q[3];
sx q[3];
rz(-3.0375286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8231507) q[0];
sx q[0];
rz(-0.57118509) q[0];
sx q[0];
rz(-0.1027064) q[0];
rz(-1.6406583) q[1];
sx q[1];
rz(-1.2702912) q[1];
sx q[1];
rz(2.695695) q[1];
rz(0.14328843) q[2];
sx q[2];
rz(-1.454151) q[2];
sx q[2];
rz(3.1232338) q[2];
rz(2.3376158) q[3];
sx q[3];
rz(-2.0686276) q[3];
sx q[3];
rz(-0.008240464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
