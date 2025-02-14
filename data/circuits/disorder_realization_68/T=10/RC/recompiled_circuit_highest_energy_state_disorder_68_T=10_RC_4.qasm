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
rz(0.23213586) q[0];
sx q[0];
rz(-0.27096662) q[0];
sx q[0];
rz(-1.1722857) q[0];
rz(1.5341893) q[1];
sx q[1];
rz(-1.6521896) q[1];
sx q[1];
rz(0.54816562) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48972691) q[0];
sx q[0];
rz(-1.681156) q[0];
sx q[0];
rz(-1.8462028) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3024523) q[2];
sx q[2];
rz(-2.4436032) q[2];
sx q[2];
rz(0.94095397) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9935659) q[1];
sx q[1];
rz(-1.8746334) q[1];
sx q[1];
rz(0.5263473) q[1];
rz(-1.1211419) q[3];
sx q[3];
rz(-1.5142875) q[3];
sx q[3];
rz(0.60349899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9008909) q[2];
sx q[2];
rz(-1.2071995) q[2];
sx q[2];
rz(-1.5473676) q[2];
rz(2.2089925) q[3];
sx q[3];
rz(-0.91351944) q[3];
sx q[3];
rz(0.80476052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7050183) q[0];
sx q[0];
rz(-0.60282928) q[0];
sx q[0];
rz(2.9233209) q[0];
rz(1.6328579) q[1];
sx q[1];
rz(-0.873133) q[1];
sx q[1];
rz(-1.1223209) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86271226) q[0];
sx q[0];
rz(-2.3718908) q[0];
sx q[0];
rz(-0.93602009) q[0];
rz(-pi) q[1];
rz(-0.74064765) q[2];
sx q[2];
rz(-1.7193931) q[2];
sx q[2];
rz(-0.94287485) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7699204) q[1];
sx q[1];
rz(-0.74410106) q[1];
sx q[1];
rz(2.4456853) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2436879) q[3];
sx q[3];
rz(-2.7460685) q[3];
sx q[3];
rz(2.4506068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1307158) q[2];
sx q[2];
rz(-1.1900095) q[2];
sx q[2];
rz(-2.7575764) q[2];
rz(-2.5859517) q[3];
sx q[3];
rz(-2.5540387) q[3];
sx q[3];
rz(-0.89239341) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71690774) q[0];
sx q[0];
rz(-0.45775828) q[0];
sx q[0];
rz(-2.6440788) q[0];
rz(-2.6566907) q[1];
sx q[1];
rz(-1.5496016) q[1];
sx q[1];
rz(2.7000694) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45831523) q[0];
sx q[0];
rz(-2.2234869) q[0];
sx q[0];
rz(2.9291332) q[0];
x q[1];
rz(1.1782568) q[2];
sx q[2];
rz(-0.52670331) q[2];
sx q[2];
rz(1.8369499) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5680629) q[1];
sx q[1];
rz(-2.2383949) q[1];
sx q[1];
rz(2.1993162) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2526337) q[3];
sx q[3];
rz(-1.3103879) q[3];
sx q[3];
rz(1.7956778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32450822) q[2];
sx q[2];
rz(-2.0980947) q[2];
sx q[2];
rz(-0.95477742) q[2];
rz(3.1214516) q[3];
sx q[3];
rz(-0.79831278) q[3];
sx q[3];
rz(-0.32625833) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37498736) q[0];
sx q[0];
rz(-0.52244455) q[0];
sx q[0];
rz(3.1368384) q[0];
rz(0.44113723) q[1];
sx q[1];
rz(-0.10103592) q[1];
sx q[1];
rz(-1.4099247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.649448) q[0];
sx q[0];
rz(-2.7810367) q[0];
sx q[0];
rz(1.1900405) q[0];
rz(-pi) q[1];
rz(0.84216046) q[2];
sx q[2];
rz(-1.7787063) q[2];
sx q[2];
rz(-2.4198857) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.46466848) q[1];
sx q[1];
rz(-2.0629971) q[1];
sx q[1];
rz(2.8850624) q[1];
x q[2];
rz(0.027214931) q[3];
sx q[3];
rz(-1.7593972) q[3];
sx q[3];
rz(-2.1299429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8498174) q[2];
sx q[2];
rz(-1.3128277) q[2];
sx q[2];
rz(0.69668359) q[2];
rz(-1.6126532) q[3];
sx q[3];
rz(-2.5051675) q[3];
sx q[3];
rz(0.67970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1729537) q[0];
sx q[0];
rz(-2.8346859) q[0];
sx q[0];
rz(-0.90721834) q[0];
rz(3.0061159) q[1];
sx q[1];
rz(-1.2716581) q[1];
sx q[1];
rz(2.4438593) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1013779) q[0];
sx q[0];
rz(-1.4411949) q[0];
sx q[0];
rz(-1.1762397) q[0];
rz(1.3759099) q[2];
sx q[2];
rz(-1.5686252) q[2];
sx q[2];
rz(1.2819829) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7373152) q[1];
sx q[1];
rz(-2.7676146) q[1];
sx q[1];
rz(2.651792) q[1];
rz(-pi) q[2];
rz(2.145153) q[3];
sx q[3];
rz(-3.0868885) q[3];
sx q[3];
rz(-0.32848725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.1873143) q[2];
sx q[2];
rz(-0.66437393) q[2];
sx q[2];
rz(-2.7962255) q[2];
rz(2.6464388) q[3];
sx q[3];
rz(-2.5457355) q[3];
sx q[3];
rz(0.12721795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0536163) q[0];
sx q[0];
rz(-1.6642267) q[0];
sx q[0];
rz(1.1548868) q[0];
rz(2.3467973) q[1];
sx q[1];
rz(-1.2966917) q[1];
sx q[1];
rz(0.48387873) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.799494) q[0];
sx q[0];
rz(-0.80518196) q[0];
sx q[0];
rz(0.99200392) q[0];
rz(-pi) q[1];
rz(0.15911289) q[2];
sx q[2];
rz(-1.0793574) q[2];
sx q[2];
rz(-0.05201498) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.72493785) q[1];
sx q[1];
rz(-2.0339801) q[1];
sx q[1];
rz(-0.2128667) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.163614) q[3];
sx q[3];
rz(-2.1107499) q[3];
sx q[3];
rz(-2.0288717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38857073) q[2];
sx q[2];
rz(-1.9066255) q[2];
sx q[2];
rz(2.0335601) q[2];
rz(-1.5293416) q[3];
sx q[3];
rz(-3.022091) q[3];
sx q[3];
rz(2.748238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.18688467) q[0];
sx q[0];
rz(-1.0857546) q[0];
sx q[0];
rz(0.67436522) q[0];
rz(-2.2652594) q[1];
sx q[1];
rz(-0.74462157) q[1];
sx q[1];
rz(0.032329917) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6276777) q[0];
sx q[0];
rz(-1.6925294) q[0];
sx q[0];
rz(-1.7072149) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.416308) q[2];
sx q[2];
rz(-1.3488028) q[2];
sx q[2];
rz(-2.7374008) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2421884) q[1];
sx q[1];
rz(-0.63540484) q[1];
sx q[1];
rz(-2.5492378) q[1];
rz(-pi) q[2];
rz(-0.19200872) q[3];
sx q[3];
rz(-2.1550094) q[3];
sx q[3];
rz(-1.8580798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9548367) q[2];
sx q[2];
rz(-2.7080471) q[2];
sx q[2];
rz(2.1046861) q[2];
rz(1.2284651) q[3];
sx q[3];
rz(-2.2322673) q[3];
sx q[3];
rz(-0.19074805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4036338) q[0];
sx q[0];
rz(-3.1138595) q[0];
sx q[0];
rz(0.24895689) q[0];
rz(0.20949334) q[1];
sx q[1];
rz(-1.8003576) q[1];
sx q[1];
rz(2.478821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4512206) q[0];
sx q[0];
rz(-1.5589542) q[0];
sx q[0];
rz(1.580834) q[0];
x q[1];
rz(0.49928645) q[2];
sx q[2];
rz(-0.82816507) q[2];
sx q[2];
rz(-1.7940831) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1472305) q[1];
sx q[1];
rz(-2.3636345) q[1];
sx q[1];
rz(1.2176355) q[1];
rz(-2.7484077) q[3];
sx q[3];
rz(-0.69614702) q[3];
sx q[3];
rz(0.37567155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2304307) q[2];
sx q[2];
rz(-1.5713567) q[2];
sx q[2];
rz(-0.39361185) q[2];
rz(0.39129928) q[3];
sx q[3];
rz(-0.53520447) q[3];
sx q[3];
rz(-0.75214255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0126208) q[0];
sx q[0];
rz(-0.18652815) q[0];
sx q[0];
rz(-0.57998002) q[0];
rz(-0.36803666) q[1];
sx q[1];
rz(-0.8152222) q[1];
sx q[1];
rz(3.0267402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84282473) q[0];
sx q[0];
rz(-1.4359763) q[0];
sx q[0];
rz(2.299286) q[0];
x q[1];
rz(-1.8813558) q[2];
sx q[2];
rz(-2.680696) q[2];
sx q[2];
rz(-2.4577934) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29909322) q[1];
sx q[1];
rz(-2.7242492) q[1];
sx q[1];
rz(-0.6562161) q[1];
rz(-pi) q[2];
rz(1.097947) q[3];
sx q[3];
rz(-2.3896296) q[3];
sx q[3];
rz(-1.0238943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4174691) q[2];
sx q[2];
rz(-0.33895156) q[2];
sx q[2];
rz(-0.69619703) q[2];
rz(2.0255069) q[3];
sx q[3];
rz(-2.3510272) q[3];
sx q[3];
rz(-0.74151403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29257947) q[0];
sx q[0];
rz(-2.1722023) q[0];
sx q[0];
rz(-2.868929) q[0];
rz(-2.4478681) q[1];
sx q[1];
rz(-1.1169746) q[1];
sx q[1];
rz(-0.56347096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16167264) q[0];
sx q[0];
rz(-1.663999) q[0];
sx q[0];
rz(0.037089238) q[0];
rz(-pi) q[1];
rz(0.54948893) q[2];
sx q[2];
rz(-0.27590431) q[2];
sx q[2];
rz(2.9018108) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2595926) q[1];
sx q[1];
rz(-1.881885) q[1];
sx q[1];
rz(-3.056219) q[1];
rz(-pi) q[2];
rz(0.70448168) q[3];
sx q[3];
rz(-0.84957963) q[3];
sx q[3];
rz(-2.3303256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9844661) q[2];
sx q[2];
rz(-2.0968585) q[2];
sx q[2];
rz(2.5926479) q[2];
rz(0.99556154) q[3];
sx q[3];
rz(-0.82058161) q[3];
sx q[3];
rz(2.5115749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41848771) q[0];
sx q[0];
rz(-0.99011078) q[0];
sx q[0];
rz(-0.44575442) q[0];
rz(0.86722974) q[1];
sx q[1];
rz(-2.0384616) q[1];
sx q[1];
rz(-1.9834317) q[1];
rz(0.42648496) q[2];
sx q[2];
rz(-2.6250962) q[2];
sx q[2];
rz(0.442183) q[2];
rz(-0.54789644) q[3];
sx q[3];
rz(-1.4461645) q[3];
sx q[3];
rz(-0.89492284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
