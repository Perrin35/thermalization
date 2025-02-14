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
rz(0.046959538) q[0];
sx q[0];
rz(-1.5488012) q[0];
sx q[0];
rz(-1.7472851) q[0];
rz(-0.46051639) q[1];
sx q[1];
rz(-1.8733652) q[1];
sx q[1];
rz(-0.4761129) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5232485) q[0];
sx q[0];
rz(-1.1914413) q[0];
sx q[0];
rz(2.1567287) q[0];
x q[1];
rz(-2.2865795) q[2];
sx q[2];
rz(-0.35021338) q[2];
sx q[2];
rz(0.24392715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.75501498) q[1];
sx q[1];
rz(-1.5881032) q[1];
sx q[1];
rz(2.9933418) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1143161) q[3];
sx q[3];
rz(-1.6347872) q[3];
sx q[3];
rz(2.3673253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8146879) q[2];
sx q[2];
rz(-1.2219595) q[2];
sx q[2];
rz(-1.8587221) q[2];
rz(0.13947105) q[3];
sx q[3];
rz(-1.8255511) q[3];
sx q[3];
rz(2.4563346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71444702) q[0];
sx q[0];
rz(-0.35816631) q[0];
sx q[0];
rz(2.8156679) q[0];
rz(-2.4320995) q[1];
sx q[1];
rz(-0.26115099) q[1];
sx q[1];
rz(-0.66169468) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70764953) q[0];
sx q[0];
rz(-3.1295332) q[0];
sx q[0];
rz(-2.2318186) q[0];
rz(-pi) q[1];
rz(2.6432132) q[2];
sx q[2];
rz(-1.0000668) q[2];
sx q[2];
rz(0.98266593) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0831415) q[1];
sx q[1];
rz(-1.3399235) q[1];
sx q[1];
rz(1.9564232) q[1];
rz(-pi) q[2];
x q[2];
rz(2.157341) q[3];
sx q[3];
rz(-0.97347608) q[3];
sx q[3];
rz(0.95109361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1208531) q[2];
sx q[2];
rz(-1.3384621) q[2];
sx q[2];
rz(0.2600812) q[2];
rz(-0.86380473) q[3];
sx q[3];
rz(-1.443202) q[3];
sx q[3];
rz(-1.1687733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23673713) q[0];
sx q[0];
rz(-0.78176347) q[0];
sx q[0];
rz(2.8535063) q[0];
rz(1.9347363) q[1];
sx q[1];
rz(-1.1129881) q[1];
sx q[1];
rz(0.77817121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.659063) q[0];
sx q[0];
rz(-0.89119833) q[0];
sx q[0];
rz(0.24404714) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6323998) q[2];
sx q[2];
rz(-0.81356102) q[2];
sx q[2];
rz(-2.1547627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6699804) q[1];
sx q[1];
rz(-1.1410574) q[1];
sx q[1];
rz(-0.55008908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0708848) q[3];
sx q[3];
rz(-0.71083655) q[3];
sx q[3];
rz(-2.5612166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8975767) q[2];
sx q[2];
rz(-1.1268758) q[2];
sx q[2];
rz(2.8934532) q[2];
rz(-1.0330307) q[3];
sx q[3];
rz(-1.4890198) q[3];
sx q[3];
rz(1.8831133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7018062) q[0];
sx q[0];
rz(-2.2070364) q[0];
sx q[0];
rz(2.5820861) q[0];
rz(1.4065546) q[1];
sx q[1];
rz(-0.94466698) q[1];
sx q[1];
rz(1.9962126) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6481768) q[0];
sx q[0];
rz(-2.6828565) q[0];
sx q[0];
rz(0.10347314) q[0];
rz(-2.1297867) q[2];
sx q[2];
rz(-0.99746639) q[2];
sx q[2];
rz(0.99775523) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57210797) q[1];
sx q[1];
rz(-2.3037801) q[1];
sx q[1];
rz(1.1694639) q[1];
rz(-pi) q[2];
rz(-1.5942176) q[3];
sx q[3];
rz(-0.56992793) q[3];
sx q[3];
rz(2.3078634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2650602) q[2];
sx q[2];
rz(-1.0141076) q[2];
sx q[2];
rz(1.2476791) q[2];
rz(1.0809336) q[3];
sx q[3];
rz(-1.388988) q[3];
sx q[3];
rz(-1.8814253) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0718229) q[0];
sx q[0];
rz(-0.48957303) q[0];
sx q[0];
rz(0.74482942) q[0];
rz(1.3818332) q[1];
sx q[1];
rz(-2.0227183) q[1];
sx q[1];
rz(-0.098085731) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52652717) q[0];
sx q[0];
rz(-1.9221109) q[0];
sx q[0];
rz(-1.8976346) q[0];
rz(-pi) q[1];
rz(0.44386835) q[2];
sx q[2];
rz(-1.265268) q[2];
sx q[2];
rz(1.6084087) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91589663) q[1];
sx q[1];
rz(-2.4476123) q[1];
sx q[1];
rz(-0.028631239) q[1];
x q[2];
rz(2.0518584) q[3];
sx q[3];
rz(-1.4816545) q[3];
sx q[3];
rz(-1.7952023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44094008) q[2];
sx q[2];
rz(-0.38267371) q[2];
sx q[2];
rz(1.0901701) q[2];
rz(-0.32554659) q[3];
sx q[3];
rz(-1.5109589) q[3];
sx q[3];
rz(0.24698273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
rz(0.22661041) q[0];
sx q[0];
rz(-2.60422) q[0];
sx q[0];
rz(1.9052624) q[0];
rz(0.63713282) q[1];
sx q[1];
rz(-2.5814711) q[1];
sx q[1];
rz(-2.2083652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97605825) q[0];
sx q[0];
rz(-1.9325496) q[0];
sx q[0];
rz(-0.975859) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1652214) q[2];
sx q[2];
rz(-2.6898807) q[2];
sx q[2];
rz(-3.0690985) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.56343397) q[1];
sx q[1];
rz(-1.8616042) q[1];
sx q[1];
rz(-1.6573424) q[1];
x q[2];
rz(2.0863462) q[3];
sx q[3];
rz(-1.7988867) q[3];
sx q[3];
rz(-0.67319621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.25202641) q[2];
sx q[2];
rz(-1.7539975) q[2];
sx q[2];
rz(-1.8760366) q[2];
rz(-1.7172074) q[3];
sx q[3];
rz(-0.5286743) q[3];
sx q[3];
rz(0.85844794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3739361) q[0];
sx q[0];
rz(-2.7601384) q[0];
sx q[0];
rz(0.23442991) q[0];
rz(-2.7081721) q[1];
sx q[1];
rz(-2.0442918) q[1];
sx q[1];
rz(-2.0030599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075468242) q[0];
sx q[0];
rz(-0.99695092) q[0];
sx q[0];
rz(0.85077758) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6509382) q[2];
sx q[2];
rz(-1.5637959) q[2];
sx q[2];
rz(2.2537083) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2039472) q[1];
sx q[1];
rz(-2.0729985) q[1];
sx q[1];
rz(-0.03003386) q[1];
rz(-0.23536162) q[3];
sx q[3];
rz(-0.70066626) q[3];
sx q[3];
rz(1.2149902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4745657) q[2];
sx q[2];
rz(-2.462025) q[2];
sx q[2];
rz(2.1843074) q[2];
rz(-0.75753093) q[3];
sx q[3];
rz(-1.6672971) q[3];
sx q[3];
rz(0.1434513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677419) q[0];
sx q[0];
rz(-2.0389281) q[0];
sx q[0];
rz(-0.17909166) q[0];
rz(2.9415019) q[1];
sx q[1];
rz(-2.4724019) q[1];
sx q[1];
rz(-0.77654138) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9574253) q[0];
sx q[0];
rz(-1.7514075) q[0];
sx q[0];
rz(2.1023048) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29751038) q[2];
sx q[2];
rz(-0.97504598) q[2];
sx q[2];
rz(1.7864625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4458865) q[1];
sx q[1];
rz(-1.4974471) q[1];
sx q[1];
rz(-0.35584764) q[1];
rz(1.7257878) q[3];
sx q[3];
rz(-0.3213045) q[3];
sx q[3];
rz(-2.6329272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17579235) q[2];
sx q[2];
rz(-2.0696023) q[2];
sx q[2];
rz(-1.779186) q[2];
rz(-1.8386748) q[3];
sx q[3];
rz(-1.6630273) q[3];
sx q[3];
rz(-1.038507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4565444) q[0];
sx q[0];
rz(-1.9467204) q[0];
sx q[0];
rz(-3.0273279) q[0];
rz(1.9396797) q[1];
sx q[1];
rz(-0.54113954) q[1];
sx q[1];
rz(-3.0466383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33126011) q[0];
sx q[0];
rz(-2.5833732) q[0];
sx q[0];
rz(1.9086188) q[0];
x q[1];
rz(2.6812115) q[2];
sx q[2];
rz(-0.91074873) q[2];
sx q[2];
rz(2.6766863) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1324117) q[1];
sx q[1];
rz(-2.5887515) q[1];
sx q[1];
rz(0.10709672) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41422959) q[3];
sx q[3];
rz(-1.0981907) q[3];
sx q[3];
rz(-2.7850658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.3007043) q[2];
sx q[2];
rz(-1.3450832) q[2];
sx q[2];
rz(-0.67187205) q[2];
rz(0.68615174) q[3];
sx q[3];
rz(-2.4785564) q[3];
sx q[3];
rz(1.9239976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9753863) q[0];
sx q[0];
rz(-0.19151846) q[0];
sx q[0];
rz(-1.5331049) q[0];
rz(2.8168822) q[1];
sx q[1];
rz(-1.8638116) q[1];
sx q[1];
rz(0.3130354) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4251809) q[0];
sx q[0];
rz(-1.4907752) q[0];
sx q[0];
rz(2.1437313) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3029477) q[2];
sx q[2];
rz(-2.9137879) q[2];
sx q[2];
rz(2.0382263) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5394143) q[1];
sx q[1];
rz(-1.0262607) q[1];
sx q[1];
rz(0.79411749) q[1];
rz(0.5996941) q[3];
sx q[3];
rz(-2.4238677) q[3];
sx q[3];
rz(-2.5120171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93840557) q[2];
sx q[2];
rz(-1.1537735) q[2];
sx q[2];
rz(1.3333092) q[2];
rz(-0.61776727) q[3];
sx q[3];
rz(-0.94527644) q[3];
sx q[3];
rz(-1.1355404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.7137322) q[0];
sx q[0];
rz(-1.8129616) q[0];
sx q[0];
rz(0.34995361) q[0];
rz(-0.89823828) q[1];
sx q[1];
rz(-1.0844834) q[1];
sx q[1];
rz(-0.35379298) q[1];
rz(-0.74358616) q[2];
sx q[2];
rz(-2.0571041) q[2];
sx q[2];
rz(-1.9509289) q[2];
rz(-1.5469503) q[3];
sx q[3];
rz(-0.84115728) q[3];
sx q[3];
rz(0.41681791) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
