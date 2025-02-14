OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0290282) q[0];
sx q[0];
rz(-1.5052786) q[0];
sx q[0];
rz(2.3583052) q[0];
rz(-2.0593491) q[1];
sx q[1];
rz(-1.6545656) q[1];
sx q[1];
rz(2.3496871) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91814268) q[0];
sx q[0];
rz(-2.8200216) q[0];
sx q[0];
rz(-0.62427036) q[0];
rz(-pi) q[1];
rz(-1.7246805) q[2];
sx q[2];
rz(-1.6072011) q[2];
sx q[2];
rz(-2.3652181) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72630771) q[1];
sx q[1];
rz(-1.0276762) q[1];
sx q[1];
rz(-1.607479) q[1];
x q[2];
rz(-2.5654985) q[3];
sx q[3];
rz(-2.2597426) q[3];
sx q[3];
rz(2.9722633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5504494) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(-3.0421416) q[2];
rz(-2.8948696) q[3];
sx q[3];
rz(-1.5244502) q[3];
sx q[3];
rz(-2.7439086) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0196911) q[0];
sx q[0];
rz(-0.97625232) q[0];
sx q[0];
rz(2.2130527) q[0];
rz(1.4506725) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(-2.4931152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.686931) q[0];
sx q[0];
rz(-1.0333916) q[0];
sx q[0];
rz(0.65673687) q[0];
rz(-pi) q[1];
rz(-1.8869867) q[2];
sx q[2];
rz(-1.0266227) q[2];
sx q[2];
rz(-0.7989102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3120756) q[1];
sx q[1];
rz(-1.5863998) q[1];
sx q[1];
rz(-0.14194686) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2821372) q[3];
sx q[3];
rz(-1.6842173) q[3];
sx q[3];
rz(0.78383776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4497946) q[2];
sx q[2];
rz(-0.75642502) q[2];
sx q[2];
rz(2.2890384) q[2];
rz(0.84613386) q[3];
sx q[3];
rz(-1.5087912) q[3];
sx q[3];
rz(-1.030863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344675) q[0];
sx q[0];
rz(-1.9159303) q[0];
sx q[0];
rz(1.0852098) q[0];
rz(0.94404864) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(-1.6403713) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1000993) q[0];
sx q[0];
rz(-1.3736808) q[0];
sx q[0];
rz(1.5310578) q[0];
rz(-0.88884647) q[2];
sx q[2];
rz(-2.8434128) q[2];
sx q[2];
rz(0.8487289) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2861917) q[1];
sx q[1];
rz(-0.58127379) q[1];
sx q[1];
rz(1.752125) q[1];
x q[2];
rz(-2.9839433) q[3];
sx q[3];
rz(-1.7846037) q[3];
sx q[3];
rz(0.42663867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7368855) q[2];
sx q[2];
rz(-0.53386226) q[2];
sx q[2];
rz(-0.85477465) q[2];
rz(0.9084304) q[3];
sx q[3];
rz(-1.5769438) q[3];
sx q[3];
rz(2.0250208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.214355) q[0];
sx q[0];
rz(-0.4011918) q[0];
sx q[0];
rz(1.1965363) q[0];
rz(-1.6664956) q[1];
sx q[1];
rz(-0.42088446) q[1];
sx q[1];
rz(1.9570785) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9036983) q[0];
sx q[0];
rz(-1.8082976) q[0];
sx q[0];
rz(2.4303275) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20303161) q[2];
sx q[2];
rz(-0.76329279) q[2];
sx q[2];
rz(2.610746) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0851497) q[1];
sx q[1];
rz(-2.5825325) q[1];
sx q[1];
rz(-1.1340301) q[1];
rz(-pi) q[2];
rz(-1.3877901) q[3];
sx q[3];
rz(-1.840045) q[3];
sx q[3];
rz(1.8257917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92675942) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(1.2614177) q[2];
rz(-2.7091806) q[3];
sx q[3];
rz(-2.0093446) q[3];
sx q[3];
rz(2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686907) q[0];
sx q[0];
rz(-1.9111159) q[0];
sx q[0];
rz(-0.029408971) q[0];
rz(-2.3853761) q[1];
sx q[1];
rz(-2.5543946) q[1];
sx q[1];
rz(-0.75278935) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9133224) q[0];
sx q[0];
rz(-3.1359657) q[0];
sx q[0];
rz(1.6371284) q[0];
rz(-pi) q[1];
rz(-3.0514206) q[2];
sx q[2];
rz(-1.8379231) q[2];
sx q[2];
rz(-1.7091441) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50548762) q[1];
sx q[1];
rz(-1.4661745) q[1];
sx q[1];
rz(-2.5701447) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54075586) q[3];
sx q[3];
rz(-1.1782559) q[3];
sx q[3];
rz(2.768571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0011562) q[2];
sx q[2];
rz(-1.493528) q[2];
sx q[2];
rz(2.746554) q[2];
rz(-0.47518528) q[3];
sx q[3];
rz(-1.3522215) q[3];
sx q[3];
rz(2.7648259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756556) q[0];
sx q[0];
rz(-0.7203311) q[0];
sx q[0];
rz(0.96250594) q[0];
rz(-2.6267701) q[1];
sx q[1];
rz(-1.6074901) q[1];
sx q[1];
rz(2.8033676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7027959) q[0];
sx q[0];
rz(-1.3745752) q[0];
sx q[0];
rz(0.81721925) q[0];
rz(-2.9043496) q[2];
sx q[2];
rz(-1.8117935) q[2];
sx q[2];
rz(0.70912305) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2185827) q[1];
sx q[1];
rz(-1.5707695) q[1];
sx q[1];
rz(-1.567651) q[1];
x q[2];
rz(1.7231483) q[3];
sx q[3];
rz(-2.7151516) q[3];
sx q[3];
rz(-1.3063198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.50011355) q[2];
sx q[2];
rz(-1.1621472) q[2];
sx q[2];
rz(2.5872453) q[2];
rz(0.85136271) q[3];
sx q[3];
rz(-0.35836372) q[3];
sx q[3];
rz(2.5135777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8165269) q[0];
sx q[0];
rz(-0.61282235) q[0];
sx q[0];
rz(0.75041962) q[0];
rz(2.5665307) q[1];
sx q[1];
rz(-1.3651747) q[1];
sx q[1];
rz(0.65779984) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1019264) q[0];
sx q[0];
rz(-0.79625477) q[0];
sx q[0];
rz(-0.58923652) q[0];
x q[1];
rz(0.051498895) q[2];
sx q[2];
rz(-2.4337075) q[2];
sx q[2];
rz(2.5619363) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9896587) q[1];
sx q[1];
rz(-0.2865782) q[1];
sx q[1];
rz(1.5119988) q[1];
x q[2];
rz(1.1893473) q[3];
sx q[3];
rz(-1.7717517) q[3];
sx q[3];
rz(-2.902346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.647992) q[2];
sx q[2];
rz(-2.1330264) q[2];
sx q[2];
rz(-1.6592525) q[2];
rz(-3.0209387) q[3];
sx q[3];
rz(-1.7574661) q[3];
sx q[3];
rz(0.83546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7928612) q[0];
sx q[0];
rz(-0.86935765) q[0];
sx q[0];
rz(2.8435775) q[0];
rz(1.0385849) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(0.15377741) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34798056) q[0];
sx q[0];
rz(-1.4958994) q[0];
sx q[0];
rz(2.8872284) q[0];
rz(0.80259097) q[2];
sx q[2];
rz(-1.0804515) q[2];
sx q[2];
rz(2.5732694) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3996622) q[1];
sx q[1];
rz(-0.35835086) q[1];
sx q[1];
rz(0.25800024) q[1];
rz(-pi) q[2];
rz(1.8337209) q[3];
sx q[3];
rz(-0.44028966) q[3];
sx q[3];
rz(1.057098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.095857233) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(-1.1886965) q[2];
rz(1.8111604) q[3];
sx q[3];
rz(-1.1878139) q[3];
sx q[3];
rz(-2.4082898) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34340149) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(-1.6424302) q[0];
rz(0.54939735) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(2.8909491) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17535398) q[0];
sx q[0];
rz(-2.4099775) q[0];
sx q[0];
rz(-1.0989972) q[0];
x q[1];
rz(1.556384) q[2];
sx q[2];
rz(-1.8317128) q[2];
sx q[2];
rz(0.018785611) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5329001) q[1];
sx q[1];
rz(-1.6402555) q[1];
sx q[1];
rz(-0.77381247) q[1];
x q[2];
rz(2.9467877) q[3];
sx q[3];
rz(-1.6945632) q[3];
sx q[3];
rz(-2.3969802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.15829076) q[2];
sx q[2];
rz(-3.0088708) q[2];
sx q[2];
rz(2.1251202) q[2];
rz(-3.0554092) q[3];
sx q[3];
rz(-1.0165756) q[3];
sx q[3];
rz(-0.76519722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30855274) q[0];
sx q[0];
rz(-0.85423952) q[0];
sx q[0];
rz(0.41900751) q[0];
rz(0.53681701) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(2.4057665) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9402855) q[0];
sx q[0];
rz(-0.92865151) q[0];
sx q[0];
rz(1.8117732) q[0];
rz(1.0334098) q[2];
sx q[2];
rz(-0.99796406) q[2];
sx q[2];
rz(-0.17102851) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.20504552) q[1];
sx q[1];
rz(-2.1959795) q[1];
sx q[1];
rz(-1.3035098) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9189776) q[3];
sx q[3];
rz(-2.0280119) q[3];
sx q[3];
rz(1.3570076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1824823) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(2.518892) q[2];
rz(-0.73838082) q[3];
sx q[3];
rz(-1.9564956) q[3];
sx q[3];
rz(0.27899376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57711346) q[0];
sx q[0];
rz(-2.8688685) q[0];
sx q[0];
rz(-2.175749) q[0];
rz(0.64361698) q[1];
sx q[1];
rz(-1.8543961) q[1];
sx q[1];
rz(-2.054945) q[1];
rz(1.2014482) q[2];
sx q[2];
rz(-1.3551753) q[2];
sx q[2];
rz(1.7418855) q[2];
rz(1.7122812) q[3];
sx q[3];
rz(-1.332193) q[3];
sx q[3];
rz(-0.79164483) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
