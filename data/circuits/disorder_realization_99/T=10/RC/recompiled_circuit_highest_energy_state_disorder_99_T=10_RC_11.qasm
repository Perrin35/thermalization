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
rz(-0.55627745) q[0];
sx q[0];
rz(3.102432) q[0];
sx q[0];
rz(10.089212) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(-1.4540949) q[1];
sx q[1];
rz(1.7323642) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7519386) q[0];
sx q[0];
rz(-1.03878) q[0];
sx q[0];
rz(1.7641032) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6970942) q[2];
sx q[2];
rz(-1.2685565) q[2];
sx q[2];
rz(-1.5950749) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4709027) q[1];
sx q[1];
rz(-1.8876612) q[1];
sx q[1];
rz(-0.58611996) q[1];
x q[2];
rz(-0.18376155) q[3];
sx q[3];
rz(-1.9270883) q[3];
sx q[3];
rz(2.3728865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0802143) q[2];
sx q[2];
rz(-2.5900216) q[2];
sx q[2];
rz(2.3870094) q[2];
rz(1.6825698) q[3];
sx q[3];
rz(-2.2857234) q[3];
sx q[3];
rz(-1.4065019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.215312) q[0];
sx q[0];
rz(-0.54406852) q[0];
sx q[0];
rz(1.1625483) q[0];
rz(2.4303719) q[1];
sx q[1];
rz(-0.7178719) q[1];
sx q[1];
rz(0.1444764) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17381771) q[0];
sx q[0];
rz(-1.4878232) q[0];
sx q[0];
rz(-2.847958) q[0];
rz(-pi) q[1];
rz(-0.87717339) q[2];
sx q[2];
rz(-0.85169221) q[2];
sx q[2];
rz(2.4404877) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4481758) q[1];
sx q[1];
rz(-0.87478335) q[1];
sx q[1];
rz(-2.059518) q[1];
x q[2];
rz(2.8322093) q[3];
sx q[3];
rz(-2.2094634) q[3];
sx q[3];
rz(-2.2974599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43313906) q[2];
sx q[2];
rz(-1.3044367) q[2];
sx q[2];
rz(-2.8348095) q[2];
rz(0.77543801) q[3];
sx q[3];
rz(-0.38475761) q[3];
sx q[3];
rz(0.14498372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6590092) q[0];
sx q[0];
rz(-0.44335303) q[0];
sx q[0];
rz(0.79202598) q[0];
rz(-1.5576942) q[1];
sx q[1];
rz(-1.3521103) q[1];
sx q[1];
rz(0.99383324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4531012) q[0];
sx q[0];
rz(-1.9217886) q[0];
sx q[0];
rz(0.61601244) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9609072) q[2];
sx q[2];
rz(-2.5441818) q[2];
sx q[2];
rz(2.1512845) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5573008) q[1];
sx q[1];
rz(-1.2879432) q[1];
sx q[1];
rz(2.0209347) q[1];
x q[2];
rz(-1.1937856) q[3];
sx q[3];
rz(-0.48156958) q[3];
sx q[3];
rz(1.4762312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9516248) q[2];
sx q[2];
rz(-0.8364532) q[2];
sx q[2];
rz(2.9659029) q[2];
rz(-2.9133993) q[3];
sx q[3];
rz(-2.5722645) q[3];
sx q[3];
rz(-0.5051676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93155414) q[0];
sx q[0];
rz(-2.2966972) q[0];
sx q[0];
rz(1.5616052) q[0];
rz(-2.2700894) q[1];
sx q[1];
rz(-1.5305287) q[1];
sx q[1];
rz(-0.16564381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4145395) q[0];
sx q[0];
rz(-3.139608) q[0];
sx q[0];
rz(-0.63568799) q[0];
rz(2.0164967) q[2];
sx q[2];
rz(-1.5786849) q[2];
sx q[2];
rz(0.6296905) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1936744) q[1];
sx q[1];
rz(-2.3478824) q[1];
sx q[1];
rz(-3.0315184) q[1];
rz(-pi) q[2];
rz(-2.7321841) q[3];
sx q[3];
rz(-0.62088517) q[3];
sx q[3];
rz(-1.3001022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6988301) q[2];
sx q[2];
rz(-2.5593968) q[2];
sx q[2];
rz(-2.3637135) q[2];
rz(-2.2569518) q[3];
sx q[3];
rz(-1.9886465) q[3];
sx q[3];
rz(-2.1625429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1409461) q[0];
sx q[0];
rz(-2.1053173) q[0];
sx q[0];
rz(-2.272814) q[0];
rz(-0.90256214) q[1];
sx q[1];
rz(-1.2259918) q[1];
sx q[1];
rz(-0.57194078) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7392699) q[0];
sx q[0];
rz(-0.86193854) q[0];
sx q[0];
rz(-2.2174382) q[0];
rz(-2.6380013) q[2];
sx q[2];
rz(-1.9707754) q[2];
sx q[2];
rz(1.3445889) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8780788) q[1];
sx q[1];
rz(-2.0058419) q[1];
sx q[1];
rz(0.33532354) q[1];
rz(1.0056313) q[3];
sx q[3];
rz(-2.701557) q[3];
sx q[3];
rz(-0.43340757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9945485) q[2];
sx q[2];
rz(-2.2037555) q[2];
sx q[2];
rz(2.5082972) q[2];
rz(-0.58081943) q[3];
sx q[3];
rz(-1.7827026) q[3];
sx q[3];
rz(2.4766428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13034114) q[0];
sx q[0];
rz(-0.3891775) q[0];
sx q[0];
rz(-1.8036386) q[0];
rz(1.7300216) q[1];
sx q[1];
rz(-2.7346225) q[1];
sx q[1];
rz(-0.79708797) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4914843) q[0];
sx q[0];
rz(-2.8239125) q[0];
sx q[0];
rz(-2.0965212) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0051959) q[2];
sx q[2];
rz(-1.505064) q[2];
sx q[2];
rz(2.1594723) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.460372) q[1];
sx q[1];
rz(-2.8275194) q[1];
sx q[1];
rz(1.0286147) q[1];
rz(-pi) q[2];
x q[2];
rz(2.09054) q[3];
sx q[3];
rz(-1.4134348) q[3];
sx q[3];
rz(1.5417772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36593124) q[2];
sx q[2];
rz(-0.66880995) q[2];
sx q[2];
rz(-2.456341) q[2];
rz(0.25964409) q[3];
sx q[3];
rz(-2.1681867) q[3];
sx q[3];
rz(-1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.892136) q[0];
sx q[0];
rz(-0.15661713) q[0];
sx q[0];
rz(-0.53949612) q[0];
rz(0.19142137) q[1];
sx q[1];
rz(-1.6554662) q[1];
sx q[1];
rz(-2.6991381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6689597) q[0];
sx q[0];
rz(-1.5704535) q[0];
sx q[0];
rz(1.5799205) q[0];
rz(-pi) q[1];
rz(-1.0422969) q[2];
sx q[2];
rz(-0.69487725) q[2];
sx q[2];
rz(2.881584) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4686582) q[1];
sx q[1];
rz(-2.4272222) q[1];
sx q[1];
rz(2.0972392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78667647) q[3];
sx q[3];
rz(-1.0287675) q[3];
sx q[3];
rz(-1.4876798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29360867) q[2];
sx q[2];
rz(-2.831735) q[2];
sx q[2];
rz(-2.4388745) q[2];
rz(0.27788776) q[3];
sx q[3];
rz(-1.0669471) q[3];
sx q[3];
rz(1.3387298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24669312) q[0];
sx q[0];
rz(-2.2791857) q[0];
sx q[0];
rz(-0.53584164) q[0];
rz(2.4193173) q[1];
sx q[1];
rz(-1.7012137) q[1];
sx q[1];
rz(-2.3586418) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0412126) q[0];
sx q[0];
rz(-1.6579275) q[0];
sx q[0];
rz(1.210733) q[0];
rz(-pi) q[1];
rz(-0.13320226) q[2];
sx q[2];
rz(-0.40051546) q[2];
sx q[2];
rz(-2.7513964) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6668876) q[1];
sx q[1];
rz(-2.0829446) q[1];
sx q[1];
rz(0.94774232) q[1];
rz(-pi) q[2];
rz(2.7045428) q[3];
sx q[3];
rz(-2.0756222) q[3];
sx q[3];
rz(0.35660353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7213514) q[2];
sx q[2];
rz(-2.6768117) q[2];
sx q[2];
rz(0.50344023) q[2];
rz(-0.33468801) q[3];
sx q[3];
rz(-1.3379593) q[3];
sx q[3];
rz(-2.9065342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47401416) q[0];
sx q[0];
rz(-2.773371) q[0];
sx q[0];
rz(-2.0817122) q[0];
rz(2.8267951) q[1];
sx q[1];
rz(-1.7960725) q[1];
sx q[1];
rz(1.559929) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81875694) q[0];
sx q[0];
rz(-2.4222932) q[0];
sx q[0];
rz(0.34157217) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8354206) q[2];
sx q[2];
rz(-1.6034063) q[2];
sx q[2];
rz(0.6219686) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25645721) q[1];
sx q[1];
rz(-2.4607435) q[1];
sx q[1];
rz(-1.9155986) q[1];
rz(-3.1195333) q[3];
sx q[3];
rz(-1.9859055) q[3];
sx q[3];
rz(-0.65655113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9554837) q[2];
sx q[2];
rz(-1.1649705) q[2];
sx q[2];
rz(2.1060409) q[2];
rz(1.995685) q[3];
sx q[3];
rz(-2.0290387) q[3];
sx q[3];
rz(0.15315332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.938852) q[0];
sx q[0];
rz(-3.0227612) q[0];
sx q[0];
rz(2.3760997) q[0];
rz(1.0368404) q[1];
sx q[1];
rz(-1.8427269) q[1];
sx q[1];
rz(0.11229215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2846887) q[0];
sx q[0];
rz(-0.56988003) q[0];
sx q[0];
rz(-2.5068552) q[0];
x q[1];
rz(-0.77941497) q[2];
sx q[2];
rz(-1.1128491) q[2];
sx q[2];
rz(-0.13546695) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.3272373) q[1];
sx q[1];
rz(-0.8443409) q[1];
sx q[1];
rz(-2.8721456) q[1];
rz(-1.5239612) q[3];
sx q[3];
rz(-0.78748393) q[3];
sx q[3];
rz(1.83444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3969193) q[2];
sx q[2];
rz(-1.8645218) q[2];
sx q[2];
rz(-3.0199158) q[2];
rz(-2.7820898) q[3];
sx q[3];
rz(-2.4183141) q[3];
sx q[3];
rz(0.017688964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.551238) q[0];
sx q[0];
rz(-1.4689162) q[0];
sx q[0];
rz(1.3015251) q[0];
rz(-1.3941258) q[1];
sx q[1];
rz(-1.9448517) q[1];
sx q[1];
rz(1.3762884) q[1];
rz(-0.57264282) q[2];
sx q[2];
rz(-2.4349458) q[2];
sx q[2];
rz(-1.4473421) q[2];
rz(1.2944503) q[3];
sx q[3];
rz(-2.69097) q[3];
sx q[3];
rz(0.16926058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
