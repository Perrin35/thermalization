OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2990155) q[0];
sx q[0];
rz(-2.9691073) q[0];
sx q[0];
rz(-0.01699288) q[0];
rz(0.51141557) q[1];
sx q[1];
rz(-0.49632448) q[1];
sx q[1];
rz(0.48547784) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4776909) q[0];
sx q[0];
rz(-1.4582275) q[0];
sx q[0];
rz(-1.4590053) q[0];
x q[1];
rz(0.67022143) q[2];
sx q[2];
rz(-0.97381684) q[2];
sx q[2];
rz(0.94979294) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6852666) q[1];
sx q[1];
rz(-2.2783372) q[1];
sx q[1];
rz(-0.68143845) q[1];
rz(1.9157682) q[3];
sx q[3];
rz(-2.8135423) q[3];
sx q[3];
rz(1.6102745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7444676) q[2];
sx q[2];
rz(-0.10245377) q[2];
sx q[2];
rz(-0.79258072) q[2];
rz(0.98627311) q[3];
sx q[3];
rz(-1.6542185) q[3];
sx q[3];
rz(-3.0084394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2546286) q[0];
sx q[0];
rz(-1.6034842) q[0];
sx q[0];
rz(-0.1057374) q[0];
rz(2.3836783) q[1];
sx q[1];
rz(-0.71280232) q[1];
sx q[1];
rz(-0.47592083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83748802) q[0];
sx q[0];
rz(-1.6679576) q[0];
sx q[0];
rz(-0.23529737) q[0];
rz(-pi) q[1];
rz(-1.4272825) q[2];
sx q[2];
rz(-1.0821618) q[2];
sx q[2];
rz(-1.8861063) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30224274) q[1];
sx q[1];
rz(-2.1363597) q[1];
sx q[1];
rz(2.9981705) q[1];
rz(-2.2146899) q[3];
sx q[3];
rz(-1.8103278) q[3];
sx q[3];
rz(2.2780971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14625749) q[2];
sx q[2];
rz(-2.6671851) q[2];
sx q[2];
rz(1.3169301) q[2];
rz(0.82777348) q[3];
sx q[3];
rz(-1.0931949) q[3];
sx q[3];
rz(-0.83724666) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4227609) q[0];
sx q[0];
rz(-1.7976924) q[0];
sx q[0];
rz(-1.9047009) q[0];
rz(1.3924567) q[1];
sx q[1];
rz(-1.4207276) q[1];
sx q[1];
rz(-0.69033355) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.515432) q[0];
sx q[0];
rz(-1.0179839) q[0];
sx q[0];
rz(2.1129069) q[0];
rz(-2.0341145) q[2];
sx q[2];
rz(-1.5672188) q[2];
sx q[2];
rz(-1.6301375) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.54385932) q[1];
sx q[1];
rz(-1.3784474) q[1];
sx q[1];
rz(-1.7406157) q[1];
rz(-pi) q[2];
rz(-2.407904) q[3];
sx q[3];
rz(-1.0881967) q[3];
sx q[3];
rz(-0.35728282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6231125) q[2];
sx q[2];
rz(-1.6192351) q[2];
sx q[2];
rz(-3.0677262) q[2];
rz(-1.1582003) q[3];
sx q[3];
rz(-1.2169714) q[3];
sx q[3];
rz(1.4069675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.912792) q[0];
sx q[0];
rz(-0.055963628) q[0];
sx q[0];
rz(2.9486616) q[0];
rz(1.8979161) q[1];
sx q[1];
rz(-1.7453777) q[1];
sx q[1];
rz(-2.9085433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6911963) q[0];
sx q[0];
rz(-1.6653227) q[0];
sx q[0];
rz(-2.9942542) q[0];
x q[1];
rz(2.2795882) q[2];
sx q[2];
rz(-1.5420015) q[2];
sx q[2];
rz(-0.51143247) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8256047) q[1];
sx q[1];
rz(-1.5145166) q[1];
sx q[1];
rz(-2.83063) q[1];
rz(2.190381) q[3];
sx q[3];
rz(-1.9635634) q[3];
sx q[3];
rz(2.5529566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9243246) q[2];
sx q[2];
rz(-1.3793719) q[2];
sx q[2];
rz(3.0775089) q[2];
rz(-0.25121769) q[3];
sx q[3];
rz(-2.678674) q[3];
sx q[3];
rz(-0.29138756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-0.083704405) q[0];
sx q[0];
rz(-1.8479481) q[0];
sx q[0];
rz(1.2842913) q[0];
rz(-0.0074370782) q[1];
sx q[1];
rz(-1.9556655) q[1];
sx q[1];
rz(-2.1844905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.590132) q[0];
sx q[0];
rz(-2.5883985) q[0];
sx q[0];
rz(2.2070918) q[0];
x q[1];
rz(-0.61422698) q[2];
sx q[2];
rz(-1.059327) q[2];
sx q[2];
rz(1.5661256) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0478555) q[1];
sx q[1];
rz(-1.4720494) q[1];
sx q[1];
rz(-2.8119823) q[1];
rz(-pi) q[2];
rz(-1.7577596) q[3];
sx q[3];
rz(-2.2174944) q[3];
sx q[3];
rz(-1.9892429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1288422) q[2];
sx q[2];
rz(-2.4350171) q[2];
sx q[2];
rz(2.1112704) q[2];
rz(-2.4230867) q[3];
sx q[3];
rz(-2.0593144) q[3];
sx q[3];
rz(2.4350186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4265823) q[0];
sx q[0];
rz(-0.74786818) q[0];
sx q[0];
rz(2.7744875) q[0];
rz(-2.7857065) q[1];
sx q[1];
rz(-2.0242736) q[1];
sx q[1];
rz(-1.5497367) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0192896) q[0];
sx q[0];
rz(-2.3117024) q[0];
sx q[0];
rz(-1.6229963) q[0];
rz(0.11271498) q[2];
sx q[2];
rz(-1.1364971) q[2];
sx q[2];
rz(0.096113228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.15799668) q[1];
sx q[1];
rz(-1.3648197) q[1];
sx q[1];
rz(-2.7218444) q[1];
rz(0.47028413) q[3];
sx q[3];
rz(-1.6374644) q[3];
sx q[3];
rz(1.0460514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1648569) q[2];
sx q[2];
rz(-1.2931436) q[2];
sx q[2];
rz(-3.1381651) q[2];
rz(0.88611832) q[3];
sx q[3];
rz(-1.0930073) q[3];
sx q[3];
rz(-1.6974983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43585983) q[0];
sx q[0];
rz(-1.6678565) q[0];
sx q[0];
rz(-2.4263897) q[0];
rz(-0.12229478) q[1];
sx q[1];
rz(-1.0194174) q[1];
sx q[1];
rz(-2.9939647) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8953596) q[0];
sx q[0];
rz(-0.48070714) q[0];
sx q[0];
rz(2.3395545) q[0];
rz(-pi) q[1];
rz(0.72509693) q[2];
sx q[2];
rz(-2.5864961) q[2];
sx q[2];
rz(1.8738418) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.25179204) q[1];
sx q[1];
rz(-1.3887059) q[1];
sx q[1];
rz(0.86343335) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8054923) q[3];
sx q[3];
rz(-0.2944335) q[3];
sx q[3];
rz(-0.89460835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3132396) q[2];
sx q[2];
rz(-2.8326663) q[2];
sx q[2];
rz(0.1114791) q[2];
rz(0.76505032) q[3];
sx q[3];
rz(-1.7181516) q[3];
sx q[3];
rz(3.1077207) q[3];
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
rz(2.2400804) q[0];
sx q[0];
rz(-0.96918786) q[0];
sx q[0];
rz(1.0028268) q[0];
rz(-2.3560933) q[1];
sx q[1];
rz(-1.5248884) q[1];
sx q[1];
rz(-1.2202107) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5883137) q[0];
sx q[0];
rz(-1.0568585) q[0];
sx q[0];
rz(1.5470424) q[0];
rz(-pi) q[1];
x q[1];
rz(0.062510291) q[2];
sx q[2];
rz(-1.4803807) q[2];
sx q[2];
rz(0.16363283) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4017851) q[1];
sx q[1];
rz(-2.5307753) q[1];
sx q[1];
rz(0.99847762) q[1];
rz(0.16943361) q[3];
sx q[3];
rz(-2.1994576) q[3];
sx q[3];
rz(-0.43600142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9057374) q[2];
sx q[2];
rz(-1.7451655) q[2];
sx q[2];
rz(0.033795707) q[2];
rz(-2.3679521) q[3];
sx q[3];
rz(-1.6126817) q[3];
sx q[3];
rz(2.0788367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3646024) q[0];
sx q[0];
rz(-0.62043014) q[0];
sx q[0];
rz(0.34238368) q[0];
rz(1.7204334) q[1];
sx q[1];
rz(-2.3183289) q[1];
sx q[1];
rz(1.294543) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50534821) q[0];
sx q[0];
rz(-2.7490135) q[0];
sx q[0];
rz(-2.7468203) q[0];
x q[1];
rz(-2.2850288) q[2];
sx q[2];
rz(-0.98443778) q[2];
sx q[2];
rz(-2.7991653) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.61664904) q[1];
sx q[1];
rz(-1.6934062) q[1];
sx q[1];
rz(-0.7864237) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3986843) q[3];
sx q[3];
rz(-1.9034608) q[3];
sx q[3];
rz(1.9486959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2932059) q[2];
sx q[2];
rz(-1.7561971) q[2];
sx q[2];
rz(-0.53655857) q[2];
rz(-0.41283354) q[3];
sx q[3];
rz(-0.53240132) q[3];
sx q[3];
rz(-1.1134953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8263016) q[0];
sx q[0];
rz(-1.2696126) q[0];
sx q[0];
rz(1.6444561) q[0];
rz(2.3313088) q[1];
sx q[1];
rz(-1.5850001) q[1];
sx q[1];
rz(0.84698814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4522471) q[0];
sx q[0];
rz(-0.23437491) q[0];
sx q[0];
rz(2.330392) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2184847) q[2];
sx q[2];
rz(-1.1121684) q[2];
sx q[2];
rz(-3.0748526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9223122) q[1];
sx q[1];
rz(-2.4488291) q[1];
sx q[1];
rz(0.42968603) q[1];
rz(-1.8448006) q[3];
sx q[3];
rz(-1.6633667) q[3];
sx q[3];
rz(-2.0333729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0392796) q[2];
sx q[2];
rz(-1.8259093) q[2];
sx q[2];
rz(0.63423356) q[2];
rz(2.5041194) q[3];
sx q[3];
rz(-1.1280779) q[3];
sx q[3];
rz(2.9243961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1524326) q[0];
sx q[0];
rz(-1.1897054) q[0];
sx q[0];
rz(-1.1783896) q[0];
rz(2.4329026) q[1];
sx q[1];
rz(-1.7955753) q[1];
sx q[1];
rz(1.2830455) q[1];
rz(0.37527714) q[2];
sx q[2];
rz(-1.2416544) q[2];
sx q[2];
rz(-1.9273368) q[2];
rz(-1.1905963) q[3];
sx q[3];
rz(-2.1722542) q[3];
sx q[3];
rz(0.62361591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
