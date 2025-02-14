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
rz(-0.79987502) q[0];
sx q[0];
rz(2.3246111) q[0];
sx q[0];
rz(9.8819879) q[0];
rz(-2.7297821) q[1];
sx q[1];
rz(-1.6219985) q[1];
sx q[1];
rz(-0.55984679) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.354411) q[0];
sx q[0];
rz(-2.5776064) q[0];
sx q[0];
rz(2.5111879) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.067808) q[2];
sx q[2];
rz(-1.8230228) q[2];
sx q[2];
rz(-2.5258738) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.96139446) q[1];
sx q[1];
rz(-1.0751659) q[1];
sx q[1];
rz(-1.2537728) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4234957) q[3];
sx q[3];
rz(-1.3527186) q[3];
sx q[3];
rz(0.17373057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37602279) q[2];
sx q[2];
rz(-2.1799808) q[2];
sx q[2];
rz(1.8454856) q[2];
rz(-0.62093312) q[3];
sx q[3];
rz(-2.4758078) q[3];
sx q[3];
rz(2.2165829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.45944443) q[0];
sx q[0];
rz(-0.58732533) q[0];
sx q[0];
rz(-2.4272954) q[0];
rz(0.722305) q[1];
sx q[1];
rz(-2.0562833) q[1];
sx q[1];
rz(-1.3246271) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52377578) q[0];
sx q[0];
rz(-1.1259698) q[0];
sx q[0];
rz(-0.33828783) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0804964) q[2];
sx q[2];
rz(-2.4078566) q[2];
sx q[2];
rz(-2.3139985) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87323071) q[1];
sx q[1];
rz(-1.2030081) q[1];
sx q[1];
rz(1.7413543) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4217266) q[3];
sx q[3];
rz(-2.4507634) q[3];
sx q[3];
rz(-0.2836424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1076727) q[2];
sx q[2];
rz(-1.718797) q[2];
sx q[2];
rz(-0.99204341) q[2];
rz(2.3658559) q[3];
sx q[3];
rz(-2.3596767) q[3];
sx q[3];
rz(-1.0824664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7139605) q[0];
sx q[0];
rz(-1.7949224) q[0];
sx q[0];
rz(0.91208518) q[0];
rz(2.7569547) q[1];
sx q[1];
rz(-0.96133989) q[1];
sx q[1];
rz(0.49547637) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3867823) q[0];
sx q[0];
rz(-1.6433924) q[0];
sx q[0];
rz(-3.0955546) q[0];
rz(-pi) q[1];
rz(-0.74929535) q[2];
sx q[2];
rz(-1.4297419) q[2];
sx q[2];
rz(-2.2749449) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4810973) q[1];
sx q[1];
rz(-0.7259136) q[1];
sx q[1];
rz(3.0654415) q[1];
rz(-pi) q[2];
rz(2.1927102) q[3];
sx q[3];
rz(-1.5601741) q[3];
sx q[3];
rz(-2.0257906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2428525) q[2];
sx q[2];
rz(-1.1246559) q[2];
sx q[2];
rz(0.43453547) q[2];
rz(3.0430326) q[3];
sx q[3];
rz(-1.7193272) q[3];
sx q[3];
rz(-2.3677473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74715215) q[0];
sx q[0];
rz(-0.23500615) q[0];
sx q[0];
rz(-0.26879841) q[0];
rz(-2.0647743) q[1];
sx q[1];
rz(-1.4067168) q[1];
sx q[1];
rz(1.1928308) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6679316) q[0];
sx q[0];
rz(-0.30465484) q[0];
sx q[0];
rz(0.61078914) q[0];
rz(2.3770726) q[2];
sx q[2];
rz(-2.9540375) q[2];
sx q[2];
rz(-0.16403596) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7320648) q[1];
sx q[1];
rz(-1.3573779) q[1];
sx q[1];
rz(-1.4103389) q[1];
rz(-pi) q[2];
rz(-0.70130879) q[3];
sx q[3];
rz(-1.4150672) q[3];
sx q[3];
rz(-2.6913672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.25527915) q[2];
sx q[2];
rz(-2.2463319) q[2];
sx q[2];
rz(-0.2992343) q[2];
rz(0.29087654) q[3];
sx q[3];
rz(-1.2521005) q[3];
sx q[3];
rz(-0.65557426) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84151477) q[0];
sx q[0];
rz(-2.1617007) q[0];
sx q[0];
rz(0.078068659) q[0];
rz(0.95371753) q[1];
sx q[1];
rz(-2.6225312) q[1];
sx q[1];
rz(-1.5740707) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46569506) q[0];
sx q[0];
rz(-0.90238304) q[0];
sx q[0];
rz(2.9107679) q[0];
rz(2.445053) q[2];
sx q[2];
rz(-1.8395506) q[2];
sx q[2];
rz(-2.0292676) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5213274) q[1];
sx q[1];
rz(-2.5286739) q[1];
sx q[1];
rz(1.6271724) q[1];
x q[2];
rz(-0.084265274) q[3];
sx q[3];
rz(-0.69484303) q[3];
sx q[3];
rz(-2.7488199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2786402) q[2];
sx q[2];
rz(-0.43206698) q[2];
sx q[2];
rz(2.8734109) q[2];
rz(2.6523318) q[3];
sx q[3];
rz(-1.093995) q[3];
sx q[3];
rz(-1.6930273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0972524) q[0];
sx q[0];
rz(-0.02359979) q[0];
sx q[0];
rz(2.6771255) q[0];
rz(-0.52344549) q[1];
sx q[1];
rz(-2.372066) q[1];
sx q[1];
rz(-0.73062599) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3364612) q[0];
sx q[0];
rz(-1.4298273) q[0];
sx q[0];
rz(-0.8223429) q[0];
rz(-pi) q[1];
rz(-2.7774657) q[2];
sx q[2];
rz(-2.4362323) q[2];
sx q[2];
rz(0.10157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1878512) q[1];
sx q[1];
rz(-1.9914728) q[1];
sx q[1];
rz(2.4286146) q[1];
rz(0.2517638) q[3];
sx q[3];
rz(-1.442896) q[3];
sx q[3];
rz(2.3826016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0475433) q[2];
sx q[2];
rz(-2.0619679) q[2];
sx q[2];
rz(-3.003982) q[2];
rz(-1.1329457) q[3];
sx q[3];
rz(-1.1762985) q[3];
sx q[3];
rz(1.2631811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15261821) q[0];
sx q[0];
rz(-1.6292097) q[0];
sx q[0];
rz(-0.033893943) q[0];
rz(-0.35941091) q[1];
sx q[1];
rz(-1.6172599) q[1];
sx q[1];
rz(0.83438897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90859014) q[0];
sx q[0];
rz(-2.6437573) q[0];
sx q[0];
rz(2.4497238) q[0];
x q[1];
rz(2.4530386) q[2];
sx q[2];
rz(-1.9594741) q[2];
sx q[2];
rz(-1.9118903) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.31454489) q[1];
sx q[1];
rz(-2.8706708) q[1];
sx q[1];
rz(-1.2071868) q[1];
rz(-pi) q[2];
rz(-1.5079751) q[3];
sx q[3];
rz(-0.90017002) q[3];
sx q[3];
rz(-0.84298979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72839165) q[2];
sx q[2];
rz(-0.30343702) q[2];
sx q[2];
rz(-2.0835853) q[2];
rz(2.6672065) q[3];
sx q[3];
rz(-2.3460903) q[3];
sx q[3];
rz(-1.0559731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0742842) q[0];
sx q[0];
rz(-1.5165167) q[0];
sx q[0];
rz(-1.3577331) q[0];
rz(-0.43074295) q[1];
sx q[1];
rz(-1.4746702) q[1];
sx q[1];
rz(-2.6745083) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3586836) q[0];
sx q[0];
rz(-3.0984146) q[0];
sx q[0];
rz(1.5367277) q[0];
rz(-1.378304) q[2];
sx q[2];
rz(-0.87710947) q[2];
sx q[2];
rz(2.326593) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.43710762) q[1];
sx q[1];
rz(-0.52981716) q[1];
sx q[1];
rz(-0.036661224) q[1];
x q[2];
rz(-2.2086618) q[3];
sx q[3];
rz(-1.5278421) q[3];
sx q[3];
rz(1.651498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0507386) q[2];
sx q[2];
rz(-2.723912) q[2];
sx q[2];
rz(1.0806855) q[2];
rz(2.4596227) q[3];
sx q[3];
rz(-0.8050279) q[3];
sx q[3];
rz(-0.69847703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4526116) q[0];
sx q[0];
rz(-0.51849759) q[0];
sx q[0];
rz(-2.3338351) q[0];
rz(-2.141433) q[1];
sx q[1];
rz(-2.2554485) q[1];
sx q[1];
rz(-0.79661405) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0611872) q[0];
sx q[0];
rz(-0.10594254) q[0];
sx q[0];
rz(1.9796014) q[0];
x q[1];
rz(2.5162637) q[2];
sx q[2];
rz(-1.3742067) q[2];
sx q[2];
rz(2.1670053) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0137009) q[1];
sx q[1];
rz(-1.7096448) q[1];
sx q[1];
rz(1.9816887) q[1];
rz(-1.2563989) q[3];
sx q[3];
rz(-1.8125497) q[3];
sx q[3];
rz(0.35255656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54958582) q[2];
sx q[2];
rz(-2.0506115) q[2];
sx q[2];
rz(2.8263212) q[2];
rz(1.5677876) q[3];
sx q[3];
rz(-1.9939634) q[3];
sx q[3];
rz(1.1228336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0386117) q[0];
sx q[0];
rz(-0.6655612) q[0];
sx q[0];
rz(0.82157201) q[0];
rz(2.6577677) q[1];
sx q[1];
rz(-1.7211823) q[1];
sx q[1];
rz(0.22463591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53157887) q[0];
sx q[0];
rz(-1.2481127) q[0];
sx q[0];
rz(0.15845297) q[0];
x q[1];
rz(-2.5098652) q[2];
sx q[2];
rz(-2.1251107) q[2];
sx q[2];
rz(-0.71507031) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1146176) q[1];
sx q[1];
rz(-1.7791505) q[1];
sx q[1];
rz(2.2592553) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3566689) q[3];
sx q[3];
rz(-1.6928634) q[3];
sx q[3];
rz(1.3510974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3839174) q[2];
sx q[2];
rz(-0.11666798) q[2];
sx q[2];
rz(-0.095001027) q[2];
rz(1.4875686) q[3];
sx q[3];
rz(-0.58949685) q[3];
sx q[3];
rz(-2.3768363) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2471531) q[0];
sx q[0];
rz(-0.40774397) q[0];
sx q[0];
rz(2.8751873) q[0];
rz(-2.5447625) q[1];
sx q[1];
rz(-1.6886371) q[1];
sx q[1];
rz(1.4773038) q[1];
rz(1.1372139) q[2];
sx q[2];
rz(-0.72535338) q[2];
sx q[2];
rz(0.85272127) q[2];
rz(1.568902) q[3];
sx q[3];
rz(-1.1952778) q[3];
sx q[3];
rz(0.3921685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
