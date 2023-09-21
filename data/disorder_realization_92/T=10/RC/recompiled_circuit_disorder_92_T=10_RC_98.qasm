OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5151514) q[0];
sx q[0];
rz(-0.03349537) q[0];
sx q[0];
rz(1.7749696) q[0];
rz(2.0904436) q[1];
sx q[1];
rz(-1.6519974) q[1];
sx q[1];
rz(2.0096013) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7819408) q[0];
sx q[0];
rz(-1.9160761) q[0];
sx q[0];
rz(-2.3495673) q[0];
x q[1];
rz(-2.3425383) q[2];
sx q[2];
rz(-2.91215) q[2];
sx q[2];
rz(2.6829751) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7992026) q[1];
sx q[1];
rz(-0.7212762) q[1];
sx q[1];
rz(-0.93625416) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44191435) q[3];
sx q[3];
rz(-0.60097296) q[3];
sx q[3];
rz(2.0244983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91036096) q[2];
sx q[2];
rz(-1.3277162) q[2];
sx q[2];
rz(-1.1532016) q[2];
rz(2.6575346) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(0.4593862) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3020878) q[0];
sx q[0];
rz(-0.33058259) q[0];
sx q[0];
rz(2.6385345) q[0];
rz(1.5548276) q[1];
sx q[1];
rz(-2.4350872) q[1];
sx q[1];
rz(0.15393004) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4707697) q[0];
sx q[0];
rz(-1.5564859) q[0];
sx q[0];
rz(2.0736681) q[0];
x q[1];
rz(0.92127992) q[2];
sx q[2];
rz(-1.3037455) q[2];
sx q[2];
rz(-2.9608179) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0347552) q[1];
sx q[1];
rz(-1.991786) q[1];
sx q[1];
rz(-1.0616598) q[1];
rz(-1.5716343) q[3];
sx q[3];
rz(-1.0610233) q[3];
sx q[3];
rz(-1.2050932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2543891) q[2];
sx q[2];
rz(-0.35787359) q[2];
sx q[2];
rz(1.8187693) q[2];
rz(-1.4860738) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(-0.71050182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(0.94674295) q[0];
sx q[0];
rz(-1.1179593) q[0];
sx q[0];
rz(-2.2316566) q[0];
rz(2.3643156) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(2.1562703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1811718) q[0];
sx q[0];
rz(-2.4665678) q[0];
sx q[0];
rz(1.8280562) q[0];
x q[1];
rz(2.2276332) q[2];
sx q[2];
rz(-1.8310391) q[2];
sx q[2];
rz(-2.7514806) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9039771) q[1];
sx q[1];
rz(-1.3217889) q[1];
sx q[1];
rz(1.6275089) q[1];
rz(-pi) q[2];
rz(0.91089532) q[3];
sx q[3];
rz(-2.1351372) q[3];
sx q[3];
rz(-0.49163715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2788006) q[2];
sx q[2];
rz(-1.1547487) q[2];
sx q[2];
rz(-1.8939691) q[2];
rz(-0.4425846) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(-2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2414395) q[0];
sx q[0];
rz(-1.0192008) q[0];
sx q[0];
rz(2.0181657) q[0];
rz(-2.5627047) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(-1.7480063) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8907341) q[0];
sx q[0];
rz(-0.77461857) q[0];
sx q[0];
rz(-0.61268341) q[0];
rz(-pi) q[1];
x q[1];
rz(2.706714) q[2];
sx q[2];
rz(-1.41845) q[2];
sx q[2];
rz(2.9762714) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.046603831) q[1];
sx q[1];
rz(-2.2294083) q[1];
sx q[1];
rz(2.1556426) q[1];
rz(2.3062069) q[3];
sx q[3];
rz(-1.3961892) q[3];
sx q[3];
rz(2.01471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0397772) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(3.139479) q[2];
rz(0.56143108) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(2.5533365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.026022) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(1.2258688) q[0];
rz(1.4670124) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(1.3668758) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05058771) q[0];
sx q[0];
rz(-2.707621) q[0];
sx q[0];
rz(2.6758053) q[0];
x q[1];
rz(-1.5405802) q[2];
sx q[2];
rz(-2.317791) q[2];
sx q[2];
rz(-3.1291762) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3102476) q[1];
sx q[1];
rz(-0.86569769) q[1];
sx q[1];
rz(3.0369333) q[1];
x q[2];
rz(2.6392691) q[3];
sx q[3];
rz(-1.886743) q[3];
sx q[3];
rz(1.5497006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.86429578) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(-0.40536353) q[2];
rz(-0.46164414) q[3];
sx q[3];
rz(-2.3159537) q[3];
sx q[3];
rz(1.5464787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49155238) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(-2.6382085) q[0];
rz(0.21884306) q[1];
sx q[1];
rz(-1.2501406) q[1];
sx q[1];
rz(-2.4898081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1056021) q[0];
sx q[0];
rz(-0.30797568) q[0];
sx q[0];
rz(-0.67291798) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4621554) q[2];
sx q[2];
rz(-2.3896653) q[2];
sx q[2];
rz(1.3736563) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5469049) q[1];
sx q[1];
rz(-1.6225796) q[1];
sx q[1];
rz(1.4361708) q[1];
rz(-pi) q[2];
rz(0.28292803) q[3];
sx q[3];
rz(-1.471721) q[3];
sx q[3];
rz(0.84474746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6289604) q[2];
sx q[2];
rz(-1.7444538) q[2];
sx q[2];
rz(-0.39166489) q[2];
rz(2.9351249) q[3];
sx q[3];
rz(-0.73409096) q[3];
sx q[3];
rz(2.4519043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3702635) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(-0.96965924) q[0];
rz(0.58352739) q[1];
sx q[1];
rz(-1.1279761) q[1];
sx q[1];
rz(-0.13024174) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7182065) q[0];
sx q[0];
rz(-1.5209746) q[0];
sx q[0];
rz(0.35276619) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98353705) q[2];
sx q[2];
rz(-0.98896719) q[2];
sx q[2];
rz(-0.26435095) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.11634593) q[1];
sx q[1];
rz(-1.3533918) q[1];
sx q[1];
rz(-0.51699637) q[1];
rz(-pi) q[2];
rz(-0.3881012) q[3];
sx q[3];
rz(-1.023205) q[3];
sx q[3];
rz(0.91528085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51817259) q[2];
sx q[2];
rz(-1.5600081) q[2];
sx q[2];
rz(-2.0557892) q[2];
rz(-0.052224934) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(-2.0558555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6770342) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(1.1664671) q[0];
rz(-2.4160066) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(2.3988147) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60177875) q[0];
sx q[0];
rz(-2.302928) q[0];
sx q[0];
rz(1.8307122) q[0];
rz(-pi) q[1];
rz(0.6767512) q[2];
sx q[2];
rz(-1.1427715) q[2];
sx q[2];
rz(-1.8523491) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.038866) q[1];
sx q[1];
rz(-0.83737774) q[1];
sx q[1];
rz(1.0086165) q[1];
rz(-pi) q[2];
rz(-2.2320896) q[3];
sx q[3];
rz(-2.0561744) q[3];
sx q[3];
rz(0.73736008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5143738) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(-0.26088866) q[2];
rz(-1.1076814) q[3];
sx q[3];
rz(-1.8543782) q[3];
sx q[3];
rz(-1.0816983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35659197) q[0];
sx q[0];
rz(-0.78105015) q[0];
sx q[0];
rz(2.2055431) q[0];
rz(0.014135663) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(-2.4900808) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0148841) q[0];
sx q[0];
rz(-1.5944905) q[0];
sx q[0];
rz(-3.1372877) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65666171) q[2];
sx q[2];
rz(-2.0920394) q[2];
sx q[2];
rz(1.9929287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9459322) q[1];
sx q[1];
rz(-1.6911427) q[1];
sx q[1];
rz(1.8912998) q[1];
rz(-0.6833965) q[3];
sx q[3];
rz(-2.5960687) q[3];
sx q[3];
rz(-2.6664536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2953879) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(-2.6182168) q[2];
rz(2.8335617) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(1.3380922) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4073407) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(-1.7379606) q[0];
rz(-0.47406468) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(1.7636991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5869851) q[0];
sx q[0];
rz(-2.4025612) q[0];
sx q[0];
rz(-2.2686195) q[0];
rz(-pi) q[1];
rz(2.2950298) q[2];
sx q[2];
rz(-0.91869527) q[2];
sx q[2];
rz(-0.50154274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.024152012) q[1];
sx q[1];
rz(-1.9630868) q[1];
sx q[1];
rz(-1.3082318) q[1];
rz(-pi) q[2];
rz(2.8836807) q[3];
sx q[3];
rz(-2.5690418) q[3];
sx q[3];
rz(-0.27975988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.750981) q[2];
sx q[2];
rz(-1.6335952) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(1.4577929) q[3];
sx q[3];
rz(-1.0554353) q[3];
sx q[3];
rz(2.2916601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4703341) q[0];
sx q[0];
rz(-2.8699734) q[0];
sx q[0];
rz(1.097453) q[0];
rz(0.63256565) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(-1.8995646) q[2];
sx q[2];
rz(-0.9267926) q[2];
sx q[2];
rz(-0.37315858) q[2];
rz(0.87631638) q[3];
sx q[3];
rz(-2.1764168) q[3];
sx q[3];
rz(-0.62302667) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
