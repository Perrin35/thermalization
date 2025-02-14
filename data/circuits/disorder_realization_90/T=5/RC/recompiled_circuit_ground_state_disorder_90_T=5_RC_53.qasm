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
rz(3.314078) q[0];
sx q[0];
rz(9.4077851) q[0];
rz(0.51141557) q[1];
sx q[1];
rz(-0.49632448) q[1];
sx q[1];
rz(-2.6561148) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89428502) q[0];
sx q[0];
rz(-1.6818769) q[0];
sx q[0];
rz(-3.0283228) q[0];
rz(-pi) q[1];
rz(-0.85632773) q[2];
sx q[2];
rz(-2.1103139) q[2];
sx q[2];
rz(-1.0401638) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.456326) q[1];
sx q[1];
rz(-0.86325544) q[1];
sx q[1];
rz(-0.68143845) q[1];
x q[2];
rz(0.11459132) q[3];
sx q[3];
rz(-1.262731) q[3];
sx q[3];
rz(-1.1685355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7444676) q[2];
sx q[2];
rz(-3.0391389) q[2];
sx q[2];
rz(-0.79258072) q[2];
rz(0.98627311) q[3];
sx q[3];
rz(-1.4873742) q[3];
sx q[3];
rz(-0.13315323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88696402) q[0];
sx q[0];
rz(-1.6034842) q[0];
sx q[0];
rz(-0.1057374) q[0];
rz(2.3836783) q[1];
sx q[1];
rz(-0.71280232) q[1];
sx q[1];
rz(-0.47592083) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3041046) q[0];
sx q[0];
rz(-1.473635) q[0];
sx q[0];
rz(0.23529737) q[0];
x q[1];
rz(-1.4272825) q[2];
sx q[2];
rz(-2.0594308) q[2];
sx q[2];
rz(1.8861063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9502752) q[1];
sx q[1];
rz(-1.4498267) q[1];
sx q[1];
rz(-2.1410393) q[1];
x q[2];
rz(-2.2146899) q[3];
sx q[3];
rz(-1.3312648) q[3];
sx q[3];
rz(0.86349559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14625749) q[2];
sx q[2];
rz(-2.6671851) q[2];
sx q[2];
rz(-1.3169301) q[2];
rz(-0.82777348) q[3];
sx q[3];
rz(-2.0483978) q[3];
sx q[3];
rz(2.304346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71883172) q[0];
sx q[0];
rz(-1.7976924) q[0];
sx q[0];
rz(-1.2368917) q[0];
rz(-1.749136) q[1];
sx q[1];
rz(-1.720865) q[1];
sx q[1];
rz(0.69033355) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25094098) q[0];
sx q[0];
rz(-2.0253882) q[0];
sx q[0];
rz(-2.5174052) q[0];
rz(-1.5788011) q[2];
sx q[2];
rz(-2.6782616) q[2];
sx q[2];
rz(-3.0750907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0818881) q[1];
sx q[1];
rz(-1.7374542) q[1];
sx q[1];
rz(2.9465066) q[1];
rz(-pi) q[2];
rz(2.1851483) q[3];
sx q[3];
rz(-2.2057475) q[3];
sx q[3];
rz(-2.3243543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51848015) q[2];
sx q[2];
rz(-1.6192351) q[2];
sx q[2];
rz(3.0677262) q[2];
rz(-1.1582003) q[3];
sx q[3];
rz(-1.9246212) q[3];
sx q[3];
rz(-1.4069675) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2288007) q[0];
sx q[0];
rz(-0.055963628) q[0];
sx q[0];
rz(0.19293109) q[0];
rz(1.2436766) q[1];
sx q[1];
rz(-1.7453777) q[1];
sx q[1];
rz(-0.23304932) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4547494) q[0];
sx q[0];
rz(-2.966724) q[0];
sx q[0];
rz(-2.5681679) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1036699) q[2];
sx q[2];
rz(-0.86235986) q[2];
sx q[2];
rz(1.0346827) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31598791) q[1];
sx q[1];
rz(-1.627076) q[1];
sx q[1];
rz(0.31096267) q[1];
rz(0.95121164) q[3];
sx q[3];
rz(-1.9635634) q[3];
sx q[3];
rz(-2.5529566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9243246) q[2];
sx q[2];
rz(-1.7622207) q[2];
sx q[2];
rz(-0.064083727) q[2];
rz(-2.890375) q[3];
sx q[3];
rz(-0.46291864) q[3];
sx q[3];
rz(-0.29138756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083704405) q[0];
sx q[0];
rz(-1.2936445) q[0];
sx q[0];
rz(-1.8573014) q[0];
rz(-3.1341556) q[1];
sx q[1];
rz(-1.1859272) q[1];
sx q[1];
rz(0.95710212) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5610661) q[0];
sx q[0];
rz(-1.8883123) q[0];
sx q[0];
rz(-1.1098212) q[0];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0937371) q[1];
sx q[1];
rz(-1.4720494) q[1];
sx q[1];
rz(-0.32961032) q[1];
x q[2];
rz(0.24139054) q[3];
sx q[3];
rz(-2.4721535) q[3];
sx q[3];
rz(-2.2934492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1288422) q[2];
sx q[2];
rz(-2.4350171) q[2];
sx q[2];
rz(1.0303222) q[2];
rz(-2.4230867) q[3];
sx q[3];
rz(-2.0593144) q[3];
sx q[3];
rz(-0.70657402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71501032) q[0];
sx q[0];
rz(-2.3937245) q[0];
sx q[0];
rz(2.7744875) q[0];
rz(2.7857065) q[1];
sx q[1];
rz(-1.117319) q[1];
sx q[1];
rz(-1.5497367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1995576) q[0];
sx q[0];
rz(-0.74239391) q[0];
sx q[0];
rz(3.0846157) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8086967) q[2];
sx q[2];
rz(-0.44778433) q[2];
sx q[2];
rz(2.9749123) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8425297) q[1];
sx q[1];
rz(-0.46485126) q[1];
sx q[1];
rz(-0.47375394) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14629062) q[3];
sx q[3];
rz(-0.47463575) q[3];
sx q[3];
rz(2.4865347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1648569) q[2];
sx q[2];
rz(-1.848449) q[2];
sx q[2];
rz(-3.1381651) q[2];
rz(0.88611832) q[3];
sx q[3];
rz(-1.0930073) q[3];
sx q[3];
rz(1.4440943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43585983) q[0];
sx q[0];
rz(-1.6678565) q[0];
sx q[0];
rz(2.4263897) q[0];
rz(-3.0192979) q[1];
sx q[1];
rz(-2.1221752) q[1];
sx q[1];
rz(0.14762793) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58260949) q[0];
sx q[0];
rz(-1.2319854) q[0];
sx q[0];
rz(0.34783439) q[0];
x q[1];
rz(-1.1806025) q[2];
sx q[2];
rz(-1.1653333) q[2];
sx q[2];
rz(2.6801339) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8898006) q[1];
sx q[1];
rz(-1.3887059) q[1];
sx q[1];
rz(-0.86343335) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8054923) q[3];
sx q[3];
rz(-2.8471591) q[3];
sx q[3];
rz(2.2469843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3132396) q[2];
sx q[2];
rz(-2.8326663) q[2];
sx q[2];
rz(0.1114791) q[2];
rz(0.76505032) q[3];
sx q[3];
rz(-1.7181516) q[3];
sx q[3];
rz(-0.033871977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2400804) q[0];
sx q[0];
rz(-0.96918786) q[0];
sx q[0];
rz(-1.0028268) q[0];
rz(0.78549939) q[1];
sx q[1];
rz(-1.6167043) q[1];
sx q[1];
rz(-1.921382) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.553279) q[0];
sx q[0];
rz(-1.0568585) q[0];
sx q[0];
rz(1.5470424) q[0];
rz(1.6613879) q[2];
sx q[2];
rz(-1.633051) q[2];
sx q[2];
rz(-1.412815) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3165908) q[1];
sx q[1];
rz(-1.8866393) q[1];
sx q[1];
rz(1.0388253) q[1];
x q[2];
rz(-2.972159) q[3];
sx q[3];
rz(-0.9421351) q[3];
sx q[3];
rz(-2.7055912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23585524) q[2];
sx q[2];
rz(-1.3964272) q[2];
sx q[2];
rz(-0.033795707) q[2];
rz(2.3679521) q[3];
sx q[3];
rz(-1.6126817) q[3];
sx q[3];
rz(-2.0788367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77699023) q[0];
sx q[0];
rz(-2.5211625) q[0];
sx q[0];
rz(0.34238368) q[0];
rz(1.7204334) q[1];
sx q[1];
rz(-0.8232638) q[1];
sx q[1];
rz(-1.294543) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50534821) q[0];
sx q[0];
rz(-2.7490135) q[0];
sx q[0];
rz(2.7468203) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72118463) q[2];
sx q[2];
rz(-0.9936665) q[2];
sx q[2];
rz(-2.3604928) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0755989) q[1];
sx q[1];
rz(-2.3477049) q[1];
sx q[1];
rz(-0.17236472) q[1];
rz(-pi) q[2];
rz(2.6813981) q[3];
sx q[3];
rz(-2.7685173) q[3];
sx q[3];
rz(2.4378547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8483868) q[2];
sx q[2];
rz(-1.3853955) q[2];
sx q[2];
rz(-2.6050341) q[2];
rz(2.7287591) q[3];
sx q[3];
rz(-0.53240132) q[3];
sx q[3];
rz(2.0280973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31529108) q[0];
sx q[0];
rz(-1.2696126) q[0];
sx q[0];
rz(1.6444561) q[0];
rz(2.3313088) q[1];
sx q[1];
rz(-1.5850001) q[1];
sx q[1];
rz(-2.2946045) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2772781) q[0];
sx q[0];
rz(-1.7314096) q[0];
sx q[0];
rz(1.3993652) q[0];
x q[1];
rz(2.2184847) q[2];
sx q[2];
rz(-2.0294242) q[2];
sx q[2];
rz(3.0748526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0125186) q[1];
sx q[1];
rz(-1.301494) q[1];
sx q[1];
rz(-2.49519) q[1];
x q[2];
rz(-1.8448006) q[3];
sx q[3];
rz(-1.6633667) q[3];
sx q[3];
rz(1.1082197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0392796) q[2];
sx q[2];
rz(-1.8259093) q[2];
sx q[2];
rz(0.63423356) q[2];
rz(-0.63747326) q[3];
sx q[3];
rz(-1.1280779) q[3];
sx q[3];
rz(-0.2171966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1524326) q[0];
sx q[0];
rz(-1.9518873) q[0];
sx q[0];
rz(1.963203) q[0];
rz(0.70869008) q[1];
sx q[1];
rz(-1.3460174) q[1];
sx q[1];
rz(-1.8585471) q[1];
rz(-1.218956) q[2];
sx q[2];
rz(-1.9250122) q[2];
sx q[2];
rz(-0.4831947) q[2];
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
