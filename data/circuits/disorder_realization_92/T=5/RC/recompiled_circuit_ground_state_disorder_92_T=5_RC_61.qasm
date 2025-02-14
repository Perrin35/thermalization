OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.3761223) q[0];
sx q[0];
rz(-1.0785311) q[0];
sx q[0];
rz(1.7662788) q[0];
rz(0.41962656) q[1];
sx q[1];
rz(-1.761275) q[1];
sx q[1];
rz(0.79545155) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8859098) q[0];
sx q[0];
rz(-1.943318) q[0];
sx q[0];
rz(1.6557978) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4515204) q[2];
sx q[2];
rz(-1.96835) q[2];
sx q[2];
rz(-1.2326604) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.011820531) q[1];
sx q[1];
rz(-1.3737965) q[1];
sx q[1];
rz(2.4596766) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.086049155) q[3];
sx q[3];
rz(-1.8294191) q[3];
sx q[3];
rz(1.8515967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2835283) q[2];
sx q[2];
rz(-2.8499446) q[2];
sx q[2];
rz(0.38823271) q[2];
rz(2.9325716) q[3];
sx q[3];
rz(-1.4744604) q[3];
sx q[3];
rz(0.89116636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6727869) q[0];
sx q[0];
rz(-2.6728632) q[0];
sx q[0];
rz(-0.76954532) q[0];
rz(2.1107213) q[1];
sx q[1];
rz(-2.2831235) q[1];
sx q[1];
rz(-1.4268202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1311817) q[0];
sx q[0];
rz(-0.37574925) q[0];
sx q[0];
rz(0.6300169) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67536036) q[2];
sx q[2];
rz(-1.5262233) q[2];
sx q[2];
rz(-1.7365484) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4974792) q[1];
sx q[1];
rz(-2.1783324) q[1];
sx q[1];
rz(1.2588905) q[1];
x q[2];
rz(0.85635045) q[3];
sx q[3];
rz(-1.6110241) q[3];
sx q[3];
rz(2.2742607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0017121) q[2];
sx q[2];
rz(-0.44846815) q[2];
sx q[2];
rz(1.2343538) q[2];
rz(2.6185696) q[3];
sx q[3];
rz(-1.1250027) q[3];
sx q[3];
rz(-0.6559059) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457542) q[0];
sx q[0];
rz(-0.0025175968) q[0];
sx q[0];
rz(2.5939831) q[0];
rz(-2.2721263) q[1];
sx q[1];
rz(-0.97620669) q[1];
sx q[1];
rz(-1.6798103) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6637409) q[0];
sx q[0];
rz(-1.9123075) q[0];
sx q[0];
rz(0.44490258) q[0];
rz(-2.8732576) q[2];
sx q[2];
rz(-0.78040022) q[2];
sx q[2];
rz(2.4805299) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.853735) q[1];
sx q[1];
rz(-0.53799858) q[1];
sx q[1];
rz(-0.96189672) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2814155) q[3];
sx q[3];
rz(-1.0499991) q[3];
sx q[3];
rz(1.2788683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4168641) q[2];
sx q[2];
rz(-0.29837307) q[2];
sx q[2];
rz(-0.63428632) q[2];
rz(1.0861081) q[3];
sx q[3];
rz(-0.48091286) q[3];
sx q[3];
rz(1.1336795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61741614) q[0];
sx q[0];
rz(-2.9845181) q[0];
sx q[0];
rz(-1.6198535) q[0];
rz(0.98328439) q[1];
sx q[1];
rz(-2.0359928) q[1];
sx q[1];
rz(0.081667893) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3043609) q[0];
sx q[0];
rz(-1.0153234) q[0];
sx q[0];
rz(0.13120478) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6105725) q[2];
sx q[2];
rz(-2.808254) q[2];
sx q[2];
rz(3.1266629) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28934902) q[1];
sx q[1];
rz(-1.0846347) q[1];
sx q[1];
rz(2.6527609) q[1];
rz(-0.053728624) q[3];
sx q[3];
rz(-1.0815291) q[3];
sx q[3];
rz(0.51988822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6433158) q[2];
sx q[2];
rz(-1.4772819) q[2];
sx q[2];
rz(2.2225883) q[2];
rz(2.8194341) q[3];
sx q[3];
rz(-1.6291658) q[3];
sx q[3];
rz(-0.34644103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084694013) q[0];
sx q[0];
rz(-1.1321122) q[0];
sx q[0];
rz(1.803501) q[0];
rz(-0.40254205) q[1];
sx q[1];
rz(-1.891338) q[1];
sx q[1];
rz(-1.58443) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4854003) q[0];
sx q[0];
rz(-1.3043405) q[0];
sx q[0];
rz(-2.0977667) q[0];
x q[1];
rz(2.3473559) q[2];
sx q[2];
rz(-1.4068687) q[2];
sx q[2];
rz(-1.2212703) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.118972) q[1];
sx q[1];
rz(-1.717219) q[1];
sx q[1];
rz(-2.6048003) q[1];
rz(-pi) q[2];
rz(-2.040287) q[3];
sx q[3];
rz(-1.8113422) q[3];
sx q[3];
rz(2.8194373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0255787) q[2];
sx q[2];
rz(-1.8639001) q[2];
sx q[2];
rz(-0.73520994) q[2];
rz(2.0871346) q[3];
sx q[3];
rz(-1.0337044) q[3];
sx q[3];
rz(1.2247156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16578199) q[0];
sx q[0];
rz(-2.7175792) q[0];
sx q[0];
rz(-1.8884416) q[0];
rz(-1.5755298) q[1];
sx q[1];
rz(-1.1352481) q[1];
sx q[1];
rz(-0.53322405) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5767007) q[0];
sx q[0];
rz(-2.3561888) q[0];
sx q[0];
rz(-1.2346953) q[0];
rz(-pi) q[1];
x q[1];
rz(0.511325) q[2];
sx q[2];
rz(-1.4361524) q[2];
sx q[2];
rz(2.0106778) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6937852) q[1];
sx q[1];
rz(-0.29305036) q[1];
sx q[1];
rz(2.4386921) q[1];
rz(-2.1142152) q[3];
sx q[3];
rz(-1.0735596) q[3];
sx q[3];
rz(1.1641722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.48625654) q[2];
sx q[2];
rz(-2.76261) q[2];
sx q[2];
rz(0.8302702) q[2];
rz(1.8298979) q[3];
sx q[3];
rz(-2.0528841) q[3];
sx q[3];
rz(-2.8970498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61520064) q[0];
sx q[0];
rz(-1.4854234) q[0];
sx q[0];
rz(-2.3080589) q[0];
rz(-2.1444881) q[1];
sx q[1];
rz(-1.3778967) q[1];
sx q[1];
rz(1.5901784) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98956052) q[0];
sx q[0];
rz(-2.6989134) q[0];
sx q[0];
rz(-0.17501207) q[0];
rz(-pi) q[1];
rz(-0.35153206) q[2];
sx q[2];
rz(-2.1334634) q[2];
sx q[2];
rz(1.2219791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4287069) q[1];
sx q[1];
rz(-0.16413153) q[1];
sx q[1];
rz(-0.81961373) q[1];
rz(1.8206014) q[3];
sx q[3];
rz(-2.0916998) q[3];
sx q[3];
rz(-1.5977719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3414063) q[2];
sx q[2];
rz(-2.1087346) q[2];
sx q[2];
rz(1.3979073) q[2];
rz(2.7684815) q[3];
sx q[3];
rz(-2.2171376) q[3];
sx q[3];
rz(2.1227409) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7031192) q[0];
sx q[0];
rz(-1.4716453) q[0];
sx q[0];
rz(-2.8224509) q[0];
rz(-1.1646264) q[1];
sx q[1];
rz(-1.171037) q[1];
sx q[1];
rz(-3.107531) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24821407) q[0];
sx q[0];
rz(-2.3015018) q[0];
sx q[0];
rz(2.7284751) q[0];
rz(-pi) q[1];
rz(2.8197779) q[2];
sx q[2];
rz(-0.79635598) q[2];
sx q[2];
rz(-2.1200695) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0553375) q[1];
sx q[1];
rz(-0.024340425) q[1];
sx q[1];
rz(-2.023995) q[1];
x q[2];
rz(1.930792) q[3];
sx q[3];
rz(-2.0742831) q[3];
sx q[3];
rz(-0.22758037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0762735) q[2];
sx q[2];
rz(-1.7460145) q[2];
sx q[2];
rz(-0.070153959) q[2];
rz(-1.3089199) q[3];
sx q[3];
rz(-2.2891243) q[3];
sx q[3];
rz(-0.13516983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081130505) q[0];
sx q[0];
rz(-1.9996996) q[0];
sx q[0];
rz(0.4749701) q[0];
rz(1.9604856) q[1];
sx q[1];
rz(-0.89473748) q[1];
sx q[1];
rz(-0.917135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4297535) q[0];
sx q[0];
rz(-1.8022402) q[0];
sx q[0];
rz(-2.5654456) q[0];
rz(-0.31355942) q[2];
sx q[2];
rz(-0.90679996) q[2];
sx q[2];
rz(-1.0194376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.57130206) q[1];
sx q[1];
rz(-1.362974) q[1];
sx q[1];
rz(-3.057722) q[1];
rz(-2.6090129) q[3];
sx q[3];
rz(-0.6709698) q[3];
sx q[3];
rz(2.4155145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8721623) q[2];
sx q[2];
rz(-1.7717125) q[2];
sx q[2];
rz(1.5577462) q[2];
rz(2.9861084) q[3];
sx q[3];
rz(-1.0289861) q[3];
sx q[3];
rz(0.99229971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6601722) q[0];
sx q[0];
rz(-1.6654797) q[0];
sx q[0];
rz(0.36556622) q[0];
rz(1.9407326) q[1];
sx q[1];
rz(-1.6114085) q[1];
sx q[1];
rz(2.0948804) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4037673) q[0];
sx q[0];
rz(-1.0407454) q[0];
sx q[0];
rz(0.60204317) q[0];
rz(0.44136845) q[2];
sx q[2];
rz(-2.0725996) q[2];
sx q[2];
rz(0.29870118) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3550008) q[1];
sx q[1];
rz(-2.2552571) q[1];
sx q[1];
rz(-0.85620086) q[1];
rz(1.0250859) q[3];
sx q[3];
rz(-2.0841597) q[3];
sx q[3];
rz(-1.8051763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6863579) q[2];
sx q[2];
rz(-2.0265667) q[2];
sx q[2];
rz(-2.4380747) q[2];
rz(-0.55758682) q[3];
sx q[3];
rz(-2.360354) q[3];
sx q[3];
rz(0.47344661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8483509) q[0];
sx q[0];
rz(-2.1108755) q[0];
sx q[0];
rz(-2.1847771) q[0];
rz(-1.035607) q[1];
sx q[1];
rz(-2.5685812) q[1];
sx q[1];
rz(-2.0249637) q[1];
rz(-0.66813095) q[2];
sx q[2];
rz(-1.1429759) q[2];
sx q[2];
rz(-0.10980724) q[2];
rz(1.5286636) q[3];
sx q[3];
rz(-1.3535454) q[3];
sx q[3];
rz(-2.7342697) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
