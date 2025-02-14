OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.056501) q[0];
sx q[0];
rz(-2.8488475) q[0];
sx q[0];
rz(-2.2226287) q[0];
rz(2.7594944) q[1];
sx q[1];
rz(-0.16799071) q[1];
sx q[1];
rz(2.0007432) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1534682) q[0];
sx q[0];
rz(-1.690295) q[0];
sx q[0];
rz(-1.7480127) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72854002) q[2];
sx q[2];
rz(-0.41245261) q[2];
sx q[2];
rz(-2.0714687) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1124494) q[1];
sx q[1];
rz(-1.9521128) q[1];
sx q[1];
rz(-1.4277532) q[1];
rz(-pi) q[2];
rz(1.1315207) q[3];
sx q[3];
rz(-1.695493) q[3];
sx q[3];
rz(-1.6623868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69316489) q[2];
sx q[2];
rz(-2.0824671) q[2];
sx q[2];
rz(1.1761752) q[2];
rz(3.0991992) q[3];
sx q[3];
rz(-1.3375125) q[3];
sx q[3];
rz(2.6068408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6723044) q[0];
sx q[0];
rz(-0.35457087) q[0];
sx q[0];
rz(-2.9571423) q[0];
rz(-3.0448659) q[1];
sx q[1];
rz(-1.6177142) q[1];
sx q[1];
rz(-1.2299889) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.205015) q[0];
sx q[0];
rz(-2.7194571) q[0];
sx q[0];
rz(2.6279215) q[0];
x q[1];
rz(1.5225836) q[2];
sx q[2];
rz(-2.1297751) q[2];
sx q[2];
rz(-2.3490482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7736771) q[1];
sx q[1];
rz(-1.333548) q[1];
sx q[1];
rz(-1.3648454) q[1];
x q[2];
rz(0.038944728) q[3];
sx q[3];
rz(-1.9488261) q[3];
sx q[3];
rz(-1.259144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.35473287) q[2];
sx q[2];
rz(-0.41149461) q[2];
sx q[2];
rz(-0.011431781) q[2];
rz(1.3139906) q[3];
sx q[3];
rz(-1.7766137) q[3];
sx q[3];
rz(-2.438681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3979724) q[0];
sx q[0];
rz(-0.13298661) q[0];
sx q[0];
rz(-0.61070329) q[0];
rz(3.1281505) q[1];
sx q[1];
rz(-1.2862658) q[1];
sx q[1];
rz(0.59649831) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0749576) q[0];
sx q[0];
rz(-0.63371113) q[0];
sx q[0];
rz(1.8541965) q[0];
rz(-0.30412401) q[2];
sx q[2];
rz(-1.0140739) q[2];
sx q[2];
rz(0.010526882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1787179) q[1];
sx q[1];
rz(-1.8824008) q[1];
sx q[1];
rz(-2.3942663) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1246931) q[3];
sx q[3];
rz(-1.5988598) q[3];
sx q[3];
rz(1.156776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4639123) q[2];
sx q[2];
rz(-1.3119421) q[2];
sx q[2];
rz(-0.00027351969) q[2];
rz(1.6655946) q[3];
sx q[3];
rz(-2.4725584) q[3];
sx q[3];
rz(-0.10703787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45604712) q[0];
sx q[0];
rz(-2.9750329) q[0];
sx q[0];
rz(-1.1628994) q[0];
rz(-0.97745013) q[1];
sx q[1];
rz(-1.6012499) q[1];
sx q[1];
rz(-1.5553442) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88240963) q[0];
sx q[0];
rz(-1.5804187) q[0];
sx q[0];
rz(-3.1338047) q[0];
x q[1];
rz(3.0968148) q[2];
sx q[2];
rz(-2.3557202) q[2];
sx q[2];
rz(1.3673683) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0229637) q[1];
sx q[1];
rz(-1.5563772) q[1];
sx q[1];
rz(2.6447536) q[1];
x q[2];
rz(-2.0647207) q[3];
sx q[3];
rz(-1.3528429) q[3];
sx q[3];
rz(0.92011425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.9419452) q[2];
sx q[2];
rz(-2.6298099) q[2];
sx q[2];
rz(-0.23435782) q[2];
rz(-2.6394081) q[3];
sx q[3];
rz(-2.1000704) q[3];
sx q[3];
rz(0.65619367) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31768826) q[0];
sx q[0];
rz(-0.81517878) q[0];
sx q[0];
rz(-1.748388) q[0];
rz(-2.4488917) q[1];
sx q[1];
rz(-1.0857948) q[1];
sx q[1];
rz(1.0903953) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2176184) q[0];
sx q[0];
rz(-1.0503142) q[0];
sx q[0];
rz(0.24407669) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1887921) q[2];
sx q[2];
rz(-2.0119466) q[2];
sx q[2];
rz(-1.6447826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0648374) q[1];
sx q[1];
rz(-1.7211815) q[1];
sx q[1];
rz(-1.1527747) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15655915) q[3];
sx q[3];
rz(-2.1053751) q[3];
sx q[3];
rz(-1.298157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56110567) q[2];
sx q[2];
rz(-1.1539536) q[2];
sx q[2];
rz(0.53606501) q[2];
rz(0.29340336) q[3];
sx q[3];
rz(-1.2111726) q[3];
sx q[3];
rz(-0.48040473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2368161) q[0];
sx q[0];
rz(-0.95229709) q[0];
sx q[0];
rz(-3.0425518) q[0];
rz(1.0320484) q[1];
sx q[1];
rz(-0.85056225) q[1];
sx q[1];
rz(-1.438407) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.563614) q[0];
sx q[0];
rz(-2.4017085) q[0];
sx q[0];
rz(-1.1881962) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8929404) q[2];
sx q[2];
rz(-0.2919582) q[2];
sx q[2];
rz(-2.0338634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3308823) q[1];
sx q[1];
rz(-1.5473978) q[1];
sx q[1];
rz(2.3512273) q[1];
rz(-pi) q[2];
rz(-2.766633) q[3];
sx q[3];
rz(-2.8228033) q[3];
sx q[3];
rz(-0.44170435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9859887) q[2];
sx q[2];
rz(-2.6376548) q[2];
sx q[2];
rz(2.3206553) q[2];
rz(1.692903) q[3];
sx q[3];
rz(-0.71805787) q[3];
sx q[3];
rz(-1.6528355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2445225) q[0];
sx q[0];
rz(-0.91461602) q[0];
sx q[0];
rz(-2.6389417) q[0];
rz(-1.2333168) q[1];
sx q[1];
rz(-2.2063467) q[1];
sx q[1];
rz(-0.9800235) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61934793) q[0];
sx q[0];
rz(-2.7947593) q[0];
sx q[0];
rz(-0.808544) q[0];
rz(-pi) q[1];
rz(-2.365391) q[2];
sx q[2];
rz(-1.3598249) q[2];
sx q[2];
rz(-2.7609776) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7777694) q[1];
sx q[1];
rz(-0.8447434) q[1];
sx q[1];
rz(-1.8822303) q[1];
rz(-pi) q[2];
rz(2.3013744) q[3];
sx q[3];
rz(-1.4158691) q[3];
sx q[3];
rz(1.8904101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9194453) q[2];
sx q[2];
rz(-1.8541226) q[2];
sx q[2];
rz(1.7670828) q[2];
rz(1.3162656) q[3];
sx q[3];
rz(-2.2856183) q[3];
sx q[3];
rz(0.98178664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027243622) q[0];
sx q[0];
rz(-0.39334941) q[0];
sx q[0];
rz(2.223176) q[0];
rz(2.9367327) q[1];
sx q[1];
rz(-1.058895) q[1];
sx q[1];
rz(-1.6580261) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.525104) q[0];
sx q[0];
rz(-1.8601928) q[0];
sx q[0];
rz(-0.15286907) q[0];
x q[1];
rz(1.7006364) q[2];
sx q[2];
rz(-2.3774638) q[2];
sx q[2];
rz(1.6418599) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0369954) q[1];
sx q[1];
rz(-2.2682087) q[1];
sx q[1];
rz(-3.0265977) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5213883) q[3];
sx q[3];
rz(-0.8522343) q[3];
sx q[3];
rz(-2.3895532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1398937) q[2];
sx q[2];
rz(-2.6504982) q[2];
sx q[2];
rz(-1.8074544) q[2];
rz(2.259518) q[3];
sx q[3];
rz(-2.4460654) q[3];
sx q[3];
rz(1.849256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081721574) q[0];
sx q[0];
rz(-0.43161714) q[0];
sx q[0];
rz(-0.01817848) q[0];
rz(2.7320618) q[1];
sx q[1];
rz(-1.7117701) q[1];
sx q[1];
rz(0.10083625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026196711) q[0];
sx q[0];
rz(-1.6053204) q[0];
sx q[0];
rz(0.63445292) q[0];
rz(-pi) q[1];
rz(-0.11282044) q[2];
sx q[2];
rz(-0.42306468) q[2];
sx q[2];
rz(-3.0213838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8162569) q[1];
sx q[1];
rz(-0.44202572) q[1];
sx q[1];
rz(-2.1620552) q[1];
x q[2];
rz(1.1222029) q[3];
sx q[3];
rz(-2.8684385) q[3];
sx q[3];
rz(2.7970683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33710256) q[2];
sx q[2];
rz(-0.10919315) q[2];
sx q[2];
rz(-1.2479372) q[2];
rz(-1.9167871) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(3.0758408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2357904) q[0];
sx q[0];
rz(-1.6658655) q[0];
sx q[0];
rz(2.8210848) q[0];
rz(0.39101741) q[1];
sx q[1];
rz(-1.1415488) q[1];
sx q[1];
rz(1.549622) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3282916) q[0];
sx q[0];
rz(-1.3718954) q[0];
sx q[0];
rz(1.7651674) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8587927) q[2];
sx q[2];
rz(-0.53164266) q[2];
sx q[2];
rz(-0.34257364) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1913073) q[1];
sx q[1];
rz(-0.37883329) q[1];
sx q[1];
rz(1.3326419) q[1];
rz(-pi) q[2];
rz(-1.5632838) q[3];
sx q[3];
rz(-1.148461) q[3];
sx q[3];
rz(-1.0653121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0512507) q[2];
sx q[2];
rz(-1.0498468) q[2];
sx q[2];
rz(-2.7412097) q[2];
rz(-0.23614899) q[3];
sx q[3];
rz(-3.0381687) q[3];
sx q[3];
rz(-1.0108488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.878933) q[0];
sx q[0];
rz(-0.50983179) q[0];
sx q[0];
rz(-1.8806993) q[0];
rz(-2.6564468) q[1];
sx q[1];
rz(-2.3427675) q[1];
sx q[1];
rz(2.6642703) q[1];
rz(-1.3624698) q[2];
sx q[2];
rz(-1.6937509) q[2];
sx q[2];
rz(-1.3954336) q[2];
rz(0.063592576) q[3];
sx q[3];
rz(-1.7983754) q[3];
sx q[3];
rz(1.6354431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
