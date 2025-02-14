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
rz(-1.1408495) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716233) q[0];
sx q[0];
rz(-2.9282018) q[0];
sx q[0];
rz(-2.1687228) q[0];
rz(1.2873285) q[2];
sx q[2];
rz(-1.8745443) q[2];
sx q[2];
rz(1.8423353) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1124494) q[1];
sx q[1];
rz(-1.1894798) q[1];
sx q[1];
rz(-1.4277532) q[1];
x q[2];
rz(0.13762044) q[3];
sx q[3];
rz(-2.0064266) q[3];
sx q[3];
rz(-0.14996687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69316489) q[2];
sx q[2];
rz(-1.0591256) q[2];
sx q[2];
rz(1.1761752) q[2];
rz(3.0991992) q[3];
sx q[3];
rz(-1.8040801) q[3];
sx q[3];
rz(0.5347518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46928826) q[0];
sx q[0];
rz(-0.35457087) q[0];
sx q[0];
rz(2.9571423) q[0];
rz(-3.0448659) q[1];
sx q[1];
rz(-1.5238785) q[1];
sx q[1];
rz(-1.9116037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.205015) q[0];
sx q[0];
rz(-0.42213556) q[0];
sx q[0];
rz(-2.6279215) q[0];
rz(-pi) q[1];
rz(1.619009) q[2];
sx q[2];
rz(-2.1297751) q[2];
sx q[2];
rz(2.3490482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0470386) q[1];
sx q[1];
rz(-2.8287005) q[1];
sx q[1];
rz(2.439586) q[1];
rz(-pi) q[2];
rz(-1.6685247) q[3];
sx q[3];
rz(-0.37993452) q[3];
sx q[3];
rz(1.1539647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7868598) q[2];
sx q[2];
rz(-2.730098) q[2];
sx q[2];
rz(3.1301609) q[2];
rz(-1.3139906) q[3];
sx q[3];
rz(-1.3649789) q[3];
sx q[3];
rz(0.7029117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74362022) q[0];
sx q[0];
rz(-3.008606) q[0];
sx q[0];
rz(0.61070329) q[0];
rz(-0.01344219) q[1];
sx q[1];
rz(-1.8553269) q[1];
sx q[1];
rz(2.5450943) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4134408) q[0];
sx q[0];
rz(-0.96609173) q[0];
sx q[0];
rz(-0.20264969) q[0];
rz(-pi) q[1];
rz(0.99278583) q[2];
sx q[2];
rz(-1.3137378) q[2];
sx q[2];
rz(1.3959194) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2141419) q[1];
sx q[1];
rz(-2.3436556) q[1];
sx q[1];
rz(-2.6990456) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1085988) q[3];
sx q[3];
rz(-1.0171431) q[3];
sx q[3];
rz(-2.7102196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67768031) q[2];
sx q[2];
rz(-1.8296506) q[2];
sx q[2];
rz(0.00027351969) q[2];
rz(1.4759981) q[3];
sx q[3];
rz(-0.6690343) q[3];
sx q[3];
rz(-0.10703787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6855455) q[0];
sx q[0];
rz(-2.9750329) q[0];
sx q[0];
rz(-1.1628994) q[0];
rz(0.97745013) q[1];
sx q[1];
rz(-1.6012499) q[1];
sx q[1];
rz(-1.5862484) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9396203) q[0];
sx q[0];
rz(-0.012379025) q[0];
sx q[0];
rz(0.89039652) q[0];
x q[1];
rz(2.3562217) q[2];
sx q[2];
rz(-1.6024688) q[2];
sx q[2];
rz(2.9065064) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6628389) q[1];
sx q[1];
rz(-2.6445619) q[1];
sx q[1];
rz(3.1113487) q[1];
x q[2];
rz(2.8951696) q[3];
sx q[3];
rz(-1.0895673) q[3];
sx q[3];
rz(-2.6068166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9419452) q[2];
sx q[2];
rz(-2.6298099) q[2];
sx q[2];
rz(-2.9072348) q[2];
rz(-0.50218454) q[3];
sx q[3];
rz(-2.1000704) q[3];
sx q[3];
rz(2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9116354) q[0];
sx q[0];
rz(-1.3595694) q[0];
sx q[0];
rz(1.0372355) q[0];
rz(-pi) q[1];
rz(0.52509016) q[2];
sx q[2];
rz(-2.1222564) q[2];
sx q[2];
rz(2.9208825) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42758917) q[1];
sx q[1];
rz(-1.98381) q[1];
sx q[1];
rz(-0.16431134) q[1];
rz(-2.9850335) q[3];
sx q[3];
rz(-2.1053751) q[3];
sx q[3];
rz(-1.298157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56110567) q[2];
sx q[2];
rz(-1.9876391) q[2];
sx q[2];
rz(2.6055276) q[2];
rz(0.29340336) q[3];
sx q[3];
rz(-1.9304201) q[3];
sx q[3];
rz(-2.6611879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2368161) q[0];
sx q[0];
rz(-2.1892956) q[0];
sx q[0];
rz(-3.0425518) q[0];
rz(1.0320484) q[1];
sx q[1];
rz(-0.85056225) q[1];
sx q[1];
rz(1.7031857) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5779786) q[0];
sx q[0];
rz(-2.4017085) q[0];
sx q[0];
rz(-1.1881962) q[0];
x q[1];
rz(-0.24865227) q[2];
sx q[2];
rz(-2.8496345) q[2];
sx q[2];
rz(2.0338634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8107103) q[1];
sx q[1];
rz(-1.5941949) q[1];
sx q[1];
rz(2.3512273) q[1];
rz(-pi) q[2];
rz(-1.6910873) q[3];
sx q[3];
rz(-1.2748537) q[3];
sx q[3];
rz(-3.0927998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9859887) q[2];
sx q[2];
rz(-0.5039379) q[2];
sx q[2];
rz(-2.3206553) q[2];
rz(-1.692903) q[3];
sx q[3];
rz(-0.71805787) q[3];
sx q[3];
rz(-1.4887571) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2445225) q[0];
sx q[0];
rz(-0.91461602) q[0];
sx q[0];
rz(-0.50265092) q[0];
rz(-1.2333168) q[1];
sx q[1];
rz(-2.2063467) q[1];
sx q[1];
rz(2.1615692) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4122881) q[0];
sx q[0];
rz(-1.8192023) q[0];
sx q[0];
rz(-2.8969943) q[0];
x q[1];
rz(1.2792311) q[2];
sx q[2];
rz(-2.3254561) q[2];
sx q[2];
rz(-1.7486435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.087638559) q[1];
sx q[1];
rz(-0.77869895) q[1];
sx q[1];
rz(-2.8092572) q[1];
rz(-0.20670047) q[3];
sx q[3];
rz(-0.85089848) q[3];
sx q[3];
rz(0.18223079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9194453) q[2];
sx q[2];
rz(-1.8541226) q[2];
sx q[2];
rz(1.3745098) q[2];
rz(-1.8253271) q[3];
sx q[3];
rz(-2.2856183) q[3];
sx q[3];
rz(-2.159806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(0.027243622) q[0];
sx q[0];
rz(-2.7482432) q[0];
sx q[0];
rz(-2.223176) q[0];
rz(-2.9367327) q[1];
sx q[1];
rz(-1.058895) q[1];
sx q[1];
rz(1.6580261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61648864) q[0];
sx q[0];
rz(-1.2813998) q[0];
sx q[0];
rz(0.15286907) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7006364) q[2];
sx q[2];
rz(-2.3774638) q[2];
sx q[2];
rz(-1.6418599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0369954) q[1];
sx q[1];
rz(-2.2682087) q[1];
sx q[1];
rz(0.11499494) q[1];
rz(-1.6202043) q[3];
sx q[3];
rz(-0.8522343) q[3];
sx q[3];
rz(0.75203943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.001699) q[2];
sx q[2];
rz(-0.49109444) q[2];
sx q[2];
rz(-1.3341382) q[2];
rz(-0.88207465) q[3];
sx q[3];
rz(-2.4460654) q[3];
sx q[3];
rz(-1.2923366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081721574) q[0];
sx q[0];
rz(-2.7099755) q[0];
sx q[0];
rz(-0.01817848) q[0];
rz(-0.40953088) q[1];
sx q[1];
rz(-1.4298226) q[1];
sx q[1];
rz(-0.10083625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6438599) q[0];
sx q[0];
rz(-2.5063305) q[0];
sx q[0];
rz(3.0833901) q[0];
x q[1];
rz(-0.11282044) q[2];
sx q[2];
rz(-2.718528) q[2];
sx q[2];
rz(3.0213838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8162569) q[1];
sx q[1];
rz(-2.6995669) q[1];
sx q[1];
rz(-2.1620552) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0193897) q[3];
sx q[3];
rz(-2.8684385) q[3];
sx q[3];
rz(0.34452439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.33710256) q[2];
sx q[2];
rz(-3.0323995) q[2];
sx q[2];
rz(1.2479372) q[2];
rz(1.9167871) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(0.065751806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2357904) q[0];
sx q[0];
rz(-1.6658655) q[0];
sx q[0];
rz(-0.32050785) q[0];
rz(-2.7505752) q[1];
sx q[1];
rz(-1.1415488) q[1];
sx q[1];
rz(-1.5919707) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029767903) q[0];
sx q[0];
rz(-2.8643908) q[0];
sx q[0];
rz(0.76407822) q[0];
x q[1];
rz(-1.7334601) q[2];
sx q[2];
rz(-2.0792336) q[2];
sx q[2];
rz(0.66772738) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2992449) q[1];
sx q[1];
rz(-1.4834373) q[1];
sx q[1];
rz(1.9398938) q[1];
rz(0.42234588) q[3];
sx q[3];
rz(-1.5639439) q[3];
sx q[3];
rz(-0.50240483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0512507) q[2];
sx q[2];
rz(-2.0917459) q[2];
sx q[2];
rz(2.7412097) q[2];
rz(-2.9054437) q[3];
sx q[3];
rz(-0.103424) q[3];
sx q[3];
rz(-1.0108488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26265963) q[0];
sx q[0];
rz(-0.50983179) q[0];
sx q[0];
rz(-1.8806993) q[0];
rz(-2.6564468) q[1];
sx q[1];
rz(-2.3427675) q[1];
sx q[1];
rz(2.6642703) q[1];
rz(-1.7791228) q[2];
sx q[2];
rz(-1.4478417) q[2];
sx q[2];
rz(1.7461591) q[2];
rz(3.0780001) q[3];
sx q[3];
rz(-1.3432172) q[3];
sx q[3];
rz(-1.5061495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
