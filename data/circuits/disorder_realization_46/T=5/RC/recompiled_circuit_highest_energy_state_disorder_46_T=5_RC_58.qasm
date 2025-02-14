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
rz(-1.3396076) q[0];
sx q[0];
rz(-0.34914246) q[0];
sx q[0];
rz(2.1139297) q[0];
rz(3.1066306) q[1];
sx q[1];
rz(-2.1244013) q[1];
sx q[1];
rz(-1.6962167) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.611269) q[0];
sx q[0];
rz(-1.230254) q[0];
sx q[0];
rz(-2.5304619) q[0];
rz(2.9541763) q[2];
sx q[2];
rz(-0.18067154) q[2];
sx q[2];
rz(0.60127163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6563182) q[1];
sx q[1];
rz(-2.7492011) q[1];
sx q[1];
rz(-2.5618784) q[1];
x q[2];
rz(0.39528109) q[3];
sx q[3];
rz(-0.78215212) q[3];
sx q[3];
rz(-0.30667337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8029636) q[2];
sx q[2];
rz(-2.6821319) q[2];
sx q[2];
rz(-2.6098693) q[2];
rz(-1.3945329) q[3];
sx q[3];
rz(-1.8683542) q[3];
sx q[3];
rz(-1.9664221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8748473) q[0];
sx q[0];
rz(-1.6371472) q[0];
sx q[0];
rz(2.3496085) q[0];
rz(-0.92164552) q[1];
sx q[1];
rz(-1.4896723) q[1];
sx q[1];
rz(-2.0335061) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8130694) q[0];
sx q[0];
rz(-1.6824772) q[0];
sx q[0];
rz(2.172676) q[0];
rz(-pi) q[1];
rz(2.1943548) q[2];
sx q[2];
rz(-1.531336) q[2];
sx q[2];
rz(-0.90496162) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11985677) q[1];
sx q[1];
rz(-0.83470063) q[1];
sx q[1];
rz(-1.0480919) q[1];
rz(-pi) q[2];
rz(1.8872126) q[3];
sx q[3];
rz(-1.5174447) q[3];
sx q[3];
rz(-2.2406468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38829982) q[2];
sx q[2];
rz(-1.029195) q[2];
sx q[2];
rz(1.8886867) q[2];
rz(-0.62475723) q[3];
sx q[3];
rz(-1.2478991) q[3];
sx q[3];
rz(-0.46318444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4598292) q[0];
sx q[0];
rz(-0.1777996) q[0];
sx q[0];
rz(2.8681927) q[0];
rz(-0.071648486) q[1];
sx q[1];
rz(-2.007808) q[1];
sx q[1];
rz(1.6260446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26165379) q[0];
sx q[0];
rz(-1.0613221) q[0];
sx q[0];
rz(-1.1306056) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6726682) q[2];
sx q[2];
rz(-1.8919626) q[2];
sx q[2];
rz(-1.4844446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3699941) q[1];
sx q[1];
rz(-0.99111667) q[1];
sx q[1];
rz(-1.7607776) q[1];
x q[2];
rz(-1.9380366) q[3];
sx q[3];
rz(-0.62588309) q[3];
sx q[3];
rz(0.82928951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7800954) q[2];
sx q[2];
rz(-2.4269035) q[2];
sx q[2];
rz(1.7717465) q[2];
rz(2.0731549) q[3];
sx q[3];
rz(-1.7504102) q[3];
sx q[3];
rz(2.1865602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2925828) q[0];
sx q[0];
rz(-2.7136901) q[0];
sx q[0];
rz(3.0191315) q[0];
rz(0.017008688) q[1];
sx q[1];
rz(-1.6288792) q[1];
sx q[1];
rz(-0.88517991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53986076) q[0];
sx q[0];
rz(-2.0575876) q[0];
sx q[0];
rz(2.9877325) q[0];
rz(-pi) q[1];
rz(-0.34220747) q[2];
sx q[2];
rz(-2.2762083) q[2];
sx q[2];
rz(-2.9660385) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7834691) q[1];
sx q[1];
rz(-1.1162288) q[1];
sx q[1];
rz(-2.697151) q[1];
rz(2.9131575) q[3];
sx q[3];
rz(-2.4266648) q[3];
sx q[3];
rz(-2.950516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.433832) q[2];
sx q[2];
rz(-1.8261352) q[2];
sx q[2];
rz(3.0641277) q[2];
rz(2.7205983) q[3];
sx q[3];
rz(-1.1869895) q[3];
sx q[3];
rz(2.8903956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0780599) q[0];
sx q[0];
rz(-0.01883004) q[0];
sx q[0];
rz(-0.68191648) q[0];
rz(-2.6497427) q[1];
sx q[1];
rz(-2.1147155) q[1];
sx q[1];
rz(-1.9815725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5321925) q[0];
sx q[0];
rz(-2.4000197) q[0];
sx q[0];
rz(1.8511008) q[0];
rz(3.0536041) q[2];
sx q[2];
rz(-0.57214979) q[2];
sx q[2];
rz(-0.83773005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3403389) q[1];
sx q[1];
rz(-1.7427708) q[1];
sx q[1];
rz(-2.1992963) q[1];
rz(1.9914636) q[3];
sx q[3];
rz(-1.3004817) q[3];
sx q[3];
rz(-0.079416954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68449768) q[2];
sx q[2];
rz(-0.86898494) q[2];
sx q[2];
rz(2.1333466) q[2];
rz(-1.6393939) q[3];
sx q[3];
rz(-2.0366171) q[3];
sx q[3];
rz(0.26386279) q[3];
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
rz(0.098792583) q[0];
sx q[0];
rz(-1.5473939) q[0];
sx q[0];
rz(-0.40476009) q[0];
rz(2.3233991) q[1];
sx q[1];
rz(-1.300756) q[1];
sx q[1];
rz(2.4249446) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7595661) q[0];
sx q[0];
rz(-1.8499825) q[0];
sx q[0];
rz(1.58103) q[0];
x q[1];
rz(2.3776182) q[2];
sx q[2];
rz(-1.171955) q[2];
sx q[2];
rz(3.0032436) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.055775) q[1];
sx q[1];
rz(-1.1113864) q[1];
sx q[1];
rz(1.8618078) q[1];
x q[2];
rz(-1.9394373) q[3];
sx q[3];
rz(-2.0799985) q[3];
sx q[3];
rz(0.64083767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8019668) q[2];
sx q[2];
rz(-2.5012987) q[2];
sx q[2];
rz(-0.65529811) q[2];
rz(2.7845434) q[3];
sx q[3];
rz(-1.3909631) q[3];
sx q[3];
rz(-0.51957875) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94012564) q[0];
sx q[0];
rz(-0.31620142) q[0];
sx q[0];
rz(1.0200208) q[0];
rz(-2.5992498) q[1];
sx q[1];
rz(-2.0261363) q[1];
sx q[1];
rz(1.7592336) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4931204) q[0];
sx q[0];
rz(-2.1628404) q[0];
sx q[0];
rz(-2.4061794) q[0];
rz(0.21127659) q[2];
sx q[2];
rz(-1.7045492) q[2];
sx q[2];
rz(-1.1115896) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0918563) q[1];
sx q[1];
rz(-1.1106792) q[1];
sx q[1];
rz(0.087398296) q[1];
rz(3.0579733) q[3];
sx q[3];
rz(-1.573054) q[3];
sx q[3];
rz(1.797054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.072784) q[2];
sx q[2];
rz(-1.6135608) q[2];
sx q[2];
rz(2.7911348) q[2];
rz(-2.6751878) q[3];
sx q[3];
rz(-0.91752183) q[3];
sx q[3];
rz(0.94356999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7826409) q[0];
sx q[0];
rz(-2.751001) q[0];
sx q[0];
rz(-0.95640957) q[0];
rz(-1.4495173) q[1];
sx q[1];
rz(-2.3608975) q[1];
sx q[1];
rz(-2.4704959) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22204493) q[0];
sx q[0];
rz(-1.4853108) q[0];
sx q[0];
rz(-1.7128904) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5726357) q[2];
sx q[2];
rz(-1.3929741) q[2];
sx q[2];
rz(-0.59019719) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.073428) q[1];
sx q[1];
rz(-0.35919093) q[1];
sx q[1];
rz(1.753162) q[1];
x q[2];
rz(1.3255886) q[3];
sx q[3];
rz(-2.2598461) q[3];
sx q[3];
rz(2.6264555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2676919) q[2];
sx q[2];
rz(-2.716422) q[2];
sx q[2];
rz(-2.0153866) q[2];
rz(-0.33024427) q[3];
sx q[3];
rz(-1.4010022) q[3];
sx q[3];
rz(0.11708524) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983593) q[0];
sx q[0];
rz(-2.5625304) q[0];
sx q[0];
rz(-0.057057127) q[0];
rz(-0.13380274) q[1];
sx q[1];
rz(-2.0963142) q[1];
sx q[1];
rz(-0.71279508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063449115) q[0];
sx q[0];
rz(-2.8196555) q[0];
sx q[0];
rz(0.19970317) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40450306) q[2];
sx q[2];
rz(-1.9484011) q[2];
sx q[2];
rz(1.9775569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5310933) q[1];
sx q[1];
rz(-1.1720017) q[1];
sx q[1];
rz(-0.51271768) q[1];
rz(-0.87656337) q[3];
sx q[3];
rz(-1.5086155) q[3];
sx q[3];
rz(2.2060195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.628525) q[2];
sx q[2];
rz(-1.2863938) q[2];
sx q[2];
rz(-3.0600424) q[2];
rz(-1.9214123) q[3];
sx q[3];
rz(-2.8704075) q[3];
sx q[3];
rz(-2.0512106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14209014) q[0];
sx q[0];
rz(-1.6195848) q[0];
sx q[0];
rz(1.5600486) q[0];
rz(2.0060495) q[1];
sx q[1];
rz(-2.0027497) q[1];
sx q[1];
rz(1.9948237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09361472) q[0];
sx q[0];
rz(-0.85913697) q[0];
sx q[0];
rz(2.467903) q[0];
rz(0.76374526) q[2];
sx q[2];
rz(-0.95852214) q[2];
sx q[2];
rz(-0.51240048) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1022328) q[1];
sx q[1];
rz(-0.22645849) q[1];
sx q[1];
rz(-0.38180085) q[1];
rz(-1.7079748) q[3];
sx q[3];
rz(-0.56932031) q[3];
sx q[3];
rz(2.5018321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5494988) q[2];
sx q[2];
rz(-1.3201951) q[2];
sx q[2];
rz(-2.7756694) q[2];
rz(-2.7538815) q[3];
sx q[3];
rz(-1.1185027) q[3];
sx q[3];
rz(-1.6754735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.53030071) q[0];
sx q[0];
rz(-0.74321754) q[0];
sx q[0];
rz(2.8331953) q[0];
rz(-2.3920234) q[1];
sx q[1];
rz(-2.0105965) q[1];
sx q[1];
rz(-1.2731193) q[1];
rz(-2.0012435) q[2];
sx q[2];
rz(-0.7067718) q[2];
sx q[2];
rz(1.5511647) q[2];
rz(2.6825873) q[3];
sx q[3];
rz(-0.68529354) q[3];
sx q[3];
rz(2.2033448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
