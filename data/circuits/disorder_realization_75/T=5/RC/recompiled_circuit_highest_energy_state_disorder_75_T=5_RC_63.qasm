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
rz(1.3156112) q[0];
sx q[0];
rz(-0.82807461) q[0];
sx q[0];
rz(-0.84870422) q[0];
rz(1.4915713) q[1];
sx q[1];
rz(-0.68869156) q[1];
sx q[1];
rz(-2.7149849) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2899982) q[0];
sx q[0];
rz(-1.6713872) q[0];
sx q[0];
rz(-1.6466584) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83580221) q[2];
sx q[2];
rz(-0.15087946) q[2];
sx q[2];
rz(-2.7250233) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6411031) q[1];
sx q[1];
rz(-1.4979384) q[1];
sx q[1];
rz(-0.44568731) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95152609) q[3];
sx q[3];
rz(-2.2689156) q[3];
sx q[3];
rz(-1.8474471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2171057) q[2];
sx q[2];
rz(-1.7542398) q[2];
sx q[2];
rz(3.1021049) q[2];
rz(1.0314137) q[3];
sx q[3];
rz(-1.02905) q[3];
sx q[3];
rz(1.2379117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51204387) q[0];
sx q[0];
rz(-1.6144253) q[0];
sx q[0];
rz(-1.7426096) q[0];
rz(1.6276739) q[1];
sx q[1];
rz(-2.2986423) q[1];
sx q[1];
rz(-1.7248076) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88089529) q[0];
sx q[0];
rz(-1.3492378) q[0];
sx q[0];
rz(-0.76789121) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8292839) q[2];
sx q[2];
rz(-0.80855364) q[2];
sx q[2];
rz(-1.1201348) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21525225) q[1];
sx q[1];
rz(-1.7071586) q[1];
sx q[1];
rz(2.7424314) q[1];
rz(0.11059304) q[3];
sx q[3];
rz(-0.40990007) q[3];
sx q[3];
rz(1.5701587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.17025718) q[2];
sx q[2];
rz(-2.038326) q[2];
sx q[2];
rz(-0.7106759) q[2];
rz(-3.1273048) q[3];
sx q[3];
rz(-2.894214) q[3];
sx q[3];
rz(1.3361196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.22064848) q[0];
sx q[0];
rz(-2.1823688) q[0];
sx q[0];
rz(0.4162108) q[0];
rz(-1.2843708) q[1];
sx q[1];
rz(-1.4764079) q[1];
sx q[1];
rz(0.31825569) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98781768) q[0];
sx q[0];
rz(-1.1889699) q[0];
sx q[0];
rz(0.29819684) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65223434) q[2];
sx q[2];
rz(-2.1822737) q[2];
sx q[2];
rz(-2.230951) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7075338) q[1];
sx q[1];
rz(-2.4726082) q[1];
sx q[1];
rz(0.21978746) q[1];
x q[2];
rz(-0.89983799) q[3];
sx q[3];
rz(-0.38649118) q[3];
sx q[3];
rz(2.6725566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2661813) q[2];
sx q[2];
rz(-0.77989945) q[2];
sx q[2];
rz(-1.2904588) q[2];
rz(-3.087888) q[3];
sx q[3];
rz(-1.8219681) q[3];
sx q[3];
rz(-1.3857589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5829492) q[0];
sx q[0];
rz(-2.4531893) q[0];
sx q[0];
rz(2.131856) q[0];
rz(-0.91521493) q[1];
sx q[1];
rz(-1.5292294) q[1];
sx q[1];
rz(0.42542747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4998326) q[0];
sx q[0];
rz(-1.4041889) q[0];
sx q[0];
rz(-0.84979041) q[0];
rz(-pi) q[1];
rz(0.073510344) q[2];
sx q[2];
rz(-0.82907721) q[2];
sx q[2];
rz(2.4501767) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.47290651) q[1];
sx q[1];
rz(-2.2232409) q[1];
sx q[1];
rz(-0.5202867) q[1];
rz(-1.7336044) q[3];
sx q[3];
rz(-1.5763487) q[3];
sx q[3];
rz(-0.83205637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2414744) q[2];
sx q[2];
rz(-1.590531) q[2];
sx q[2];
rz(-0.11108622) q[2];
rz(2.9491718) q[3];
sx q[3];
rz(-2.7399053) q[3];
sx q[3];
rz(-0.69453159) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38839328) q[0];
sx q[0];
rz(-0.73047262) q[0];
sx q[0];
rz(2.2764192) q[0];
rz(-2.6490037) q[1];
sx q[1];
rz(-1.0961696) q[1];
sx q[1];
rz(1.385484) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139112) q[0];
sx q[0];
rz(-1.7075141) q[0];
sx q[0];
rz(2.3522931) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0086914) q[2];
sx q[2];
rz(-2.3302493) q[2];
sx q[2];
rz(-1.0824301) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3802478) q[1];
sx q[1];
rz(-1.963539) q[1];
sx q[1];
rz(-0.56403807) q[1];
rz(-pi) q[2];
rz(0.82204865) q[3];
sx q[3];
rz(-2.1072142) q[3];
sx q[3];
rz(2.3182403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8869141) q[2];
sx q[2];
rz(-2.6723537) q[2];
sx q[2];
rz(2.4349507) q[2];
rz(-0.32736579) q[3];
sx q[3];
rz(-1.2923765) q[3];
sx q[3];
rz(-2.4842026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3051598) q[0];
sx q[0];
rz(-1.7921472) q[0];
sx q[0];
rz(2.7850372) q[0];
rz(-0.56218475) q[1];
sx q[1];
rz(-1.1327344) q[1];
sx q[1];
rz(2.1551989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0092457744) q[0];
sx q[0];
rz(-0.86902394) q[0];
sx q[0];
rz(2.470507) q[0];
rz(2.3190034) q[2];
sx q[2];
rz(-2.6554567) q[2];
sx q[2];
rz(-1.8784472) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74979177) q[1];
sx q[1];
rz(-1.5448678) q[1];
sx q[1];
rz(0.75849979) q[1];
rz(-pi) q[2];
rz(-0.18137698) q[3];
sx q[3];
rz(-1.3311609) q[3];
sx q[3];
rz(-0.58663034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3761313) q[2];
sx q[2];
rz(-0.71530801) q[2];
sx q[2];
rz(-2.4410655) q[2];
rz(-0.96945196) q[3];
sx q[3];
rz(-2.1078096) q[3];
sx q[3];
rz(-0.9303003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9483865) q[0];
sx q[0];
rz(-0.32873118) q[0];
sx q[0];
rz(-2.9902003) q[0];
rz(1.5104177) q[1];
sx q[1];
rz(-2.0163586) q[1];
sx q[1];
rz(1.6815394) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6918644) q[0];
sx q[0];
rz(-2.1546493) q[0];
sx q[0];
rz(-2.1644781) q[0];
rz(-pi) q[1];
rz(-2.0978155) q[2];
sx q[2];
rz(-1.8879381) q[2];
sx q[2];
rz(-2.4316367) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5826539) q[1];
sx q[1];
rz(-13/(6*pi)) q[1];
sx q[1];
rz(0.93780545) q[1];
rz(-1.3760819) q[3];
sx q[3];
rz(-1.7140183) q[3];
sx q[3];
rz(0.16370521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2363362) q[2];
sx q[2];
rz(-2.4589804) q[2];
sx q[2];
rz(0.35161099) q[2];
rz(2.6998399) q[3];
sx q[3];
rz(-0.92919246) q[3];
sx q[3];
rz(2.9898804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.041895954) q[0];
sx q[0];
rz(-1.8459039) q[0];
sx q[0];
rz(-1.1610485) q[0];
rz(-2.4940122) q[1];
sx q[1];
rz(-1.3787965) q[1];
sx q[1];
rz(-2.3191998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1402991) q[0];
sx q[0];
rz(-1.1736794) q[0];
sx q[0];
rz(-2.1236093) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0089125) q[2];
sx q[2];
rz(-0.88399502) q[2];
sx q[2];
rz(-2.7504454) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6127738) q[1];
sx q[1];
rz(-2.3102323) q[1];
sx q[1];
rz(-0.39503132) q[1];
rz(-pi) q[2];
x q[2];
rz(1.137072) q[3];
sx q[3];
rz(-1.2081283) q[3];
sx q[3];
rz(-2.7421302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9245727) q[2];
sx q[2];
rz(-1.9968888) q[2];
sx q[2];
rz(-1.8355969) q[2];
rz(-0.71803391) q[3];
sx q[3];
rz(-0.82457232) q[3];
sx q[3];
rz(-0.048862783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43593916) q[0];
sx q[0];
rz(-1.0834563) q[0];
sx q[0];
rz(0.017024592) q[0];
rz(-1.9450933) q[1];
sx q[1];
rz(-0.39533177) q[1];
sx q[1];
rz(-0.33504018) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8340613) q[0];
sx q[0];
rz(-2.1387707) q[0];
sx q[0];
rz(-0.012004367) q[0];
x q[1];
rz(2.4990988) q[2];
sx q[2];
rz(-1.2351002) q[2];
sx q[2];
rz(3.1106126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0373203) q[1];
sx q[1];
rz(-1.3824711) q[1];
sx q[1];
rz(0.89481797) q[1];
rz(1.0862971) q[3];
sx q[3];
rz(-1.1346176) q[3];
sx q[3];
rz(-1.0624113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7177141) q[2];
sx q[2];
rz(-0.7462036) q[2];
sx q[2];
rz(-0.27900532) q[2];
rz(-1.3689857) q[3];
sx q[3];
rz(-1.6266581) q[3];
sx q[3];
rz(-0.92931187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245215) q[0];
sx q[0];
rz(-2.9947424) q[0];
sx q[0];
rz(0.76706925) q[0];
rz(-0.37297878) q[1];
sx q[1];
rz(-1.6423128) q[1];
sx q[1];
rz(-2.9708718) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6304008) q[0];
sx q[0];
rz(-2.3053279) q[0];
sx q[0];
rz(-1.6096576) q[0];
rz(-2.7880048) q[2];
sx q[2];
rz(-0.43904009) q[2];
sx q[2];
rz(-3.0692284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.062791477) q[1];
sx q[1];
rz(-1.2578221) q[1];
sx q[1];
rz(2.880135) q[1];
rz(-pi) q[2];
rz(0.3091273) q[3];
sx q[3];
rz(-2.6539475) q[3];
sx q[3];
rz(-1.6548827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79147044) q[2];
sx q[2];
rz(-0.56362027) q[2];
sx q[2];
rz(0.0027837022) q[2];
rz(2.2400098) q[3];
sx q[3];
rz(-2.1193347) q[3];
sx q[3];
rz(-2.3367052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0881385) q[0];
sx q[0];
rz(-2.8202941) q[0];
sx q[0];
rz(1.1711076) q[0];
rz(2.3114655) q[1];
sx q[1];
rz(-1.3273888) q[1];
sx q[1];
rz(-1.8524016) q[1];
rz(2.6302677) q[2];
sx q[2];
rz(-2.1187388) q[2];
sx q[2];
rz(-0.50730898) q[2];
rz(-0.8030025) q[3];
sx q[3];
rz(-2.0661125) q[3];
sx q[3];
rz(-2.6281602) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
