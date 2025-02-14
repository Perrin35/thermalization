OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.390653) q[0];
sx q[0];
rz(-0.57354623) q[0];
sx q[0];
rz(-0.030800495) q[0];
rz(-2.7544694) q[1];
sx q[1];
rz(-0.20063278) q[1];
sx q[1];
rz(-2.2320342) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46900374) q[0];
sx q[0];
rz(-1.0116972) q[0];
sx q[0];
rz(2.3536273) q[0];
rz(-pi) q[1];
rz(2.163226) q[2];
sx q[2];
rz(-2.1804951) q[2];
sx q[2];
rz(-1.1014706) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8583044) q[1];
sx q[1];
rz(-1.4333144) q[1];
sx q[1];
rz(-1.3686352) q[1];
x q[2];
rz(-1.2113282) q[3];
sx q[3];
rz(-0.28924527) q[3];
sx q[3];
rz(0.91212439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6261141) q[2];
sx q[2];
rz(-1.4592183) q[2];
sx q[2];
rz(-2.1558351) q[2];
rz(-1.9043026) q[3];
sx q[3];
rz(-0.83043778) q[3];
sx q[3];
rz(1.5498836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80376959) q[0];
sx q[0];
rz(-1.060312) q[0];
sx q[0];
rz(2.7983303) q[0];
rz(-2.901851) q[1];
sx q[1];
rz(-2.7570351) q[1];
sx q[1];
rz(-3.1125617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6343347) q[0];
sx q[0];
rz(-2.0478978) q[0];
sx q[0];
rz(1.7768589) q[0];
x q[1];
rz(0.72783651) q[2];
sx q[2];
rz(-1.542188) q[2];
sx q[2];
rz(1.7012353) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5146755) q[1];
sx q[1];
rz(-0.59319741) q[1];
sx q[1];
rz(3.1137756) q[1];
rz(1.9594479) q[3];
sx q[3];
rz(-2.209348) q[3];
sx q[3];
rz(-0.030872542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79641882) q[2];
sx q[2];
rz(-0.56946218) q[2];
sx q[2];
rz(-0.84863895) q[2];
rz(0.38550115) q[3];
sx q[3];
rz(-1.6698839) q[3];
sx q[3];
rz(2.8482385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.8141396) q[0];
sx q[0];
rz(-1.3804945) q[0];
sx q[0];
rz(0.13963786) q[0];
rz(-0.26593727) q[1];
sx q[1];
rz(-1.9640924) q[1];
sx q[1];
rz(-1.0136484) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33260298) q[0];
sx q[0];
rz(-2.4255891) q[0];
sx q[0];
rz(0.96366514) q[0];
x q[1];
rz(1.9825473) q[2];
sx q[2];
rz(-1.1687741) q[2];
sx q[2];
rz(0.92317239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2340548) q[1];
sx q[1];
rz(-0.50633303) q[1];
sx q[1];
rz(0.34409006) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.069310686) q[3];
sx q[3];
rz(-1.4678174) q[3];
sx q[3];
rz(-1.5147095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93918148) q[2];
sx q[2];
rz(-0.95042578) q[2];
sx q[2];
rz(2.6964296) q[2];
rz(-0.72994453) q[3];
sx q[3];
rz(-1.688262) q[3];
sx q[3];
rz(1.2030407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79778033) q[0];
sx q[0];
rz(-2.7773363) q[0];
sx q[0];
rz(1.181418) q[0];
rz(1.5161318) q[1];
sx q[1];
rz(-1.6226035) q[1];
sx q[1];
rz(2.6536062) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5582191) q[0];
sx q[0];
rz(-1.7285556) q[0];
sx q[0];
rz(2.8789916) q[0];
x q[1];
rz(2.4662211) q[2];
sx q[2];
rz(-0.84091917) q[2];
sx q[2];
rz(-2.5414064) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9727461) q[1];
sx q[1];
rz(-2.1543145) q[1];
sx q[1];
rz(0.31202392) q[1];
rz(-2.7236688) q[3];
sx q[3];
rz(-1.897614) q[3];
sx q[3];
rz(2.5824997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.92859149) q[2];
sx q[2];
rz(-1.7276305) q[2];
sx q[2];
rz(-0.50814) q[2];
rz(-0.43934923) q[3];
sx q[3];
rz(-0.63826799) q[3];
sx q[3];
rz(0.58258575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64684922) q[0];
sx q[0];
rz(-1.2127533) q[0];
sx q[0];
rz(0.61432046) q[0];
rz(-0.99614227) q[1];
sx q[1];
rz(-0.83729815) q[1];
sx q[1];
rz(3.0706629) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3230266) q[0];
sx q[0];
rz(-0.97296087) q[0];
sx q[0];
rz(-0.37436178) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1251758) q[2];
sx q[2];
rz(-0.51018894) q[2];
sx q[2];
rz(-1.4025127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54878753) q[1];
sx q[1];
rz(-2.1419464) q[1];
sx q[1];
rz(2.3086045) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89226858) q[3];
sx q[3];
rz(-2.1315977) q[3];
sx q[3];
rz(-0.73169198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7190242) q[2];
sx q[2];
rz(-1.7302128) q[2];
sx q[2];
rz(-3.0651429) q[2];
rz(3.0494173) q[3];
sx q[3];
rz(-2.01391) q[3];
sx q[3];
rz(2.5341212) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80716625) q[0];
sx q[0];
rz(-0.18944117) q[0];
sx q[0];
rz(2.6763951) q[0];
rz(3.0472962) q[1];
sx q[1];
rz(-2.1262719) q[1];
sx q[1];
rz(-1.8781072) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2282361) q[0];
sx q[0];
rz(-2.7688518) q[0];
sx q[0];
rz(1.9035702) q[0];
x q[1];
rz(0.63226065) q[2];
sx q[2];
rz(-1.5723197) q[2];
sx q[2];
rz(-1.8830118) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9003927) q[1];
sx q[1];
rz(-2.3896895) q[1];
sx q[1];
rz(0.65691595) q[1];
x q[2];
rz(-0.46390905) q[3];
sx q[3];
rz(-1.9134484) q[3];
sx q[3];
rz(3.1009348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59244001) q[2];
sx q[2];
rz(-1.1391888) q[2];
sx q[2];
rz(-2.450954) q[2];
rz(-1.45951) q[3];
sx q[3];
rz(-0.032460902) q[3];
sx q[3];
rz(2.8800817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14677793) q[0];
sx q[0];
rz(-2.1187466) q[0];
sx q[0];
rz(-2.6652375) q[0];
rz(-1.9626544) q[1];
sx q[1];
rz(-2.4504688) q[1];
sx q[1];
rz(-1.3271416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0020719) q[0];
sx q[0];
rz(-0.61113531) q[0];
sx q[0];
rz(0.72211653) q[0];
rz(-pi) q[1];
rz(1.671538) q[2];
sx q[2];
rz(-2.6664554) q[2];
sx q[2];
rz(0.4630928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8031297) q[1];
sx q[1];
rz(-1.6410562) q[1];
sx q[1];
rz(-0.52609013) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61481878) q[3];
sx q[3];
rz(-1.58883) q[3];
sx q[3];
rz(-0.65476894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8644766) q[2];
sx q[2];
rz(-1.4863374) q[2];
sx q[2];
rz(-0.76249301) q[2];
rz(-1.3660733) q[3];
sx q[3];
rz(-1.450918) q[3];
sx q[3];
rz(-2.8097927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6890126) q[0];
sx q[0];
rz(-0.068050139) q[0];
sx q[0];
rz(0.77823234) q[0];
rz(-1.7315841) q[1];
sx q[1];
rz(-0.86763132) q[1];
sx q[1];
rz(2.5708139) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1335412) q[0];
sx q[0];
rz(-2.2184209) q[0];
sx q[0];
rz(-2.6106541) q[0];
x q[1];
rz(-1.8536751) q[2];
sx q[2];
rz(-1.520051) q[2];
sx q[2];
rz(-0.03374781) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.135317) q[1];
sx q[1];
rz(-2.0915739) q[1];
sx q[1];
rz(-0.09002491) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2064785) q[3];
sx q[3];
rz(-1.5908656) q[3];
sx q[3];
rz(2.3754313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1125696) q[2];
sx q[2];
rz(-1.7971635) q[2];
sx q[2];
rz(-2.6018108) q[2];
rz(-0.44627407) q[3];
sx q[3];
rz(-0.73791426) q[3];
sx q[3];
rz(-2.3257183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.0453294) q[0];
sx q[0];
rz(-0.49880609) q[0];
sx q[0];
rz(-0.86736429) q[0];
rz(-1.7034886) q[1];
sx q[1];
rz(-0.73355621) q[1];
sx q[1];
rz(2.6677456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6561013) q[0];
sx q[0];
rz(-2.7655619) q[0];
sx q[0];
rz(2.4185527) q[0];
x q[1];
rz(-2.176193) q[2];
sx q[2];
rz(-1.9022446) q[2];
sx q[2];
rz(-1.4476349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6240166) q[1];
sx q[1];
rz(-1.712447) q[1];
sx q[1];
rz(2.4017548) q[1];
rz(-0.97029347) q[3];
sx q[3];
rz(-0.78509319) q[3];
sx q[3];
rz(0.11150211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69502407) q[2];
sx q[2];
rz(-2.5769672) q[2];
sx q[2];
rz(-1.1308283) q[2];
rz(1.9628149) q[3];
sx q[3];
rz(-1.3935573) q[3];
sx q[3];
rz(0.47275561) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21993615) q[0];
sx q[0];
rz(-0.49880767) q[0];
sx q[0];
rz(2.1529799) q[0];
rz(1.0722718) q[1];
sx q[1];
rz(-1.6055454) q[1];
sx q[1];
rz(1.4354717) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0804744) q[0];
sx q[0];
rz(-0.99493626) q[0];
sx q[0];
rz(-0.030035069) q[0];
x q[1];
rz(2.9133116) q[2];
sx q[2];
rz(-1.5782701) q[2];
sx q[2];
rz(0.26075903) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.86462532) q[1];
sx q[1];
rz(-1.5709779) q[1];
sx q[1];
rz(1.6305883) q[1];
rz(-2.7434968) q[3];
sx q[3];
rz(-2.4090066) q[3];
sx q[3];
rz(1.4040412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2995305) q[2];
sx q[2];
rz(-0.64075035) q[2];
sx q[2];
rz(-2.2068742) q[2];
rz(2.0958021) q[3];
sx q[3];
rz(-2.9691634) q[3];
sx q[3];
rz(-0.24833965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0228148) q[0];
sx q[0];
rz(-1.7775443) q[0];
sx q[0];
rz(-1.4015629) q[0];
rz(0.10764311) q[1];
sx q[1];
rz(-1.6247152) q[1];
sx q[1];
rz(-2.7693137) q[1];
rz(-2.3533271) q[2];
sx q[2];
rz(-1.3769281) q[2];
sx q[2];
rz(1.4396254) q[2];
rz(2.204328) q[3];
sx q[3];
rz(-2.4087003) q[3];
sx q[3];
rz(2.0743662) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
