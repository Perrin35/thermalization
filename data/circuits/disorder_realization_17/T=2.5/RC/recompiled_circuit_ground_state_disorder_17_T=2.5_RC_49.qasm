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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5296171) q[0];
sx q[0];
rz(-0.92607609) q[0];
sx q[0];
rz(2.2964649) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67439331) q[2];
sx q[2];
rz(-0.82288089) q[2];
sx q[2];
rz(-2.9059913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8583044) q[1];
sx q[1];
rz(-1.4333144) q[1];
sx q[1];
rz(1.7729574) q[1];
x q[2];
rz(1.2991139) q[3];
sx q[3];
rz(-1.6713023) q[3];
sx q[3];
rz(2.1372014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51547852) q[2];
sx q[2];
rz(-1.4592183) q[2];
sx q[2];
rz(-0.98575753) q[2];
rz(1.23729) q[3];
sx q[3];
rz(-2.3111549) q[3];
sx q[3];
rz(1.591709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80376959) q[0];
sx q[0];
rz(-1.060312) q[0];
sx q[0];
rz(-2.7983303) q[0];
rz(-0.23974165) q[1];
sx q[1];
rz(-0.38455757) q[1];
sx q[1];
rz(0.029031001) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15923289) q[0];
sx q[0];
rz(-1.3880215) q[0];
sx q[0];
rz(-0.48587783) q[0];
x q[1];
rz(-1.5324872) q[2];
sx q[2];
rz(-0.84332436) q[2];
sx q[2];
rz(3.0366355) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5482201) q[1];
sx q[1];
rz(-2.1637329) q[1];
sx q[1];
rz(-1.5520468) q[1];
x q[2];
rz(2.6695811) q[3];
sx q[3];
rz(-0.73314141) q[3];
sx q[3];
rz(2.5705119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3451738) q[2];
sx q[2];
rz(-2.5721305) q[2];
sx q[2];
rz(-0.84863895) q[2];
rz(2.7560915) q[3];
sx q[3];
rz(-1.4717088) q[3];
sx q[3];
rz(-0.29335415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.327453) q[0];
sx q[0];
rz(-1.3804945) q[0];
sx q[0];
rz(3.0019548) q[0];
rz(-2.8756554) q[1];
sx q[1];
rz(-1.1775002) q[1];
sx q[1];
rz(-1.0136484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.06482) q[0];
sx q[0];
rz(-1.0014604) q[0];
sx q[0];
rz(0.4607309) q[0];
rz(0.75516197) q[2];
sx q[2];
rz(-0.56729672) q[2];
sx q[2];
rz(3.0586363) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.35956906) q[1];
sx q[1];
rz(-1.7351362) q[1];
sx q[1];
rz(2.6604628) q[1];
rz(-pi) q[2];
rz(3.072282) q[3];
sx q[3];
rz(-1.4678174) q[3];
sx q[3];
rz(-1.5147095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2024112) q[2];
sx q[2];
rz(-2.1911669) q[2];
sx q[2];
rz(0.44516304) q[2];
rz(-0.72994453) q[3];
sx q[3];
rz(-1.688262) q[3];
sx q[3];
rz(1.2030407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79778033) q[0];
sx q[0];
rz(-2.7773363) q[0];
sx q[0];
rz(-1.181418) q[0];
rz(-1.6254609) q[1];
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
rz(-0.51633126) q[0];
sx q[0];
rz(-0.30540403) q[0];
sx q[0];
rz(2.5918079) q[0];
rz(-pi) q[1];
rz(0.67537157) q[2];
sx q[2];
rz(-0.84091917) q[2];
sx q[2];
rz(-0.60018626) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1688465) q[1];
sx q[1];
rz(-0.98727814) q[1];
sx q[1];
rz(-0.31202392) q[1];
x q[2];
rz(-0.41792385) q[3];
sx q[3];
rz(-1.897614) q[3];
sx q[3];
rz(0.55909294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2130012) q[2];
sx q[2];
rz(-1.7276305) q[2];
sx q[2];
rz(-2.6334527) q[2];
rz(-2.7022434) q[3];
sx q[3];
rz(-0.63826799) q[3];
sx q[3];
rz(-0.58258575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4947434) q[0];
sx q[0];
rz(-1.2127533) q[0];
sx q[0];
rz(-0.61432046) q[0];
rz(0.99614227) q[1];
sx q[1];
rz(-0.83729815) q[1];
sx q[1];
rz(0.070929758) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3230266) q[0];
sx q[0];
rz(-2.1686318) q[0];
sx q[0];
rz(0.37436178) q[0];
rz(-pi) q[1];
rz(-2.1251758) q[2];
sx q[2];
rz(-2.6314037) q[2];
sx q[2];
rz(-1.73908) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6560737) q[1];
sx q[1];
rz(-2.2427301) q[1];
sx q[1];
rz(0.80826064) q[1];
rz(-pi) q[2];
rz(-0.89226858) q[3];
sx q[3];
rz(-2.1315977) q[3];
sx q[3];
rz(0.73169198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4225685) q[2];
sx q[2];
rz(-1.7302128) q[2];
sx q[2];
rz(3.0651429) q[2];
rz(-3.0494173) q[3];
sx q[3];
rz(-2.01391) q[3];
sx q[3];
rz(0.60747147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344264) q[0];
sx q[0];
rz(-0.18944117) q[0];
sx q[0];
rz(-2.6763951) q[0];
rz(-3.0472962) q[1];
sx q[1];
rz(-2.1262719) q[1];
sx q[1];
rz(1.8781072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.87288) q[0];
sx q[0];
rz(-1.9221734) q[0];
sx q[0];
rz(-0.12704487) q[0];
rz(-pi) q[1];
rz(-1.5689079) q[2];
sx q[2];
rz(-0.93853653) q[2];
sx q[2];
rz(2.8304932) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81663075) q[1];
sx q[1];
rz(-1.1405319) q[1];
sx q[1];
rz(-0.63741405) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4685361) q[3];
sx q[3];
rz(-0.56922072) q[3];
sx q[3];
rz(-1.0200899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59244001) q[2];
sx q[2];
rz(-1.1391888) q[2];
sx q[2];
rz(-2.450954) q[2];
rz(1.6820827) q[3];
sx q[3];
rz(-0.032460902) q[3];
sx q[3];
rz(2.8800817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14677793) q[0];
sx q[0];
rz(-2.1187466) q[0];
sx q[0];
rz(0.47635517) q[0];
rz(-1.9626544) q[1];
sx q[1];
rz(-0.6911239) q[1];
sx q[1];
rz(-1.8144511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0020719) q[0];
sx q[0];
rz(-0.61113531) q[0];
sx q[0];
rz(2.4194761) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05169241) q[2];
sx q[2];
rz(-2.043327) q[2];
sx q[2];
rz(-0.57628265) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8031297) q[1];
sx q[1];
rz(-1.5005365) q[1];
sx q[1];
rz(-0.52609013) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61481878) q[3];
sx q[3];
rz(-1.5527627) q[3];
sx q[3];
rz(0.65476894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2771161) q[2];
sx q[2];
rz(-1.6552552) q[2];
sx q[2];
rz(0.76249301) q[2];
rz(-1.7755194) q[3];
sx q[3];
rz(-1.450918) q[3];
sx q[3];
rz(-0.33179992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45258006) q[0];
sx q[0];
rz(-3.0735425) q[0];
sx q[0];
rz(2.3633603) q[0];
rz(-1.4100086) q[1];
sx q[1];
rz(-0.86763132) q[1];
sx q[1];
rz(0.57077879) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0080514114) q[0];
sx q[0];
rz(-2.2184209) q[0];
sx q[0];
rz(-2.6106541) q[0];
rz(-1.7507872) q[2];
sx q[2];
rz(-2.8543185) q[2];
sx q[2];
rz(1.7773079) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.135317) q[1];
sx q[1];
rz(-2.0915739) q[1];
sx q[1];
rz(-3.0515677) q[1];
rz(-pi) q[2];
rz(0.024939288) q[3];
sx q[3];
rz(-0.93526269) q[3];
sx q[3];
rz(0.81944114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0290231) q[2];
sx q[2];
rz(-1.7971635) q[2];
sx q[2];
rz(0.53978187) q[2];
rz(2.6953186) q[3];
sx q[3];
rz(-0.73791426) q[3];
sx q[3];
rz(0.81587434) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096263252) q[0];
sx q[0];
rz(-0.49880609) q[0];
sx q[0];
rz(-0.86736429) q[0];
rz(1.438104) q[1];
sx q[1];
rz(-0.73355621) q[1];
sx q[1];
rz(-0.47384706) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2445969) q[0];
sx q[0];
rz(-1.8497494) q[0];
sx q[0];
rz(1.8263234) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.176193) q[2];
sx q[2];
rz(-1.2393481) q[2];
sx q[2];
rz(-1.6939577) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.20666981) q[1];
sx q[1];
rz(-2.3908472) q[1];
sx q[1];
rz(0.20845558) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51406887) q[3];
sx q[3];
rz(-0.94810728) q[3];
sx q[3];
rz(-2.4838187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4465686) q[2];
sx q[2];
rz(-2.5769672) q[2];
sx q[2];
rz(1.1308283) q[2];
rz(-1.9628149) q[3];
sx q[3];
rz(-1.3935573) q[3];
sx q[3];
rz(2.668837) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21993615) q[0];
sx q[0];
rz(-2.642785) q[0];
sx q[0];
rz(0.98861277) q[0];
rz(-2.0693208) q[1];
sx q[1];
rz(-1.5360473) q[1];
sx q[1];
rz(-1.4354717) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0804744) q[0];
sx q[0];
rz(-2.1466564) q[0];
sx q[0];
rz(0.030035069) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5631235) q[2];
sx q[2];
rz(-1.3425217) q[2];
sx q[2];
rz(1.3117737) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2769673) q[1];
sx q[1];
rz(-1.5706148) q[1];
sx q[1];
rz(-1.6305883) q[1];
rz(-2.4492757) q[3];
sx q[3];
rz(-1.8330599) q[3];
sx q[3];
rz(0.13627288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2995305) q[2];
sx q[2];
rz(-2.5008423) q[2];
sx q[2];
rz(0.93471849) q[2];
rz(-1.0457906) q[3];
sx q[3];
rz(-2.9691634) q[3];
sx q[3];
rz(-0.24833965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0228148) q[0];
sx q[0];
rz(-1.7775443) q[0];
sx q[0];
rz(-1.4015629) q[0];
rz(3.0339495) q[1];
sx q[1];
rz(-1.5168774) q[1];
sx q[1];
rz(0.37227896) q[1];
rz(1.2992181) q[2];
sx q[2];
rz(-0.80119802) q[2];
sx q[2];
rz(-0.3225633) q[2];
rz(2.6519898) q[3];
sx q[3];
rz(-2.1402749) q[3];
sx q[3];
rz(-0.28771852) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
