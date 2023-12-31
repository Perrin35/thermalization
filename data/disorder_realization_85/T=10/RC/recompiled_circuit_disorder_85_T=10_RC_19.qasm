OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.934259) q[0];
sx q[0];
rz(-0.59036314) q[0];
sx q[0];
rz(-2.7705749) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(1.376027) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5698619) q[0];
sx q[0];
rz(-0.96643448) q[0];
sx q[0];
rz(-2.5992924) q[0];
rz(-1.0797834) q[2];
sx q[2];
rz(-2.2248189) q[2];
sx q[2];
rz(-2.430928) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.260533) q[1];
sx q[1];
rz(-0.36306371) q[1];
sx q[1];
rz(-1.1313603) q[1];
rz(-pi) q[2];
rz(-0.81935482) q[3];
sx q[3];
rz(-1.5324394) q[3];
sx q[3];
rz(-1.0591782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0573037) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-0.83077103) q[3];
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
rz(-2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(-1.1799312) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(-0.72431272) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3121719) q[0];
sx q[0];
rz(-2.2342626) q[0];
sx q[0];
rz(2.7396766) q[0];
rz(-pi) q[1];
rz(-2.7472277) q[2];
sx q[2];
rz(-0.21557237) q[2];
sx q[2];
rz(2.8087316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.69246768) q[1];
sx q[1];
rz(-1.3652703) q[1];
sx q[1];
rz(-1.9301901) q[1];
x q[2];
rz(2.4507387) q[3];
sx q[3];
rz(-0.30403954) q[3];
sx q[3];
rz(-3.0512878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90536845) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-2.779707) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1419462) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(1.9225072) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8273979) q[0];
sx q[0];
rz(-2.6959246) q[0];
sx q[0];
rz(-2.218194) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9169541) q[2];
sx q[2];
rz(-2.4097754) q[2];
sx q[2];
rz(-0.099230448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29875007) q[1];
sx q[1];
rz(-2.0992273) q[1];
sx q[1];
rz(0.28038402) q[1];
x q[2];
rz(-0.015467042) q[3];
sx q[3];
rz(-1.3028212) q[3];
sx q[3];
rz(-2.1249352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7200155) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.4455618) q[2];
rz(0.56882632) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(-3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22359426) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(-2.4261684) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(-2.3148361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2430192) q[0];
sx q[0];
rz(-1.1613701) q[0];
sx q[0];
rz(-0.55222521) q[0];
rz(1.7051758) q[2];
sx q[2];
rz(-2.3339286) q[2];
sx q[2];
rz(2.813051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3863694) q[1];
sx q[1];
rz(-1.9519023) q[1];
sx q[1];
rz(-1.2632881) q[1];
x q[2];
rz(0.34927807) q[3];
sx q[3];
rz(-0.75540245) q[3];
sx q[3];
rz(-0.7790156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.15239079) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(1.7626804) q[2];
rz(-3.0692696) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(-1.4962083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82692659) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(-2.0157053) q[0];
rz(-2.2391438) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(-0.28516969) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42442214) q[0];
sx q[0];
rz(-1.5283913) q[0];
sx q[0];
rz(-1.0958584) q[0];
x q[1];
rz(-0.6638078) q[2];
sx q[2];
rz(-1.8823349) q[2];
sx q[2];
rz(-2.2500492) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73357108) q[1];
sx q[1];
rz(-1.5703778) q[1];
sx q[1];
rz(1.8838521) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22484803) q[3];
sx q[3];
rz(-2.7307011) q[3];
sx q[3];
rz(1.2332682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29331648) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(2.6021393) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-0.15701292) q[0];
rz(0.69333386) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(-1.3670115) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6002944) q[0];
sx q[0];
rz(-1.5026662) q[0];
sx q[0];
rz(3.1085204) q[0];
rz(0.81142601) q[2];
sx q[2];
rz(-1.1277414) q[2];
sx q[2];
rz(-1.4366988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1937618) q[1];
sx q[1];
rz(-1.281922) q[1];
sx q[1];
rz(0.08430251) q[1];
x q[2];
rz(-1.4671765) q[3];
sx q[3];
rz(-1.9220256) q[3];
sx q[3];
rz(2.2675089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.29629016) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(-0.40346754) q[2];
rz(2.6599595) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(-0.51923716) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24213174) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(2.2677299) q[0];
rz(-0.44772398) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(1.1706932) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255575) q[0];
sx q[0];
rz(-3.1177525) q[0];
sx q[0];
rz(2.8799879) q[0];
x q[1];
rz(-2.7295693) q[2];
sx q[2];
rz(-0.12672666) q[2];
sx q[2];
rz(0.7574946) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.88041828) q[1];
sx q[1];
rz(-2.2194127) q[1];
sx q[1];
rz(-2.3748114) q[1];
x q[2];
rz(-0.53415926) q[3];
sx q[3];
rz(-1.8338406) q[3];
sx q[3];
rz(1.1155038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0968904) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(-0.61075413) q[2];
rz(0.47510535) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(-2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(2.8444667) q[0];
rz(-1.3946474) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.811759) q[0];
sx q[0];
rz(-2.9114897) q[0];
sx q[0];
rz(2.0287201) q[0];
rz(-pi) q[1];
rz(1.4023151) q[2];
sx q[2];
rz(-2.0373166) q[2];
sx q[2];
rz(2.434935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4375953) q[1];
sx q[1];
rz(-0.62347177) q[1];
sx q[1];
rz(2.5295579) q[1];
rz(2.2046702) q[3];
sx q[3];
rz(-0.86808944) q[3];
sx q[3];
rz(-3.1261409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(0.70139766) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4421473) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(-2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(-0.10841766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.425697) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(0.9066559) q[0];
rz(1.8234532) q[2];
sx q[2];
rz(-0.18441072) q[2];
sx q[2];
rz(0.84469634) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0260967) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(-0.91067578) q[1];
x q[2];
rz(-1.7897723) q[3];
sx q[3];
rz(-0.7162381) q[3];
sx q[3];
rz(2.4718474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(0.33995315) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-0.72559124) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(2.7774096) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-0.055158786) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.508651) q[0];
sx q[0];
rz(-1.635449) q[0];
sx q[0];
rz(-2.7642194) q[0];
x q[1];
rz(1.7807547) q[2];
sx q[2];
rz(-1.0210438) q[2];
sx q[2];
rz(1.2325665) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.058051) q[1];
sx q[1];
rz(-2.1122879) q[1];
sx q[1];
rz(2.4400649) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1147898) q[3];
sx q[3];
rz(-1.5378693) q[3];
sx q[3];
rz(-1.1820716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1577592) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(0.13752078) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72538439) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(-1.4738884) q[2];
sx q[2];
rz(-1.8816392) q[2];
sx q[2];
rz(0.30602602) q[2];
rz(1.6155852) q[3];
sx q[3];
rz(-1.8106034) q[3];
sx q[3];
rz(-2.8869224) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
