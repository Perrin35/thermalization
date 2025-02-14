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
rz(-0.55627745) q[0];
sx q[0];
rz(-0.039160691) q[0];
sx q[0];
rz(-0.66443366) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(-1.4540949) q[1];
sx q[1];
rz(1.7323642) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.389654) q[0];
sx q[0];
rz(-1.03878) q[0];
sx q[0];
rz(1.3774894) q[0];
rz(-pi) q[1];
rz(-1.9033236) q[2];
sx q[2];
rz(-1.1477787) q[2];
sx q[2];
rz(3.0250383) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8027026) q[1];
sx q[1];
rz(-2.4842508) q[1];
sx q[1];
rz(2.606462) q[1];
rz(-1.9326747) q[3];
sx q[3];
rz(-1.7428977) q[3];
sx q[3];
rz(-0.73735305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0802143) q[2];
sx q[2];
rz(-0.55157101) q[2];
sx q[2];
rz(-0.75458327) q[2];
rz(1.4590229) q[3];
sx q[3];
rz(-0.85586923) q[3];
sx q[3];
rz(-1.4065019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.215312) q[0];
sx q[0];
rz(-0.54406852) q[0];
sx q[0];
rz(-1.9790443) q[0];
rz(-2.4303719) q[1];
sx q[1];
rz(-0.7178719) q[1];
sx q[1];
rz(-0.1444764) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1294589) q[0];
sx q[0];
rz(-2.8367865) q[0];
sx q[0];
rz(0.27979677) q[0];
rz(-pi) q[1];
rz(2.5108482) q[2];
sx q[2];
rz(-0.95385433) q[2];
sx q[2];
rz(-1.601905) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.45119845) q[1];
sx q[1];
rz(-1.9393801) q[1];
sx q[1];
rz(0.75779961) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9599592) q[3];
sx q[3];
rz(-0.70011052) q[3];
sx q[3];
rz(2.7896529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7084536) q[2];
sx q[2];
rz(-1.837156) q[2];
sx q[2];
rz(-2.8348095) q[2];
rz(-0.77543801) q[3];
sx q[3];
rz(-2.756835) q[3];
sx q[3];
rz(0.14498372) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48258346) q[0];
sx q[0];
rz(-2.6982396) q[0];
sx q[0];
rz(2.3495667) q[0];
rz(-1.5838985) q[1];
sx q[1];
rz(-1.7894824) q[1];
sx q[1];
rz(0.99383324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4980443) q[0];
sx q[0];
rz(-0.99731612) q[0];
sx q[0];
rz(-1.9924966) q[0];
x q[1];
rz(1.4491399) q[2];
sx q[2];
rz(-0.98441974) q[2];
sx q[2];
rz(-2.3687378) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0210304) q[1];
sx q[1];
rz(-1.1397727) q[1];
sx q[1];
rz(-2.8293508) q[1];
rz(-pi) q[2];
rz(1.1184887) q[3];
sx q[3];
rz(-1.399446) q[3];
sx q[3];
rz(0.43205827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9516248) q[2];
sx q[2];
rz(-0.8364532) q[2];
sx q[2];
rz(-2.9659029) q[2];
rz(2.9133993) q[3];
sx q[3];
rz(-0.56932813) q[3];
sx q[3];
rz(-0.5051676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93155414) q[0];
sx q[0];
rz(-0.84489548) q[0];
sx q[0];
rz(1.5616052) q[0];
rz(2.2700894) q[1];
sx q[1];
rz(-1.5305287) q[1];
sx q[1];
rz(0.16564381) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7270532) q[0];
sx q[0];
rz(-0.001984607) q[0];
sx q[0];
rz(2.5059047) q[0];
x q[1];
rz(3.13285) q[2];
sx q[2];
rz(-2.0164818) q[2];
sx q[2];
rz(2.2042556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1936744) q[1];
sx q[1];
rz(-2.3478824) q[1];
sx q[1];
rz(-3.0315184) q[1];
rz(-pi) q[2];
rz(1.8481726) q[3];
sx q[3];
rz(-2.1337389) q[3];
sx q[3];
rz(1.3514047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6988301) q[2];
sx q[2];
rz(-2.5593968) q[2];
sx q[2];
rz(0.77787918) q[2];
rz(-2.2569518) q[3];
sx q[3];
rz(-1.1529461) q[3];
sx q[3];
rz(2.1625429) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1409461) q[0];
sx q[0];
rz(-1.0362754) q[0];
sx q[0];
rz(0.86877862) q[0];
rz(0.90256214) q[1];
sx q[1];
rz(-1.2259918) q[1];
sx q[1];
rz(0.57194078) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4298852) q[0];
sx q[0];
rz(-1.0957624) q[0];
sx q[0];
rz(2.320313) q[0];
x q[1];
rz(-2.6380013) q[2];
sx q[2];
rz(-1.9707754) q[2];
sx q[2];
rz(-1.7970038) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16145591) q[1];
sx q[1];
rz(-1.2677578) q[1];
sx q[1];
rz(-1.113446) q[1];
rz(-pi) q[2];
rz(1.9492416) q[3];
sx q[3];
rz(-1.3406383) q[3];
sx q[3];
rz(-2.5251021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9945485) q[2];
sx q[2];
rz(-0.93783718) q[2];
sx q[2];
rz(0.63329548) q[2];
rz(2.5607732) q[3];
sx q[3];
rz(-1.3588901) q[3];
sx q[3];
rz(-2.4766428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13034114) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(-1.3379541) q[0];
rz(1.7300216) q[1];
sx q[1];
rz(-0.4069702) q[1];
sx q[1];
rz(-2.3445047) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6501084) q[0];
sx q[0];
rz(-2.8239125) q[0];
sx q[0];
rz(1.0450715) q[0];
rz(1.448579) q[2];
sx q[2];
rz(-2.5725992) q[2];
sx q[2];
rz(2.6560266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2458249) q[1];
sx q[1];
rz(-1.8386158) q[1];
sx q[1];
rz(2.9755249) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2615629) q[3];
sx q[3];
rz(-0.54094523) q[3];
sx q[3];
rz(0.23829392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7756614) q[2];
sx q[2];
rz(-0.66880995) q[2];
sx q[2];
rz(2.456341) q[2];
rz(0.25964409) q[3];
sx q[3];
rz(-0.97340596) q[3];
sx q[3];
rz(1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24945666) q[0];
sx q[0];
rz(-2.9849755) q[0];
sx q[0];
rz(0.53949612) q[0];
rz(0.19142137) q[1];
sx q[1];
rz(-1.4861264) q[1];
sx q[1];
rz(2.6991381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4726329) q[0];
sx q[0];
rz(-1.5711391) q[0];
sx q[0];
rz(-1.5799205) q[0];
x q[1];
rz(0.39789756) q[2];
sx q[2];
rz(-0.98491633) q[2];
sx q[2];
rz(-2.2316124) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.328622) q[1];
sx q[1];
rz(-2.1729706) q[1];
sx q[1];
rz(-2.7306836) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2768657) q[3];
sx q[3];
rz(-0.9191117) q[3];
sx q[3];
rz(0.39419277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.847984) q[2];
sx q[2];
rz(-2.831735) q[2];
sx q[2];
rz(-2.4388745) q[2];
rz(-0.27788776) q[3];
sx q[3];
rz(-1.0669471) q[3];
sx q[3];
rz(-1.3387298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24669312) q[0];
sx q[0];
rz(-2.2791857) q[0];
sx q[0];
rz(0.53584164) q[0];
rz(0.72227532) q[1];
sx q[1];
rz(-1.7012137) q[1];
sx q[1];
rz(-0.78295082) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7039258) q[0];
sx q[0];
rz(-1.9294318) q[0];
sx q[0];
rz(-3.0485247) q[0];
rz(-1.6269685) q[2];
sx q[2];
rz(-1.9675641) q[2];
sx q[2];
rz(0.53469354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4346501) q[1];
sx q[1];
rz(-2.1044113) q[1];
sx q[1];
rz(-0.60551079) q[1];
rz(-1.0231185) q[3];
sx q[3];
rz(-1.9503106) q[3];
sx q[3];
rz(-1.4364157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4202412) q[2];
sx q[2];
rz(-0.46478096) q[2];
sx q[2];
rz(2.6381524) q[2];
rz(-0.33468801) q[3];
sx q[3];
rz(-1.3379593) q[3];
sx q[3];
rz(-2.9065342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6675785) q[0];
sx q[0];
rz(-0.36822167) q[0];
sx q[0];
rz(2.0817122) q[0];
rz(0.31479752) q[1];
sx q[1];
rz(-1.7960725) q[1];
sx q[1];
rz(1.5816636) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49071872) q[0];
sx q[0];
rz(-1.3482674) q[0];
sx q[0];
rz(2.4516546) q[0];
rz(-pi) q[1];
x q[1];
rz(0.033785162) q[2];
sx q[2];
rz(-1.3063161) q[2];
sx q[2];
rz(0.95766256) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25645721) q[1];
sx q[1];
rz(-2.4607435) q[1];
sx q[1];
rz(-1.9155986) q[1];
rz(1.6208036) q[3];
sx q[3];
rz(-0.41566089) q[3];
sx q[3];
rz(-0.60189825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.18610893) q[2];
sx q[2];
rz(-1.1649705) q[2];
sx q[2];
rz(-2.1060409) q[2];
rz(1.995685) q[3];
sx q[3];
rz(-2.0290387) q[3];
sx q[3];
rz(0.15315332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.938852) q[0];
sx q[0];
rz(-0.11883141) q[0];
sx q[0];
rz(-0.76549292) q[0];
rz(2.1047523) q[1];
sx q[1];
rz(-1.8427269) q[1];
sx q[1];
rz(-0.11229215) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85690391) q[0];
sx q[0];
rz(-0.56988003) q[0];
sx q[0];
rz(2.5068552) q[0];
rz(-pi) q[1];
rz(2.5300171) q[2];
sx q[2];
rz(-0.87885746) q[2];
sx q[2];
rz(1.0146146) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72126216) q[1];
sx q[1];
rz(-2.3754076) q[1];
sx q[1];
rz(-1.2797194) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.046979172) q[3];
sx q[3];
rz(-2.3571797) q[3];
sx q[3];
rz(1.3734773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3969193) q[2];
sx q[2];
rz(-1.8645218) q[2];
sx q[2];
rz(0.12167682) q[2];
rz(-0.35950288) q[3];
sx q[3];
rz(-2.4183141) q[3];
sx q[3];
rz(3.1239037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.551238) q[0];
sx q[0];
rz(-1.4689162) q[0];
sx q[0];
rz(1.3015251) q[0];
rz(1.7474668) q[1];
sx q[1];
rz(-1.9448517) q[1];
sx q[1];
rz(1.3762884) q[1];
rz(-2.0040705) q[2];
sx q[2];
rz(-0.99356298) q[2];
sx q[2];
rz(-2.1504924) q[2];
rz(1.1351552) q[3];
sx q[3];
rz(-1.6899077) q[3];
sx q[3];
rz(-1.6514889) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
