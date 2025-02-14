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
rz(1.5975098) q[0];
sx q[0];
rz(-1.37356) q[0];
sx q[0];
rz(1.0184259) q[0];
rz(-0.51813689) q[1];
sx q[1];
rz(-2.3845446) q[1];
sx q[1];
rz(0.63176027) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3903094) q[0];
sx q[0];
rz(-1.7770808) q[0];
sx q[0];
rz(-3.0136316) q[0];
x q[1];
rz(-0.024876923) q[2];
sx q[2];
rz(-1.8734249) q[2];
sx q[2];
rz(0.51371117) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.077191) q[1];
sx q[1];
rz(-2.4418412) q[1];
sx q[1];
rz(1.1279068) q[1];
rz(1.2445035) q[3];
sx q[3];
rz(-1.4832116) q[3];
sx q[3];
rz(-1.0911986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.53039256) q[2];
sx q[2];
rz(-2.7860614) q[2];
sx q[2];
rz(-1.4872888) q[2];
rz(-2.2330331) q[3];
sx q[3];
rz(-2.9049951) q[3];
sx q[3];
rz(-2.1366185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1210043) q[0];
sx q[0];
rz(-1.4326743) q[0];
sx q[0];
rz(2.9793136) q[0];
rz(-2.8970215) q[1];
sx q[1];
rz(-1.9385612) q[1];
sx q[1];
rz(0.29168209) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26839089) q[0];
sx q[0];
rz(-1.8531728) q[0];
sx q[0];
rz(-1.3703521) q[0];
x q[1];
rz(-0.22352085) q[2];
sx q[2];
rz(-1.2940027) q[2];
sx q[2];
rz(-1.2621438) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1202787) q[1];
sx q[1];
rz(-2.3244384) q[1];
sx q[1];
rz(-1.2517125) q[1];
x q[2];
rz(0.41145153) q[3];
sx q[3];
rz(-1.9570159) q[3];
sx q[3];
rz(-2.7038108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4633999) q[2];
sx q[2];
rz(-0.86898154) q[2];
sx q[2];
rz(2.0657067) q[2];
rz(1.894527) q[3];
sx q[3];
rz(-1.5445292) q[3];
sx q[3];
rz(-1.4472848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4397044) q[0];
sx q[0];
rz(-2.0151558) q[0];
sx q[0];
rz(2.9578748) q[0];
rz(1.5048997) q[1];
sx q[1];
rz(-2.3972062) q[1];
sx q[1];
rz(-0.16477975) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8863109) q[0];
sx q[0];
rz(-1.5613397) q[0];
sx q[0];
rz(0.9849809) q[0];
x q[1];
rz(1.8833227) q[2];
sx q[2];
rz(-2.222568) q[2];
sx q[2];
rz(1.5787909) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91080571) q[1];
sx q[1];
rz(-1.9525098) q[1];
sx q[1];
rz(-2.4934105) q[1];
rz(-pi) q[2];
rz(2.1042473) q[3];
sx q[3];
rz(-1.2743534) q[3];
sx q[3];
rz(2.8840349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.034721) q[2];
sx q[2];
rz(-1.160459) q[2];
sx q[2];
rz(-1.1064233) q[2];
rz(0.70118457) q[3];
sx q[3];
rz(-1.3283575) q[3];
sx q[3];
rz(0.4655233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3727386) q[0];
sx q[0];
rz(-2.0059858) q[0];
sx q[0];
rz(-2.1970774) q[0];
rz(1.5185897) q[1];
sx q[1];
rz(-0.98097643) q[1];
sx q[1];
rz(-1.6015582) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9409143) q[0];
sx q[0];
rz(-0.8884065) q[0];
sx q[0];
rz(2.3290578) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11941274) q[2];
sx q[2];
rz(-1.4541199) q[2];
sx q[2];
rz(-1.2774955) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51241261) q[1];
sx q[1];
rz(-0.80958074) q[1];
sx q[1];
rz(2.0960686) q[1];
x q[2];
rz(0.75042689) q[3];
sx q[3];
rz(-1.5380757) q[3];
sx q[3];
rz(-2.1145543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1314142) q[2];
sx q[2];
rz(-0.62128908) q[2];
sx q[2];
rz(2.9739001) q[2];
rz(3.0883279) q[3];
sx q[3];
rz(-1.1054509) q[3];
sx q[3];
rz(1.7648599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57529706) q[0];
sx q[0];
rz(-0.10558858) q[0];
sx q[0];
rz(2.8422624) q[0];
rz(-2.7686367) q[1];
sx q[1];
rz(-1.8920218) q[1];
sx q[1];
rz(1.6711055) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0051992) q[0];
sx q[0];
rz(-1.0233425) q[0];
sx q[0];
rz(-1.0639079) q[0];
rz(-pi) q[1];
rz(0.46385455) q[2];
sx q[2];
rz(-1.6783444) q[2];
sx q[2];
rz(0.35075089) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.3186444) q[1];
sx q[1];
rz(-0.80893436) q[1];
sx q[1];
rz(0.72973324) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8138795) q[3];
sx q[3];
rz(-2.5057) q[3];
sx q[3];
rz(-0.85195527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9386998) q[2];
sx q[2];
rz(-2.4788269) q[2];
sx q[2];
rz(-1.1104442) q[2];
rz(-0.86197305) q[3];
sx q[3];
rz(-1.5098666) q[3];
sx q[3];
rz(1.2601669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92457572) q[0];
sx q[0];
rz(-0.68704263) q[0];
sx q[0];
rz(-0.4050912) q[0];
rz(-0.336126) q[1];
sx q[1];
rz(-1.7205709) q[1];
sx q[1];
rz(-2.1902671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2971056) q[0];
sx q[0];
rz(-1.1702303) q[0];
sx q[0];
rz(1.8482918) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23299322) q[2];
sx q[2];
rz(-2.39175) q[2];
sx q[2];
rz(2.1342333) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8321633) q[1];
sx q[1];
rz(-2.1346666) q[1];
sx q[1];
rz(-1.840074) q[1];
rz(1.3780932) q[3];
sx q[3];
rz(-1.8451705) q[3];
sx q[3];
rz(1.5128795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0508017) q[2];
sx q[2];
rz(-1.4046706) q[2];
sx q[2];
rz(-0.23078272) q[2];
rz(-0.18203059) q[3];
sx q[3];
rz(-0.56806505) q[3];
sx q[3];
rz(-2.6128984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35293216) q[0];
sx q[0];
rz(-0.039529888) q[0];
sx q[0];
rz(-2.2094862) q[0];
rz(0.63356361) q[1];
sx q[1];
rz(-1.1195868) q[1];
sx q[1];
rz(1.8036141) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18620488) q[0];
sx q[0];
rz(-3.0078631) q[0];
sx q[0];
rz(-0.9132847) q[0];
x q[1];
rz(2.7510178) q[2];
sx q[2];
rz(-0.29478595) q[2];
sx q[2];
rz(0.24402555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89796126) q[1];
sx q[1];
rz(-2.1719031) q[1];
sx q[1];
rz(-2.6578085) q[1];
rz(-pi) q[2];
rz(-2.6432645) q[3];
sx q[3];
rz(-1.5428203) q[3];
sx q[3];
rz(-1.8701815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6771217) q[2];
sx q[2];
rz(-1.9147583) q[2];
sx q[2];
rz(1.202549) q[2];
rz(-0.36457148) q[3];
sx q[3];
rz(-0.69260827) q[3];
sx q[3];
rz(-0.83474368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4778336) q[0];
sx q[0];
rz(-0.57357016) q[0];
sx q[0];
rz(2.9649576) q[0];
rz(0.44003507) q[1];
sx q[1];
rz(-1.1853848) q[1];
sx q[1];
rz(2.1766591) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0394536) q[0];
sx q[0];
rz(-2.9055465) q[0];
sx q[0];
rz(2.5419767) q[0];
rz(-pi) q[1];
rz(-2.4201323) q[2];
sx q[2];
rz(-2.4000492) q[2];
sx q[2];
rz(-0.78765819) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7495887) q[1];
sx q[1];
rz(-0.67396213) q[1];
sx q[1];
rz(-2.1753009) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8550526) q[3];
sx q[3];
rz(-2.8769917) q[3];
sx q[3];
rz(-0.79597571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9289916) q[2];
sx q[2];
rz(-1.6182199) q[2];
sx q[2];
rz(-3.1316481) q[2];
rz(0.11659226) q[3];
sx q[3];
rz(-0.3392342) q[3];
sx q[3];
rz(-2.5998083) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5803439) q[0];
sx q[0];
rz(-1.9191701) q[0];
sx q[0];
rz(-2.5196581) q[0];
rz(-1.4382582) q[1];
sx q[1];
rz(-2.56918) q[1];
sx q[1];
rz(-0.66351801) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85304442) q[0];
sx q[0];
rz(-2.2877734) q[0];
sx q[0];
rz(1.2794446) q[0];
rz(0.12309317) q[2];
sx q[2];
rz(-1.8522339) q[2];
sx q[2];
rz(-1.4299453) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.864526) q[1];
sx q[1];
rz(-0.73523318) q[1];
sx q[1];
rz(1.0864054) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6648977) q[3];
sx q[3];
rz(-1.851436) q[3];
sx q[3];
rz(-0.68583127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5180987) q[2];
sx q[2];
rz(-1.3891209) q[2];
sx q[2];
rz(-2.7790879) q[2];
rz(0.47719657) q[3];
sx q[3];
rz(-2.0878744) q[3];
sx q[3];
rz(-1.1227192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.057673205) q[0];
sx q[0];
rz(-2.4102983) q[0];
sx q[0];
rz(-2.2398563) q[0];
rz(2.4665191) q[1];
sx q[1];
rz(-0.98552862) q[1];
sx q[1];
rz(2.9973082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77567277) q[0];
sx q[0];
rz(-2.7808041) q[0];
sx q[0];
rz(-1.8065288) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7117908) q[2];
sx q[2];
rz(-2.9007705) q[2];
sx q[2];
rz(-2.4483829) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.54507885) q[1];
sx q[1];
rz(-2.1987913) q[1];
sx q[1];
rz(-1.2596115) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4422632) q[3];
sx q[3];
rz(-1.398842) q[3];
sx q[3];
rz(-2.8124335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5054063) q[2];
sx q[2];
rz(-1.9289086) q[2];
sx q[2];
rz(-0.20067659) q[2];
rz(-2.0319273) q[3];
sx q[3];
rz(-2.6789013) q[3];
sx q[3];
rz(-0.5717352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5526445) q[0];
sx q[0];
rz(-2.1727967) q[0];
sx q[0];
rz(-1.1217242) q[0];
rz(0.48925346) q[1];
sx q[1];
rz(-1.4514634) q[1];
sx q[1];
rz(-1.0101752) q[1];
rz(1.7792173) q[2];
sx q[2];
rz(-2.7187612) q[2];
sx q[2];
rz(0.446212) q[2];
rz(-0.44224593) q[3];
sx q[3];
rz(-1.5968235) q[3];
sx q[3];
rz(2.1376283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
