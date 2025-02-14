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
rz(-2.812204) q[0];
sx q[0];
rz(-0.63151276) q[0];
sx q[0];
rz(3.1407177) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(5.2823245) q[1];
sx q[1];
rz(9.7699788) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787155) q[0];
sx q[0];
rz(-1.5460617) q[0];
sx q[0];
rz(1.6616761) q[0];
rz(-pi) q[1];
rz(-1.7183185) q[2];
sx q[2];
rz(-1.194724) q[2];
sx q[2];
rz(-0.12898117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.4243489) q[1];
sx q[1];
rz(-2.1645438) q[1];
sx q[1];
rz(0.51148606) q[1];
rz(-pi) q[2];
rz(-1.7854858) q[3];
sx q[3];
rz(-0.076138894) q[3];
sx q[3];
rz(-2.3977181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(0.20797569) q[2];
rz(-0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(-1.1408898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8107373) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(-2.148707) q[0];
rz(1.8244686) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(2.3670926) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3178783) q[0];
sx q[0];
rz(-1.5597938) q[0];
sx q[0];
rz(1.3925793) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6235145) q[2];
sx q[2];
rz(-2.2642914) q[2];
sx q[2];
rz(0.77502807) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5238873) q[1];
sx q[1];
rz(-1.7951709) q[1];
sx q[1];
rz(-2.7588506) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22529545) q[3];
sx q[3];
rz(-1.1114128) q[3];
sx q[3];
rz(2.6090906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2449067) q[2];
sx q[2];
rz(-1.8550355) q[2];
sx q[2];
rz(-2.1265105) q[2];
rz(0.24909881) q[3];
sx q[3];
rz(-0.86779147) q[3];
sx q[3];
rz(0.36756137) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86431137) q[0];
sx q[0];
rz(-2.4503777) q[0];
sx q[0];
rz(2.9840898) q[0];
rz(-2.1306254) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(-0.13883042) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15181825) q[0];
sx q[0];
rz(-0.4954557) q[0];
sx q[0];
rz(-0.83924967) q[0];
rz(-pi) q[1];
rz(-1.1558258) q[2];
sx q[2];
rz(-2.2425644) q[2];
sx q[2];
rz(-2.6443554) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5032387) q[1];
sx q[1];
rz(-1.8625096) q[1];
sx q[1];
rz(2.9364763) q[1];
rz(-pi) q[2];
rz(1.3957455) q[3];
sx q[3];
rz(-1.850046) q[3];
sx q[3];
rz(-1.6498914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46131721) q[2];
sx q[2];
rz(-2.2291144) q[2];
sx q[2];
rz(0.3581363) q[2];
rz(-0.67874587) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(1.0391191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0963652) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(-2.8367693) q[0];
rz(0.40799704) q[1];
sx q[1];
rz(-1.7034737) q[1];
sx q[1];
rz(2.2136484) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6154895) q[0];
sx q[0];
rz(-0.23357059) q[0];
sx q[0];
rz(-2.156267) q[0];
rz(-pi) q[1];
rz(0.81176968) q[2];
sx q[2];
rz(-0.78321811) q[2];
sx q[2];
rz(-2.4350172) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4907704) q[1];
sx q[1];
rz(-0.70413744) q[1];
sx q[1];
rz(-2.9068374) q[1];
rz(-1.0545066) q[3];
sx q[3];
rz(-2.1355503) q[3];
sx q[3];
rz(-1.3206583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1912332) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(-2.9534269) q[2];
rz(1.865271) q[3];
sx q[3];
rz(-1.3151582) q[3];
sx q[3];
rz(1.7108542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6998049) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(3.1222043) q[0];
rz(-2.0806606) q[1];
sx q[1];
rz(-1.7273936) q[1];
sx q[1];
rz(1.9020938) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1000308) q[0];
sx q[0];
rz(-1.8503555) q[0];
sx q[0];
rz(2.1651405) q[0];
x q[1];
rz(-0.6940191) q[2];
sx q[2];
rz(-0.87022129) q[2];
sx q[2];
rz(-1.7050336) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.46249786) q[1];
sx q[1];
rz(-2.0789008) q[1];
sx q[1];
rz(1.8930045) q[1];
x q[2];
rz(2.7321759) q[3];
sx q[3];
rz(-2.0338661) q[3];
sx q[3];
rz(-0.45209979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1218607) q[2];
sx q[2];
rz(-1.3132361) q[2];
sx q[2];
rz(1.3266374) q[2];
rz(-2.895368) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(-0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5354079) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(-0.73915172) q[0];
rz(-2.7958561) q[1];
sx q[1];
rz(-1.812499) q[1];
sx q[1];
rz(-2.2703222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5369878) q[0];
sx q[0];
rz(-2.3571627) q[0];
sx q[0];
rz(-3.0434199) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61001444) q[2];
sx q[2];
rz(-2.2918309) q[2];
sx q[2];
rz(-0.80243669) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4193486) q[1];
sx q[1];
rz(-0.39545317) q[1];
sx q[1];
rz(1.761318) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63251791) q[3];
sx q[3];
rz(-2.4261203) q[3];
sx q[3];
rz(-2.4533437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.31200108) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(0.87257067) q[2];
rz(1.5729337) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(-0.93305552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0912112) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(-2.3799489) q[0];
rz(-0.42539445) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(0.92686191) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0480014) q[0];
sx q[0];
rz(-0.33080745) q[0];
sx q[0];
rz(-0.76865102) q[0];
rz(-pi) q[1];
rz(2.2144187) q[2];
sx q[2];
rz(-2.1774051) q[2];
sx q[2];
rz(-1.3869029) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95185223) q[1];
sx q[1];
rz(-1.2092672) q[1];
sx q[1];
rz(0.23747634) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6377444) q[3];
sx q[3];
rz(-2.5173325) q[3];
sx q[3];
rz(0.74893307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1130134) q[2];
sx q[2];
rz(-0.69602746) q[2];
sx q[2];
rz(-2.2704303) q[2];
rz(-0.29843676) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(-1.8106073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.4153862) q[0];
sx q[0];
rz(-0.50006777) q[0];
sx q[0];
rz(0.4183847) q[0];
rz(-2.8437974) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(-0.13430886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9878085) q[0];
sx q[0];
rz(-2.1001108) q[0];
sx q[0];
rz(2.3321926) q[0];
x q[1];
rz(1.9842582) q[2];
sx q[2];
rz(-0.43817876) q[2];
sx q[2];
rz(-2.7587121) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.011439104) q[1];
sx q[1];
rz(-1.4362037) q[1];
sx q[1];
rz(0.18421872) q[1];
x q[2];
rz(2.8577096) q[3];
sx q[3];
rz(-1.1177269) q[3];
sx q[3];
rz(2.4456152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(0.64201391) q[2];
rz(2.0783453) q[3];
sx q[3];
rz(-2.4898873) q[3];
sx q[3];
rz(-2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5409656) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(-3.0294321) q[0];
rz(1.3345831) q[1];
sx q[1];
rz(-0.59252512) q[1];
sx q[1];
rz(2.1122011) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063113041) q[0];
sx q[0];
rz(-1.6029583) q[0];
sx q[0];
rz(0.86634791) q[0];
rz(-pi) q[1];
rz(-2.8376828) q[2];
sx q[2];
rz(-1.5329156) q[2];
sx q[2];
rz(1.607995) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9255213) q[1];
sx q[1];
rz(-0.23434429) q[1];
sx q[1];
rz(-1.7611124) q[1];
rz(-pi) q[2];
rz(-0.46799) q[3];
sx q[3];
rz(-1.6525998) q[3];
sx q[3];
rz(-1.8379267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7182497) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(-3.1035799) q[2];
rz(2.5293317) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(-0.14122252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24899471) q[0];
sx q[0];
rz(-1.6706415) q[0];
sx q[0];
rz(-2.7227962) q[0];
rz(-1.8839802) q[1];
sx q[1];
rz(-1.5092756) q[1];
sx q[1];
rz(0.14258252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2217956) q[0];
sx q[0];
rz(-2.9807973) q[0];
sx q[0];
rz(0.17531403) q[0];
x q[1];
rz(1.8480186) q[2];
sx q[2];
rz(-0.96053329) q[2];
sx q[2];
rz(-2.0137613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8293162) q[1];
sx q[1];
rz(-1.8354776) q[1];
sx q[1];
rz(0.97977248) q[1];
rz(-pi) q[2];
rz(-2.0311277) q[3];
sx q[3];
rz(-1.2602196) q[3];
sx q[3];
rz(-2.6643857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7964145) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(2.8734015) q[2];
rz(-1.4043407) q[3];
sx q[3];
rz(-2.0495448) q[3];
sx q[3];
rz(1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0161229) q[0];
sx q[0];
rz(-2.7653427) q[0];
sx q[0];
rz(0.34633037) q[0];
rz(-3.012433) q[1];
sx q[1];
rz(-1.2551413) q[1];
sx q[1];
rz(-1.7358949) q[1];
rz(-1.2186981) q[2];
sx q[2];
rz(-2.0205971) q[2];
sx q[2];
rz(-1.481075) q[2];
rz(1.3553452) q[3];
sx q[3];
rz(-0.30134311) q[3];
sx q[3];
rz(-2.338196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
