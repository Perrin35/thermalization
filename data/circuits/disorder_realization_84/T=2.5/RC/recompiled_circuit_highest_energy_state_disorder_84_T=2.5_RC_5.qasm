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
rz(2.5166729) q[0];
sx q[0];
rz(-1.807037) q[0];
sx q[0];
rz(-3.0883489) q[0];
rz(0.48625311) q[1];
sx q[1];
rz(6.1558131) q[1];
sx q[1];
rz(11.131412) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4006179) q[0];
sx q[0];
rz(-1.8549998) q[0];
sx q[0];
rz(1.4301197) q[0];
x q[1];
rz(2.8743319) q[2];
sx q[2];
rz(-1.6858941) q[2];
sx q[2];
rz(2.5989344) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1827774) q[1];
sx q[1];
rz(-2.7441886) q[1];
sx q[1];
rz(3.0031877) q[1];
x q[2];
rz(0.49519914) q[3];
sx q[3];
rz(-0.99118587) q[3];
sx q[3];
rz(-0.51900348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3367553) q[2];
sx q[2];
rz(-1.3433604) q[2];
sx q[2];
rz(-2.8678144) q[2];
rz(-1.2182419) q[3];
sx q[3];
rz(-0.67353407) q[3];
sx q[3];
rz(-2.8555253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022920595) q[0];
sx q[0];
rz(-0.76347041) q[0];
sx q[0];
rz(2.3619695) q[0];
rz(2.4769705) q[1];
sx q[1];
rz(-1.2088935) q[1];
sx q[1];
rz(1.2453311) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98545233) q[0];
sx q[0];
rz(-0.13053556) q[0];
sx q[0];
rz(-2.0640316) q[0];
rz(-pi) q[1];
x q[1];
rz(0.092463569) q[2];
sx q[2];
rz(-2.2801054) q[2];
sx q[2];
rz(2.6576328) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43269952) q[1];
sx q[1];
rz(-1.3586167) q[1];
sx q[1];
rz(3.0794705) q[1];
rz(-pi) q[2];
rz(1.2492485) q[3];
sx q[3];
rz(-2.5587497) q[3];
sx q[3];
rz(3.1414523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43210426) q[2];
sx q[2];
rz(-0.69807845) q[2];
sx q[2];
rz(-3.091605) q[2];
rz(-1.4738119) q[3];
sx q[3];
rz(-1.4155017) q[3];
sx q[3];
rz(-0.93562359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525928) q[0];
sx q[0];
rz(-0.86243668) q[0];
sx q[0];
rz(-1.9631901) q[0];
rz(1.1475457) q[1];
sx q[1];
rz(-0.18745628) q[1];
sx q[1];
rz(0.60290927) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7258747) q[0];
sx q[0];
rz(-0.19042579) q[0];
sx q[0];
rz(-2.1384396) q[0];
x q[1];
rz(2.2608537) q[2];
sx q[2];
rz(-1.9866441) q[2];
sx q[2];
rz(-1.8772454) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29290043) q[1];
sx q[1];
rz(-2.22977) q[1];
sx q[1];
rz(1.0408569) q[1];
rz(-0.52899482) q[3];
sx q[3];
rz(-0.96312614) q[3];
sx q[3];
rz(2.0207456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9511562) q[2];
sx q[2];
rz(-0.78593212) q[2];
sx q[2];
rz(0.3832761) q[2];
rz(-1.8303309) q[3];
sx q[3];
rz(-1.0878599) q[3];
sx q[3];
rz(-0.12152984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7168147) q[0];
sx q[0];
rz(-2.1489693) q[0];
sx q[0];
rz(-0.92856652) q[0];
rz(2.3021452) q[1];
sx q[1];
rz(-1.3176368) q[1];
sx q[1];
rz(0.80926698) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033484785) q[0];
sx q[0];
rz(-1.9376457) q[0];
sx q[0];
rz(0.97292329) q[0];
rz(0.35688422) q[2];
sx q[2];
rz(-2.6864671) q[2];
sx q[2];
rz(-0.17562107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7245771) q[1];
sx q[1];
rz(-0.25522403) q[1];
sx q[1];
rz(-1.8468813) q[1];
x q[2];
rz(1.0471116) q[3];
sx q[3];
rz(-2.4274656) q[3];
sx q[3];
rz(-1.3830101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3889918) q[2];
sx q[2];
rz(-1.6099124) q[2];
sx q[2];
rz(-1.6947702) q[2];
rz(2.0145448) q[3];
sx q[3];
rz(-2.4067252) q[3];
sx q[3];
rz(-2.5126422) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9160354) q[0];
sx q[0];
rz(-1.1981523) q[0];
sx q[0];
rz(1.4300562) q[0];
rz(-2.9043708) q[1];
sx q[1];
rz(-2.1784541) q[1];
sx q[1];
rz(2.846834) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65564102) q[0];
sx q[0];
rz(-1.427934) q[0];
sx q[0];
rz(-0.44012196) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2331401) q[2];
sx q[2];
rz(-1.8690495) q[2];
sx q[2];
rz(2.3191888) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39043104) q[1];
sx q[1];
rz(-2.0195578) q[1];
sx q[1];
rz(-2.6780891) q[1];
rz(-pi) q[2];
rz(-2.9137524) q[3];
sx q[3];
rz(-0.58594847) q[3];
sx q[3];
rz(-0.70825746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99990591) q[2];
sx q[2];
rz(-0.20076951) q[2];
sx q[2];
rz(-3.0086009) q[2];
rz(2.3717132) q[3];
sx q[3];
rz(-1.8498288) q[3];
sx q[3];
rz(2.8225115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28984508) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(1.7293365) q[0];
rz(3.0899561) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(-2.9488865) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3782327) q[0];
sx q[0];
rz(-0.51412778) q[0];
sx q[0];
rz(-0.038551081) q[0];
rz(-pi) q[1];
rz(2.1837854) q[2];
sx q[2];
rz(-2.2039206) q[2];
sx q[2];
rz(-1.1867961) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6289857) q[1];
sx q[1];
rz(-0.64388212) q[1];
sx q[1];
rz(-2.3888247) q[1];
rz(1.9514378) q[3];
sx q[3];
rz(-2.1771113) q[3];
sx q[3];
rz(-2.9553138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.073033832) q[2];
sx q[2];
rz(-1.8517588) q[2];
sx q[2];
rz(-1.3812836) q[2];
rz(2.6265788) q[3];
sx q[3];
rz(-2.2541855) q[3];
sx q[3];
rz(-1.380434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9800022) q[0];
sx q[0];
rz(-1.9965633) q[0];
sx q[0];
rz(0.34647754) q[0];
rz(-2.1222291) q[1];
sx q[1];
rz(-2.4788224) q[1];
sx q[1];
rz(-1.4541218) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8603191) q[0];
sx q[0];
rz(-0.77206445) q[0];
sx q[0];
rz(2.4738929) q[0];
rz(-pi) q[1];
rz(-0.54155751) q[2];
sx q[2];
rz(-2.1514153) q[2];
sx q[2];
rz(-1.1714747) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30502637) q[1];
sx q[1];
rz(-1.8469715) q[1];
sx q[1];
rz(-2.1146875) q[1];
x q[2];
rz(3.0410513) q[3];
sx q[3];
rz(-1.7243598) q[3];
sx q[3];
rz(-2.6996149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.60160294) q[2];
sx q[2];
rz(-1.3301962) q[2];
sx q[2];
rz(2.9023564) q[2];
rz(2.7573977) q[3];
sx q[3];
rz(-2.4592063) q[3];
sx q[3];
rz(-0.95170963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20188986) q[0];
sx q[0];
rz(-1.6738142) q[0];
sx q[0];
rz(-1.3772759) q[0];
rz(1.2205623) q[1];
sx q[1];
rz(-0.86253291) q[1];
sx q[1];
rz(0.16622226) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9181053) q[0];
sx q[0];
rz(-2.0265731) q[0];
sx q[0];
rz(2.8546643) q[0];
rz(-pi) q[1];
rz(2.407309) q[2];
sx q[2];
rz(-1.3314651) q[2];
sx q[2];
rz(0.83877968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5095013) q[1];
sx q[1];
rz(-2.9281514) q[1];
sx q[1];
rz(-0.14519338) q[1];
x q[2];
rz(2.1020426) q[3];
sx q[3];
rz(-0.47787468) q[3];
sx q[3];
rz(0.22051624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7938457) q[2];
sx q[2];
rz(-2.2913427) q[2];
sx q[2];
rz(2.094685) q[2];
rz(-1.8868123) q[3];
sx q[3];
rz(-1.0898277) q[3];
sx q[3];
rz(1.2966398) q[3];
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
rz(-2.5407402) q[0];
sx q[0];
rz(-0.37373251) q[0];
sx q[0];
rz(2.0284213) q[0];
rz(2.3241849) q[1];
sx q[1];
rz(-1.4459041) q[1];
sx q[1];
rz(-2.7730952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2042646) q[0];
sx q[0];
rz(-0.44631347) q[0];
sx q[0];
rz(-0.89881368) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4442367) q[2];
sx q[2];
rz(-0.83240763) q[2];
sx q[2];
rz(2.0703933) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.011960192) q[1];
sx q[1];
rz(-2.3756177) q[1];
sx q[1];
rz(0.38928826) q[1];
x q[2];
rz(2.542607) q[3];
sx q[3];
rz(-2.3040207) q[3];
sx q[3];
rz(-0.70928228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14785279) q[2];
sx q[2];
rz(-0.61823121) q[2];
sx q[2];
rz(-1.7842133) q[2];
rz(-0.74710685) q[3];
sx q[3];
rz(-0.96304572) q[3];
sx q[3];
rz(2.4660306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6152182) q[0];
sx q[0];
rz(-0.45658657) q[0];
sx q[0];
rz(1.0427465) q[0];
rz(0.77049795) q[1];
sx q[1];
rz(-1.0105402) q[1];
sx q[1];
rz(0.76046336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7211518) q[0];
sx q[0];
rz(-1.3320702) q[0];
sx q[0];
rz(-1.9107642) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3622401) q[2];
sx q[2];
rz(-2.957666) q[2];
sx q[2];
rz(2.6435801) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4573496) q[1];
sx q[1];
rz(-1.4919229) q[1];
sx q[1];
rz(-0.42399391) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4861927) q[3];
sx q[3];
rz(-1.8701347) q[3];
sx q[3];
rz(-1.4531524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.20327917) q[2];
sx q[2];
rz(-0.61050296) q[2];
sx q[2];
rz(0.40943134) q[2];
rz(-1.5583386) q[3];
sx q[3];
rz(-2.6810472) q[3];
sx q[3];
rz(2.515365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64774491) q[0];
sx q[0];
rz(-1.0418325) q[0];
sx q[0];
rz(-2.7000725) q[0];
rz(1.582675) q[1];
sx q[1];
rz(-1.1987004) q[1];
sx q[1];
rz(-0.62335062) q[1];
rz(-1.5123488) q[2];
sx q[2];
rz(-1.6736021) q[2];
sx q[2];
rz(-2.1291669) q[2];
rz(-0.24258976) q[3];
sx q[3];
rz(-1.442083) q[3];
sx q[3];
rz(0.51221893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
