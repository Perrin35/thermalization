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
rz(-0.62491971) q[0];
sx q[0];
rz(-1.3345557) q[0];
sx q[0];
rz(-0.053243756) q[0];
rz(0.48625311) q[1];
sx q[1];
rz(6.1558131) q[1];
sx q[1];
rz(11.131412) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4006179) q[0];
sx q[0];
rz(-1.8549998) q[0];
sx q[0];
rz(1.7114729) q[0];
rz(1.6900914) q[2];
sx q[2];
rz(-1.8362459) q[2];
sx q[2];
rz(-2.1448898) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8088805) q[1];
sx q[1];
rz(-1.96419) q[1];
sx q[1];
rz(1.5129526) q[1];
x q[2];
rz(-0.49519914) q[3];
sx q[3];
rz(-2.1504068) q[3];
sx q[3];
rz(2.6225892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80483738) q[2];
sx q[2];
rz(-1.7982322) q[2];
sx q[2];
rz(2.8678144) q[2];
rz(1.9233507) q[3];
sx q[3];
rz(-2.4680586) q[3];
sx q[3];
rz(2.8555253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.022920595) q[0];
sx q[0];
rz(-2.3781222) q[0];
sx q[0];
rz(2.3619695) q[0];
rz(-0.66462213) q[1];
sx q[1];
rz(-1.2088935) q[1];
sx q[1];
rz(-1.8962616) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6593334) q[0];
sx q[0];
rz(-1.6856992) q[0];
sx q[0];
rz(3.0795133) q[0];
x q[1];
rz(-0.092463569) q[2];
sx q[2];
rz(-2.2801054) q[2];
sx q[2];
rz(-2.6576328) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14549637) q[1];
sx q[1];
rz(-2.9206373) q[1];
sx q[1];
rz(1.851382) q[1];
rz(-pi) q[2];
rz(-0.20540463) q[3];
sx q[3];
rz(-2.1202135) q[3];
sx q[3];
rz(0.37976625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7094884) q[2];
sx q[2];
rz(-2.4435142) q[2];
sx q[2];
rz(-0.049987642) q[2];
rz(1.4738119) q[3];
sx q[3];
rz(-1.4155017) q[3];
sx q[3];
rz(-2.2059691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1525928) q[0];
sx q[0];
rz(-0.86243668) q[0];
sx q[0];
rz(-1.1784026) q[0];
rz(1.1475457) q[1];
sx q[1];
rz(-0.18745628) q[1];
sx q[1];
rz(0.60290927) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5459203) q[0];
sx q[0];
rz(-1.4688558) q[0];
sx q[0];
rz(1.7319182) q[0];
rz(-pi) q[1];
rz(2.2608537) q[2];
sx q[2];
rz(-1.1549486) q[2];
sx q[2];
rz(1.8772454) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29290043) q[1];
sx q[1];
rz(-0.91182263) q[1];
sx q[1];
rz(2.1007358) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94305493) q[3];
sx q[3];
rz(-2.3585417) q[3];
sx q[3];
rz(2.8179226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19043645) q[2];
sx q[2];
rz(-2.3556605) q[2];
sx q[2];
rz(2.7583165) q[2];
rz(1.3112618) q[3];
sx q[3];
rz(-2.0537328) q[3];
sx q[3];
rz(-3.0200628) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424778) q[0];
sx q[0];
rz(-0.99262339) q[0];
sx q[0];
rz(-2.2130261) q[0];
rz(-2.3021452) q[1];
sx q[1];
rz(-1.3176368) q[1];
sx q[1];
rz(2.3323257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8438551) q[0];
sx q[0];
rz(-1.017573) q[0];
sx q[0];
rz(-0.4351625) q[0];
rz(-0.35688422) q[2];
sx q[2];
rz(-0.45512558) q[2];
sx q[2];
rz(2.9659716) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41701554) q[1];
sx q[1];
rz(-2.8863686) q[1];
sx q[1];
rz(-1.8468813) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2146642) q[3];
sx q[3];
rz(-1.2371105) q[3];
sx q[3];
rz(2.9179171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75260085) q[2];
sx q[2];
rz(-1.5316803) q[2];
sx q[2];
rz(-1.4468225) q[2];
rz(2.0145448) q[3];
sx q[3];
rz(-0.73486745) q[3];
sx q[3];
rz(-0.62895044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22555722) q[0];
sx q[0];
rz(-1.9434403) q[0];
sx q[0];
rz(-1.7115364) q[0];
rz(0.23722181) q[1];
sx q[1];
rz(-2.1784541) q[1];
sx q[1];
rz(-0.29475862) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65564102) q[0];
sx q[0];
rz(-1.427934) q[0];
sx q[0];
rz(2.7014707) q[0];
x q[1];
rz(2.7698603) q[2];
sx q[2];
rz(-0.9424389) q[2];
sx q[2];
rz(-2.6184751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39043104) q[1];
sx q[1];
rz(-2.0195578) q[1];
sx q[1];
rz(0.46350355) q[1];
rz(1.421991) q[3];
sx q[3];
rz(-1.0019046) q[3];
sx q[3];
rz(-0.97967813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99990591) q[2];
sx q[2];
rz(-2.9408231) q[2];
sx q[2];
rz(-0.13299175) q[2];
rz(-0.76987949) q[3];
sx q[3];
rz(-1.2917638) q[3];
sx q[3];
rz(-2.8225115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28984508) q[0];
sx q[0];
rz(-2.802749) q[0];
sx q[0];
rz(-1.7293365) q[0];
rz(-3.0899561) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(2.9488865) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9826003) q[0];
sx q[0];
rz(-1.5518414) q[0];
sx q[0];
rz(-0.51380957) q[0];
x q[1];
rz(-0.66484837) q[2];
sx q[2];
rz(-0.85077219) q[2];
sx q[2];
rz(0.31535044) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51260692) q[1];
sx q[1];
rz(-0.64388212) q[1];
sx q[1];
rz(2.3888247) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5000743) q[3];
sx q[3];
rz(-1.2605485) q[3];
sx q[3];
rz(-1.1603242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0685588) q[2];
sx q[2];
rz(-1.2898338) q[2];
sx q[2];
rz(1.3812836) q[2];
rz(-0.51501385) q[3];
sx q[3];
rz(-2.2541855) q[3];
sx q[3];
rz(-1.380434) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16159049) q[0];
sx q[0];
rz(-1.9965633) q[0];
sx q[0];
rz(0.34647754) q[0];
rz(-2.1222291) q[1];
sx q[1];
rz(-0.66277021) q[1];
sx q[1];
rz(1.4541218) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0270868) q[0];
sx q[0];
rz(-0.99106228) q[0];
sx q[0];
rz(-2.113335) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2367661) q[2];
sx q[2];
rz(-0.77205333) q[2];
sx q[2];
rz(1.1386629) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84216181) q[1];
sx q[1];
rz(-0.60363942) q[1];
sx q[1];
rz(-1.0697212) q[1];
x q[2];
rz(3.0410513) q[3];
sx q[3];
rz(-1.7243598) q[3];
sx q[3];
rz(0.44197772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5399897) q[2];
sx q[2];
rz(-1.3301962) q[2];
sx q[2];
rz(-0.23923624) q[2];
rz(2.7573977) q[3];
sx q[3];
rz(-2.4592063) q[3];
sx q[3];
rz(2.189883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.20188986) q[0];
sx q[0];
rz(-1.6738142) q[0];
sx q[0];
rz(-1.7643167) q[0];
rz(1.2205623) q[1];
sx q[1];
rz(-2.2790597) q[1];
sx q[1];
rz(-0.16622226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3275422) q[0];
sx q[0];
rz(-0.53314994) q[0];
sx q[0];
rz(1.0472121) q[0];
rz(-pi) q[1];
rz(1.8883838) q[2];
sx q[2];
rz(-0.86195213) q[2];
sx q[2];
rz(2.6203757) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6580171) q[1];
sx q[1];
rz(-1.7819575) q[1];
sx q[1];
rz(-1.6021452) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1020426) q[3];
sx q[3];
rz(-0.47787468) q[3];
sx q[3];
rz(2.9210764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.347747) q[2];
sx q[2];
rz(-0.85024992) q[2];
sx q[2];
rz(-1.0469077) q[2];
rz(-1.2547803) q[3];
sx q[3];
rz(-2.051765) q[3];
sx q[3];
rz(1.2966398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5407402) q[0];
sx q[0];
rz(-0.37373251) q[0];
sx q[0];
rz(1.1131713) q[0];
rz(0.81740776) q[1];
sx q[1];
rz(-1.4459041) q[1];
sx q[1];
rz(-0.36849749) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21462003) q[0];
sx q[0];
rz(-1.2262188) q[0];
sx q[0];
rz(-0.2895245) q[0];
x q[1];
rz(-3.0037874) q[2];
sx q[2];
rz(-2.3944582) q[2];
sx q[2];
rz(-1.2580308) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8702183) q[1];
sx q[1];
rz(-1.8370359) q[1];
sx q[1];
rz(-0.72721796) q[1];
x q[2];
rz(0.59898563) q[3];
sx q[3];
rz(-2.3040207) q[3];
sx q[3];
rz(-2.4323104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9937399) q[2];
sx q[2];
rz(-0.61823121) q[2];
sx q[2];
rz(-1.3573793) q[2];
rz(-2.3944858) q[3];
sx q[3];
rz(-0.96304572) q[3];
sx q[3];
rz(0.67556206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6152182) q[0];
sx q[0];
rz(-2.6850061) q[0];
sx q[0];
rz(1.0427465) q[0];
rz(0.77049795) q[1];
sx q[1];
rz(-1.0105402) q[1];
sx q[1];
rz(0.76046336) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805394) q[0];
sx q[0];
rz(-0.4127316) q[0];
sx q[0];
rz(0.94032918) q[0];
rz(-pi) q[1];
rz(-1.7793525) q[2];
sx q[2];
rz(-0.18392662) q[2];
sx q[2];
rz(2.6435801) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14899602) q[1];
sx q[1];
rz(-1.1482052) q[1];
sx q[1];
rz(-1.6572957) q[1];
rz(-pi) q[2];
rz(-1.1995777) q[3];
sx q[3];
rz(-0.94915024) q[3];
sx q[3];
rz(0.10528639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20327917) q[2];
sx q[2];
rz(-0.61050296) q[2];
sx q[2];
rz(-2.7321613) q[2];
rz(-1.5832541) q[3];
sx q[3];
rz(-0.46054545) q[3];
sx q[3];
rz(2.515365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4938477) q[0];
sx q[0];
rz(-2.0997601) q[0];
sx q[0];
rz(0.44152015) q[0];
rz(-1.582675) q[1];
sx q[1];
rz(-1.9428923) q[1];
sx q[1];
rz(2.518242) q[1];
rz(0.10298038) q[2];
sx q[2];
rz(-1.5126577) q[2];
sx q[2];
rz(2.589227) q[2];
rz(-2.6473896) q[3];
sx q[3];
rz(-2.8675666) q[3];
sx q[3];
rz(1.6045229) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
