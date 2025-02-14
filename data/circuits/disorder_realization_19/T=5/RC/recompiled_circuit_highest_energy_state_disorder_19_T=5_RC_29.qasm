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
rz(-2.6416566) q[0];
sx q[0];
rz(-1.9498107) q[0];
sx q[0];
rz(-1.7697822) q[0];
rz(-1.491188) q[1];
sx q[1];
rz(-2.2743382) q[1];
sx q[1];
rz(-0.39630085) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7157519) q[0];
sx q[0];
rz(-1.6501199) q[0];
sx q[0];
rz(-2.9315445) q[0];
rz(0.17163817) q[2];
sx q[2];
rz(-1.2002103) q[2];
sx q[2];
rz(2.5570392) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5057136) q[1];
sx q[1];
rz(-1.4766271) q[1];
sx q[1];
rz(-2.0999184) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39537786) q[3];
sx q[3];
rz(-1.8413183) q[3];
sx q[3];
rz(0.45365712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1402011) q[2];
sx q[2];
rz(-0.80696693) q[2];
sx q[2];
rz(-1.8001455) q[2];
rz(1.0691102) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(-1.7599546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5694201) q[0];
sx q[0];
rz(-0.91745806) q[0];
sx q[0];
rz(0.65895748) q[0];
rz(-0.072703687) q[1];
sx q[1];
rz(-1.0560938) q[1];
sx q[1];
rz(1.2603849) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35048553) q[0];
sx q[0];
rz(-2.022615) q[0];
sx q[0];
rz(0.090403948) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86671202) q[2];
sx q[2];
rz(-0.63642353) q[2];
sx q[2];
rz(-1.1918024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5301828) q[1];
sx q[1];
rz(-1.5059222) q[1];
sx q[1];
rz(-2.8758509) q[1];
rz(-pi) q[2];
rz(0.96003344) q[3];
sx q[3];
rz(-1.6690147) q[3];
sx q[3];
rz(-1.7383505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2468804) q[2];
sx q[2];
rz(-2.3983045) q[2];
sx q[2];
rz(-1.7903719) q[2];
rz(-2.5249935) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(0.29115796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63224822) q[0];
sx q[0];
rz(-0.46724874) q[0];
sx q[0];
rz(2.9314991) q[0];
rz(-3.0585152) q[1];
sx q[1];
rz(-1.9030842) q[1];
sx q[1];
rz(-3.074379) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34850304) q[0];
sx q[0];
rz(-1.5572892) q[0];
sx q[0];
rz(1.5587139) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8930882) q[2];
sx q[2];
rz(-1.9502769) q[2];
sx q[2];
rz(-2.7082682) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3785079) q[1];
sx q[1];
rz(-0.77662599) q[1];
sx q[1];
rz(2.0586527) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4806858) q[3];
sx q[3];
rz(-0.86719524) q[3];
sx q[3];
rz(-1.4058324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3136966) q[2];
sx q[2];
rz(-1.6500429) q[2];
sx q[2];
rz(-1.3649887) q[2];
rz(-0.40417534) q[3];
sx q[3];
rz(-1.3173236) q[3];
sx q[3];
rz(1.9425758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86554027) q[0];
sx q[0];
rz(-1.4006389) q[0];
sx q[0];
rz(-2.6633967) q[0];
rz(0.95282355) q[1];
sx q[1];
rz(-2.9915504) q[1];
sx q[1];
rz(0.40963867) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.206868) q[0];
sx q[0];
rz(-1.544401) q[0];
sx q[0];
rz(-0.65576632) q[0];
rz(1.4729926) q[2];
sx q[2];
rz(-1.1162283) q[2];
sx q[2];
rz(-1.3196047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9497021) q[1];
sx q[1];
rz(-0.97964761) q[1];
sx q[1];
rz(-0.77456919) q[1];
x q[2];
rz(-2.9982938) q[3];
sx q[3];
rz(-1.2526726) q[3];
sx q[3];
rz(-2.1963726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5459368) q[2];
sx q[2];
rz(-1.4384392) q[2];
sx q[2];
rz(1.1379918) q[2];
rz(-0.99500895) q[3];
sx q[3];
rz(-0.55750877) q[3];
sx q[3];
rz(-3.0549684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74395162) q[0];
sx q[0];
rz(-1.3084363) q[0];
sx q[0];
rz(-0.27807903) q[0];
rz(1.2644348) q[1];
sx q[1];
rz(-1.2828628) q[1];
sx q[1];
rz(-1.5637195) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57347689) q[0];
sx q[0];
rz(-1.9327123) q[0];
sx q[0];
rz(-0.135735) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78863849) q[2];
sx q[2];
rz(-1.83162) q[2];
sx q[2];
rz(-2.3344085) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.63999004) q[1];
sx q[1];
rz(-1.6899365) q[1];
sx q[1];
rz(-0.57699012) q[1];
x q[2];
rz(-3.0487537) q[3];
sx q[3];
rz(-1.7291837) q[3];
sx q[3];
rz(-0.024933727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2084674) q[2];
sx q[2];
rz(-1.5310023) q[2];
sx q[2];
rz(0.36766407) q[2];
rz(1.095088) q[3];
sx q[3];
rz(-2.118066) q[3];
sx q[3];
rz(-2.9663405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4472189) q[0];
sx q[0];
rz(-0.29964888) q[0];
sx q[0];
rz(-1.8023941) q[0];
rz(0.12204349) q[1];
sx q[1];
rz(-1.782878) q[1];
sx q[1];
rz(0.36062127) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084753732) q[0];
sx q[0];
rz(-2.2331616) q[0];
sx q[0];
rz(1.0785036) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1507658) q[2];
sx q[2];
rz(-1.6070865) q[2];
sx q[2];
rz(-2.0535713) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0039978) q[1];
sx q[1];
rz(-1.2892525) q[1];
sx q[1];
rz(1.1124193) q[1];
rz(-pi) q[2];
rz(2.8887094) q[3];
sx q[3];
rz(-1.5759625) q[3];
sx q[3];
rz(-0.21765003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3660761) q[2];
sx q[2];
rz(-1.4630432) q[2];
sx q[2];
rz(-1.1386846) q[2];
rz(0.25389296) q[3];
sx q[3];
rz(-2.4963899) q[3];
sx q[3];
rz(-2.3937288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.89895407) q[0];
sx q[0];
rz(-0.88976088) q[0];
sx q[0];
rz(2.1658072) q[0];
rz(-2.0121393) q[1];
sx q[1];
rz(-0.44810805) q[1];
sx q[1];
rz(1.7879558) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8944446) q[0];
sx q[0];
rz(-1.2784875) q[0];
sx q[0];
rz(-0.72299374) q[0];
x q[1];
rz(-1.1466402) q[2];
sx q[2];
rz(-1.4978906) q[2];
sx q[2];
rz(2.6404501) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28747121) q[1];
sx q[1];
rz(-1.7023037) q[1];
sx q[1];
rz(-2.8033957) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9951172) q[3];
sx q[3];
rz(-2.0034308) q[3];
sx q[3];
rz(2.603765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5740616) q[2];
sx q[2];
rz(-2.3109544) q[2];
sx q[2];
rz(2.2243824) q[2];
rz(-1.5926825) q[3];
sx q[3];
rz(-0.33532381) q[3];
sx q[3];
rz(-0.77939051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9354644) q[0];
sx q[0];
rz(-1.6213106) q[0];
sx q[0];
rz(-3.1167378) q[0];
rz(-1.8553597) q[1];
sx q[1];
rz(-1.7846466) q[1];
sx q[1];
rz(-1.53055) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0162189) q[0];
sx q[0];
rz(-1.1825996) q[0];
sx q[0];
rz(-2.6632705) q[0];
rz(-pi) q[1];
rz(-2.4553005) q[2];
sx q[2];
rz(-1.7721906) q[2];
sx q[2];
rz(-1.9274595) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3250585) q[1];
sx q[1];
rz(-1.5360502) q[1];
sx q[1];
rz(-2.860022) q[1];
rz(1.7814662) q[3];
sx q[3];
rz(-1.7066741) q[3];
sx q[3];
rz(-1.6257515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.020869104) q[2];
sx q[2];
rz(-2.1376231) q[2];
sx q[2];
rz(-1.2845767) q[2];
rz(0.32127109) q[3];
sx q[3];
rz(-0.64906859) q[3];
sx q[3];
rz(0.2215213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1247509) q[0];
sx q[0];
rz(-1.9945972) q[0];
sx q[0];
rz(0.91957134) q[0];
rz(0.59182566) q[1];
sx q[1];
rz(-1.0706736) q[1];
sx q[1];
rz(-0.33669534) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3172835) q[0];
sx q[0];
rz(-1.705184) q[0];
sx q[0];
rz(1.8407673) q[0];
rz(-0.30981242) q[2];
sx q[2];
rz(-0.37978077) q[2];
sx q[2];
rz(-1.3889695) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8996657) q[1];
sx q[1];
rz(-1.4414235) q[1];
sx q[1];
rz(2.632377) q[1];
rz(-pi) q[2];
x q[2];
rz(1.82545) q[3];
sx q[3];
rz(-1.7249722) q[3];
sx q[3];
rz(-1.3772937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0261953) q[2];
sx q[2];
rz(-0.43033174) q[2];
sx q[2];
rz(-2.3587312) q[2];
rz(-0.10910263) q[3];
sx q[3];
rz(-1.0166054) q[3];
sx q[3];
rz(1.2469863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1054909) q[0];
sx q[0];
rz(-1.4694659) q[0];
sx q[0];
rz(2.2176657) q[0];
rz(2.6149514) q[1];
sx q[1];
rz(-1.275332) q[1];
sx q[1];
rz(0.99871666) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0629421) q[0];
sx q[0];
rz(-1.5178158) q[0];
sx q[0];
rz(-0.15356252) q[0];
rz(1.2998215) q[2];
sx q[2];
rz(-0.30108115) q[2];
sx q[2];
rz(-1.6903413) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6560581) q[1];
sx q[1];
rz(-1.6053992) q[1];
sx q[1];
rz(0.98283474) q[1];
rz(-pi) q[2];
rz(-1.9400334) q[3];
sx q[3];
rz(-2.2155361) q[3];
sx q[3];
rz(0.45756868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0052789) q[2];
sx q[2];
rz(-2.2515991) q[2];
sx q[2];
rz(-2.8078553) q[2];
rz(0.8606832) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(3.1349643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6776047) q[0];
sx q[0];
rz(-1.1900359) q[0];
sx q[0];
rz(0.66587454) q[0];
rz(0.55466501) q[1];
sx q[1];
rz(-1.8634836) q[1];
sx q[1];
rz(-2.9992933) q[1];
rz(1.8205504) q[2];
sx q[2];
rz(-1.142923) q[2];
sx q[2];
rz(-0.10071071) q[2];
rz(2.4246115) q[3];
sx q[3];
rz(-0.76631279) q[3];
sx q[3];
rz(-2.8240639) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
