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
rz(0.49993604) q[0];
sx q[0];
rz(-1.191782) q[0];
sx q[0];
rz(-1.3718104) q[0];
rz(-1.491188) q[1];
sx q[1];
rz(-2.2743382) q[1];
sx q[1];
rz(-0.39630085) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4258408) q[0];
sx q[0];
rz(-1.4914728) q[0];
sx q[0];
rz(-2.9315445) q[0];
rz(-1.9849586) q[2];
sx q[2];
rz(-0.4067308) q[2];
sx q[2];
rz(2.1106281) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1314376) q[1];
sx q[1];
rz(-2.0973296) q[1];
sx q[1];
rz(3.0326157) q[1];
rz(-1.8627152) q[3];
sx q[3];
rz(-1.1905498) q[3];
sx q[3];
rz(-2.1355262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0013915) q[2];
sx q[2];
rz(-0.80696693) q[2];
sx q[2];
rz(-1.8001455) q[2];
rz(2.0724824) q[3];
sx q[3];
rz(-2.9954698) q[3];
sx q[3];
rz(-1.7599546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57217252) q[0];
sx q[0];
rz(-2.2241346) q[0];
sx q[0];
rz(2.4826352) q[0];
rz(-3.068889) q[1];
sx q[1];
rz(-1.0560938) q[1];
sx q[1];
rz(-1.2603849) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9958268) q[0];
sx q[0];
rz(-0.46016177) q[0];
sx q[0];
rz(1.7547248) q[0];
rz(-pi) q[1];
rz(0.86671202) q[2];
sx q[2];
rz(-0.63642353) q[2];
sx q[2];
rz(1.1918024) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1998493) q[1];
sx q[1];
rz(-1.3056271) q[1];
sx q[1];
rz(-1.5035692) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96003344) q[3];
sx q[3];
rz(-1.6690147) q[3];
sx q[3];
rz(-1.4032422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2468804) q[2];
sx q[2];
rz(-2.3983045) q[2];
sx q[2];
rz(1.7903719) q[2];
rz(2.5249935) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(2.8504347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.63224822) q[0];
sx q[0];
rz(-2.6743439) q[0];
sx q[0];
rz(2.9314991) q[0];
rz(-0.083077438) q[1];
sx q[1];
rz(-1.9030842) q[1];
sx q[1];
rz(-0.067213623) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7930896) q[0];
sx q[0];
rz(-1.5843035) q[0];
sx q[0];
rz(-1.5828787) q[0];
rz(-0.24850444) q[2];
sx q[2];
rz(-1.1913158) q[2];
sx q[2];
rz(-2.7082682) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9720502) q[1];
sx q[1];
rz(-1.2360555) q[1];
sx q[1];
rz(-0.85595815) q[1];
rz(-pi) q[2];
rz(-0.10566575) q[3];
sx q[3];
rz(-0.70836954) q[3];
sx q[3];
rz(-1.8745223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.82789603) q[2];
sx q[2];
rz(-1.6500429) q[2];
sx q[2];
rz(1.776604) q[2];
rz(-2.7374173) q[3];
sx q[3];
rz(-1.3173236) q[3];
sx q[3];
rz(1.1990168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2760524) q[0];
sx q[0];
rz(-1.7409538) q[0];
sx q[0];
rz(-0.47819594) q[0];
rz(2.1887691) q[1];
sx q[1];
rz(-2.9915504) q[1];
sx q[1];
rz(2.731954) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8119558) q[0];
sx q[0];
rz(-2.4853737) q[0];
sx q[0];
rz(0.043270525) q[0];
x q[1];
rz(1.6686001) q[2];
sx q[2];
rz(-2.0253643) q[2];
sx q[2];
rz(-1.3196047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9497021) q[1];
sx q[1];
rz(-2.161945) q[1];
sx q[1];
rz(-2.3670235) q[1];
rz(-1.2495991) q[3];
sx q[3];
rz(-1.4347335) q[3];
sx q[3];
rz(0.67067671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59565583) q[2];
sx q[2];
rz(-1.7031534) q[2];
sx q[2];
rz(-1.1379918) q[2];
rz(-2.1465837) q[3];
sx q[3];
rz(-0.55750877) q[3];
sx q[3];
rz(-0.086624302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74395162) q[0];
sx q[0];
rz(-1.8331563) q[0];
sx q[0];
rz(0.27807903) q[0];
rz(1.8771578) q[1];
sx q[1];
rz(-1.8587298) q[1];
sx q[1];
rz(-1.5637195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0456384) q[0];
sx q[0];
rz(-1.4439034) q[0];
sx q[0];
rz(-1.2058099) q[0];
rz(0.78863849) q[2];
sx q[2];
rz(-1.83162) q[2];
sx q[2];
rz(-0.80718416) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.63999004) q[1];
sx q[1];
rz(-1.6899365) q[1];
sx q[1];
rz(2.5646025) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.092838959) q[3];
sx q[3];
rz(-1.412409) q[3];
sx q[3];
rz(-0.024933727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2084674) q[2];
sx q[2];
rz(-1.5310023) q[2];
sx q[2];
rz(2.7739286) q[2];
rz(-2.0465046) q[3];
sx q[3];
rz(-1.0235267) q[3];
sx q[3];
rz(-0.17525214) q[3];
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
rz(1.4472189) q[0];
sx q[0];
rz(-2.8419438) q[0];
sx q[0];
rz(-1.8023941) q[0];
rz(3.0195492) q[1];
sx q[1];
rz(-1.3587147) q[1];
sx q[1];
rz(0.36062127) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80194762) q[0];
sx q[0];
rz(-2.3390798) q[0];
sx q[0];
rz(0.54484493) q[0];
x q[1];
rz(-3.1018512) q[2];
sx q[2];
rz(-1.9905328) q[2];
sx q[2];
rz(-2.6750203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0039978) q[1];
sx q[1];
rz(-1.2892525) q[1];
sx q[1];
rz(2.0291734) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5654605) q[3];
sx q[3];
rz(-1.3179165) q[3];
sx q[3];
rz(1.3518113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.7755166) q[2];
sx q[2];
rz(-1.4630432) q[2];
sx q[2];
rz(-1.1386846) q[2];
rz(0.25389296) q[3];
sx q[3];
rz(-0.64520276) q[3];
sx q[3];
rz(2.3937288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.89895407) q[0];
sx q[0];
rz(-0.88976088) q[0];
sx q[0];
rz(2.1658072) q[0];
rz(2.0121393) q[1];
sx q[1];
rz(-0.44810805) q[1];
sx q[1];
rz(-1.7879558) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8944446) q[0];
sx q[0];
rz(-1.8631051) q[0];
sx q[0];
rz(2.4185989) q[0];
rz(-pi) q[1];
rz(1.7464306) q[2];
sx q[2];
rz(-2.7115916) q[2];
sx q[2];
rz(-1.2295837) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9043541) q[1];
sx q[1];
rz(-1.2356346) q[1];
sx q[1];
rz(-1.7101014) q[1];
rz(-pi) q[2];
rz(-0.72809345) q[3];
sx q[3];
rz(-2.545176) q[3];
sx q[3];
rz(-2.8567258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5740616) q[2];
sx q[2];
rz(-0.83063829) q[2];
sx q[2];
rz(0.91721025) q[2];
rz(1.5489102) q[3];
sx q[3];
rz(-2.8062688) q[3];
sx q[3];
rz(0.77939051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(0.2061283) q[0];
sx q[0];
rz(-1.6213106) q[0];
sx q[0];
rz(-0.024854831) q[0];
rz(-1.286233) q[1];
sx q[1];
rz(-1.7846466) q[1];
sx q[1];
rz(-1.6110427) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0162189) q[0];
sx q[0];
rz(-1.1825996) q[0];
sx q[0];
rz(-2.6632705) q[0];
rz(-2.8299061) q[2];
sx q[2];
rz(-2.4309845) q[2];
sx q[2];
rz(-0.59618261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3250585) q[1];
sx q[1];
rz(-1.5360502) q[1];
sx q[1];
rz(2.860022) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1389109) q[3];
sx q[3];
rz(-1.362097) q[3];
sx q[3];
rz(3.1155966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1207235) q[2];
sx q[2];
rz(-1.0039696) q[2];
sx q[2];
rz(1.8570159) q[2];
rz(0.32127109) q[3];
sx q[3];
rz(-2.4925241) q[3];
sx q[3];
rz(2.9200714) q[3];
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
rz(-pi/2) q[3];
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
rz(-0.016841737) q[0];
sx q[0];
rz(-1.9945972) q[0];
sx q[0];
rz(2.2220213) q[0];
rz(0.59182566) q[1];
sx q[1];
rz(-2.070919) q[1];
sx q[1];
rz(-2.8048973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3172835) q[0];
sx q[0];
rz(-1.705184) q[0];
sx q[0];
rz(1.8407673) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6918963) q[2];
sx q[2];
rz(-1.9316439) q[2];
sx q[2];
rz(2.0845513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8996657) q[1];
sx q[1];
rz(-1.4414235) q[1];
sx q[1];
rz(2.632377) q[1];
rz(2.1235663) q[3];
sx q[3];
rz(-2.8447731) q[3];
sx q[3];
rz(2.8021013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0261953) q[2];
sx q[2];
rz(-2.7112609) q[2];
sx q[2];
rz(-2.3587312) q[2];
rz(3.03249) q[3];
sx q[3];
rz(-2.1249873) q[3];
sx q[3];
rz(1.8946064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03610177) q[0];
sx q[0];
rz(-1.6721268) q[0];
sx q[0];
rz(0.92392695) q[0];
rz(-2.6149514) q[1];
sx q[1];
rz(-1.275332) q[1];
sx q[1];
rz(-0.99871666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.963388) q[0];
sx q[0];
rz(-2.9792157) q[0];
sx q[0];
rz(0.33372648) q[0];
rz(-pi) q[1];
rz(-0.082926875) q[2];
sx q[2];
rz(-1.2810263) q[2];
sx q[2];
rz(-1.7343327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6560581) q[1];
sx q[1];
rz(-1.5361934) q[1];
sx q[1];
rz(2.1587579) q[1];
x q[2];
rz(0.67855723) q[3];
sx q[3];
rz(-1.2781837) q[3];
sx q[3];
rz(-1.7998526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1363137) q[2];
sx q[2];
rz(-0.88999358) q[2];
sx q[2];
rz(-2.8078553) q[2];
rz(2.2809095) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(0.0066283289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6776047) q[0];
sx q[0];
rz(-1.1900359) q[0];
sx q[0];
rz(0.66587454) q[0];
rz(2.5869276) q[1];
sx q[1];
rz(-1.278109) q[1];
sx q[1];
rz(0.14229933) q[1];
rz(0.49665819) q[2];
sx q[2];
rz(-2.6500812) q[2];
sx q[2];
rz(-0.65190114) q[2];
rz(0.71698112) q[3];
sx q[3];
rz(-2.3752799) q[3];
sx q[3];
rz(0.3175288) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
