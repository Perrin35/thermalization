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
rz(0.82290736) q[0];
sx q[0];
rz(-0.35879254) q[0];
sx q[0];
rz(0.86831492) q[0];
rz(-1.4153642) q[1];
sx q[1];
rz(-1.1338898) q[1];
sx q[1];
rz(-0.99376065) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446237) q[0];
sx q[0];
rz(-2.7236189) q[0];
sx q[0];
rz(-2.860642) q[0];
x q[1];
rz(-0.70902087) q[2];
sx q[2];
rz(-2.5160501) q[2];
sx q[2];
rz(2.4838205) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.704761) q[1];
sx q[1];
rz(-1.1306354) q[1];
sx q[1];
rz(-2.6635567) q[1];
rz(-pi) q[2];
rz(1.078061) q[3];
sx q[3];
rz(-0.48731128) q[3];
sx q[3];
rz(-1.1307956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.913784) q[2];
sx q[2];
rz(-1.3334393) q[2];
sx q[2];
rz(-2.559973) q[2];
rz(-0.73451129) q[3];
sx q[3];
rz(-1.651265) q[3];
sx q[3];
rz(-2.7233126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.89473474) q[0];
sx q[0];
rz(-1.7623836) q[0];
sx q[0];
rz(-2.6439164) q[0];
rz(-1.0408164) q[1];
sx q[1];
rz(-2.7383995) q[1];
sx q[1];
rz(1.9042447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5631249) q[0];
sx q[0];
rz(-1.6131365) q[0];
sx q[0];
rz(2.4189831) q[0];
x q[1];
rz(-0.18816973) q[2];
sx q[2];
rz(-1.2154748) q[2];
sx q[2];
rz(-0.73588348) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9927382) q[1];
sx q[1];
rz(-1.142456) q[1];
sx q[1];
rz(2.9866339) q[1];
rz(-pi) q[2];
rz(1.7877616) q[3];
sx q[3];
rz(-1.8304018) q[3];
sx q[3];
rz(1.6093773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.441338) q[2];
sx q[2];
rz(-0.73107084) q[2];
sx q[2];
rz(2.8311484) q[2];
rz(-1.9849518) q[3];
sx q[3];
rz(-1.1247331) q[3];
sx q[3];
rz(-1.1667075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0079086) q[0];
sx q[0];
rz(-2.5569361) q[0];
sx q[0];
rz(2.7957918) q[0];
rz(2.8969104) q[1];
sx q[1];
rz(-0.9318277) q[1];
sx q[1];
rz(1.4366879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5253752) q[0];
sx q[0];
rz(-1.2572917) q[0];
sx q[0];
rz(-2.0140735) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9635003) q[2];
sx q[2];
rz(-2.1059193) q[2];
sx q[2];
rz(1.8491883) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50592283) q[1];
sx q[1];
rz(-1.3590006) q[1];
sx q[1];
rz(-0.79896547) q[1];
x q[2];
rz(-1.1967902) q[3];
sx q[3];
rz(-2.0197649) q[3];
sx q[3];
rz(-1.506492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1338542) q[2];
sx q[2];
rz(-1.1275007) q[2];
sx q[2];
rz(-0.63068843) q[2];
rz(-0.33356365) q[3];
sx q[3];
rz(-2.1400698) q[3];
sx q[3];
rz(2.1400616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0665322) q[0];
sx q[0];
rz(-2.1933031) q[0];
sx q[0];
rz(1.4045658) q[0];
rz(-2.7724077) q[1];
sx q[1];
rz(-1.7351979) q[1];
sx q[1];
rz(0.085748347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4082216) q[0];
sx q[0];
rz(-1.8693755) q[0];
sx q[0];
rz(-2.7135506) q[0];
rz(-0.46041885) q[2];
sx q[2];
rz(-1.2017837) q[2];
sx q[2];
rz(-0.93840124) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2601187) q[1];
sx q[1];
rz(-2.6131995) q[1];
sx q[1];
rz(1.0851938) q[1];
rz(-pi) q[2];
rz(-0.39601456) q[3];
sx q[3];
rz(-1.8654483) q[3];
sx q[3];
rz(-0.174232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3021476) q[2];
sx q[2];
rz(-1.2371233) q[2];
sx q[2];
rz(-2.6117924) q[2];
rz(-0.038330404) q[3];
sx q[3];
rz(-0.72961346) q[3];
sx q[3];
rz(1.5995601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2235276) q[0];
sx q[0];
rz(-2.6730838) q[0];
sx q[0];
rz(-1.202762) q[0];
rz(-2.862152) q[1];
sx q[1];
rz(-1.025082) q[1];
sx q[1];
rz(-0.7739982) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5764172) q[0];
sx q[0];
rz(-0.87067184) q[0];
sx q[0];
rz(-2.7267674) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4882795) q[2];
sx q[2];
rz(-1.8029717) q[2];
sx q[2];
rz(1.2007932) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.30602905) q[1];
sx q[1];
rz(-1.7098428) q[1];
sx q[1];
rz(-3.1285888) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1388009) q[3];
sx q[3];
rz(-1.0353902) q[3];
sx q[3];
rz(-1.1100779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1191795) q[2];
sx q[2];
rz(-0.1410307) q[2];
sx q[2];
rz(0.5640344) q[2];
rz(-2.4885079) q[3];
sx q[3];
rz(-1.1354732) q[3];
sx q[3];
rz(-2.691332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3611203) q[0];
sx q[0];
rz(-0.55279624) q[0];
sx q[0];
rz(0.64315382) q[0];
rz(1.9505352) q[1];
sx q[1];
rz(-1.4570313) q[1];
sx q[1];
rz(-2.4868884) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7793286) q[0];
sx q[0];
rz(-1.6342499) q[0];
sx q[0];
rz(0.16880798) q[0];
rz(0.046682552) q[2];
sx q[2];
rz(-2.4571223) q[2];
sx q[2];
rz(-2.4425263) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5712329) q[1];
sx q[1];
rz(-1.028299) q[1];
sx q[1];
rz(-3.065055) q[1];
rz(-pi) q[2];
rz(0.07043802) q[3];
sx q[3];
rz(-1.1813191) q[3];
sx q[3];
rz(-0.73274437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5074629) q[2];
sx q[2];
rz(-2.7644988) q[2];
sx q[2];
rz(1.4748352) q[2];
rz(0.48464388) q[3];
sx q[3];
rz(-0.99769297) q[3];
sx q[3];
rz(-1.8528329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6396879) q[0];
sx q[0];
rz(-0.26308331) q[0];
sx q[0];
rz(-2.4581773) q[0];
rz(-3.0112093) q[1];
sx q[1];
rz(-1.5905292) q[1];
sx q[1];
rz(3.108976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6382173) q[0];
sx q[0];
rz(-1.0116568) q[0];
sx q[0];
rz(2.9342164) q[0];
x q[1];
rz(-2.4559385) q[2];
sx q[2];
rz(-1.4156121) q[2];
sx q[2];
rz(2.6557166) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72717818) q[1];
sx q[1];
rz(-1.7373996) q[1];
sx q[1];
rz(1.3018621) q[1];
x q[2];
rz(-2.8835589) q[3];
sx q[3];
rz(-2.2562108) q[3];
sx q[3];
rz(-0.097921927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7183097) q[2];
sx q[2];
rz(-2.0330567) q[2];
sx q[2];
rz(-0.56524593) q[2];
rz(1.3609173) q[3];
sx q[3];
rz(-2.8748685) q[3];
sx q[3];
rz(-2.3274073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77124202) q[0];
sx q[0];
rz(-3.0841565) q[0];
sx q[0];
rz(-2.8420319) q[0];
rz(-1.6948505) q[1];
sx q[1];
rz(-0.87211496) q[1];
sx q[1];
rz(0.95132336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46630105) q[0];
sx q[0];
rz(-1.6032752) q[0];
sx q[0];
rz(1.3656653) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1323593) q[2];
sx q[2];
rz(-1.3592741) q[2];
sx q[2];
rz(2.4516425) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.044477) q[1];
sx q[1];
rz(-1.5504147) q[1];
sx q[1];
rz(2.7173032) q[1];
rz(-pi) q[2];
rz(1.2532611) q[3];
sx q[3];
rz(-2.866689) q[3];
sx q[3];
rz(1.8567228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99913725) q[2];
sx q[2];
rz(-2.3951525) q[2];
sx q[2];
rz(2.8738521) q[2];
rz(-2.1382051) q[3];
sx q[3];
rz(-2.181874) q[3];
sx q[3];
rz(0.80823922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35850152) q[0];
sx q[0];
rz(-2.6434904) q[0];
sx q[0];
rz(0.95712334) q[0];
rz(-2.762291) q[1];
sx q[1];
rz(-0.4117659) q[1];
sx q[1];
rz(-0.1651102) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50277621) q[0];
sx q[0];
rz(-0.90327016) q[0];
sx q[0];
rz(1.1291885) q[0];
x q[1];
rz(1.6696641) q[2];
sx q[2];
rz(-2.6672088) q[2];
sx q[2];
rz(2.6883467) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58453023) q[1];
sx q[1];
rz(-1.3898661) q[1];
sx q[1];
rz(-2.1541938) q[1];
rz(-pi) q[2];
rz(2.9028724) q[3];
sx q[3];
rz(-0.98084274) q[3];
sx q[3];
rz(-1.9815529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.86753201) q[2];
sx q[2];
rz(-1.910285) q[2];
sx q[2];
rz(0.77862281) q[2];
rz(0.33729956) q[3];
sx q[3];
rz(-2.1179347) q[3];
sx q[3];
rz(1.5571099) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.857665) q[0];
sx q[0];
rz(-1.090467) q[0];
sx q[0];
rz(-2.0528059) q[0];
rz(0.42824832) q[1];
sx q[1];
rz(-2.1183522) q[1];
sx q[1];
rz(0.64839378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8729432) q[0];
sx q[0];
rz(-2.6336484) q[0];
sx q[0];
rz(-0.56468876) q[0];
rz(-pi) q[1];
rz(-1.4444541) q[2];
sx q[2];
rz(-0.34046945) q[2];
sx q[2];
rz(2.7012417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2164444) q[1];
sx q[1];
rz(-0.92324644) q[1];
sx q[1];
rz(-1.5179894) q[1];
rz(0.17953486) q[3];
sx q[3];
rz(-0.57099062) q[3];
sx q[3];
rz(1.6761615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8757561) q[2];
sx q[2];
rz(-2.0529604) q[2];
sx q[2];
rz(2.9618373) q[2];
rz(-1.9836551) q[3];
sx q[3];
rz(-1.6854743) q[3];
sx q[3];
rz(1.9013083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614048) q[0];
sx q[0];
rz(-0.15459939) q[0];
sx q[0];
rz(0.79994487) q[0];
rz(-1.1663306) q[1];
sx q[1];
rz(-0.98465289) q[1];
sx q[1];
rz(-2.226895) q[1];
rz(1.8682754) q[2];
sx q[2];
rz(-1.8061721) q[2];
sx q[2];
rz(2.8190148) q[2];
rz(-0.33954444) q[3];
sx q[3];
rz(-1.3673269) q[3];
sx q[3];
rz(0.23585503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
