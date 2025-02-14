OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.360541) q[0];
sx q[0];
rz(-0.6147576) q[0];
sx q[0];
rz(2.0701261) q[0];
rz(2.7331424) q[1];
sx q[1];
rz(-0.93745679) q[1];
sx q[1];
rz(1.4484922) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1030185) q[0];
sx q[0];
rz(-1.9644613) q[0];
sx q[0];
rz(-2.066924) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.229847) q[2];
sx q[2];
rz(-1.023479) q[2];
sx q[2];
rz(2.7173619) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9706124) q[1];
sx q[1];
rz(-0.17380789) q[1];
sx q[1];
rz(2.5998678) q[1];
x q[2];
rz(-2.2812649) q[3];
sx q[3];
rz(-1.2324047) q[3];
sx q[3];
rz(-0.30358728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.94413269) q[2];
sx q[2];
rz(-2.2977273) q[2];
sx q[2];
rz(0.79305631) q[2];
rz(-0.26886764) q[3];
sx q[3];
rz(-1.8157418) q[3];
sx q[3];
rz(-2.1813006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8307513) q[0];
sx q[0];
rz(-0.38250592) q[0];
sx q[0];
rz(2.261396) q[0];
rz(-0.63086069) q[1];
sx q[1];
rz(-2.3289101) q[1];
sx q[1];
rz(2.0426483) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5901075) q[0];
sx q[0];
rz(-1.5758568) q[0];
sx q[0];
rz(0.17781114) q[0];
rz(-pi) q[1];
rz(-2.1889792) q[2];
sx q[2];
rz(-1.435155) q[2];
sx q[2];
rz(2.408825) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4447282) q[1];
sx q[1];
rz(-2.676739) q[1];
sx q[1];
rz(1.9397199) q[1];
rz(-2.2078832) q[3];
sx q[3];
rz(-1.9572581) q[3];
sx q[3];
rz(-2.4368047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7546996) q[2];
sx q[2];
rz(-1.6563481) q[2];
sx q[2];
rz(-2.5999542) q[2];
rz(2.4382639) q[3];
sx q[3];
rz(-0.43916217) q[3];
sx q[3];
rz(-2.7642803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90848732) q[0];
sx q[0];
rz(-2.2623514) q[0];
sx q[0];
rz(-0.096906699) q[0];
rz(-2.5540409) q[1];
sx q[1];
rz(-2.2511626) q[1];
sx q[1];
rz(0.60120916) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22569779) q[0];
sx q[0];
rz(-0.89381274) q[0];
sx q[0];
rz(-1.3381132) q[0];
x q[1];
rz(1.5069836) q[2];
sx q[2];
rz(-1.4087965) q[2];
sx q[2];
rz(-0.29159233) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86089462) q[1];
sx q[1];
rz(-1.1105683) q[1];
sx q[1];
rz(-1.5794419) q[1];
rz(-pi) q[2];
rz(2.2279068) q[3];
sx q[3];
rz(-1.2877712) q[3];
sx q[3];
rz(1.0237657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84843695) q[2];
sx q[2];
rz(-1.6907254) q[2];
sx q[2];
rz(0.71451521) q[2];
rz(-1.4826639) q[3];
sx q[3];
rz(-0.27479333) q[3];
sx q[3];
rz(-2.7627435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42245427) q[0];
sx q[0];
rz(-1.9062573) q[0];
sx q[0];
rz(1.5041014) q[0];
rz(-2.0488886) q[1];
sx q[1];
rz(-1.8605109) q[1];
sx q[1];
rz(-0.5853931) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8849759) q[0];
sx q[0];
rz(-0.69640771) q[0];
sx q[0];
rz(-0.43088669) q[0];
rz(1.4665514) q[2];
sx q[2];
rz(-1.5603258) q[2];
sx q[2];
rz(2.143303) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3964557) q[1];
sx q[1];
rz(-1.0978699) q[1];
sx q[1];
rz(-0.84348444) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62111698) q[3];
sx q[3];
rz(-0.74291544) q[3];
sx q[3];
rz(2.7662504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66712159) q[2];
sx q[2];
rz(-2.4690364) q[2];
sx q[2];
rz(-2.7453864) q[2];
rz(2.951156) q[3];
sx q[3];
rz(-2.2021144) q[3];
sx q[3];
rz(0.52085352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025295479) q[0];
sx q[0];
rz(-2.681356) q[0];
sx q[0];
rz(-2.7468371) q[0];
rz(-2.1341628) q[1];
sx q[1];
rz(-1.1578683) q[1];
sx q[1];
rz(1.6068858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0368113) q[0];
sx q[0];
rz(-2.3320873) q[0];
sx q[0];
rz(-0.39489371) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2053648) q[2];
sx q[2];
rz(-2.2534568) q[2];
sx q[2];
rz(2.3548369) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.04567178) q[1];
sx q[1];
rz(-1.878698) q[1];
sx q[1];
rz(-0.20521693) q[1];
rz(-2.0232361) q[3];
sx q[3];
rz(-1.1471738) q[3];
sx q[3];
rz(0.37003368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69636238) q[2];
sx q[2];
rz(-0.69914377) q[2];
sx q[2];
rz(2.9635079) q[2];
rz(0.94545025) q[3];
sx q[3];
rz(-0.33792308) q[3];
sx q[3];
rz(0.43615714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6984542) q[0];
sx q[0];
rz(-0.66265023) q[0];
sx q[0];
rz(-0.15289256) q[0];
rz(-1.0180417) q[1];
sx q[1];
rz(-2.1986304) q[1];
sx q[1];
rz(-0.61000383) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40326443) q[0];
sx q[0];
rz(-2.3685507) q[0];
sx q[0];
rz(-2.9570262) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1334396) q[2];
sx q[2];
rz(-2.1446531) q[2];
sx q[2];
rz(1.2795841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78298104) q[1];
sx q[1];
rz(-0.84221101) q[1];
sx q[1];
rz(1.2549574) q[1];
x q[2];
rz(3.043706) q[3];
sx q[3];
rz(-1.6523419) q[3];
sx q[3];
rz(1.3697325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6792004) q[2];
sx q[2];
rz(-1.8588763) q[2];
sx q[2];
rz(-0.25897762) q[2];
rz(-2.9955043) q[3];
sx q[3];
rz(-0.77272213) q[3];
sx q[3];
rz(2.5025388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(1.6530957) q[0];
sx q[0];
rz(-0.86576068) q[0];
sx q[0];
rz(-0.75188941) q[0];
rz(2.409626) q[1];
sx q[1];
rz(-1.7611793) q[1];
sx q[1];
rz(-0.52925777) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2471591) q[0];
sx q[0];
rz(-0.34810796) q[0];
sx q[0];
rz(2.0624119) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8445149) q[2];
sx q[2];
rz(-2.6714795) q[2];
sx q[2];
rz(-2.672247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.098245278) q[1];
sx q[1];
rz(-1.6520513) q[1];
sx q[1];
rz(-3.1286865) q[1];
rz(2.8236709) q[3];
sx q[3];
rz(-2.0057851) q[3];
sx q[3];
rz(0.11162139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8818714) q[2];
sx q[2];
rz(-1.0292116) q[2];
sx q[2];
rz(-2.3465346) q[2];
rz(1.9184387) q[3];
sx q[3];
rz(-1.7984248) q[3];
sx q[3];
rz(2.3651626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.935598) q[0];
sx q[0];
rz(-1.5933651) q[0];
sx q[0];
rz(0.11496168) q[0];
rz(3.0523172) q[1];
sx q[1];
rz(-2.142579) q[1];
sx q[1];
rz(-0.1549391) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6892825) q[0];
sx q[0];
rz(-1.6712821) q[0];
sx q[0];
rz(1.7013676) q[0];
x q[1];
rz(0.9047382) q[2];
sx q[2];
rz(-0.78287941) q[2];
sx q[2];
rz(-2.7017186) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8304883) q[1];
sx q[1];
rz(-1.1837256) q[1];
sx q[1];
rz(1.8482659) q[1];
x q[2];
rz(-3.0691419) q[3];
sx q[3];
rz(-0.84699291) q[3];
sx q[3];
rz(-0.75954306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0838919) q[2];
sx q[2];
rz(-1.8996779) q[2];
sx q[2];
rz(-0.6883626) q[2];
rz(-0.56811959) q[3];
sx q[3];
rz(-1.0024242) q[3];
sx q[3];
rz(-2.4585371) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883009) q[0];
sx q[0];
rz(-0.81128565) q[0];
sx q[0];
rz(1.9360833) q[0];
rz(1.0909117) q[1];
sx q[1];
rz(-0.58735192) q[1];
sx q[1];
rz(-1.6039414) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6610049) q[0];
sx q[0];
rz(-1.2086226) q[0];
sx q[0];
rz(1.1042547) q[0];
x q[1];
rz(1.5700266) q[2];
sx q[2];
rz(-2.5736897) q[2];
sx q[2];
rz(-2.9132194) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.007215167) q[1];
sx q[1];
rz(-2.341801) q[1];
sx q[1];
rz(-0.15534955) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5439369) q[3];
sx q[3];
rz(-1.171805) q[3];
sx q[3];
rz(1.1433218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.18498147) q[2];
sx q[2];
rz(-1.1143755) q[2];
sx q[2];
rz(-1.9967009) q[2];
rz(1.6635118) q[3];
sx q[3];
rz(-2.3749115) q[3];
sx q[3];
rz(-1.112282) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0612653) q[0];
sx q[0];
rz(-0.63617951) q[0];
sx q[0];
rz(-1.7058477) q[0];
rz(-1.1371293) q[1];
sx q[1];
rz(-2.4500193) q[1];
sx q[1];
rz(-0.93708509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6244753) q[0];
sx q[0];
rz(-1.269425) q[0];
sx q[0];
rz(1.6624438) q[0];
x q[1];
rz(0.50855277) q[2];
sx q[2];
rz(-1.9483231) q[2];
sx q[2];
rz(3.0005531) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.62793362) q[1];
sx q[1];
rz(-0.80027831) q[1];
sx q[1];
rz(-1.5350128) q[1];
rz(-pi) q[2];
rz(-0.61599515) q[3];
sx q[3];
rz(-1.8503891) q[3];
sx q[3];
rz(3.0699025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3252141) q[2];
sx q[2];
rz(-0.64322317) q[2];
sx q[2];
rz(-0.31036672) q[2];
rz(0.54272932) q[3];
sx q[3];
rz(-2.2118745) q[3];
sx q[3];
rz(-0.31699666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19001374) q[0];
sx q[0];
rz(-2.0484476) q[0];
sx q[0];
rz(0.53566326) q[0];
rz(-0.71470064) q[1];
sx q[1];
rz(-1.8938046) q[1];
sx q[1];
rz(-1.6751777) q[1];
rz(-2.705639) q[2];
sx q[2];
rz(-1.9212848) q[2];
sx q[2];
rz(-2.3876847) q[2];
rz(1.5101931) q[3];
sx q[3];
rz(-0.84280673) q[3];
sx q[3];
rz(-1.6731586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
