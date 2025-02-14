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
rz(2.6440808) q[0];
sx q[0];
rz(-1.3389791) q[0];
sx q[0];
rz(2.8259377) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(2.7538731) q[1];
sx q[1];
rz(8.8207689) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9785288) q[0];
sx q[0];
rz(-1.8730358) q[0];
sx q[0];
rz(0.053044293) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9004132) q[2];
sx q[2];
rz(-1.9758103) q[2];
sx q[2];
rz(0.23460282) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0993201) q[1];
sx q[1];
rz(-1.0612773) q[1];
sx q[1];
rz(-2.7541942) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2617926) q[3];
sx q[3];
rz(-1.4674392) q[3];
sx q[3];
rz(-2.7637533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.151256) q[2];
sx q[2];
rz(-1.1948816) q[2];
sx q[2];
rz(-1.7830431) q[2];
rz(-1.547706) q[3];
sx q[3];
rz(-0.93021506) q[3];
sx q[3];
rz(-1.5096629) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58594054) q[0];
sx q[0];
rz(-2.5037615) q[0];
sx q[0];
rz(-1.1974539) q[0];
rz(-0.50152957) q[1];
sx q[1];
rz(-1.5634368) q[1];
sx q[1];
rz(1.4362358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89051188) q[0];
sx q[0];
rz(-0.62411004) q[0];
sx q[0];
rz(1.9735543) q[0];
rz(-pi) q[1];
rz(1.6701489) q[2];
sx q[2];
rz(-2.2713285) q[2];
sx q[2];
rz(1.6068899) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7926678) q[1];
sx q[1];
rz(-1.3745921) q[1];
sx q[1];
rz(-1.6755107) q[1];
rz(-pi) q[2];
rz(-2.7327939) q[3];
sx q[3];
rz(-2.650947) q[3];
sx q[3];
rz(-1.4178432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2531835) q[2];
sx q[2];
rz(-1.9056355) q[2];
sx q[2];
rz(-3.1186228) q[2];
rz(-2.2394771) q[3];
sx q[3];
rz(-1.2227367) q[3];
sx q[3];
rz(-2.7149916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70812923) q[0];
sx q[0];
rz(-1.3490278) q[0];
sx q[0];
rz(0.9915114) q[0];
rz(-1.0333215) q[1];
sx q[1];
rz(-2.8585377) q[1];
sx q[1];
rz(-1.2352157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18444321) q[0];
sx q[0];
rz(-0.98688025) q[0];
sx q[0];
rz(1.8690442) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63629432) q[2];
sx q[2];
rz(-2.4264196) q[2];
sx q[2];
rz(1.4610804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0491525) q[1];
sx q[1];
rz(-1.8180646) q[1];
sx q[1];
rz(2.33806) q[1];
x q[2];
rz(-2.1622873) q[3];
sx q[3];
rz(-1.1011754) q[3];
sx q[3];
rz(-0.33794935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.027355) q[2];
sx q[2];
rz(-2.3894252) q[2];
sx q[2];
rz(-0.99119622) q[2];
rz(1.8441955) q[3];
sx q[3];
rz(-1.8810561) q[3];
sx q[3];
rz(1.4170925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39795136) q[0];
sx q[0];
rz(-2.209111) q[0];
sx q[0];
rz(-0.81746307) q[0];
rz(-2.815333) q[1];
sx q[1];
rz(-1.1810818) q[1];
sx q[1];
rz(0.78525966) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9049587) q[0];
sx q[0];
rz(-1.8371305) q[0];
sx q[0];
rz(-1.2575218) q[0];
x q[1];
rz(0.26008028) q[2];
sx q[2];
rz(-0.65979119) q[2];
sx q[2];
rz(1.545797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1568415) q[1];
sx q[1];
rz(-1.9923216) q[1];
sx q[1];
rz(2.0578029) q[1];
x q[2];
rz(2.9296285) q[3];
sx q[3];
rz(-0.3623687) q[3];
sx q[3];
rz(0.22892287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.12104812) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(-0.56309593) q[2];
rz(-0.8935039) q[3];
sx q[3];
rz(-1.9681135) q[3];
sx q[3];
rz(1.9068498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23984443) q[0];
sx q[0];
rz(-0.32308602) q[0];
sx q[0];
rz(-1.5455986) q[0];
rz(-1.0218989) q[1];
sx q[1];
rz(-2.0183759) q[1];
sx q[1];
rz(-0.18128577) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2342628) q[0];
sx q[0];
rz(-2.0600304) q[0];
sx q[0];
rz(0.2184043) q[0];
rz(-pi) q[1];
rz(1.5230122) q[2];
sx q[2];
rz(-1.0956229) q[2];
sx q[2];
rz(0.20296861) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1039155) q[1];
sx q[1];
rz(-1.5427151) q[1];
sx q[1];
rz(-1.3485391) q[1];
rz(-pi) q[2];
rz(0.44293176) q[3];
sx q[3];
rz(-1.7128403) q[3];
sx q[3];
rz(1.0012817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1229317) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(1.8956511) q[2];
rz(0.59257007) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(0.84111253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.5331921) q[0];
sx q[0];
rz(-0.99219457) q[0];
sx q[0];
rz(-2.8327508) q[0];
rz(-2.8187075) q[1];
sx q[1];
rz(-0.37728089) q[1];
sx q[1];
rz(3.0016532) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9788743) q[0];
sx q[0];
rz(-0.97467438) q[0];
sx q[0];
rz(-0.30124248) q[0];
rz(-pi) q[1];
rz(-1.762319) q[2];
sx q[2];
rz(-2.6783248) q[2];
sx q[2];
rz(-0.5199711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7227712) q[1];
sx q[1];
rz(-2.3496637) q[1];
sx q[1];
rz(-1.140711) q[1];
rz(-0.91173429) q[3];
sx q[3];
rz(-1.3305802) q[3];
sx q[3];
rz(1.6068589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66190019) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(-2.8996186) q[2];
rz(0.74609977) q[3];
sx q[3];
rz(-1.3662162) q[3];
sx q[3];
rz(1.411875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.022843) q[0];
sx q[0];
rz(-2.0976837) q[0];
sx q[0];
rz(-1.8705077) q[0];
rz(2.5785043) q[1];
sx q[1];
rz(-2.2124898) q[1];
sx q[1];
rz(-1.7005327) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6800425) q[0];
sx q[0];
rz(-0.6667887) q[0];
sx q[0];
rz(2.0186485) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3344603) q[2];
sx q[2];
rz(-2.7401667) q[2];
sx q[2];
rz(2.4106467) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5477834) q[1];
sx q[1];
rz(-1.9873706) q[1];
sx q[1];
rz(1.814247) q[1];
x q[2];
rz(1.3541) q[3];
sx q[3];
rz(-2.2842477) q[3];
sx q[3];
rz(0.59123033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9219804) q[2];
sx q[2];
rz(-1.2954804) q[2];
sx q[2];
rz(1.8249576) q[2];
rz(-2.5804139) q[3];
sx q[3];
rz(-0.96446529) q[3];
sx q[3];
rz(1.5285899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0372666) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(0.88687801) q[0];
rz(-0.2001702) q[1];
sx q[1];
rz(-1.472241) q[1];
sx q[1];
rz(-1.0221457) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3322395) q[0];
sx q[0];
rz(-0.13357559) q[0];
sx q[0];
rz(-1.3788401) q[0];
rz(-0.81424197) q[2];
sx q[2];
rz(-2.1541284) q[2];
sx q[2];
rz(-0.63177201) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31769842) q[1];
sx q[1];
rz(-0.824747) q[1];
sx q[1];
rz(-0.0024248799) q[1];
rz(-pi) q[2];
rz(-1.2194013) q[3];
sx q[3];
rz(-2.1179652) q[3];
sx q[3];
rz(1.6148293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6518121) q[2];
sx q[2];
rz(-2.688789) q[2];
sx q[2];
rz(-2.32302) q[2];
rz(-2.1583648) q[3];
sx q[3];
rz(-2.1849617) q[3];
sx q[3];
rz(0.53808588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035456903) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(1.5552833) q[0];
rz(0.69681329) q[1];
sx q[1];
rz(-0.15334829) q[1];
sx q[1];
rz(-0.78479016) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1589543) q[0];
sx q[0];
rz(-2.1110015) q[0];
sx q[0];
rz(-1.0544974) q[0];
x q[1];
rz(1.773514) q[2];
sx q[2];
rz(-1.7200816) q[2];
sx q[2];
rz(-1.210807) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6482249) q[1];
sx q[1];
rz(-0.43637143) q[1];
sx q[1];
rz(3.0619951) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56702964) q[3];
sx q[3];
rz(-2.4483557) q[3];
sx q[3];
rz(-2.6854613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67184225) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(-2.3308241) q[2];
rz(1.2189216) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(-2.0764652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3588381) q[0];
sx q[0];
rz(-2.8420119) q[0];
sx q[0];
rz(1.4240356) q[0];
rz(-2.3639823) q[1];
sx q[1];
rz(-2.168226) q[1];
sx q[1];
rz(0.75540677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1799602) q[0];
sx q[0];
rz(-2.9316219) q[0];
sx q[0];
rz(1.5927954) q[0];
rz(-pi) q[1];
rz(-0.39647409) q[2];
sx q[2];
rz(-1.586953) q[2];
sx q[2];
rz(-2.0153869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.81138071) q[1];
sx q[1];
rz(-0.63618681) q[1];
sx q[1];
rz(-0.61417555) q[1];
rz(-pi) q[2];
rz(-0.40047365) q[3];
sx q[3];
rz(-0.71610928) q[3];
sx q[3];
rz(-0.082883714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5648254) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(-2.9887065) q[2];
rz(1.8704002) q[3];
sx q[3];
rz(-1.252424) q[3];
sx q[3];
rz(2.474474) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2818903) q[0];
sx q[0];
rz(-2.6489881) q[0];
sx q[0];
rz(-0.76982605) q[0];
rz(1.9181171) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(0.6003004) q[2];
sx q[2];
rz(-2.4523198) q[2];
sx q[2];
rz(-0.92858989) q[2];
rz(-1.9742692) q[3];
sx q[3];
rz(-1.2444166) q[3];
sx q[3];
rz(1.9593596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
