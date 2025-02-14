OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.40795657) q[0];
sx q[0];
rz(-2.9562558) q[0];
sx q[0];
rz(1.6428525) q[0];
rz(-0.19417956) q[1];
sx q[1];
rz(-0.3436389) q[1];
sx q[1];
rz(2.763881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35750264) q[0];
sx q[0];
rz(-1.329485) q[0];
sx q[0];
rz(-1.3426128) q[0];
x q[1];
rz(0.80773662) q[2];
sx q[2];
rz(-2.5662072) q[2];
sx q[2];
rz(-2.6859716) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9895674) q[1];
sx q[1];
rz(-1.5366239) q[1];
sx q[1];
rz(2.4844643) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1123766) q[3];
sx q[3];
rz(-1.3061285) q[3];
sx q[3];
rz(2.9801369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3898042) q[2];
sx q[2];
rz(-1.7916388) q[2];
sx q[2];
rz(-0.26503116) q[2];
rz(0.42936471) q[3];
sx q[3];
rz(-1.8931484) q[3];
sx q[3];
rz(-0.00092367729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65381831) q[0];
sx q[0];
rz(-2.9980897) q[0];
sx q[0];
rz(1.300746) q[0];
rz(0.067151345) q[1];
sx q[1];
rz(-0.90922272) q[1];
sx q[1];
rz(0.86004177) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22500998) q[0];
sx q[0];
rz(-2.2448178) q[0];
sx q[0];
rz(-2.3300578) q[0];
rz(-pi) q[1];
rz(2.4323787) q[2];
sx q[2];
rz(-2.1265891) q[2];
sx q[2];
rz(1.1467939) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7262267) q[1];
sx q[1];
rz(-2.4176855) q[1];
sx q[1];
rz(0.93127802) q[1];
x q[2];
rz(-1.2705401) q[3];
sx q[3];
rz(-1.7466063) q[3];
sx q[3];
rz(2.6838357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5106875) q[2];
sx q[2];
rz(-0.61669934) q[2];
sx q[2];
rz(2.5751233) q[2];
rz(-2.412292) q[3];
sx q[3];
rz(-0.84435487) q[3];
sx q[3];
rz(-0.20310371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1931964) q[0];
sx q[0];
rz(-2.0553135) q[0];
sx q[0];
rz(-2.7753944) q[0];
rz(-1.5858448) q[1];
sx q[1];
rz(-1.0428753) q[1];
sx q[1];
rz(3.0911456) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18264601) q[0];
sx q[0];
rz(-3.052127) q[0];
sx q[0];
rz(-2.0561809) q[0];
rz(-pi) q[1];
rz(-2.0642199) q[2];
sx q[2];
rz(-2.1972547) q[2];
sx q[2];
rz(-2.3979417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55712426) q[1];
sx q[1];
rz(-1.4038175) q[1];
sx q[1];
rz(-2.9191769) q[1];
rz(-pi) q[2];
rz(1.8837711) q[3];
sx q[3];
rz(-0.44383263) q[3];
sx q[3];
rz(-2.7601506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0744276) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(-1.830843) q[2];
rz(-1.7835167) q[3];
sx q[3];
rz(-0.57099968) q[3];
sx q[3];
rz(1.3091492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686789) q[0];
sx q[0];
rz(-2.9153115) q[0];
sx q[0];
rz(1.3847466) q[0];
rz(-1.2387431) q[1];
sx q[1];
rz(-1.9107995) q[1];
sx q[1];
rz(1.1741656) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0525111) q[0];
sx q[0];
rz(-1.3818372) q[0];
sx q[0];
rz(1.1367537) q[0];
rz(-pi) q[1];
rz(-2.410589) q[2];
sx q[2];
rz(-2.7766697) q[2];
sx q[2];
rz(2.475955) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2432855) q[1];
sx q[1];
rz(-1.2445868) q[1];
sx q[1];
rz(1.0350569) q[1];
rz(-0.45278182) q[3];
sx q[3];
rz(-2.2119129) q[3];
sx q[3];
rz(2.6124817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7063286) q[2];
sx q[2];
rz(-1.3202983) q[2];
sx q[2];
rz(-0.038979385) q[2];
rz(-2.4987761) q[3];
sx q[3];
rz(-2.1233852) q[3];
sx q[3];
rz(2.5207998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83127999) q[0];
sx q[0];
rz(-1.0142925) q[0];
sx q[0];
rz(-3.1211299) q[0];
rz(0.19275716) q[1];
sx q[1];
rz(-0.61734504) q[1];
sx q[1];
rz(-1.1654759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11949018) q[0];
sx q[0];
rz(-2.179232) q[0];
sx q[0];
rz(-1.0761989) q[0];
rz(-pi) q[1];
rz(0.67365174) q[2];
sx q[2];
rz(-1.21017) q[2];
sx q[2];
rz(1.1471105) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9679035) q[1];
sx q[1];
rz(-2.6375131) q[1];
sx q[1];
rz(1.6126812) q[1];
rz(1.3806512) q[3];
sx q[3];
rz(-1.4961637) q[3];
sx q[3];
rz(-0.41342218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3141025) q[2];
sx q[2];
rz(-0.9919439) q[2];
sx q[2];
rz(0.52725434) q[2];
rz(-2.5984247) q[3];
sx q[3];
rz(-2.4542464) q[3];
sx q[3];
rz(2.7988722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0813893) q[0];
sx q[0];
rz(-2.9855766) q[0];
sx q[0];
rz(-0.40670893) q[0];
rz(1.1609062) q[1];
sx q[1];
rz(-2.2229767) q[1];
sx q[1];
rz(-1.0103753) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.564931) q[0];
sx q[0];
rz(-2.1157678) q[0];
sx q[0];
rz(3.0002563) q[0];
x q[1];
rz(-2.7356304) q[2];
sx q[2];
rz(-0.34892198) q[2];
sx q[2];
rz(2.5418848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95376172) q[1];
sx q[1];
rz(-1.7825025) q[1];
sx q[1];
rz(-2.7642438) q[1];
x q[2];
rz(-0.80548894) q[3];
sx q[3];
rz(-1.7542766) q[3];
sx q[3];
rz(-3.0940521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0747718) q[2];
sx q[2];
rz(-1.8596884) q[2];
sx q[2];
rz(-1.034896) q[2];
rz(-1.3567989) q[3];
sx q[3];
rz(-0.59485888) q[3];
sx q[3];
rz(2.518173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9676301) q[0];
sx q[0];
rz(-1.5246464) q[0];
sx q[0];
rz(-0.3279283) q[0];
rz(0.11218849) q[1];
sx q[1];
rz(-1.9393549) q[1];
sx q[1];
rz(-0.97253886) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16614322) q[0];
sx q[0];
rz(-1.219141) q[0];
sx q[0];
rz(-2.1696198) q[0];
rz(0.29914029) q[2];
sx q[2];
rz(-1.9950331) q[2];
sx q[2];
rz(-1.5787293) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2518721) q[1];
sx q[1];
rz(-2.1600284) q[1];
sx q[1];
rz(-1.5934029) q[1];
rz(1.2054382) q[3];
sx q[3];
rz(-2.1653173) q[3];
sx q[3];
rz(0.48156092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5358676) q[2];
sx q[2];
rz(-0.87656993) q[2];
sx q[2];
rz(0.65315872) q[2];
rz(2.7086835) q[3];
sx q[3];
rz(-0.50986367) q[3];
sx q[3];
rz(0.21231095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2046278) q[0];
sx q[0];
rz(-2.300394) q[0];
sx q[0];
rz(2.9168108) q[0];
rz(2.3460491) q[1];
sx q[1];
rz(-1.3139775) q[1];
sx q[1];
rz(1.068211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5607308) q[0];
sx q[0];
rz(-0.95968548) q[0];
sx q[0];
rz(-0.9471425) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62309391) q[2];
sx q[2];
rz(-0.87298191) q[2];
sx q[2];
rz(2.0562003) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2676937) q[1];
sx q[1];
rz(-1.5452478) q[1];
sx q[1];
rz(2.4539095) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3146888) q[3];
sx q[3];
rz(-2.912622) q[3];
sx q[3];
rz(-3.1310905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.75605679) q[2];
sx q[2];
rz(-1.7588561) q[2];
sx q[2];
rz(-0.63560152) q[2];
rz(2.8086737) q[3];
sx q[3];
rz(-1.0115441) q[3];
sx q[3];
rz(2.6788768) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3391089) q[0];
sx q[0];
rz(-0.026263069) q[0];
sx q[0];
rz(-0.12301692) q[0];
rz(1.3975551) q[1];
sx q[1];
rz(-1.7536283) q[1];
sx q[1];
rz(2.3618598) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5194979) q[0];
sx q[0];
rz(-2.4197398) q[0];
sx q[0];
rz(1.1546385) q[0];
rz(-pi) q[1];
rz(-1.7177203) q[2];
sx q[2];
rz(-0.89040745) q[2];
sx q[2];
rz(2.9522487) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.349682) q[1];
sx q[1];
rz(-1.9432707) q[1];
sx q[1];
rz(2.9279686) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0340683) q[3];
sx q[3];
rz(-1.4443099) q[3];
sx q[3];
rz(-1.3832987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4745549) q[2];
sx q[2];
rz(-1.8755308) q[2];
sx q[2];
rz(3.0511268) q[2];
rz(0.41680923) q[3];
sx q[3];
rz(-0.79206812) q[3];
sx q[3];
rz(-2.1478103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
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
rz(1.6592634) q[0];
sx q[0];
rz(-1.3220795) q[0];
sx q[0];
rz(-3.0521159) q[0];
rz(0.80884519) q[1];
sx q[1];
rz(-2.6409812) q[1];
sx q[1];
rz(-2.083875) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19729511) q[0];
sx q[0];
rz(-0.57387251) q[0];
sx q[0];
rz(-2.3969335) q[0];
rz(-pi) q[1];
rz(1.5799149) q[2];
sx q[2];
rz(-0.55202019) q[2];
sx q[2];
rz(2.7087351) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8174929) q[1];
sx q[1];
rz(-1.1760532) q[1];
sx q[1];
rz(2.7604719) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5954893) q[3];
sx q[3];
rz(-2.3273483) q[3];
sx q[3];
rz(2.7443742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40314254) q[2];
sx q[2];
rz(-0.53692997) q[2];
sx q[2];
rz(0.37330791) q[2];
rz(-0.25660723) q[3];
sx q[3];
rz(-1.4213057) q[3];
sx q[3];
rz(-3.0671425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858717) q[0];
sx q[0];
rz(-1.5789565) q[0];
sx q[0];
rz(-1.4174905) q[0];
rz(1.7120842) q[1];
sx q[1];
rz(-2.7919339) q[1];
sx q[1];
rz(-2.3333593) q[1];
rz(-2.6827742) q[2];
sx q[2];
rz(-1.5622258) q[2];
sx q[2];
rz(-0.020389204) q[2];
rz(3.1028845) q[3];
sx q[3];
rz(-0.9699655) q[3];
sx q[3];
rz(-1.7326596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
