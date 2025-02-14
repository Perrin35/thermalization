OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91035834) q[0];
sx q[0];
rz(-2.2725821) q[0];
sx q[0];
rz(2.056871) q[0];
rz(1.9864858) q[1];
sx q[1];
rz(-2.3218563) q[1];
sx q[1];
rz(-2.2302332) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93854967) q[0];
sx q[0];
rz(-1.8126376) q[0];
sx q[0];
rz(-0.77745243) q[0];
rz(1.4248104) q[2];
sx q[2];
rz(-2.2772191) q[2];
sx q[2];
rz(-1.8391158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6877796) q[1];
sx q[1];
rz(-2.0168843) q[1];
sx q[1];
rz(-0.58961726) q[1];
rz(-pi) q[2];
rz(1.7098268) q[3];
sx q[3];
rz(-2.7043501) q[3];
sx q[3];
rz(-2.4831877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.800941) q[2];
sx q[2];
rz(-2.030535) q[2];
sx q[2];
rz(3.0621373) q[2];
rz(0.60845145) q[3];
sx q[3];
rz(-1.918101) q[3];
sx q[3];
rz(-0.89103812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-0.65421739) q[0];
sx q[0];
rz(-0.86367622) q[0];
sx q[0];
rz(1.7864216) q[0];
rz(2.1036509) q[1];
sx q[1];
rz(-2.4619921) q[1];
sx q[1];
rz(-0.80744809) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8184745) q[0];
sx q[0];
rz(-1.8606631) q[0];
sx q[0];
rz(-2.1978343) q[0];
rz(-pi) q[1];
rz(2.1218929) q[2];
sx q[2];
rz(-2.4906213) q[2];
sx q[2];
rz(0.62007444) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5894306) q[1];
sx q[1];
rz(-2.8928601) q[1];
sx q[1];
rz(2.8692607) q[1];
rz(-0.89887606) q[3];
sx q[3];
rz(-2.9714317) q[3];
sx q[3];
rz(-0.32246043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9281533) q[2];
sx q[2];
rz(-1.5779053) q[2];
sx q[2];
rz(-1.4595002) q[2];
rz(-2.6835119) q[3];
sx q[3];
rz(-0.57772485) q[3];
sx q[3];
rz(-2.1817082) q[3];
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
rz(-2.9816575) q[0];
sx q[0];
rz(-2.0913048) q[0];
sx q[0];
rz(-2.4666069) q[0];
rz(-2.5534897) q[1];
sx q[1];
rz(-2.2876078) q[1];
sx q[1];
rz(1.810422) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0818401) q[0];
sx q[0];
rz(-0.65048238) q[0];
sx q[0];
rz(-0.68997927) q[0];
x q[1];
rz(-0.66548062) q[2];
sx q[2];
rz(-1.988171) q[2];
sx q[2];
rz(-0.62705428) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.59440295) q[1];
sx q[1];
rz(-1.8246027) q[1];
sx q[1];
rz(0.82271346) q[1];
rz(1.9310476) q[3];
sx q[3];
rz(-0.52740806) q[3];
sx q[3];
rz(1.0077734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62446928) q[2];
sx q[2];
rz(-0.34999592) q[2];
sx q[2];
rz(2.9208753) q[2];
rz(0.80495009) q[3];
sx q[3];
rz(-1.5996108) q[3];
sx q[3];
rz(-1.8055003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.540156) q[0];
sx q[0];
rz(-2.8995081) q[0];
sx q[0];
rz(0.618774) q[0];
rz(-3.0586808) q[1];
sx q[1];
rz(-0.34268788) q[1];
sx q[1];
rz(1.6212911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5280209) q[0];
sx q[0];
rz(-2.6074886) q[0];
sx q[0];
rz(-1.7604802) q[0];
rz(1.8074715) q[2];
sx q[2];
rz(-2.1370721) q[2];
sx q[2];
rz(-0.90619722) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.027923) q[1];
sx q[1];
rz(-1.4884559) q[1];
sx q[1];
rz(-2.0827977) q[1];
rz(-1.3906562) q[3];
sx q[3];
rz(-1.9400175) q[3];
sx q[3];
rz(0.1843017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3380022) q[2];
sx q[2];
rz(-3.0323196) q[2];
sx q[2];
rz(-2.9962311) q[2];
rz(-0.27211443) q[3];
sx q[3];
rz(-1.4067255) q[3];
sx q[3];
rz(-1.9947778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014932545) q[0];
sx q[0];
rz(-1.3260051) q[0];
sx q[0];
rz(0.10619157) q[0];
rz(-0.95942489) q[1];
sx q[1];
rz(-0.54490772) q[1];
sx q[1];
rz(0.19439654) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0036164408) q[0];
sx q[0];
rz(-0.73775333) q[0];
sx q[0];
rz(-3.0135148) q[0];
rz(-pi) q[1];
rz(3.1367231) q[2];
sx q[2];
rz(-1.4559573) q[2];
sx q[2];
rz(-1.8094935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0982617) q[1];
sx q[1];
rz(-0.98214591) q[1];
sx q[1];
rz(1.5215988) q[1];
rz(-1.164489) q[3];
sx q[3];
rz(-1.4999564) q[3];
sx q[3];
rz(2.336647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.306119) q[2];
sx q[2];
rz(-1.6872311) q[2];
sx q[2];
rz(-0.56841889) q[2];
rz(0.052637188) q[3];
sx q[3];
rz(-1.5024374) q[3];
sx q[3];
rz(0.13681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0282054) q[0];
sx q[0];
rz(-2.5882692) q[0];
sx q[0];
rz(0.18948874) q[0];
rz(2.1151309) q[1];
sx q[1];
rz(-0.84954134) q[1];
sx q[1];
rz(-2.386327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8403977) q[0];
sx q[0];
rz(-0.9009255) q[0];
sx q[0];
rz(-2.2706881) q[0];
rz(2.3113475) q[2];
sx q[2];
rz(-0.70607042) q[2];
sx q[2];
rz(2.490807) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2456296) q[1];
sx q[1];
rz(-1.4019483) q[1];
sx q[1];
rz(-0.16014512) q[1];
rz(0.35474687) q[3];
sx q[3];
rz(-1.2257396) q[3];
sx q[3];
rz(1.3416895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.75817627) q[2];
sx q[2];
rz(-1.5551609) q[2];
sx q[2];
rz(1.4413393) q[2];
rz(1.7656743) q[3];
sx q[3];
rz(-2.223189) q[3];
sx q[3];
rz(-0.70022303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7021084) q[0];
sx q[0];
rz(-1.0478042) q[0];
sx q[0];
rz(-1.1442319) q[0];
rz(-2.8817835) q[1];
sx q[1];
rz(-0.83994284) q[1];
sx q[1];
rz(2.1536486) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0596302) q[0];
sx q[0];
rz(-3.0407627) q[0];
sx q[0];
rz(-2.9743845) q[0];
rz(-pi) q[1];
rz(1.1455215) q[2];
sx q[2];
rz(-2.2605611) q[2];
sx q[2];
rz(2.9381616) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14079796) q[1];
sx q[1];
rz(-2.4355445) q[1];
sx q[1];
rz(0.39800675) q[1];
rz(0.3026721) q[3];
sx q[3];
rz(-2.1244123) q[3];
sx q[3];
rz(0.68676567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7320431) q[2];
sx q[2];
rz(-1.6403551) q[2];
sx q[2];
rz(0.6130971) q[2];
rz(2.4260855) q[3];
sx q[3];
rz(-2.3698273) q[3];
sx q[3];
rz(2.7806921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9382984) q[0];
sx q[0];
rz(-2.2977915) q[0];
sx q[0];
rz(-2.8734558) q[0];
rz(-2.9109491) q[1];
sx q[1];
rz(-1.5985039) q[1];
sx q[1];
rz(0.014009744) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5174311) q[0];
sx q[0];
rz(-0.78963477) q[0];
sx q[0];
rz(2.8454418) q[0];
x q[1];
rz(-1.2410937) q[2];
sx q[2];
rz(-2.6791141) q[2];
sx q[2];
rz(1.3792737) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18302984) q[1];
sx q[1];
rz(-1.5261478) q[1];
sx q[1];
rz(0.53739287) q[1];
rz(-1.713578) q[3];
sx q[3];
rz(-1.5413949) q[3];
sx q[3];
rz(-2.2063696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1676499) q[2];
sx q[2];
rz(-1.8029982) q[2];
sx q[2];
rz(2.8323284) q[2];
rz(-0.15657982) q[3];
sx q[3];
rz(-1.4045818) q[3];
sx q[3];
rz(1.0369161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2579047) q[0];
sx q[0];
rz(-2.4990999) q[0];
sx q[0];
rz(2.8908492) q[0];
rz(-2.5550487) q[1];
sx q[1];
rz(-1.5637014) q[1];
sx q[1];
rz(-0.7647382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4630554) q[0];
sx q[0];
rz(-1.8856115) q[0];
sx q[0];
rz(-0.30588715) q[0];
rz(-pi) q[1];
x q[1];
rz(2.938835) q[2];
sx q[2];
rz(-2.360095) q[2];
sx q[2];
rz(2.4818713) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8426399) q[1];
sx q[1];
rz(-3.0184426) q[1];
sx q[1];
rz(3.0663436) q[1];
rz(3.0050817) q[3];
sx q[3];
rz(-0.38204604) q[3];
sx q[3];
rz(2.7979224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5137198) q[2];
sx q[2];
rz(-0.79078117) q[2];
sx q[2];
rz(-2.9676843) q[2];
rz(0.74226132) q[3];
sx q[3];
rz(-2.3895013) q[3];
sx q[3];
rz(-0.043005634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.1327508) q[0];
sx q[0];
rz(-0.77021563) q[0];
sx q[0];
rz(2.8736864) q[0];
rz(1.335089) q[1];
sx q[1];
rz(-2.2309512) q[1];
sx q[1];
rz(1.7693899) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5822179) q[0];
sx q[0];
rz(-1.1211044) q[0];
sx q[0];
rz(0.26828464) q[0];
rz(-0.29628654) q[2];
sx q[2];
rz(-1.4742659) q[2];
sx q[2];
rz(-1.5636843) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42555537) q[1];
sx q[1];
rz(-2.4074984) q[1];
sx q[1];
rz(-2.3922763) q[1];
x q[2];
rz(0.80268152) q[3];
sx q[3];
rz(-2.6978328) q[3];
sx q[3];
rz(0.26324031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.043896) q[2];
sx q[2];
rz(-2.1052269) q[2];
sx q[2];
rz(1.3807266) q[2];
rz(-1.1287639) q[3];
sx q[3];
rz(-0.94996101) q[3];
sx q[3];
rz(1.5855764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783405) q[0];
sx q[0];
rz(-1.6077519) q[0];
sx q[0];
rz(2.0335249) q[0];
rz(-0.68645984) q[1];
sx q[1];
rz(-1.3931128) q[1];
sx q[1];
rz(-1.211094) q[1];
rz(-2.6348719) q[2];
sx q[2];
rz(-0.55247775) q[2];
sx q[2];
rz(-1.1341118) q[2];
rz(-0.41175089) q[3];
sx q[3];
rz(-1.0706378) q[3];
sx q[3];
rz(2.0720458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
