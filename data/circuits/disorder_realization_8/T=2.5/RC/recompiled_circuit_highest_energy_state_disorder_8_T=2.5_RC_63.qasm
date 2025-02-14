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
rz(1.3266069) q[0];
sx q[0];
rz(-2.9171483) q[0];
sx q[0];
rz(1.8087968) q[0];
rz(-0.83238554) q[1];
sx q[1];
rz(-1.4406942) q[1];
sx q[1];
rz(1.1103777) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.719037) q[0];
sx q[0];
rz(-2.9326322) q[0];
sx q[0];
rz(2.4709769) q[0];
x q[1];
rz(0.30958561) q[2];
sx q[2];
rz(-0.88304115) q[2];
sx q[2];
rz(-1.0716455) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77423461) q[1];
sx q[1];
rz(-1.6239161) q[1];
sx q[1];
rz(1.5989283) q[1];
rz(-pi) q[2];
rz(-0.63817231) q[3];
sx q[3];
rz(-1.8211357) q[3];
sx q[3];
rz(-1.7145715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9139468) q[2];
sx q[2];
rz(-3.1326742) q[2];
sx q[2];
rz(1.439636) q[2];
rz(1.7281744) q[3];
sx q[3];
rz(-0.011912502) q[3];
sx q[3];
rz(-0.2695151) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.444376) q[0];
sx q[0];
rz(-1.603729) q[0];
sx q[0];
rz(0.61475301) q[0];
rz(-2.6055824) q[1];
sx q[1];
rz(-0.025608048) q[1];
sx q[1];
rz(2.8047309) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5081072) q[0];
sx q[0];
rz(-2.9683873) q[0];
sx q[0];
rz(-2.2499535) q[0];
rz(-3.1104149) q[2];
sx q[2];
rz(-2.665069) q[2];
sx q[2];
rz(2.2533803) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.555111) q[1];
sx q[1];
rz(-3.0947436) q[1];
sx q[1];
rz(1.2270801) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2516465) q[3];
sx q[3];
rz(-1.7515514) q[3];
sx q[3];
rz(1.1584566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2577995) q[2];
sx q[2];
rz(-3.1287441) q[2];
sx q[2];
rz(0.69965714) q[2];
rz(-2.9720225) q[3];
sx q[3];
rz(-3.1278059) q[3];
sx q[3];
rz(0.56210303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4366142) q[0];
sx q[0];
rz(-0.545937) q[0];
sx q[0];
rz(-2.8563232) q[0];
rz(2.8778695) q[1];
sx q[1];
rz(-0.00058760651) q[1];
sx q[1];
rz(0.79424167) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6503822) q[0];
sx q[0];
rz(-0.69612078) q[0];
sx q[0];
rz(-2.9589423) q[0];
rz(-pi) q[1];
rz(1.2474485) q[2];
sx q[2];
rz(-1.4686606) q[2];
sx q[2];
rz(2.7030526) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2546187) q[1];
sx q[1];
rz(-1.6052263) q[1];
sx q[1];
rz(-1.5772343) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91762791) q[3];
sx q[3];
rz(-2.1289187) q[3];
sx q[3];
rz(2.2150326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11702015) q[2];
sx q[2];
rz(-0.046155013) q[2];
sx q[2];
rz(-0.864492) q[2];
rz(2.3965059) q[3];
sx q[3];
rz(-2.2741208) q[3];
sx q[3];
rz(0.043070506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87557781) q[0];
sx q[0];
rz(-3.0953396) q[0];
sx q[0];
rz(0.86638802) q[0];
rz(0.59151793) q[1];
sx q[1];
rz(-0.94647995) q[1];
sx q[1];
rz(-1.0513069) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0795938) q[0];
sx q[0];
rz(-1.5709086) q[0];
sx q[0];
rz(-3.1395509) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7655623) q[2];
sx q[2];
rz(-3.1392619) q[2];
sx q[2];
rz(-0.19007209) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4128139) q[1];
sx q[1];
rz(-2.3024913) q[1];
sx q[1];
rz(1.0207291) q[1];
rz(-pi) q[2];
rz(2.5539033) q[3];
sx q[3];
rz(-1.1212548) q[3];
sx q[3];
rz(1.12236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5279348) q[2];
sx q[2];
rz(-3.077007) q[2];
sx q[2];
rz(-2.2876372) q[2];
rz(-2.4438786) q[3];
sx q[3];
rz(-1.8055975) q[3];
sx q[3];
rz(0.31049389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7455604) q[0];
sx q[0];
rz(-0.10265352) q[0];
sx q[0];
rz(0.41203099) q[0];
rz(1.283006) q[1];
sx q[1];
rz(-2.367815) q[1];
sx q[1];
rz(2.3903019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.368741) q[0];
sx q[0];
rz(-1.6882467) q[0];
sx q[0];
rz(2.9020578) q[0];
rz(-2.9137026) q[2];
sx q[2];
rz(-2.2376559) q[2];
sx q[2];
rz(1.6713072) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4163782) q[1];
sx q[1];
rz(-2.0592923) q[1];
sx q[1];
rz(-0.67711551) q[1];
rz(-3.0194332) q[3];
sx q[3];
rz(-1.7137104) q[3];
sx q[3];
rz(2.0412594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46386197) q[2];
sx q[2];
rz(-0.025721392) q[2];
sx q[2];
rz(0.99979293) q[2];
rz(-0.60100466) q[3];
sx q[3];
rz(-3.0809564) q[3];
sx q[3];
rz(0.84990466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.8782225) q[0];
sx q[0];
rz(-0.090494089) q[0];
sx q[0];
rz(2.7409842) q[0];
rz(1.406631) q[1];
sx q[1];
rz(-2.7901283) q[1];
sx q[1];
rz(1.7469143) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18086223) q[0];
sx q[0];
rz(-0.020190857) q[0];
sx q[0];
rz(-1.8024615) q[0];
rz(-pi) q[1];
rz(2.8408634) q[2];
sx q[2];
rz(-0.82737714) q[2];
sx q[2];
rz(-2.9234795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1876389) q[1];
sx q[1];
rz(-1.5752183) q[1];
sx q[1];
rz(-0.19051801) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9506545) q[3];
sx q[3];
rz(-0.74511601) q[3];
sx q[3];
rz(2.7332954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3957735) q[2];
sx q[2];
rz(-0.58818156) q[2];
sx q[2];
rz(2.0664717) q[2];
rz(0.51946688) q[3];
sx q[3];
rz(-2.963701) q[3];
sx q[3];
rz(1.5943257) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2125242) q[0];
sx q[0];
rz(-1.2146177) q[0];
sx q[0];
rz(-1.6125096) q[0];
rz(0.56135881) q[1];
sx q[1];
rz(-3.1293271) q[1];
sx q[1];
rz(2.5711109) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5382696) q[0];
sx q[0];
rz(-3.0099359) q[0];
sx q[0];
rz(1.1696474) q[0];
rz(2.1856603) q[2];
sx q[2];
rz(-0.38262767) q[2];
sx q[2];
rz(-2.0418233) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7509527) q[1];
sx q[1];
rz(-1.5719853) q[1];
sx q[1];
rz(-1.5863064) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6030406) q[3];
sx q[3];
rz(-0.90044124) q[3];
sx q[3];
rz(0.30645257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0824215) q[2];
sx q[2];
rz(-0.26796451) q[2];
sx q[2];
rz(1.9783665) q[2];
rz(1.6022812) q[3];
sx q[3];
rz(-0.10207615) q[3];
sx q[3];
rz(0.99360895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.38692835) q[0];
sx q[0];
rz(-2.8480777) q[0];
sx q[0];
rz(-2.5013404) q[0];
rz(2.2786268) q[1];
sx q[1];
rz(-2.9997365) q[1];
sx q[1];
rz(-2.364025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95381584) q[0];
sx q[0];
rz(-1.9627155) q[0];
sx q[0];
rz(-1.1338393) q[0];
x q[1];
rz(-2.8993494) q[2];
sx q[2];
rz(-2.3628334) q[2];
sx q[2];
rz(2.3295516) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.41987644) q[1];
sx q[1];
rz(-1.6437746) q[1];
sx q[1];
rz(3.130129) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4005017) q[3];
sx q[3];
rz(-0.6395542) q[3];
sx q[3];
rz(-0.049413817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5295279) q[2];
sx q[2];
rz(-2.7175588) q[2];
sx q[2];
rz(2.6112134) q[2];
rz(0.027675962) q[3];
sx q[3];
rz(-0.045462463) q[3];
sx q[3];
rz(1.8224705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46488047) q[0];
sx q[0];
rz(-2.9806529) q[0];
sx q[0];
rz(2.2093534) q[0];
rz(1.543401) q[1];
sx q[1];
rz(-0.80359572) q[1];
sx q[1];
rz(2.7995321) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8810976) q[0];
sx q[0];
rz(-1.6552345) q[0];
sx q[0];
rz(1.5556094) q[0];
rz(-pi) q[1];
rz(0.73782302) q[2];
sx q[2];
rz(-1.617811) q[2];
sx q[2];
rz(1.1931057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4885713) q[1];
sx q[1];
rz(-2.584051) q[1];
sx q[1];
rz(0.36399059) q[1];
x q[2];
rz(-1.341712) q[3];
sx q[3];
rz(-1.5543665) q[3];
sx q[3];
rz(-1.5180902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0590234) q[2];
sx q[2];
rz(-0.00073585357) q[2];
sx q[2];
rz(-1.9334582) q[2];
rz(-0.9515323) q[3];
sx q[3];
rz(-0.0080684302) q[3];
sx q[3];
rz(0.72037303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35219881) q[0];
sx q[0];
rz(-2.3645526) q[0];
sx q[0];
rz(-3.0598031) q[0];
rz(-2.8261322) q[1];
sx q[1];
rz(-3.0888562) q[1];
sx q[1];
rz(-1.2299406) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8982809) q[0];
sx q[0];
rz(-1.7276242) q[0];
sx q[0];
rz(2.9812198) q[0];
x q[1];
rz(0.66218485) q[2];
sx q[2];
rz(-2.9176468) q[2];
sx q[2];
rz(-2.0020773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9435159) q[1];
sx q[1];
rz(-1.6606037) q[1];
sx q[1];
rz(1.4526558) q[1];
x q[2];
rz(-3.0026912) q[3];
sx q[3];
rz(-1.0437696) q[3];
sx q[3];
rz(-2.4580818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6007467) q[2];
sx q[2];
rz(-3.1223065) q[2];
sx q[2];
rz(-2.02796) q[2];
rz(-3.1019548) q[3];
sx q[3];
rz(-3.1318635) q[3];
sx q[3];
rz(-0.59687328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6780846) q[0];
sx q[0];
rz(-1.9263374) q[0];
sx q[0];
rz(1.3081464) q[0];
rz(-0.61617638) q[1];
sx q[1];
rz(-2.0694852) q[1];
sx q[1];
rz(-2.9328666) q[1];
rz(-2.2439416) q[2];
sx q[2];
rz(-2.8406526) q[2];
sx q[2];
rz(0.1546897) q[2];
rz(1.7442654) q[3];
sx q[3];
rz(-2.9174505) q[3];
sx q[3];
rz(2.9944921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
