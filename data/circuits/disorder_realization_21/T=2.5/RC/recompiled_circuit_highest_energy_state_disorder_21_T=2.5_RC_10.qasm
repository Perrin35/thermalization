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
rz(-1.4323956) q[0];
sx q[0];
rz(-0.30269912) q[0];
sx q[0];
rz(-2.1085289) q[0];
rz(-1.2248224) q[1];
sx q[1];
rz(-1.4900102) q[1];
sx q[1];
rz(2.9472247) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8914999) q[0];
sx q[0];
rz(-2.19133) q[0];
sx q[0];
rz(-1.7583855) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35440234) q[2];
sx q[2];
rz(-0.72847073) q[2];
sx q[2];
rz(0.013106339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74339515) q[1];
sx q[1];
rz(-1.6170073) q[1];
sx q[1];
rz(3.1070437) q[1];
x q[2];
rz(-1.1367873) q[3];
sx q[3];
rz(-2.1510486) q[3];
sx q[3];
rz(-2.2541719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.73244652) q[2];
sx q[2];
rz(-2.0398085) q[2];
sx q[2];
rz(0.5644325) q[2];
rz(1.8730646) q[3];
sx q[3];
rz(-0.0081491834) q[3];
sx q[3];
rz(-2.234999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0821575) q[0];
sx q[0];
rz(-3.1340288) q[0];
sx q[0];
rz(-0.91307688) q[0];
rz(-2.4572241) q[1];
sx q[1];
rz(-3.1412536) q[1];
sx q[1];
rz(-2.2025462) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20439273) q[0];
sx q[0];
rz(-2.4449722) q[0];
sx q[0];
rz(-0.60782363) q[0];
x q[1];
rz(-2.6776161) q[2];
sx q[2];
rz(-1.5590406) q[2];
sx q[2];
rz(1.8092312) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1176743) q[1];
sx q[1];
rz(-0.54421762) q[1];
sx q[1];
rz(-0.88440374) q[1];
rz(-1.9775852) q[3];
sx q[3];
rz(-1.6322246) q[3];
sx q[3];
rz(2.8451426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82274503) q[2];
sx q[2];
rz(-3.1343967) q[2];
sx q[2];
rz(2.2491271) q[2];
rz(-1.5626296) q[3];
sx q[3];
rz(-0.024450863) q[3];
sx q[3];
rz(-1.8306556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398294) q[0];
sx q[0];
rz(-2.6391397) q[0];
sx q[0];
rz(2.7318562) q[0];
rz(3.1306664) q[1];
sx q[1];
rz(-0.22053638) q[1];
sx q[1];
rz(-1.3440557) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4554268) q[0];
sx q[0];
rz(-0.51138216) q[0];
sx q[0];
rz(-1.6126942) q[0];
rz(-1.0667886) q[2];
sx q[2];
rz(-3.0098999) q[2];
sx q[2];
rz(-1.7273336) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5907659) q[1];
sx q[1];
rz(-1.3129586) q[1];
sx q[1];
rz(-3.0745561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6556754) q[3];
sx q[3];
rz(-0.1156919) q[3];
sx q[3];
rz(1.8574238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4695796) q[2];
sx q[2];
rz(-1.4521705) q[2];
sx q[2];
rz(3.0984042) q[2];
rz(-2.0016661) q[3];
sx q[3];
rz(-0.15634263) q[3];
sx q[3];
rz(-1.9388916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7286872) q[0];
sx q[0];
rz(-0.40189704) q[0];
sx q[0];
rz(-1.3380949) q[0];
rz(-2.4463704) q[1];
sx q[1];
rz(-0.1278563) q[1];
sx q[1];
rz(-2.7241838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25487568) q[0];
sx q[0];
rz(-0.69807295) q[0];
sx q[0];
rz(1.3017734) q[0];
rz(-pi) q[1];
rz(-0.92441332) q[2];
sx q[2];
rz(-1.7862537) q[2];
sx q[2];
rz(0.59981031) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.771811) q[1];
sx q[1];
rz(-2.2038748) q[1];
sx q[1];
rz(-2.8903583) q[1];
x q[2];
rz(-1.2715075) q[3];
sx q[3];
rz(-2.6911754) q[3];
sx q[3];
rz(-1.3189486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7030316) q[2];
sx q[2];
rz(-0.87623864) q[2];
sx q[2];
rz(1.7523127) q[2];
rz(-0.82052463) q[3];
sx q[3];
rz(-1.5503784) q[3];
sx q[3];
rz(2.4337721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.163212) q[0];
sx q[0];
rz(-3.1040525) q[0];
sx q[0];
rz(0.96330825) q[0];
rz(-2.9523201) q[1];
sx q[1];
rz(-3.1258588) q[1];
sx q[1];
rz(-2.9761369) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2236639) q[0];
sx q[0];
rz(-0.05089137) q[0];
sx q[0];
rz(2.5975157) q[0];
rz(2.2654183) q[2];
sx q[2];
rz(-1.8037705) q[2];
sx q[2];
rz(2.7829952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7041118) q[1];
sx q[1];
rz(-0.81161122) q[1];
sx q[1];
rz(0.1603756) q[1];
x q[2];
rz(2.1894375) q[3];
sx q[3];
rz(-2.6252271) q[3];
sx q[3];
rz(-2.0455893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8339707) q[2];
sx q[2];
rz(-2.622719) q[2];
sx q[2];
rz(-1.0597672) q[2];
rz(-2.4506532) q[3];
sx q[3];
rz(-2.2773404) q[3];
sx q[3];
rz(1.8701657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6178013) q[0];
sx q[0];
rz(-3.0484564) q[0];
sx q[0];
rz(1.5366489) q[0];
rz(-0.45283428) q[1];
sx q[1];
rz(-3.1335399) q[1];
sx q[1];
rz(1.4401999) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33537597) q[0];
sx q[0];
rz(-1.7378238) q[0];
sx q[0];
rz(1.7352292) q[0];
rz(-pi) q[1];
rz(-0.22162205) q[2];
sx q[2];
rz(-1.3563507) q[2];
sx q[2];
rz(-0.947244) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.837484) q[1];
sx q[1];
rz(-2.9867801) q[1];
sx q[1];
rz(-0.047708851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0945419) q[3];
sx q[3];
rz(-0.19652995) q[3];
sx q[3];
rz(2.9158338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3692533) q[2];
sx q[2];
rz(-0.99268308) q[2];
sx q[2];
rz(3.0625647) q[2];
rz(2.2760271) q[3];
sx q[3];
rz(-0.95439684) q[3];
sx q[3];
rz(1.0237833) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9158151) q[0];
sx q[0];
rz(-0.0053891698) q[0];
sx q[0];
rz(-2.916577) q[0];
rz(0.30613884) q[1];
sx q[1];
rz(-3.1248416) q[1];
sx q[1];
rz(0.87047815) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7891312) q[0];
sx q[0];
rz(-1.6183743) q[0];
sx q[0];
rz(-0.0713047) q[0];
x q[1];
rz(-0.71484675) q[2];
sx q[2];
rz(-1.1230334) q[2];
sx q[2];
rz(-0.078850672) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0902033) q[1];
sx q[1];
rz(-1.6312113) q[1];
sx q[1];
rz(-2.3714972) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2629396) q[3];
sx q[3];
rz(-2.5402398) q[3];
sx q[3];
rz(2.2215312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47591448) q[2];
sx q[2];
rz(-3.0147538) q[2];
sx q[2];
rz(-2.1062984) q[2];
rz(0.78694844) q[3];
sx q[3];
rz(-3.0988155) q[3];
sx q[3];
rz(-0.27534819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03054522) q[0];
sx q[0];
rz(-0.011140911) q[0];
sx q[0];
rz(-0.025644843) q[0];
rz(1.8772839) q[1];
sx q[1];
rz(-0.023921078) q[1];
sx q[1];
rz(-0.69902507) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049533904) q[0];
sx q[0];
rz(-1.7932685) q[0];
sx q[0];
rz(0.27357863) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8139171) q[2];
sx q[2];
rz(-0.49105308) q[2];
sx q[2];
rz(-1.9126522) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5311218) q[1];
sx q[1];
rz(-1.4969331) q[1];
sx q[1];
rz(-2.7570037) q[1];
rz(-0.54663901) q[3];
sx q[3];
rz(-2.1501503) q[3];
sx q[3];
rz(-2.1648615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.25253025) q[2];
sx q[2];
rz(-0.74873304) q[2];
sx q[2];
rz(0.32386455) q[2];
rz(2.9845386) q[3];
sx q[3];
rz(-1.9040949) q[3];
sx q[3];
rz(-2.9878555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57482982) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(-2.5515442) q[0];
rz(-2.4216962) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(-0.86729008) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9218719) q[0];
sx q[0];
rz(-0.94014535) q[0];
sx q[0];
rz(-1.7008855) q[0];
x q[1];
rz(1.2008823) q[2];
sx q[2];
rz(-2.5061766) q[2];
sx q[2];
rz(2.4398092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5979851) q[1];
sx q[1];
rz(-1.6389009) q[1];
sx q[1];
rz(1.4582514) q[1];
x q[2];
rz(1.7884364) q[3];
sx q[3];
rz(-3*pi/13) q[3];
sx q[3];
rz(2.9293037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94816339) q[2];
sx q[2];
rz(-0.89274222) q[2];
sx q[2];
rz(0.11290045) q[2];
rz(-1.0490949) q[3];
sx q[3];
rz(-1.5821404) q[3];
sx q[3];
rz(0.73672867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0650487) q[0];
sx q[0];
rz(-1.5815409) q[0];
sx q[0];
rz(-1.5034058) q[0];
rz(0.63649559) q[1];
sx q[1];
rz(-2.681585) q[1];
sx q[1];
rz(-1.5696625) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51278567) q[0];
sx q[0];
rz(-0.11462002) q[0];
sx q[0];
rz(-0.8091457) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9365065) q[2];
sx q[2];
rz(-3.1383927) q[2];
sx q[2];
rz(-2.5815258) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1366982) q[1];
sx q[1];
rz(-1.5678798) q[1];
sx q[1];
rz(1.5705622) q[1];
rz(0.51497634) q[3];
sx q[3];
rz(-1.4287717) q[3];
sx q[3];
rz(1.2339301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.87127176) q[2];
sx q[2];
rz(-1.511829) q[2];
sx q[2];
rz(-2.9809269) q[2];
rz(-1.2284944) q[3];
sx q[3];
rz(-0.03385032) q[3];
sx q[3];
rz(-0.20139774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092030839) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(1.6231712) q[1];
sx q[1];
rz(-0.36133125) q[1];
sx q[1];
rz(0.28892118) q[1];
rz(-1.8636017) q[2];
sx q[2];
rz(-1.5488727) q[2];
sx q[2];
rz(1.9643754) q[2];
rz(-1.2613983) q[3];
sx q[3];
rz(-0.45844504) q[3];
sx q[3];
rz(-2.2506056) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
