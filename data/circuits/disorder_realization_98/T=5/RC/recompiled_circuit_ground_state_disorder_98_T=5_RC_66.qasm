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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41357562) q[0];
sx q[0];
rz(-2.8110162) q[0];
sx q[0];
rz(-2.3982993) q[0];
rz(-pi) q[1];
rz(-2.0091363) q[2];
sx q[2];
rz(-1.1852263) q[2];
sx q[2];
rz(-1.350268) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39241782) q[1];
sx q[1];
rz(-2.2274744) q[1];
sx q[1];
rz(-1.6139469) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6781647) q[3];
sx q[3];
rz(-0.26623785) q[3];
sx q[3];
rz(0.050198089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3898042) q[2];
sx q[2];
rz(-1.3499539) q[2];
sx q[2];
rz(-2.8765615) q[2];
rz(2.7122279) q[3];
sx q[3];
rz(-1.8931484) q[3];
sx q[3];
rz(-3.140669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65381831) q[0];
sx q[0];
rz(-0.14350292) q[0];
sx q[0];
rz(1.300746) q[0];
rz(3.0744413) q[1];
sx q[1];
rz(-0.90922272) q[1];
sx q[1];
rz(-0.86004177) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81102449) q[0];
sx q[0];
rz(-1.0029135) q[0];
sx q[0];
rz(0.83356838) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76171909) q[2];
sx q[2];
rz(-2.2713619) q[2];
sx q[2];
rz(-0.127244) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7262267) q[1];
sx q[1];
rz(-0.72390717) q[1];
sx q[1];
rz(0.93127802) q[1];
rz(-2.9577291) q[3];
sx q[3];
rz(-1.2753092) q[3];
sx q[3];
rz(1.1671403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5106875) q[2];
sx q[2];
rz(-0.61669934) q[2];
sx q[2];
rz(-0.56646937) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94839621) q[0];
sx q[0];
rz(-1.0862792) q[0];
sx q[0];
rz(2.7753944) q[0];
rz(1.5858448) q[1];
sx q[1];
rz(-2.0987174) q[1];
sx q[1];
rz(-0.050447024) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18264601) q[0];
sx q[0];
rz(-0.0894657) q[0];
sx q[0];
rz(1.0854118) q[0];
x q[1];
rz(1.0773727) q[2];
sx q[2];
rz(-0.9443379) q[2];
sx q[2];
rz(-0.74365091) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5844684) q[1];
sx q[1];
rz(-1.4038175) q[1];
sx q[1];
rz(0.22241576) q[1];
rz(-pi) q[2];
rz(-1.1459703) q[3];
sx q[3];
rz(-1.4381988) q[3];
sx q[3];
rz(-1.4736922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.067165) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(-1.3107497) q[2];
rz(1.358076) q[3];
sx q[3];
rz(-0.57099968) q[3];
sx q[3];
rz(1.3091492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27291372) q[0];
sx q[0];
rz(-2.9153115) q[0];
sx q[0];
rz(1.7568461) q[0];
rz(1.2387431) q[1];
sx q[1];
rz(-1.2307931) q[1];
sx q[1];
rz(-1.967427) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096701972) q[0];
sx q[0];
rz(-2.6706123) q[0];
sx q[0];
rz(1.1440008) q[0];
x q[1];
rz(0.27710813) q[2];
sx q[2];
rz(-1.8113675) q[2];
sx q[2];
rz(-1.5391369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8983072) q[1];
sx q[1];
rz(-1.8970058) q[1];
sx q[1];
rz(-2.1065358) q[1];
rz(-pi) q[2];
rz(-0.87808319) q[3];
sx q[3];
rz(-1.2125848) q[3];
sx q[3];
rz(-0.75853759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7063286) q[2];
sx q[2];
rz(-1.8212943) q[2];
sx q[2];
rz(-0.038979385) q[2];
rz(-0.6428166) q[3];
sx q[3];
rz(-2.1233852) q[3];
sx q[3];
rz(-2.5207998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83127999) q[0];
sx q[0];
rz(-2.1273002) q[0];
sx q[0];
rz(-0.020462791) q[0];
rz(0.19275716) q[1];
sx q[1];
rz(-0.61734504) q[1];
sx q[1];
rz(-1.1654759) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63686164) q[0];
sx q[0];
rz(-0.76380542) q[0];
sx q[0];
rz(-0.59815191) q[0];
rz(-2.5978659) q[2];
sx q[2];
rz(-2.3909878) q[2];
sx q[2];
rz(-3.1341396) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.36042696) q[1];
sx q[1];
rz(-1.5505704) q[1];
sx q[1];
rz(-1.0670877) q[1];
x q[2];
rz(1.9475157) q[3];
sx q[3];
rz(-2.93749) q[3];
sx q[3];
rz(-2.3538176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8274902) q[2];
sx q[2];
rz(-2.1496488) q[2];
sx q[2];
rz(-0.52725434) q[2];
rz(-2.5984247) q[3];
sx q[3];
rz(-0.68734622) q[3];
sx q[3];
rz(-2.7988722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0813893) q[0];
sx q[0];
rz(-2.9855766) q[0];
sx q[0];
rz(0.40670893) q[0];
rz(1.1609062) q[1];
sx q[1];
rz(-0.91861594) q[1];
sx q[1];
rz(-2.1312174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079491678) q[0];
sx q[0];
rz(-1.4500424) q[0];
sx q[0];
rz(1.0213721) q[0];
x q[1];
rz(1.7134885) q[2];
sx q[2];
rz(-1.2512904) q[2];
sx q[2];
rz(-2.112893) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1045253) q[1];
sx q[1];
rz(-2.7114093) q[1];
sx q[1];
rz(-2.6135315) q[1];
rz(2.8897417) q[3];
sx q[3];
rz(-0.82149071) q[3];
sx q[3];
rz(1.3497533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0668209) q[2];
sx q[2];
rz(-1.8596884) q[2];
sx q[2];
rz(1.034896) q[2];
rz(1.3567989) q[3];
sx q[3];
rz(-2.5467338) q[3];
sx q[3];
rz(-0.6234197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1739625) q[0];
sx q[0];
rz(-1.6169463) q[0];
sx q[0];
rz(2.8136643) q[0];
rz(0.11218849) q[1];
sx q[1];
rz(-1.9393549) q[1];
sx q[1];
rz(-0.97253886) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5060726) q[0];
sx q[0];
rz(-2.1284261) q[0];
sx q[0];
rz(2.7235759) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0123291) q[2];
sx q[2];
rz(-1.8427094) q[2];
sx q[2];
rz(-0.11833469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2518721) q[1];
sx q[1];
rz(-0.98156428) q[1];
sx q[1];
rz(1.5934029) q[1];
x q[2];
rz(-0.48613207) q[3];
sx q[3];
rz(-2.4555169) q[3];
sx q[3];
rz(0.11763517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5358676) q[2];
sx q[2];
rz(-2.2650227) q[2];
sx q[2];
rz(-2.4884339) q[2];
rz(-2.7086835) q[3];
sx q[3];
rz(-2.631729) q[3];
sx q[3];
rz(-2.9292817) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9369649) q[0];
sx q[0];
rz(-2.300394) q[0];
sx q[0];
rz(0.2247819) q[0];
rz(0.79554355) q[1];
sx q[1];
rz(-1.3139775) q[1];
sx q[1];
rz(2.0733817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5607308) q[0];
sx q[0];
rz(-0.95968548) q[0];
sx q[0];
rz(2.1944502) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3722367) q[2];
sx q[2];
rz(-1.1072259) q[2];
sx q[2];
rz(3.0887512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33419106) q[1];
sx q[1];
rz(-0.68808031) q[1];
sx q[1];
rz(-0.040236878) q[1];
rz(-1.3146888) q[3];
sx q[3];
rz(-0.22897069) q[3];
sx q[3];
rz(3.1310905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3855359) q[2];
sx q[2];
rz(-1.3827366) q[2];
sx q[2];
rz(2.5059911) q[2];
rz(-2.8086737) q[3];
sx q[3];
rz(-1.0115441) q[3];
sx q[3];
rz(-2.6788768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8024837) q[0];
sx q[0];
rz(-0.026263069) q[0];
sx q[0];
rz(-3.0185757) q[0];
rz(-1.7440375) q[1];
sx q[1];
rz(-1.3879644) q[1];
sx q[1];
rz(0.77973286) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132039) q[0];
sx q[0];
rz(-1.3003948) q[0];
sx q[0];
rz(2.2487075) q[0];
x q[1];
rz(0.68569195) q[2];
sx q[2];
rz(-1.6848411) q[2];
sx q[2];
rz(1.6673078) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9537264) q[1];
sx q[1];
rz(-0.426891) q[1];
sx q[1];
rz(-1.0736476) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86984313) q[3];
sx q[3];
rz(-0.16582684) q[3];
sx q[3];
rz(-1.0505249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4745549) q[2];
sx q[2];
rz(-1.8755308) q[2];
sx q[2];
rz(3.0511268) q[2];
rz(-0.41680923) q[3];
sx q[3];
rz(-0.79206812) q[3];
sx q[3];
rz(2.1478103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4823293) q[0];
sx q[0];
rz(-1.8195131) q[0];
sx q[0];
rz(0.089476712) q[0];
rz(0.80884519) q[1];
sx q[1];
rz(-0.50061148) q[1];
sx q[1];
rz(2.083875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0322005) q[0];
sx q[0];
rz(-1.9475749) q[0];
sx q[0];
rz(-2.6978639) q[0];
x q[1];
rz(3.1359768) q[2];
sx q[2];
rz(-2.1227909) q[2];
sx q[2];
rz(0.44356669) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.32409975) q[1];
sx q[1];
rz(-1.9655394) q[1];
sx q[1];
rz(2.7604719) q[1];
rz(-pi) q[2];
rz(3.1154409) q[3];
sx q[3];
rz(-0.75687486) q[3];
sx q[3];
rz(-0.43318403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7384501) q[2];
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
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2557209) q[0];
sx q[0];
rz(-1.5626361) q[0];
sx q[0];
rz(1.7241021) q[0];
rz(1.7120842) q[1];
sx q[1];
rz(-2.7919339) q[1];
sx q[1];
rz(-2.3333593) q[1];
rz(-1.5612372) q[2];
sx q[2];
rz(-2.0295967) q[2];
sx q[2];
rz(1.5461736) q[2];
rz(-0.038708121) q[3];
sx q[3];
rz(-0.9699655) q[3];
sx q[3];
rz(-1.7326596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
