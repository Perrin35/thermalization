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
rz(-0.62491971) q[0];
sx q[0];
rz(-1.3345557) q[0];
sx q[0];
rz(-0.053243756) q[0];
rz(0.48625311) q[1];
sx q[1];
rz(6.1558131) q[1];
sx q[1];
rz(11.131412) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7409748) q[0];
sx q[0];
rz(-1.2865928) q[0];
sx q[0];
rz(-1.7114729) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8743319) q[2];
sx q[2];
rz(-1.4556985) q[2];
sx q[2];
rz(0.54265825) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8813144) q[1];
sx q[1];
rz(-1.5173755) q[1];
sx q[1];
rz(2.747606) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49519914) q[3];
sx q[3];
rz(-0.99118587) q[3];
sx q[3];
rz(-0.51900348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80483738) q[2];
sx q[2];
rz(-1.3433604) q[2];
sx q[2];
rz(2.8678144) q[2];
rz(-1.2182419) q[3];
sx q[3];
rz(-0.67353407) q[3];
sx q[3];
rz(-2.8555253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1186721) q[0];
sx q[0];
rz(-2.3781222) q[0];
sx q[0];
rz(2.3619695) q[0];
rz(0.66462213) q[1];
sx q[1];
rz(-1.9326991) q[1];
sx q[1];
rz(-1.8962616) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1561403) q[0];
sx q[0];
rz(-0.13053556) q[0];
sx q[0];
rz(-2.0640316) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0491291) q[2];
sx q[2];
rz(-0.86148724) q[2];
sx q[2];
rz(0.48395983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1511953) q[1];
sx q[1];
rz(-1.510069) q[1];
sx q[1];
rz(1.3582188) q[1];
x q[2];
rz(-2.1297203) q[3];
sx q[3];
rz(-1.3959612) q[3];
sx q[3];
rz(1.8421941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43210426) q[2];
sx q[2];
rz(-0.69807845) q[2];
sx q[2];
rz(3.091605) q[2];
rz(-1.6677808) q[3];
sx q[3];
rz(-1.4155017) q[3];
sx q[3];
rz(0.93562359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.98899984) q[0];
sx q[0];
rz(-0.86243668) q[0];
sx q[0];
rz(-1.1784026) q[0];
rz(1.9940469) q[1];
sx q[1];
rz(-0.18745628) q[1];
sx q[1];
rz(-0.60290927) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99166223) q[0];
sx q[0];
rz(-1.4105182) q[0];
sx q[0];
rz(-3.0383238) q[0];
x q[1];
rz(-0.88073894) q[2];
sx q[2];
rz(-1.9866441) q[2];
sx q[2];
rz(-1.8772454) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47040874) q[1];
sx q[1];
rz(-0.8200596) q[1];
sx q[1];
rz(0.57826184) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.892671) q[3];
sx q[3];
rz(-1.1436) q[3];
sx q[3];
rz(0.77208608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19043645) q[2];
sx q[2];
rz(-0.78593212) q[2];
sx q[2];
rz(2.7583165) q[2];
rz(-1.8303309) q[3];
sx q[3];
rz(-2.0537328) q[3];
sx q[3];
rz(0.12152984) q[3];
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
rz(2.424778) q[0];
sx q[0];
rz(-0.99262339) q[0];
sx q[0];
rz(-0.92856652) q[0];
rz(2.3021452) q[1];
sx q[1];
rz(-1.3176368) q[1];
sx q[1];
rz(0.80926698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033484785) q[0];
sx q[0];
rz(-1.2039469) q[0];
sx q[0];
rz(2.1686694) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7116466) q[2];
sx q[2];
rz(-1.4166178) q[2];
sx q[2];
rz(-1.7183314) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2553808) q[1];
sx q[1];
rz(-1.6396697) q[1];
sx q[1];
rz(-1.3248454) q[1];
x q[2];
rz(2.094481) q[3];
sx q[3];
rz(-0.71412702) q[3];
sx q[3];
rz(1.7585825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3889918) q[2];
sx q[2];
rz(-1.5316803) q[2];
sx q[2];
rz(1.4468225) q[2];
rz(2.0145448) q[3];
sx q[3];
rz(-2.4067252) q[3];
sx q[3];
rz(-2.5126422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22555722) q[0];
sx q[0];
rz(-1.1981523) q[0];
sx q[0];
rz(1.7115364) q[0];
rz(2.9043708) q[1];
sx q[1];
rz(-0.96313852) q[1];
sx q[1];
rz(-0.29475862) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65564102) q[0];
sx q[0];
rz(-1.7136586) q[0];
sx q[0];
rz(2.7014707) q[0];
rz(-2.0343658) q[2];
sx q[2];
rz(-2.4245302) q[2];
sx q[2];
rz(-2.0328558) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7511616) q[1];
sx q[1];
rz(-1.1220349) q[1];
sx q[1];
rz(0.46350355) q[1];
rz(-pi) q[2];
x q[2];
rz(1.421991) q[3];
sx q[3];
rz(-2.1396881) q[3];
sx q[3];
rz(-2.1619145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99990591) q[2];
sx q[2];
rz(-2.9408231) q[2];
sx q[2];
rz(3.0086009) q[2];
rz(-2.3717132) q[3];
sx q[3];
rz(-1.8498288) q[3];
sx q[3];
rz(-2.8225115) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28984508) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(-1.4122562) q[0];
rz(-3.0899561) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(2.9488865) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.76336) q[0];
sx q[0];
rz(-2.6274649) q[0];
sx q[0];
rz(-0.038551081) q[0];
rz(0.95780722) q[2];
sx q[2];
rz(-2.2039206) q[2];
sx q[2];
rz(-1.9547966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.440408) q[1];
sx q[1];
rz(-1.9936947) q[1];
sx q[1];
rz(-0.50114034) q[1];
rz(-pi) q[2];
rz(0.49183853) q[3];
sx q[3];
rz(-0.70295111) q[3];
sx q[3];
rz(0.79849488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0685588) q[2];
sx q[2];
rz(-1.8517588) q[2];
sx q[2];
rz(-1.7603091) q[2];
rz(2.6265788) q[3];
sx q[3];
rz(-0.88740715) q[3];
sx q[3];
rz(-1.7611586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16159049) q[0];
sx q[0];
rz(-1.9965633) q[0];
sx q[0];
rz(0.34647754) q[0];
rz(1.0193635) q[1];
sx q[1];
rz(-2.4788224) q[1];
sx q[1];
rz(-1.4541218) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2812735) q[0];
sx q[0];
rz(-2.3695282) q[0];
sx q[0];
rz(-2.4738929) q[0];
rz(-0.54155751) q[2];
sx q[2];
rz(-2.1514153) q[2];
sx q[2];
rz(1.970118) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4291961) q[1];
sx q[1];
rz(-1.0496666) q[1];
sx q[1];
rz(2.8217535) q[1];
rz(-pi) q[2];
rz(-1.4164657) q[3];
sx q[3];
rz(-1.4714421) q[3];
sx q[3];
rz(1.1442483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5399897) q[2];
sx q[2];
rz(-1.3301962) q[2];
sx q[2];
rz(-0.23923624) q[2];
rz(0.38419497) q[3];
sx q[3];
rz(-2.4592063) q[3];
sx q[3];
rz(-2.189883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20188986) q[0];
sx q[0];
rz(-1.4677784) q[0];
sx q[0];
rz(1.3772759) q[0];
rz(-1.2205623) q[1];
sx q[1];
rz(-2.2790597) q[1];
sx q[1];
rz(-2.9753704) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3275422) q[0];
sx q[0];
rz(-2.6084427) q[0];
sx q[0];
rz(-2.0943805) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34923415) q[2];
sx q[2];
rz(-2.376261) q[2];
sx q[2];
rz(0.98883307) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4835755) q[1];
sx q[1];
rz(-1.7819575) q[1];
sx q[1];
rz(-1.5394475) q[1];
rz(-pi) q[2];
rz(2.1020426) q[3];
sx q[3];
rz(-2.663718) q[3];
sx q[3];
rz(-0.22051624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7938457) q[2];
sx q[2];
rz(-2.2913427) q[2];
sx q[2];
rz(1.0469077) q[2];
rz(1.8868123) q[3];
sx q[3];
rz(-1.0898277) q[3];
sx q[3];
rz(-1.2966398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60085249) q[0];
sx q[0];
rz(-0.37373251) q[0];
sx q[0];
rz(2.0284213) q[0];
rz(-0.81740776) q[1];
sx q[1];
rz(-1.6956885) q[1];
sx q[1];
rz(2.7730952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93732801) q[0];
sx q[0];
rz(-0.44631347) q[0];
sx q[0];
rz(-2.242779) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7423875) q[2];
sx q[2];
rz(-1.4773122) q[2];
sx q[2];
rz(0.41416083) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.011960192) q[1];
sx q[1];
rz(-0.76597491) q[1];
sx q[1];
rz(0.38928826) q[1];
rz(0.59898563) q[3];
sx q[3];
rz(-2.3040207) q[3];
sx q[3];
rz(-2.4323104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14785279) q[2];
sx q[2];
rz(-0.61823121) q[2];
sx q[2];
rz(-1.7842133) q[2];
rz(-2.3944858) q[3];
sx q[3];
rz(-0.96304572) q[3];
sx q[3];
rz(0.67556206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52637446) q[0];
sx q[0];
rz(-0.45658657) q[0];
sx q[0];
rz(-1.0427465) q[0];
rz(-2.3710947) q[1];
sx q[1];
rz(-1.0105402) q[1];
sx q[1];
rz(-2.3811293) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5805394) q[0];
sx q[0];
rz(-2.7288611) q[0];
sx q[0];
rz(-2.2012635) q[0];
rz(-pi) q[1];
rz(-0.038497849) q[2];
sx q[2];
rz(-1.7506934) q[2];
sx q[2];
rz(2.8556089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.684243) q[1];
sx q[1];
rz(-1.4919229) q[1];
sx q[1];
rz(0.42399391) q[1];
rz(2.4861927) q[3];
sx q[3];
rz(-1.2714579) q[3];
sx q[3];
rz(-1.6884402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9383135) q[2];
sx q[2];
rz(-2.5310897) q[2];
sx q[2];
rz(2.7321613) q[2];
rz(-1.5832541) q[3];
sx q[3];
rz(-0.46054545) q[3];
sx q[3];
rz(-0.62622768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4938477) q[0];
sx q[0];
rz(-1.0418325) q[0];
sx q[0];
rz(-2.7000725) q[0];
rz(-1.5589177) q[1];
sx q[1];
rz(-1.1987004) q[1];
sx q[1];
rz(-0.62335062) q[1];
rz(-1.5123488) q[2];
sx q[2];
rz(-1.6736021) q[2];
sx q[2];
rz(-2.1291669) q[2];
rz(-0.24258976) q[3];
sx q[3];
rz(-1.442083) q[3];
sx q[3];
rz(0.51221893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
