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
rz(-2.3772216) q[0];
sx q[0];
rz(-1.8005014) q[0];
sx q[0];
rz(-2.2361225) q[0];
rz(1.5068997) q[1];
sx q[1];
rz(-2.1807179) q[1];
sx q[1];
rz(1.5009872) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.35957) q[0];
sx q[0];
rz(-1.3728314) q[0];
sx q[0];
rz(2.2819166) q[0];
rz(-1.9328961) q[2];
sx q[2];
rz(-2.675867) q[2];
sx q[2];
rz(-3.1052239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.562362) q[1];
sx q[1];
rz(-0.048431245) q[1];
sx q[1];
rz(-0.50244759) q[1];
rz(-pi) q[2];
rz(-2.7952173) q[3];
sx q[3];
rz(-2.1057671) q[3];
sx q[3];
rz(0.76081027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2291439) q[2];
sx q[2];
rz(-1.392776) q[2];
sx q[2];
rz(2.8196715) q[2];
rz(2.1866482) q[3];
sx q[3];
rz(-0.82986444) q[3];
sx q[3];
rz(1.1436852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8973273) q[0];
sx q[0];
rz(-2.0877512) q[0];
sx q[0];
rz(2.6296997) q[0];
rz(-1.5882675) q[1];
sx q[1];
rz(-0.48738185) q[1];
sx q[1];
rz(2.0194676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.666709) q[0];
sx q[0];
rz(-2.8132952) q[0];
sx q[0];
rz(-2.8216381) q[0];
rz(-pi) q[1];
rz(-1.4830515) q[2];
sx q[2];
rz(-2.7105687) q[2];
sx q[2];
rz(1.4813678) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1074688) q[1];
sx q[1];
rz(-1.9975047) q[1];
sx q[1];
rz(2.772985) q[1];
rz(-pi) q[2];
rz(1.0997222) q[3];
sx q[3];
rz(-2.2212914) q[3];
sx q[3];
rz(-0.6944523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.7684218) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(0.46305099) q[2];
rz(-0.57524663) q[3];
sx q[3];
rz(-1.3653711) q[3];
sx q[3];
rz(0.099253207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.7333882) q[0];
sx q[0];
rz(-2.2938804) q[0];
sx q[0];
rz(-2.7114765) q[0];
rz(-2.6759713) q[1];
sx q[1];
rz(-2.4211113) q[1];
sx q[1];
rz(-0.63708416) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2475605) q[0];
sx q[0];
rz(-1.5848418) q[0];
sx q[0];
rz(1.6117023) q[0];
rz(-0.026224296) q[2];
sx q[2];
rz(-0.78556873) q[2];
sx q[2];
rz(-0.81119591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3780669) q[1];
sx q[1];
rz(-0.77516205) q[1];
sx q[1];
rz(2.7833392) q[1];
x q[2];
rz(-2.5175321) q[3];
sx q[3];
rz(-2.5284323) q[3];
sx q[3];
rz(-1.1414736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35885262) q[2];
sx q[2];
rz(-2.095486) q[2];
sx q[2];
rz(0.72506881) q[2];
rz(-1.3310165) q[3];
sx q[3];
rz(-2.1264117) q[3];
sx q[3];
rz(2.5031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26938874) q[0];
sx q[0];
rz(-1.6870455) q[0];
sx q[0];
rz(-1.4440906) q[0];
rz(-3.0929502) q[1];
sx q[1];
rz(-1.8154058) q[1];
sx q[1];
rz(2.8526502) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5015857) q[0];
sx q[0];
rz(-1.3292392) q[0];
sx q[0];
rz(1.1904741) q[0];
rz(-pi) q[1];
rz(-1.7213351) q[2];
sx q[2];
rz(-1.5427914) q[2];
sx q[2];
rz(1.3893407) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98499456) q[1];
sx q[1];
rz(-0.20594507) q[1];
sx q[1];
rz(1.0099645) q[1];
rz(-pi) q[2];
rz(0.66529243) q[3];
sx q[3];
rz(-1.6136618) q[3];
sx q[3];
rz(-2.8957518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85663969) q[2];
sx q[2];
rz(-0.53264561) q[2];
sx q[2];
rz(-0.3802158) q[2];
rz(0.63816655) q[3];
sx q[3];
rz(-1.9117982) q[3];
sx q[3];
rz(1.1285454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3087092) q[0];
sx q[0];
rz(-0.55279151) q[0];
sx q[0];
rz(2.047245) q[0];
rz(-2.265918) q[1];
sx q[1];
rz(-0.62962571) q[1];
sx q[1];
rz(2.0424776) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5364781) q[0];
sx q[0];
rz(-0.41145675) q[0];
sx q[0];
rz(0.34939975) q[0];
rz(-pi) q[1];
rz(-2.3540456) q[2];
sx q[2];
rz(-2.6733477) q[2];
sx q[2];
rz(0.59216532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.588394) q[1];
sx q[1];
rz(-2.1171682) q[1];
sx q[1];
rz(0.1853892) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91644561) q[3];
sx q[3];
rz(-1.4438085) q[3];
sx q[3];
rz(-0.45230745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8984453) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(-0.97664991) q[2];
rz(-0.53168932) q[3];
sx q[3];
rz(-2.6952126) q[3];
sx q[3];
rz(-0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0584745) q[0];
sx q[0];
rz(-2.0396905) q[0];
sx q[0];
rz(-1.8699159) q[0];
rz(-2.7187128) q[1];
sx q[1];
rz(-1.6115178) q[1];
sx q[1];
rz(-1.5333102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23841183) q[0];
sx q[0];
rz(-0.067771284) q[0];
sx q[0];
rz(1.3890024) q[0];
rz(-pi) q[1];
rz(1.6688231) q[2];
sx q[2];
rz(-1.5801801) q[2];
sx q[2];
rz(2.2414152) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21580885) q[1];
sx q[1];
rz(-1.6263576) q[1];
sx q[1];
rz(-0.30345281) q[1];
rz(-pi) q[2];
rz(-2.3284376) q[3];
sx q[3];
rz(-1.7405563) q[3];
sx q[3];
rz(0.90857279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61234683) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(0.38522729) q[2];
rz(-0.084608229) q[3];
sx q[3];
rz(-0.43360964) q[3];
sx q[3];
rz(-0.87578526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92152921) q[0];
sx q[0];
rz(-1.055287) q[0];
sx q[0];
rz(-2.7440942) q[0];
rz(-0.53783224) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(-0.0040815512) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57848141) q[0];
sx q[0];
rz(-0.77881506) q[0];
sx q[0];
rz(-1.388873) q[0];
rz(-pi) q[1];
rz(-2.053431) q[2];
sx q[2];
rz(-2.3950999) q[2];
sx q[2];
rz(-1.5160402) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7023125) q[1];
sx q[1];
rz(-1.1902307) q[1];
sx q[1];
rz(0.61813942) q[1];
rz(-2.2479731) q[3];
sx q[3];
rz(-1.0056408) q[3];
sx q[3];
rz(0.86097417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6323382) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(-0.21305591) q[2];
rz(0.80900711) q[3];
sx q[3];
rz(-3.1000948) q[3];
sx q[3];
rz(-1.2808778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0179366) q[0];
sx q[0];
rz(-1.2024711) q[0];
sx q[0];
rz(-1.9785471) q[0];
rz(2.842438) q[1];
sx q[1];
rz(-1.9367633) q[1];
sx q[1];
rz(2.1655653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3523063) q[0];
sx q[0];
rz(-1.3255766) q[0];
sx q[0];
rz(-1.6847436) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9891122) q[2];
sx q[2];
rz(-1.3950384) q[2];
sx q[2];
rz(-0.7712785) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4055172) q[1];
sx q[1];
rz(-1.8092522) q[1];
sx q[1];
rz(1.9333436) q[1];
rz(1.8194514) q[3];
sx q[3];
rz(-1.1583405) q[3];
sx q[3];
rz(-1.4506884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8396847) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(-0.1304661) q[2];
rz(1.8179551) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(-2.8470993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94006222) q[0];
sx q[0];
rz(-0.37684965) q[0];
sx q[0];
rz(1.6424204) q[0];
rz(2.6014853) q[1];
sx q[1];
rz(-0.79743782) q[1];
sx q[1];
rz(-1.7971719) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3657065) q[0];
sx q[0];
rz(-1.2822064) q[0];
sx q[0];
rz(-1.1613174) q[0];
rz(-pi) q[1];
rz(-0.70602472) q[2];
sx q[2];
rz(-2.386552) q[2];
sx q[2];
rz(1.4817099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0186179) q[1];
sx q[1];
rz(-0.6280762) q[1];
sx q[1];
rz(1.1372034) q[1];
rz(-pi) q[2];
rz(0.88560652) q[3];
sx q[3];
rz(-0.84784283) q[3];
sx q[3];
rz(2.1849439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53238791) q[2];
sx q[2];
rz(-0.33668533) q[2];
sx q[2];
rz(1.4538291) q[2];
rz(0.42803556) q[3];
sx q[3];
rz(-1.5513709) q[3];
sx q[3];
rz(-1.154703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353272) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(1.4916627) q[0];
rz(0.55117575) q[1];
sx q[1];
rz(-2.1291514) q[1];
sx q[1];
rz(0.59250441) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3640219) q[0];
sx q[0];
rz(-2.7710272) q[0];
sx q[0];
rz(-2.4353566) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5985322) q[2];
sx q[2];
rz(-1.1064227) q[2];
sx q[2];
rz(-2.4891702) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7043276) q[1];
sx q[1];
rz(-2.6710837) q[1];
sx q[1];
rz(-2.5321042) q[1];
rz(1.6208526) q[3];
sx q[3];
rz(-2.069469) q[3];
sx q[3];
rz(-3.0704481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42980117) q[2];
sx q[2];
rz(-1.0415123) q[2];
sx q[2];
rz(-0.05833021) q[2];
rz(0.37060261) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(3.1378194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.8534828) q[0];
sx q[0];
rz(-1.7848889) q[0];
sx q[0];
rz(-3.0369192) q[0];
rz(-1.3660322) q[1];
sx q[1];
rz(-0.80148253) q[1];
sx q[1];
rz(-2.186224) q[1];
rz(-0.38717196) q[2];
sx q[2];
rz(-0.6091112) q[2];
sx q[2];
rz(2.1412639) q[2];
rz(-1.0287063) q[3];
sx q[3];
rz(-2.0301314) q[3];
sx q[3];
rz(0.83740656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
