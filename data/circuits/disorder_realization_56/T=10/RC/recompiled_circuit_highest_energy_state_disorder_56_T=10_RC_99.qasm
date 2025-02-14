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
rz(-0.40060488) q[0];
sx q[0];
rz(-2.7317943) q[0];
sx q[0];
rz(2.0706489) q[0];
rz(0.61869705) q[1];
sx q[1];
rz(3.7096042) q[1];
sx q[1];
rz(8.5811442) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5045568) q[0];
sx q[0];
rz(-1.3063653) q[0];
sx q[0];
rz(-0.030585551) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4272408) q[2];
sx q[2];
rz(-1.2782405) q[2];
sx q[2];
rz(-1.0247599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2045793) q[1];
sx q[1];
rz(-2.3790303) q[1];
sx q[1];
rz(2.049905) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0655965) q[3];
sx q[3];
rz(-1.9952979) q[3];
sx q[3];
rz(-0.44874661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3598651) q[2];
sx q[2];
rz(-0.83911506) q[2];
sx q[2];
rz(-3.1245933) q[2];
rz(1.7403691) q[3];
sx q[3];
rz(-2.5798116) q[3];
sx q[3];
rz(-1.9248272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725937) q[0];
sx q[0];
rz(-1.8091135) q[0];
sx q[0];
rz(-2.3542985) q[0];
rz(-0.17838082) q[1];
sx q[1];
rz(-1.2063824) q[1];
sx q[1];
rz(-1.23752) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4263692) q[0];
sx q[0];
rz(-2.5551717) q[0];
sx q[0];
rz(-0.1118025) q[0];
rz(-pi) q[1];
rz(-1.9745047) q[2];
sx q[2];
rz(-2.4879527) q[2];
sx q[2];
rz(-2.6622651) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.54531714) q[1];
sx q[1];
rz(-1.8870237) q[1];
sx q[1];
rz(-1.113722) q[1];
rz(-pi) q[2];
rz(-0.44053733) q[3];
sx q[3];
rz(-1.8228662) q[3];
sx q[3];
rz(2.1217176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7201207) q[2];
sx q[2];
rz(-1.6933491) q[2];
sx q[2];
rz(-0.063701542) q[2];
rz(-1.2740159) q[3];
sx q[3];
rz(-2.2985022) q[3];
sx q[3];
rz(-1.564285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.098123) q[0];
sx q[0];
rz(-1.948057) q[0];
sx q[0];
rz(-0.7435588) q[0];
rz(-2.6877563) q[1];
sx q[1];
rz(-1.791714) q[1];
sx q[1];
rz(-2.7226864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48506698) q[0];
sx q[0];
rz(-0.48975268) q[0];
sx q[0];
rz(1.144125) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9100283) q[2];
sx q[2];
rz(-2.3776109) q[2];
sx q[2];
rz(2.6126249) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53330227) q[1];
sx q[1];
rz(-2.4455297) q[1];
sx q[1];
rz(-1.1185557) q[1];
rz(1.4019074) q[3];
sx q[3];
rz(-2.4554376) q[3];
sx q[3];
rz(-2.9952733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0539187) q[2];
sx q[2];
rz(-2.3232465) q[2];
sx q[2];
rz(1.9107001) q[2];
rz(2.6026717) q[3];
sx q[3];
rz(-1.475622) q[3];
sx q[3];
rz(-1.6787136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8777799) q[0];
sx q[0];
rz(-2.3872264) q[0];
sx q[0];
rz(-0.60212773) q[0];
rz(-1.4471794) q[1];
sx q[1];
rz(-2.4592631) q[1];
sx q[1];
rz(-2.2785861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8558288) q[0];
sx q[0];
rz(-1.4572954) q[0];
sx q[0];
rz(2.4185926) q[0];
x q[1];
rz(-0.90410467) q[2];
sx q[2];
rz(-0.46510425) q[2];
sx q[2];
rz(0.99921528) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.132024) q[1];
sx q[1];
rz(-1.3088262) q[1];
sx q[1];
rz(0.97673194) q[1];
rz(-pi) q[2];
rz(1.9163777) q[3];
sx q[3];
rz(-2.3732819) q[3];
sx q[3];
rz(-2.5490724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0756695) q[2];
sx q[2];
rz(-1.9663591) q[2];
sx q[2];
rz(-2.2389331) q[2];
rz(0.81501189) q[3];
sx q[3];
rz(-2.5383526) q[3];
sx q[3];
rz(-2.7284315) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0478504) q[0];
sx q[0];
rz(-3.0939026) q[0];
sx q[0];
rz(-0.93511859) q[0];
rz(1.5330261) q[1];
sx q[1];
rz(-0.32564274) q[1];
sx q[1];
rz(2.8757222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8489129) q[0];
sx q[0];
rz(-2.8124708) q[0];
sx q[0];
rz(-1.9051244) q[0];
rz(-pi) q[1];
rz(-2.3912848) q[2];
sx q[2];
rz(-1.3641832) q[2];
sx q[2];
rz(-0.37918249) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26755493) q[1];
sx q[1];
rz(-1.7312164) q[1];
sx q[1];
rz(2.7880413) q[1];
rz(-pi) q[2];
rz(1.9343801) q[3];
sx q[3];
rz(-0.28237469) q[3];
sx q[3];
rz(2.7350712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68044668) q[2];
sx q[2];
rz(-1.5004044) q[2];
sx q[2];
rz(0.45073304) q[2];
rz(1.1905253) q[3];
sx q[3];
rz(-0.7014941) q[3];
sx q[3];
rz(0.43959555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.180069) q[0];
sx q[0];
rz(-2.838205) q[0];
sx q[0];
rz(-3.0971089) q[0];
rz(-2.8728409) q[1];
sx q[1];
rz(-2.2046397) q[1];
sx q[1];
rz(-0.43089795) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92545729) q[0];
sx q[0];
rz(-1.3089942) q[0];
sx q[0];
rz(0.63762224) q[0];
x q[1];
rz(1.0111647) q[2];
sx q[2];
rz(-0.75087386) q[2];
sx q[2];
rz(2.4023285) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9932672) q[1];
sx q[1];
rz(-1.1701411) q[1];
sx q[1];
rz(2.2877778) q[1];
x q[2];
rz(-2.481302) q[3];
sx q[3];
rz(-1.5753645) q[3];
sx q[3];
rz(-0.82483236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1162794) q[2];
sx q[2];
rz(-0.37386027) q[2];
sx q[2];
rz(-0.6244134) q[2];
rz(-2.0404909) q[3];
sx q[3];
rz(-0.81239429) q[3];
sx q[3];
rz(-2.1541434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2295912) q[0];
sx q[0];
rz(-2.0822552) q[0];
sx q[0];
rz(0.68369317) q[0];
rz(1.4423485) q[1];
sx q[1];
rz(-2.3577299) q[1];
sx q[1];
rz(2.9396465) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6091222) q[0];
sx q[0];
rz(-1.2426071) q[0];
sx q[0];
rz(-0.0013602982) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82345404) q[2];
sx q[2];
rz(-1.5353893) q[2];
sx q[2];
rz(-0.096114352) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8117042) q[1];
sx q[1];
rz(-2.0913908) q[1];
sx q[1];
rz(1.64774) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0544087) q[3];
sx q[3];
rz(-2.6823061) q[3];
sx q[3];
rz(2.313314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0553637) q[2];
sx q[2];
rz(-1.4198885) q[2];
sx q[2];
rz(-1.8602547) q[2];
rz(0.66644871) q[3];
sx q[3];
rz(-2.0446348) q[3];
sx q[3];
rz(-0.54653978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619038) q[0];
sx q[0];
rz(-2.5894916) q[0];
sx q[0];
rz(-2.6276278) q[0];
rz(-2.1886096) q[1];
sx q[1];
rz(-2.516808) q[1];
sx q[1];
rz(-2.0228588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054487451) q[0];
sx q[0];
rz(-1.4892254) q[0];
sx q[0];
rz(-2.7177627) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4301901) q[2];
sx q[2];
rz(-0.83900827) q[2];
sx q[2];
rz(0.48626394) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8427022) q[1];
sx q[1];
rz(-1.7849079) q[1];
sx q[1];
rz(1.1529117) q[1];
rz(2.5947078) q[3];
sx q[3];
rz(-0.62084475) q[3];
sx q[3];
rz(-1.3399233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.755456) q[2];
sx q[2];
rz(-0.35217199) q[2];
sx q[2];
rz(0.90292162) q[2];
rz(1.6821945) q[3];
sx q[3];
rz(-1.7995588) q[3];
sx q[3];
rz(-2.3193147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67895401) q[0];
sx q[0];
rz(-2.8215388) q[0];
sx q[0];
rz(-1.1902887) q[0];
rz(0.32084385) q[1];
sx q[1];
rz(-1.4535934) q[1];
sx q[1];
rz(1.2215337) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7090593) q[0];
sx q[0];
rz(-1.7072868) q[0];
sx q[0];
rz(0.20705072) q[0];
x q[1];
rz(1.874711) q[2];
sx q[2];
rz(-2.3836145) q[2];
sx q[2];
rz(2.9044819) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.62432837) q[1];
sx q[1];
rz(-2.359694) q[1];
sx q[1];
rz(-2.8502591) q[1];
x q[2];
rz(-0.19467312) q[3];
sx q[3];
rz(-1.3910196) q[3];
sx q[3];
rz(-0.67357066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2373576) q[2];
sx q[2];
rz(-2.9866437) q[2];
sx q[2];
rz(-2.8083189) q[2];
rz(2.1244369) q[3];
sx q[3];
rz(-2.1867496) q[3];
sx q[3];
rz(2.8216383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.3453813) q[0];
sx q[0];
rz(-0.51207241) q[0];
sx q[0];
rz(-0.91878015) q[0];
rz(2.9504919) q[1];
sx q[1];
rz(-2.7156576) q[1];
sx q[1];
rz(0.79413116) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16852233) q[0];
sx q[0];
rz(-2.2380073) q[0];
sx q[0];
rz(-2.4486827) q[0];
rz(-2.625611) q[2];
sx q[2];
rz(-2.7323639) q[2];
sx q[2];
rz(0.87142631) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9136466) q[1];
sx q[1];
rz(-1.8381048) q[1];
sx q[1];
rz(-2.2993907) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32148949) q[3];
sx q[3];
rz(-1.7712542) q[3];
sx q[3];
rz(3.0908302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0309151) q[2];
sx q[2];
rz(-0.6157178) q[2];
sx q[2];
rz(-2.3827629) q[2];
rz(-2.2714254) q[3];
sx q[3];
rz(-2.3535959) q[3];
sx q[3];
rz(-2.0738535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(3.1398685) q[0];
sx q[0];
rz(-1.4586466) q[0];
sx q[0];
rz(1.9139105) q[0];
rz(-0.43791804) q[1];
sx q[1];
rz(-1.7057849) q[1];
sx q[1];
rz(1.5954856) q[1];
rz(2.6545637) q[2];
sx q[2];
rz(-1.1985881) q[2];
sx q[2];
rz(2.3771622) q[2];
rz(-1.036676) q[3];
sx q[3];
rz(-2.277959) q[3];
sx q[3];
rz(-2.1716519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
