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
rz(-0.10997009) q[0];
sx q[0];
rz(0.87183824) q[0];
sx q[0];
rz(11.750615) q[0];
rz(1.323913) q[1];
sx q[1];
rz(-1.3988928) q[1];
sx q[1];
rz(-0.89371347) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1146671) q[0];
sx q[0];
rz(-2.6485778) q[0];
sx q[0];
rz(0.42401029) q[0];
rz(-pi) q[1];
rz(-2.9479974) q[2];
sx q[2];
rz(-1.6565588) q[2];
sx q[2];
rz(-0.69995689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47706931) q[1];
sx q[1];
rz(-2.0023195) q[1];
sx q[1];
rz(1.7791788) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2170691) q[3];
sx q[3];
rz(-0.24532977) q[3];
sx q[3];
rz(-2.7360578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90014234) q[2];
sx q[2];
rz(-0.19108471) q[2];
sx q[2];
rz(-3.0639263) q[2];
rz(0.11134722) q[3];
sx q[3];
rz(-1.0166549) q[3];
sx q[3];
rz(0.98285037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10649189) q[0];
sx q[0];
rz(-0.74420559) q[0];
sx q[0];
rz(-0.27150723) q[0];
rz(-2.5356281) q[1];
sx q[1];
rz(-2.7022336) q[1];
sx q[1];
rz(-0.6449759) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6124961) q[0];
sx q[0];
rz(-1.2550651) q[0];
sx q[0];
rz(1.015618) q[0];
rz(-2.4753202) q[2];
sx q[2];
rz(-1.4372908) q[2];
sx q[2];
rz(1.324151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6771705) q[1];
sx q[1];
rz(-1.7126516) q[1];
sx q[1];
rz(-2.9173098) q[1];
rz(-2.2521001) q[3];
sx q[3];
rz(-2.9078662) q[3];
sx q[3];
rz(0.97433486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6509387) q[2];
sx q[2];
rz(-1.9927315) q[2];
sx q[2];
rz(-2.8059354) q[2];
rz(-3.007174) q[3];
sx q[3];
rz(-0.69177827) q[3];
sx q[3];
rz(2.5366096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.26948872) q[0];
sx q[0];
rz(-1.9600927) q[0];
sx q[0];
rz(-0.14542018) q[0];
rz(0.91436404) q[1];
sx q[1];
rz(-1.7040355) q[1];
sx q[1];
rz(-0.61633715) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37164524) q[0];
sx q[0];
rz(-2.0326737) q[0];
sx q[0];
rz(0.99528639) q[0];
rz(-1.6583865) q[2];
sx q[2];
rz(-2.6283205) q[2];
sx q[2];
rz(0.22944726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7246252) q[1];
sx q[1];
rz(-0.47096241) q[1];
sx q[1];
rz(2.4344786) q[1];
rz(-pi) q[2];
rz(1.3389204) q[3];
sx q[3];
rz(-2.1926741) q[3];
sx q[3];
rz(-2.9987914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99438897) q[2];
sx q[2];
rz(-2.2115967) q[2];
sx q[2];
rz(-0.27929107) q[2];
rz(0.3768557) q[3];
sx q[3];
rz(-2.2241204) q[3];
sx q[3];
rz(0.74523029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.299861) q[0];
sx q[0];
rz(-2.0722516) q[0];
sx q[0];
rz(1.8002864) q[0];
rz(-0.74814859) q[1];
sx q[1];
rz(-0.52296269) q[1];
sx q[1];
rz(1.062692) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11618075) q[0];
sx q[0];
rz(-2.5635911) q[0];
sx q[0];
rz(1.7837753) q[0];
x q[1];
rz(0.80047582) q[2];
sx q[2];
rz(-1.8006341) q[2];
sx q[2];
rz(2.2198913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.53790435) q[1];
sx q[1];
rz(-2.474437) q[1];
sx q[1];
rz(1.0546507) q[1];
x q[2];
rz(2.371754) q[3];
sx q[3];
rz(-2.0645118) q[3];
sx q[3];
rz(-0.57800402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5144389) q[2];
sx q[2];
rz(-1.2320765) q[2];
sx q[2];
rz(-2.3717234) q[2];
rz(-1.7047966) q[3];
sx q[3];
rz(-1.0154279) q[3];
sx q[3];
rz(-2.3427486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-1.1771667) q[0];
sx q[0];
rz(-1.3097958) q[0];
sx q[0];
rz(0.31846309) q[0];
rz(0.98694363) q[1];
sx q[1];
rz(-1.3401597) q[1];
sx q[1];
rz(2.9950704) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28956911) q[0];
sx q[0];
rz(-0.98678723) q[0];
sx q[0];
rz(-2.7147074) q[0];
x q[1];
rz(-2.4618838) q[2];
sx q[2];
rz(-2.3041937) q[2];
sx q[2];
rz(-1.8200761) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9695786) q[1];
sx q[1];
rz(-2.7272228) q[1];
sx q[1];
rz(-2.5408391) q[1];
x q[2];
rz(0.73934619) q[3];
sx q[3];
rz(-1.2962747) q[3];
sx q[3];
rz(1.8806461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5557308) q[2];
sx q[2];
rz(-0.18438688) q[2];
sx q[2];
rz(-1.2510074) q[2];
rz(-2.3986554) q[3];
sx q[3];
rz(-1.5115967) q[3];
sx q[3];
rz(1.9353297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10155216) q[0];
sx q[0];
rz(-1.0572301) q[0];
sx q[0];
rz(1.2054766) q[0];
rz(-0.34791738) q[1];
sx q[1];
rz(-1.1707062) q[1];
sx q[1];
rz(-1.9084557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41878149) q[0];
sx q[0];
rz(-1.6555047) q[0];
sx q[0];
rz(-2.8891488) q[0];
rz(2.3382931) q[2];
sx q[2];
rz(-1.9484919) q[2];
sx q[2];
rz(2.3944254) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4193488) q[1];
sx q[1];
rz(-2.1170685) q[1];
sx q[1];
rz(-2.6180079) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5830501) q[3];
sx q[3];
rz(-2.5457446) q[3];
sx q[3];
rz(-2.5003025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7463688) q[2];
sx q[2];
rz(-0.93337983) q[2];
sx q[2];
rz(0.4371492) q[2];
rz(0.44132597) q[3];
sx q[3];
rz(-2.2289942) q[3];
sx q[3];
rz(0.92379409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31170347) q[0];
sx q[0];
rz(-1.7376816) q[0];
sx q[0];
rz(0.29790685) q[0];
rz(-0.20826134) q[1];
sx q[1];
rz(-0.88349897) q[1];
sx q[1];
rz(-2.7347402) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132459) q[0];
sx q[0];
rz(-1.6016895) q[0];
sx q[0];
rz(1.6591765) q[0];
rz(-0.24123206) q[2];
sx q[2];
rz(-2.2451984) q[2];
sx q[2];
rz(-2.8087552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.033483) q[1];
sx q[1];
rz(-0.86520608) q[1];
sx q[1];
rz(-2.8005141) q[1];
rz(2.6471944) q[3];
sx q[3];
rz(-1.0596078) q[3];
sx q[3];
rz(-2.7123812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75140658) q[2];
sx q[2];
rz(-2.4441256) q[2];
sx q[2];
rz(-0.81525272) q[2];
rz(-0.32663545) q[3];
sx q[3];
rz(-1.1465466) q[3];
sx q[3];
rz(-0.013464125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99454749) q[0];
sx q[0];
rz(-2.4319686) q[0];
sx q[0];
rz(0.56353322) q[0];
rz(-1.4348449) q[1];
sx q[1];
rz(-0.98054612) q[1];
sx q[1];
rz(2.837406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72918159) q[0];
sx q[0];
rz(-0.92386073) q[0];
sx q[0];
rz(1.6135733) q[0];
rz(-pi) q[1];
rz(2.2320692) q[2];
sx q[2];
rz(-0.82445383) q[2];
sx q[2];
rz(-0.72469358) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4106928) q[1];
sx q[1];
rz(-0.58847702) q[1];
sx q[1];
rz(1.9517577) q[1];
rz(-pi) q[2];
rz(-1.1564141) q[3];
sx q[3];
rz(-1.3638921) q[3];
sx q[3];
rz(0.072553886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7131416) q[2];
sx q[2];
rz(-0.17927543) q[2];
sx q[2];
rz(-2.156554) q[2];
rz(1.3215514) q[3];
sx q[3];
rz(-1.225084) q[3];
sx q[3];
rz(2.9086746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.709885) q[0];
sx q[0];
rz(-0.49254492) q[0];
sx q[0];
rz(-2.7370969) q[0];
rz(-2.3996023) q[1];
sx q[1];
rz(-1.9335577) q[1];
sx q[1];
rz(-2.580339) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7405734) q[0];
sx q[0];
rz(-2.3971746) q[0];
sx q[0];
rz(-2.1932426) q[0];
rz(-pi) q[1];
rz(2.77131) q[2];
sx q[2];
rz(-2.0136626) q[2];
sx q[2];
rz(0.40391573) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8617647) q[1];
sx q[1];
rz(-0.4273673) q[1];
sx q[1];
rz(-2.159987) q[1];
rz(1.9434378) q[3];
sx q[3];
rz(-2.2868881) q[3];
sx q[3];
rz(1.3195697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.81885091) q[2];
sx q[2];
rz(-2.782414) q[2];
sx q[2];
rz(0.90510577) q[2];
rz(-2.1985066) q[3];
sx q[3];
rz(-1.9102996) q[3];
sx q[3];
rz(-2.4020307) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1484225) q[0];
sx q[0];
rz(-0.44023308) q[0];
sx q[0];
rz(0.37047186) q[0];
rz(-2.2263893) q[1];
sx q[1];
rz(-1.4280495) q[1];
sx q[1];
rz(-2.3781093) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5452427) q[0];
sx q[0];
rz(-2.204329) q[0];
sx q[0];
rz(2.0526396) q[0];
rz(-2.6652226) q[2];
sx q[2];
rz(-2.5231276) q[2];
sx q[2];
rz(-2.9805984) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6333446) q[1];
sx q[1];
rz(-2.2142383) q[1];
sx q[1];
rz(0.28336759) q[1];
rz(-pi) q[2];
rz(1.6053469) q[3];
sx q[3];
rz(-1.5175382) q[3];
sx q[3];
rz(-2.0617503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8107599) q[2];
sx q[2];
rz(-1.5786889) q[2];
sx q[2];
rz(-1.2609437) q[2];
rz(1.2817945) q[3];
sx q[3];
rz(-2.1112879) q[3];
sx q[3];
rz(1.9542046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24487615) q[0];
sx q[0];
rz(-1.6112994) q[0];
sx q[0];
rz(-1.9095008) q[0];
rz(0.92944073) q[1];
sx q[1];
rz(-0.96439958) q[1];
sx q[1];
rz(0.31698116) q[1];
rz(-1.2435883) q[2];
sx q[2];
rz(-1.6251593) q[2];
sx q[2];
rz(-1.9030824) q[2];
rz(-1.9396325) q[3];
sx q[3];
rz(-2.0963674) q[3];
sx q[3];
rz(2.5921303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
