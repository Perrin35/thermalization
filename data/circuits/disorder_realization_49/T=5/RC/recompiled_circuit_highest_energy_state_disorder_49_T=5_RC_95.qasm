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
rz(3.0316226) q[0];
sx q[0];
rz(-0.87183824) q[0];
sx q[0];
rz(-0.81575552) q[0];
rz(-1.8176796) q[1];
sx q[1];
rz(-1.7426999) q[1];
sx q[1];
rz(0.89371347) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92233673) q[0];
sx q[0];
rz(-1.3748264) q[0];
sx q[0];
rz(-2.6862269) q[0];
rz(-pi) q[1];
rz(-1.4834095) q[2];
sx q[2];
rz(-1.3779216) q[2];
sx q[2];
rz(0.88763104) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1964394) q[1];
sx q[1];
rz(-2.6652539) q[1];
sx q[1];
rz(0.4222541) q[1];
rz(0.23979322) q[3];
sx q[3];
rz(-1.6231281) q[3];
sx q[3];
rz(2.1870942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2414503) q[2];
sx q[2];
rz(-2.9505079) q[2];
sx q[2];
rz(0.07766635) q[2];
rz(-3.0302454) q[3];
sx q[3];
rz(-1.0166549) q[3];
sx q[3];
rz(0.98285037) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0351008) q[0];
sx q[0];
rz(-0.74420559) q[0];
sx q[0];
rz(-2.8700854) q[0];
rz(0.60596451) q[1];
sx q[1];
rz(-2.7022336) q[1];
sx q[1];
rz(-0.6449759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2319658) q[0];
sx q[0];
rz(-1.0459959) q[0];
sx q[0];
rz(2.7746137) q[0];
x q[1];
rz(-2.4753202) q[2];
sx q[2];
rz(-1.4372908) q[2];
sx q[2];
rz(1.324151) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6771705) q[1];
sx q[1];
rz(-1.4289411) q[1];
sx q[1];
rz(-2.9173098) q[1];
x q[2];
rz(-2.9927587) q[3];
sx q[3];
rz(-1.7516802) q[3];
sx q[3];
rz(0.27950865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.490654) q[2];
sx q[2];
rz(-1.9927315) q[2];
sx q[2];
rz(2.8059354) q[2];
rz(0.13441864) q[3];
sx q[3];
rz(-0.69177827) q[3];
sx q[3];
rz(2.5366096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26948872) q[0];
sx q[0];
rz(-1.9600927) q[0];
sx q[0];
rz(0.14542018) q[0];
rz(-0.91436404) q[1];
sx q[1];
rz(-1.7040355) q[1];
sx q[1];
rz(-2.5252555) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59728226) q[0];
sx q[0];
rz(-2.4203886) q[0];
sx q[0];
rz(0.82996675) q[0];
rz(-pi) q[1];
rz(-3.0923244) q[2];
sx q[2];
rz(-1.0596837) q[2];
sx q[2];
rz(3.0126103) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80464493) q[1];
sx q[1];
rz(-1.8700127) q[1];
sx q[1];
rz(0.3693337) q[1];
x q[2];
rz(2.8313274) q[3];
sx q[3];
rz(-2.4832926) q[3];
sx q[3];
rz(-2.8993115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99438897) q[2];
sx q[2];
rz(-2.2115967) q[2];
sx q[2];
rz(0.27929107) q[2];
rz(2.764737) q[3];
sx q[3];
rz(-0.91747228) q[3];
sx q[3];
rz(0.74523029) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.299861) q[0];
sx q[0];
rz(-2.0722516) q[0];
sx q[0];
rz(1.8002864) q[0];
rz(0.74814859) q[1];
sx q[1];
rz(-2.61863) q[1];
sx q[1];
rz(1.062692) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6338017) q[0];
sx q[0];
rz(-1.6865382) q[0];
sx q[0];
rz(2.1383889) q[0];
x q[1];
rz(0.80047582) q[2];
sx q[2];
rz(-1.3409586) q[2];
sx q[2];
rz(0.92170131) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53790435) q[1];
sx q[1];
rz(-2.474437) q[1];
sx q[1];
rz(2.0869419) q[1];
x q[2];
rz(2.483401) q[3];
sx q[3];
rz(-0.88630967) q[3];
sx q[3];
rz(1.6940862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5144389) q[2];
sx q[2];
rz(-1.2320765) q[2];
sx q[2];
rz(0.76986924) q[2];
rz(1.7047966) q[3];
sx q[3];
rz(-1.0154279) q[3];
sx q[3];
rz(2.3427486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1771667) q[0];
sx q[0];
rz(-1.3097958) q[0];
sx q[0];
rz(-0.31846309) q[0];
rz(0.98694363) q[1];
sx q[1];
rz(-1.3401597) q[1];
sx q[1];
rz(-0.14652227) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5269546) q[0];
sx q[0];
rz(-1.2181158) q[0];
sx q[0];
rz(2.1988064) q[0];
x q[1];
rz(2.179903) q[2];
sx q[2];
rz(-2.1868248) q[2];
sx q[2];
rz(-2.2005657) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95899907) q[1];
sx q[1];
rz(-1.800391) q[1];
sx q[1];
rz(2.793538) q[1];
rz(-pi) q[2];
rz(-1.2066483) q[3];
sx q[3];
rz(-0.86508646) q[3];
sx q[3];
rz(-0.067506703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5557308) q[2];
sx q[2];
rz(-2.9572058) q[2];
sx q[2];
rz(1.8905852) q[2];
rz(2.3986554) q[3];
sx q[3];
rz(-1.5115967) q[3];
sx q[3];
rz(1.2062629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.10155216) q[0];
sx q[0];
rz(-2.0843625) q[0];
sx q[0];
rz(-1.2054766) q[0];
rz(2.7936753) q[1];
sx q[1];
rz(-1.9708865) q[1];
sx q[1];
rz(-1.233137) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83506993) q[0];
sx q[0];
rz(-0.26598922) q[0];
sx q[0];
rz(2.8138922) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80329956) q[2];
sx q[2];
rz(-1.1931007) q[2];
sx q[2];
rz(-0.74716728) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5812787) q[1];
sx q[1];
rz(-0.73773161) q[1];
sx q[1];
rz(-0.88256617) q[1];
rz(-pi) q[2];
rz(2.5830501) q[3];
sx q[3];
rz(-2.5457446) q[3];
sx q[3];
rz(2.5003025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7463688) q[2];
sx q[2];
rz(-0.93337983) q[2];
sx q[2];
rz(-2.7044435) q[2];
rz(2.7002667) q[3];
sx q[3];
rz(-0.91259846) q[3];
sx q[3];
rz(0.92379409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8298892) q[0];
sx q[0];
rz(-1.7376816) q[0];
sx q[0];
rz(-0.29790685) q[0];
rz(0.20826134) q[1];
sx q[1];
rz(-0.88349897) q[1];
sx q[1];
rz(-0.40685245) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1778616) q[0];
sx q[0];
rz(-3.047982) q[0];
sx q[0];
rz(-1.2340182) q[0];
rz(2.9003606) q[2];
sx q[2];
rz(-2.2451984) q[2];
sx q[2];
rz(0.33283744) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6089203) q[1];
sx q[1];
rz(-0.77076503) q[1];
sx q[1];
rz(-1.1966329) q[1];
x q[2];
rz(-1.0034542) q[3];
sx q[3];
rz(-1.1441244) q[3];
sx q[3];
rz(2.2578491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3901861) q[2];
sx q[2];
rz(-2.4441256) q[2];
sx q[2];
rz(-2.3263399) q[2];
rz(-0.32663545) q[3];
sx q[3];
rz(-1.9950461) q[3];
sx q[3];
rz(-3.1281285) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1470452) q[0];
sx q[0];
rz(-0.70962405) q[0];
sx q[0];
rz(-2.5780594) q[0];
rz(1.4348449) q[1];
sx q[1];
rz(-0.98054612) q[1];
sx q[1];
rz(0.30418667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4833058) q[0];
sx q[0];
rz(-2.493447) q[0];
sx q[0];
rz(3.0850406) q[0];
rz(-pi) q[1];
rz(-0.86433347) q[2];
sx q[2];
rz(-1.1030518) q[2];
sx q[2];
rz(0.36004972) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4106928) q[1];
sx q[1];
rz(-2.5531156) q[1];
sx q[1];
rz(1.9517577) q[1];
rz(-pi) q[2];
rz(-2.9161737) q[3];
sx q[3];
rz(-1.9758164) q[3];
sx q[3];
rz(1.5532358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7131416) q[2];
sx q[2];
rz(-2.9623172) q[2];
sx q[2];
rz(-0.98503867) q[2];
rz(-1.8200412) q[3];
sx q[3];
rz(-1.225084) q[3];
sx q[3];
rz(-0.2329181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.709885) q[0];
sx q[0];
rz(-0.49254492) q[0];
sx q[0];
rz(-2.7370969) q[0];
rz(2.3996023) q[1];
sx q[1];
rz(-1.208035) q[1];
sx q[1];
rz(-2.580339) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1741176) q[0];
sx q[0];
rz(-2.1537279) q[0];
sx q[0];
rz(-2.6487104) q[0];
rz(-pi) q[1];
rz(1.1001585) q[2];
sx q[2];
rz(-1.2377035) q[2];
sx q[2];
rz(-1.3317219) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2798279) q[1];
sx q[1];
rz(-0.4273673) q[1];
sx q[1];
rz(0.98160569) q[1];
rz(-pi) q[2];
rz(-0.75144479) q[3];
sx q[3];
rz(-1.8490233) q[3];
sx q[3];
rz(3.1415526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.81885091) q[2];
sx q[2];
rz(-0.35917869) q[2];
sx q[2];
rz(2.2364869) q[2];
rz(0.94308606) q[3];
sx q[3];
rz(-1.9102996) q[3];
sx q[3];
rz(-2.4020307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(-0.99317011) q[0];
sx q[0];
rz(-0.44023308) q[0];
sx q[0];
rz(0.37047186) q[0];
rz(0.91520339) q[1];
sx q[1];
rz(-1.4280495) q[1];
sx q[1];
rz(-2.3781093) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3199056) q[0];
sx q[0];
rz(-2.3662462) q[0];
sx q[0];
rz(-0.5628235) q[0];
rz(1.8862091) q[2];
sx q[2];
rz(-2.1120781) q[2];
sx q[2];
rz(-0.72555055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6333446) q[1];
sx q[1];
rz(-0.92735433) q[1];
sx q[1];
rz(0.28336759) q[1];
rz(-pi) q[2];
rz(0.053289847) q[3];
sx q[3];
rz(-1.6052979) q[3];
sx q[3];
rz(-0.49279398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8107599) q[2];
sx q[2];
rz(-1.5786889) q[2];
sx q[2];
rz(-1.8806489) q[2];
rz(-1.8597982) q[3];
sx q[3];
rz(-2.1112879) q[3];
sx q[3];
rz(-1.1873881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8967165) q[0];
sx q[0];
rz(-1.5302932) q[0];
sx q[0];
rz(1.2320919) q[0];
rz(-2.2121519) q[1];
sx q[1];
rz(-0.96439958) q[1];
sx q[1];
rz(0.31698116) q[1];
rz(1.4030761) q[2];
sx q[2];
rz(-0.3315331) q[2];
sx q[2];
rz(-0.49102993) q[2];
rz(1.9396325) q[3];
sx q[3];
rz(-1.0452253) q[3];
sx q[3];
rz(-0.54946231) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
