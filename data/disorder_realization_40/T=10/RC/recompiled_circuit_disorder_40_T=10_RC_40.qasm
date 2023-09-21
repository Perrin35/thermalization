OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(-3.0298046) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(4.7987727) q[1];
sx q[1];
rz(6.12943) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.789334) q[0];
sx q[0];
rz(-0.52868045) q[0];
sx q[0];
rz(2.5369011) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0646677) q[2];
sx q[2];
rz(-1.3300606) q[2];
sx q[2];
rz(1.392729) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.44126025) q[1];
sx q[1];
rz(-1.889519) q[1];
sx q[1];
rz(-0.031949921) q[1];
x q[2];
rz(0.33711707) q[3];
sx q[3];
rz(-1.1532591) q[3];
sx q[3];
rz(1.9330213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6979606) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(1.4367746) q[2];
rz(-2.4076961) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88965082) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(0.92457986) q[0];
rz(-2.1444767) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(-1.8181713) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91081496) q[0];
sx q[0];
rz(-1.9998904) q[0];
sx q[0];
rz(-0.50165117) q[0];
x q[1];
rz(0.3496062) q[2];
sx q[2];
rz(-1.5290302) q[2];
sx q[2];
rz(1.8984399) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1393226) q[1];
sx q[1];
rz(-1.3011258) q[1];
sx q[1];
rz(-2.3766999) q[1];
x q[2];
rz(-2.1341392) q[3];
sx q[3];
rz(-0.65348071) q[3];
sx q[3];
rz(-0.37936488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20415846) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(-0.75817529) q[2];
rz(-0.6289064) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(-1.1531856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73308289) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(-2.4531903) q[0];
rz(0.06772659) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(-0.53007954) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353697) q[0];
sx q[0];
rz(-2.831499) q[0];
sx q[0];
rz(1.9276516) q[0];
rz(-pi) q[1];
rz(2.2377551) q[2];
sx q[2];
rz(-0.63637667) q[2];
sx q[2];
rz(0.19665502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9024132) q[1];
sx q[1];
rz(-1.2352408) q[1];
sx q[1];
rz(-2.2526342) q[1];
rz(-pi) q[2];
rz(1.9534555) q[3];
sx q[3];
rz(-2.8842162) q[3];
sx q[3];
rz(-2.5352258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80785859) q[2];
sx q[2];
rz(-3.1299751) q[2];
sx q[2];
rz(-2.2401436) q[2];
rz(2.3060913) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(-1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19105844) q[0];
sx q[0];
rz(-1.5427417) q[0];
sx q[0];
rz(0.78432551) q[0];
rz(0.061231881) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(-0.13664666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2211321) q[0];
sx q[0];
rz(-1.4401299) q[0];
sx q[0];
rz(-2.0544102) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5606974) q[2];
sx q[2];
rz(-0.53933203) q[2];
sx q[2];
rz(-0.5667333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20679684) q[1];
sx q[1];
rz(-3.0691642) q[1];
sx q[1];
rz(-1.2345242) q[1];
rz(0.99478787) q[3];
sx q[3];
rz(-1.4928865) q[3];
sx q[3];
rz(2.5324477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1293929) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(0.56048918) q[2];
rz(-3.1292606) q[3];
sx q[3];
rz(-2.2380232) q[3];
sx q[3];
rz(-2.0509317) q[3];
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
rz(-0.43301582) q[0];
sx q[0];
rz(-0.60537678) q[0];
sx q[0];
rz(-2.3204455) q[0];
rz(-0.87617809) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(1.4076153) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6038937) q[0];
sx q[0];
rz(-1.961381) q[0];
sx q[0];
rz(-1.1917398) q[0];
x q[1];
rz(-2.0166964) q[2];
sx q[2];
rz(-0.6302399) q[2];
sx q[2];
rz(1.557204) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40286139) q[1];
sx q[1];
rz(-2.1364473) q[1];
sx q[1];
rz(2.9823142) q[1];
rz(-1.2110662) q[3];
sx q[3];
rz(-1.779428) q[3];
sx q[3];
rz(-1.2071351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.451482) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(3.0991128) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-2.5094331) q[3];
sx q[3];
rz(0.49155864) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38934389) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(0.4831627) q[0];
rz(-2.0893611) q[1];
sx q[1];
rz(-1.1455043) q[1];
sx q[1];
rz(-2.5767456) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0676346) q[0];
sx q[0];
rz(-1.4356151) q[0];
sx q[0];
rz(2.1002752) q[0];
rz(-1.1499314) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(-1.7771306) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40183345) q[1];
sx q[1];
rz(-1.9976915) q[1];
sx q[1];
rz(0.065211936) q[1];
rz(-pi) q[2];
rz(-0.58638339) q[3];
sx q[3];
rz(-0.27554232) q[3];
sx q[3];
rz(-1.1875718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5715282) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(1.338039) q[2];
rz(1.3048874) q[3];
sx q[3];
rz(-2.0740502) q[3];
sx q[3];
rz(2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682482) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(2.5937953) q[0];
rz(-2.3563747) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(-2.8731667) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9831532) q[0];
sx q[0];
rz(-1.8162677) q[0];
sx q[0];
rz(1.4863187) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7873017) q[2];
sx q[2];
rz(-0.95677081) q[2];
sx q[2];
rz(-1.3348483) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.06936) q[1];
sx q[1];
rz(-1.4727122) q[1];
sx q[1];
rz(0.98274883) q[1];
rz(-pi) q[2];
rz(-0.3018474) q[3];
sx q[3];
rz(-2.0904623) q[3];
sx q[3];
rz(2.2083851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5238374) q[2];
sx q[2];
rz(-0.80344168) q[2];
sx q[2];
rz(-2.3274373) q[2];
rz(2.7653149) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(-0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5557264) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(1.0193753) q[0];
rz(2.2881919) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(-2.6928435) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42496029) q[0];
sx q[0];
rz(-0.53629959) q[0];
sx q[0];
rz(1.1015571) q[0];
rz(1.2690527) q[2];
sx q[2];
rz(-1.8174968) q[2];
sx q[2];
rz(1.4898642) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1467421) q[1];
sx q[1];
rz(-0.59616201) q[1];
sx q[1];
rz(2.6776828) q[1];
x q[2];
rz(-2.429871) q[3];
sx q[3];
rz(-0.91152836) q[3];
sx q[3];
rz(-3.1275415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2723508) q[2];
sx q[2];
rz(-1.7607471) q[2];
sx q[2];
rz(-0.31420079) q[2];
rz(0.82434404) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(2.3468988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070351275) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(1.2605793) q[0];
rz(-2.4977327) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(3.1226645) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2429598) q[0];
sx q[0];
rz(-1.8343933) q[0];
sx q[0];
rz(-1.5522869) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9855965) q[2];
sx q[2];
rz(-1.6678572) q[2];
sx q[2];
rz(0.43924606) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.70871204) q[1];
sx q[1];
rz(-1.5203262) q[1];
sx q[1];
rz(-2.9570079) q[1];
x q[2];
rz(0.67846672) q[3];
sx q[3];
rz(-0.36558357) q[3];
sx q[3];
rz(-1.6775223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54946047) q[2];
sx q[2];
rz(-2.7533054) q[2];
sx q[2];
rz(0.67031676) q[2];
rz(-2.629225) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(1.4204773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5877514) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(2.642139) q[0];
rz(1.5746501) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(-2.0589028) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2912746) q[0];
sx q[0];
rz(-2.5524338) q[0];
sx q[0];
rz(-0.088081443) q[0];
rz(-pi) q[1];
rz(1.6062276) q[2];
sx q[2];
rz(-2.3897768) q[2];
sx q[2];
rz(-2.5978136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4808828) q[1];
sx q[1];
rz(-2.1143267) q[1];
sx q[1];
rz(1.5488312) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4389078) q[3];
sx q[3];
rz(-2.6583238) q[3];
sx q[3];
rz(-1.0815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3988951) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(0.51188525) q[2];
rz(0.39294696) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(-1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95505161) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(-0.74116771) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(2.1096061) q[2];
sx q[2];
rz(-2.451755) q[2];
sx q[2];
rz(-2.3103726) q[2];
rz(2.9560271) q[3];
sx q[3];
rz(-2.3864828) q[3];
sx q[3];
rz(-1.3226487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
