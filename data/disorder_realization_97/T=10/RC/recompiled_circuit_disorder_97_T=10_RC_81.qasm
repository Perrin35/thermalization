OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(-1.6516049) q[0];
sx q[0];
rz(-2.2111501) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(4.2760744) q[1];
sx q[1];
rz(8.3174336) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21270277) q[0];
sx q[0];
rz(-1.2455997) q[0];
sx q[0];
rz(-2.4341499) q[0];
x q[1];
rz(2.2912575) q[2];
sx q[2];
rz(-1.9413661) q[2];
sx q[2];
rz(2.4239899) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2335637) q[1];
sx q[1];
rz(-1.4587914) q[1];
sx q[1];
rz(-1.2909375) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1125621) q[3];
sx q[3];
rz(-0.68800612) q[3];
sx q[3];
rz(-1.0109166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.063623108) q[2];
sx q[2];
rz(-2.412553) q[2];
sx q[2];
rz(-1.8135653) q[2];
rz(2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.48224738) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(-0.88062084) q[0];
rz(1.847514) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(2.3243288) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9416695) q[0];
sx q[0];
rz(-0.37460735) q[0];
sx q[0];
rz(1.1767715) q[0];
rz(-1.1772637) q[2];
sx q[2];
rz(-0.52429188) q[2];
sx q[2];
rz(0.92698586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0751794) q[1];
sx q[1];
rz(-1.315409) q[1];
sx q[1];
rz(-0.65402072) q[1];
rz(1.9278139) q[3];
sx q[3];
rz(-0.5001874) q[3];
sx q[3];
rz(-0.57750765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91784224) q[2];
sx q[2];
rz(-0.68085256) q[2];
sx q[2];
rz(-0.36402738) q[2];
rz(-2.1552127) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(1.4646437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7746975) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(2.1752775) q[0];
rz(0.19293383) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(1.4470709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23611785) q[0];
sx q[0];
rz(-0.21572278) q[0];
sx q[0];
rz(-0.082213684) q[0];
x q[1];
rz(-0.082712163) q[2];
sx q[2];
rz(-1.2809922) q[2];
sx q[2];
rz(-2.2207584) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9008873) q[1];
sx q[1];
rz(-1.9449688) q[1];
sx q[1];
rz(-2.4918873) q[1];
rz(-1.9871739) q[3];
sx q[3];
rz(-1.795711) q[3];
sx q[3];
rz(-0.51605663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43859279) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(2.036371) q[2];
rz(0.74622074) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(-0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.1352017) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(1.6756469) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-2.7526061) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8303191) q[0];
sx q[0];
rz(-1.8000326) q[0];
sx q[0];
rz(-2.8515408) q[0];
rz(-pi) q[1];
rz(-0.19548266) q[2];
sx q[2];
rz(-1.0921548) q[2];
sx q[2];
rz(0.10749707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7547238) q[1];
sx q[1];
rz(-1.9935973) q[1];
sx q[1];
rz(2.5368284) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9205434) q[3];
sx q[3];
rz(-2.3813558) q[3];
sx q[3];
rz(-1.5787214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-0.45670613) q[2];
rz(-1.6263973) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91963768) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(2.7815681) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(-0.49450758) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7617103) q[0];
sx q[0];
rz(-2.8344791) q[0];
sx q[0];
rz(-2.7193927) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.102166) q[2];
sx q[2];
rz(-1.4354424) q[2];
sx q[2];
rz(2.6089422) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8620017) q[1];
sx q[1];
rz(-1.630736) q[1];
sx q[1];
rz(-0.63043352) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0881808) q[3];
sx q[3];
rz(-0.44465372) q[3];
sx q[3];
rz(-0.85944552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4313724) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(2.8430856) q[2];
rz(-0.31202894) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(-0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1453778) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(-0.19590713) q[0];
rz(-0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(1.235199) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93341953) q[0];
sx q[0];
rz(-2.6801077) q[0];
sx q[0];
rz(0.74605201) q[0];
x q[1];
rz(0.39288315) q[2];
sx q[2];
rz(-0.74379951) q[2];
sx q[2];
rz(1.5647897) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27855733) q[1];
sx q[1];
rz(-2.1457991) q[1];
sx q[1];
rz(-2.0726191) q[1];
x q[2];
rz(-2.9052832) q[3];
sx q[3];
rz(-1.6522539) q[3];
sx q[3];
rz(0.27796516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49332508) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(-2.7098999) q[2];
rz(-1.4124983) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(-3.0055962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39111185) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(0.11211638) q[0];
rz(-2.926459) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(-2.0281866) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8047785) q[0];
sx q[0];
rz(-1.3114197) q[0];
sx q[0];
rz(2.74385) q[0];
rz(-pi) q[1];
rz(-2.8261975) q[2];
sx q[2];
rz(-1.8176259) q[2];
sx q[2];
rz(-2.18404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0943162) q[1];
sx q[1];
rz(-2.2502406) q[1];
sx q[1];
rz(-1.1761155) q[1];
rz(-1.1464305) q[3];
sx q[3];
rz(-2.6570435) q[3];
sx q[3];
rz(-2.4534006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1404184) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(3.0155638) q[2];
rz(-1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(0.70820156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.66013181) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(-1.460176) q[0];
rz(1.2449645) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(-1.9326899) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7891114) q[0];
sx q[0];
rz(-1.4915823) q[0];
sx q[0];
rz(-1.3080025) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65931321) q[2];
sx q[2];
rz(-2.1954143) q[2];
sx q[2];
rz(-2.5059932) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3371256) q[1];
sx q[1];
rz(-2.2837688) q[1];
sx q[1];
rz(-1.1285524) q[1];
rz(-pi) q[2];
rz(0.94001694) q[3];
sx q[3];
rz(-1.7490083) q[3];
sx q[3];
rz(-0.34561397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37844354) q[2];
sx q[2];
rz(-1.2377137) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(0.59213263) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(-3.106451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.2329344) q[0];
sx q[0];
rz(-0.52353752) q[0];
sx q[0];
rz(-1.3611025) q[0];
rz(-1.9305485) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(0.39168721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91613382) q[0];
sx q[0];
rz(-0.43457169) q[0];
sx q[0];
rz(-0.01390121) q[0];
rz(0.39536706) q[2];
sx q[2];
rz(-2.6569416) q[2];
sx q[2];
rz(-1.5026827) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3561331) q[1];
sx q[1];
rz(-1.9003938) q[1];
sx q[1];
rz(-0.81570028) q[1];
rz(-pi) q[2];
rz(0.52569315) q[3];
sx q[3];
rz(-1.634053) q[3];
sx q[3];
rz(2.4478108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52035511) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(1.3367782) q[2];
rz(-0.76198602) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(-2.3412162) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(15/(14*pi)) q[0];
sx q[0];
rz(-0.30277345) q[0];
sx q[0];
rz(0.57089943) q[0];
rz(-1.7123429) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(-2.9796519) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8853332) q[0];
sx q[0];
rz(-0.67978871) q[0];
sx q[0];
rz(-1.1023561) q[0];
rz(1.5230721) q[2];
sx q[2];
rz(-2.952791) q[2];
sx q[2];
rz(1.4023086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83878126) q[1];
sx q[1];
rz(-1.4421842) q[1];
sx q[1];
rz(1.3978811) q[1];
rz(-pi) q[2];
rz(-1.6041683) q[3];
sx q[3];
rz(-1.0479234) q[3];
sx q[3];
rz(-3.0748622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1841715) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(1.7133678) q[2];
rz(1.9421633) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(-1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.4409055) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.7715001) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(0.11907555) q[2];
sx q[2];
rz(-1.7172114) q[2];
sx q[2];
rz(-1.7211771) q[2];
rz(0.45392848) q[3];
sx q[3];
rz(-0.47149999) q[3];
sx q[3];
rz(-0.75054689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
