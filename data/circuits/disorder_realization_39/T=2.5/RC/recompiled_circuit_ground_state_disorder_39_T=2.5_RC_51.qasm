OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7118536) q[0];
sx q[0];
rz(5.138152) q[0];
sx q[0];
rz(10.44147) q[0];
rz(-0.84438762) q[1];
sx q[1];
rz(-0.8165741) q[1];
sx q[1];
rz(-2.5752697) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53439683) q[0];
sx q[0];
rz(-1.6677987) q[0];
sx q[0];
rz(0.96235458) q[0];
x q[1];
rz(3.064134) q[2];
sx q[2];
rz(-1.3734387) q[2];
sx q[2];
rz(-3.0826621) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8992246) q[1];
sx q[1];
rz(-1.4903233) q[1];
sx q[1];
rz(-2.7182275) q[1];
rz(0.30632354) q[3];
sx q[3];
rz(-0.60738436) q[3];
sx q[3];
rz(-2.4709783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5113968) q[2];
sx q[2];
rz(-2.2679195) q[2];
sx q[2];
rz(-2.013618) q[2];
rz(1.5026211) q[3];
sx q[3];
rz(-2.2962544) q[3];
sx q[3];
rz(0.42475167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.86654919) q[0];
sx q[0];
rz(-2.6847222) q[0];
sx q[0];
rz(0.04091111) q[0];
rz(2.5919137) q[1];
sx q[1];
rz(-2.7626541) q[1];
sx q[1];
rz(3.0612225) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.550404) q[0];
sx q[0];
rz(-0.72232095) q[0];
sx q[0];
rz(-3.1352311) q[0];
rz(-pi) q[1];
rz(2.2533312) q[2];
sx q[2];
rz(-0.76768926) q[2];
sx q[2];
rz(-1.7861799) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.077420635) q[1];
sx q[1];
rz(-1.6257833) q[1];
sx q[1];
rz(-1.085678) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96719043) q[3];
sx q[3];
rz(-1.3466102) q[3];
sx q[3];
rz(-2.8160915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72767672) q[2];
sx q[2];
rz(-0.84299403) q[2];
sx q[2];
rz(-0.28398871) q[2];
rz(1.6207638) q[3];
sx q[3];
rz(-1.5903383) q[3];
sx q[3];
rz(-0.76333299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82336998) q[0];
sx q[0];
rz(-0.1544054) q[0];
sx q[0];
rz(-2.7247317) q[0];
rz(-0.93337026) q[1];
sx q[1];
rz(-0.88637543) q[1];
sx q[1];
rz(-2.0534168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8899208) q[0];
sx q[0];
rz(-1.8150738) q[0];
sx q[0];
rz(-2.1273462) q[0];
rz(1.5303262) q[2];
sx q[2];
rz(-2.1310316) q[2];
sx q[2];
rz(0.76228415) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8510906) q[1];
sx q[1];
rz(-0.41176418) q[1];
sx q[1];
rz(2.3985833) q[1];
x q[2];
rz(1.635664) q[3];
sx q[3];
rz(-1.8577575) q[3];
sx q[3];
rz(-2.502041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0205959) q[2];
sx q[2];
rz(-2.4211297) q[2];
sx q[2];
rz(-0.32568112) q[2];
rz(2.7459512) q[3];
sx q[3];
rz(-2.6511657) q[3];
sx q[3];
rz(2.4237848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91738236) q[0];
sx q[0];
rz(-0.93467394) q[0];
sx q[0];
rz(3.000946) q[0];
rz(2.7463101) q[1];
sx q[1];
rz(-1.544516) q[1];
sx q[1];
rz(-2.7093754) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0790063) q[0];
sx q[0];
rz(-2.5012928) q[0];
sx q[0];
rz(-0.85377924) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0555058) q[2];
sx q[2];
rz(-2.0755092) q[2];
sx q[2];
rz(-1.9534257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85104698) q[1];
sx q[1];
rz(-1.6509027) q[1];
sx q[1];
rz(0.27388957) q[1];
x q[2];
rz(-2.5514546) q[3];
sx q[3];
rz(-1.034684) q[3];
sx q[3];
rz(2.9394583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29493368) q[2];
sx q[2];
rz(-2.0294971) q[2];
sx q[2];
rz(-0.46889949) q[2];
rz(-2.9851959) q[3];
sx q[3];
rz(-2.1726435) q[3];
sx q[3];
rz(-1.2676574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0676607) q[0];
sx q[0];
rz(-2.0079375) q[0];
sx q[0];
rz(2.3611948) q[0];
rz(2.228915) q[1];
sx q[1];
rz(-2.3527805) q[1];
sx q[1];
rz(1.7524293) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75994223) q[0];
sx q[0];
rz(-1.8382605) q[0];
sx q[0];
rz(2.777369) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19540968) q[2];
sx q[2];
rz(-0.59048803) q[2];
sx q[2];
rz(2.5974642) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4215736) q[1];
sx q[1];
rz(-1.6813845) q[1];
sx q[1];
rz(-2.5441267) q[1];
rz(-pi) q[2];
rz(0.88646226) q[3];
sx q[3];
rz(-1.9406978) q[3];
sx q[3];
rz(-2.1641638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3461561) q[2];
sx q[2];
rz(-2.4861591) q[2];
sx q[2];
rz(-0.1795086) q[2];
rz(-1.6620212) q[3];
sx q[3];
rz(-0.69412762) q[3];
sx q[3];
rz(3.1365385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805304) q[0];
sx q[0];
rz(-1.5164277) q[0];
sx q[0];
rz(-1.6777212) q[0];
rz(0.013710984) q[1];
sx q[1];
rz(-1.4669712) q[1];
sx q[1];
rz(-0.94508583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6701407) q[0];
sx q[0];
rz(-1.2218804) q[0];
sx q[0];
rz(-2.3405911) q[0];
rz(1.7795947) q[2];
sx q[2];
rz(-1.8755302) q[2];
sx q[2];
rz(-1.3688251) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42691222) q[1];
sx q[1];
rz(-1.1537997) q[1];
sx q[1];
rz(-0.51583536) q[1];
x q[2];
rz(2.6078647) q[3];
sx q[3];
rz(-2.0954359) q[3];
sx q[3];
rz(2.644705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1549687) q[2];
sx q[2];
rz(-1.0633609) q[2];
sx q[2];
rz(-2.092579) q[2];
rz(-1.6074041) q[3];
sx q[3];
rz(-1.0020071) q[3];
sx q[3];
rz(1.3656507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.5904215) q[0];
sx q[0];
rz(-0.48749247) q[0];
sx q[0];
rz(2.3792939) q[0];
rz(2.6679692) q[1];
sx q[1];
rz(-1.7601687) q[1];
sx q[1];
rz(-0.90829888) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98971924) q[0];
sx q[0];
rz(-0.79209581) q[0];
sx q[0];
rz(3.1342952) q[0];
rz(-pi) q[1];
rz(1.7920831) q[2];
sx q[2];
rz(-1.7402116) q[2];
sx q[2];
rz(-2.7986023) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6544853) q[1];
sx q[1];
rz(-1.5853549) q[1];
sx q[1];
rz(1.8076987) q[1];
rz(-3.129035) q[3];
sx q[3];
rz(-1.5987442) q[3];
sx q[3];
rz(1.8699153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.121754) q[2];
sx q[2];
rz(-0.25889954) q[2];
sx q[2];
rz(2.1799977) q[2];
rz(-0.52729765) q[3];
sx q[3];
rz(-1.3798102) q[3];
sx q[3];
rz(0.56078792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.90802646) q[0];
sx q[0];
rz(-0.59995025) q[0];
sx q[0];
rz(2.7954234) q[0];
rz(1.8442122) q[1];
sx q[1];
rz(-0.66645122) q[1];
sx q[1];
rz(-2.5328439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.611871) q[0];
sx q[0];
rz(-2.4606315) q[0];
sx q[0];
rz(-1.7044071) q[0];
rz(-pi) q[1];
rz(0.66242156) q[2];
sx q[2];
rz(-2.8067744) q[2];
sx q[2];
rz(-1.7551369) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1947555) q[1];
sx q[1];
rz(-1.5018592) q[1];
sx q[1];
rz(-0.26534715) q[1];
rz(-pi) q[2];
rz(-2.1431646) q[3];
sx q[3];
rz(-0.86248904) q[3];
sx q[3];
rz(-0.64841333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79494563) q[2];
sx q[2];
rz(-2.5865159) q[2];
sx q[2];
rz(-0.65183276) q[2];
rz(-0.73976222) q[3];
sx q[3];
rz(-2.0756105) q[3];
sx q[3];
rz(-0.8968001) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5935434) q[0];
sx q[0];
rz(-0.47320047) q[0];
sx q[0];
rz(2.5439673) q[0];
rz(-1.8400486) q[1];
sx q[1];
rz(-1.5751244) q[1];
sx q[1];
rz(0.32599932) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2151011) q[0];
sx q[0];
rz(-0.85868663) q[0];
sx q[0];
rz(2.7551921) q[0];
rz(0.65970274) q[2];
sx q[2];
rz(-2.7355425) q[2];
sx q[2];
rz(-0.34161196) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.9218775) q[1];
sx q[1];
rz(-1.3120717) q[1];
sx q[1];
rz(0.90176505) q[1];
rz(0.80504412) q[3];
sx q[3];
rz(-2.9784103) q[3];
sx q[3];
rz(-1.9022579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7035383) q[2];
sx q[2];
rz(-1.2039098) q[2];
sx q[2];
rz(0.34029141) q[2];
rz(-1.641364) q[3];
sx q[3];
rz(-2.5989125) q[3];
sx q[3];
rz(2.5435508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.930645) q[0];
sx q[0];
rz(-1.1177381) q[0];
sx q[0];
rz(-0.043638226) q[0];
rz(-0.71815193) q[1];
sx q[1];
rz(-1.3190045) q[1];
sx q[1];
rz(-1.4290379) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6788504) q[0];
sx q[0];
rz(-1.7587816) q[0];
sx q[0];
rz(-1.3975271) q[0];
rz(-pi) q[1];
rz(-0.67769717) q[2];
sx q[2];
rz(-1.4645828) q[2];
sx q[2];
rz(2.8863557) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2344115) q[1];
sx q[1];
rz(-2.6446807) q[1];
sx q[1];
rz(2.9610544) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7765482) q[3];
sx q[3];
rz(-1.2170346) q[3];
sx q[3];
rz(2.7937074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67937294) q[2];
sx q[2];
rz(-2.2023109) q[2];
sx q[2];
rz(2.3636554) q[2];
rz(2.098162) q[3];
sx q[3];
rz(-1.611004) q[3];
sx q[3];
rz(1.9061609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048024561) q[0];
sx q[0];
rz(-1.0014191) q[0];
sx q[0];
rz(0.77145664) q[0];
rz(-1.7780766) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(-0.72259283) q[2];
sx q[2];
rz(-1.1579787) q[2];
sx q[2];
rz(-2.9415393) q[2];
rz(-2.6778775) q[3];
sx q[3];
rz(-0.9526997) q[3];
sx q[3];
rz(0.37982527) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
