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
rz(1.4899878) q[0];
sx q[0];
rz(8.4943354) q[0];
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(-2.0342483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6247834) q[0];
sx q[0];
rz(-2.2342355) q[0];
sx q[0];
rz(1.9883363) q[0];
rz(0.47714699) q[2];
sx q[2];
rz(-2.2331182) q[2];
sx q[2];
rz(1.9805816) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5109374) q[1];
sx q[1];
rz(-1.8488548) q[1];
sx q[1];
rz(-0.11649881) q[1];
x q[2];
rz(2.2060478) q[3];
sx q[3];
rz(-1.2860635) q[3];
sx q[3];
rz(-0.92393827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.063623108) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(1.8135653) q[2];
rz(0.32087457) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(3.0096171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6593453) q[0];
sx q[0];
rz(-3.0292065) q[0];
sx q[0];
rz(-2.2609718) q[0];
rz(-1.847514) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(-0.81726384) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1999232) q[0];
sx q[0];
rz(-2.7669853) q[0];
sx q[0];
rz(1.1767715) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1772637) q[2];
sx q[2];
rz(-2.6173008) q[2];
sx q[2];
rz(2.2146068) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4546928) q[1];
sx q[1];
rz(-2.2001839) q[1];
sx q[1];
rz(-1.8886186) q[1];
rz(-pi) q[2];
rz(-2.0440631) q[3];
sx q[3];
rz(-1.739199) q[3];
sx q[3];
rz(0.67697224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91784224) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(-0.36402738) q[2];
rz(0.98637995) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7746975) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(-0.96631518) q[0];
rz(0.19293383) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(1.6945217) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15196249) q[0];
sx q[0];
rz(-1.3558136) q[0];
sx q[0];
rz(1.5528029) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8615396) q[2];
sx q[2];
rz(-1.4915407) q[2];
sx q[2];
rz(0.62627625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6009532) q[1];
sx q[1];
rz(-0.97266957) q[1];
sx q[1];
rz(2.0289434) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2451285) q[3];
sx q[3];
rz(-1.1655302) q[3];
sx q[3];
rz(2.1851636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43859279) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(2.036371) q[2];
rz(-0.74622074) q[3];
sx q[3];
rz(-1.6522224) q[3];
sx q[3];
rz(-0.97572774) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352017) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(1.6756469) q[0];
rz(2.8566467) q[1];
sx q[1];
rz(-1.0703215) q[1];
sx q[1];
rz(2.7526061) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3272414) q[0];
sx q[0];
rz(-1.2885433) q[0];
sx q[0];
rz(-1.809657) q[0];
x q[1];
rz(1.0842501) q[2];
sx q[2];
rz(-1.3975189) q[2];
sx q[2];
rz(1.372352) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.092408808) q[1];
sx q[1];
rz(-1.0256983) q[1];
sx q[1];
rz(-1.070302) q[1];
x q[2];
rz(0.22104927) q[3];
sx q[3];
rz(-2.3813558) q[3];
sx q[3];
rz(1.5628712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(0.45670613) q[2];
rz(-1.6263973) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.91963768) q[0];
sx q[0];
rz(-2.0018405) q[0];
sx q[0];
rz(0.36002457) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.6123687) q[1];
sx q[1];
rz(-0.49450758) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8202782) q[0];
sx q[0];
rz(-1.2914133) q[0];
sx q[0];
rz(-1.7000291) q[0];
x q[1];
rz(1.7062543) q[2];
sx q[2];
rz(-1.6098621) q[2];
sx q[2];
rz(-2.1087697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33489409) q[1];
sx q[1];
rz(-2.1999199) q[1];
sx q[1];
rz(1.6449528) q[1];
rz(-1.5962283) q[3];
sx q[3];
rz(-1.126822) q[3];
sx q[3];
rz(-2.3412995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4313724) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(2.8430856) q[2];
rz(0.31202894) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1453778) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(-0.19590713) q[0];
rz(0.021082489) q[1];
sx q[1];
rz(-1.7430051) q[1];
sx q[1];
rz(1.9063937) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7346749) q[0];
sx q[0];
rz(-1.9039246) q[0];
sx q[0];
rz(-1.245265) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4371106) q[2];
sx q[2];
rz(-1.3085758) q[2];
sx q[2];
rz(-2.8395677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1393309) q[1];
sx q[1];
rz(-1.1552703) q[1];
sx q[1];
rz(-0.63654391) q[1];
x q[2];
rz(-1.487021) q[3];
sx q[3];
rz(-1.8063074) q[3];
sx q[3];
rz(-1.868353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6482676) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(-2.7098999) q[2];
rz(1.4124983) q[3];
sx q[3];
rz(-2.4192211) q[3];
sx q[3];
rz(-3.0055962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39111185) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(3.0294763) q[0];
rz(0.21513367) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(1.1134061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3413234) q[0];
sx q[0];
rz(-1.1870664) q[0];
sx q[0];
rz(1.2905489) q[0];
rz(-pi) q[1];
rz(-2.8261975) q[2];
sx q[2];
rz(-1.8176259) q[2];
sx q[2];
rz(-2.18404) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4091195) q[1];
sx q[1];
rz(-1.8745683) q[1];
sx q[1];
rz(2.4227546) q[1];
rz(1.1464305) q[3];
sx q[3];
rz(-0.48454912) q[3];
sx q[3];
rz(0.68819203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.001174288) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(-0.12602885) q[2];
rz(2.0942988) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(-2.4333911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66013181) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(1.2449645) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(1.2089027) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35248127) q[0];
sx q[0];
rz(-1.6500104) q[0];
sx q[0];
rz(-1.8335901) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4822794) q[2];
sx q[2];
rz(-2.1954143) q[2];
sx q[2];
rz(0.63559947) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3371256) q[1];
sx q[1];
rz(-0.85782385) q[1];
sx q[1];
rz(-2.0130403) q[1];
x q[2];
rz(-1.8672089) q[3];
sx q[3];
rz(-2.489438) q[3];
sx q[3];
rz(-1.6782827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7631491) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(-2.54946) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(3.106451) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9086583) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(1.7804902) q[0];
rz(1.9305485) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(-0.39168721) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2254588) q[0];
sx q[0];
rz(-0.43457169) q[0];
sx q[0];
rz(0.01390121) q[0];
rz(-2.6892745) q[2];
sx q[2];
rz(-1.3903793) q[2];
sx q[2];
rz(-0.42186055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0605676) q[1];
sx q[1];
rz(-0.86522663) q[1];
sx q[1];
rz(-2.702436) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4976981) q[3];
sx q[3];
rz(-2.0953296) q[3];
sx q[3];
rz(0.91367164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6212375) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(-1.8048145) q[2];
rz(-0.76198602) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(15/(14*pi)) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(-0.57089943) q[0];
rz(1.4292498) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(-0.16194078) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83308342) q[0];
sx q[0];
rz(-0.97531318) q[0];
sx q[0];
rz(-0.34992976) q[0];
rz(1.5230721) q[2];
sx q[2];
rz(-2.952791) q[2];
sx q[2];
rz(-1.7392841) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0429749) q[1];
sx q[1];
rz(-0.21511714) q[1];
sx q[1];
rz(2.2153562) q[1];
rz(-pi) q[2];
rz(-1.6041683) q[3];
sx q[3];
rz(-1.0479234) q[3];
sx q[3];
rz(-3.0748622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.95742115) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(-1.7133678) q[2];
rz(1.9421633) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(1.7709581) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4409055) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.3700925) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(1.7182405) q[2];
sx q[2];
rz(-1.4530008) q[2];
sx q[2];
rz(3.0086649) q[2];
rz(2.6876642) q[3];
sx q[3];
rz(-2.6700927) q[3];
sx q[3];
rz(2.3910458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
