OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(-0.49180254) q[0];
sx q[0];
rz(-2.9536182) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(-1.6245276) q[1];
sx q[1];
rz(-2.7741073) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6715235) q[0];
sx q[0];
rz(-1.0363665) q[0];
sx q[0];
rz(-0.1202017) q[0];
rz(1.8194524) q[2];
sx q[2];
rz(-2.6373632) q[2];
sx q[2];
rz(2.3766975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0109947) q[1];
sx q[1];
rz(-1.3291318) q[1];
sx q[1];
rz(2.9690837) q[1];
rz(-0.46842694) q[3];
sx q[3];
rz(-2.763063) q[3];
sx q[3];
rz(0.98640501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1774896) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(2.5906079) q[2];
rz(-1.3059113) q[3];
sx q[3];
rz(-1.4923613) q[3];
sx q[3];
rz(-1.3163542) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630163) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(-0.4719032) q[0];
rz(-0.42981237) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(-2.205251) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443003) q[0];
sx q[0];
rz(-2.8951277) q[0];
sx q[0];
rz(0.36578567) q[0];
rz(-pi) q[1];
rz(1.9905375) q[2];
sx q[2];
rz(-0.83101666) q[2];
sx q[2];
rz(1.3295528) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71948235) q[1];
sx q[1];
rz(-1.6750095) q[1];
sx q[1];
rz(2.4476493) q[1];
x q[2];
rz(1.068088) q[3];
sx q[3];
rz(-1.5912676) q[3];
sx q[3];
rz(2.3308144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3669746) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(0.42638391) q[2];
rz(1.9042227) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24580978) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(-2.202503) q[0];
rz(-0.89871961) q[1];
sx q[1];
rz(-2.6627314) q[1];
sx q[1];
rz(-2.5476707) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.827841) q[0];
sx q[0];
rz(-1.8694287) q[0];
sx q[0];
rz(1.9065501) q[0];
rz(-0.8823231) q[2];
sx q[2];
rz(-1.5261298) q[2];
sx q[2];
rz(2.4633212) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9582639) q[1];
sx q[1];
rz(-1.7839583) q[1];
sx q[1];
rz(-2.0002685) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5251758) q[3];
sx q[3];
rz(-2.5865002) q[3];
sx q[3];
rz(0.68716955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64017355) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(-1.7017986) q[2];
rz(0.38763186) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(-0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3751635) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(-2.638812) q[0];
rz(2.373383) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(0.75685135) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216953) q[0];
sx q[0];
rz(-1.9031525) q[0];
sx q[0];
rz(-0.38145782) q[0];
rz(-pi) q[1];
x q[1];
rz(0.012979522) q[2];
sx q[2];
rz(-2.0364967) q[2];
sx q[2];
rz(-0.73060689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5731569) q[1];
sx q[1];
rz(-0.70306289) q[1];
sx q[1];
rz(0.99848024) q[1];
rz(-pi) q[2];
rz(0.73369153) q[3];
sx q[3];
rz(-1.1847704) q[3];
sx q[3];
rz(1.2804077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7148774) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(-1.4871917) q[2];
rz(2.5590844) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(-2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(2.3838682) q[0];
rz(1.2879397) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(1.0505189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61705631) q[0];
sx q[0];
rz(-2.7334088) q[0];
sx q[0];
rz(0.88390669) q[0];
x q[1];
rz(-2.1804785) q[2];
sx q[2];
rz(-2.4498307) q[2];
sx q[2];
rz(-1.320653) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.037462385) q[1];
sx q[1];
rz(-0.81172746) q[1];
sx q[1];
rz(-2.3292259) q[1];
x q[2];
rz(2.1220783) q[3];
sx q[3];
rz(-0.84421221) q[3];
sx q[3];
rz(2.7725078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.918255) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(0.63344947) q[2];
rz(1.9472306) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(-0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.441992) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(-0.25318405) q[0];
rz(-1.6075915) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(-1.4621428) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1519449) q[0];
sx q[0];
rz(-1.4914574) q[0];
sx q[0];
rz(0.21110714) q[0];
rz(-0.92163779) q[2];
sx q[2];
rz(-1.3772794) q[2];
sx q[2];
rz(0.18093872) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51965442) q[1];
sx q[1];
rz(-2.7192273) q[1];
sx q[1];
rz(-0.42610355) q[1];
rz(0.31520321) q[3];
sx q[3];
rz(-1.8990371) q[3];
sx q[3];
rz(-0.70716508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9399461) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(-0.80491006) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(-1.0837519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.2729623) q[0];
sx q[0];
rz(-2.0721764) q[0];
sx q[0];
rz(0.92765635) q[0];
rz(-2.1169128) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(-2.129508) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9575858) q[0];
sx q[0];
rz(-2.036096) q[0];
sx q[0];
rz(2.1561949) q[0];
x q[1];
rz(1.7094678) q[2];
sx q[2];
rz(-1.1853293) q[2];
sx q[2];
rz(-0.0027545714) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.95841366) q[1];
sx q[1];
rz(-1.603754) q[1];
sx q[1];
rz(-0.54380137) q[1];
x q[2];
rz(1.4885694) q[3];
sx q[3];
rz(-2.1279018) q[3];
sx q[3];
rz(2.0092056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53081375) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5291418) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(0.21417831) q[0];
rz(1.051349) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(0.28373757) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8003214) q[0];
sx q[0];
rz(-2.8386136) q[0];
sx q[0];
rz(-3.0269701) q[0];
rz(-1.0546513) q[2];
sx q[2];
rz(-2.7661341) q[2];
sx q[2];
rz(-1.6106538) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2092065) q[1];
sx q[1];
rz(-1.7665518) q[1];
sx q[1];
rz(-0.50435658) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65225668) q[3];
sx q[3];
rz(-1.0271003) q[3];
sx q[3];
rz(-3.101845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.45067898) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(1.696375) q[2];
rz(1.5971659) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(-0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504869) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(3.0723363) q[0];
rz(1.6537369) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(1.5690631) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42800316) q[0];
sx q[0];
rz(-2.028095) q[0];
sx q[0];
rz(-2.2767115) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4885159) q[2];
sx q[2];
rz(-1.2459718) q[2];
sx q[2];
rz(2.5282853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8005341) q[1];
sx q[1];
rz(-1.3375999) q[1];
sx q[1];
rz(1.484616) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8435555) q[3];
sx q[3];
rz(-2.4747304) q[3];
sx q[3];
rz(0.71344261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1853603) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(-2.774003) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(2.1504413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-0.39127025) q[0];
sx q[0];
rz(2.4998253) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.1522013) q[1];
sx q[1];
rz(0.26783255) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5689241) q[0];
sx q[0];
rz(-0.84416443) q[0];
sx q[0];
rz(-2.5323244) q[0];
rz(-2.1688813) q[2];
sx q[2];
rz(-0.22062606) q[2];
sx q[2];
rz(2.0711183) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6707582) q[1];
sx q[1];
rz(-1.2467614) q[1];
sx q[1];
rz(-1.2653989) q[1];
rz(-pi) q[2];
rz(2.855741) q[3];
sx q[3];
rz(-2.0945858) q[3];
sx q[3];
rz(-1.7459735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3315167) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(1.4769185) q[2];
rz(-0.26633513) q[3];
sx q[3];
rz(-2.895152) q[3];
sx q[3];
rz(0.5464856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(-0.77990445) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(-0.65200381) q[2];
sx q[2];
rz(-2.3163788) q[2];
sx q[2];
rz(-0.083995081) q[2];
rz(-0.6775425) q[3];
sx q[3];
rz(-1.7384221) q[3];
sx q[3];
rz(1.3330028) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
