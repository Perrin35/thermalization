OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.36800185) q[0];
sx q[0];
rz(-0.79080963) q[0];
sx q[0];
rz(0.33413449) q[0];
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(-1.2184719) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96555644) q[0];
sx q[0];
rz(-0.82057014) q[0];
sx q[0];
rz(1.5216989) q[0];
rz(2.676744) q[2];
sx q[2];
rz(-1.6072304) q[2];
sx q[2];
rz(0.51793098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6005046) q[1];
sx q[1];
rz(-2.0588015) q[1];
sx q[1];
rz(1.0469251) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2741332) q[3];
sx q[3];
rz(-0.28775035) q[3];
sx q[3];
rz(-0.61144637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5228287) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(2.5640326) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(-2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98786551) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(-2.7541449) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.739025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57099045) q[0];
sx q[0];
rz(-1.5413069) q[0];
sx q[0];
rz(1.5748613) q[0];
x q[1];
rz(-1.7895133) q[2];
sx q[2];
rz(-1.0220851) q[2];
sx q[2];
rz(-2.1825841) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1075588) q[1];
sx q[1];
rz(-0.65980136) q[1];
sx q[1];
rz(2.4504513) q[1];
x q[2];
rz(-0.12947793) q[3];
sx q[3];
rz(-2.5054512) q[3];
sx q[3];
rz(-2.5667218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7188321) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(2.823901) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(-2.3247705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(-1.4136219) q[0];
rz(-0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(0.40107045) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.385752) q[0];
sx q[0];
rz(-2.7999561) q[0];
sx q[0];
rz(1.2582448) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6443278) q[2];
sx q[2];
rz(-1.7228848) q[2];
sx q[2];
rz(2.6459141) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2420826) q[1];
sx q[1];
rz(-2.3009355) q[1];
sx q[1];
rz(1.5793369) q[1];
rz(1.1776351) q[3];
sx q[3];
rz(-1.3893681) q[3];
sx q[3];
rz(-1.3822615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59427375) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(0.55580124) q[2];
rz(-2.1650971) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(1.697631) q[0];
rz(1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(-2.8881883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2052106) q[0];
sx q[0];
rz(-2.10996) q[0];
sx q[0];
rz(-3.0217231) q[0];
rz(-pi) q[1];
rz(0.01732145) q[2];
sx q[2];
rz(-1.0850731) q[2];
sx q[2];
rz(2.3522365) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6654012) q[1];
sx q[1];
rz(-2.5051077) q[1];
sx q[1];
rz(-0.16237662) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13417379) q[3];
sx q[3];
rz(-2.4878256) q[3];
sx q[3];
rz(1.8133481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(-0.78732642) q[2];
rz(0.91286719) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(-2.1319938) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(0.062967904) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(0.45809349) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9034018) q[0];
sx q[0];
rz(-1.5731249) q[0];
sx q[0];
rz(2.6989614) q[0];
x q[1];
rz(0.45256726) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(0.82816154) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.304368) q[1];
sx q[1];
rz(-0.75290426) q[1];
sx q[1];
rz(-1.3992845) q[1];
x q[2];
rz(3.0466988) q[3];
sx q[3];
rz(-2.0695544) q[3];
sx q[3];
rz(-0.14548485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1725585) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(-1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(-0.91947412) q[0];
rz(3.0793076) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(1.2671635) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125268) q[0];
sx q[0];
rz(-0.8493087) q[0];
sx q[0];
rz(-0.52961403) q[0];
x q[1];
rz(0.15684776) q[2];
sx q[2];
rz(-2.3914797) q[2];
sx q[2];
rz(-1.4175121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25176469) q[1];
sx q[1];
rz(-1.3740731) q[1];
sx q[1];
rz(-1.9038492) q[1];
rz(-pi) q[2];
x q[2];
rz(1.091286) q[3];
sx q[3];
rz(-1.7877842) q[3];
sx q[3];
rz(0.12245164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.9617812) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(0.58376694) q[2];
rz(2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07847438) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(2.069058) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(1.1118719) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3277153) q[0];
sx q[0];
rz(-1.0881256) q[0];
sx q[0];
rz(-3.0359603) q[0];
rz(-pi) q[1];
rz(-0.56356168) q[2];
sx q[2];
rz(-0.60699082) q[2];
sx q[2];
rz(2.3129472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.48549451) q[1];
sx q[1];
rz(-1.494325) q[1];
sx q[1];
rz(-0.2904201) q[1];
x q[2];
rz(0.92915793) q[3];
sx q[3];
rz(-2.540179) q[3];
sx q[3];
rz(-1.8638368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(-1.973935) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(0.73079601) q[0];
rz(0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(0.75497595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46144339) q[0];
sx q[0];
rz(-1.7168683) q[0];
sx q[0];
rz(2.2938674) q[0];
rz(-pi) q[1];
rz(2.1515498) q[2];
sx q[2];
rz(-2.2120737) q[2];
sx q[2];
rz(0.61604797) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6983812) q[1];
sx q[1];
rz(-1.6834007) q[1];
sx q[1];
rz(-2.5618782) q[1];
rz(-0.49236492) q[3];
sx q[3];
rz(-1.0400606) q[3];
sx q[3];
rz(0.20519557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3770611) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(0.72717845) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(2.0027347) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(0.019502217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560023) q[0];
sx q[0];
rz(-1.7770355) q[0];
sx q[0];
rz(3.0598559) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7712014) q[2];
sx q[2];
rz(-1.4223756) q[2];
sx q[2];
rz(-0.60147775) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1972678) q[1];
sx q[1];
rz(-1.7515344) q[1];
sx q[1];
rz(-0.9428057) q[1];
rz(-pi) q[2];
rz(-0.1653413) q[3];
sx q[3];
rz(-1.2623566) q[3];
sx q[3];
rz(0.43499085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3432817) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(-0.48450255) q[0];
rz(1.7548521) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(1.1482931) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0120221) q[0];
sx q[0];
rz(-0.14229933) q[0];
sx q[0];
rz(-2.9905031) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5103288) q[2];
sx q[2];
rz(-1.3901276) q[2];
sx q[2];
rz(2.0618912) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53459586) q[1];
sx q[1];
rz(-1.7417223) q[1];
sx q[1];
rz(-2.2439438) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6230574) q[3];
sx q[3];
rz(-1.9196379) q[3];
sx q[3];
rz(-0.40961743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(-0.45583367) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1098332) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-0.96314349) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(-0.49114901) q[2];
sx q[2];
rz(-0.72577234) q[2];
sx q[2];
rz(2.5189248) q[2];
rz(-2.1309489) q[3];
sx q[3];
rz(-1.5394566) q[3];
sx q[3];
rz(1.4395366) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
