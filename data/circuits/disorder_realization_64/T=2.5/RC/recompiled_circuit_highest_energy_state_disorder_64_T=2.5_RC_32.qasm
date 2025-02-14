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
rz(-0.023802726) q[0];
sx q[0];
rz(-2.0355712) q[0];
sx q[0];
rz(0.7769146) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(-0.97631747) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1084749) q[0];
sx q[0];
rz(-0.43370789) q[0];
sx q[0];
rz(-1.8932305) q[0];
rz(-pi) q[1];
x q[1];
rz(1.637865) q[2];
sx q[2];
rz(-2.133495) q[2];
sx q[2];
rz(0.85516155) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7853422) q[1];
sx q[1];
rz(-1.1467198) q[1];
sx q[1];
rz(-2.8431358) q[1];
x q[2];
rz(-2.1438445) q[3];
sx q[3];
rz(-1.5437417) q[3];
sx q[3];
rz(1.388035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6355847) q[2];
sx q[2];
rz(-3.0878461) q[2];
sx q[2];
rz(-2.7347943) q[2];
rz(-2.9721416) q[3];
sx q[3];
rz(-2.6120766) q[3];
sx q[3];
rz(2.0690401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2605543) q[0];
sx q[0];
rz(-2.9091703) q[0];
sx q[0];
rz(-0.0037923092) q[0];
rz(3.0637528) q[1];
sx q[1];
rz(-2.4816315) q[1];
sx q[1];
rz(-0.30581623) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2394442) q[0];
sx q[0];
rz(-1.7551219) q[0];
sx q[0];
rz(2.091616) q[0];
rz(-pi) q[1];
rz(-1.7051093) q[2];
sx q[2];
rz(-2.3181653) q[2];
sx q[2];
rz(0.01625492) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7485436) q[1];
sx q[1];
rz(-2.3094588) q[1];
sx q[1];
rz(-1.0219177) q[1];
rz(-2.3313794) q[3];
sx q[3];
rz(-0.088220291) q[3];
sx q[3];
rz(2.719413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1514312) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(-0.64275098) q[2];
rz(-0.07240545) q[3];
sx q[3];
rz(-1.0591155) q[3];
sx q[3];
rz(-1.36093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0533503) q[0];
sx q[0];
rz(-0.14277661) q[0];
sx q[0];
rz(-0.15637583) q[0];
rz(-0.0414255) q[1];
sx q[1];
rz(-2.5138469) q[1];
sx q[1];
rz(-1.590439) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32937059) q[0];
sx q[0];
rz(-2.5539264) q[0];
sx q[0];
rz(0.37826041) q[0];
rz(-pi) q[1];
rz(-0.088038283) q[2];
sx q[2];
rz(-2.4836342) q[2];
sx q[2];
rz(-0.41706271) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.98387466) q[1];
sx q[1];
rz(-1.0441171) q[1];
sx q[1];
rz(-1.7009363) q[1];
rz(0.66371347) q[3];
sx q[3];
rz(-0.98582375) q[3];
sx q[3];
rz(-3.0970517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40858832) q[2];
sx q[2];
rz(-0.53885794) q[2];
sx q[2];
rz(-3.1206701) q[2];
rz(-0.18713348) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(-3.0278897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9631831) q[0];
sx q[0];
rz(-0.33247501) q[0];
sx q[0];
rz(0.43854976) q[0];
rz(1.5248388) q[1];
sx q[1];
rz(-2.8068145) q[1];
sx q[1];
rz(-2.8964892) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1634451) q[0];
sx q[0];
rz(-0.72448271) q[0];
sx q[0];
rz(1.6969157) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7736362) q[2];
sx q[2];
rz(-1.7178255) q[2];
sx q[2];
rz(-1.8859175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4192701) q[1];
sx q[1];
rz(-2.5825204) q[1];
sx q[1];
rz(1.30617) q[1];
rz(-pi) q[2];
rz(-0.33632261) q[3];
sx q[3];
rz(-0.78236474) q[3];
sx q[3];
rz(-2.4049644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(-2.8098246) q[2];
rz(-2.6541384) q[3];
sx q[3];
rz(-1.0468227) q[3];
sx q[3];
rz(-2.148518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29193923) q[0];
sx q[0];
rz(-1.6878457) q[0];
sx q[0];
rz(-2.3680903) q[0];
rz(-1.1812814) q[1];
sx q[1];
rz(-0.14207323) q[1];
sx q[1];
rz(1.389651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6365125) q[0];
sx q[0];
rz(-0.55547041) q[0];
sx q[0];
rz(-0.73969264) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85948617) q[2];
sx q[2];
rz(-1.0946678) q[2];
sx q[2];
rz(-1.3422333) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4331095) q[1];
sx q[1];
rz(-2.6603087) q[1];
sx q[1];
rz(0.17848707) q[1];
rz(-pi) q[2];
rz(-0.64719154) q[3];
sx q[3];
rz(-0.29509896) q[3];
sx q[3];
rz(1.5403252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2191849) q[2];
sx q[2];
rz(-1.2097404) q[2];
sx q[2];
rz(2.6694471) q[2];
rz(-1.8426497) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(-0.81056547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9206813) q[0];
sx q[0];
rz(-2.8784316) q[0];
sx q[0];
rz(-2.8780908) q[0];
rz(-1.1031411) q[1];
sx q[1];
rz(-1.8270854) q[1];
sx q[1];
rz(0.37364328) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6275245) q[0];
sx q[0];
rz(-1.4980982) q[0];
sx q[0];
rz(-3.0811653) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7079855) q[2];
sx q[2];
rz(-0.73167668) q[2];
sx q[2];
rz(-2.8525994) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8220351) q[1];
sx q[1];
rz(-2.4801835) q[1];
sx q[1];
rz(-2.103785) q[1];
x q[2];
rz(-2.5534036) q[3];
sx q[3];
rz(-1.08687) q[3];
sx q[3];
rz(1.5744792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11334795) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(0.55383468) q[2];
rz(1.3977741) q[3];
sx q[3];
rz(-0.60269409) q[3];
sx q[3];
rz(-0.3127313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
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
rz(2.5572307) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(-1.2297909) q[0];
rz(-2.9025485) q[1];
sx q[1];
rz(-1.5115279) q[1];
sx q[1];
rz(-2.8412433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25324437) q[0];
sx q[0];
rz(-1.408256) q[0];
sx q[0];
rz(-1.4534611) q[0];
rz(-pi) q[1];
rz(-2.5744252) q[2];
sx q[2];
rz(-2.3544899) q[2];
sx q[2];
rz(-1.2445104) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3546037) q[1];
sx q[1];
rz(-1.5704078) q[1];
sx q[1];
rz(1.5860737) q[1];
rz(-2.2836779) q[3];
sx q[3];
rz(-1.5522209) q[3];
sx q[3];
rz(1.1934848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0099237) q[2];
sx q[2];
rz(-1.6841623) q[2];
sx q[2];
rz(-0.24492502) q[2];
rz(-0.51982546) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(2.4533217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35172611) q[0];
sx q[0];
rz(-2.7930197) q[0];
sx q[0];
rz(1.9785731) q[0];
rz(3.0746958) q[1];
sx q[1];
rz(-1.6480548) q[1];
sx q[1];
rz(-2.1287207) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6352309) q[0];
sx q[0];
rz(-1.073045) q[0];
sx q[0];
rz(2.228581) q[0];
x q[1];
rz(-1.8469454) q[2];
sx q[2];
rz(-2.5564402) q[2];
sx q[2];
rz(-0.25979751) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0433181) q[1];
sx q[1];
rz(-1.2407082) q[1];
sx q[1];
rz(1.2374452) q[1];
x q[2];
rz(-0.49533923) q[3];
sx q[3];
rz(-0.74995774) q[3];
sx q[3];
rz(2.0928252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0248727) q[2];
sx q[2];
rz(-2.1634384) q[2];
sx q[2];
rz(-0.32279521) q[2];
rz(2.5358477) q[3];
sx q[3];
rz(-2.3492458) q[3];
sx q[3];
rz(0.34887031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.087273) q[0];
sx q[0];
rz(-2.9948586) q[0];
sx q[0];
rz(0.012454575) q[0];
rz(-2.3948578) q[1];
sx q[1];
rz(-2.2181999) q[1];
sx q[1];
rz(2.8616203) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0664205) q[0];
sx q[0];
rz(-0.1233347) q[0];
sx q[0];
rz(1.9955817) q[0];
rz(1.1367646) q[2];
sx q[2];
rz(-1.869259) q[2];
sx q[2];
rz(3.0500183) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38281967) q[1];
sx q[1];
rz(-2.3291322) q[1];
sx q[1];
rz(2.7121905) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98484184) q[3];
sx q[3];
rz(-2.7142314) q[3];
sx q[3];
rz(2.7387184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81194699) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(0.21198708) q[2];
rz(2.3181465) q[3];
sx q[3];
rz(-1.4444838) q[3];
sx q[3];
rz(-2.8823891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14281808) q[0];
sx q[0];
rz(-3.0840254) q[0];
sx q[0];
rz(2.4488191) q[0];
rz(-2.5686) q[1];
sx q[1];
rz(-1.3447821) q[1];
sx q[1];
rz(-0.43100345) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26769022) q[0];
sx q[0];
rz(-0.55492102) q[0];
sx q[0];
rz(-3.0303427) q[0];
rz(-pi) q[1];
rz(-0.60811483) q[2];
sx q[2];
rz(-0.69497847) q[2];
sx q[2];
rz(-2.8335477) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0242053) q[1];
sx q[1];
rz(-2.3820602) q[1];
sx q[1];
rz(2.7717436) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0447356) q[3];
sx q[3];
rz(-0.60322475) q[3];
sx q[3];
rz(-2.4639377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7196322) q[2];
sx q[2];
rz(-0.27406359) q[2];
sx q[2];
rz(0.56023041) q[2];
rz(-2.6719921) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(0.71389055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8568759) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(-2.3003385) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(0.31568676) q[2];
sx q[2];
rz(-1.2178979) q[2];
sx q[2];
rz(0.60565368) q[2];
rz(-1.5359405) q[3];
sx q[3];
rz(-0.65924725) q[3];
sx q[3];
rz(1.02871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
