OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8371589) q[0];
sx q[0];
rz(4.9871939) q[0];
sx q[0];
rz(11.468588) q[0];
rz(-0.9552362) q[1];
sx q[1];
rz(7.0248338) q[1];
sx q[1];
rz(6.4344814) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3442952) q[0];
sx q[0];
rz(-1.5452641) q[0];
sx q[0];
rz(1.6961369) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7593616) q[2];
sx q[2];
rz(-2.5841568) q[2];
sx q[2];
rz(0.093890015) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9469493) q[1];
sx q[1];
rz(-2.2990755) q[1];
sx q[1];
rz(-0.76232736) q[1];
rz(1.2829078) q[3];
sx q[3];
rz(-0.79018738) q[3];
sx q[3];
rz(-1.1367281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0155045) q[2];
sx q[2];
rz(-1.2903004) q[2];
sx q[2];
rz(-0.31910953) q[2];
rz(1.1110405) q[3];
sx q[3];
rz(-0.55912656) q[3];
sx q[3];
rz(1.8266953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1183209) q[0];
sx q[0];
rz(-0.51902223) q[0];
sx q[0];
rz(-0.58854377) q[0];
rz(2.5449246) q[1];
sx q[1];
rz(-1.3304973) q[1];
sx q[1];
rz(-2.9002424) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5198181) q[0];
sx q[0];
rz(-2.1960917) q[0];
sx q[0];
rz(0.53859512) q[0];
rz(1.0701837) q[2];
sx q[2];
rz(-1.6334849) q[2];
sx q[2];
rz(2.3586065) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0343098) q[1];
sx q[1];
rz(-1.4713381) q[1];
sx q[1];
rz(-2.0080272) q[1];
rz(-1.4292795) q[3];
sx q[3];
rz(-0.63014537) q[3];
sx q[3];
rz(-0.41434789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1179463) q[2];
sx q[2];
rz(-1.6829374) q[2];
sx q[2];
rz(-0.092197593) q[2];
rz(-0.89208952) q[3];
sx q[3];
rz(-0.8725608) q[3];
sx q[3];
rz(1.0649072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43111619) q[0];
sx q[0];
rz(-1.6350063) q[0];
sx q[0];
rz(-3.0431252) q[0];
rz(0.59421986) q[1];
sx q[1];
rz(-1.0585982) q[1];
sx q[1];
rz(2.1509511) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7567609) q[0];
sx q[0];
rz(-0.38479003) q[0];
sx q[0];
rz(1.8033474) q[0];
rz(2.5712396) q[2];
sx q[2];
rz(-2.2945171) q[2];
sx q[2];
rz(-0.090426771) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2222683) q[1];
sx q[1];
rz(-1.3263371) q[1];
sx q[1];
rz(0.9217086) q[1];
x q[2];
rz(-2.9143067) q[3];
sx q[3];
rz(-1.1133988) q[3];
sx q[3];
rz(1.5279087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6777665) q[2];
sx q[2];
rz(-0.78654424) q[2];
sx q[2];
rz(-2.3393935) q[2];
rz(-1.5471316) q[3];
sx q[3];
rz(-1.0651257) q[3];
sx q[3];
rz(-2.2772363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7441854) q[0];
sx q[0];
rz(-0.35905251) q[0];
sx q[0];
rz(-1.8205951) q[0];
rz(-1.6162704) q[1];
sx q[1];
rz(-2.1883712) q[1];
sx q[1];
rz(0.1604518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1421776) q[0];
sx q[0];
rz(-0.83926187) q[0];
sx q[0];
rz(-2.8150303) q[0];
rz(2.552659) q[2];
sx q[2];
rz(-2.9351165) q[2];
sx q[2];
rz(-1.4331417) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2962118) q[1];
sx q[1];
rz(-1.7239162) q[1];
sx q[1];
rz(0.91836849) q[1];
rz(2.6053107) q[3];
sx q[3];
rz(-1.9347526) q[3];
sx q[3];
rz(-2.3126147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82115951) q[2];
sx q[2];
rz(-1.6511788) q[2];
sx q[2];
rz(-0.56905812) q[2];
rz(1.9197561) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(-3.0088185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4959167) q[0];
sx q[0];
rz(-0.78302947) q[0];
sx q[0];
rz(2.3107279) q[0];
rz(1.9310541) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(-2.1499706) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.338991) q[0];
sx q[0];
rz(-2.3590238) q[0];
sx q[0];
rz(2.0156142) q[0];
rz(2.2575602) q[2];
sx q[2];
rz(-1.3515389) q[2];
sx q[2];
rz(2.925556) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2048671) q[1];
sx q[1];
rz(-0.74581742) q[1];
sx q[1];
rz(1.3406483) q[1];
rz(-2.6001106) q[3];
sx q[3];
rz(-2.1953744) q[3];
sx q[3];
rz(-1.2219714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7603989) q[2];
sx q[2];
rz(-0.92279592) q[2];
sx q[2];
rz(2.7823616) q[2];
rz(0.71581101) q[3];
sx q[3];
rz(-0.78407136) q[3];
sx q[3];
rz(-1.951096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0773709) q[0];
sx q[0];
rz(-1.8697898) q[0];
sx q[0];
rz(0.52128681) q[0];
rz(-0.84292665) q[1];
sx q[1];
rz(-1.1784252) q[1];
sx q[1];
rz(0.11046031) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0960707) q[0];
sx q[0];
rz(-1.318249) q[0];
sx q[0];
rz(0.70742328) q[0];
rz(-pi) q[1];
rz(-0.92824061) q[2];
sx q[2];
rz(-1.5945598) q[2];
sx q[2];
rz(1.0500963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0065436) q[1];
sx q[1];
rz(-1.1056756) q[1];
sx q[1];
rz(-2.9198482) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.13917) q[3];
sx q[3];
rz(-0.56474287) q[3];
sx q[3];
rz(-0.016591681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20737401) q[2];
sx q[2];
rz(-1.5418345) q[2];
sx q[2];
rz(0.99384394) q[2];
rz(0.55073109) q[3];
sx q[3];
rz(-2.6271074) q[3];
sx q[3];
rz(-1.6186835) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2391424) q[0];
sx q[0];
rz(-0.37343326) q[0];
sx q[0];
rz(2.9504839) q[0];
rz(0.36901078) q[1];
sx q[1];
rz(-1.399682) q[1];
sx q[1];
rz(2.5083127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8275237) q[0];
sx q[0];
rz(-1.6681328) q[0];
sx q[0];
rz(2.901696) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0379637) q[2];
sx q[2];
rz(-2.2351801) q[2];
sx q[2];
rz(-0.54069041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.889588) q[1];
sx q[1];
rz(-1.7637296) q[1];
sx q[1];
rz(-1.7798406) q[1];
rz(0.19509372) q[3];
sx q[3];
rz(-2.4332402) q[3];
sx q[3];
rz(-0.1699902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9621027) q[2];
sx q[2];
rz(-1.7682163) q[2];
sx q[2];
rz(2.7607259) q[2];
rz(-1.3880091) q[3];
sx q[3];
rz(-2.1151147) q[3];
sx q[3];
rz(1.7109722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4717344) q[0];
sx q[0];
rz(-0.40100455) q[0];
sx q[0];
rz(1.5420472) q[0];
rz(1.3767287) q[1];
sx q[1];
rz(-1.3841261) q[1];
sx q[1];
rz(0.9309887) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.72351) q[0];
sx q[0];
rz(-0.93659725) q[0];
sx q[0];
rz(-1.67976) q[0];
rz(-1.8225708) q[2];
sx q[2];
rz(-1.1295756) q[2];
sx q[2];
rz(3.0721498) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37759763) q[1];
sx q[1];
rz(-2.3990409) q[1];
sx q[1];
rz(1.5729891) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2263636) q[3];
sx q[3];
rz(-1.1438362) q[3];
sx q[3];
rz(2.2674436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.57637438) q[2];
sx q[2];
rz(-0.56052506) q[2];
sx q[2];
rz(3.0726688) q[2];
rz(-1.8779514) q[3];
sx q[3];
rz(-0.37216035) q[3];
sx q[3];
rz(-2.4431958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7428335) q[0];
sx q[0];
rz(-2.1837406) q[0];
sx q[0];
rz(0.051890705) q[0];
rz(0.17008153) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(-2.7896519) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1387685) q[0];
sx q[0];
rz(-1.616298) q[0];
sx q[0];
rz(0.36624927) q[0];
rz(2.3805915) q[2];
sx q[2];
rz(-0.78520757) q[2];
sx q[2];
rz(2.1987555) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4008949) q[1];
sx q[1];
rz(-2.4603421) q[1];
sx q[1];
rz(-2.2227312) q[1];
rz(-pi) q[2];
rz(1.7991583) q[3];
sx q[3];
rz(-0.31410892) q[3];
sx q[3];
rz(-0.64947646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6841782) q[2];
sx q[2];
rz(-0.64932051) q[2];
sx q[2];
rz(2.7049098) q[2];
rz(2.7501578) q[3];
sx q[3];
rz(-1.6792363) q[3];
sx q[3];
rz(2.2552538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6201685) q[0];
sx q[0];
rz(-0.7069718) q[0];
sx q[0];
rz(0.11908764) q[0];
rz(-1.2991615) q[1];
sx q[1];
rz(-1.3928587) q[1];
sx q[1];
rz(-1.7838759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085310809) q[0];
sx q[0];
rz(-1.9765903) q[0];
sx q[0];
rz(-0.88078518) q[0];
rz(-2.7804271) q[2];
sx q[2];
rz(-2.2884011) q[2];
sx q[2];
rz(-1.4360365) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5636041) q[1];
sx q[1];
rz(-1.9050026) q[1];
sx q[1];
rz(-0.6329221) q[1];
x q[2];
rz(-2.3230053) q[3];
sx q[3];
rz(-1.9789816) q[3];
sx q[3];
rz(3.0455923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.208821) q[2];
sx q[2];
rz(-1.5893156) q[2];
sx q[2];
rz(0.61441747) q[2];
rz(1.5116073) q[3];
sx q[3];
rz(-2.4698518) q[3];
sx q[3];
rz(-0.083812788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.3096302) q[0];
sx q[0];
rz(-0.93292581) q[0];
sx q[0];
rz(0.25136872) q[0];
rz(2.108719) q[1];
sx q[1];
rz(-1.947247) q[1];
sx q[1];
rz(-1.4364545) q[1];
rz(-1.4658374) q[2];
sx q[2];
rz(-1.5009673) q[2];
sx q[2];
rz(-3.0143723) q[2];
rz(-2.4182416) q[3];
sx q[3];
rz(-2.1925329) q[3];
sx q[3];
rz(-0.88919269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
