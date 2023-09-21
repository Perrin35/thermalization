OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(3.3189964) q[0];
sx q[0];
rz(11.431974) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(0.66361767) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463319) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(2.6803826) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2167591) q[2];
sx q[2];
rz(-2.3166222) q[2];
sx q[2];
rz(-1.1302469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.692688) q[1];
sx q[1];
rz(-1.2346039) q[1];
sx q[1];
rz(-0.2775788) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.033028) q[3];
sx q[3];
rz(-1.3285471) q[3];
sx q[3];
rz(0.066075174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87876451) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(-1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(0.59659514) q[0];
rz(2.3157628) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(1.2260431) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85428836) q[0];
sx q[0];
rz(-2.0137631) q[0];
sx q[0];
rz(1.3501549) q[0];
rz(-2.3184899) q[2];
sx q[2];
rz(-2.1777993) q[2];
sx q[2];
rz(-0.58128231) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5920168) q[1];
sx q[1];
rz(-1.5703652) q[1];
sx q[1];
rz(-1.3509343) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1321208) q[3];
sx q[3];
rz(-0.99820271) q[3];
sx q[3];
rz(-3.115311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(-0.80667574) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9817292) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(-1.8925517) q[0];
rz(-0.088009134) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(-1.0294611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80699608) q[0];
sx q[0];
rz(-2.1801729) q[0];
sx q[0];
rz(0.2683123) q[0];
rz(-2.7715893) q[2];
sx q[2];
rz(-0.68379935) q[2];
sx q[2];
rz(-1.5918658) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0975115) q[1];
sx q[1];
rz(-1.4315726) q[1];
sx q[1];
rz(0.15686762) q[1];
rz(-pi) q[2];
rz(-1.6920407) q[3];
sx q[3];
rz(-0.83449927) q[3];
sx q[3];
rz(-2.8672225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.007894667) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(2.6339445) q[2];
rz(-1.3890022) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65748173) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(-1.5959651) q[0];
rz(-1.0428628) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(1.625659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0632616) q[0];
sx q[0];
rz(-2.445979) q[0];
sx q[0];
rz(2.9177833) q[0];
x q[1];
rz(-2.1298725) q[2];
sx q[2];
rz(-1.6435197) q[2];
sx q[2];
rz(2.5474472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2559291) q[1];
sx q[1];
rz(-0.71223488) q[1];
sx q[1];
rz(-0.42339143) q[1];
x q[2];
rz(-2.5528615) q[3];
sx q[3];
rz(-1.7541459) q[3];
sx q[3];
rz(0.4962894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76379124) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(-3.0002248) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(-0.98658371) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35448733) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(-1.6217344) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(-0.89486665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2762404) q[0];
sx q[0];
rz(-1.9812752) q[0];
sx q[0];
rz(-0.68666896) q[0];
x q[1];
rz(1.0084444) q[2];
sx q[2];
rz(-2.214606) q[2];
sx q[2];
rz(-1.84562) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.85256165) q[1];
sx q[1];
rz(-2.3173996) q[1];
sx q[1];
rz(-0.95552163) q[1];
rz(-pi) q[2];
rz(-0.28646333) q[3];
sx q[3];
rz(-1.5634007) q[3];
sx q[3];
rz(-1.3163819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.616509) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(1.9990702) q[2];
rz(-0.74674314) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(-2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58364761) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(-0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(-2.6766052) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87100055) q[0];
sx q[0];
rz(-0.56068476) q[0];
sx q[0];
rz(2.2339348) q[0];
rz(-pi) q[1];
rz(-0.34157413) q[2];
sx q[2];
rz(-1.6401059) q[2];
sx q[2];
rz(2.9937033) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7558414) q[1];
sx q[1];
rz(-2.5045536) q[1];
sx q[1];
rz(-1.4904651) q[1];
x q[2];
rz(-1.7808077) q[3];
sx q[3];
rz(-0.87267733) q[3];
sx q[3];
rz(-0.39892808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.650699) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(0.4450376) q[2];
rz(-2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12748195) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(-1.7215464) q[0];
rz(3.1177915) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(-2.9856317) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10931817) q[0];
sx q[0];
rz(-0.30800691) q[0];
sx q[0];
rz(-0.59662915) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67955741) q[2];
sx q[2];
rz(-2.5153) q[2];
sx q[2];
rz(1.2223513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8487726) q[1];
sx q[1];
rz(-0.18310586) q[1];
sx q[1];
rz(-1.2700901) q[1];
rz(-pi) q[2];
rz(0.74219269) q[3];
sx q[3];
rz(-0.67897292) q[3];
sx q[3];
rz(-1.2364482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(-1.0726661) q[2];
rz(-0.32564751) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(-1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1919365) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(-2.8038213) q[0];
rz(2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(0.24857323) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6876858) q[0];
sx q[0];
rz(-2.0810063) q[0];
sx q[0];
rz(1.2841671) q[0];
rz(-pi) q[1];
rz(-2.0806662) q[2];
sx q[2];
rz(-1.14398) q[2];
sx q[2];
rz(0.61444297) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8057115) q[1];
sx q[1];
rz(-2.1599249) q[1];
sx q[1];
rz(-0.59405234) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0646938) q[3];
sx q[3];
rz(-1.7401164) q[3];
sx q[3];
rz(2.3016735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(-3.0548813) q[2];
rz(0.48197204) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.988525) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865737) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(2.8588262) q[0];
rz(-0.70156082) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(1.823002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222262) q[0];
sx q[0];
rz(-2.0113693) q[0];
sx q[0];
rz(-0.088935436) q[0];
rz(-2.4852072) q[2];
sx q[2];
rz(-1.0644541) q[2];
sx q[2];
rz(2.2040747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9642155) q[1];
sx q[1];
rz(-2.4417158) q[1];
sx q[1];
rz(-1.7978976) q[1];
rz(-pi) q[2];
rz(2.3913279) q[3];
sx q[3];
rz(-1.8599469) q[3];
sx q[3];
rz(-2.9671448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7302154) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(-2.7424157) q[2];
rz(-2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76319641) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(1.5989074) q[0];
rz(1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(-1.8803966) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48860088) q[0];
sx q[0];
rz(-1.6761259) q[0];
sx q[0];
rz(2.7306042) q[0];
rz(-pi) q[1];
rz(0.032632685) q[2];
sx q[2];
rz(-0.59851461) q[2];
sx q[2];
rz(1.061071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2640472) q[1];
sx q[1];
rz(-2.0685158) q[1];
sx q[1];
rz(1.643814) q[1];
rz(-pi) q[2];
rz(1.4114755) q[3];
sx q[3];
rz(-2.2485002) q[3];
sx q[3];
rz(-2.2417828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.5579055) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(0.56979257) q[2];
rz(-1.2184881) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31310836) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(0.60824153) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(-2.1424978) q[2];
sx q[2];
rz(-1.5928762) q[2];
sx q[2];
rz(1.5230509) q[2];
rz(-2.9363587) q[3];
sx q[3];
rz(-2.5406465) q[3];
sx q[3];
rz(3.061486) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];