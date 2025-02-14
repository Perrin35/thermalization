OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1528435) q[0];
sx q[0];
rz(3.710521) q[0];
sx q[0];
rz(13.518128) q[0];
rz(3.9495502) q[1];
sx q[1];
rz(6.1679975) q[1];
sx q[1];
rz(8.4313784) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4848413) q[0];
sx q[0];
rz(-2.1661421) q[0];
sx q[0];
rz(-2.8369321) q[0];
rz(-2.1615218) q[2];
sx q[2];
rz(-2.9108742) q[2];
sx q[2];
rz(-1.7144817) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7639072) q[1];
sx q[1];
rz(-2.0983258) q[1];
sx q[1];
rz(-2.1930013) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4798054) q[3];
sx q[3];
rz(-1.6787623) q[3];
sx q[3];
rz(1.770883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6318165) q[2];
sx q[2];
rz(-2.5413373) q[2];
sx q[2];
rz(0.14524761) q[2];
rz(0.97233573) q[3];
sx q[3];
rz(-1.626868) q[3];
sx q[3];
rz(-1.8632896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4165322) q[0];
sx q[0];
rz(-2.5711377) q[0];
sx q[0];
rz(2.3544627) q[0];
rz(2.9540673) q[1];
sx q[1];
rz(-1.7555321) q[1];
sx q[1];
rz(3.131391) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.762741) q[0];
sx q[0];
rz(-1.6798899) q[0];
sx q[0];
rz(-0.0082882546) q[0];
rz(-pi) q[1];
x q[1];
rz(0.065930837) q[2];
sx q[2];
rz(-1.0869972) q[2];
sx q[2];
rz(0.20434665) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82380159) q[1];
sx q[1];
rz(-1.4483869) q[1];
sx q[1];
rz(-2.9106004) q[1];
rz(-pi) q[2];
rz(-0.13957406) q[3];
sx q[3];
rz(-0.8349542) q[3];
sx q[3];
rz(-2.8680108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0804704) q[2];
sx q[2];
rz(-3.0219813) q[2];
sx q[2];
rz(-1.1629026) q[2];
rz(1.2969147) q[3];
sx q[3];
rz(-1.7087405) q[3];
sx q[3];
rz(-0.0059303693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6892683) q[0];
sx q[0];
rz(-2.5732714) q[0];
sx q[0];
rz(0.34533056) q[0];
rz(0.055222424) q[1];
sx q[1];
rz(-0.60631141) q[1];
sx q[1];
rz(2.9323554) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064406618) q[0];
sx q[0];
rz(-1.6492873) q[0];
sx q[0];
rz(-2.2839943) q[0];
rz(1.8183858) q[2];
sx q[2];
rz(-0.074562975) q[2];
sx q[2];
rz(-0.32081315) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0303231) q[1];
sx q[1];
rz(-1.7220338) q[1];
sx q[1];
rz(0.72152414) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8322631) q[3];
sx q[3];
rz(-2.3785354) q[3];
sx q[3];
rz(0.3857358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.050134) q[2];
sx q[2];
rz(-1.4860934) q[2];
sx q[2];
rz(-0.44198188) q[2];
rz(0.4778536) q[3];
sx q[3];
rz(-1.7171532) q[3];
sx q[3];
rz(1.2598239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0100937) q[0];
sx q[0];
rz(-2.5690014) q[0];
sx q[0];
rz(-1.2226489) q[0];
rz(-1.6652416) q[1];
sx q[1];
rz(-1.0226701) q[1];
sx q[1];
rz(2.8388035) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3633805) q[0];
sx q[0];
rz(-1.1882093) q[0];
sx q[0];
rz(1.873111) q[0];
x q[1];
rz(0.035149375) q[2];
sx q[2];
rz(-1.4917231) q[2];
sx q[2];
rz(-3.1052232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.041639665) q[1];
sx q[1];
rz(-1.2570459) q[1];
sx q[1];
rz(-3.0932752) q[1];
rz(-2.5141671) q[3];
sx q[3];
rz(-2.1221106) q[3];
sx q[3];
rz(-2.3105636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0786232) q[2];
sx q[2];
rz(-1.4457694) q[2];
sx q[2];
rz(1.3485738) q[2];
rz(-0.013269987) q[3];
sx q[3];
rz(-0.66418663) q[3];
sx q[3];
rz(-1.9224904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.9297946) q[0];
sx q[0];
rz(-3.0199265) q[0];
sx q[0];
rz(-1.1965055) q[0];
rz(-0.78367805) q[1];
sx q[1];
rz(-1.0107026) q[1];
sx q[1];
rz(2.3616621) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257192) q[0];
sx q[0];
rz(-0.78479973) q[0];
sx q[0];
rz(0.37338202) q[0];
rz(-3.0230396) q[2];
sx q[2];
rz(-0.44008068) q[2];
sx q[2];
rz(-1.5755744) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.096141013) q[1];
sx q[1];
rz(-2.2287031) q[1];
sx q[1];
rz(-2.0631353) q[1];
rz(-pi) q[2];
rz(-1.8891009) q[3];
sx q[3];
rz(-1.0530143) q[3];
sx q[3];
rz(-0.095967231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1448867) q[2];
sx q[2];
rz(-3.0250664) q[2];
sx q[2];
rz(3.0987926) q[2];
rz(1.608611) q[3];
sx q[3];
rz(-1.3442842) q[3];
sx q[3];
rz(0.90095055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75339371) q[0];
sx q[0];
rz(-2.2008984) q[0];
sx q[0];
rz(2.4401376) q[0];
rz(0.88297168) q[1];
sx q[1];
rz(-1.7794304) q[1];
sx q[1];
rz(1.0900452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1397247) q[0];
sx q[0];
rz(-1.1764948) q[0];
sx q[0];
rz(0.35625881) q[0];
rz(2.6323363) q[2];
sx q[2];
rz(-2.351077) q[2];
sx q[2];
rz(0.59206405) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2631677) q[1];
sx q[1];
rz(-0.15053883) q[1];
sx q[1];
rz(1.9012039) q[1];
rz(-0.16441508) q[3];
sx q[3];
rz(-1.8071756) q[3];
sx q[3];
rz(0.5099834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8604454) q[2];
sx q[2];
rz(-0.88938418) q[2];
sx q[2];
rz(-1.7640198) q[2];
rz(-3.0456165) q[3];
sx q[3];
rz(-0.70982659) q[3];
sx q[3];
rz(0.53203741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9607361) q[0];
sx q[0];
rz(-2.2611698) q[0];
sx q[0];
rz(-1.2197422) q[0];
rz(2.5324054) q[1];
sx q[1];
rz(-1.5193308) q[1];
sx q[1];
rz(-2.6763197) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1022961) q[0];
sx q[0];
rz(-1.0900153) q[0];
sx q[0];
rz(-2.395438) q[0];
rz(-2.1805167) q[2];
sx q[2];
rz(-1.3111918) q[2];
sx q[2];
rz(2.0849092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6890258) q[1];
sx q[1];
rz(-1.6151607) q[1];
sx q[1];
rz(-1.1242799) q[1];
rz(-pi) q[2];
rz(1.7666529) q[3];
sx q[3];
rz(-0.80482414) q[3];
sx q[3];
rz(-2.5087961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6294127) q[2];
sx q[2];
rz(-2.8359154) q[2];
sx q[2];
rz(-2.7194887) q[2];
rz(-0.01853881) q[3];
sx q[3];
rz(-1.5471349) q[3];
sx q[3];
rz(-0.6137994) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2043532) q[0];
sx q[0];
rz(-0.48870191) q[0];
sx q[0];
rz(2.3204284) q[0];
rz(-2.4429854) q[1];
sx q[1];
rz(-0.76118529) q[1];
sx q[1];
rz(3.1331114) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5831006) q[0];
sx q[0];
rz(-1.8752397) q[0];
sx q[0];
rz(2.4598756) q[0];
x q[1];
rz(2.7445265) q[2];
sx q[2];
rz(-1.0806568) q[2];
sx q[2];
rz(-1.2063783) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.20336452) q[1];
sx q[1];
rz(-1.3044323) q[1];
sx q[1];
rz(-0.78856923) q[1];
rz(-0.04106122) q[3];
sx q[3];
rz(-0.92721894) q[3];
sx q[3];
rz(0.82541556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8734133) q[2];
sx q[2];
rz(-0.11205967) q[2];
sx q[2];
rz(-2.6123135) q[2];
rz(1.4389634) q[3];
sx q[3];
rz(-1.8301423) q[3];
sx q[3];
rz(2.8992991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.5731803) q[0];
sx q[0];
rz(-0.72887623) q[0];
sx q[0];
rz(-0.53701425) q[0];
rz(-2.2611387) q[1];
sx q[1];
rz(-1.302779) q[1];
sx q[1];
rz(0.74768487) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76067257) q[0];
sx q[0];
rz(-2.2539646) q[0];
sx q[0];
rz(1.1364514) q[0];
rz(-pi) q[1];
rz(-0.40966655) q[2];
sx q[2];
rz(-2.1958399) q[2];
sx q[2];
rz(1.1898439) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1701053) q[1];
sx q[1];
rz(-2.7344646) q[1];
sx q[1];
rz(-2.219063) q[1];
rz(-pi) q[2];
rz(-0.35028065) q[3];
sx q[3];
rz(-2.4995188) q[3];
sx q[3];
rz(0.16360276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.65413862) q[2];
sx q[2];
rz(-2.4345001) q[2];
sx q[2];
rz(-2.1716165) q[2];
rz(-1.4966494) q[3];
sx q[3];
rz(-2.1275438) q[3];
sx q[3];
rz(1.0223201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8146166) q[0];
sx q[0];
rz(-1.6642445) q[0];
sx q[0];
rz(-2.5388663) q[0];
rz(-3.1372435) q[1];
sx q[1];
rz(-1.3163722) q[1];
sx q[1];
rz(2.0306921) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2839226) q[0];
sx q[0];
rz(-0.94597781) q[0];
sx q[0];
rz(-2.3032673) q[0];
rz(-pi) q[1];
rz(-2.1793785) q[2];
sx q[2];
rz(-2.3446744) q[2];
sx q[2];
rz(-1.6730409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0266089) q[1];
sx q[1];
rz(-0.98140684) q[1];
sx q[1];
rz(1.4242054) q[1];
rz(1.2405715) q[3];
sx q[3];
rz(-2.0194411) q[3];
sx q[3];
rz(-0.90106746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6415619) q[2];
sx q[2];
rz(-1.2189453) q[2];
sx q[2];
rz(-2.5489589) q[2];
rz(-0.57989132) q[3];
sx q[3];
rz(-0.65120828) q[3];
sx q[3];
rz(-1.7013223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92689571) q[0];
sx q[0];
rz(-1.0746645) q[0];
sx q[0];
rz(-0.97434531) q[0];
rz(0.018085619) q[1];
sx q[1];
rz(-2.7985202) q[1];
sx q[1];
rz(0.090029686) q[1];
rz(-1.3152348) q[2];
sx q[2];
rz(-1.8696052) q[2];
sx q[2];
rz(2.2754696) q[2];
rz(-2.6206489) q[3];
sx q[3];
rz(-2.3250865) q[3];
sx q[3];
rz(-0.10100828) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
