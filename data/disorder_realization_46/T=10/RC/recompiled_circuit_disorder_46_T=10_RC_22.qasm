OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5912136) q[0];
sx q[0];
rz(-0.0033291078) q[0];
sx q[0];
rz(-0.34531265) q[0];
rz(-1.3357063) q[1];
sx q[1];
rz(-0.3392646) q[1];
sx q[1];
rz(0.27944922) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1169259) q[0];
sx q[0];
rz(-1.3660087) q[0];
sx q[0];
rz(-1.2555528) q[0];
rz(-pi) q[1];
rz(-0.34502132) q[2];
sx q[2];
rz(-0.81232386) q[2];
sx q[2];
rz(-1.617384) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2024723) q[1];
sx q[1];
rz(-1.3807553) q[1];
sx q[1];
rz(-3.128016) q[1];
rz(-2.4996098) q[3];
sx q[3];
rz(-1.5390453) q[3];
sx q[3];
rz(2.8905405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7001069) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(-2.1842365) q[2];
rz(2.8422614) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(-2.7295952) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9807724) q[0];
sx q[0];
rz(-2.1919724) q[0];
sx q[0];
rz(2.7278996) q[0];
rz(-1.7970239) q[1];
sx q[1];
rz(-2.3584056) q[1];
sx q[1];
rz(0.63562524) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9522063) q[0];
sx q[0];
rz(-1.4385907) q[0];
sx q[0];
rz(-3.0968451) q[0];
rz(1.7027249) q[2];
sx q[2];
rz(-1.6172098) q[2];
sx q[2];
rz(-3.1293491) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.41308584) q[1];
sx q[1];
rz(-1.8131078) q[1];
sx q[1];
rz(0.26279022) q[1];
x q[2];
rz(-0.91931822) q[3];
sx q[3];
rz(-2.161536) q[3];
sx q[3];
rz(0.35349333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3893434) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(0.52865571) q[2];
rz(1.2403437) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(-2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.56025958) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(0.2581968) q[0];
rz(-1.5064346) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(0.43446508) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6601352) q[0];
sx q[0];
rz(-2.0148206) q[0];
sx q[0];
rz(-0.40159479) q[0];
rz(-pi) q[1];
rz(-1.6824109) q[2];
sx q[2];
rz(-2.6093405) q[2];
sx q[2];
rz(0.75512952) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9745969) q[1];
sx q[1];
rz(-1.7803065) q[1];
sx q[1];
rz(-0.4619044) q[1];
x q[2];
rz(0.64819077) q[3];
sx q[3];
rz(-2.18581) q[3];
sx q[3];
rz(-2.0730413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42052856) q[2];
sx q[2];
rz(-1.3238182) q[2];
sx q[2];
rz(-0.055796441) q[2];
rz(2.2037286) q[3];
sx q[3];
rz(-2.8580229) q[3];
sx q[3];
rz(2.8985033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043512251) q[0];
sx q[0];
rz(-2.1976017) q[0];
sx q[0];
rz(-3.0010624) q[0];
rz(0.17164104) q[1];
sx q[1];
rz(-1.3069897) q[1];
sx q[1];
rz(-2.8780639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4123654) q[0];
sx q[0];
rz(-2.1365039) q[0];
sx q[0];
rz(-1.4501146) q[0];
rz(1.0657004) q[2];
sx q[2];
rz(-0.96955883) q[2];
sx q[2];
rz(0.87841735) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0553953) q[1];
sx q[1];
rz(-0.89790422) q[1];
sx q[1];
rz(1.124568) q[1];
rz(-2.8670026) q[3];
sx q[3];
rz(-1.2488741) q[3];
sx q[3];
rz(1.6695175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0984829) q[2];
sx q[2];
rz(-2.7973599) q[2];
sx q[2];
rz(1.2496703) q[2];
rz(1.194687) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(2.4627114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93976218) q[0];
sx q[0];
rz(-1.285346) q[0];
sx q[0];
rz(0.029065954) q[0];
rz(-1.7395696) q[1];
sx q[1];
rz(-0.35750917) q[1];
sx q[1];
rz(-2.8129541) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0120481) q[0];
sx q[0];
rz(-2.64376) q[0];
sx q[0];
rz(1.356202) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57428898) q[2];
sx q[2];
rz(-0.31775489) q[2];
sx q[2];
rz(2.6025835) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29341104) q[1];
sx q[1];
rz(-1.1383346) q[1];
sx q[1];
rz(1.5007988) q[1];
x q[2];
rz(-2.0629115) q[3];
sx q[3];
rz(-2.2552367) q[3];
sx q[3];
rz(-1.9198315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0236726) q[2];
sx q[2];
rz(-0.4549883) q[2];
sx q[2];
rz(2.9809791) q[2];
rz(-1.9953856) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(-0.62079352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6495431) q[0];
sx q[0];
rz(-0.8594802) q[0];
sx q[0];
rz(2.3550526) q[0];
rz(-0.37711626) q[1];
sx q[1];
rz(-0.84723324) q[1];
sx q[1];
rz(0.4424817) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0241644) q[0];
sx q[0];
rz(-2.6167343) q[0];
sx q[0];
rz(0.79458046) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2636678) q[2];
sx q[2];
rz(-2.1233635) q[2];
sx q[2];
rz(-0.27586684) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1861371) q[1];
sx q[1];
rz(-1.8147787) q[1];
sx q[1];
rz(-1.9013491) q[1];
rz(-pi) q[2];
rz(-1.8507067) q[3];
sx q[3];
rz(-2.8140682) q[3];
sx q[3];
rz(-2.1675828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4601712) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(-0.9712514) q[2];
rz(0.69532895) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(2.7142081) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63715315) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(2.1321645) q[0];
rz(2.5908453) q[1];
sx q[1];
rz(-1.2754722) q[1];
sx q[1];
rz(-1.1632464) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79272072) q[0];
sx q[0];
rz(-0.36076818) q[0];
sx q[0];
rz(1.9804721) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8613569) q[2];
sx q[2];
rz(-2.4121948) q[2];
sx q[2];
rz(1.1682208) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35112652) q[1];
sx q[1];
rz(-0.27517056) q[1];
sx q[1];
rz(0.64063425) q[1];
rz(1.936124) q[3];
sx q[3];
rz(-0.72971535) q[3];
sx q[3];
rz(-0.038289379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.78272468) q[2];
sx q[2];
rz(-1.9903368) q[2];
sx q[2];
rz(0.7981832) q[2];
rz(2.362137) q[3];
sx q[3];
rz(-0.53648406) q[3];
sx q[3];
rz(1.1727758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9629795) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(-0.12284199) q[0];
rz(3.0154862) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(1.2164446) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1232237) q[0];
sx q[0];
rz(-0.95802486) q[0];
sx q[0];
rz(-2.0809253) q[0];
rz(-pi) q[1];
rz(0.33151303) q[2];
sx q[2];
rz(-0.47795948) q[2];
sx q[2];
rz(0.31994672) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0876906) q[1];
sx q[1];
rz(-1.9551827) q[1];
sx q[1];
rz(-1.7840506) q[1];
rz(-pi) q[2];
rz(1.3835195) q[3];
sx q[3];
rz(-0.98687275) q[3];
sx q[3];
rz(2.5919979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(3.0984666) q[2];
rz(2.96636) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(-1.5568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6373428) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(2.9328226) q[0];
rz(-1.6437221) q[1];
sx q[1];
rz(-1.4893963) q[1];
sx q[1];
rz(-2.0170905) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61125206) q[0];
sx q[0];
rz(-1.7065305) q[0];
sx q[0];
rz(1.432857) q[0];
rz(-pi) q[1];
rz(-1.3051946) q[2];
sx q[2];
rz(-2.2679066) q[2];
sx q[2];
rz(3.0150989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75301718) q[1];
sx q[1];
rz(-1.57197) q[1];
sx q[1];
rz(1.7041824) q[1];
rz(-pi) q[2];
rz(-0.87499683) q[3];
sx q[3];
rz(-1.1831302) q[3];
sx q[3];
rz(-0.33867237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0008529) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(0.17803426) q[2];
rz(1.5464276) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(0.49062887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7889325) q[0];
sx q[0];
rz(-1.9910318) q[0];
sx q[0];
rz(-2.1066522) q[0];
rz(0.79822284) q[1];
sx q[1];
rz(-2.1612576) q[1];
sx q[1];
rz(2.9916874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56652727) q[0];
sx q[0];
rz(-1.2034104) q[0];
sx q[0];
rz(-2.4627986) q[0];
rz(-0.28658861) q[2];
sx q[2];
rz(-1.5550685) q[2];
sx q[2];
rz(-0.77428267) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18969892) q[1];
sx q[1];
rz(-0.67283291) q[1];
sx q[1];
rz(-0.876902) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3224162) q[3];
sx q[3];
rz(-0.90667778) q[3];
sx q[3];
rz(0.76457232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7297111) q[2];
sx q[2];
rz(-1.231266) q[2];
sx q[2];
rz(0.40851545) q[2];
rz(-2.216693) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(-2.6598721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9020486) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(-1.7655903) q[1];
sx q[1];
rz(-1.9648432) q[1];
sx q[1];
rz(1.2479938) q[1];
rz(2.3464936) q[2];
sx q[2];
rz(-1.8021402) q[2];
sx q[2];
rz(-2.6593593) q[2];
rz(-2.8077447) q[3];
sx q[3];
rz(-2.4945989) q[3];
sx q[3];
rz(0.37024959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];