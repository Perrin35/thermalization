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
rz(1.5481663) q[0];
sx q[0];
rz(3.610008) q[0];
sx q[0];
rz(9.2863291) q[0];
rz(1.8847213) q[1];
sx q[1];
rz(4.2121834) q[1];
sx q[1];
rz(9.4033006) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81651743) q[0];
sx q[0];
rz(-1.7735574) q[0];
sx q[0];
rz(1.9295503) q[0];
rz(-pi) q[1];
rz(-2.3214705) q[2];
sx q[2];
rz(-2.409365) q[2];
sx q[2];
rz(-0.9257462) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5925046) q[1];
sx q[1];
rz(-0.23356423) q[1];
sx q[1];
rz(-0.84850581) q[1];
rz(-1.2091694) q[3];
sx q[3];
rz(-2.1361793) q[3];
sx q[3];
rz(-2.9361182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5452177) q[2];
sx q[2];
rz(-0.45520982) q[2];
sx q[2];
rz(-1.5111766) q[2];
rz(-0.049985416) q[3];
sx q[3];
rz(-1.9129916) q[3];
sx q[3];
rz(-0.25240067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91994691) q[0];
sx q[0];
rz(-1.9939461) q[0];
sx q[0];
rz(-0.037121437) q[0];
rz(-0.014017398) q[1];
sx q[1];
rz(-2.5633096) q[1];
sx q[1];
rz(-0.23606539) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35821298) q[0];
sx q[0];
rz(-1.5851846) q[0];
sx q[0];
rz(0.074812513) q[0];
rz(-pi) q[1];
rz(2.5949941) q[2];
sx q[2];
rz(-0.79397092) q[2];
sx q[2];
rz(-3.1304718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29480793) q[1];
sx q[1];
rz(-2.4681512) q[1];
sx q[1];
rz(0.93017857) q[1];
x q[2];
rz(2.8909952) q[3];
sx q[3];
rz(-1.9204233) q[3];
sx q[3];
rz(2.6299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1290805) q[2];
sx q[2];
rz(-1.7725638) q[2];
sx q[2];
rz(3.1079666) q[2];
rz(2.8513837) q[3];
sx q[3];
rz(-2.2955194) q[3];
sx q[3];
rz(0.39053759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6212807) q[0];
sx q[0];
rz(-0.97977591) q[0];
sx q[0];
rz(-2.8424971) q[0];
rz(1.6200804) q[1];
sx q[1];
rz(-0.027093096) q[1];
sx q[1];
rz(3.1153968) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25616249) q[0];
sx q[0];
rz(-3.1195398) q[0];
sx q[0];
rz(1.8697478) q[0];
rz(-1.2904335) q[2];
sx q[2];
rz(-1.7086662) q[2];
sx q[2];
rz(3.1119135) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.67945665) q[1];
sx q[1];
rz(-0.95210451) q[1];
sx q[1];
rz(2.2652744) q[1];
rz(-pi) q[2];
rz(-1.800816) q[3];
sx q[3];
rz(-1.134308) q[3];
sx q[3];
rz(0.088584049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.68411487) q[2];
sx q[2];
rz(-0.44197765) q[2];
sx q[2];
rz(-0.6074062) q[2];
rz(0.36951798) q[3];
sx q[3];
rz(-1.4990467) q[3];
sx q[3];
rz(-3.1336237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5615416) q[0];
sx q[0];
rz(-2.9396368) q[0];
sx q[0];
rz(1.1308905) q[0];
rz(-2.789403) q[1];
sx q[1];
rz(-2.6512841) q[1];
sx q[1];
rz(2.9255548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7894831) q[0];
sx q[0];
rz(-2.0228068) q[0];
sx q[0];
rz(1.3003967) q[0];
rz(-pi) q[1];
rz(-0.44682002) q[2];
sx q[2];
rz(-1.6761314) q[2];
sx q[2];
rz(0.0086431816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6747848) q[1];
sx q[1];
rz(-0.66716107) q[1];
sx q[1];
rz(2.2391367) q[1];
rz(-pi) q[2];
rz(-1.0918255) q[3];
sx q[3];
rz(-0.19102016) q[3];
sx q[3];
rz(-0.52680069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.617368) q[2];
sx q[2];
rz(-2.2286131) q[2];
sx q[2];
rz(-0.53140223) q[2];
rz(1.4517387) q[3];
sx q[3];
rz(-2.6302591) q[3];
sx q[3];
rz(2.1112198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3568929) q[0];
sx q[0];
rz(-2.0578616) q[0];
sx q[0];
rz(0.30237958) q[0];
rz(2.0284292) q[1];
sx q[1];
rz(-2.1402363) q[1];
sx q[1];
rz(1.0611634) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7776124) q[0];
sx q[0];
rz(-1.5060695) q[0];
sx q[0];
rz(1.510468) q[0];
rz(1.6689014) q[2];
sx q[2];
rz(-1.3472234) q[2];
sx q[2];
rz(-1.443271) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9446774) q[1];
sx q[1];
rz(-0.79087559) q[1];
sx q[1];
rz(0.18038919) q[1];
rz(-pi) q[2];
rz(-0.69740422) q[3];
sx q[3];
rz(-2.7602969) q[3];
sx q[3];
rz(-0.21496102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82659668) q[2];
sx q[2];
rz(-1.7642517) q[2];
sx q[2];
rz(-2.946741) q[2];
rz(-2.3991614) q[3];
sx q[3];
rz(-2.4104379) q[3];
sx q[3];
rz(-0.2754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18405296) q[0];
sx q[0];
rz(-0.66127151) q[0];
sx q[0];
rz(1.1141962) q[0];
rz(-2.8320352) q[1];
sx q[1];
rz(-2.20859) q[1];
sx q[1];
rz(-1.385744) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4969869) q[0];
sx q[0];
rz(-0.33503767) q[0];
sx q[0];
rz(-1.0100288) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3676332) q[2];
sx q[2];
rz(-0.51649714) q[2];
sx q[2];
rz(1.9488283) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0506865) q[1];
sx q[1];
rz(-0.93877568) q[1];
sx q[1];
rz(-1.2463831) q[1];
rz(-0.77499109) q[3];
sx q[3];
rz(-2.0609566) q[3];
sx q[3];
rz(-0.22966269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.679057) q[2];
sx q[2];
rz(-2.6947196) q[2];
sx q[2];
rz(2.7317969) q[2];
rz(-3.0637686) q[3];
sx q[3];
rz(-1.2160559) q[3];
sx q[3];
rz(2.3791544) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78793144) q[0];
sx q[0];
rz(-0.15058148) q[0];
sx q[0];
rz(-1.6702363) q[0];
rz(2.3279066) q[1];
sx q[1];
rz(-0.7209456) q[1];
sx q[1];
rz(0.90027666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0896331) q[0];
sx q[0];
rz(-1.51043) q[0];
sx q[0];
rz(-3.0230762) q[0];
rz(3.056053) q[2];
sx q[2];
rz(-1.2844049) q[2];
sx q[2];
rz(-2.6154892) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1623692) q[1];
sx q[1];
rz(-1.3047393) q[1];
sx q[1];
rz(0.0084657808) q[1];
rz(-1.6358709) q[3];
sx q[3];
rz(-1.3378007) q[3];
sx q[3];
rz(-0.20291337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.5208931) q[2];
sx q[2];
rz(-0.40105477) q[2];
sx q[2];
rz(2.9996784) q[2];
rz(-2.9698931) q[3];
sx q[3];
rz(-1.9008235) q[3];
sx q[3];
rz(-2.5920674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5080344) q[0];
sx q[0];
rz(-1.9867851) q[0];
sx q[0];
rz(-1.5800193) q[0];
rz(-0.70478565) q[1];
sx q[1];
rz(-0.90122688) q[1];
sx q[1];
rz(-2.8660692) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4954105) q[0];
sx q[0];
rz(-0.14178628) q[0];
sx q[0];
rz(-2.0210193) q[0];
x q[1];
rz(0.24532206) q[2];
sx q[2];
rz(-1.5333912) q[2];
sx q[2];
rz(0.87329162) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.91159484) q[1];
sx q[1];
rz(-1.8640638) q[1];
sx q[1];
rz(-1.0469584) q[1];
rz(-pi) q[2];
rz(1.9884692) q[3];
sx q[3];
rz(-2.1942007) q[3];
sx q[3];
rz(-0.88043326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0602818) q[2];
sx q[2];
rz(-1.5014481) q[2];
sx q[2];
rz(-2.8509129) q[2];
rz(-0.79689133) q[3];
sx q[3];
rz(-2.6698038) q[3];
sx q[3];
rz(-0.98381388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.31442916) q[0];
sx q[0];
rz(-1.6570579) q[0];
sx q[0];
rz(-2.2737801) q[0];
rz(-2.9293291) q[1];
sx q[1];
rz(-0.8808732) q[1];
sx q[1];
rz(0.27264047) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6295687) q[0];
sx q[0];
rz(-1.6457412) q[0];
sx q[0];
rz(-1.3550979) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10137239) q[2];
sx q[2];
rz(-2.3194312) q[2];
sx q[2];
rz(0.062495141) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.18913774) q[1];
sx q[1];
rz(-2.5057372) q[1];
sx q[1];
rz(-1.5067504) q[1];
x q[2];
rz(-0.41447422) q[3];
sx q[3];
rz(-1.2948841) q[3];
sx q[3];
rz(1.7981573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.049204443) q[2];
sx q[2];
rz(-0.35085756) q[2];
sx q[2];
rz(1.997088) q[2];
rz(-0.62659621) q[3];
sx q[3];
rz(-0.75855362) q[3];
sx q[3];
rz(2.8265317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7476244) q[0];
sx q[0];
rz(-2.4691041) q[0];
sx q[0];
rz(2.5482063) q[0];
rz(2.4001832) q[1];
sx q[1];
rz(-0.64677042) q[1];
sx q[1];
rz(1.5945565) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56549673) q[0];
sx q[0];
rz(-1.4617697) q[0];
sx q[0];
rz(1.7897357) q[0];
x q[1];
rz(-2.7266961) q[2];
sx q[2];
rz(-1.3438359) q[2];
sx q[2];
rz(1.0912053) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50362557) q[1];
sx q[1];
rz(-2.4071998) q[1];
sx q[1];
rz(-3.0375558) q[1];
rz(-0.029742777) q[3];
sx q[3];
rz(-1.3689201) q[3];
sx q[3];
rz(-1.0591398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0934304) q[2];
sx q[2];
rz(-0.54485816) q[2];
sx q[2];
rz(0.94259924) q[2];
rz(-1.7132828) q[3];
sx q[3];
rz(-1.2497679) q[3];
sx q[3];
rz(1.1552756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9485332) q[0];
sx q[0];
rz(-1.5848703) q[0];
sx q[0];
rz(1.7901044) q[0];
rz(1.9652741) q[1];
sx q[1];
rz(-0.90570025) q[1];
sx q[1];
rz(-0.66014231) q[1];
rz(1.0084441) q[2];
sx q[2];
rz(-0.32505072) q[2];
sx q[2];
rz(1.6735003) q[2];
rz(-1.3073714) q[3];
sx q[3];
rz(-2.0888804) q[3];
sx q[3];
rz(-2.9874939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
