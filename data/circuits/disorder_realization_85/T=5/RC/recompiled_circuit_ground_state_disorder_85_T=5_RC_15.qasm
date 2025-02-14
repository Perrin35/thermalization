OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0394734) q[0];
sx q[0];
rz(-1.4730299) q[0];
sx q[0];
rz(-3.0102475) q[0];
rz(0.036845358) q[1];
sx q[1];
rz(-0.52739066) q[1];
sx q[1];
rz(1.8324469) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2764171) q[0];
sx q[0];
rz(-0.88224471) q[0];
sx q[0];
rz(0.49746969) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6845161) q[2];
sx q[2];
rz(-1.0731878) q[2];
sx q[2];
rz(-0.08809419) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6908474) q[1];
sx q[1];
rz(-1.8826797) q[1];
sx q[1];
rz(1.3203095) q[1];
rz(0.18351002) q[3];
sx q[3];
rz(-1.7391296) q[3];
sx q[3];
rz(0.77372293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6068136) q[2];
sx q[2];
rz(-2.7853192) q[2];
sx q[2];
rz(2.479539) q[2];
rz(2.3946297) q[3];
sx q[3];
rz(-2.2253939) q[3];
sx q[3];
rz(-1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58716431) q[0];
sx q[0];
rz(-2.9968408) q[0];
sx q[0];
rz(1.1821049) q[0];
rz(0.74554044) q[1];
sx q[1];
rz(-0.59919557) q[1];
sx q[1];
rz(3.0381957) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4415084) q[0];
sx q[0];
rz(-1.9943976) q[0];
sx q[0];
rz(-1.3157428) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1076043) q[2];
sx q[2];
rz(-2.4721382) q[2];
sx q[2];
rz(-0.34399271) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7999801) q[1];
sx q[1];
rz(-2.3368521) q[1];
sx q[1];
rz(2.0384203) q[1];
rz(-2.96805) q[3];
sx q[3];
rz(-2.3056917) q[3];
sx q[3];
rz(2.4835582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5891002) q[2];
sx q[2];
rz(-1.2673667) q[2];
sx q[2];
rz(0.44431552) q[2];
rz(1.0537423) q[3];
sx q[3];
rz(-3.0709303) q[3];
sx q[3];
rz(-1.1963371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0788197) q[0];
sx q[0];
rz(-1.2508996) q[0];
sx q[0];
rz(-2.479082) q[0];
rz(-1.6569116) q[1];
sx q[1];
rz(-1.0228446) q[1];
sx q[1];
rz(0.2167162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22610006) q[0];
sx q[0];
rz(-0.88577548) q[0];
sx q[0];
rz(-2.1846376) q[0];
rz(-pi) q[1];
rz(2.0092874) q[2];
sx q[2];
rz(-0.88738933) q[2];
sx q[2];
rz(0.33861288) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.89284929) q[1];
sx q[1];
rz(-0.67393747) q[1];
sx q[1];
rz(-1.6350063) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29189887) q[3];
sx q[3];
rz(-2.6287098) q[3];
sx q[3];
rz(1.8512902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4776939) q[2];
sx q[2];
rz(-0.52637664) q[2];
sx q[2];
rz(-2.9065175) q[2];
rz(-0.37522405) q[3];
sx q[3];
rz(-0.90685654) q[3];
sx q[3];
rz(0.71845976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.5594056) q[0];
sx q[0];
rz(-3.0538054) q[0];
sx q[0];
rz(2.1413595) q[0];
rz(1.2031215) q[1];
sx q[1];
rz(-1.5088046) q[1];
sx q[1];
rz(-2.5968754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2208015) q[0];
sx q[0];
rz(-1.615633) q[0];
sx q[0];
rz(0.71323552) q[0];
x q[1];
rz(-0.83196173) q[2];
sx q[2];
rz(-2.1429981) q[2];
sx q[2];
rz(-1.2158647) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96918584) q[1];
sx q[1];
rz(-1.6891857) q[1];
sx q[1];
rz(-1.7612996) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29981837) q[3];
sx q[3];
rz(-0.97917367) q[3];
sx q[3];
rz(-0.064379582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41877052) q[2];
sx q[2];
rz(-0.79221574) q[2];
sx q[2];
rz(-0.79696068) q[2];
rz(0.289251) q[3];
sx q[3];
rz(-2.1713493) q[3];
sx q[3];
rz(-0.42118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5274984) q[0];
sx q[0];
rz(-0.088936381) q[0];
sx q[0];
rz(1.8726789) q[0];
rz(2.1753963) q[1];
sx q[1];
rz(-1.7485917) q[1];
sx q[1];
rz(-2.6752245) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20579977) q[0];
sx q[0];
rz(-1.9917492) q[0];
sx q[0];
rz(-2.5942786) q[0];
rz(0.46196533) q[2];
sx q[2];
rz(-1.2547387) q[2];
sx q[2];
rz(0.3332999) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6154229) q[1];
sx q[1];
rz(-2.1716482) q[1];
sx q[1];
rz(-2.7533965) q[1];
rz(-2.0556695) q[3];
sx q[3];
rz(-2.4182662) q[3];
sx q[3];
rz(-1.5679899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40631488) q[2];
sx q[2];
rz(-2.3390528) q[2];
sx q[2];
rz(0.80424133) q[2];
rz(2.6546226) q[3];
sx q[3];
rz(-2.099497) q[3];
sx q[3];
rz(-0.0074726661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.2020579) q[0];
sx q[0];
rz(-0.29629961) q[0];
sx q[0];
rz(-0.95788389) q[0];
rz(1.6288039) q[1];
sx q[1];
rz(-2.084338) q[1];
sx q[1];
rz(-1.588795) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4544425) q[0];
sx q[0];
rz(-1.5936562) q[0];
sx q[0];
rz(1.4909362) q[0];
rz(2.8634809) q[2];
sx q[2];
rz(-1.8459003) q[2];
sx q[2];
rz(0.29488568) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.342345) q[1];
sx q[1];
rz(-2.4589897) q[1];
sx q[1];
rz(-0.28283878) q[1];
x q[2];
rz(1.1670052) q[3];
sx q[3];
rz(-0.98282114) q[3];
sx q[3];
rz(-0.53596067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2019041) q[2];
sx q[2];
rz(-1.0558015) q[2];
sx q[2];
rz(-2.4064257) q[2];
rz(0.9489263) q[3];
sx q[3];
rz(-0.85845033) q[3];
sx q[3];
rz(-0.86400509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(2.3801124) q[0];
sx q[0];
rz(-2.1367456) q[0];
sx q[0];
rz(-0.6066221) q[0];
rz(0.78423777) q[1];
sx q[1];
rz(-1.3332858) q[1];
sx q[1];
rz(1.3444791) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3494251) q[0];
sx q[0];
rz(-1.0972404) q[0];
sx q[0];
rz(0.72832941) q[0];
rz(-0.49656252) q[2];
sx q[2];
rz(-1.7933868) q[2];
sx q[2];
rz(0.45997657) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91025464) q[1];
sx q[1];
rz(-1.7926072) q[1];
sx q[1];
rz(2.7342058) q[1];
rz(-pi) q[2];
rz(1.7502039) q[3];
sx q[3];
rz(-2.6699237) q[3];
sx q[3];
rz(2.6149132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.532423) q[2];
sx q[2];
rz(-2.6340941) q[2];
sx q[2];
rz(1.7519105) q[2];
rz(0.17942795) q[3];
sx q[3];
rz(-0.98589412) q[3];
sx q[3];
rz(-0.40670407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0591902) q[0];
sx q[0];
rz(-2.817509) q[0];
sx q[0];
rz(-2.3865336) q[0];
rz(-2.6853307) q[1];
sx q[1];
rz(-1.7165963) q[1];
sx q[1];
rz(2.0645352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9108601) q[0];
sx q[0];
rz(-2.5344779) q[0];
sx q[0];
rz(-2.7348628) q[0];
x q[1];
rz(0.12142373) q[2];
sx q[2];
rz(-0.7875663) q[2];
sx q[2];
rz(-1.1932258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67959784) q[1];
sx q[1];
rz(-2.5933929) q[1];
sx q[1];
rz(-2.5735709) q[1];
x q[2];
rz(-2.1113339) q[3];
sx q[3];
rz(-1.4344627) q[3];
sx q[3];
rz(-0.34878525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4130752) q[2];
sx q[2];
rz(-2.04144) q[2];
sx q[2];
rz(2.4033578) q[2];
rz(1.5089367) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(-1.9807321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.64410925) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(2.4926376) q[0];
rz(-2.451918) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(0.41742691) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48106347) q[0];
sx q[0];
rz(-1.821035) q[0];
sx q[0];
rz(0.60926837) q[0];
x q[1];
rz(1.4944878) q[2];
sx q[2];
rz(-1.2718437) q[2];
sx q[2];
rz(2.0336917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0829186) q[1];
sx q[1];
rz(-1.4665468) q[1];
sx q[1];
rz(-2.7306795) q[1];
rz(-pi) q[2];
rz(-2.7415258) q[3];
sx q[3];
rz(-1.2037945) q[3];
sx q[3];
rz(0.36959592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6680341) q[2];
sx q[2];
rz(-1.4213976) q[2];
sx q[2];
rz(0.045698015) q[2];
rz(2.8081196) q[3];
sx q[3];
rz(-0.57531753) q[3];
sx q[3];
rz(3.0157109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4561653) q[0];
sx q[0];
rz(-0.5394772) q[0];
sx q[0];
rz(1.217655) q[0];
rz(-2.894891) q[1];
sx q[1];
rz(-1.8496937) q[1];
sx q[1];
rz(0.47952476) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8586774) q[0];
sx q[0];
rz(-2.0020131) q[0];
sx q[0];
rz(-2.2515576) q[0];
rz(-0.99202435) q[2];
sx q[2];
rz(-2.357748) q[2];
sx q[2];
rz(1.0933361) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.61726924) q[1];
sx q[1];
rz(-1.2795382) q[1];
sx q[1];
rz(-2.6074383) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6761618) q[3];
sx q[3];
rz(-1.7836535) q[3];
sx q[3];
rz(1.2133017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55629998) q[2];
sx q[2];
rz(-1.8712529) q[2];
sx q[2];
rz(1.7362107) q[2];
rz(-0.65226883) q[3];
sx q[3];
rz(-2.8758719) q[3];
sx q[3];
rz(-0.71776596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2300867) q[0];
sx q[0];
rz(-1.7510887) q[0];
sx q[0];
rz(-0.44816309) q[0];
rz(-2.3975092) q[1];
sx q[1];
rz(-2.4241445) q[1];
sx q[1];
rz(-1.4345899) q[1];
rz(2.5581735) q[2];
sx q[2];
rz(-1.4392244) q[2];
sx q[2];
rz(0.34860162) q[2];
rz(2.7884095) q[3];
sx q[3];
rz(-2.1174869) q[3];
sx q[3];
rz(-0.01275851) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
