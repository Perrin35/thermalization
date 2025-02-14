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
rz(0.13134512) q[0];
rz(-3.1047473) q[1];
sx q[1];
rz(3.6689833) q[1];
sx q[1];
rz(13.875516) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1037169) q[0];
sx q[0];
rz(-1.193422) q[0];
sx q[0];
rz(2.3233633) q[0];
x q[1];
rz(0.45707656) q[2];
sx q[2];
rz(-1.0731878) q[2];
sx q[2];
rz(3.0534985) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.943191) q[1];
sx q[1];
rz(-1.8089589) q[1];
sx q[1];
rz(-0.32126255) q[1];
x q[2];
rz(-2.9580826) q[3];
sx q[3];
rz(-1.402463) q[3];
sx q[3];
rz(2.3678697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5347791) q[2];
sx q[2];
rz(-0.35627347) q[2];
sx q[2];
rz(-2.479539) q[2];
rz(0.74696294) q[3];
sx q[3];
rz(-0.91619879) q[3];
sx q[3];
rz(-1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5544283) q[0];
sx q[0];
rz(-2.9968408) q[0];
sx q[0];
rz(1.1821049) q[0];
rz(-0.74554044) q[1];
sx q[1];
rz(-0.59919557) q[1];
sx q[1];
rz(-3.0381957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022522634) q[0];
sx q[0];
rz(-1.8028717) q[0];
sx q[0];
rz(0.43605767) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1076043) q[2];
sx q[2];
rz(-0.66945449) q[2];
sx q[2];
rz(0.34399271) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1704639) q[1];
sx q[1];
rz(-2.269575) q[1];
sx q[1];
rz(2.7034195) q[1];
rz(-0.17354266) q[3];
sx q[3];
rz(-2.3056917) q[3];
sx q[3];
rz(0.65803448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55249247) q[2];
sx q[2];
rz(-1.2673667) q[2];
sx q[2];
rz(0.44431552) q[2];
rz(2.0878504) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(1.9452555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.0788197) q[0];
sx q[0];
rz(-1.8906931) q[0];
sx q[0];
rz(-2.479082) q[0];
rz(1.484681) q[1];
sx q[1];
rz(-1.0228446) q[1];
sx q[1];
rz(-2.9248765) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92531918) q[0];
sx q[0];
rz(-1.1084304) q[0];
sx q[0];
rz(-2.3564649) q[0];
rz(0.7324842) q[2];
sx q[2];
rz(-1.9062796) q[2];
sx q[2];
rz(1.6215289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3308476) q[1];
sx q[1];
rz(-2.243089) q[1];
sx q[1];
rz(-0.051203392) q[1];
rz(1.4101348) q[3];
sx q[3];
rz(-1.0815797) q[3];
sx q[3];
rz(-2.1833724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.66389877) q[2];
sx q[2];
rz(-2.615216) q[2];
sx q[2];
rz(2.9065175) q[2];
rz(-2.7663686) q[3];
sx q[3];
rz(-2.2347361) q[3];
sx q[3];
rz(0.71845976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5594056) q[0];
sx q[0];
rz(-0.087787293) q[0];
sx q[0];
rz(-1.0002332) q[0];
rz(-1.2031215) q[1];
sx q[1];
rz(-1.5088046) q[1];
sx q[1];
rz(2.5968754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398374) q[0];
sx q[0];
rz(-2.4271963) q[0];
sx q[0];
rz(3.073126) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3339231) q[2];
sx q[2];
rz(-2.2414506) q[2];
sx q[2];
rz(0.18137056) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.57883731) q[1];
sx q[1];
rz(-1.3816427) q[1];
sx q[1];
rz(3.021043) q[1];
rz(-2.8417743) q[3];
sx q[3];
rz(-0.97917367) q[3];
sx q[3];
rz(-3.0772131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41877052) q[2];
sx q[2];
rz(-0.79221574) q[2];
sx q[2];
rz(-2.344632) q[2];
rz(2.8523417) q[3];
sx q[3];
rz(-0.97024337) q[3];
sx q[3];
rz(2.7204035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5274984) q[0];
sx q[0];
rz(-3.0526563) q[0];
sx q[0];
rz(1.2689137) q[0];
rz(-2.1753963) q[1];
sx q[1];
rz(-1.7485917) q[1];
sx q[1];
rz(-0.46636811) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0206576) q[0];
sx q[0];
rz(-1.0758721) q[0];
sx q[0];
rz(-1.0878956) q[0];
x q[1];
rz(2.6796273) q[2];
sx q[2];
rz(-1.8868539) q[2];
sx q[2];
rz(-2.8082928) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1524017) q[1];
sx q[1];
rz(-0.70211239) q[1];
sx q[1];
rz(2.0753808) q[1];
rz(-2.7511699) q[3];
sx q[3];
rz(-0.94518928) q[3];
sx q[3];
rz(2.1805891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40631488) q[2];
sx q[2];
rz(-2.3390528) q[2];
sx q[2];
rz(2.3373513) q[2];
rz(0.4869701) q[3];
sx q[3];
rz(-2.099497) q[3];
sx q[3];
rz(-3.13412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2020579) q[0];
sx q[0];
rz(-0.29629961) q[0];
sx q[0];
rz(2.1837088) q[0];
rz(1.5127888) q[1];
sx q[1];
rz(-2.084338) q[1];
sx q[1];
rz(-1.5527976) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2561172) q[0];
sx q[0];
rz(-1.6506356) q[0];
sx q[0];
rz(3.1186597) q[0];
rz(-pi) q[1];
rz(-0.79926305) q[2];
sx q[2];
rz(-0.38866079) q[2];
sx q[2];
rz(-2.6262019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9839638) q[1];
sx q[1];
rz(-0.92008725) q[1];
sx q[1];
rz(1.3476861) q[1];
rz(-pi) q[2];
rz(2.6090066) q[3];
sx q[3];
rz(-0.69952337) q[3];
sx q[3];
rz(-0.12040779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9396886) q[2];
sx q[2];
rz(-2.0857911) q[2];
sx q[2];
rz(-2.4064257) q[2];
rz(-0.9489263) q[3];
sx q[3];
rz(-0.85845033) q[3];
sx q[3];
rz(-2.2775876) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3801124) q[0];
sx q[0];
rz(-1.004847) q[0];
sx q[0];
rz(-0.6066221) q[0];
rz(-2.3573549) q[1];
sx q[1];
rz(-1.8083068) q[1];
sx q[1];
rz(-1.3444791) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356268) q[0];
sx q[0];
rz(-2.2971662) q[0];
sx q[0];
rz(-2.4854922) q[0];
rz(-1.8227578) q[2];
sx q[2];
rz(-2.0540385) q[2];
sx q[2];
rz(-1.9117219) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5658907) q[1];
sx q[1];
rz(-1.9676349) q[1];
sx q[1];
rz(1.8116519) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7502039) q[3];
sx q[3];
rz(-0.47166892) q[3];
sx q[3];
rz(0.52667945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6091696) q[2];
sx q[2];
rz(-0.50749856) q[2];
sx q[2];
rz(1.7519105) q[2];
rz(-0.17942795) q[3];
sx q[3];
rz(-2.1556985) q[3];
sx q[3];
rz(-0.40670407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.0824025) q[0];
sx q[0];
rz(-0.32408369) q[0];
sx q[0];
rz(-0.75505906) q[0];
rz(0.45626196) q[1];
sx q[1];
rz(-1.7165963) q[1];
sx q[1];
rz(2.0645352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74769831) q[0];
sx q[0];
rz(-2.1222881) q[0];
sx q[0];
rz(1.8389788) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6918534) q[2];
sx q[2];
rz(-2.350995) q[2];
sx q[2];
rz(-1.3644753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3216599) q[1];
sx q[1];
rz(-1.1159619) q[1];
sx q[1];
rz(-1.2534035) q[1];
x q[2];
rz(1.0302587) q[3];
sx q[3];
rz(-1.7071299) q[3];
sx q[3];
rz(-2.7928074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72851744) q[2];
sx q[2];
rz(-2.04144) q[2];
sx q[2];
rz(-2.4033578) q[2];
rz(-1.5089367) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(-1.1608605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4974834) q[0];
sx q[0];
rz(-1.6832385) q[0];
sx q[0];
rz(0.64895502) q[0];
rz(2.451918) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(2.7241657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6605292) q[0];
sx q[0];
rz(-1.821035) q[0];
sx q[0];
rz(0.60926837) q[0];
x q[1];
rz(2.8990977) q[2];
sx q[2];
rz(-0.30825492) q[2];
sx q[2];
rz(0.85390845) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6747804) q[1];
sx q[1];
rz(-1.1622475) q[1];
sx q[1];
rz(1.4571587) q[1];
x q[2];
rz(-2.7415258) q[3];
sx q[3];
rz(-1.9377982) q[3];
sx q[3];
rz(2.7719967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6680341) q[2];
sx q[2];
rz(-1.4213976) q[2];
sx q[2];
rz(3.0958946) q[2];
rz(0.33347305) q[3];
sx q[3];
rz(-0.57531753) q[3];
sx q[3];
rz(0.12588178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4561653) q[0];
sx q[0];
rz(-2.6021155) q[0];
sx q[0];
rz(-1.9239377) q[0];
rz(-2.894891) q[1];
sx q[1];
rz(-1.8496937) q[1];
sx q[1];
rz(-2.6620679) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2829153) q[0];
sx q[0];
rz(-1.1395795) q[0];
sx q[0];
rz(-0.89003508) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2662451) q[2];
sx q[2];
rz(-1.1743059) q[2];
sx q[2];
rz(3.0975395) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.61726924) q[1];
sx q[1];
rz(-1.8620544) q[1];
sx q[1];
rz(0.53415438) q[1];
x q[2];
rz(-1.3334951) q[3];
sx q[3];
rz(-2.0249244) q[3];
sx q[3];
rz(0.25179201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.55629998) q[2];
sx q[2];
rz(-1.2703398) q[2];
sx q[2];
rz(1.7362107) q[2];
rz(0.65226883) q[3];
sx q[3];
rz(-0.26572078) q[3];
sx q[3];
rz(2.4238267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
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
rz(-2.5581735) q[2];
sx q[2];
rz(-1.7023682) q[2];
sx q[2];
rz(-2.792991) q[2];
rz(-0.99540972) q[3];
sx q[3];
rz(-1.8707471) q[3];
sx q[3];
rz(-1.7729014) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
