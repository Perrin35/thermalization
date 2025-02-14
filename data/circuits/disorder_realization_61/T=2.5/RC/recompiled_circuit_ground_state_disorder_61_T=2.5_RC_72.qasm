OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9899848) q[0];
sx q[0];
rz(-2.2014872) q[0];
sx q[0];
rz(-0.54036933) q[0];
rz(-2.6481533) q[1];
sx q[1];
rz(-0.72055888) q[1];
sx q[1];
rz(-2.5277353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2342398) q[0];
sx q[0];
rz(-0.94830238) q[0];
sx q[0];
rz(1.3433775) q[0];
rz(-pi) q[1];
rz(1.5510606) q[2];
sx q[2];
rz(-0.29183772) q[2];
sx q[2];
rz(-2.8449051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.39345523) q[1];
sx q[1];
rz(-2.2917213) q[1];
sx q[1];
rz(-1.9490521) q[1];
rz(1.5924388) q[3];
sx q[3];
rz(-1.8179245) q[3];
sx q[3];
rz(0.46268845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2962239) q[2];
sx q[2];
rz(-1.2426528) q[2];
sx q[2];
rz(2.5073012) q[2];
rz(2.2058709) q[3];
sx q[3];
rz(-0.24565419) q[3];
sx q[3];
rz(0.74265695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82035404) q[0];
sx q[0];
rz(-1.6004434) q[0];
sx q[0];
rz(-0.4441922) q[0];
rz(-2.3815637) q[1];
sx q[1];
rz(-1.1385695) q[1];
sx q[1];
rz(-0.98145032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6985921) q[0];
sx q[0];
rz(-0.60054296) q[0];
sx q[0];
rz(-0.53437676) q[0];
x q[1];
rz(-1.3992127) q[2];
sx q[2];
rz(-1.7497471) q[2];
sx q[2];
rz(1.6323665) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3722575) q[1];
sx q[1];
rz(-2.6943639) q[1];
sx q[1];
rz(2.9370537) q[1];
rz(0.53175521) q[3];
sx q[3];
rz(-0.49177836) q[3];
sx q[3];
rz(-2.7454706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0280219) q[2];
sx q[2];
rz(-0.61442033) q[2];
sx q[2];
rz(-1.2051955) q[2];
rz(0.53660721) q[3];
sx q[3];
rz(-1.9034932) q[3];
sx q[3];
rz(-1.7369778) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9179012) q[0];
sx q[0];
rz(-1.8603928) q[0];
sx q[0];
rz(-2.3028288) q[0];
rz(-0.63610786) q[1];
sx q[1];
rz(-1.4898224) q[1];
sx q[1];
rz(3.1210693) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3588226) q[0];
sx q[0];
rz(-2.2789445) q[0];
sx q[0];
rz(-2.3282611) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0383561) q[2];
sx q[2];
rz(-2.1439432) q[2];
sx q[2];
rz(1.8351042) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5387303) q[1];
sx q[1];
rz(-2.7285721) q[1];
sx q[1];
rz(-3.1228423) q[1];
rz(-pi) q[2];
rz(2.0666615) q[3];
sx q[3];
rz(-2.819546) q[3];
sx q[3];
rz(-2.5888458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8013578) q[2];
sx q[2];
rz(-1.5284208) q[2];
sx q[2];
rz(-0.94432962) q[2];
rz(-1.0229735) q[3];
sx q[3];
rz(-1.2228271) q[3];
sx q[3];
rz(2.003722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706547) q[0];
sx q[0];
rz(-1.174467) q[0];
sx q[0];
rz(-1.9842072) q[0];
rz(-1.6429139) q[1];
sx q[1];
rz(-1.6694371) q[1];
sx q[1];
rz(1.1882943) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0152215) q[0];
sx q[0];
rz(-1.0980716) q[0];
sx q[0];
rz(-1.0593828) q[0];
x q[1];
rz(-1.525773) q[2];
sx q[2];
rz(-1.0226215) q[2];
sx q[2];
rz(-0.5612824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7674347) q[1];
sx q[1];
rz(-1.0113455) q[1];
sx q[1];
rz(2.4591695) q[1];
rz(-pi) q[2];
rz(-1.1340302) q[3];
sx q[3];
rz(-0.371007) q[3];
sx q[3];
rz(1.8414258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1184065) q[2];
sx q[2];
rz(-1.1772757) q[2];
sx q[2];
rz(0.6905306) q[2];
rz(-2.0265419) q[3];
sx q[3];
rz(-2.3670022) q[3];
sx q[3];
rz(-1.640865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046415064) q[0];
sx q[0];
rz(-1.7061808) q[0];
sx q[0];
rz(2.114356) q[0];
rz(-1.3014303) q[1];
sx q[1];
rz(-0.65168989) q[1];
sx q[1];
rz(-0.32381907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1808609) q[0];
sx q[0];
rz(-1.9692076) q[0];
sx q[0];
rz(-0.84503998) q[0];
rz(-0.01417966) q[2];
sx q[2];
rz(-1.018033) q[2];
sx q[2];
rz(-1.9270037) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71621543) q[1];
sx q[1];
rz(-1.0310804) q[1];
sx q[1];
rz(-1.9464689) q[1];
x q[2];
rz(0.15669723) q[3];
sx q[3];
rz(-2.0771871) q[3];
sx q[3];
rz(-1.8739669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7258437) q[2];
sx q[2];
rz(-2.9314633) q[2];
sx q[2];
rz(-0.48247639) q[2];
rz(-1.1680565) q[3];
sx q[3];
rz(-1.9246293) q[3];
sx q[3];
rz(-2.4647958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2060858) q[0];
sx q[0];
rz(-2.2914903) q[0];
sx q[0];
rz(-2.2555943) q[0];
rz(2.50792) q[1];
sx q[1];
rz(-0.88992563) q[1];
sx q[1];
rz(2.450313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9871498) q[0];
sx q[0];
rz(-1.5342752) q[0];
sx q[0];
rz(-0.64572389) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5070314) q[2];
sx q[2];
rz(-0.87997961) q[2];
sx q[2];
rz(1.4379759) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9554841) q[1];
sx q[1];
rz(-2.3209474) q[1];
sx q[1];
rz(-0.11891831) q[1];
rz(-pi) q[2];
rz(0.78929928) q[3];
sx q[3];
rz(-1.7501481) q[3];
sx q[3];
rz(-1.8478249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1401691) q[2];
sx q[2];
rz(-0.836335) q[2];
sx q[2];
rz(2.8720065) q[2];
rz(-0.94830281) q[3];
sx q[3];
rz(-1.6092665) q[3];
sx q[3];
rz(1.74291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36169323) q[0];
sx q[0];
rz(-2.7044856) q[0];
sx q[0];
rz(-2.1345188) q[0];
rz(-2.7359447) q[1];
sx q[1];
rz(-2.5463153) q[1];
sx q[1];
rz(-0.55353037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29964511) q[0];
sx q[0];
rz(-1.3163693) q[0];
sx q[0];
rz(3.0823452) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78041665) q[2];
sx q[2];
rz(-2.4274106) q[2];
sx q[2];
rz(1.2250021) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.036052536) q[1];
sx q[1];
rz(-0.80784384) q[1];
sx q[1];
rz(-2.4921662) q[1];
x q[2];
rz(-3.0789525) q[3];
sx q[3];
rz(-0.38068145) q[3];
sx q[3];
rz(-0.21089396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3125399) q[2];
sx q[2];
rz(-2.0487787) q[2];
sx q[2];
rz(-1.9276169) q[2];
rz(2.35516) q[3];
sx q[3];
rz(-1.395547) q[3];
sx q[3];
rz(0.023155183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9455652) q[0];
sx q[0];
rz(-0.29226154) q[0];
sx q[0];
rz(-1.3319525) q[0];
rz(2.5281483) q[1];
sx q[1];
rz(-2.1211801) q[1];
sx q[1];
rz(2.6920998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025279538) q[0];
sx q[0];
rz(-0.49508245) q[0];
sx q[0];
rz(1.7465003) q[0];
rz(3.055213) q[2];
sx q[2];
rz(-1.9055718) q[2];
sx q[2];
rz(2.2074062) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66287884) q[1];
sx q[1];
rz(-2.3509074) q[1];
sx q[1];
rz(1.7639169) q[1];
rz(-pi) q[2];
rz(0.2576377) q[3];
sx q[3];
rz(-1.6272568) q[3];
sx q[3];
rz(-2.8650521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8248262) q[2];
sx q[2];
rz(-2.4591441) q[2];
sx q[2];
rz(2.5332434) q[2];
rz(1.9781205) q[3];
sx q[3];
rz(-1.3694265) q[3];
sx q[3];
rz(2.5782862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0332396) q[0];
sx q[0];
rz(-0.9110564) q[0];
sx q[0];
rz(0.60229993) q[0];
rz(0.96744084) q[1];
sx q[1];
rz(-0.93092218) q[1];
sx q[1];
rz(-1.1036576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7491584) q[0];
sx q[0];
rz(-1.9604248) q[0];
sx q[0];
rz(-2.1632739) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5561364) q[2];
sx q[2];
rz(-0.93428627) q[2];
sx q[2];
rz(2.9755862) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98075726) q[1];
sx q[1];
rz(-2.6904562) q[1];
sx q[1];
rz(1.5734929) q[1];
rz(-2.1488701) q[3];
sx q[3];
rz(-1.5141271) q[3];
sx q[3];
rz(-0.68777675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2063107) q[2];
sx q[2];
rz(-1.3056359) q[2];
sx q[2];
rz(2.804011) q[2];
rz(0.28389367) q[3];
sx q[3];
rz(-2.4604535) q[3];
sx q[3];
rz(-0.18686992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5934481) q[0];
sx q[0];
rz(-1.5080867) q[0];
sx q[0];
rz(-0.067807587) q[0];
rz(-1.9920805) q[1];
sx q[1];
rz(-1.5510473) q[1];
sx q[1];
rz(0.4745208) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10013469) q[0];
sx q[0];
rz(-0.23915072) q[0];
sx q[0];
rz(-1.5140223) q[0];
rz(-2.6340747) q[2];
sx q[2];
rz(-0.27406494) q[2];
sx q[2];
rz(0.12390359) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4700807) q[1];
sx q[1];
rz(-1.5776411) q[1];
sx q[1];
rz(-1.4840625) q[1];
rz(2.3274997) q[3];
sx q[3];
rz(-1.0225227) q[3];
sx q[3];
rz(0.89422885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6600251) q[2];
sx q[2];
rz(-0.93067545) q[2];
sx q[2];
rz(2.068326) q[2];
rz(-2.977071) q[3];
sx q[3];
rz(-1.8142895) q[3];
sx q[3];
rz(-2.8209414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5545223) q[0];
sx q[0];
rz(-1.5759435) q[0];
sx q[0];
rz(-1.6389621) q[0];
rz(1.5578237) q[1];
sx q[1];
rz(-2.0638034) q[1];
sx q[1];
rz(2.6149909) q[1];
rz(0.46109328) q[2];
sx q[2];
rz(-1.8965707) q[2];
sx q[2];
rz(-2.0988219) q[2];
rz(0.14231331) q[3];
sx q[3];
rz(-1.6518136) q[3];
sx q[3];
rz(1.470737) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
