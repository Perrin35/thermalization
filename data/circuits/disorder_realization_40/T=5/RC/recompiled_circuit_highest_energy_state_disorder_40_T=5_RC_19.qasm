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
rz(1.7215913) q[0];
sx q[0];
rz(-2.1821332) q[0];
sx q[0];
rz(1.0406915) q[0];
rz(-0.33848441) q[1];
sx q[1];
rz(-1.4893293) q[1];
sx q[1];
rz(-1.4298061) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0105699) q[0];
sx q[0];
rz(-1.4522465) q[0];
sx q[0];
rz(0.27310102) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0645039) q[2];
sx q[2];
rz(-2.4790451) q[2];
sx q[2];
rz(2.8665989) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0104645) q[1];
sx q[1];
rz(-0.87112037) q[1];
sx q[1];
rz(-1.9424428) q[1];
x q[2];
rz(2.5761371) q[3];
sx q[3];
rz(-1.3420336) q[3];
sx q[3];
rz(2.8086587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0040943) q[2];
sx q[2];
rz(-1.4278922) q[2];
sx q[2];
rz(0.049169866) q[2];
rz(2.5947425) q[3];
sx q[3];
rz(-2.2873736) q[3];
sx q[3];
rz(-1.2949519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8787254) q[0];
sx q[0];
rz(-2.290949) q[0];
sx q[0];
rz(-3.1021297) q[0];
rz(-2.4354758) q[1];
sx q[1];
rz(-1.1742274) q[1];
sx q[1];
rz(-1.2605234) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8055819) q[0];
sx q[0];
rz(-1.1534497) q[0];
sx q[0];
rz(-0.063800617) q[0];
rz(-pi) q[1];
rz(-2.4751365) q[2];
sx q[2];
rz(-1.1315695) q[2];
sx q[2];
rz(-2.0298438) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.305734) q[1];
sx q[1];
rz(-1.5286338) q[1];
sx q[1];
rz(-0.034138676) q[1];
x q[2];
rz(-2.8144263) q[3];
sx q[3];
rz(-0.20450243) q[3];
sx q[3];
rz(-0.12302264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4827106) q[2];
sx q[2];
rz(-2.1363246) q[2];
sx q[2];
rz(-0.91747326) q[2];
rz(0.42258036) q[3];
sx q[3];
rz(-1.8924507) q[3];
sx q[3];
rz(2.1036072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68040401) q[0];
sx q[0];
rz(-0.55484158) q[0];
sx q[0];
rz(2.19221) q[0];
rz(1.2481015) q[1];
sx q[1];
rz(-2.5578942) q[1];
sx q[1];
rz(1.5138352) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4405524) q[0];
sx q[0];
rz(-1.3258385) q[0];
sx q[0];
rz(-2.7666758) q[0];
x q[1];
rz(-2.0545839) q[2];
sx q[2];
rz(-1.3788333) q[2];
sx q[2];
rz(1.0707945) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.48958923) q[1];
sx q[1];
rz(-1.6449274) q[1];
sx q[1];
rz(-0.28848047) q[1];
rz(-pi) q[2];
rz(1.85905) q[3];
sx q[3];
rz(-1.6950399) q[3];
sx q[3];
rz(-1.8831203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44024399) q[2];
sx q[2];
rz(-1.1260208) q[2];
sx q[2];
rz(-1.6181642) q[2];
rz(-0.31420389) q[3];
sx q[3];
rz(-2.115963) q[3];
sx q[3];
rz(-1.6395125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0358129) q[0];
sx q[0];
rz(-0.19136763) q[0];
sx q[0];
rz(-2.9845003) q[0];
rz(-3.0063903) q[1];
sx q[1];
rz(-0.78510761) q[1];
sx q[1];
rz(-0.33321998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5409637) q[0];
sx q[0];
rz(-2.0218533) q[0];
sx q[0];
rz(-0.55229295) q[0];
rz(1.7945788) q[2];
sx q[2];
rz(-0.74088851) q[2];
sx q[2];
rz(0.2303094) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95851719) q[1];
sx q[1];
rz(-2.6974899) q[1];
sx q[1];
rz(-0.91981319) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2832421) q[3];
sx q[3];
rz(-1.3305404) q[3];
sx q[3];
rz(2.1879856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.80804431) q[2];
sx q[2];
rz(-0.6490038) q[2];
sx q[2];
rz(-1.6928847) q[2];
rz(-2.4088805) q[3];
sx q[3];
rz(-1.9219739) q[3];
sx q[3];
rz(-1.5318416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6919959) q[0];
sx q[0];
rz(-1.6549598) q[0];
sx q[0];
rz(2.9100371) q[0];
rz(-0.72136503) q[1];
sx q[1];
rz(-0.8911348) q[1];
sx q[1];
rz(0.84842938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68616406) q[0];
sx q[0];
rz(-1.9047184) q[0];
sx q[0];
rz(-1.0592106) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5987342) q[2];
sx q[2];
rz(-1.9221537) q[2];
sx q[2];
rz(1.5962102) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0368122) q[1];
sx q[1];
rz(-1.9580578) q[1];
sx q[1];
rz(-1.6464064) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23691688) q[3];
sx q[3];
rz(-1.6730783) q[3];
sx q[3];
rz(-1.8558242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.007175) q[2];
sx q[2];
rz(-2.538372) q[2];
sx q[2];
rz(1.2009386) q[2];
rz(2.4522771) q[3];
sx q[3];
rz(-1.9715693) q[3];
sx q[3];
rz(-2.5758666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1597964) q[0];
sx q[0];
rz(-1.5317651) q[0];
sx q[0];
rz(-1.7933886) q[0];
rz(2.4234405) q[1];
sx q[1];
rz(-0.57128692) q[1];
sx q[1];
rz(-1.2917554) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6245218) q[0];
sx q[0];
rz(-0.41247955) q[0];
sx q[0];
rz(-0.58858354) q[0];
rz(-pi) q[1];
rz(-0.37448762) q[2];
sx q[2];
rz(-2.1364223) q[2];
sx q[2];
rz(-0.35978064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8706269) q[1];
sx q[1];
rz(-2.9521204) q[1];
sx q[1];
rz(0.47196526) q[1];
rz(-pi) q[2];
x q[2];
rz(0.060543493) q[3];
sx q[3];
rz(-0.39602867) q[3];
sx q[3];
rz(-0.76914364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1527839) q[2];
sx q[2];
rz(-0.96502105) q[2];
sx q[2];
rz(-2.7394845) q[2];
rz(-2.1357338) q[3];
sx q[3];
rz(-2.918225) q[3];
sx q[3];
rz(-1.7042101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.724029) q[0];
sx q[0];
rz(-0.14823866) q[0];
sx q[0];
rz(-1.8120026) q[0];
rz(-1.3279462) q[1];
sx q[1];
rz(-1.4168394) q[1];
sx q[1];
rz(0.45164576) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490135) q[0];
sx q[0];
rz(-2.012741) q[0];
sx q[0];
rz(1.4699292) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7164104) q[2];
sx q[2];
rz(-1.4431134) q[2];
sx q[2];
rz(0.38056254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0473808) q[1];
sx q[1];
rz(-1.3806731) q[1];
sx q[1];
rz(-3.090108) q[1];
x q[2];
rz(-0.65589738) q[3];
sx q[3];
rz(-2.128388) q[3];
sx q[3];
rz(-2.2849071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9391276) q[2];
sx q[2];
rz(-1.7697325) q[2];
sx q[2];
rz(0.28787127) q[2];
rz(-1.6599844) q[3];
sx q[3];
rz(-2.2899254) q[3];
sx q[3];
rz(1.456858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57956368) q[0];
sx q[0];
rz(-0.96051878) q[0];
sx q[0];
rz(-0.55291837) q[0];
rz(1.8909854) q[1];
sx q[1];
rz(-0.41025531) q[1];
sx q[1];
rz(1.3804573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4794248) q[0];
sx q[0];
rz(-1.2081376) q[0];
sx q[0];
rz(2.5890686) q[0];
x q[1];
rz(-0.98357551) q[2];
sx q[2];
rz(-1.4693575) q[2];
sx q[2];
rz(-1.3124397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1565754) q[1];
sx q[1];
rz(-0.24057287) q[1];
sx q[1];
rz(-0.17372082) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7431511) q[3];
sx q[3];
rz(-0.77756778) q[3];
sx q[3];
rz(0.48157495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55014253) q[2];
sx q[2];
rz(-1.3478792) q[2];
sx q[2];
rz(2.4264917) q[2];
rz(-3.0969369) q[3];
sx q[3];
rz(-2.5584593) q[3];
sx q[3];
rz(0.62838069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1735246) q[0];
sx q[0];
rz(-0.76924291) q[0];
sx q[0];
rz(0.64088696) q[0];
rz(1.4961273) q[1];
sx q[1];
rz(-2.75664) q[1];
sx q[1];
rz(0.65151757) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1524044) q[0];
sx q[0];
rz(-1.8274283) q[0];
sx q[0];
rz(2.2874831) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7213621) q[2];
sx q[2];
rz(-2.6131008) q[2];
sx q[2];
rz(-0.53208379) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1478737) q[1];
sx q[1];
rz(-2.2145956) q[1];
sx q[1];
rz(2.6753475) q[1];
rz(1.9189758) q[3];
sx q[3];
rz(-1.7507675) q[3];
sx q[3];
rz(0.5157541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7351825) q[2];
sx q[2];
rz(-2.0969022) q[2];
sx q[2];
rz(2.6835105) q[2];
rz(0.011890751) q[3];
sx q[3];
rz(-1.0606822) q[3];
sx q[3];
rz(-0.021765821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40282014) q[0];
sx q[0];
rz(-2.9467376) q[0];
sx q[0];
rz(-0.028976945) q[0];
rz(-0.92652357) q[1];
sx q[1];
rz(-1.4484826) q[1];
sx q[1];
rz(0.40625939) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6299141) q[0];
sx q[0];
rz(-1.2907012) q[0];
sx q[0];
rz(1.5797432) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.026097) q[2];
sx q[2];
rz(-2.4433854) q[2];
sx q[2];
rz(0.9309665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.528842) q[1];
sx q[1];
rz(-1.539673) q[1];
sx q[1];
rz(1.2791388) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3876312) q[3];
sx q[3];
rz(-0.9545325) q[3];
sx q[3];
rz(-2.1631416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77068344) q[2];
sx q[2];
rz(-1.6434881) q[2];
sx q[2];
rz(2.5468199) q[2];
rz(-2.4316783) q[3];
sx q[3];
rz(-0.9570595) q[3];
sx q[3];
rz(2.1962568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5668673) q[0];
sx q[0];
rz(-2.18676) q[0];
sx q[0];
rz(-2.2796897) q[0];
rz(-2.8361539) q[1];
sx q[1];
rz(-1.8157235) q[1];
sx q[1];
rz(-2.9095412) q[1];
rz(2.7965056) q[2];
sx q[2];
rz(-1.9660334) q[2];
sx q[2];
rz(2.6031969) q[2];
rz(-2.6628982) q[3];
sx q[3];
rz(-2.6810418) q[3];
sx q[3];
rz(-2.2543805) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
