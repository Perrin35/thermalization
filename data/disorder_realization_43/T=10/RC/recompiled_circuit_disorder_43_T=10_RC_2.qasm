OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2551978) q[0];
sx q[0];
rz(-1.891341) q[0];
sx q[0];
rz(1.3347081) q[0];
rz(2.788738) q[1];
sx q[1];
rz(-2.9810413) q[1];
sx q[1];
rz(-0.97595739) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2655576) q[0];
sx q[0];
rz(-1.3296488) q[0];
sx q[0];
rz(-2.5493456) q[0];
rz(-pi) q[1];
rz(2.8085254) q[2];
sx q[2];
rz(-2.1085848) q[2];
sx q[2];
rz(-1.2531467) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.576697) q[1];
sx q[1];
rz(-1.4370059) q[1];
sx q[1];
rz(0.96058515) q[1];
rz(2.7482412) q[3];
sx q[3];
rz(-0.8439807) q[3];
sx q[3];
rz(-0.8918744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7261937) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(-2.3577918) q[2];
rz(0.33256724) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48288229) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-2.9630307) q[0];
rz(1.3372955) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(-0.006342412) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9316677) q[0];
sx q[0];
rz(-1.4405662) q[0];
sx q[0];
rz(1.3757214) q[0];
rz(-pi) q[1];
rz(0.86654051) q[2];
sx q[2];
rz(-0.73359493) q[2];
sx q[2];
rz(-1.7689592) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9176863) q[1];
sx q[1];
rz(-0.31653857) q[1];
sx q[1];
rz(-0.19885893) q[1];
x q[2];
rz(0.56224058) q[3];
sx q[3];
rz(-1.3573109) q[3];
sx q[3];
rz(1.436304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.94770849) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(0.87810278) q[2];
rz(0.86205035) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5851615) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(-0.01097824) q[0];
rz(0.36704656) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(3.045851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84155267) q[0];
sx q[0];
rz(-1.3944355) q[0];
sx q[0];
rz(-2.9862613) q[0];
rz(2.4143366) q[2];
sx q[2];
rz(-0.63823344) q[2];
sx q[2];
rz(0.22932316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1167736) q[1];
sx q[1];
rz(-1.0305335) q[1];
sx q[1];
rz(0.54622548) q[1];
rz(-pi) q[2];
rz(0.5511958) q[3];
sx q[3];
rz(-1.6008198) q[3];
sx q[3];
rz(2.4195645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0720955) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(0.83479184) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(-2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.9469706) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(-0.34657493) q[0];
rz(0.52571458) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(2.1077164) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82499408) q[0];
sx q[0];
rz(-0.98485095) q[0];
sx q[0];
rz(-0.98866776) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1448574) q[2];
sx q[2];
rz(-2.0248932) q[2];
sx q[2];
rz(3.0388289) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32767195) q[1];
sx q[1];
rz(-1.76183) q[1];
sx q[1];
rz(-3.0753067) q[1];
rz(-pi) q[2];
rz(-1.8978118) q[3];
sx q[3];
rz(-1.838284) q[3];
sx q[3];
rz(-0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(1.3467849) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088257) q[0];
sx q[0];
rz(-0.79341745) q[0];
sx q[0];
rz(2.72686) q[0];
rz(-1.746009) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(-0.57410747) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60526472) q[0];
sx q[0];
rz(-2.5308373) q[0];
sx q[0];
rz(1.7954134) q[0];
rz(-1.5933883) q[2];
sx q[2];
rz(-0.85561692) q[2];
sx q[2];
rz(-1.1720282) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.638006) q[1];
sx q[1];
rz(-1.9203016) q[1];
sx q[1];
rz(-1.867678) q[1];
x q[2];
rz(0.90548924) q[3];
sx q[3];
rz(-1.724616) q[3];
sx q[3];
rz(-2.820435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2888912) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(-1.5138907) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0742652) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(2.4434027) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(0.73227698) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66757827) q[0];
sx q[0];
rz(-1.7828373) q[0];
sx q[0];
rz(1.408512) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4562624) q[2];
sx q[2];
rz(-0.88980674) q[2];
sx q[2];
rz(0.53390098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60019833) q[1];
sx q[1];
rz(-1.3458034) q[1];
sx q[1];
rz(0.66585983) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5703339) q[3];
sx q[3];
rz(-1.3864281) q[3];
sx q[3];
rz(-1.0384699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7531062) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(-0.30900624) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(-1.1674315) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56918615) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(-2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(2.3715473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2419287) q[0];
sx q[0];
rz(-2.1673492) q[0];
sx q[0];
rz(1.7799737) q[0];
rz(-pi) q[1];
rz(-2.7355843) q[2];
sx q[2];
rz(-1.0205262) q[2];
sx q[2];
rz(-1.6183491) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8416482) q[1];
sx q[1];
rz(-2.3023459) q[1];
sx q[1];
rz(-2.4323835) q[1];
x q[2];
rz(-2.5823334) q[3];
sx q[3];
rz(-0.57563215) q[3];
sx q[3];
rz(-1.5534793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4975171) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(-1.7049568) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62676936) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(-0.24169895) q[0];
rz(2.4027951) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-2.7752005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05224932) q[0];
sx q[0];
rz(-2.0848668) q[0];
sx q[0];
rz(-2.423717) q[0];
x q[1];
rz(2.512152) q[2];
sx q[2];
rz(-1.2120314) q[2];
sx q[2];
rz(1.7164873) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3185127) q[1];
sx q[1];
rz(-1.7365343) q[1];
sx q[1];
rz(-2.4351099) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4922769) q[3];
sx q[3];
rz(-2.588387) q[3];
sx q[3];
rz(1.5403403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1239132) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(-2.8477342) q[2];
rz(0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(-0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095187) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(0.88395399) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(0.79137897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.489483) q[0];
sx q[0];
rz(-1.8882505) q[0];
sx q[0];
rz(0.86274685) q[0];
x q[1];
rz(-0.80957885) q[2];
sx q[2];
rz(-0.62925816) q[2];
sx q[2];
rz(-1.4776243) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3190805) q[1];
sx q[1];
rz(-1.1573536) q[1];
sx q[1];
rz(-2.6507069) q[1];
x q[2];
rz(1.8459686) q[3];
sx q[3];
rz(-0.55579805) q[3];
sx q[3];
rz(2.9142227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.066594921) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(1.0212612) q[2];
rz(-2.7630473) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(0.59797257) q[3];
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
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026697712) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(-0.60780203) q[0];
rz(-2.9027477) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.7609319) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35683435) q[0];
sx q[0];
rz(-2.1958302) q[0];
sx q[0];
rz(2.9805095) q[0];
x q[1];
rz(-2.0062749) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(-1.6053111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8163562) q[1];
sx q[1];
rz(-1.9776077) q[1];
sx q[1];
rz(2.7951294) q[1];
x q[2];
rz(2.4139666) q[3];
sx q[3];
rz(-0.25087038) q[3];
sx q[3];
rz(-0.48654702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(2.805368) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0040141) q[0];
sx q[0];
rz(-1.4083569) q[0];
sx q[0];
rz(1.2271723) q[0];
rz(2.5972988) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(0.82472807) q[2];
sx q[2];
rz(-1.3168954) q[2];
sx q[2];
rz(-1.6864824) q[2];
rz(1.929677) q[3];
sx q[3];
rz(-2.5419895) q[3];
sx q[3];
rz(0.06751577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
