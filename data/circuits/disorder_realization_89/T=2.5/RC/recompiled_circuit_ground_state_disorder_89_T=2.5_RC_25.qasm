OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1998347) q[0];
sx q[0];
rz(-0.69154843) q[0];
sx q[0];
rz(-0.77342311) q[0];
rz(3.050488) q[1];
sx q[1];
rz(-0.79371047) q[1];
sx q[1];
rz(1.8082126) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1379194) q[0];
sx q[0];
rz(-0.9124476) q[0];
sx q[0];
rz(-1.5594894) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23024586) q[2];
sx q[2];
rz(-1.0827999) q[2];
sx q[2];
rz(-2.5813318) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62480085) q[1];
sx q[1];
rz(-1.1784823) q[1];
sx q[1];
rz(-0.90992163) q[1];
rz(1.8190906) q[3];
sx q[3];
rz(-0.069337519) q[3];
sx q[3];
rz(1.2420036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.29204631) q[2];
sx q[2];
rz(-1.8201733) q[2];
sx q[2];
rz(0.41479659) q[2];
rz(-1.241812) q[3];
sx q[3];
rz(-1.1153699) q[3];
sx q[3];
rz(-2.950086) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6280129) q[0];
sx q[0];
rz(-0.11592557) q[0];
sx q[0];
rz(2.7017748) q[0];
rz(3.0385333) q[1];
sx q[1];
rz(-0.28918806) q[1];
sx q[1];
rz(1.1691079) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98840731) q[0];
sx q[0];
rz(-0.64278945) q[0];
sx q[0];
rz(1.8039186) q[0];
rz(-pi) q[1];
x q[1];
rz(1.493426) q[2];
sx q[2];
rz(-2.5499509) q[2];
sx q[2];
rz(-0.40654069) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4354187) q[1];
sx q[1];
rz(-1.7507554) q[1];
sx q[1];
rz(-2.2246996) q[1];
rz(-2.2803047) q[3];
sx q[3];
rz(-1.420974) q[3];
sx q[3];
rz(3.0784208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.64572) q[2];
sx q[2];
rz(-1.7503259) q[2];
sx q[2];
rz(0.24620852) q[2];
rz(1.5496893) q[3];
sx q[3];
rz(-2.527745) q[3];
sx q[3];
rz(-1.7483819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98591268) q[0];
sx q[0];
rz(-2.0014626) q[0];
sx q[0];
rz(2.929856) q[0];
rz(3.1302997) q[1];
sx q[1];
rz(-2.3432422) q[1];
sx q[1];
rz(-1.6976154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.885545) q[0];
sx q[0];
rz(-0.035158947) q[0];
sx q[0];
rz(-0.0020744046) q[0];
rz(1.774774) q[2];
sx q[2];
rz(-0.24572554) q[2];
sx q[2];
rz(-1.3947006) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1423751) q[1];
sx q[1];
rz(-2.1793167) q[1];
sx q[1];
rz(-1.1097679) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6614603) q[3];
sx q[3];
rz(-1.1216402) q[3];
sx q[3];
rz(-0.36270519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4766562) q[2];
sx q[2];
rz(-1.0756476) q[2];
sx q[2];
rz(-2.5962489) q[2];
rz(2.2319345) q[3];
sx q[3];
rz(-0.46234149) q[3];
sx q[3];
rz(1.9285412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42105168) q[0];
sx q[0];
rz(-0.67890972) q[0];
sx q[0];
rz(0.82756591) q[0];
rz(1.1830117) q[1];
sx q[1];
rz(-1.2748101) q[1];
sx q[1];
rz(1.8756728) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0166548) q[0];
sx q[0];
rz(-1.6724574) q[0];
sx q[0];
rz(-0.09135017) q[0];
rz(-1.5180196) q[2];
sx q[2];
rz(-0.20773331) q[2];
sx q[2];
rz(-0.29333255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.75531653) q[1];
sx q[1];
rz(-2.1960253) q[1];
sx q[1];
rz(-3.1023854) q[1];
x q[2];
rz(-1.3807348) q[3];
sx q[3];
rz(-1.7281088) q[3];
sx q[3];
rz(0.33657956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8852692) q[2];
sx q[2];
rz(-1.4726535) q[2];
sx q[2];
rz(-0.081776865) q[2];
rz(2.8335588) q[3];
sx q[3];
rz(-2.6576198) q[3];
sx q[3];
rz(-1.5380194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772188) q[0];
sx q[0];
rz(-2.866221) q[0];
sx q[0];
rz(-3.0710869) q[0];
rz(2.8071857) q[1];
sx q[1];
rz(-0.53135482) q[1];
sx q[1];
rz(2.6883584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8716259) q[0];
sx q[0];
rz(-2.0701417) q[0];
sx q[0];
rz(0.64887394) q[0];
x q[1];
rz(-1.0595991) q[2];
sx q[2];
rz(-0.44253749) q[2];
sx q[2];
rz(2.0997186) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2259648) q[1];
sx q[1];
rz(-2.8702822) q[1];
sx q[1];
rz(-2.7167999) q[1];
rz(-pi) q[2];
rz(-2.2599561) q[3];
sx q[3];
rz(-0.83587468) q[3];
sx q[3];
rz(2.3595032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8411023) q[2];
sx q[2];
rz(-1.5056242) q[2];
sx q[2];
rz(-0.21031586) q[2];
rz(-2.7503843) q[3];
sx q[3];
rz(-0.20520964) q[3];
sx q[3];
rz(2.6989663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.6415569) q[0];
sx q[0];
rz(-2.0099202) q[0];
sx q[0];
rz(2.9697707) q[0];
rz(1.7956644) q[1];
sx q[1];
rz(-0.95564061) q[1];
sx q[1];
rz(-1.9926434) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7846735) q[0];
sx q[0];
rz(-1.5789512) q[0];
sx q[0];
rz(3.0300333) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8549598) q[2];
sx q[2];
rz(-1.3475085) q[2];
sx q[2];
rz(-2.5866825) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.51275245) q[1];
sx q[1];
rz(-1.0003371) q[1];
sx q[1];
rz(2.9243584) q[1];
rz(-2.9385185) q[3];
sx q[3];
rz(-2.0690898) q[3];
sx q[3];
rz(-0.6607252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.020784) q[2];
sx q[2];
rz(-0.9641996) q[2];
sx q[2];
rz(2.7210893) q[2];
rz(1.6546107) q[3];
sx q[3];
rz(-1.9872811) q[3];
sx q[3];
rz(-1.1544863) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3625951) q[0];
sx q[0];
rz(-1.2263466) q[0];
sx q[0];
rz(2.763789) q[0];
rz(-1.3899577) q[1];
sx q[1];
rz(-1.4327587) q[1];
sx q[1];
rz(-0.43209824) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8209131) q[0];
sx q[0];
rz(-0.28962505) q[0];
sx q[0];
rz(-2.4885213) q[0];
rz(-0.092707002) q[2];
sx q[2];
rz(-1.5376629) q[2];
sx q[2];
rz(-2.9404145) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3258265) q[1];
sx q[1];
rz(-1.5349421) q[1];
sx q[1];
rz(0.46200606) q[1];
x q[2];
rz(-1.0749519) q[3];
sx q[3];
rz(-0.73710873) q[3];
sx q[3];
rz(-2.6219311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7695339) q[2];
sx q[2];
rz(-2.5199315) q[2];
sx q[2];
rz(2.9080234) q[2];
rz(1.6893049) q[3];
sx q[3];
rz(-1.4315616) q[3];
sx q[3];
rz(-1.5610032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772783) q[0];
sx q[0];
rz(-1.0294585) q[0];
sx q[0];
rz(-2.8218063) q[0];
rz(0.86209595) q[1];
sx q[1];
rz(-1.8902706) q[1];
sx q[1];
rz(1.7822942) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089238374) q[0];
sx q[0];
rz(-2.1349094) q[0];
sx q[0];
rz(2.4648239) q[0];
rz(3.0002757) q[2];
sx q[2];
rz(-1.3719146) q[2];
sx q[2];
rz(2.3945253) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8492125) q[1];
sx q[1];
rz(-1.74382) q[1];
sx q[1];
rz(-0.72489691) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8639761) q[3];
sx q[3];
rz(-0.95562387) q[3];
sx q[3];
rz(0.53078101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84997815) q[2];
sx q[2];
rz(-2.7225967) q[2];
sx q[2];
rz(-0.46763793) q[2];
rz(-2.6575798) q[3];
sx q[3];
rz(-1.7734807) q[3];
sx q[3];
rz(-1.2776933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9676232) q[0];
sx q[0];
rz(-1.6521709) q[0];
sx q[0];
rz(2.0066579) q[0];
rz(-0.060215503) q[1];
sx q[1];
rz(-2.4048012) q[1];
sx q[1];
rz(-1.6103475) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0031207) q[0];
sx q[0];
rz(-1.1116189) q[0];
sx q[0];
rz(-0.91158406) q[0];
rz(-2.4545099) q[2];
sx q[2];
rz(-1.9959046) q[2];
sx q[2];
rz(1.2883505) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5694014) q[1];
sx q[1];
rz(-1.4735771) q[1];
sx q[1];
rz(-0.26406322) q[1];
rz(-pi) q[2];
rz(1.5315191) q[3];
sx q[3];
rz(-2.352592) q[3];
sx q[3];
rz(-1.6226744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.81426364) q[2];
sx q[2];
rz(-1.2171823) q[2];
sx q[2];
rz(-0.72861707) q[2];
rz(-2.9260351) q[3];
sx q[3];
rz(-1.39648) q[3];
sx q[3];
rz(2.047915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88290596) q[0];
sx q[0];
rz(-3.1104493) q[0];
sx q[0];
rz(0.31797847) q[0];
rz(1.7440354) q[1];
sx q[1];
rz(-1.8122858) q[1];
sx q[1];
rz(2.8584282) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69071373) q[0];
sx q[0];
rz(-2.9895999) q[0];
sx q[0];
rz(0.99894036) q[0];
rz(1.4941169) q[2];
sx q[2];
rz(-0.81816219) q[2];
sx q[2];
rz(0.61329182) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6826815) q[1];
sx q[1];
rz(-1.1106756) q[1];
sx q[1];
rz(-2.3691872) q[1];
rz(1.3628325) q[3];
sx q[3];
rz(-1.4803518) q[3];
sx q[3];
rz(-1.2048723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91750034) q[2];
sx q[2];
rz(-1.6250236) q[2];
sx q[2];
rz(-0.089281233) q[2];
rz(1.3034405) q[3];
sx q[3];
rz(-0.83804122) q[3];
sx q[3];
rz(-2.1681521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3598809) q[0];
sx q[0];
rz(-2.4926873) q[0];
sx q[0];
rz(0.98095184) q[0];
rz(0.5858865) q[1];
sx q[1];
rz(-1.2955019) q[1];
sx q[1];
rz(-1.6265709) q[1];
rz(-0.21296756) q[2];
sx q[2];
rz(-1.3942547) q[2];
sx q[2];
rz(-0.19993776) q[2];
rz(-2.8845434) q[3];
sx q[3];
rz(-1.201073) q[3];
sx q[3];
rz(1.4504688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
