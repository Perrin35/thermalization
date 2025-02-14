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
rz(1.9795228) q[0];
sx q[0];
rz(-2.906769) q[0];
sx q[0];
rz(1.809037) q[0];
rz(1.3136343) q[1];
sx q[1];
rz(-1.87275) q[1];
sx q[1];
rz(-2.3743) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3911661) q[0];
sx q[0];
rz(-0.39505233) q[0];
sx q[0];
rz(-3.1373298) q[0];
x q[1];
rz(0.69756223) q[2];
sx q[2];
rz(-2.5346906) q[2];
sx q[2];
rz(-1.4542435) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.060081944) q[1];
sx q[1];
rz(-0.36735555) q[1];
sx q[1];
rz(1.0698331) q[1];
rz(-pi) q[2];
rz(0.6758322) q[3];
sx q[3];
rz(-2.5321372) q[3];
sx q[3];
rz(-0.5389365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8511429) q[2];
sx q[2];
rz(-2.3252611) q[2];
sx q[2];
rz(-1.8905224) q[2];
rz(-0.81392455) q[3];
sx q[3];
rz(-1.1000752) q[3];
sx q[3];
rz(1.4936911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55452764) q[0];
sx q[0];
rz(-1.4545414) q[0];
sx q[0];
rz(2.5982507) q[0];
rz(0.73703611) q[1];
sx q[1];
rz(-2.4529424) q[1];
sx q[1];
rz(-0.063311689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86156002) q[0];
sx q[0];
rz(-0.82516042) q[0];
sx q[0];
rz(0.78365032) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28191992) q[2];
sx q[2];
rz(-1.8056944) q[2];
sx q[2];
rz(0.24215936) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.649401) q[1];
sx q[1];
rz(-2.1604558) q[1];
sx q[1];
rz(2.7604514) q[1];
rz(-pi) q[2];
rz(1.7821941) q[3];
sx q[3];
rz(-1.711004) q[3];
sx q[3];
rz(1.141215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43725499) q[2];
sx q[2];
rz(-2.0180118) q[2];
sx q[2];
rz(2.227318) q[2];
rz(2.8546913) q[3];
sx q[3];
rz(-2.5723781) q[3];
sx q[3];
rz(0.60681075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8617926) q[0];
sx q[0];
rz(-1.8656116) q[0];
sx q[0];
rz(1.3943425) q[0];
rz(-2.2718248) q[1];
sx q[1];
rz(-0.52647796) q[1];
sx q[1];
rz(2.8960752) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.190827) q[0];
sx q[0];
rz(-2.6689235) q[0];
sx q[0];
rz(-2.5373775) q[0];
rz(1.6310224) q[2];
sx q[2];
rz(-2.5909419) q[2];
sx q[2];
rz(-0.27397317) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80563021) q[1];
sx q[1];
rz(-1.9850678) q[1];
sx q[1];
rz(1.2719243) q[1];
x q[2];
rz(2.3342553) q[3];
sx q[3];
rz(-1.5669979) q[3];
sx q[3];
rz(-1.2646874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3123582) q[2];
sx q[2];
rz(-1.1623397) q[2];
sx q[2];
rz(-1.8355628) q[2];
rz(0.28451377) q[3];
sx q[3];
rz(-1.6407216) q[3];
sx q[3];
rz(1.3792926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6156085) q[0];
sx q[0];
rz(-2.5475976) q[0];
sx q[0];
rz(2.3748412) q[0];
rz(1.426567) q[1];
sx q[1];
rz(-1.7568935) q[1];
sx q[1];
rz(2.0791304) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3037231) q[0];
sx q[0];
rz(-0.024191054) q[0];
sx q[0];
rz(2.2538857) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97924657) q[2];
sx q[2];
rz(-1.4580168) q[2];
sx q[2];
rz(1.2119669) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6973138) q[1];
sx q[1];
rz(-2.0525888) q[1];
sx q[1];
rz(-2.8413111) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6312319) q[3];
sx q[3];
rz(-1.7157367) q[3];
sx q[3];
rz(1.0355024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1876424) q[2];
sx q[2];
rz(-3.1270449) q[2];
sx q[2];
rz(-3.0237831) q[2];
rz(-2.5976962) q[3];
sx q[3];
rz(-2.1102648) q[3];
sx q[3];
rz(0.73345524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2619005) q[0];
sx q[0];
rz(-1.6629135) q[0];
sx q[0];
rz(-0.53801584) q[0];
rz(0.063848786) q[1];
sx q[1];
rz(-1.0212746) q[1];
sx q[1];
rz(1.341238) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2573451) q[0];
sx q[0];
rz(-2.1187415) q[0];
sx q[0];
rz(-2.1764285) q[0];
rz(-pi) q[1];
rz(1.8939429) q[2];
sx q[2];
rz(-2.1003336) q[2];
sx q[2];
rz(0.33796899) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76439124) q[1];
sx q[1];
rz(-2.1294597) q[1];
sx q[1];
rz(-1.6425981) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43301626) q[3];
sx q[3];
rz(-1.6891494) q[3];
sx q[3];
rz(1.0326715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4591879) q[2];
sx q[2];
rz(-0.29067278) q[2];
sx q[2];
rz(0.9064557) q[2];
rz(1.9393548) q[3];
sx q[3];
rz(-1.5952483) q[3];
sx q[3];
rz(-1.7513587) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15659139) q[0];
sx q[0];
rz(-2.5522794) q[0];
sx q[0];
rz(-2.6949448) q[0];
rz(2.7062972) q[1];
sx q[1];
rz(-2.5050102) q[1];
sx q[1];
rz(-0.84904233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34693107) q[0];
sx q[0];
rz(-1.4022572) q[0];
sx q[0];
rz(-3.0656673) q[0];
x q[1];
rz(-2.8809636) q[2];
sx q[2];
rz(-0.5472551) q[2];
sx q[2];
rz(-1.7420285) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6852279) q[1];
sx q[1];
rz(-1.5816018) q[1];
sx q[1];
rz(-0.2266591) q[1];
x q[2];
rz(-2.7924323) q[3];
sx q[3];
rz(-1.4085692) q[3];
sx q[3];
rz(2.8256299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1242096) q[2];
sx q[2];
rz(-2.1653192) q[2];
sx q[2];
rz(1.3503831) q[2];
rz(-0.24556686) q[3];
sx q[3];
rz(-2.1731302) q[3];
sx q[3];
rz(2.811331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6733112) q[0];
sx q[0];
rz(-0.89235965) q[0];
sx q[0];
rz(-2.6924676) q[0];
rz(1.4324073) q[1];
sx q[1];
rz(-0.9811554) q[1];
sx q[1];
rz(-1.2264576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1311199) q[0];
sx q[0];
rz(-2.1852699) q[0];
sx q[0];
rz(0.09297405) q[0];
rz(-1.9513556) q[2];
sx q[2];
rz(-2.042843) q[2];
sx q[2];
rz(-2.7573331) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12237182) q[1];
sx q[1];
rz(-1.4667515) q[1];
sx q[1];
rz(-1.4833443) q[1];
rz(-pi) q[2];
rz(0.38568078) q[3];
sx q[3];
rz(-1.2880518) q[3];
sx q[3];
rz(-1.4644506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8760828) q[2];
sx q[2];
rz(-2.2777057) q[2];
sx q[2];
rz(-0.49501219) q[2];
rz(2.072295) q[3];
sx q[3];
rz(-2.4102231) q[3];
sx q[3];
rz(2.4643912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4070364) q[0];
sx q[0];
rz(-0.41663909) q[0];
sx q[0];
rz(-0.46105841) q[0];
rz(-0.1289434) q[1];
sx q[1];
rz(-2.415633) q[1];
sx q[1];
rz(0.9022378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8811097) q[0];
sx q[0];
rz(-0.04180464) q[0];
sx q[0];
rz(2.9737472) q[0];
x q[1];
rz(-1.938615) q[2];
sx q[2];
rz(-1.6253643) q[2];
sx q[2];
rz(-0.56625783) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5758491) q[1];
sx q[1];
rz(-1.1593124) q[1];
sx q[1];
rz(-0.38994148) q[1];
rz(-1.447313) q[3];
sx q[3];
rz(-2.4816645) q[3];
sx q[3];
rz(0.29327682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43877131) q[2];
sx q[2];
rz(-2.2192025) q[2];
sx q[2];
rz(1.927446) q[2];
rz(-0.37902784) q[3];
sx q[3];
rz(-1.4641848) q[3];
sx q[3];
rz(-1.7927143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2269065) q[0];
sx q[0];
rz(-1.5808957) q[0];
sx q[0];
rz(-3.1353986) q[0];
rz(2.6365623) q[1];
sx q[1];
rz(-2.6513702) q[1];
sx q[1];
rz(-0.44713155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8861388) q[0];
sx q[0];
rz(-2.2199174) q[0];
sx q[0];
rz(-2.4711424) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6225409) q[2];
sx q[2];
rz(-1.468892) q[2];
sx q[2];
rz(-0.48619147) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3382692) q[1];
sx q[1];
rz(-1.4191966) q[1];
sx q[1];
rz(-3.1123509) q[1];
x q[2];
rz(-1.379844) q[3];
sx q[3];
rz(-1.1292283) q[3];
sx q[3];
rz(0.34363765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7827683) q[2];
sx q[2];
rz(-2.0147822) q[2];
sx q[2];
rz(-0.44754851) q[2];
rz(1.1131845) q[3];
sx q[3];
rz(-1.0380849) q[3];
sx q[3];
rz(-0.32570496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.177445) q[0];
sx q[0];
rz(-3.0772458) q[0];
sx q[0];
rz(2.294975) q[0];
rz(2.8586491) q[1];
sx q[1];
rz(-1.6694262) q[1];
sx q[1];
rz(-0.67323321) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1812912) q[0];
sx q[0];
rz(-0.01894572) q[0];
sx q[0];
rz(0.51769729) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77263562) q[2];
sx q[2];
rz(-1.7353206) q[2];
sx q[2];
rz(-2.5009837) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9615508) q[1];
sx q[1];
rz(-0.93856877) q[1];
sx q[1];
rz(0.2747196) q[1];
rz(-pi) q[2];
rz(0.85906927) q[3];
sx q[3];
rz(-2.7112656) q[3];
sx q[3];
rz(2.8266065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8474569) q[2];
sx q[2];
rz(-0.86926502) q[2];
sx q[2];
rz(2.4693303) q[2];
rz(-1.2761891) q[3];
sx q[3];
rz(-1.0381235) q[3];
sx q[3];
rz(-0.41260317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9857674) q[0];
sx q[0];
rz(-1.1341996) q[0];
sx q[0];
rz(-1.2962935) q[0];
rz(0.014421944) q[1];
sx q[1];
rz(-1.8564693) q[1];
sx q[1];
rz(1.2546702) q[1];
rz(2.4449271) q[2];
sx q[2];
rz(-1.9370542) q[2];
sx q[2];
rz(2.5123971) q[2];
rz(0.23019595) q[3];
sx q[3];
rz(-1.7741813) q[3];
sx q[3];
rz(0.81188191) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
