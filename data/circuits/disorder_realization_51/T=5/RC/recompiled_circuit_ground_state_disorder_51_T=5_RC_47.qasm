OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.12282523) q[0];
sx q[0];
rz(3.1376165) q[0];
sx q[0];
rz(9.3038179) q[0];
rz(0.51789415) q[1];
sx q[1];
rz(2.8105812) q[1];
sx q[1];
rz(10.612648) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56138229) q[0];
sx q[0];
rz(-2.8299709) q[0];
sx q[0];
rz(0.28579692) q[0];
rz(0.14622525) q[2];
sx q[2];
rz(-1.0492965) q[2];
sx q[2];
rz(-2.5086049) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.814683) q[1];
sx q[1];
rz(-1.9724477) q[1];
sx q[1];
rz(-0.99237751) q[1];
x q[2];
rz(0.63159794) q[3];
sx q[3];
rz(-1.609841) q[3];
sx q[3];
rz(1.6717048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0821705) q[2];
sx q[2];
rz(-0.94258451) q[2];
sx q[2];
rz(1.7131294) q[2];
rz(1.7510471) q[3];
sx q[3];
rz(-2.3640552) q[3];
sx q[3];
rz(2.9318504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9222337) q[0];
sx q[0];
rz(-1.5809504) q[0];
sx q[0];
rz(1.6844164) q[0];
rz(-0.68151418) q[1];
sx q[1];
rz(-2.4512955) q[1];
sx q[1];
rz(-0.95726475) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0326704) q[0];
sx q[0];
rz(-1.3650948) q[0];
sx q[0];
rz(2.6673596) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8419015) q[2];
sx q[2];
rz(-2.3699591) q[2];
sx q[2];
rz(-2.6234805) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73814002) q[1];
sx q[1];
rz(-2.7305718) q[1];
sx q[1];
rz(-1.6938126) q[1];
rz(-pi) q[2];
rz(1.4232457) q[3];
sx q[3];
rz(-0.06961623) q[3];
sx q[3];
rz(0.92629647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7822781) q[2];
sx q[2];
rz(-1.0474397) q[2];
sx q[2];
rz(2.2347343) q[2];
rz(-2.4685229) q[3];
sx q[3];
rz(-0.38452092) q[3];
sx q[3];
rz(-1.2955906) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4562255) q[0];
sx q[0];
rz(-0.65995589) q[0];
sx q[0];
rz(-0.54687706) q[0];
rz(-1.3705672) q[1];
sx q[1];
rz(-1.1616881) q[1];
sx q[1];
rz(-0.10017698) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4456383) q[0];
sx q[0];
rz(-1.4379408) q[0];
sx q[0];
rz(-1.3039051) q[0];
x q[1];
rz(-2.6351026) q[2];
sx q[2];
rz(-1.2236986) q[2];
sx q[2];
rz(-0.60817761) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.924739) q[1];
sx q[1];
rz(-1.3432055) q[1];
sx q[1];
rz(-2.4000969) q[1];
x q[2];
rz(1.1421575) q[3];
sx q[3];
rz(-1.996962) q[3];
sx q[3];
rz(1.4166735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3935516) q[2];
sx q[2];
rz(-1.9358205) q[2];
sx q[2];
rz(2.098341) q[2];
rz(0.68759632) q[3];
sx q[3];
rz(-2.6396773) q[3];
sx q[3];
rz(-0.25575328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2954751) q[0];
sx q[0];
rz(-2.6961374) q[0];
sx q[0];
rz(-1.0365781) q[0];
rz(-1.2713894) q[1];
sx q[1];
rz(-2.3174353) q[1];
sx q[1];
rz(1.2870671) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6223479) q[0];
sx q[0];
rz(-2.17832) q[0];
sx q[0];
rz(2.1871479) q[0];
rz(-pi) q[1];
x q[1];
rz(1.14187) q[2];
sx q[2];
rz(-0.95488906) q[2];
sx q[2];
rz(-2.3162637) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.27582622) q[1];
sx q[1];
rz(-0.681502) q[1];
sx q[1];
rz(-1.1033415) q[1];
rz(0.98388328) q[3];
sx q[3];
rz(-2.623284) q[3];
sx q[3];
rz(-2.2494418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53671545) q[2];
sx q[2];
rz(-0.36467364) q[2];
sx q[2];
rz(-0.019901179) q[2];
rz(0.59602916) q[3];
sx q[3];
rz(-0.28087956) q[3];
sx q[3];
rz(1.178044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48536456) q[0];
sx q[0];
rz(-1.8940268) q[0];
sx q[0];
rz(-1.1744936) q[0];
rz(2.8354722) q[1];
sx q[1];
rz(-2.0968292) q[1];
sx q[1];
rz(-1.3694481) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5590483) q[0];
sx q[0];
rz(-1.567548) q[0];
sx q[0];
rz(1.5772468) q[0];
rz(-pi) q[1];
rz(1.3197863) q[2];
sx q[2];
rz(-0.6133056) q[2];
sx q[2];
rz(-1.6005188) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3432797) q[1];
sx q[1];
rz(-1.7213744) q[1];
sx q[1];
rz(1.9184789) q[1];
rz(-pi) q[2];
rz(-1.9846652) q[3];
sx q[3];
rz(-2.0479432) q[3];
sx q[3];
rz(-0.26354313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29641637) q[2];
sx q[2];
rz(-0.74411074) q[2];
sx q[2];
rz(0.80294341) q[2];
rz(2.0261649) q[3];
sx q[3];
rz(-0.69474703) q[3];
sx q[3];
rz(1.3692793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95336103) q[0];
sx q[0];
rz(-0.53934923) q[0];
sx q[0];
rz(0.11235919) q[0];
rz(0.27680963) q[1];
sx q[1];
rz(-1.9748297) q[1];
sx q[1];
rz(-0.64754957) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.707412) q[0];
sx q[0];
rz(-2.4828504) q[0];
sx q[0];
rz(-1.0007658) q[0];
rz(-pi) q[1];
rz(-2.4729257) q[2];
sx q[2];
rz(-1.200182) q[2];
sx q[2];
rz(-2.6636555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6426443) q[1];
sx q[1];
rz(-2.0181839) q[1];
sx q[1];
rz(-1.4397463) q[1];
rz(-1.3350491) q[3];
sx q[3];
rz(-1.4000579) q[3];
sx q[3];
rz(-0.52908191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.67927805) q[2];
sx q[2];
rz(-2.9913112) q[2];
sx q[2];
rz(1.3432937) q[2];
rz(-0.51370931) q[3];
sx q[3];
rz(-1.4880344) q[3];
sx q[3];
rz(0.59521365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3419471) q[0];
sx q[0];
rz(-1.3222313) q[0];
sx q[0];
rz(2.7213726) q[0];
rz(1.7814024) q[1];
sx q[1];
rz(-1.9748634) q[1];
sx q[1];
rz(2.3374048) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4046832) q[0];
sx q[0];
rz(-2.1054287) q[0];
sx q[0];
rz(2.695309) q[0];
x q[1];
rz(-1.8558414) q[2];
sx q[2];
rz(-0.80696304) q[2];
sx q[2];
rz(2.3318044) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0252991) q[1];
sx q[1];
rz(-2.5466047) q[1];
sx q[1];
rz(-1.3061252) q[1];
rz(-2.1523989) q[3];
sx q[3];
rz(-0.19671725) q[3];
sx q[3];
rz(-2.2328053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5899902) q[2];
sx q[2];
rz(-2.1253822) q[2];
sx q[2];
rz(1.972398) q[2];
rz(2.9679838) q[3];
sx q[3];
rz(-2.2326525) q[3];
sx q[3];
rz(1.199523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6731897) q[0];
sx q[0];
rz(-1.0262187) q[0];
sx q[0];
rz(1.9386559) q[0];
rz(0.24083336) q[1];
sx q[1];
rz(-2.2151561) q[1];
sx q[1];
rz(-0.9333207) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4881977) q[0];
sx q[0];
rz(-1.1943333) q[0];
sx q[0];
rz(0.36247902) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1432518) q[2];
sx q[2];
rz(-1.7939374) q[2];
sx q[2];
rz(0.23912341) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.49737469) q[1];
sx q[1];
rz(-1.3709732) q[1];
sx q[1];
rz(2.7544061) q[1];
x q[2];
rz(-2.9894514) q[3];
sx q[3];
rz(-2.1295432) q[3];
sx q[3];
rz(1.5882176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11747083) q[2];
sx q[2];
rz(-1.6880219) q[2];
sx q[2];
rz(3.0061099) q[2];
rz(2.1425715) q[3];
sx q[3];
rz(-2.8935367) q[3];
sx q[3];
rz(-2.5854056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4713772) q[0];
sx q[0];
rz(-1.1286292) q[0];
sx q[0];
rz(1.7673329) q[0];
rz(-1.0820092) q[1];
sx q[1];
rz(-2.3550985) q[1];
sx q[1];
rz(-1.5622004) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.608273) q[0];
sx q[0];
rz(-1.0315622) q[0];
sx q[0];
rz(3.1031552) q[0];
rz(1.9553929) q[2];
sx q[2];
rz(-2.0973258) q[2];
sx q[2];
rz(1.0004472) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1154707) q[1];
sx q[1];
rz(-1.9023832) q[1];
sx q[1];
rz(-2.0381171) q[1];
x q[2];
rz(2.0252941) q[3];
sx q[3];
rz(-0.70279944) q[3];
sx q[3];
rz(1.1490319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18206) q[2];
sx q[2];
rz(-1.5740732) q[2];
sx q[2];
rz(2.5173397) q[2];
rz(2.4652081) q[3];
sx q[3];
rz(-1.9703777) q[3];
sx q[3];
rz(-0.83522767) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4312209) q[0];
sx q[0];
rz(-1.2483163) q[0];
sx q[0];
rz(-1.7927908) q[0];
rz(-2.8377332) q[1];
sx q[1];
rz(-1.5307348) q[1];
sx q[1];
rz(-2.2644728) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0766834) q[0];
sx q[0];
rz(-0.8479894) q[0];
sx q[0];
rz(1.377489) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2459846) q[2];
sx q[2];
rz(-1.6742252) q[2];
sx q[2];
rz(1.6425486) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.408016) q[1];
sx q[1];
rz(-2.004679) q[1];
sx q[1];
rz(-2.7201319) q[1];
x q[2];
rz(-2.5780664) q[3];
sx q[3];
rz(-1.0358264) q[3];
sx q[3];
rz(0.68171147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0861686) q[2];
sx q[2];
rz(-0.40579) q[2];
sx q[2];
rz(-1.0322734) q[2];
rz(-2.0530733) q[3];
sx q[3];
rz(-1.4884596) q[3];
sx q[3];
rz(0.77435875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1879723) q[0];
sx q[0];
rz(-2.3379876) q[0];
sx q[0];
rz(0.048820989) q[0];
rz(1.8232518) q[1];
sx q[1];
rz(-1.7346458) q[1];
sx q[1];
rz(0.55191747) q[1];
rz(-1.3011408) q[2];
sx q[2];
rz(-1.9741304) q[2];
sx q[2];
rz(-1.3000549) q[2];
rz(-2.0069176) q[3];
sx q[3];
rz(-1.5831309) q[3];
sx q[3];
rz(2.3203608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
