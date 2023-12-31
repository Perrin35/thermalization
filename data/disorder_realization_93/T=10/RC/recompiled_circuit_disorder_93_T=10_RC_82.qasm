OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(-1.260489) q[0];
rz(2.1029544) q[1];
sx q[1];
rz(-1.3488052) q[1];
sx q[1];
rz(0.92372149) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0215065) q[0];
sx q[0];
rz(-1.9061631) q[0];
sx q[0];
rz(-1.0943227) q[0];
x q[1];
rz(1.107723) q[2];
sx q[2];
rz(-2.827364) q[2];
sx q[2];
rz(-2.3067834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94581374) q[1];
sx q[1];
rz(-2.0683156) q[1];
sx q[1];
rz(-1.0832018) q[1];
rz(-pi) q[2];
x q[2];
rz(1.190891) q[3];
sx q[3];
rz(-2.0348747) q[3];
sx q[3];
rz(-0.78117785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91360056) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(0.16201924) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(2.4285765) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0682003) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(-1.1967999) q[0];
rz(0.67990047) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.686036) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3128132) q[0];
sx q[0];
rz(-1.4038741) q[0];
sx q[0];
rz(-2.1524327) q[0];
rz(-pi) q[1];
rz(1.4886841) q[2];
sx q[2];
rz(-1.1973235) q[2];
sx q[2];
rz(0.79370802) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62555186) q[1];
sx q[1];
rz(-1.5572773) q[1];
sx q[1];
rz(1.0479755) q[1];
rz(-1.5973813) q[3];
sx q[3];
rz(-1.1592602) q[3];
sx q[3];
rz(0.4707903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(1.7896174) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(2.341111) q[0];
rz(0.02877409) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(1.9690537) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81329936) q[0];
sx q[0];
rz(-1.0070224) q[0];
sx q[0];
rz(-1.9378807) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8981947) q[2];
sx q[2];
rz(-2.2990169) q[2];
sx q[2];
rz(0.74795216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2492003) q[1];
sx q[1];
rz(-1.8141659) q[1];
sx q[1];
rz(1.0358441) q[1];
rz(-pi) q[2];
rz(1.8215239) q[3];
sx q[3];
rz(-1.6276976) q[3];
sx q[3];
rz(-1.0292605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0744434) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(-0.91119901) q[2];
rz(0.95101142) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(-2.2495911) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38055414) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(-3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(2.6180843) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2142221) q[0];
sx q[0];
rz(-0.71338755) q[0];
sx q[0];
rz(0.58332304) q[0];
rz(-pi) q[1];
rz(1.6558311) q[2];
sx q[2];
rz(-1.8997314) q[2];
sx q[2];
rz(-2.0563682) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47141155) q[1];
sx q[1];
rz(-1.8243316) q[1];
sx q[1];
rz(2.2815435) q[1];
x q[2];
rz(3.0047699) q[3];
sx q[3];
rz(-0.98383437) q[3];
sx q[3];
rz(0.50859355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(-2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(-2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(-0.87961698) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-2.1496444) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6742453) q[0];
sx q[0];
rz(-1.464198) q[0];
sx q[0];
rz(-0.17515134) q[0];
x q[1];
rz(-1.5422103) q[2];
sx q[2];
rz(-2.8386142) q[2];
sx q[2];
rz(-2.6945393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76584133) q[1];
sx q[1];
rz(-1.5305133) q[1];
sx q[1];
rz(0.30893107) q[1];
x q[2];
rz(2.3283142) q[3];
sx q[3];
rz(-1.4762029) q[3];
sx q[3];
rz(2.7300342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-0.27080718) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72702423) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(-1.3487934) q[0];
rz(-0.38189608) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(1.4250925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2857367) q[0];
sx q[0];
rz(-2.5897313) q[0];
sx q[0];
rz(0.43666552) q[0];
rz(-3.0503057) q[2];
sx q[2];
rz(-1.316615) q[2];
sx q[2];
rz(2.9448178) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6188366) q[1];
sx q[1];
rz(-0.4546051) q[1];
sx q[1];
rz(-1.418581) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85961996) q[3];
sx q[3];
rz(-2.4058127) q[3];
sx q[3];
rz(-1.3501292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1340593) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(-2.9610736) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-0.40294161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5086223) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-0.41931835) q[0];
rz(-1.5527027) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(-2.3197876) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1497027) q[0];
sx q[0];
rz(-0.8677965) q[0];
sx q[0];
rz(2.4292612) q[0];
rz(0.84341151) q[2];
sx q[2];
rz(-1.1864098) q[2];
sx q[2];
rz(-0.20993983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8824749) q[1];
sx q[1];
rz(-1.3438517) q[1];
sx q[1];
rz(-2.4005753) q[1];
rz(0.58737289) q[3];
sx q[3];
rz(-1.2578739) q[3];
sx q[3];
rz(-1.8250993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3372779) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(2.896893) q[2];
rz(-3.0120567) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(-0.82292557) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(1.3051422) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6179498) q[0];
sx q[0];
rz(-2.4542913) q[0];
sx q[0];
rz(0.53074093) q[0];
rz(-pi) q[1];
rz(-0.42758503) q[2];
sx q[2];
rz(-0.75140778) q[2];
sx q[2];
rz(2.4497355) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.23755632) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(-1.1035641) q[1];
x q[2];
rz(-0.023530258) q[3];
sx q[3];
rz(-1.2315893) q[3];
sx q[3];
rz(-0.48285218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.4902327) q[2];
rz(-2.0643318) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867144) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(-0.3219147) q[0];
rz(-1.6053258) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(0.70294356) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212595) q[0];
sx q[0];
rz(-2.5998305) q[0];
sx q[0];
rz(0.74777491) q[0];
rz(-1.4998352) q[2];
sx q[2];
rz(-0.5103726) q[2];
sx q[2];
rz(-0.77529782) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8150836) q[1];
sx q[1];
rz(-2.106296) q[1];
sx q[1];
rz(-2.9498847) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1861107) q[3];
sx q[3];
rz(-1.8468879) q[3];
sx q[3];
rz(2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90074173) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(1.3396324) q[2];
rz(0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16383485) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(-1.9158069) q[0];
rz(-2.2380791) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43384957) q[0];
sx q[0];
rz(-1.187547) q[0];
sx q[0];
rz(-1.3462523) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2655728) q[2];
sx q[2];
rz(-2.6555736) q[2];
sx q[2];
rz(-1.3181869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74433078) q[1];
sx q[1];
rz(-1.8363734) q[1];
sx q[1];
rz(0.070752146) q[1];
rz(0.81810276) q[3];
sx q[3];
rz(-0.59554282) q[3];
sx q[3];
rz(-2.9593352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6580711) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(-1.998385) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951915) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(-2.519683) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(2.4516104) q[2];
sx q[2];
rz(-2.1524515) q[2];
sx q[2];
rz(-0.099302789) q[2];
rz(2.2223496) q[3];
sx q[3];
rz(-1.5083434) q[3];
sx q[3];
rz(-1.2283243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
