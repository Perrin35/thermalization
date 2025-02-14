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
rz(-1.2920657) q[0];
sx q[0];
rz(-2.3176365) q[0];
sx q[0];
rz(0.92744654) q[0];
rz(-1.0499586) q[1];
sx q[1];
rz(-2.1702622) q[1];
sx q[1];
rz(-2.3271022) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2353219) q[0];
sx q[0];
rz(-2.5320279) q[0];
sx q[0];
rz(0.8602575) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79618228) q[2];
sx q[2];
rz(-2.0416235) q[2];
sx q[2];
rz(1.7288961) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51920096) q[1];
sx q[1];
rz(-1.0990874) q[1];
sx q[1];
rz(0.5640819) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3130034) q[3];
sx q[3];
rz(-1.0416721) q[3];
sx q[3];
rz(0.55381004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77259511) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(2.7102846) q[2];
rz(1.6469693) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(-2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5539219) q[0];
sx q[0];
rz(-2.9449154) q[0];
sx q[0];
rz(1.824463) q[0];
rz(1.0493086) q[1];
sx q[1];
rz(-0.53593719) q[1];
sx q[1];
rz(2.5228693) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65370377) q[0];
sx q[0];
rz(-1.7152707) q[0];
sx q[0];
rz(-2.5313951) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0929865) q[2];
sx q[2];
rz(-0.75121643) q[2];
sx q[2];
rz(-2.4260245) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6438457) q[1];
sx q[1];
rz(-2.7676555) q[1];
sx q[1];
rz(-1.9356329) q[1];
x q[2];
rz(1.9369164) q[3];
sx q[3];
rz(-1.06166) q[3];
sx q[3];
rz(-2.3222271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.416136) q[2];
sx q[2];
rz(-1.340103) q[2];
sx q[2];
rz(0.49187342) q[2];
rz(-1.5791996) q[3];
sx q[3];
rz(-1.4926566) q[3];
sx q[3];
rz(2.5621342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270776) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(0.28011093) q[0];
rz(-0.07490553) q[1];
sx q[1];
rz(-0.43498755) q[1];
sx q[1];
rz(-1.3704971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460187) q[0];
sx q[0];
rz(-1.7798806) q[0];
sx q[0];
rz(-0.81935291) q[0];
rz(0.89981516) q[2];
sx q[2];
rz(-2.0986989) q[2];
sx q[2];
rz(-1.4769276) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0947511) q[1];
sx q[1];
rz(-1.0221204) q[1];
sx q[1];
rz(-1.8803174) q[1];
rz(-pi) q[2];
rz(1.5896443) q[3];
sx q[3];
rz(-1.4385838) q[3];
sx q[3];
rz(2.978505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8568153) q[2];
sx q[2];
rz(-1.4692401) q[2];
sx q[2];
rz(2.3397297) q[2];
rz(-2.2297468) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(-0.97889939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3367679) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(-0.62070745) q[0];
rz(0.2785109) q[1];
sx q[1];
rz(-1.4418437) q[1];
sx q[1];
rz(2.1381569) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0898363) q[0];
sx q[0];
rz(-0.47430719) q[0];
sx q[0];
rz(-0.82292212) q[0];
x q[1];
rz(3.0714967) q[2];
sx q[2];
rz(-1.9906487) q[2];
sx q[2];
rz(1.1481089) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.675337) q[1];
sx q[1];
rz(-1.3020483) q[1];
sx q[1];
rz(-1.6273725) q[1];
x q[2];
rz(-0.33903709) q[3];
sx q[3];
rz(-2.5154696) q[3];
sx q[3];
rz(-0.98527415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3044546) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(-0.066085286) q[2];
rz(1.1262013) q[3];
sx q[3];
rz(-0.76032138) q[3];
sx q[3];
rz(-2.5662305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5758301) q[0];
sx q[0];
rz(-0.081871651) q[0];
sx q[0];
rz(0.78731147) q[0];
rz(2.2832504) q[1];
sx q[1];
rz(-2.1635677) q[1];
sx q[1];
rz(-1.0816921) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3560548) q[0];
sx q[0];
rz(-2.303146) q[0];
sx q[0];
rz(0.69489702) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41795798) q[2];
sx q[2];
rz(-2.0047028) q[2];
sx q[2];
rz(1.5669296) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3216721) q[1];
sx q[1];
rz(-1.5645488) q[1];
sx q[1];
rz(2.856918) q[1];
x q[2];
rz(-0.95329625) q[3];
sx q[3];
rz(-0.82103182) q[3];
sx q[3];
rz(-1.2987069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16296884) q[2];
sx q[2];
rz(-1.9917515) q[2];
sx q[2];
rz(0.84257379) q[2];
rz(0.27635559) q[3];
sx q[3];
rz(-2.1601951) q[3];
sx q[3];
rz(1.792255) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2862947) q[0];
sx q[0];
rz(-1.8305625) q[0];
sx q[0];
rz(1.1899765) q[0];
rz(1.2225993) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(0.58293265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.615743) q[0];
sx q[0];
rz(-1.9186363) q[0];
sx q[0];
rz(-2.6418532) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1276064) q[2];
sx q[2];
rz(-0.99726935) q[2];
sx q[2];
rz(-1.397246) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8376298) q[1];
sx q[1];
rz(-0.25165855) q[1];
sx q[1];
rz(0.32984094) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0058026) q[3];
sx q[3];
rz(-0.9851176) q[3];
sx q[3];
rz(-3.1068561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0827834) q[2];
sx q[2];
rz(-1.8913816) q[2];
sx q[2];
rz(1.0397376) q[2];
rz(2.5522363) q[3];
sx q[3];
rz(-1.6683234) q[3];
sx q[3];
rz(-2.1350433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650918) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(-1.3366706) q[0];
rz(-0.22557766) q[1];
sx q[1];
rz(-2.070319) q[1];
sx q[1];
rz(-0.76990661) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5107489) q[0];
sx q[0];
rz(-2.0389557) q[0];
sx q[0];
rz(1.7450784) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82775292) q[2];
sx q[2];
rz(-2.0619287) q[2];
sx q[2];
rz(2.998005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47178956) q[1];
sx q[1];
rz(-0.21179971) q[1];
sx q[1];
rz(-2.7752911) q[1];
rz(-pi) q[2];
rz(0.1128686) q[3];
sx q[3];
rz(-0.60823554) q[3];
sx q[3];
rz(3.1361767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.3930438) q[2];
sx q[2];
rz(1.0153655) q[2];
rz(-0.20396248) q[3];
sx q[3];
rz(-1.5927619) q[3];
sx q[3];
rz(1.7502194) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2400951) q[0];
sx q[0];
rz(-2.5101341) q[0];
sx q[0];
rz(-0.39328662) q[0];
rz(-2.6914864) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(0.030287655) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5022885) q[0];
sx q[0];
rz(-0.24950108) q[0];
sx q[0];
rz(-0.12608237) q[0];
x q[1];
rz(-0.25595902) q[2];
sx q[2];
rz(-0.82590196) q[2];
sx q[2];
rz(-1.9108552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5133678) q[1];
sx q[1];
rz(-1.0252044) q[1];
sx q[1];
rz(-1.5978769) q[1];
x q[2];
rz(2.8775086) q[3];
sx q[3];
rz(-1.698016) q[3];
sx q[3];
rz(0.41632392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0344737) q[2];
sx q[2];
rz(-1.5566885) q[2];
sx q[2];
rz(-1.1220453) q[2];
rz(1.0551039) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(1.8900227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17290641) q[0];
sx q[0];
rz(-0.98772573) q[0];
sx q[0];
rz(-0.77447844) q[0];
rz(-0.6340181) q[1];
sx q[1];
rz(-2.1114383) q[1];
sx q[1];
rz(2.074362) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7726583) q[0];
sx q[0];
rz(-0.83266363) q[0];
sx q[0];
rz(-1.5475818) q[0];
rz(1.2634694) q[2];
sx q[2];
rz(-1.5241703) q[2];
sx q[2];
rz(-1.7152552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6744819) q[1];
sx q[1];
rz(-1.8201105) q[1];
sx q[1];
rz(1.8915908) q[1];
rz(1.2379856) q[3];
sx q[3];
rz(-0.5026256) q[3];
sx q[3];
rz(-2.7183333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5344703) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(1.2657451) q[2];
rz(-3.1014465) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(-3.1025187) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48113111) q[0];
sx q[0];
rz(-2.2534695) q[0];
sx q[0];
rz(1.2598502) q[0];
rz(-1.6549567) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(3.1390417) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1387864) q[0];
sx q[0];
rz(-1.696179) q[0];
sx q[0];
rz(3.0210397) q[0];
rz(1.4223406) q[2];
sx q[2];
rz(-3.0768148) q[2];
sx q[2];
rz(-1.5671519) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.59460179) q[1];
sx q[1];
rz(-2.3560215) q[1];
sx q[1];
rz(-2.5681096) q[1];
x q[2];
rz(1.5137767) q[3];
sx q[3];
rz(-2.7064763) q[3];
sx q[3];
rz(-0.42729898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5551787) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(-0.67443887) q[2];
rz(1.1174348) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(-2.8286772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085070327) q[0];
sx q[0];
rz(-2.0105579) q[0];
sx q[0];
rz(2.6149268) q[0];
rz(2.7936735) q[1];
sx q[1];
rz(-1.2384474) q[1];
sx q[1];
rz(1.7927982) q[1];
rz(1.3580218) q[2];
sx q[2];
rz(-1.2322896) q[2];
sx q[2];
rz(-0.88862669) q[2];
rz(-1.450235) q[3];
sx q[3];
rz(-0.14738722) q[3];
sx q[3];
rz(1.9167388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
