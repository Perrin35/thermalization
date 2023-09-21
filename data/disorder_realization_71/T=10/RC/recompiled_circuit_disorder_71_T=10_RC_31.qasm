OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(0.82984501) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(-0.87632626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012955879) q[0];
sx q[0];
rz(-0.46015938) q[0];
sx q[0];
rz(-0.53928661) q[0];
rz(-2.949602) q[2];
sx q[2];
rz(-0.99388323) q[2];
sx q[2];
rz(2.752395) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10325601) q[1];
sx q[1];
rz(-1.5557319) q[1];
sx q[1];
rz(1.9008093) q[1];
x q[2];
rz(-1.417744) q[3];
sx q[3];
rz(-1.9733841) q[3];
sx q[3];
rz(-3.0229085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(0.7286287) q[2];
rz(0.5209926) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8347297) q[0];
sx q[0];
rz(-1.9711718) q[0];
sx q[0];
rz(2.0200502) q[0];
rz(-2.8858378) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(2.2671525) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53423184) q[0];
sx q[0];
rz(-0.53593862) q[0];
sx q[0];
rz(0.94632728) q[0];
rz(-pi) q[1];
rz(2.1463257) q[2];
sx q[2];
rz(-1.4412291) q[2];
sx q[2];
rz(1.0075943) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.16972152) q[1];
sx q[1];
rz(-2.019878) q[1];
sx q[1];
rz(0.68193087) q[1];
rz(-3.087895) q[3];
sx q[3];
rz(-1.3018381) q[3];
sx q[3];
rz(0.69124903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4008537) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(-2.7056616) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(0.40288231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9044559) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(-1.0748192) q[0];
rz(0.83956051) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(-2.7456465) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0229189) q[0];
sx q[0];
rz(-0.6624822) q[0];
sx q[0];
rz(-2.2454717) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6171574) q[2];
sx q[2];
rz(-2.8472387) q[2];
sx q[2];
rz(-2.4183395) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0242651) q[1];
sx q[1];
rz(-1.5803442) q[1];
sx q[1];
rz(-0.68318232) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.08926908) q[3];
sx q[3];
rz(-2.8850728) q[3];
sx q[3];
rz(1.4005816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5376771) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(2.9023857) q[2];
rz(3.0662597) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72162119) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(-2.3216632) q[0];
rz(-0.48768249) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(-2.908169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4963835) q[0];
sx q[0];
rz(-1.5102981) q[0];
sx q[0];
rz(1.4738826) q[0];
rz(-1.3350305) q[2];
sx q[2];
rz(-2.1996017) q[2];
sx q[2];
rz(2.6218888) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.0049131752) q[1];
sx q[1];
rz(-0.79489743) q[1];
sx q[1];
rz(2.6585048) q[1];
rz(-2.7531284) q[3];
sx q[3];
rz(-0.82287517) q[3];
sx q[3];
rz(-0.32271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0507811) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(-0.74742571) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915879) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(3.0773556) q[0];
rz(-2.1977987) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(-2.3805526) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31537406) q[0];
sx q[0];
rz(-1.5407019) q[0];
sx q[0];
rz(-1.8316359) q[0];
rz(-pi) q[1];
rz(0.20247395) q[2];
sx q[2];
rz(-1.379181) q[2];
sx q[2];
rz(-1.3445878) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0423454) q[1];
sx q[1];
rz(-1.5464916) q[1];
sx q[1];
rz(0.32056067) q[1];
rz(-pi) q[2];
rz(2.8588572) q[3];
sx q[3];
rz(-1.8394107) q[3];
sx q[3];
rz(0.80073592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3349907) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(-3.0495194) q[2];
rz(0.66172415) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(-1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-2.6732366) q[0];
sx q[0];
rz(-1.8259003) q[0];
sx q[0];
rz(0.026542149) q[0];
rz(-2.2684855) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(0.32593265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1074166) q[0];
sx q[0];
rz(-1.0506442) q[0];
sx q[0];
rz(1.1362032) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58165254) q[2];
sx q[2];
rz(-2.622421) q[2];
sx q[2];
rz(0.66228629) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0210452) q[1];
sx q[1];
rz(-1.435934) q[1];
sx q[1];
rz(-0.30121505) q[1];
rz(-pi) q[2];
rz(0.10109191) q[3];
sx q[3];
rz(-1.5875845) q[3];
sx q[3];
rz(-1.8198131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5439593) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(-1.7116961) q[3];
sx q[3];
rz(-1.1681898) q[3];
sx q[3];
rz(2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8687246) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(2.4196999) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(3.022335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8598547) q[0];
sx q[0];
rz(-0.87346948) q[0];
sx q[0];
rz(1.9182693) q[0];
rz(-pi) q[1];
rz(2.1642045) q[2];
sx q[2];
rz(-0.19492976) q[2];
sx q[2];
rz(-0.56602636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6357928) q[1];
sx q[1];
rz(-2.0675106) q[1];
sx q[1];
rz(-2.990681) q[1];
rz(-2.7151832) q[3];
sx q[3];
rz(-0.75424131) q[3];
sx q[3];
rz(-2.756556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(-1.7822441) q[2];
rz(2.3826777) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(-2.6045077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(2.0781562) q[0];
rz(-2.8670782) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(2.2559821) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0622334) q[0];
sx q[0];
rz(-1.0842807) q[0];
sx q[0];
rz(-1.833672) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9912668) q[2];
sx q[2];
rz(-1.5577003) q[2];
sx q[2];
rz(3.1377369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1427666) q[1];
sx q[1];
rz(-0.66122675) q[1];
sx q[1];
rz(2.2426474) q[1];
rz(-pi) q[2];
rz(2.2742911) q[3];
sx q[3];
rz(-2.2059545) q[3];
sx q[3];
rz(-2.9150073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4132335) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(-2.0098861) q[2];
rz(1.0845832) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(-1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(-2.0595179) q[0];
rz(-1.2754296) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(2.0057604) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95490676) q[0];
sx q[0];
rz(-1.0404772) q[0];
sx q[0];
rz(-2.9493939) q[0];
rz(-pi) q[1];
x q[1];
rz(0.047770569) q[2];
sx q[2];
rz(-2.7063745) q[2];
sx q[2];
rz(-1.4620632) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3867934) q[1];
sx q[1];
rz(-0.42897412) q[1];
sx q[1];
rz(3.0241443) q[1];
rz(-1.0335835) q[3];
sx q[3];
rz(-0.72519231) q[3];
sx q[3];
rz(-2.4560526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.11848005) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(2.2980799) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(-0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89235598) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(1.2783485) q[0];
rz(1.0247914) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(1.1970253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14180627) q[0];
sx q[0];
rz(-1.0021266) q[0];
sx q[0];
rz(-0.52642157) q[0];
rz(-1.808666) q[2];
sx q[2];
rz(-1.6878205) q[2];
sx q[2];
rz(-2.7931917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5843643) q[1];
sx q[1];
rz(-0.96712501) q[1];
sx q[1];
rz(2.5283458) q[1];
x q[2];
rz(-0.50935575) q[3];
sx q[3];
rz(-1.9037814) q[3];
sx q[3];
rz(1.8978564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3060351) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-2.3416134) q[2];
rz(1.9647313) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(-2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4326614) q[0];
sx q[0];
rz(-2.9928757) q[0];
sx q[0];
rz(-2.3401674) q[0];
rz(-2.6196383) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(0.37408806) q[2];
sx q[2];
rz(-1.903152) q[2];
sx q[2];
rz(-2.8852035) q[2];
rz(-0.63511499) q[3];
sx q[3];
rz(-1.0233581) q[3];
sx q[3];
rz(-0.81627853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
