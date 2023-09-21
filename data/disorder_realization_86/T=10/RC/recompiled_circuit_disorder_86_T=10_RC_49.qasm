OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(1.4120742) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(-0.59564367) q[1];
sx q[1];
rz(-1.6593978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94443653) q[0];
sx q[0];
rz(-1.7261788) q[0];
sx q[0];
rz(0.42874254) q[0];
rz(-pi) q[1];
rz(2.8843845) q[2];
sx q[2];
rz(-1.6991985) q[2];
sx q[2];
rz(-0.53127015) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.23956242) q[1];
sx q[1];
rz(-0.67443496) q[1];
sx q[1];
rz(0.85427888) q[1];
rz(-1.4115303) q[3];
sx q[3];
rz(-2.5730238) q[3];
sx q[3];
rz(-2.0025314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98510629) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(2.1872897) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(-1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(3.1153733) q[0];
rz(-1.5401309) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-0.96347934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.11461) q[0];
sx q[0];
rz(-0.61404213) q[0];
sx q[0];
rz(-1.5668037) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33032592) q[2];
sx q[2];
rz(-1.0268372) q[2];
sx q[2];
rz(-2.3842173) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0986833) q[1];
sx q[1];
rz(-2.5334362) q[1];
sx q[1];
rz(0.98867464) q[1];
rz(1.9217334) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(1.2082781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5144689) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(3.0070686) q[2];
rz(2.3965805) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-0.79743687) q[0];
rz(2.0939317) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(-0.55999666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3968351) q[0];
sx q[0];
rz(-2.7042537) q[0];
sx q[0];
rz(1.0768946) q[0];
rz(-1.5494924) q[2];
sx q[2];
rz(-1.2676123) q[2];
sx q[2];
rz(1.4959178) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3126038) q[1];
sx q[1];
rz(-1.7692411) q[1];
sx q[1];
rz(0.36638422) q[1];
rz(2.3476944) q[3];
sx q[3];
rz(-0.63459914) q[3];
sx q[3];
rz(1.8397699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3893163) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(-0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(2.0836232) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383706) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(2.8248887) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(-1.2987312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7553058) q[0];
sx q[0];
rz(-0.55314976) q[0];
sx q[0];
rz(1.132071) q[0];
rz(-2.3372075) q[2];
sx q[2];
rz(-1.569869) q[2];
sx q[2];
rz(-1.5915807) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1620996) q[1];
sx q[1];
rz(-2.5433308) q[1];
sx q[1];
rz(1.23566) q[1];
x q[2];
rz(1.8989765) q[3];
sx q[3];
rz(-1.4644074) q[3];
sx q[3];
rz(-2.1387517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5359042) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(-2.7569125) q[2];
rz(-0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(1.4543021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.6032747) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(-1.7549365) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(0.2968266) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41243991) q[0];
sx q[0];
rz(-1.6099596) q[0];
sx q[0];
rz(-0.57106437) q[0];
rz(-pi) q[1];
rz(1.7237687) q[2];
sx q[2];
rz(-0.83380552) q[2];
sx q[2];
rz(2.134915) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8991124) q[1];
sx q[1];
rz(-0.63226262) q[1];
sx q[1];
rz(-2.2805023) q[1];
rz(-pi) q[2];
rz(1.4818707) q[3];
sx q[3];
rz(-2.2825135) q[3];
sx q[3];
rz(-2.5951648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.37830535) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(-2.7491167) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.0157938) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(-1.8136576) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(-0.60633916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7327001) q[0];
sx q[0];
rz(-0.048763976) q[0];
sx q[0];
rz(-2.7932037) q[0];
rz(2.7251284) q[2];
sx q[2];
rz(-0.93566862) q[2];
sx q[2];
rz(3.0481899) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.11583466) q[1];
sx q[1];
rz(-2.5587213) q[1];
sx q[1];
rz(-1.9028266) q[1];
rz(-pi) q[2];
x q[2];
rz(1.243152) q[3];
sx q[3];
rz(-0.22938211) q[3];
sx q[3];
rz(-2.4250507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5027344) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(-2.0022557) q[2];
rz(1.6566488) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6095603) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(-2.6334921) q[0];
rz(1.5628901) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(-2.3513444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9433141) q[0];
sx q[0];
rz(-0.6753079) q[0];
sx q[0];
rz(2.6576463) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21364084) q[2];
sx q[2];
rz(-1.5593312) q[2];
sx q[2];
rz(-1.2917047) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8040647) q[1];
sx q[1];
rz(-1.84066) q[1];
sx q[1];
rz(-2.7359664) q[1];
x q[2];
rz(1.1931476) q[3];
sx q[3];
rz(-0.95818633) q[3];
sx q[3];
rz(2.2539504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.053085176) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(-1.6648071) q[3];
sx q[3];
rz(-2.1063185) q[3];
sx q[3];
rz(-0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4247894) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(-2.5324902) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(1.3887127) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4318651) q[0];
sx q[0];
rz(-1.0114397) q[0];
sx q[0];
rz(3.1064242) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1461166) q[2];
sx q[2];
rz(-1.0815902) q[2];
sx q[2];
rz(-0.81791544) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70506239) q[1];
sx q[1];
rz(-0.092518004) q[1];
sx q[1];
rz(-3.0180879) q[1];
rz(-pi) q[2];
rz(-1.1561398) q[3];
sx q[3];
rz(-2.4932043) q[3];
sx q[3];
rz(0.63601953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1887112) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(2.2593373) q[2];
rz(-1.7404209) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8700478) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(1.4260938) q[0];
rz(3.0601314) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(2.5833599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8919864) q[0];
sx q[0];
rz(-2.5476646) q[0];
sx q[0];
rz(2.3114165) q[0];
x q[1];
rz(1.6782645) q[2];
sx q[2];
rz(-1.5170013) q[2];
sx q[2];
rz(1.6966284) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0135632) q[1];
sx q[1];
rz(-2.3131144) q[1];
sx q[1];
rz(-3.0548884) q[1];
rz(-2.6551412) q[3];
sx q[3];
rz(-1.7351741) q[3];
sx q[3];
rz(-2.1742976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(-0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50487173) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(1.6037534) q[0];
rz(2.3161855) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-2.6182981) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6494203) q[0];
sx q[0];
rz(-1.5086552) q[0];
sx q[0];
rz(-0.6092351) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.043838219) q[2];
sx q[2];
rz(-0.99928108) q[2];
sx q[2];
rz(2.8457355) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1112422) q[1];
sx q[1];
rz(-2.6551464) q[1];
sx q[1];
rz(-2.4187947) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8966122) q[3];
sx q[3];
rz(-1.3462726) q[3];
sx q[3];
rz(0.45104879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.837073) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(-2.8722897) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8158648) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(0.011209839) q[2];
sx q[2];
rz(-1.8335473) q[2];
sx q[2];
rz(2.0469472) q[2];
rz(-0.79694637) q[3];
sx q[3];
rz(-1.7399825) q[3];
sx q[3];
rz(0.83132838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];