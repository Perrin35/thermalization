OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6381792) q[0];
sx q[0];
rz(4.4210202) q[0];
sx q[0];
rz(11.790334) q[0];
rz(0.70932055) q[1];
sx q[1];
rz(-1.6242937) q[1];
sx q[1];
rz(-0.59901839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4871019) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(-1.0949277) q[0];
rz(-pi) q[1];
rz(2.8593596) q[2];
sx q[2];
rz(-2.2405409) q[2];
sx q[2];
rz(-1.9053659) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4091275) q[1];
sx q[1];
rz(-2.5149087) q[1];
sx q[1];
rz(-2.5123358) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3612411) q[3];
sx q[3];
rz(-1.2557185) q[3];
sx q[3];
rz(1.9596069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1385931) q[2];
sx q[2];
rz(-2.004576) q[2];
sx q[2];
rz(-2.9453759) q[2];
rz(2.0972706) q[3];
sx q[3];
rz(-2.315867) q[3];
sx q[3];
rz(2.830982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1021378) q[0];
sx q[0];
rz(-0.65334833) q[0];
sx q[0];
rz(2.2838604) q[0];
rz(-2.6013382) q[1];
sx q[1];
rz(-1.0539571) q[1];
sx q[1];
rz(-0.78261715) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40561134) q[0];
sx q[0];
rz(-0.95658871) q[0];
sx q[0];
rz(0.10615291) q[0];
rz(-pi) q[1];
rz(2.906053) q[2];
sx q[2];
rz(-0.96220926) q[2];
sx q[2];
rz(-2.4289301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2007717) q[1];
sx q[1];
rz(-1.7634974) q[1];
sx q[1];
rz(0.60639834) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88528648) q[3];
sx q[3];
rz(-1.1210223) q[3];
sx q[3];
rz(-1.7974082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9767849) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(-2.5505193) q[2];
rz(2.1137386) q[3];
sx q[3];
rz(-2.5058392) q[3];
sx q[3];
rz(-0.13636057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9117821) q[0];
sx q[0];
rz(-2.6214143) q[0];
sx q[0];
rz(2.0853364) q[0];
rz(-0.99110574) q[1];
sx q[1];
rz(-1.3737563) q[1];
sx q[1];
rz(0.80449218) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0206385) q[0];
sx q[0];
rz(-0.94392969) q[0];
sx q[0];
rz(1.9009717) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.095011906) q[2];
sx q[2];
rz(-2.1974051) q[2];
sx q[2];
rz(2.1206405) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.04136297) q[1];
sx q[1];
rz(-1.8112665) q[1];
sx q[1];
rz(2.9036456) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7630173) q[3];
sx q[3];
rz(-2.2092735) q[3];
sx q[3];
rz(0.22787991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6911917) q[2];
sx q[2];
rz(-0.85019008) q[2];
sx q[2];
rz(-0.56785339) q[2];
rz(-1.9495226) q[3];
sx q[3];
rz(-0.53693938) q[3];
sx q[3];
rz(-0.15538628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464226) q[0];
sx q[0];
rz(-1.3986873) q[0];
sx q[0];
rz(1.9544741) q[0];
rz(-0.25539708) q[1];
sx q[1];
rz(-1.7385769) q[1];
sx q[1];
rz(-0.38937169) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02006836) q[0];
sx q[0];
rz(-0.25857718) q[0];
sx q[0];
rz(0.54988396) q[0];
x q[1];
rz(-3.1115628) q[2];
sx q[2];
rz(-2.2976934) q[2];
sx q[2];
rz(2.6384357) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8592718) q[1];
sx q[1];
rz(-0.92073694) q[1];
sx q[1];
rz(2.9031624) q[1];
x q[2];
rz(-3.0087972) q[3];
sx q[3];
rz(-1.3758052) q[3];
sx q[3];
rz(2.8094069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4046459) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(-1.3673937) q[2];
rz(2.6020452) q[3];
sx q[3];
rz(-1.5893785) q[3];
sx q[3];
rz(-1.3246983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0668199) q[0];
sx q[0];
rz(-0.19817752) q[0];
sx q[0];
rz(2.5812126) q[0];
rz(-2.9226774) q[1];
sx q[1];
rz(-2.0603265) q[1];
sx q[1];
rz(-0.17094368) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4881011) q[0];
sx q[0];
rz(-2.2224036) q[0];
sx q[0];
rz(2.1733858) q[0];
x q[1];
rz(0.60336171) q[2];
sx q[2];
rz(-2.8953279) q[2];
sx q[2];
rz(-2.6880996) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5569607) q[1];
sx q[1];
rz(-2.0804394) q[1];
sx q[1];
rz(-0.72344785) q[1];
rz(-pi) q[2];
rz(0.070329026) q[3];
sx q[3];
rz(-2.189895) q[3];
sx q[3];
rz(-0.17507565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7818452) q[2];
sx q[2];
rz(-1.1283987) q[2];
sx q[2];
rz(-2.183059) q[2];
rz(-0.16252276) q[3];
sx q[3];
rz(-0.84238094) q[3];
sx q[3];
rz(1.5490279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.573134) q[0];
sx q[0];
rz(-0.064366654) q[0];
sx q[0];
rz(2.3934613) q[0];
rz(2.023078) q[1];
sx q[1];
rz(-2.7225814) q[1];
sx q[1];
rz(-2.8245139) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1823163) q[0];
sx q[0];
rz(-1.8166564) q[0];
sx q[0];
rz(1.6410752) q[0];
x q[1];
rz(1.6197312) q[2];
sx q[2];
rz(-1.3746007) q[2];
sx q[2];
rz(2.115415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4273964) q[1];
sx q[1];
rz(-2.3347989) q[1];
sx q[1];
rz(0.098735672) q[1];
rz(-pi) q[2];
rz(-1.5311144) q[3];
sx q[3];
rz(-2.3473266) q[3];
sx q[3];
rz(-3.119885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76724425) q[2];
sx q[2];
rz(-0.92864645) q[2];
sx q[2];
rz(1.7283776) q[2];
rz(0.080538571) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(2.9197689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(0.11248511) q[0];
sx q[0];
rz(-1.2487829) q[0];
sx q[0];
rz(-1.8674194) q[0];
rz(1.3451276) q[1];
sx q[1];
rz(-0.74256623) q[1];
sx q[1];
rz(-1.3412195) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8250834) q[0];
sx q[0];
rz(-2.1232455) q[0];
sx q[0];
rz(2.1567206) q[0];
x q[1];
rz(2.6599081) q[2];
sx q[2];
rz(-1.8284214) q[2];
sx q[2];
rz(3.1236908) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0643031) q[1];
sx q[1];
rz(-1.7619362) q[1];
sx q[1];
rz(-2.9851476) q[1];
rz(1.7246095) q[3];
sx q[3];
rz(-0.65585583) q[3];
sx q[3];
rz(0.7747618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93366569) q[2];
sx q[2];
rz(-1.9542481) q[2];
sx q[2];
rz(-2.6742317) q[2];
rz(0.76006877) q[3];
sx q[3];
rz(-1.6642539) q[3];
sx q[3];
rz(1.8925586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6811328) q[0];
sx q[0];
rz(-2.12119) q[0];
sx q[0];
rz(-1.6946174) q[0];
rz(-2.815822) q[1];
sx q[1];
rz(-0.21369801) q[1];
sx q[1];
rz(-1.6206585) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2163154) q[0];
sx q[0];
rz(-1.0641353) q[0];
sx q[0];
rz(1.8761047) q[0];
rz(-pi) q[1];
x q[1];
rz(0.09163945) q[2];
sx q[2];
rz(-1.5711725) q[2];
sx q[2];
rz(2.3918652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9190678) q[1];
sx q[1];
rz(-1.7456858) q[1];
sx q[1];
rz(1.2122077) q[1];
x q[2];
rz(2.9321947) q[3];
sx q[3];
rz(-0.87762524) q[3];
sx q[3];
rz(1.587876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.96523607) q[2];
sx q[2];
rz(-2.3945645) q[2];
sx q[2];
rz(-0.61526862) q[2];
rz(-1.2728914) q[3];
sx q[3];
rz(-0.26510173) q[3];
sx q[3];
rz(0.42241514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1212921) q[0];
sx q[0];
rz(-1.2739807) q[0];
sx q[0];
rz(-0.85025382) q[0];
rz(2.0203159) q[1];
sx q[1];
rz(-1.774615) q[1];
sx q[1];
rz(2.3698295) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5284755) q[0];
sx q[0];
rz(-0.69418797) q[0];
sx q[0];
rz(-1.608169) q[0];
x q[1];
rz(1.6408368) q[2];
sx q[2];
rz(-1.3904499) q[2];
sx q[2];
rz(-1.5926306) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.38852019) q[1];
sx q[1];
rz(-0.24986146) q[1];
sx q[1];
rz(-1.3461963) q[1];
x q[2];
rz(0.87987514) q[3];
sx q[3];
rz(-2.7304683) q[3];
sx q[3];
rz(3.0539235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11543342) q[2];
sx q[2];
rz(-1.5449646) q[2];
sx q[2];
rz(0.8405295) q[2];
rz(3.0606411) q[3];
sx q[3];
rz(-2.5637124) q[3];
sx q[3];
rz(2.8988083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42846546) q[0];
sx q[0];
rz(-0.19793887) q[0];
sx q[0];
rz(-1.9848829) q[0];
rz(2.1139862) q[1];
sx q[1];
rz(-1.3778069) q[1];
sx q[1];
rz(-2.3311232) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5549703) q[0];
sx q[0];
rz(-1.3548618) q[0];
sx q[0];
rz(0.76540375) q[0];
rz(-pi) q[1];
rz(0.52105302) q[2];
sx q[2];
rz(-2.5312662) q[2];
sx q[2];
rz(-1.2674171) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15511423) q[1];
sx q[1];
rz(-0.41448516) q[1];
sx q[1];
rz(2.3398967) q[1];
rz(1.719172) q[3];
sx q[3];
rz(-2.1273489) q[3];
sx q[3];
rz(-0.30239964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.345574) q[2];
sx q[2];
rz(-1.8859665) q[2];
sx q[2];
rz(1.5239117) q[2];
rz(2.0412622) q[3];
sx q[3];
rz(-1.1805781) q[3];
sx q[3];
rz(0.17980096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0236459) q[0];
sx q[0];
rz(-1.7272341) q[0];
sx q[0];
rz(-0.9247307) q[0];
rz(-3.0991411) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(1.4412075) q[2];
sx q[2];
rz(-0.40934206) q[2];
sx q[2];
rz(-0.48624292) q[2];
rz(1.6067947) q[3];
sx q[3];
rz(-1.1136342) q[3];
sx q[3];
rz(-1.3200091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
